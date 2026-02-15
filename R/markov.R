#' Create a Markov model for `nlmixr2`
#'
#' @param ... Passed to `createMarkovTransitionMatrix()`
#' @param ignoreProbLt Do not estimate probabilities less than the given probability. Set to 0 to estimate all probabilities with any transition chance, use `estimateZeroTransitions` to give an initial estimate when there is zero probability of transition into the state.
#' @param ignoreProbGt Do not estimate probabilities on a row where any probability is greater than the given probability (treating it as a collector row). Set to 1 to estimate all probabilities with any transition chance away from the state, use `estimateZeroTransitions` to give an initial estimate when there is zero probability of transition from the state.
#' @param transitions Use this manually-created transition matrix rather than an automatically-creating one
#' @returns A template `nlmixr2` model as a character string with `ini()` and `model()` blocks for the Markov model
#' @family Markov models
#' @export
createMarkovModel <- function(..., ignoreProbLt = 0, ignoreProbGt = 1, transitions = NULL) {
  # Check inputs
  checkmate::assert_number(ignoreProbLt, lower = 0, upper = 0.99, na.ok = FALSE)
  checkmate::assert_number(ignoreProbGt, lower = 0.01, upper = 1, na.ok = FALSE)
  if (is.null(transitions)) {
    transitions <- createMarkovTransitionMatrix(...)
  }
  checkmate::assert_matrix(transitions, mode = "double", any.missing = FALSE)
  if (nrow(transitions) != ncol(transitions)) {
    # This is not necessarily required; rectangular transition matrices are possible
    stop("`transitions` matrix must be square")
  } else if (!all(rownames(transitions) == colnames(transitions))) {
    # This is not necessarily required in the case of rectangular transition matrices
    stop("row and column names of `transitions` matrix must be the same")
  }

  allStates <- stats::setNames(rownames(transitions), nm = make.names(rownames(transitions), unique = TRUE))
  modelCallStr <-
    c(
      "# Create the following one-hot encoded columns for previous and current Markov states (this can be done with `createMarkovModelDataset()`)",
      paste0("# For state ", unname(allStates), ": prev", names(allStates), ", cur", names(allStates))
    )

  transitionsCumSum <- lapply(X = stats::setNames(nm = rownames(transitions)), FUN = function(x) transitions[x, ])
  # Drop collector rows
  maskDropCollector <- vapply(X = transitionsCumSum, FUN = function(x) any(x >= ignoreProbGt), FUN.VALUE = NA)
  if (any(maskDropCollector)) {
    modelCallStr <-
      c(
        modelCallStr,
        paste("# No transitions created for collector state,", names(maskDropCollector)[maskDropCollector])
      )
    transitionsCumSum <- transitionsCumSum[!maskDropCollector]
  }
  # Drop ignorably small probabilities
  for (nm in names(transitionsCumSum)) {
    maskKeepHighProb <- transitionsCumSum[[nm]] > ignoreProbLt
    transitionsCumSum[[nm]] <- transitionsCumSum[[nm]][maskKeepHighProb]
  }
  # Validate that each non-collector state still has at least two transitions
  insufficientTransitions <-
    vapply(
      X = transitionsCumSum,
      FUN = function(x) length(x) <= 1L,
      FUN.VALUE = logical(1L)
    )
  if (any(insufficientTransitions)) {
    stop(
      sprintf(
        "After filtering with ignoreProbLt = %g, the following states have insufficient outgoing transitions (<= 1): %s",
        ignoreProbLt,
        paste(names(transitionsCumSum)[insufficientTransitions], collapse = ", ")
      ),
      call. = FALSE
    )
  }
  uiPartsPerState <-
    lapply(
      X = stats::setNames(nm = names(transitionsCumSum)),
      FUN = function(nm) createMarkovModelFromSingleState(transitionsCumSum[nm], stateNames = allStates)
    )

  iniParts <- unlist(sapply(X = uiPartsPerState, FUN = "[", "ini"), use.names = FALSE)
  modelParts <- unlist(sapply(X = uiPartsPerState, FUN = "[", "model"), use.names = FALSE)

  # Build the full log-likelihood
  transitioningStates <- allStates[!maskDropCollector]
  modelParts[length(modelParts) + 1] <- "# Overall Markov model log-likelihood"
  modelParts[length(modelParts) + 1] <-
    sprintf(
      "llMarkov <- %s",
       paste(sprintf(
        "ll%s", names(transitioningStates)
       ), collapse = " + ")
    )
  modelParts[length(modelParts) + 1] <- "ll(err) ~ llMarkov"

  iniFull <- c("ini({", paste(" ", iniParts), "})")
  modelFull <- c("model({", paste(" ", c(modelCallStr, modelParts)), "})")
  metaMarkovStates <- sprintf("markovStates <- %s", paste(deparse(allStates), collapse = " "))
  paste(
    c("function() {", paste(" ", c(metaMarkovStates, iniFull, modelFull)), "}"),
    collapse = "\n"
  )
}

#' Create the parts of a Markov model for transitioning from a single state
#' @param transitionRow a single-element named list with a named vector of all transitions
#' @param stateNames a named vector of state names where the name is the name for use in the model parameters
#' @returns A list with two elements, "ini" and "model", where each element is a character vector of lines of code for the model
#' @keywords Internal
#' @family Markov models
createMarkovModelFromSingleState <- function(transitionRow, stateNames) {
  checkmate::assert_list(transitionRow, len = 1)
  checkmate::assert_named(transitionRow)
  fromState <- names(transitionRow)
  fromStateName <- names(stateNames[stateNames == fromState])
  # The cumulative probability of transitioning between each state
  cumsumRow <- cumsum(transitionRow[[1]])
  stopifnot(length(cumsumRow) > 1)

  # Setup the ini() part

  # The initial value that will go into the ini() block
  iniValue <- rep(NA_real_, length(cumsumRow) - 1)
  ## The first initial value is simply the logit
  iniValue[1] <- logit(cumsumRow[1])
  ## Subsequent initial values are the log of the difference in cumulative
  ## probability from the previous value
  if (length(iniValue) > 1) {
    indices <- 2:length(iniValue)
    iniValue[indices] <- log(diff(cumsumRow[-length(cumsumRow)]))
  }
  toState <- names(cumsumRow)
  toStateName <- names(stateNames[stateNames %in% names(cumsumRow)])
  # Many Markov operations are different for the final state. Give a simple way
  # to index that out of the operations.
  notLastState <- -length(cumsumRow)

  iniParamPrefix <- ifelse(seq_along(iniValue) == 1, "logit", "log")
  iniParams <- sprintf("%s%sto%s", iniParamPrefix, fromStateName, toStateName[notLastState])
  iniUnit <- ifelse(seq_along(iniValue) == 1, "logit probability", "log-logit link difference from prior state")

  # Safely escape state names for use in string literals
  fromStateEscaped <- gsub('"', '\\\\"', fromState[notLastState], fixed = TRUE)
  toStateEscaped <- gsub('"', '\\\\"', toState[notLastState], fixed = TRUE)

  retIni <-
    sprintf(
      '%s <- %g; label("Probability of transition from state %s to %s (%s)")',
      iniParams, signif(iniValue, digits = 4),
      fromStateEscaped, toStateEscaped, iniUnit
    )

  # Setup the model() part
  modelParams <- sprintf("%sto%s", fromStateName, toStateName)
  modelLinkParams <- paste0("link", modelParams) # link function
  modelCumProbParams <- paste0("cumpr", modelParams) # cumulative probability
  modelProbParams <- paste0("pr", modelParams) # probability
  retModel <- character() # The code for the model (as a character vector)

  # Do not estimate a value for the final probability as it's the difference with 100%
  for (idx in seq_len(length(cumsumRow) - 1)) {
    retModel[length(retModel) + 1] <- sprintf('# transition from state "%s" to state "%s"', fromState, toState[idx])
    if (idx == 1) {
      retModel[length(retModel) + 1] <- sprintf("%s <- %s", modelLinkParams[idx], iniParams[idx])
    } else {
      retModel[length(retModel) + 1] <- sprintf("%s <- %s + exp(%s)", modelLinkParams[idx], modelLinkParams[idx-1], iniParams[idx])
    }
    retModel[length(retModel) + 1] <- sprintf('%s <- expit(%s)', modelCumProbParams[idx], modelLinkParams[idx])
  }

  # Probability of each state transition
  retModel[length(retModel) + 1] <- "# Probability of each state transition"
  for (idx in seq_along(cumsumRow)) {
    if (idx == 1) {
      # The initial state
      retModel[length(retModel) + 1] <- sprintf('%s <- %s', modelProbParams[idx], modelCumProbParams[[idx]])
    } else if (idx == length(cumsumRow)) {
      # The final state
      retModel[length(retModel) + 1] <- sprintf('%s <- 1 - %s', modelProbParams[idx], modelCumProbParams[[idx - 1]])
    } else {
      retModel[length(retModel) + 1] <- sprintf("%s <- %s - %s", modelProbParams[idx], modelCumProbParams[[idx]], modelCumProbParams[[idx - 1]])
    }
    retModel[length(retModel)] <-
      sprintf(
        "%s # Probability of transition from state %s to %s",
        retModel[length(retModel)], fromState, toState[idx]
      )
  }

  # Overall probability of state transition when from the current state
  retModel[length(retModel) + 1] <- sprintf("# log-likelihood of any transition from state %s", fromState)
  retModel[length(retModel) + 1] <-
    sprintf(
      'll%s <- prev%s*(%s)',
      fromStateName, fromStateName,
      paste(
        sprintf("cur%s*log(%s)", toStateName, modelProbParams),
        collapse = " + "
      )
    )

  list(
    ini = retIni,
    model = retModel
  )
}

#' Create a Markov transition matrix with probabilities of transitioning between every state
#'
#' @inheritParams createMarkovModelDataset
#' @param estimateZeroTransitions Should transitions that have zero occurrences be estimated? This is done by setting the state to have a single transition.
#' @param estimateZeroTransitionsInitial Should transitions that are only initial states be estimated (ignored if `estimateZeroTransitions = FALSE`)
#' @param ... Ignored
#' @returns A square matrix with row and column names for each state where rows are the prior state and columns are the current state.
#' @family Markov models
#' @export
createMarkovTransitionMatrix <- function(colPrev, colCur, estimateZeroTransitions = FALSE, estimateZeroTransitionsInitial = FALSE, ...) {
  # Create the transition matrix
  if (any(is.na(colPrev))) {
    stop("`colPrev` cannot be `NA`")
  } else if (any(is.na(colCur))) {
    stop("`colCur` cannot be `NA`")
  }
  # find all states in the data
  allStates <- sort(unique(c(colPrev, colCur)))
  if (length(allStates) < 2) {
    stop("Only one state detected, cannot create a nontrivial Markov model")
  }
  transitionCount <-
    matrix(
      NA_real_,
      nrow = length(allStates), ncol = length(allStates),
      dimnames = list(as.character(allStates), as.character(allStates))
    )
  # From state is the row, to state is the column
  for (idx1 in seq_along(allStates)) {
    mask1 <- colPrev %in% allStates[idx1]
    for (idx2 in seq_along(allStates)) {
      mask2 <- colCur %in% allStates[idx2]
      transitionCount[idx1, idx2] <- sum(mask1 & mask2)
    }
  }

  initialStateOnly <- colSums(transitionCount) == 0
  initialStateNames <- names(initialStateOnly)[initialStateOnly]
  collectingStateOnly <- rowSums(transitionCount != 0) == 1
  collectingStateNames <- names(collectingStateOnly)[collectingStateOnly]
  if (any(initialStateOnly)) {
    message("The following state(s) appear only to be initial states: ", paste(initialStateNames, collapse = ", "))
  }
  if (any(collectingStateOnly)) {
    message("The following appear only to be collecting state(s): ", paste(collectingStateNames, collapse = ", "))
  }
  if (estimateZeroTransitions && any(transitionCount %in% 0)) {
    nmStates <- as.character(allStates)
    if (!estimateZeroTransitionsInitial) {
      nmStates <- setdiff(nmStates, initialStateNames)
    }
    for (currentNm in nmStates) {
      currentCount <- transitionCount[, currentNm]
      currentCount[currentCount %in% 0] <- 1
      transitionCount[, currentNm] <- currentCount
    }
  }
  transitionProbs <- transitionCount
  transitionProbsSums <- rowSums(transitionCount)
  for (idx in seq_len(ncol(transitionProbs))) {
    if (transitionProbsSums[[idx]] > 0) {
      transitionProbs[idx, ] <- transitionProbs[idx, ] / transitionProbsSums[[idx]]
    }
  }
  transitionProbs
}

#' Create a Markov model dataset for use in estimating Markov models
#'
#' @param x The object to create a Markov dataset within
#' @param ... Passed to methods
#' @family Markov models
#' @export
createMarkovModelDataset <- function(x, ...) {
  UseMethod("createMarkovModelDataset")
}

#' @describeIn createMarkovModelDataset Create a Markov dataset from a pair of
#'   vectors
#' @param colPrev,colCur Vectors for the previous and current states (may be any
#'   R object that can be coerced to a character string, ordered objects are
#'   usually preferred)
#' @param prefixPrev,prefixCur Column name prefixes for the previous and current
#'   states
#' @returns The data.frame modified by adding one-hot columns for the previous
#'   and current states
#' @export
#' @examples
#' createMarkovModelDataset(c(1, 2, 1), colCur = c(1, 2, 3))
createMarkovModelDataset.default <- function(x, colCur, colPrev = x, prefixPrev = "prev", prefixCur = "cur", ...) {
  checkmate::assert_vector(colPrev)
  checkmate::assert_vector(colCur)
  stopifnot(length(colPrev) == length(colCur))
  # sort() will automatically omit NA values
  allStates <- sort(unique(c(colPrev, colCur)), na.last = NA)
  allStates <- stats::setNames(allStates, nm = make.names(allStates, unique = TRUE))
  ret <- data.frame(prev = colPrev, cur = colCur)
  for (stateIdx in seq_along(allStates)) {
    stateVal <- allStates[[stateIdx]]
    stateNm <- names(allStates[stateIdx])
    stateColPrev <- paste0(prefixPrev, stateNm)
    stateColCur <- paste0(prefixCur, stateNm)
    ret[[stateColPrev]] <- colPrev %in% stateVal
    ret[[stateColCur]] <- colCur %in% stateVal
  }
  ret$prev <- ret$cur <- NULL
  ret
}

#' @export
createMarkovModelDataset.factor <- function(x, colCur, colPrev=x, ...) {
  checkmate::assert_factor(x)
  checkmate::assert_factor(colPrev)
  checkmate::assert_factor(colCur)
  # Factor levels must match (consider handling if one's levels are a strict
  # superset of the other's)
  checkmate::assert_factor(colCur, levels = levels(colPrev))
  createMarkovModelDataset.default(colPrev = x, colCur = colCur, ...)
}

#' @describeIn createMarkovModelDataset Create a Markov dataset from a
#'   data.frame
#' @param colPrev,colCur Column names in the dataset for the previous and
#'   current states (may be any R object that can be coerced to a character
#'   string, ordered objects are usually preferred)
#' @returns The data.frame modified by adding one-hot columns for the previous
#'   and current states
#' @export
createMarkovModelDataset.data.frame <- function(x, colPrev, colCur, ...) {
  checkmate::assert_data_frame(x)
  checkmate::assert_names(names(x), must.include = c(colPrev, colCur))
  cbind(x, createMarkovModelDataset.default(colPrev = x[[colPrev]], colCur = x[[colCur]], ...))
}

#' Convert nlmixr2 simulations to Markov states
#'
#' Simulations are performed per sim.id (if applicable) and per subject id
#' within the simulation. The `sim.id` column must be named exactly that, and
#' the `id` column is detected case-insensitvely.
#'
#' See the Markov vignette (\code{vignette("markov", package = "nlmixr2lib")})
#' for an example.
#'
#' @param ui a model fit or simulation object to simulate from or a data.frame
#'   containing the transition probability columns
#' @param initialState The initial Markov state
#' @param states All states in the model as a character vector
#' @inheritParams createMarkovModelDataset
#' @returns A data.frame with the original data along with the new state columns
#'   (`colPrev` and `colCur`) added
#' @family Markov models
#' @export
simMarkov <- function(ui, initialState, states, colPrev = "previous", colCur = "current", ...) {
  checkmate::assert_choice(initialState, choices = states, null.ok = FALSE)
  checkmate::assert_character(states, min.len = 2, any.missing = FALSE)
  if (is.null(names(states))) {
    names(states) <- make.names(states, unique = TRUE)
  }
  ret <- as.data.frame(ui)
  ret[[colPrev]] <- ret[[colCur]] <- NA
  idCol <- names(ret)[tolower(names(ret)) == "id"]
  if (length(idCol) == 0L) {
    stop("Could not find an ID column in `ui`; expected a column named 'id' (case-insensitive).")
  }
  addSimId <- FALSE
  if (!("sim.id" %in% names(ret))) {
    ret$sim.id <- 1
    addSimId <- TRUE
  }

  transitionPrCols <-
    lapply(
      X = stats::setNames(nm = names(states)),
      FUN = function(x) {
        allPrCols <- stats::setNames(paste0("pr", x, "to", names(states)), nm = names(states))
        # Don't use intersect() to preserve the names
        allPrCols[allPrCols %in% names(ret)]
      }
    )

  for (simCur in unique(ret$sim.id)) {
    # Handle multiple simulations
    maskSim <- ret$sim.id %in% simCur
    for (idCur in unique(ret[[idCol]])) {
      maskId <- maskSim & (ret[[idCol]] %in% idCur)
      initialStateName <- names(states)[states == initialState]
      retId <- simMarkovId(ret[maskId, ], initialState = initialStateName, prCols = transitionPrCols)
      ret[[colPrev]][maskId] <- retId$prev
      ret[[colCur]][maskId] <- retId$cur
    }
  }
  ret[[colPrev]] <- factor(ret[[colPrev]], levels = states)
  ret[[colCur]] <- factor(ret[[colCur]], levels = states)
  if (addSimId) {
    ret$sim.id <- NULL
  }
  ret
}

#' Simulate a Markov model broken down by ID
#' @param data A data.frame for the individual with columns for each of the
#'   state transition probabilities
#' @inheritParams simMarkov
#' @param prCols A named list of probability columns. List names are the
#'   previous state, and list values are a character vector of probability
#'   columns.
#' @returns A data.frame with two columns named "prev" and "cur"
#' @keywords Internal
#' @family Markov models
simMarkovId <- function(data, initialState, prCols) {
  checkmate::assert_data_frame(data)
  # Make sure that the data has all of the columns used as probability columns
  checkmate::assert_names(names(data), must.include = unlist(prCols))
  checkmate::assert_list(prCols, names = "unique")
  collectingStates <- character()
  for (nm in names(prCols)) {
    if (length(prCols[[nm]]) > 0) {
      varName <- sprintf("prCols[['%s']]", nm)
      checkmate::assert_character(prCols[[nm]], names = "unique", .var.name = varName)
      # Verify that each of the probability column specifications is named for a
      # state
      checkmate::assert_names(names(prCols[[nm]]), subset.of = names(prCols), .var.name = varName)
    } else {
      collectingStates <- c(collectingStates, nm)
    }
  }
  ret <- data.frame(prev = rep(NA, nrow(data)), cur = NA)
  prevState <- initialState
  randNums <- stats::runif(n = nrow(data))
  # The state that goes from the final percentage to 100% is the last one in the
  # list
  for (idx in seq_len(nrow(ret))) {
    ret$prev[idx] <- prevState
    # Create the transition matrix row of cumulative sums
    # cumulative probability columns
    availableCols <- prCols[[ret$prev[idx]]]
    transitionVec <- data[idx, availableCols, drop = FALSE]
    # Find the first cumulative probability greater than the randomly selected
    # probability. Protect from floating point issues where the sum could add to
    # a value <1 by always having the final category as a default, last value.
    curStateIdx <- c(which(cumsum(unlist(transitionVec)) > randNums[idx]), length(availableCols))[1]
    ret$cur[idx] <- prevState <- names(availableCols[curStateIdx])
    if ((prevState %in% collectingStates) && (idx < nrow(ret))) {
      # If the simulation enters a collecting state, the rest of the results
      # stay in that state
      ret$prev[(idx+1):nrow(ret)] <- ret$cur[(idx+1):nrow(ret)] <- prevState
      break
    }
  }
  ret
}
