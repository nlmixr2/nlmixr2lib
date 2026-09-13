# Machinery for solving every model in the registry at typical values.
#
# Used by test-modeldb-solve-gate.R.  It exists as a helper rather than
# inline in the test file so that the covariate-derivation rules -- the part
# most likely to need extending as new models land -- are documented in one
# place and can be exercised on their own.
#
# THE PROBE MUST NOT BE DEGENERATE.  Setting every covariate to 0 is the
# obvious shortcut and it is wrong: nearly every model in this registry
# carries an allometric term of the form `(WT / 70)^0.75`, so a zero
# covariate collapses the denominator and manufactures BOTH false zeros and
# absurd non-zeros (a probe that did exactly this reported max(Cc) = 1.35e18
# for one model and a solve failure for another; see
# refharvest/reports/oasweep_PMC12648368.md).  Every covariate value below is
# therefore derived, and the rule that produced it is recorded so a surprising
# probe result can be traced back to the value that caused it.

# Fallback typical values by unit string, used ONLY when neither the model's
# covariateData nor its model() body supplies a reference.  These are
# order-of-magnitude physiological placeholders, not claims about any
# particular paper: the gate they feed asks whether a dose reaches the system
# at all, which is insensitive to the exact value but very sensitive to a
# zero in a denominator.  Nothing here may be 0.
.probeUnitDefaults <- c(
  "kg" = 70, "g" = 70000, "lb" = 154,
  "kg/m^2" = 25, "cm" = 170, "m" = 1.7, "m^2" = 1.73,
  "years" = 40, "year" = 40, "y" = 40, "months" = 480, "month" = 480,
  "weeks" = 2080, "week" = 2080, "days" = 14600, "day" = 14600, "h" = 350000,
  "mL/min" = 100, "mL/min/1.73 m^2" = 90, "L/h" = 6, "mL/h" = 6000,
  "g/L" = 40, "g/dL" = 4, "mg/dL" = 1, "mg/L" = 10,
  "umol/L" = 80, "mmol/L" = 5, "nmol/L" = 100, "nM" = 100, "uM" = 1,
  "U/L" = 30, "IU/L" = 30, "U" = 30,
  "ng/mL" = 10, "ug/mL" = 1, "pg/mL" = 100,
  "mg" = 100, "ug" = 100, "g " = 1, "nmol" = 1, "umol" = 1, "mmol" = 1,
  "%" = 50, "fraction" = 0.5,
  "mg/day" = 100, "mg/kg" = 1, "mg/m^2" = 100,
  "(count)" = 1, "(binary)" = 0, "(categorical)" = 0
)

# Read the leading number out of a reference_category, which is written either
# as a bare number (`reference_category = 0`) or as an annotated string
# (`reference_category = "0 (male)"`).  Returns NA when there is no leading
# number, which is the signal to fall through to the next rule.
.probeLeadingNumber <- function(x) {
  if (is.null(x) || length(x) != 1L) return(NA_real_)
  if (is.numeric(x)) return(if (is.finite(x)) as.numeric(x) else NA_real_)
  if (!is.character(x)) return(NA_real_)
  m <- regmatches(x, regexpr("^\\s*[-+]?[0-9]*\\.?[0-9]+", x))
  if (!length(m)) return(NA_real_)
  suppressWarnings(as.numeric(m))
}

# Walk a model() block looking for the constant a covariate is normalised
# against.  Population-PK models almost always write their reference value
# literally into the equation -- `(WT / 70)^0.75`, `CRCL / 68`, `AGE - 40` --
# so the author's own typical value is recoverable from the source rather than
# guessed at.  Returns NA when no such constant appears.
.probeNormalisingWalk <- function(e, covName, hit) {
  if (is.finite(hit$value) || !is.call(e)) return(invisible(NULL))
  if (length(e) == 3L &&
      (identical(e[[1]], as.name("/")) || identical(e[[1]], as.name("-")) ||
       identical(e[[1]], as.name("+")))) {
    lhs <- e[[2]]
    rhs <- e[[3]]
    if (is.name(lhs) && as.character(lhs) == covName && is.numeric(rhs) &&
        length(rhs) == 1L && is.finite(rhs) && rhs != 0) {
      hit$value <- as.numeric(rhs)
    } else if (identical(e[[1]], as.name("/")) && is.name(rhs) &&
               as.character(rhs) == covName && is.numeric(lhs) &&
               length(lhs) == 1L && is.finite(lhs) && lhs != 0) {
      # `70 / WT` -- the reference sits on the other side
      hit$value <- as.numeric(lhs)
    }
  }
  for (i in seq_along(e)) {
    if (!is.null(e[[i]])) .probeNormalisingWalk(e[[i]], covName, hit)
  }
  invisible(NULL)
}

.probeNormalisingConstant <- function(expr, covName) {
  # An environment rather than a closure variable: reference semantics mean the
  # recursion can record its hit without a superassignment, and the assignment
  # target is named right where the assignment happens.
  hit <- new.env(parent = emptyenv())
  hit$value <- NA_real_
  .probeNormalisingWalk(expr, covName, hit)
  hit$value
}

# Pull the model() block out of an rxUi so the normalising-constant search has
# something to walk.
.probeModelExpr <- function(ui) {
  lst <- try(ui$lstExpr, silent = TRUE)
  if (inherits(lst, "try-error") || is.null(lst)) return(NULL)
  as.call(c(list(as.name("{")), lst))
}

#' Derive a probe value for every covariate a model requires
#'
#' @param ui A decompressed rxUi.
#' @return A list with `values` (named numeric, one per covariate in
#'   `ui$allCovs`) and `why` (named character, the rule that produced each).
#' @noRd
probeCovariateValues <- function(ui) {
  covs <- ui$allCovs
  if (is.null(covs) || !length(covs)) {
    return(list(values = numeric(0), why = character(0)))
  }
  meta <- as.list(ui$meta)
  cd <- meta$covariateData
  body <- .probeModelExpr(ui)
  values <- stats::setNames(rep(NA_real_, length(covs)), covs)
  why <- stats::setNames(rep(NA_character_, length(covs)), covs)
  for (nm in covs) {
    entry <- if (is.list(cd) && !is.null(cd[[nm]]) && is.list(cd[[nm]])) cd[[nm]] else list()

    # 1. An explicit reference category.  This is the right answer for every
    #    binary and categorical covariate: it is the value the paper's typical
    #    subject has.
    v <- .probeLeadingNumber(entry$reference_category)
    if (is.finite(v)) {
      values[[nm]] <- v
      why[[nm]] <- paste0("covariateData$reference_category (",
                          as.character(entry$reference_category)[1], ")")
      next
    }

    # 2. The constant the model itself normalises against.  For a continuous
    #    covariate this is strictly better than any table: it is the
    #    population reference the authors fitted to, read out of their own
    #    equation.
    if (!is.null(body)) {
      v <- .probeNormalisingConstant(body, nm)
      if (is.finite(v)) {
        values[[nm]] <- v
        why[[nm]] <- paste0("normalising constant in model() (", format(v), ")")
        next
      }
    }

    # 3. Any other numeric reference the metadata happens to carry.
    for (field in c("reference_value", "reference", "median", "typical_value")) {
      v <- .probeLeadingNumber(entry[[field]])
      if (is.finite(v)) break
    }
    if (is.finite(v)) {
      values[[nm]] <- v
      why[[nm]] <- paste0("covariateData$", field)
      next
    }

    # 4. A physiological placeholder chosen from the declared units.
    u <- entry$units
    if (is.character(u) && length(u) == 1L && u %in% names(.probeUnitDefaults)) {
      values[[nm]] <- unname(.probeUnitDefaults[[u]])
      why[[nm]] <- paste0("unit default for '", u, "'")
      next
    }

    # 5. Last resort, by declared type.  1 rather than 0 for anything
    #    continuous, because a continuous covariate is far more likely to sit
    #    in a denominator than a binary one.
    ty <- entry$type
    ty <- if (is.character(ty) && length(ty) == 1L) ty else "unknown"
    values[[nm]] <- if (ty %in% c("binary", "categorical")) 0 else 1
    why[[nm]] <- paste0("type default for '", ty, "' (no reference available)")
  }
  list(values = values, why = why)
}

# Time grid, in the model's own time units.  A monoclonal antibody on a
# 0-48 grid looks indistinguishable from a dead model if the units are days.
probeTimeGrid <- function(ui) {
  u <- as.list(ui$meta)$units
  tu <- if (is.list(u) && is.character(u$time) && length(u$time) == 1L) tolower(u$time) else ""
  end <- switch(tu,
                "min" = 720, "minute" = 720, "minutes" = 720,
                "h" = 48, "hr" = 48, "hour" = 48, "hours" = 48,
                "day" = 56, "days" = 56, "d" = 56,
                "week" = 24, "weeks" = 24,
                "month" = 12, "months" = 12,
                48)
  seq(0, end, length.out = 97L)
}

# The compartment a model expects to be dosed into, taken from the registry
# rather than guessed.  Models that declare none are not probeable.
probeDoseCmt <- function(name, db) {
  d <- db$dosing[db$name == name]
  if (!length(d) || is.na(d) || !nzchar(d)) return(NA_character_)
  strsplit(d, ",")[[1]][1]
}

# Build the event table.  Observation records need a compartment or dvid once
# a model has more than one endpoint, otherwise rxode2 raises
# "'dvid'->'cmt' or 'cmt' on observation record" -- which is an event-table
# error, not a model defect, and must not be reported as one.
#
# `form` selects how multi-endpoint observations are tagged.  No single form
# works for every model in the registry, so probeSolveModel() tries them in
# turn: naming the endpoint compartment is the most readable and works for the
# large majority, `dvid` indices cover models whose endpoint names are not
# compartment names, and untagged observations cover the rest.
probeEvents <- function(ui, doseCmt, amt, times, form = "cmt") {
  ev <- rxode2::et(amt = amt, cmt = doseCmt)
  eps <- ui$predDf$cond
  if (length(eps) > 1L && !identical(form, "plain")) {
    if (identical(form, "dvid")) {
      for (i in seq_along(eps)) ev <- rxode2::et(ev, times, dvid = i)
    } else {
      for (ep in eps) ev <- rxode2::et(ev, times, cmt = ep)
    }
  } else {
    ev <- rxode2::et(ev, times)
  }
  ev
}

.probeEventForms <- c("cmt", "dvid", "plain")

# Every column worth comparing between two solves of the same model: the ODE
# states and the model's own computed quantities, minus the stochastic
# columns, which differ run to run and say nothing about whether the dose
# arrived.
.probeValueCols <- function(df) {
  drop <- c("id", "time", "evid", "amt", "cmt", "dvid", "ii", "addl", "dur",
            "rate", "ss", "sim", "ipredSim")
  setdiff(names(df), drop)
}

#' Largest relative disagreement between two solves of the same model
#'
#' Scaled per column so a 1e-4 threshold means the same thing for a
#' concentration in ng/mL and an amount in mg. Columns present in only one of
#' the two solves are ignored: the point is to compare the numbers the two
#' paths both claim to produce.
#'
#' @param a,b Results of `probeSolveModel()`, both with `status == "ok"`.
#' @return A single non-negative number; 0 when the two solves agree.
#' @noRd
probeMaxRelDiff <- function(a, b) {
  cols <- intersect(a$cols, b$cols)
  if (!length(cols) || nrow(a$dosed) != nrow(b$dosed)) return(0)
  rel <- vapply(cols, function(cc) {
    x <- a$dosed[[cc]]
    y <- b$dosed[[cc]]
    sc <- suppressWarnings(max(abs(c(x, y)), na.rm = TRUE))
    if (!is.finite(sc) || sc == 0) return(0)
    d <- suppressWarnings(max(abs(x - y), na.rm = TRUE))
    if (!is.finite(d)) return(0)
    d / sc
  }, numeric(1))
  if (!length(rel) || all(is.na(rel))) return(0)
  max(rel, na.rm = TRUE)
}

# Screening for models at risk of the ODE-to-linCmt conversion defect ------
#
# rxode2 converts an ODE system to `linCmt()` only when it matches a
# depot/central/peripheral topology, and the conversion only loses information
# when some right-hand side carries a term that is not proportional to a
# state. Both conditions are visible in the source, so the candidate set can be
# computed from the model files in a couple of seconds -- no compiling, no
# solving. The screen deliberately errs towards including models the converter
# would reject (it checks the cheap structural preconditions, not the full
# topology), because a false positive costs two solves and a false negative
# costs a silently broken model.

.probeModelBlockOfFile <- function(path) {
  env <- new.env(parent = globalenv())
  fname <- sub("[.]R$", "", basename(path))
  sys.source(path, envir = env)
  if (!exists(fname, envir = env, inherits = FALSE)) return(NULL)
  b <- body(get(fname, envir = env))
  for (i in seq_along(b)) {
    node <- b[[i]]
    if (is.call(node) && identical(node[[1]], as.name("model"))) return(node[[2]])
  }
  NULL
}

# `d/dt(x)` parses as `/`(d, dt(x))
.probeIsDdt <- function(target) {
  is.call(target) && identical(target[[1]], as.name("/")) &&
    identical(target[[2]], as.name("d")) && is.call(target[[3]]) &&
    identical(target[[3]][[1]], as.name("dt"))
}

# Split an expression into signed additive terms, the same decomposition
# rxode2's `.collectAddTerms()` performs.
.probeAddTerms <- function(e, sign = 1L, acc = list()) {
  if (is.call(e) && length(e) == 3L && identical(e[[1]], as.name("+"))) {
    return(.probeAddTerms(e[[3]], sign, .probeAddTerms(e[[2]], sign, acc)))
  }
  if (is.call(e) && length(e) == 3L && identical(e[[1]], as.name("-"))) {
    return(.probeAddTerms(e[[3]], -sign, .probeAddTerms(e[[2]], sign, acc)))
  }
  if (is.call(e) && length(e) == 2L && identical(e[[1]], as.name("-"))) {
    return(.probeAddTerms(e[[2]], -sign, acc))
  }
  if (is.call(e) && length(e) == 2L && identical(e[[1]], as.name("("))) {
    return(.probeAddTerms(e[[2]], sign, acc))
  }
  c(acc, list(list(sign = sign, expr = e)))
}

.probeScreenOneFile <- function(path) {
  mb <- tryCatch(.probeModelBlockOfFile(path), error = function(e) NULL)
  if (is.null(mb)) return(FALSE)
  odes <- list()
  assigns <- list()
  for (i in seq_along(mb)) {
    st <- mb[[i]]
    isAssign <- is.call(st) && length(st) >= 3 &&
      (identical(st[[1]], as.name("<-")) || identical(st[[1]], as.name("=")))
    if (!isAssign) next
    if (.probeIsDdt(st[[2]])) {
      odes[[length(odes) + 1L]] <- list(cmt = as.character(st[[2]][[3]][[2]]), rhs = st[[3]])
    } else if (is.name(st[[2]])) {
      assigns[[length(assigns) + 1L]] <- st[[3]]
    }
  }
  # The converter handles at most four compartments.
  if (!length(odes) || length(odes) > 4L) return(FALSE)
  states <- vapply(odes, function(o) o$cmt, character(1))
  # It also needs an output line of the form `var <- <state> / <expr>`.
  hasOut <- FALSE
  for (r in assigns) {
    if (is.call(r) && length(r) == 3L && identical(r[[1]], as.name("/")) &&
        is.name(r[[2]]) && as.character(r[[2]]) %in% states) {
      hasOut <- TRUE
      break
    }
  }
  if (!hasOut) return(FALSE)
  # Finally, something on a right-hand side that references no state at all:
  # an exogenous input the analytical solution cannot carry.
  for (o in odes) {
    for (tm in .probeAddTerms(o$rhs)) {
      if (is.numeric(tm$expr) && length(tm$expr) == 1L && tm$expr == 0) next
      if (!any(all.vars(tm$expr) %in% states)) return(TRUE)
    }
  }
  FALSE
}

#' Registry models structurally at risk of losing a term to `linCmt()`
#'
#' Enumerates every model file shipped in `inst/modeldb`; nothing is
#' hand-listed, so a model added later is screened without touching this code.
#'
#' @return Character vector of model names.
#' @noRd
linCmtRiskCandidates <- function() {
  root <- system.file("modeldb", package = "nlmixr2lib")
  if (!nzchar(root)) return(character(0))
  files <- list.files(root, pattern = "[.]R$", recursive = TRUE, full.names = TRUE)
  keep <- vapply(files, function(f) isTRUE(tryCatch(.probeScreenOneFile(f),
                                                    error = function(e) FALSE)),
                 logical(1))
  unname(sub("[.]R$", "", basename(files[keep])))
}

# One solve of an already-prepared model. `amt = 0` gives the no-dose control,
# which keeps the two solves identical in every other respect -- same event
# grid, same covariates, same solver settings.
.probeSolveOnce <- function(ui, doseCmt, amt, times, covValues, useLinCmt, form) {
  ev <- probeEvents(ui, doseCmt, amt, times, form)
  args <- list(object = ui, events = ev, returnType = "data.frame",
               addDosing = FALSE, useLinCmt = useLinCmt)
  if (length(covValues)) args$params <- covValues
  suppressWarnings(suppressMessages(do.call(rxode2::rxSolve, args)))
}

# Find an observation-record form this model accepts, and return the dosed
# solve made with it.  Returns a `try-error` carrying the FIRST form's message
# when none works, because that message names the real obstacle.
.probeSolveWithAnyForm <- function(ui, doseCmt, amt, times, covValues, useLinCmt) {
  first <- NULL
  for (form in .probeEventForms) {
    r <- try(.probeSolveOnce(ui, doseCmt, amt, times, covValues, useLinCmt, form),
             silent = TRUE)
    if (!inherits(r, "try-error")) return(list(result = r, form = form))
    if (is.null(first)) first <- r
    # Only the multi-endpoint tagging varies between forms; for a
    # single-endpoint model all three build the same table, so retrying is
    # pointless.
    if (nrow(ui$predDf) <= 1L) break
  }
  list(result = first, form = NA_character_)
}

#' Solve one registry model at typical values, with and without a dose
#'
#' The no-dose solve is the control.  Asking "is the solve identically zero"
#' misses a model whose PK pathway died but whose PD baseline is non-zero;
#' asking "did the dose change anything" cannot.  It is also a difference
#' rather than a ratio, so it cannot pass by dividing one zero by another.
#'
#' @param name Registry model name.
#' @param db The model registry (`modeldb`), passed in so the caller can
#'   subset it.
#' @param useLinCmt Passed through to `rxSolve()`; `FALSE` disables rxode2's
#'   ODE-to-`linCmt()` optimisation.
#' @return A list with `status` (one of "ok", "skip", "error"), and on success
#'   `doseEffect` (the largest absolute difference the dose made to any
#'   model quantity), `dosed`/`undosed` data frames, and the covariate
#'   provenance.
#' @noRd
probeSolveModel <- function(name, db, useLinCmt = TRUE, amt = 100) {
  out <- list(name = name, status = "error", message = NA_character_)
  ui <- try(rxode2::rxUiDecompress(rxode2::as.rxUi(nlmixr2lib::readModelDb(name))),
            silent = TRUE)
  if (inherits(ui, "try-error")) {
    out$message <- paste("could not load:", conditionMessage(attr(ui, "condition")))
    return(out)
  }
  doseCmt <- probeDoseCmt(name, db)
  if (is.na(doseCmt)) {
    out$status <- "skip"
    out$message <- "no dosing compartment declared in the registry"
    return(out)
  }
  ui <- try(rxode2::zeroRe(ui), silent = TRUE)
  if (inherits(ui, "try-error")) {
    out$message <- paste("zeroRe() failed:", conditionMessage(attr(ui, "condition")))
    return(out)
  }
  cov <- probeCovariateValues(ui)
  times <- probeTimeGrid(ui)
  attempt <- .probeSolveWithAnyForm(ui, doseCmt, amt, times, cov$values, useLinCmt)
  dosed <- attempt$result
  if (inherits(dosed, "try-error")) {
    out$message <- paste("solve failed:", conditionMessage(attr(dosed, "condition")))
    out$covWhy <- cov$why
    return(out)
  }
  out$eventForm <- attempt$form
  undosed <- try(.probeSolveOnce(ui, doseCmt, 0, times, cov$values, useLinCmt,
                                 attempt$form),
                 silent = TRUE)
  if (inherits(undosed, "try-error")) {
    out$message <- paste("control solve failed:",
                         conditionMessage(attr(undosed, "condition")))
    out$covWhy <- cov$why
    return(out)
  }
  cols <- intersect(.probeValueCols(dosed), .probeValueCols(undosed))
  cols <- cols[vapply(cols, function(cc) is.numeric(dosed[[cc]]), logical(1))]
  if (!length(cols) || nrow(dosed) != nrow(undosed)) {
    out$status <- "skip"
    out$message <- "no comparable numeric output columns"
    return(out)
  }
  diffs <- vapply(cols, function(cc) {
    d <- abs(dosed[[cc]] - undosed[[cc]])
    if (all(is.na(d))) return(NA_real_)
    max(d, na.rm = TRUE)
  }, numeric(1))
  out$status <- "ok"
  out$doseCmt <- doseCmt
  out$doseEffect <- if (all(is.na(diffs))) NA_real_ else max(diffs, na.rm = TRUE)
  out$nonFinite <- any(!is.finite(as.matrix(dosed[, cols, drop = FALSE])))
  out$dosed <- dosed
  out$cols <- cols
  out$covWhy <- cov$why
  out$covValues <- cov$values
  out
}
