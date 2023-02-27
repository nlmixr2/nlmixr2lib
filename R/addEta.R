.getVarLhs <- function(model) {
  if (!inherits(model, "rxUi")) {
    .ui <- nlmixr2est::nlmixr2(model)
  } else {
    .ui <- model
  }
  .varLhs <- .ui$varLhs
  if (is.null(.varLhs)) .varLhs <- .ui$getSplitMuModel$pureMuRef
  .varLhs
}

.dropEnv <- c("nonmemData", "etaData", "ipredAtol", "ipredRtol",
              "ipredCompare", "predAtol", "predRtol", "predCompare",
              "thetaMat", "dfSub", "dfObs")
# keep sigma

#' Add random effects to a model
#'
#' @param model The model as a function
#' @param eta vector with the parameters to add random effects (sometimes
#'   referred to as inter-individual variability, IIV) on
#' @return The model with eta added to the requested parameters
#' @examples
#' readModelDb("PK_1cmt") %>% addEta("ka")
#' @export
addEta <- function(model, eta) {
  if (is.character(eta)) {
    # Assign a default value
    eta <- stats::setNames(rep(0.1, length(eta)), eta)
  }
  checkmate::assert_numeric(eta, lower = 0, null.ok = FALSE, min.len = 1)
  .ls <- NULL
  .cls <- NULL
  if (inherits(model, "rxUi")) {
    .model <- rxode2::rxUiDecompress(model)
    .ls <- ls(envir=.model)
    .cls <- class(model)
  }
  # Get the mu-referenced parameter names
  # getSplitMuModel requires nlmixr2est, so the model is parsed from there...
  # This will add the S3 method to allow $getSplitModel to work
  murefNames <- .getVarLhs(model)
  for (currentEta in names(eta)) {
    if (currentEta %in% names(murefNames)) {
      # do nothing
    } else if (currentEta %in% murefNames) {
      # Remap the parameter to the mu-referenced value for modification
      priorEta <- currentEta
      currentEta <- names(murefNames)[murefNames %in% currentEta]
      names(eta)[names(eta) == priorEta] <- currentEta
      cli::cli_alert(sprintf("Adding eta to %s instead of %s due to mu-referencing", currentEta, priorEta))
    }
    model <-
      searchReplace(
        object = model,
        find = currentEta,
        replace = sprintf("%s + eta%s", currentEta, currentEta)
      )
  }
  etaIni <- lapply(X = paste0("eta", names(eta), "~", eta), FUN = base::str2lang)
  iniArgs <-
    append(
      list(model), etaIni
    )
  # Work around rxode2 issue #277
  lotri <- rxode2::lotri
  # Return the parsed ui itself
  .ret <- do.call(rxode2::ini, iniArgs)$fun()
  if (!is.null(.ls)) {
    .ret <- rxode2::rxUiDecompress(.ret)
    .ls2 <- ls(.ret)
    .env <- new.env(parent=emptyenv())
    .env$ignore <- NULL
    .env$keep <- NULL
    lapply(setdiff(.ls, .ls2), function(v) {
      if (v %in% .dropEnv) {
        .env$ignore <- c(.env$ignore, v)
        return(NULL)
      }
      .env$keep <- c(.env$keep, v)
      assign(v, get(v, envir=.model), envir=.ret)
      NULL
    })
    if (length(.env$keep) > 0) {
      cli::cli_alert(sprintf("Kept in model: '%s'",
                             paste(paste0("$",.env$keep), collapse="', '")))
    }
    if (length(.env$ignore) > 0) {
      cli::cli_alert(sprintf("Removed from model: '%s'",
                             paste(paste0("$", .env$ignore), collapse="', '")))
    }
    if (inherits(.model, "raw")) {
      .ret <- rxode2::rxUiCompress(.ret)
    }
    class(.ret) <- .cls
  }
  return(.ret)
}
