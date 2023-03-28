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

#' Add random effects to a model
#'
#' @param model The model as a function
#' @param eta vector with the parameters to add random effects (sometimes
#'   referred to as inter-individual variability, IIV) on
#' @return The model with eta added to the requested parameters
#' @examples
#' library(rxode2)
#' readModelDb("PK_1cmt") %>% addEta("ka")
#' @export
addEta <- function(model, eta) {
  mod <- model # save to apply everything later
  if (is.character(eta)) {
    # Assign a default value
    eta <- stats::setNames(rep(0.1, length(eta)), eta)
  }
  checkmate::assert_numeric(eta, lower = 0, null.ok = FALSE, min.len = 1)
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
  # Return the function itself or the updated ui
  fun <- do.call(rxode2::ini, iniArgs)$fun
  rxode2::rxode2(mod) <- body(fun)
  mod
}
