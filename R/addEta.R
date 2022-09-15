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
    eta <- setNames(rep(0.1, length(eta)), eta)
  }
  checkmate::expect_numeric(eta, lower = 0, null.ok = FALSE, min.len = 1)
  # Get the mu-referenced parameter names
  murefNames <- rxode2::rxode2(model)$getSplitMuModel$pureMuRef
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
  etaIni <- lapply(X = paste0("eta", names(eta), "~", eta), FUN = as.formula)
  iniArgs <-
    append(
      list(model), etaIni
    )
  # Work around rxode2 issue #277
  lotri <- rxode2::lotri
  # Return the function itself
  do.call(rxode2::ini, iniArgs)$fun
}
