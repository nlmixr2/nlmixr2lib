#' @title Get the left handed of the model based on defined initial conditions
#' @param model rxode2 model
#' @return The variable to left handed side of the model return
#' @noRd
#' @author Matthew L. Fidler
.getVarLhs <- function(model) {
  if (!inherits(model, "rxUi")) {
    .ui <- nlmixr2est::nlmixr2(model)
  } else {
    .ui <- model
  }
  .varLhs <- .ui$varLhs
  # getSplitMuModel requires nlmixr2est, so the model is parsed from there...
  # This will add the S3 method to allow $getSplitModel to work
  if (is.null(.varLhs)) {
    .varLhs <- .ui$getSplitMuModel$pureMuRef
  }
  .varLhs
}
