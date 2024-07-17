#' Verify that a value is a valid nlmixr2 compartment name
#'
#' @param x The compartment name to test
#' @returns The compartment name or an error
#' @export
assertCompartmentName <- function(x) {
  checkmate::assertCharacter(
    x,
    pattern = "^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$",
    len = 1,
    any.missing = FALSE,
    min.chars = 1,
    .var.name = paste0(deparse(eval.parent(substitute(substitute(x))), width.cutoff = 500L), collapse = "\n")
  )
}
