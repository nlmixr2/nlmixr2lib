#' Verify that a value is a valid nlmixr2 compartment name
#'
#' @param x The value to test
#' @returns The value or an error
#' @family Assertions
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

#' @describeIn assertCompartmentName Verify that a value is a valid nlmixr2 variable name
#' @export
assertVariableName <- assertCompartmentName

#' @describeIn assertCompartmentName Verify that a value is a valid nlmixr2 parameter value
#' @export
assertParameterValue <- function(x) {
  checkmate::assertNumeric(
    x,
    len=1,
    any.missing=FALSE,
    finite = TRUE,
    .var.name = paste0(deparse(eval.parent(substitute(substitute(x))), width.cutoff = 500L), collapse = "\n")
  )
}
