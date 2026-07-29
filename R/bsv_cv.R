#' Convert coefficient of variation (percent) to between-subject variance
#'
#' `cvToBsv()` converts a log-normal coefficient of variation (expressed as a
#' percent) to the variance on the log scale, `omega^2 = log((CV/100)^2 + 1)`.
#' `bsvToCv()` is the inverse transform.
#'
#' @param cv The coefficient of variation, expressed as a percent.
#' @returns The between-subject variability on the variance (omega^2) scale.
#' @export
cvToBsv <- function(cv) {
  log((cv / 100)^2 + 1)
}

#' @describeIn cvToBsv Convert log-scale between-subject variance back to a
#'   percent coefficient of variation.
#' @param bsv The between-subject variability on the variance (omega^2) scale.
#' @returns The coefficient of variation on the percent scale.
#' @export
bsvToCv <- function(bsv) {
  sqrt(exp(bsv) - 1) * 100
}
