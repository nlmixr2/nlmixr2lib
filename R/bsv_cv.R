#' Convert coefficient of variation (percent) to standard deviation
#'
#' @param cv The coefficient of variation (%)
#' @returns The between-subject variabilty on the standard deviation scale
#' @export
cvToBsv <- function(cv) {
  log((cv/100)^2 + 1)
}

#' @describeIn cvToBsv Convert standard deviation of between-subject variability to coefficient of variation
#' @returns The coefficient of variation on the percent scale
#' @export
bsvToCv <- function(bsv) {
  sqrt(exp(bsv) - 1) * 100
}
