#' Get the MD5 hash of the model library
#' 
#' @export
#' @examples
#' nlmix2libMd5()
nlmix2libMd5 <- function() {
  .Call(`_nlmixr2libMd5`)
}
