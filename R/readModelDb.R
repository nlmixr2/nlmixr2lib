#' Read a model from the nlmixr2 model database
#'
#' @param name The name of the model (must be one of \code{modeldb$name})
#' @return The model as a function
#' @export
readModelDb <- function(name) {
  if (!(name %in% modeldb$name)) {
    stop("'name' not in database")
  } else {
    ret <- eval(parse(file = modeldb$filename[modeldb$name == name]))
  }
  ret
}
