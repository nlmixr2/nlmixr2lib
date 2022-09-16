#' Read a model from the nlmixr2 model database
#'
#' @param name The name of the model (must be one of \code{modeldb$name})
#' @return The model as a function
#' @export
readModelDb <- function(name) {
  if (!(name %in% modeldb$name)) {
    stop("'name' not in database")
  } else {
    .fileName <- modeldb$filename[modeldb$name == name]
    if (!file.exists(.fileName)) {
      if (file.exists(system(.fileName, package="nlmixr2lib"))) {
        .fileName <- system(.fileName, package="nlmixr2lib")
      }
    }
    ret <- eval(parse(file = .fileName, keep.source=TRUE))
  }
  ret
}
