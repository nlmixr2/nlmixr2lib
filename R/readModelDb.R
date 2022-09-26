#' Read a model from the nlmixr2 model database
#'
#' @param name The name of the model (must be one of \code{modeldb$name})
#' @return The model as a function
#' @export
#' @examples
#' readModelDb("PK_1cmt")
readModelDb <- function(name) {
  if (!(name %in% modeldb$name)) {
    stop("'name' not in database")
  } else {
    # Check for the base file path
    .fileName <- modeldb$filename[modeldb$name == name]
    if (!file.exists(.fileName)) {
      # Check within the package
      .fileName <- gsub("inst[/\\]", "", .fileName)
      if (file.exists(system.file(.fileName, package="nlmixr2lib"))) {
        .fileName <- system.file(.fileName, package="nlmixr2lib")
      }
    }
    ret <- eval(parse(file = .fileName, keep.source=TRUE))
  }
  ret
}
