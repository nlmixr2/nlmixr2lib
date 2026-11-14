#' Read a model from the nlmixr2 model database
#'
#' @param name The name of the model (must be one of \code{modeldb$name})
#' @return The model as a function
#' @export
#' @examples
#' readModelDb("PK_1cmt")
readModelDb <- function(name) {
  .fileName <- .modelDbFile(name)
  ret <- eval(parse(file = .fileName, keep.source = TRUE))
  msg <- attr(ret, "message")
  if (!is.null(msg)) {
    cli::cli_alert_info(msg, class = "nlmixr2lib_model_message")
  }
  ret
}

# The path of the file of the registered model `name`: the registry's filename
# when it exists as given, else that file within the package.
.modelDbFile <- function(name) {
  .modeldb <- readRDS(system.file("modeldb.rds", package = "nlmixr2lib"))
  if (!(name %in% .modeldb$name)) {
    stop("'name' not in database")
  }
  .fileName <- .modeldb$filename[.modeldb$name == name]
  if (!file.exists(.fileName)) {
    .fileName <- system.file(file.path("modeldb", .fileName), package = "nlmixr2lib")
  }
  .fileName
}
