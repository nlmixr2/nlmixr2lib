#' Read a model from the nlmixr2 model database
#'
#' @param name The name of the model (must be one of \code{modeldb$name})
#' @return The model as a function
#' @export
#' @examples
#' readModelDb("PK_1cmt")
readModelDb <- function(name) {
  .modeldb <- qs::qread(system.file("modeldb.qs", package="nlmixr2lib"))
  if (!(name %in% .modeldb$name)) {
    stop("'name' not in database")
  } else {
    # Check for the base file path
    .fileName <- .modeldb$filename[.modeldb$name == name]
    if (!file.exists(.fileName)) {
      # Check within the package
      .fileName <- system.file(file.path("modeldb", .fileName), package = "nlmixr2lib")
    }
    ret <- eval(parse(file = .fileName, keep.source = TRUE))
    msg <- attr(ret, "message")
    if (!is.null(msg)) {
      cli::cli_alert_info(msg, class = "nlmixr2lib_model_message")
    }
  }
  ret
}
