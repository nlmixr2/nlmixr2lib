#' Get the model from the model library
#'
#' This function gets a model from the available model library
#'
#' @param name character with the name of the model to load (if \code{NULL},
#'   lists all available base models)
#' @inheritParams addEta
#' @inheritParams addResErr
#'
#' @details This is a very first draft just to look at the proof of concept
#'
#' @return The function returns a function the model code (or \code{NULL} if the
#'   \code{model = NULL})
#'
#' @export
#' @examples
#'
#' \dontrun{
#'   modellib(name="PK_1cmt")
#'   modellib(name="PK_1cmt", eta = c("ka", "vc"), reserr = "add")
#'   modellib(name="PK_1cmt", reserr = "add")
#' }
modellib <- function(name=NULL, eta=NULL, reserr=NULL){
  if (is.null(name)) {
    # List available models
    cat(paste0(modeldb$name," (",modeldb$description,")"),sep="\n")
    return(invisible(NULL))
  }
  modr <- readModelDb(name = name)
  if (!is.null(eta)) {
    modr <- addEta(modr, eta = eta)
  }
  if (!is.null(reserr)) {
    modr <- addResErr(modr, reserr = reserr)
  }
  # currently a vector with the model function is returned
  # we could easily extend this, e.g.
  # to console: deparse(modr)
  # to file: writeLines(modr,paste0(model,".r"))
  # directly evaluated: eval(parse(text=modr))
  # to Rstudio's current script: rstudioapi::insertText(paste(modr,"\n"))
  return(modr)
}
