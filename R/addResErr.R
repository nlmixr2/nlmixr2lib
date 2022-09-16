#' Add residual error to a model
#'
#' @param model The model as a function
#' @param reserr character with the type of residual error (currently "add",
#'   "prop" and "add+prop" are accepted)
#' @return The model with residual error modified
#' @examples
#' readModelDb("PK_1cmt") %>% addResErr("add")
#' @export
addResErr <- function(model, reserr) {
  modelUi <- rxode2::rxode(model)
  paramErr <- modelUi$predDf$cond
  if ("rxLinCmt" %in% paramErr) {
    paramErr[paramErr %in% "rxLinCmt"] <- "linCmt()"
  }
  checkmate::expect_character(paramErr, len = 1, min.chars = 1)

  if (is.character(reserr)) {
    if (reserr == "add") {
      newErrLine <- sprintf("%s ~ add(add.err)", paramErr)
      modelUi <- do.call(rxode2::model, list(modelUi, str2lang(newErrLine)))
      modelUi <- rxode2::ini(modelUi, add.err=1)
    } else if (reserr == "prop") {
      newErrLine <- sprintf("%s ~ prop(prop.err)", paramErr)
      modelUi <- do.call(rxode2::model, list(modelUi, str2lang(newErrLine)))
      modelUi <- rxode2::ini(modelUi, prop.err=0.5)
    } else if (reserr == "add+prop") {
      newErrLine <- sprintf("%s ~ add(add.err) + prop(prop.err)", paramErr)
      modelUi <- do.call(rxode2::model, list(modelUi, str2lang(newErrLine)))
      modelUi <- rxode2::ini(modelUi, add.err=1, prop.err=0.5)
    } else {
      cli::cli_abort("unknown residual error")
    }
  } else {
    cli::cli_abort("reserr must be character")
  }
  modelUi$fun
}
