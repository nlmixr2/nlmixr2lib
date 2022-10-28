#' Add residual error to a model
#'
#' @param model The model as a function
#' @param reserr character with the type of residual error (currently "addSd",
#'   "propSd" and "addSd+propSd" are accepted)
#' @return The model with residual error modified
#' @examples
#' readModelDb("PK_1cmt") %>% addResErr("addSd")
#' @export
addResErr <- function(model, reserr) {
  modelUi <- rxode2::rxode(model)
  paramErr <- modelUi$predDf$cond
  if ("rxLinCmt" %in% paramErr) {
    paramErr[paramErr %in% "rxLinCmt"] <- "linCmt()"
  }
  checkmate::expect_character(paramErr, len = 1, min.chars = 1)

  if (is.character(reserr)) {
    if (reserr == "addSd") {
      newErrLine <- sprintf("%s ~ add(addSd)", paramErr)
      modelUi <- do.call(rxode2::model, list(modelUi, str2lang(newErrLine)))
      modelUi <- rxode2::ini(modelUi, addSd=1)
    } else if (reserr == "propSd") {
      newErrLine <- sprintf("%s ~ prop(propSd)", paramErr)
      modelUi <- do.call(rxode2::model, list(modelUi, str2lang(newErrLine)))
      modelUi <- rxode2::ini(modelUi, propSd=0.5)
    } else if (reserr == "addSd+propSd") {
      newErrLine <- sprintf("%s ~ add(addSd) + prop(propSd)", paramErr)
      modelUi <- do.call(rxode2::model, list(modelUi, str2lang(newErrLine)))
      modelUi <- rxode2::ini(modelUi, addSd=1, propSd=0.5)
    } else {
      cli::cli_abort("unknown residual error")
    }
  } else {
    cli::cli_abort("reserr must be character")
  }
  modelUi$fun
}
