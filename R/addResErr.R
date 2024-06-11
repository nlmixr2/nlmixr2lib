#' Add residual error to a model
#'
#' @details For \code{reserr}, the parameter will be named with the dependent
#'   variable from the model as a prefix.  For example, if the dependent
#'   variable in the model is \code{Cc}, the parameter name for \code{propSd}
#'   will become \code{CcpropSd}.
#'
#' @param model The model as a function
#' @param reserr The type or types of residual error (currently \code{"addSd"},
#'   \code{"propSd"}, and \code{"lnormSd"} are accepted)
#' @return The model with residual error modified
#' @examples
#' library(rxode2)
#' readModelDb("PK_1cmt") |> addResErr("addSd")
#' readModelDb("PK_1cmt") |> addResErr("lnormSd")
#' readModelDb("PK_1cmt") |> addResErr(c("addSd", "propSd"))
#' @export
addResErr <- function(model, reserr) {
  mod <- model
  modelUi <- rxode2::rxode(model)
  paramErr <- modelUi$predDf$cond
  if ("rxLinCmt" %in% paramErr) {
    paramErr[paramErr %in% "rxLinCmt"] <- "linCmt()"
  }
  checkmate::assert_character(paramErr, len = 1, min.chars = 1)
  errFunMap <-
    c(
      addSd = "add(%s)",
      propSd = "prop(%s)",
      lnormSd = "lnorm(%s)"
    )
  defaultIniEst <- c(addSd = 1, propSd = 0.5, lnormSd = 0.5)
  # Confirm that the code is up to date
  stopifnot(length(errFunMap) == length(defaultIniEst))
  stopifnot(all(names(errFunMap) %in% names(defaultIniEst)))
  if (is.character(reserr)) {
    checkmate::assert_subset(reserr, choices = names(defaultIniEst), empty.ok = FALSE)
    reserr <- defaultIniEst[reserr]
  } else if (!is.numeric(reserr)) {
    cli::cli_abort("reserr must be a character string or a named numeric vector")
  }
  checkmate::assert_numeric(reserr, min.len = 1, lower = 0, finite = TRUE, any.missing = FALSE)
  checkmate::assert_names(names(reserr), subset.of = names(errFunMap))

  newErrLineRhs <-
    paste(
      sprintf(errFunMap[names(reserr)], paste0(paramErr, names(reserr))),
      collapse = " + "
    )
  newErrLine <- sprintf("%s ~ %s", paramErr, newErrLineRhs)
  newIniEst <- reserr
  names(newIniEst) <- paste0(paramErr, names(reserr))

  # Update the model with the new residual error line and the new initial
  # estimates
  modelUi <- do.call(rxode2::model, list(modelUi, str2lang(newErrLine)))
  modelUi <- rxode2::ini(modelUi, newIniEst)
  # Return the model function or ui with props
  rxode2::rxode2(mod) <- modelUi$fun
  mod
}
