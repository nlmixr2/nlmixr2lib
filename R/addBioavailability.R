#' Add bioavailability to a model
#'
#' The parameter name in the model will be `paste0("lf", compartment)` if
#' estimated on the log scale or `paste0("tvf", compartment)` if estimated on the
#' linear scale, for example, "lfdepot".  The parameter used in the model will be named
#'
#' @inheritParams addDepot
#' @param compartment The compartment to add bioavailability to (typically
#'   "depot")
#' @param iniValue The initial value to set (on the linear scale)
#' @param logScale Should estimation be performed on the log scale?
#' @returns The model with bioavailability added to the `compartment`
#' @examples
#' readModelDb("PK_1cmt_des") |>
#'   addBioavailability(iniValue = 1)
#' @export
addBioavailability <- function(model, compartment = "depot", iniValue, logScale = TRUE) {
  innerModel <- rxode2::assertRxUi(model)
  assertCompartmentName(compartment)
  assertParameterValue(iniValue)
  checkmate::assertLogical(logScale, any.missing = FALSE, len = 1)

  if (logScale) {
    fparam <- paste0("lf", compartment)
    fcompartmentModel <- paste0("f", compartment, " <- exp(", fparam, ")")
    iniValue <- log(iniValue)
  } else {
    fparam <- paste0("tvf", compartment)
    fcompartmentModel <- paste0("f", compartment, " <- ", fparam, "")
  }

  # Find the line to add the f(compartment) after
  lincmtLine <-
    vapply(
      X = rxode2::model(innerModel)[[2]],
      FUN = rxode2::.matchesLangTemplate,
      template = str2lang("linCmt()"),
      FUN.VALUE = TRUE
    )
  odeLine <-
    vapply(
      X = rxode2::model(innerModel)[[2]],
      FUN = rxode2::.matchesLangTemplate,
      template = str2lang(sprintf("d/dt(%s) <- .", compartment)),
      FUN.VALUE = TRUE
    )

  browser()
  stop()

  modelExtract <- rxode2::model(innerModel)

  if (any(lincmtLine)) {
    appendLine <- max(which(lincmtLine))
  } else if (any(odeLine)) {
    appendLine <- max(which(odeLine))
  } else {
    # If unsure, go to the end
    appendLine <- length(modelExtract[[2]])
  }
  model_mod <-

  rxode2::model(innerModel) <-
    c(
      str2lang(fcompartmentModel),
      modelExtract[[2]][1:appendLine],
      str2lang(paste0("f(", compartment, ") <- f", compartment)),
      modelExtract[[2]][(appendLine + 1):length(modelExtract[[2]])]
    )

  fIni <- str2lang(paste0(fparam, " <-", iniValue))
  modelUpdate <- rxode2::ini(modelUpdate, fIni)

  modelUpdate
}
