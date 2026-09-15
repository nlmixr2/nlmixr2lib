#' Converts first order absorption model to Weibull absorption model
#'
#' @inheritParams addTransit
#' @param wa Weibull alpha parameter name
#' @param wb Weibull beta parameter name
#' @return model where first order absorption is changed to Weibull absorption model
#' @family absorption
#' @export
#' @author Matthew L. Fidler
#' @examples
#' readModelDb("PK_1cmt_des") |>
#'   addWeibullAbs()
addWeibullAbs <- function(
  ui,
  ntransit,
  central = "central",
  depot = "depot",
  transit = "transit",
  wa = "wa",
  wb = "wb",
  ka = "ka",
  ktr = "ktr"
) {
  .ui <- rxode2::assertRxUi(ui)
  central <- rxode2::assertCompartmentExists(.ui, central)
  if (!rxode2::testCompartmentExists(.ui, depot)) {
    .ui <- addDepot(.ui, central = central, depot = depot, ka = ka)
    warning("'", depot, "' added to model for Weibull absorption model", call. = FALSE)
  } else if (rxode2::testCompartmentExists(.ui, paste0(transit, "1"))) {
    .ui <- removeTransit(.ui, central = central, depot = depot, transit = transit, ktr = ktr, ka = ka)
    warning("transit compartments removed for Weibull absorption model", call. = FALSE)
  }
  ka <- rxode2::assertVariableExists(.ui, ka)
  rxode2::assertVariableNew(.ui, wa)
  rxode2::assertVariableNew(.ui, wb)

  # replace ka*depot with Weibull absorption model in central
  .modelLines <- .ui$lstExpr
  .wb <- paste0("(", wb, "/", wa, ")*(tad0(", depot, ")/", wa, ")^(", wb, "-1)*", depot)
  .w <- .whichDdt(.modelLines, central)
  .tmp <- .extractModelLinesAtW(.modelLines, .w)
  .modelLines <- c(.tmp$pre, .replaceMult(.tmp$w, v1 = depot, v2 = ka, ret = .wb), .tmp$post)

  # replace ka*depot with Weibull absorption model in depot, and add
  # initial estimates
  .w <- .whichDdt(.modelLines, depot)
  .tmp <- .extractModelLinesAtW(.modelLines, .w)
  .modelLines <- c(
    list(str2lang(paste0(wa, " <- exp(l", wa, ")")), str2lang(paste0(wb, " <- exp(l", wb, ")"))),
    .tmp$pre,
    .replaceMult(.tmp$w, v1 = depot, v2 = ka, ret = .wb),
    .tmp$post
  )
  # add parameter estimates
  .tmp <- .getEtaTheta(.ui)
  .iniDf <- .tmp$iniDf
  .theta <- .tmp$theta
  .eta <- .tmp$eta
  .tmp <- .dropLines(.ui, .modelLines, .theta, .eta, ka)
  .modelLines <- .tmp$modelLines
  .theta <- .tmp$theta
  .eta <- .tmp$eta
  .ntheta <- .tmp$ntheta

  .ui <- rxode2::rxUiDecompress(.ui)
  # .theta/.eta are subsets of this model's own iniDf (ka dropped), so binding
  # them assumes nothing about its columns; the two new parameters go in
  # through ini().
  .ui$iniDf <- rbind(.theta, .eta)
  if (exists("description", envir = .ui$meta)) {
    rm("description", envir = .ui$meta)
  }

  # modify model block
  rxode2::model(.ui) <- .modelLines
  .ui <- .iniAddTheta(.ui, paste0("l", wa), label = paste0("Weibull absorption alpha (", wa, ")"))
  .iniAddTheta(.ui, paste0("l", wb), label = paste0("Weibull absorption beta (", wb, ")"))
}
