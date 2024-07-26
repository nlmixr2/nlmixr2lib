#' Converts first order absorption model to Weibull absorption model
#'
#' @inheritParams addTransit
#' @param wa weibull alpha parameter name
#' @param wb weibull beta parameter name
#' @return model where first order absorption is changed to weibull absorption model
#' @family absorption
#' @export
#' @author Matthew L. Fidler
#' @examples
#' readModelDb("PK_1cmt_des") |>
#'   addWeibullAbs()
addWeibullAbs <- function(ui, ntransit, central = "central",
                          depot = "depot",
                          transit="transit",
                          wa="wa",
                          wb="wb",
                          ka="ka",
                          ktr="ktr") {
  .ui <- rxode2::assertRxUi(ui)
  central <- rxode2::assertCompartmentExists(.ui, central)
  if (!rxode2::testCompartmentExists(.ui, depot)) {
    .ui <- addDepot(.ui, central=central, depot=depot, ka=ka)
    .mv <- rxode2::rxModelVars(.ui)
    warning("'", depot, "' added to model for weibull absroption model", call.=FALSE)
  } else if (rxode2::testCompartmentExists(.ui, paste0(transit, "1"))) {
    .ui <- removeTransit(.ui,
                         central = central,
                         depot = depot, transit=transit,
                         ktr = ktr,
                         ka=ka)
    warning("transit compartments removed for weibull absroption model", call.=FALSE)
  }
  ka <- rxode2::assertVariableExists(.ui, ka)
  rxode2::assertVariableNew(.ui, wa)
  rxode2::assertVariableNew(.ui, wb)
  rxode2::assertVariableNew(.ui, wbk)

  # replace ka*depot with Weibull absorption model in central
  .modelLines <- .ui$lstExpr
  .wb <- paste0("(", wb, "/", wa, ")*(tad0(", depot, ")/", wa,
                ")^(", wb, "-1)*", depot)
  .w <- .whichDdt(.modelLines, central)
  .tmp <- .extractModelLinesAtW(.modelLines, .w)
  .modelLines <- c(.tmp$pre,
                   .replaceMult(.tmp$w,
                           v1=depot, v2=ka,
                           ret=.wb),
                   .tmp$post)

  # replace ka*depot with Weibull absorption model in depot, and add
  # initial estimates
  .w <- .whichDdt(.modelLines, depot)
  .tmp <- .extractModelLinesAtW(.modelLines, .w)
  .modelLines <- c(list(str2lang(paste0(wa, " <- exp(l", wa, ")")),
                        str2lang(paste0(wb, " <- exp(l", wb, ")"))),
                   .tmp$pre,
                   .replaceMult(.tmp$w,
                                v1=depot, v2=ka,
                                ret=.wb),
                   .tmp$post)
  # add parameter estimates
  .tmp <- .getEtaThetaTheta1(.ui)
  .iniDf <- .tmp$iniDf
  .theta <- .tmp$theta
  .theta1 <- .tmp$theta1
  .eta <- .tmp$eta

  if (length(.theta$name) == 0L) {
    .ntheta <- 0
  } else {
    .ntheta <- max(.theta$ntheta)
  }
  .thetawa <- .get1theta(wa, .theta1, .ntheta,
                          label=paste0("Weibull absorption alpha (", wa, ")"))
  .ntheta <- .ntheta + 1

  .thetawb <- .get1theta(wb, .theta1, .ntheta,
                         label=paste0("Weibull absorption beta (", wa, ")"))
  .ntheta <- .ntheta + 1
  .ui <- rxode2::rxUiDecompress(.ui)
  .ui$iniDf <- rbind(.theta,
                     .thetawa,
                     .thetawb,
                     .eta)
  if (exists("description", envir=.ui$meta)) {
    rm("description", envir=.ui$meta)
  }

  # modify model block
  rxode2::model(.ui) <- .modelLines
  .ui
}
