#' Convert a kin/kout indirect response model to R0 and kout
#'
#' This replaces the kin/kout parametrization to the R0 and kout parametrization
#'
#' @param ui a rxode2 user function
#' @param kin the kin variable (by default is "kin")
#' @param kout the kout variable (by default is "kout")
#' @param R the compartment variable (by default is "R")
#' @param R0 the R0 variable (by default is "R0")
#' @return a model where the estimated kin is changed to the estimated R0
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' addIndirect(stim="in") |> convertKinR0()
convertKinR0 <- function(ui,
                         kin="kin",
                         kout="kout",
                         R="R",
                         R0="R0") {
  .ui <- rxode2::assertRxUi(ui)
  kin <- rxode2::assertVariableExists(.ui, kin)
  rxode2::assertVariableNew(.ui, R0)
  R <- rxode2::assertCompartmentExists(.ui, R)
  .tmp <- .getEtaThetaTheta1(.ui)
  .iniDf <- .tmp$iniDf
  .theta <- .tmp$theta
  .theta1 <- .tmp$theta1
  .eta <- .tmp$eta
  if(length(.theta$ntheta) == 0) {
    .ntheta <- 0
  } else {
    .ntheta <- max(.theta$ntheta)
  }
  .modelLines <- .ui$lstExpr
  .w <- .whichDdt(.modelLines, R, start="", end="(0)")
  if (!length(.w)) {
    stop(paste0("the model does not have the expected ",
                R, "(0) expression"),
         call.=FALSE)
  }
  .tmp <- .extractModelLinesAtW(.modelLines, .w)
  if (!identical(.tmp$w[[3]], str2lang(paste0(kin, "/", kout)))) {
    stop(paste0("the model does not have the expected ",
                R, "(0) <- ", kin,"/",kout, " expression"),
         call.=FALSE)
  }
  .modelLines <- c(
    str2lang(paste0(R0, "<- u", R0)),
    .tmp$pre,
    list(str2lang(paste0(R, "(0) <- ", R0))),
    .tmp$post)

  .w <- .whichDdt(.modelLines, R)
  if (length(.w) != 1L) {
    stop(paste0("the model does not have the expected d/dt(",
                R, ") expression"),
         call.=FALSE)
  }
  .tmp <- .extractModelLinesAtW(.modelLines, .w)
  .tmp$w <- searchReplaceHelper(.tmp$w, str2lang(kin), str2lang(paste0(kout, "*", R0)))
  .modelLines <- c(.tmp$pre, list(.tmp$w), .tmp$post)

  .tmp <- .dropLines(.ui, .modelLines, .theta, .eta, kin)
  .modelLines <- .tmp$modelLines
  .theta <- .tmp$theta
  .eta <- .tmp$eta

  .thetaR0 <- .get1theta(R0,  .theta1, .ntheta,
                         name=paste0("u", R0),
                         label=paste0("untransformed baseline (",
                                      R0, ")"))
  .ui <- rxode2::rxUiDecompress(.ui)
  .ui$iniDf <- rbind(.theta,
                     .thetaR0,
                     .eta)
  if (exists("description", envir=.ui$meta)) {
    rm("description", envir=.ui$meta)
  }
  rxode2::model(.ui) <- .modelLines
  .ui
}
