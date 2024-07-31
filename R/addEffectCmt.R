#' Add effect compartment
#'
#' @inheritParams addIndirect
#' @param ke0 This is the effect compartment keo rate
#' @param ce This is the concentration in the effect compartment
#' @return a model with an effect compartment attached
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' readModelDb("PK_2cmt_no_depot") |>
#'    addEffectCmtLin()
#'
#' # Can also be changed to the more typical Emax with constant (estimated) baselie
#' readModelDb("PK_2cmt_no_depot") |>
#'   addEffectCmtLin() |>
#'   convertEmaxHill() |>
#'   addBaselineConst()
#'
addEffectCmtLin <- function(ui,
                            ke0="ke0",
                            cc="Cc",
                            ce="Ce",
                            ek="Ek",
                            effect="effect") {
  if (missing(ui)) {
    return(fakeCc(addEffectCmtLin, ke0=ke0, cc=cc, ce=ce, ek=ek, effect=effect))
  }
  .effectSd <- defaultCombine(effect, "sd")
  .ui <- rxode2::assertRxUi(ui)
  .ui <- rxode2::rxUiDecompress(.ui)
  rxode2::assertExists(.ui, cc)
  rxode2::assertCompartmentNew(.ui, ce)
  rxode2::assertCompartmentNew(.ui, ek)
  rxode2::assertVariableNew(.ui, effect)
  rxode2::assertVariableNew(.ui, .effectSd)
  .ce <- str2lang(paste0("d/dt(", ce, ") <- ", ke0, "*(", cc, "-", ce, ")"))
  .ef <- str2lang(paste0("effect <- ", ce, "*", ek))
  .err <- str2lang(paste0("effect ~ add(", .effectSd, ")"))
  .modelLines <- c(list(str2lang(paste0(ke0, "<- exp(l", ke0, ")")),
                        str2lang(paste0(ek, "<- u", ek))),
                   .ui$lstExpr,
                   list(.ce,
                        .ef,
                        .err))

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

  .thetaKe0 <- .get1theta(ke0, .theta1, .ntheta,
                          label=paste0("effect compartment rate (", ke0 , ")"))
  .ntheta <- .ntheta + 1

  .thetaEk <- .get1theta(ek, .theta1, .ntheta,
                         name=paste0("u", ek),
                         label=paste0("untransformed linear slope (",
                                      ek, ")"))
  .ntheta <- .ntheta + 1

  .thetaErr <-  .get1theta(.effectSd, .theta1, .ntheta,
                           lower=0,
                           label=paste0("additive error for ", effect),
                           name=.effectSd)
  .thetaErr$condition <- effect
  .thetaErr$err <- "add"
  .ui$iniDf <- rbind(.theta,
                     .thetaKe0,
                     .thetaEk,
                     .thetaErr,
                     .eta)
  if (exists("description", envir=.ui$meta)) {
    rm("description", envir=.ui$meta)
  }
  rxode2::model(.ui) <- .modelLines
  .ui
}
