#' Add direct linear effect with baseline=0
#'
#' @inheritParams addIndirectLin
#' @family PD
#' @return model with direct linear effect added (baseline=0)
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#'
#' # Direct linear model
#' readModelDb("PK_2cmt_no_depot") |>
#'   addDirectLin()
#'
#' # Direct emax model
#' readModelDb("PK_2cmt_no_depot") |>
#'   addDirectLin() |>
#'   convertEmax()
#'
addDirectLin <- function(ui,
                         ek="Ek",
                         cc=c("Ce", "Cc"),
                         effect="effect") {
  if (missing(ui)) {
    return(fakeCc(addDirectLin,
                  ek=ek, cc=cc, effect=effect))
  }
  .ui <- rxode2::assertRxUi(ui)
  cc <- rxode2::assertExists(.ui, cc)
  .effectSd <- paste0(effect, "Sd")
  rxode2::assertVariableNew(.ui, ek)
  rxode2::assertVariableNew(.ui, effect)
  rxode2::assertVariableNew(.ui, .effectSd)

  .eff <- str2lang(paste0(effect, " <- ", ek, "*", cc))
  .modelLines <- c(list(paste0(ek, " <- u", ek)),
                   .ui$lstExpr,
                   .eff,
                   str2lang(paste0(effect, " ~ add(", .effectSd, ")")))

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
  .thetaEk <- .get1theta(ek, .theta1, .ntheta,
                         name=paste0("u", ek),
                         label=paste0("untransformed slope (", ek, ")"))
  .ntheta <- .ntheta + 1
  .thetaErr <-  .get1theta(.effectSd, .theta1, .ntheta,
                           lower=0,
                           label=paste0("additive error for ", effect),
                           name=.effectSd)
  .thetaErr$condition <- effect
  .thetaErr$err <- "add"

  .ui <- rxode2::rxUiDecompress(.ui)
  .ui$iniDf <- rbind(.theta,
                     .thetaEk,
                     .thetaErr,
                     .eta)
  if (exists("description", envir=.ui$meta)) {
    rm("description", envir=.ui$meta)
  }
  rxode2::model(.ui) <- .modelLines
  .ui
}
