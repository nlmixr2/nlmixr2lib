#' Converts a linear effect to a log-linear effect
#'
#' @inheritParams addIndirect
#' @family PD
#' @return model converted from linear to log-linear effect
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' readModelDb("PK_2cmt_no_depot") |>
#'   addDirectLin() |>
#'   convertLogLin()
#'
#' readModelDb("PK_2cmt_no_depot") |>
#'   addIndirectLin(stim="out") |>
#'   convertLogLin()
convertLogLin <- function(ui, ek=c("Ik", "Ek"), cc=c("Ce", "Cc")) {
  .ui <- rxode2::assertRxUi(ui)
  cc <- rxode2::assertExists(.ui, cc)
  ek <- rxode2::assertVariableExists(.ui, ek)
  .modelLines <- .replaceMult(.ui$lstExpr,
                              v1=ek, v2=cc,
                              ret=paste0(ek, "*log(", cc, ")"))
  .ui <- rxode2::rxUiDecompress(.ui)
  ## .ui$iniDf <- rbind(.theta,
  ##                    .thetaEk,
  ##                    .thetaErr,
  ##                    .eta)
  if (exists("description", envir=.ui$meta)) {
    rm("description", envir=.ui$meta)
  }
  rxode2::model(.ui) <- .modelLines
  .ui
}
