#' Convert linear elimination to quadratic elimination
#'
#' @inheritParams addIndirect
#' @param ek2 quadratic coefficient
#' @family PD
#' @return model with linear effect converted to quadratic effect
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' readModelDb("PK_2cmt_no_depot") |>
#'   addIndirectLin(stim = "out") |>
#'   convertQuad()
#'
#' readModelDb("PK_2cmt_no_depot") |>
#'   addDirectLin() |>
#'   convertQuad()
#'
#' readModelDb("PK_2cmt_no_depot") |>
#'   addEffectCmtLin() |>
#'   convertQuad()
convertQuad <- function(ui, ek = c("Ik", "Ek"), cc = c("Ce", "Cc"), ek2 = "Ek2") {
  .ui <- rxode2::assertRxUi(ui)
  .ui <- rxode2::rxUiDecompress(.ui)
  cc <- rxode2::assertExists(.ui, cc)
  ek <- rxode2::assertVariableExists(.ui, ek)
  .modelLines <- c(list(str2lang(paste0(ek2, " <- u", ek2))),
    .replaceMult(.ui$lstExpr,
      v1 = ek, v2 = cc,
      ret = paste0(ek, "*", cc, "+", ek2, "*", cc, "^2")))
  if (exists("description", envir = .ui$meta)) {
    rm("description", envir = .ui$meta)
  }
  rxode2::model(.ui) <- .modelLines
  .iniAddTheta(.ui, paste0("u", ek2),
    label = paste0("untransformed quadratic slope (", ek2, ")"))
}
