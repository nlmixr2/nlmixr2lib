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
addDirectLin <- function(ui, ek = "Ek", cc = c("Ce", "Cc"), effect = "effect") {
  if (missing(ui)) {
    return(fakeCc(addDirectLin, ek = ek, cc = cc, effect = effect))
  }
  .ui <- rxode2::assertRxUi(ui)
  cc <- rxode2::assertExists(.ui, cc)
  .effectSd <- defaultCombine(effect, "sd")
  rxode2::assertVariableNew(.ui, ek)
  rxode2::assertVariableNew(.ui, effect)
  rxode2::assertVariableNew(.ui, .effectSd)

  .eff <- str2lang(paste0(effect, " <- ", ek, "*", cc))
  .modelLines <- c(list(paste0(ek, " <- u", ek)), .ui$lstExpr, .eff)
  .errLine <- str2lang(paste0(effect, " ~ add(", .effectSd, ")"))

  .ui <- rxode2::rxUiDecompress(.ui)
  if (exists("description", envir = .ui$meta)) {
    rm("description", envir = .ui$meta)
  }
  rxode2::model(.ui) <- .modelLines
  .ui <- .iniAddTheta(.ui, paste0("u", ek), label = paste0("untransformed slope (", ek, ")"))
  # the endpoint goes on by itself so rxode2 creates the residual parameter
  # and decides its condition, err and lower bound
  .ui <- .modelAppend(.ui, list(.errLine))
  .iniAddTheta(.ui, .effectSd, label = paste0("additive error for ", effect))
}
