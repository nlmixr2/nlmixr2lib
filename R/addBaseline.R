#' add the baseline information to model lines
#'
#' @param ui rxode2 ui
#' @param effect the variable to add the baseline to
#' @param eb baseline string expression
#' @return model lines
#' @family PD
#' @noRd
#' @author Matthew L. Fidler
.addBaseline <- function(ui, effect = "effect",
                         eb = "Eb") {
  .modelLines <- ui$lstExpr
  .w <- .whichDdt(.modelLines, effect, start = "", end = "")
  .tmp <- .extractModelLinesAtW(.modelLines, .w)
  .tmp$w <- list(str2lang(paste0(deparse1(.tmp$w), "+", eb)))
  c(.tmp$pre,
    .tmp$w,
    .tmp$post)
}

#' Add an estimated baseline constant
#'
#' @inheritParams addIndirect
#' @param eb baseline constant parameter
#' @return model with baseline constant
#' @family PD
#' @export
#' @author Matthew L. Fidler
#' @examples
#' readModelDb("PK_2cmt_no_depot") |>
#'   addDirectLin() |>
#'   convertQuad() |>
#'   addBaselineConst()
addBaselineConst <- function(ui, effect = "effect", eb = "Eb") {
  .ui <- rxode2::assertRxUi(ui)
  .ui <- rxode2::rxUiDecompress(.ui)
  rxode2::assertVariableExists(.ui, effect)
  rxode2::assertVariableNew(.ui, eb)
  .modelLines <- c(list(str2lang(paste0(eb, "<- u", eb))),
    .addBaseline(.ui, effect = effect, eb = eb))

  if (exists("description", envir = .ui$meta)) {
    rm("description", envir = .ui$meta)
  }
  rxode2::model(.ui) <- .modelLines
  .iniAddTheta(.ui, paste0("u", eb),
    label = paste0("untransformed constant baseline (", eb, ")"))
}
#' Add baseline that decays exponential with time
#'
#' @inheritParams addIndirect
#' @inheritParams addBaselineConst
#' @param time the time or other variable used for baseline decay
#' @param kb the first order baseline decay constant
#' @return model with baseline constant
#' @export
#' @family PD
#' @author Matthew L. Fidler
#' @examples
#' readModelDb("PK_2cmt_no_depot") |>
#'   addDirectLin() |>
#'   convertQuad() |>
#'   addBaselineExp()
addBaselineExp <- function(ui, effect = "effect", eb = "Eb",
                           time = "time", kb = "kb") {
  .ui <- rxode2::assertRxUi(ui)
  .ui <- rxode2::rxUiDecompress(.ui)
  rxode2::assertVariableExists(.ui, effect)
  rxode2::assertVariableNew(.ui, eb)
  rxode2::assertVariableNew(.ui, kb)
  .modelLines <- c(list(str2lang(paste0(eb, "<- u", eb)),
    str2lang(paste0(kb, "<- exp(l", kb, ")"))),
  .addBaseline(.ui, effect = effect,
    eb = paste0(eb, "*exp(-", kb, "*",
      time,
      ")")))

  if (exists("description", envir = .ui$meta)) {
    rm("description", envir = .ui$meta)
  }
  rxode2::model(.ui) <- .modelLines
  .ui <- .iniAddTheta(.ui, paste0("u", eb),
    label = paste0("untransformed constant baseline (", eb, ")"))
  .iniAddTheta(.ui, paste0("l", kb),
    label = paste0("baseline time-decay constant (", kb, ")"))
}

#' Add baseline that decays exponential with time
#'
#' @inheritParams addIndirect
#' @inheritParams addBaselineExp
#' @param time the time or other variable used for baseline decay
#' @param kb the first order baseline decay constant
#' @return model with baseline constant
#' @export
#' @family PD
#' @author Matthew L. Fidler
#' @examples
#' readModelDb("PK_2cmt_no_depot") |>
#'   addDirectLin() |>
#'   convertQuad() |>
#'   addBaseline1exp()
addBaseline1exp <- function(ui, effect = "effect", eb = "Eb",
                            time = "time", kb = "kb") {
  .ui <- rxode2::assertRxUi(ui)
  .ui <- rxode2::rxUiDecompress(.ui)
  rxode2::assertVariableExists(.ui, effect)
  rxode2::assertVariableNew(.ui, eb)
  rxode2::assertVariableNew(.ui, kb)
  .modelLines <- c(list(str2lang(paste0(eb, "<- u", eb)),
    str2lang(paste0(kb, "<- exp(l", kb, ")"))),
  .addBaseline(.ui, effect = effect,
    eb = paste0(eb, "*(1-exp(-", kb, "*",
      time,
      "))")))

  if (exists("description", envir = .ui$meta)) {
    rm("description", envir = .ui$meta)
  }
  rxode2::model(.ui) <- .modelLines
  .ui <- .iniAddTheta(.ui, paste0("u", eb),
    label = paste0("untransformed constant baseline (", eb, ")"))
  .iniAddTheta(.ui, paste0("l", kb),
    label = paste0("baseline time-decay constant (", kb, ")"))
}

#' Add an estimated baseline linear constant
#'
#' @inheritParams addBaselineExp
#' @return model with baseline linear constant
#' @family PD
#' @export
#' @author Matthew L. Fidler
#' @examples
#' readModelDb("PK_2cmt_no_depot") |>
#'   addDirectLin() |>
#'   convertQuad() |>
#'   addBaselineLin()
addBaselineLin <- function(ui, effect = "effect", eb = "Eb",
                           time = "time") {
  .ui <- rxode2::assertRxUi(ui)
  .ui <- rxode2::rxUiDecompress(.ui)
  rxode2::assertVariableExists(.ui, effect)
  rxode2::assertVariableNew(.ui, eb)
  .modelLines <- c(list(str2lang(paste0(eb, "<- u", eb))),
    .addBaseline(.ui, effect = effect, eb = paste0(eb, "*", time)))
  if (exists("description", envir = .ui$meta)) {
    rm("description", envir = .ui$meta)
  }
  rxode2::model(.ui) <- .modelLines
  .iniAddTheta(.ui, paste0("u", eb),
    label = paste0("untransformed constant baseline (", eb, ")"))
}
