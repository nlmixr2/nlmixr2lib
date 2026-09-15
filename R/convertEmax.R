#' Convert linear effect to Emax effect
#'
#' @param ui rxode2 model
#'
#' @param emax Emax parameter
#' @param ec50 EC50 parameter
#' @param imax Imax parameter used when input model contains "Ik"
#'   instead of "Ek"
#' @param ic50 IC50 parameter used when input model contains "Ik"
#'   instead of "Ek"
#' @inheritParams addIndirectLin
#'
#' @family PD
#'
#' @return Model with the linear effect converted to an Emax effect
#'
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' readModelDb("PK_2cmt_no_depot") |>
#'   addIndirectLin(stim = "in") |>
#'   convertEmax()
#'
#' # When emax=1
#' readModelDb("PK_2cmt_no_depot") |>
#'   addIndirectLin(stim = "in") |>
#'   convertEmax(emax = 1)
#'
convertEmax <- function(
  ui,
  emax = "Emax",
  ec50 = "EC50",
  imax = "Imax",
  ic50 = "IC50",
  ek = c("Ik", "Ek"),
  cc = c("Ec", "Cc")
) {
  .ui <- rxode2::assertRxUi(ui)
  cc <- rxode2::assertExists(.ui, cc)
  ek <- rxode2::assertVariableExists(.ui, ek)
  if (ek == "Ik") {
    emax <- imax
    ec50 <- ic50
  }
  if (inherits(emax, "character")) {
    rxode2::assertVariableNew(.ui, emax)
    .emaxMult <- paste0(emax, "*")
  } else if (is.numeric(emax) && emax == 1.0) {
    .emaxMult <- ""
  } else {
    stop("'", emax, "' not specified correctly", call. = FALSE)
  }
  rxode2::assertVariableNew(.ui, ec50)
  .modelLines <- .replaceMult(.ui$lstExpr, v1 = ek, v2 = cc, ret = paste0(.emaxMult, cc, "/(", cc, "+", ec50, ")"))
  .tmp <- .getEtaTheta(.ui)
  .iniDf <- .tmp$iniDf
  .theta <- .tmp$theta
  .eta <- .tmp$eta
  .tmp <- .dropLines(.ui, .modelLines, .theta, .eta, ek)
  .modelLines <- .tmp$modelLines
  .theta <- .tmp$theta
  .eta <- .tmp$eta

  if (length(.theta$ntheta) == 0) {
    .ntheta <- 0
  } else {
    .ntheta <- max(.theta$ntheta)
  }

  .emaxLine <- NULL
  if (.emaxMult != "") {
    .emaxLine <- list(str2lang(paste0(emax, "<- exp(l", emax, ")")))
  }

  .ui <- rxode2::rxUiDecompress(.ui)
  # .theta/.eta are subsets of this model's own iniDf (ek dropped), so binding
  # them assumes nothing about its columns; the new parameters go in via ini().
  .ui$iniDf <- rbind(.theta, .eta)
  if (exists("description", envir = .ui$meta)) {
    rm("description", envir = .ui$meta)
  }
  rxode2::model(.ui) <- c(.emaxLine, list(str2lang(paste0(ec50, "<- exp(l", ec50, ")"))), .modelLines)
  if (!is.null(.emaxLine)) {
    .ui <- .iniAddTheta(.ui, paste0("l", emax), label = paste0("Maximum effect (", emax, ")"))
  }
  .iniAddTheta(.ui, paste0("l", ec50), label = paste0("Concentration of 50% ", emax, " (", emax, ")"))
}

#'  Convert linear effect to Emax-Hill effect
#'
#' @inheritParams convertEmax
#' @param g hill coefficient
#' @return Model with the linear effect converted to an Emax effect
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' readModelDb("PK_2cmt_no_depot") |>
#'   addIndirectLin(stim = "in") |>
#'   convertEmaxHill()
#'
#' # can also specify as emax=1
#'
#' readModelDb("PK_2cmt_no_depot") |>
#'   addIndirectLin(stim = "in") |>
#'   convertEmaxHill(emax = 1)
#'
convertEmaxHill <- function(
  ui,
  emax = "Emax",
  ec50 = "EC50",
  g = "g",
  imax = "Imax",
  ic50 = "IC50",
  ek = c("Ik", "Ek"),
  cc = c("Ec", "Cc")
) {
  .ui <- rxode2::assertRxUi(ui)
  cc <- rxode2::assertExists(.ui, cc)
  ek <- rxode2::assertVariableExists(.ui, ek)
  if (ek == "Ik") {
    emax <- imax
    ec50 <- ic50
  }
  .emaxMult <- paste0(emax, "*")
  if (inherits(emax, "character")) {
    rxode2::assertVariableNew(.ui, emax)
  } else if (is.numeric(emax) && emax == 1.0) {
    .emaxMult <- ""
  } else {
    stop("'", emax, "' not specified correctly", call. = FALSE)
  }
  rxode2::assertVariableNew(.ui, ec50)

  .modelLines <- .replaceMult(
    .ui$lstExpr,
    v1 = ek,
    v2 = cc,
    ret = paste0(.emaxMult, cc, "^", g, "/(", cc, "^", g, "+", ec50, "^", g, ")")
  )
  .tmp <- .getEtaTheta(.ui)
  .iniDf <- .tmp$iniDf
  .theta <- .tmp$theta
  .eta <- .tmp$eta
  .tmp <- .dropLines(.ui, .modelLines, .theta, .eta, ek)
  .modelLines <- .tmp$modelLines
  .theta <- .tmp$theta
  .eta <- .tmp$eta

  if (length(.theta$ntheta) == 0) {
    .ntheta <- 0
  } else {
    .ntheta <- max(.theta$ntheta)
  }

  .emaxLine <- NULL
  if (.emaxMult != "") {
    .emaxLine <- list(str2lang(paste0(emax, "<- exp(l", emax, ")")))
  }

  .ui <- rxode2::rxUiDecompress(.ui)
  # .theta/.eta are subsets of this model's own iniDf (ek dropped), so binding
  # them assumes nothing about its columns; the new parameters go in via ini().
  .ui$iniDf <- rbind(.theta, .eta)
  if (exists("description", envir = .ui$meta)) {
    rm("description", envir = .ui$meta)
  }
  rxode2::model(.ui) <- c(
    .emaxLine,
    list(str2lang(paste0(ec50, "<- exp(l", ec50, ")")), str2lang(paste0(g, "<- expit(lg", g, ", 0.1, 10)"))),
    .modelLines
  )
  if (!is.null(.emaxLine)) {
    .ui <- .iniAddTheta(.ui, paste0("l", emax), label = paste0("Maximum effect (", emax, ")"))
  }
  .ui <- .iniAddTheta(.ui, paste0("l", ec50), label = paste0("Concentration of 50% ", emax, " (", emax, ")"))
  .iniAddTheta(.ui, paste0("lg", g), est = logit(1, 0.1, 10), label = paste0("logit-constrained Hill coefficient ", g))
}
