#' Add a property to a compartment
#'
#' @param ui rxode2 ui object
#' @param prop property to add to a compartment:
#'
#'  - \code{F}: bioavailability
#'
#'  - \code{lag}: absorption lag time
#'
#' - \code{dur}: modeled duration of infusion
#'
#' - \code{rate}: modeled infusion rate
#'
#' - \code{ini}: initial value of the compartment
#'
#' @param cmt compartment to apply the property to
#'
#' @return rxode2 ui object with property applied
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' readModelDb("PK_3cmt_des") |> addCmtProp("f", "depot")
#'
#' readModelDb("PK_3cmt_des") |> addBioavailability(depot)
#'
#' readModelDb("PK_3cmt_des") |> addLag(depot)
#'
#' readModelDb("PK_3cmt_des") |> addDur(depot)
#'
#' readModelDb("PK_3cmt_des") |> addRate(depot)
#'
#' readModelDb("PK_3cmt_des") |> addIni(depot)
#'
addCmtProp <- function(ui, prop = c("f", "lag", "dur", "rate", "ini"),
                       cmt) {
  .ui <- rxode2::assertRxUi(ui)
  .cmt <- as.character(substitute(cmt))
  cmt <- try(force(cmt), silent = TRUE)
  if (inherits(cmt, "try-error")) {
    cmt <- .cmt
  }
  prop <- match.arg(prop)
  if (rxode2::testCompartmentExists(.ui, cmt)) {
    .cmt <- rxode2::assertCompartmentExists(.ui, cmt)
    .modelLines <- .ui$lstExpr
    .var <- defaultCombine(prop, .cmt)
    .ui <- addLogEstimates(.ui, .var)
    .modelLines <- .ui$lstExpr
    .w <- .whichDdt(.modelLines, .cmt)
    .tmp <- .extractModelLinesAtW(.modelLines, .w)

    .ui <- rxode2::rxUiDecompress(.ui)
    if (exists("description", envir = .ui$meta)) {
      rm("description", envir = .ui$meta)
    }
    if (prop == "ini") {
      rxode2::model(.ui) <- c(.tmp$pre,
        .tmp$w,
        str2lang(paste0(.cmt, "(0) <- ", .var)),
        .tmp$post)
    } else {
      rxode2::model(.ui) <- c(.tmp$pre,
        .tmp$w,
        str2lang(paste0(prop, "(", .cmt, ") <- ", .var)),
        .tmp$post)
    }
    .ui
  } else {
    stop("Compartment ", cmt, " does not exist")
  }
}

#' @describeIn addCmtProp Adds the bioavailability to a compartment in the model
#' @export
addBioavailability <- function(ui, cmt) {
  .cmt <- as.character(substitute(cmt))
  cmt <- try(force(cmt), silent = TRUE)
  if (inherits(cmt, "try-error")) {
    cmt <- .cmt
  }
  addCmtProp(ui, prop = "f", cmt = cmt)
}

#' Add a logit-parameterized bioavailability to a compartment
#'
#' Adds a bioavailability fraction constrained to (0,1) with a
#' logit/expit parameterization, matching Monolix-style \code{F}
#' (oral, subcutaneous) and the \code{F1}/\code{Fd} dose-apportionment
#' fractions.  When \code{cmt2} is supplied, a single fraction is
#' estimated for \code{cmt} and \code{f(cmt2)} is set to \code{1 -
#' f(cmt)} (the double-absorption \code{F1} dose split).
#'
#' Unlike [addBioavailability()], which log-parameterizes \code{F}
#' unboundedly, this keeps \code{F} in (0,1): the initial estimate is
#' stored on the logit scale and the model block defines `fDepot <-
#' expit(lgfDepot)`, which keeps the mu-referenced bounds.
#'
#' @param ui rxode2 ui object
#' @param cmt compartment to apply the bioavailability to
#' @param cmt2 optional second compartment; when supplied its
#'   bioavailability is set to `1 - f(cmt)` (dose split)
#' @param f initial bioavailability fraction, in (0,1)
#' @return rxode2 ui object with the logit bioavailability applied
#' @export
#' @family absorption
#' @author Matthew L. Fidler
#' @examples
#'
#' readModelDb("PK_1cmt_des") |> addLogitBioavailability(depot)
#'
#' # dose split between two absorption paths (F1, 1-F1)
#' readModelDb("PK_1cmt_des") |>
#'   addDepot(depot = "depot2", ka = "ka2") |>
#'   addLogitBioavailability(depot, depot2)
#'
addLogitBioavailability <- function(ui, cmt, cmt2 = NULL, f = 0.8) {
  checkmate::assertNumeric(f, len = 1L, any.missing = FALSE)
  if (f <= 0 || f >= 1) {
    stop("f must be in (0, 1)", call. = FALSE)
  }
  .cmt <- as.character(substitute(cmt))
  cmt <- try(force(cmt), silent = TRUE)
  if (inherits(cmt, "try-error")) {
    cmt <- .cmt
  }
  .cmt2 <- as.character(substitute(cmt2))
  cmt2 <- try(force(cmt2), silent = TRUE)
  if (inherits(cmt2, "try-error")) {
    cmt2 <- .cmt2
  }
  if (is.character(cmt2) && identical(cmt2, "NULL")) {
    cmt2 <- NULL
  }
  .ui <- rxode2::assertRxUi(ui)
  .cmt1 <- rxode2::assertCompartmentExists(.ui, cmt)
  .split <- !is.null(cmt2)
  if (.split) {
    .cmt2 <- rxode2::assertCompartmentExists(.ui, cmt2)
    if (identical(.cmt1, .cmt2)) {
      stop("cmt and cmt2 must be different compartments", call. = FALSE)
    }
  }
  .var <- defaultCombine("f", .cmt1)
  .lvar <- paste0("lg", .var)
  .label <- if (.split) {
    paste0("Fraction of dose absorbed from ", .cmt1, " (", .var, ")")
  } else {
    paste0("Bioavailability fraction (", .var, ")")
  }
  .tmp <- .getEtaThetaTheta1(.ui)
  .theta <- .tmp$theta
  .theta1 <- .tmp$theta1
  .eta <- .tmp$eta
  if (length(.theta$ntheta) == 0) {
    .ntheta <- 0
  } else {
    .ntheta <- max(.theta$ntheta)
  }
  .theta <- rbind(.theta,
    .get1theta(.var, .theta1, .ntheta,
      est = logit(f),
      label = .label,
      name = .lvar
    ))
  .modelLines <- .ui$lstExpr
  .w1 <- .whichDdt(.modelLines, .cmt1)
  .f1 <- str2lang(paste0("f(", .cmt1, ") <- ", .var))
  if (.split) {
    .w2 <- .whichDdt(.modelLines, .cmt2)
    .f2 <- str2lang(paste0("f(", .cmt2, ") <- 1-", .var))
  }
  .out <- list(str2lang(paste0(.var, " <- expit(", .lvar, ")")))
  for (.i in seq_along(.modelLines)) {
    .out[[length(.out) + 1L]] <- .modelLines[[.i]]
    if (.i == .w1) {
      .out[[length(.out) + 1L]] <- .f1
    }
    if (.split && .i == .w2) {
      .out[[length(.out) + 1L]] <- .f2
    }
  }
  .ui <- rxode2::rxUiDecompress(.ui)
  .ui$iniDf <- rbind(.theta, .eta)
  if (exists("description", envir = .ui$meta)) {
    rm("description", envir = .ui$meta)
  }
  rxode2::model(.ui) <- .out
  rxode2::rxUiCompress(.ui)
}

#' @describeIn addCmtProp Adds the lag-time to a compartment in the model
#' @export
addLag <- function(ui, cmt) {
  .cmt <- as.character(substitute(cmt))
  cmt <- try(force(cmt), silent = TRUE)
  if (inherits(cmt, "try-error")) {
    cmt <- .cmt
  }
  addCmtProp(ui, prop = "lag", cmt = cmt)
}

#' @describeIn addCmtProp Adds the modeled duration to a compartment in the model
#' @export
addDur <- function(ui, cmt) {
  .cmt <- as.character(substitute(cmt))
  cmt <- try(force(cmt), silent = TRUE)
  if (inherits(cmt, "try-error")) {
    cmt <- .cmt
  }
  addCmtProp(ui, prop = "dur", cmt = cmt)
}
#' @describeIn addCmtProp Adds the modeled rate to a compartment in the model
#' @export
addRate <- function(ui, cmt) {
  .cmt <- as.character(substitute(cmt))
  cmt <- try(force(cmt), silent = TRUE)
  if (inherits(cmt, "try-error")) {
    cmt <- .cmt
  }
  addCmtProp(ui, prop = "rate", cmt = cmt)
}

#' @describeIn addCmtProp Adds the initial value to the compartment
#' @export
addIni <- function(ui, cmt) {
  .cmt <- as.character(substitute(cmt))
  cmt <- try(force(cmt), silent = TRUE)
  if (inherits(cmt, "try-error")) {
    cmt <- .cmt
  }
  addCmtProp(ui, prop = "ini", cmt = cmt)
}
