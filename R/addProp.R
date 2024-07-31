#' Add a property to a compartment
#'
#' @param ui rxode2 ui object
#' @param prop property to add to a compartment:
#'
#'  - \code{F}: bioavailability
#'
#'  - \code{alag}: absorption lag time
#'
#' - \code{dur}: modeled duration of infusion
#'
#' - \code{rate}: modeled infusion rate
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
addCmtProp <- function(ui, prop=c("f", "alag", "dur", "rate", "ini"),
                       cmt) {
  .ui <- rxode2::assertRxUi(ui)
  .cmt <- as.character(substitute(cmt))
  cmt <-try(force(cmt), silent=TRUE)
  if (inherits(cmt, "try-error")) {
    cmt <- .cmt
  }
  if (rxode2::testCompartmentExists(.ui, cmt)) {
    .cmt <- rxode2::assertCompartmentExists(.ui, cmt)
    .modelLines <- .ui$lstExpr
    .var <- defaultCombine(prop, .cmt)
    .ui <- addLogEstimates(.ui, .var)
    .modelLines <- .ui$lstExpr
    .w <- .whichDdt(.modelLines, .cmt)
    .tmp <- .extractModelLinesAtW(.modelLines, .w)

    .ui <- rxode2::rxUiDecompress(.ui)
    if (exists("description", envir=.ui$meta)) {
      rm("description", envir=.ui$meta)
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
  cmt <-try(force(cmt), silent=TRUE)
  if (inherits(cmt, "try-error")) {
    cmt <- .cmt
  }
  addCmtProp(ui, prop="f", cmt=cmt)
}

#' @describeIn addCmtProp Adds the lag-time to a compartment in the model
#' @export
addLag <- function(ui, cmt) {
  .cmt <- as.character(substitute(cmt))
  cmt <-try(force(cmt), silent=TRUE)
  if (inherits(cmt, "try-error")) {
    cmt <- .cmt
  }
  addCmtProp(ui, prop="alag", cmt=cmt)
}

#' @describeIn addCmtProp Adds the modeled duration to a compartment in the model
#' @export
addDur <- function(ui, cmt) {
  .cmt <- as.character(substitute(cmt))
  cmt <-try(force(cmt), silent=TRUE)
  if (inherits(cmt, "try-error")) {
    cmt <- .cmt
  }
  addCmtProp(ui, prop="dur", cmt=cmt)
}
#' @describeIn addCmtProp Adds the modeled rate to a compartment in the model
#' @export
addRate <- function(ui, cmt) {
  .cmt <- as.character(substitute(cmt))
  cmt <-try(force(cmt), silent=TRUE)
  if (inherits(cmt, "try-error")) {
    cmt <- .cmt
  }
  addCmtProp(ui, prop="rate", cmt=cmt)
}
