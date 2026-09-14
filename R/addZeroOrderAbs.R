#' Convert a model to zero-order absorption (Monolix \code{Tk0})
#'
#' Zero-order absorption releases the dose into the central compartment
#' at a constant rate over a modeled duration \code{Tk0}, matching
#' Monolix's \code{absorption(type=2, Tk0)} oral route.  Following
#' monolix2rx's translation, this is expressed as a modeled duration on
#' the central compartment (\code{dur(central) <- tk0}), so a depot
#' compartment (and its \code{ka}) is removed first.  Dosing records
#' for the zero-order route must request the modeled duration in the
#' event table (\code{RATE = -2}; see \code{et(rate=-2)}).
#'
#' Lag time and bioavailability can be combined with [addLag()] and
#' [addLogitBioavailability()].
#'
#' @inheritParams addDepot
#' @inheritParams addTransit
#' @param tk0 zero-order absorption duration parameter name (Monolix's
#'   \code{Tk0})
#' @return a model with zero-order absorption
#' @family absorption
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' readModelDb("PK_1cmt_des") |> addZeroOrderAbs()
#'
#' # a model without a depot gets the same modeled duration input
#' readModelDb("PK_2cmt_no_depot") |> addZeroOrderAbs()
addZeroOrderAbs <- function(ui, central = "central", depot = "depot",
                            transit = "transit", ktr = "ktr", ka = "ka",
                            tk0 = "tk0") {
  .ui <- rxode2::assertRxUi(ui)
  central <- rxode2::assertCompartmentExists(.ui, central)
  assertCompartmentName(depot)
  assertVariableName(tk0)
  rxode2::assertVariableNew(.ui, tk0)
  .cp <- .ui$props$cmtProp
  if (!is.null(.cp) &&
        any(.cp$Compartment == central & .cp$Property == "dur")) {
    stop("modeled duration already present for compartment '", central, "'",
         call. = FALSE)
  }
  if (rxode2::testCompartmentExists(.ui, paste0(transit, "1"))) {
    .ui <- removeTransit(.ui,
                         central = central, depot = depot,
                         transit = transit, ktr = ktr, ka = ka)
    warning("transit compartments removed for zero-order absorption model",
            call. = FALSE)
  }
  if (rxode2::testCompartmentExists(.ui, depot)) {
    .ui <- removeDepot(.ui, central = central, depot = depot, ka = ka)
    warning("'", depot, "' removed for zero-order absorption model",
            call. = FALSE)
  }
  .ui <- addLogEstimates(.ui,
                         stats::setNames("Zero-order absorption duration (Tk0)",
                                         tk0))
  .modelLines <- .ui$lstExpr
  .w <- .whichDdt(.modelLines, central)
  .tmp <- .extractModelLinesAtW(.modelLines, .w)
  .ui <- rxode2::rxUiDecompress(.ui)
  if (exists("description", envir = .ui$meta)) {
    rm("description", envir = .ui$meta)
  }
  rxode2::model(.ui) <- c(.tmp$pre,
                          .tmp$w,
                          str2lang(paste0("dur(", central, ") <- ", tk0)),
                          .tmp$post)
  rxode2::rxUiCompress(.ui)
}

#' Remove zero-order absorption from a model
#'
#' Removes the modeled duration added by [addZeroOrderAbs()] and the
#' \code{tk0} parameter, leaving the dose as an intravenous bolus into
#' the central compartment.  First-order absorption can be restored
#' with [addDepot()].
#'
#' @inheritParams addDepot
#' @param tk0 zero-order absorption duration parameter name
#' @return a model where the zero-order absorption is removed
#' @family absorption
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' readModelDb("PK_1cmt_des") |> addZeroOrderAbs() |> removeZeroOrderAbs()
removeZeroOrderAbs <- function(ui, central = "central", tk0 = "tk0") {
  .ui <- rxode2::assertRxUi(ui)
  central <- rxode2::assertCompartmentExists(.ui, central)
  assertVariableName(tk0)
  .modelLines <- .ui$lstExpr
  .w <- .whichDdt(.modelLines, central, start = "dur(", end = ")")
  .modelLines <- .modelLines[-.w]
  .ui <- rxode2::rxUiDecompress(.ui)
  if (exists("description", envir = .ui$meta)) {
    rm("description", envir = .ui$meta)
  }
  rxode2::model(.ui) <- .modelLines
  .ui <- rxode2::rxUiCompress(.ui)
  if (rxode2::testVariableExists(.ui, tk0)) {
    .ui <- removeLinesAndInis(.ui, tk0)
  }
  rxode2::rxUiCompress(.ui)
}
