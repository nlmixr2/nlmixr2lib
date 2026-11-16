#' Reparameterize a transit chain by its mean transit time
#'
#' Replaces the transit rate `ktr <- exp(lktr)` of a transit chain with
#' the mean-transit-time form Monolix uses: `Mtt <- exp(lMtt)` and
#' `ktr <- (n+1)/Mtt`, where `n` is the number of transit compartments.
#' The estimated parameter becomes the mean transit time, and the
#' initial estimate is converted so the model is unchanged at the
#' starting point.  The chain rate expression keeps `ktr` as its
#' name, so the rest of the model is untouched.
#'
#' The number of transit compartments is counted from the model, so
#' the chain must already be present (see [addTransit()] and
#' [addSecondAbsorption()]).
#'
#' @inheritParams addDepot
#' @param transit the transit compartment prefix (`"transit"` for the
#'   first path, `"transitB"` for a second path added by
#'   [addSecondAbsorption()])
#' @param ktr the transit rate variable name
#' @param mtt the mean transit time variable name
#' @return a model with the transit chain parameterized by `mtt`
#' @family absorption
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' readModelDb("PK_1cmt_des") |> addTransit(3) |> convertTransitMtt()
#'
convertTransitMtt <- function(ui, transit = "transit", ktr = "ktr", mtt = "Mtt") {
  .useModelAsUi()
  .ui <- rxode2::assertRxUi(ui)
  rxode2::assertCompartmentName(transit)
  assertVariableName(ktr)
  assertVariableName(mtt)
  # count the chain: transit1 .. transitN must be consecutive
  .n <- 0L
  while (rxode2::testCompartmentExists(.ui, paste0(transit, .n + 1L))) {
    .n <- .n + 1L
  }
  if (.n == 0L) {
    stop("no transit chain '", transit, "1' is present", call. = FALSE)
  }
  # the chain must be exactly transit1..N: a gap, or a longer chain
  # beyond the loop, would otherwise be converted with the wrong n
  .chain <- rxode2::rxModelVars(.ui)$state
  .chain <- .chain[grepl(paste0("^", transit, "[0-9]+$"), .chain)]
  if (length(.chain) != .n) {
    stop("transit chain '", transit, "' is not consecutive (found ", length(.chain),
      " compartments, expected ", .n, ")", call. = FALSE)
  }
  .lktr <- paste0("l", ktr)
  .lmtt <- paste0("l", mtt)
  rxode2::assertVariableNew(.ui, mtt)
  rxode2::assertVariableNew(.ui, .lmtt)
  if (!rxode2::testVariableExists(.ui, .lktr)) {
    stop("'", ktr, "' is not a log-scale transit rate in the model", call. = FALSE)
  }
  # the chain's entry rate must be this chain's rate: a transit line
  # driven by another rate would be silently re-parameterized
  .entryPat <- paste0("^d/dt\\(", transit, "1\\)")
  .entry <- Filter(function(l) grepl(.entryPat, deparse1(l)), .ui$lstExpr)
  if (length(.entry) != 1L || !grepl(paste0("\\b", ktr, "\\b"), deparse1(.entry[[1L]]))) {
    stop("'d/dt(", transit, "1)' does not use the rate '", ktr, "'", call. = FALSE)
  }
  .modelLines <- .ui$lstExpr
  .exprs <- list(str2lang(paste0(ktr, " <- exp(", .lktr, ")")))
  .w <- which(vapply(
    .modelLines,
    function(l) {
      rxode2::.matchesLangTemplate(l, .exprs[[1L]])
    },
    logical(1),
    USE.NAMES = FALSE
  ))
  if (length(.w) != 1L) {
    stop("expected one '", ktr, " <- exp(", .lktr, ")' line", call. = FALSE)
  }
  # the estimate moves to the mean transit time: ktr = (n+1)/Mtt
  .tmp <- .getEtaTheta(.ui)
  .est <- .tmp$theta$est[.tmp$theta$name == .lktr]
  .mttEst <- (.n + 1) / exp(.est)
  .modelLines[[.w]] <- str2lang(paste0(mtt, " <- exp(", .lmtt, ")"))
  .modelLines <- c(
    .modelLines[seq_len(.w)],
    list(str2lang(paste0(ktr, " <- (", .n, " + 1)/", mtt))),
    .modelLines[-seq_len(.w)]
  )
  .theta <- .dropTheta(.tmp$theta, .lktr)
  .eta <- .dropEta(.tmp$eta, .lktr)
  .ui <- rxode2::rxUiDecompress(.ui)
  .ui$iniDf <- rbind(.theta, .eta)
  if (exists("description", envir = .ui$meta)) {
    rm("description", envir = .ui$meta)
  }
  rxode2::model(.ui) <- .modelLines
  .ui <- .iniAddTheta(.ui, .lmtt, est = log(.mttEst), label = paste0("Mean transit time (", mtt, ")"))
  rxode2::rxUiCompress(rxode2::as.rxUi(.ui))
}
