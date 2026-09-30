#' Add a peripheral compartment to a model
#'
#' Adds a peripheral distribution compartment (`peripheral1`, then
#' `peripheral2`) to a central compartment, following the
#' `PK_2cmt_des`/`PK_3cmt_des` seed conventions: intercompartmental
#' clearance `q` (`q2` for the second) and peripheral volume `vp`
#' (`vp2`), with `k12 <- q/vc`, `k21 <- q/vp` (and `k13`/`k31` for the
#' second).  At most two peripheral compartments are supported,
#' matching the seeds.
#'
#' @inheritParams addDepot
#' @param n which peripheral to add: 1 (`peripheral1`, `q`, `vp`) or 2
#'   (`peripheral2`, `q2`, `vp2`); defaults to the first missing one
#' @return a model with the peripheral compartment added
#' @family distribution
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' readModelDb("PK_1cmt_des") |> addPeriph()
#'
#' readModelDb("PK_2cmt_des") |> addPeriph()
#'
addPeriph <- function(ui, n = NULL, central = "central", model) {
  .useModelAsUi()
  .ui <- rxode2::assertRxUi(ui)
  .cmt <- as.character(substitute(central))
  central <- try(force(central), silent = TRUE)
  if (inherits(central, "try-error")) {
    central <- .cmt
  }
  central <- rxode2::assertCompartmentExists(.ui, central)
  .have1 <- rxode2::testCompartmentExists(.ui, "peripheral1")
  .have2 <- rxode2::testCompartmentExists(.ui, "peripheral2")
  if (is.null(n)) {
    n <- if (!.have1) {
      1L
    } else if (!.have2) {
      2L
    } else {
      stop("both peripheral compartments are already present", call. = FALSE)
    }
  }
  checkmate::assertIntegerish(n, lower = 1L, upper = 2L, len = 1L, any.missing = FALSE)
  n <- as.integer(n)
  if (n == 1L && .have1) {
    stop("'peripheral1' is already present", call. = FALSE)
  }
  if (n == 2L && .have2) {
    stop("'peripheral2' is already present", call. = FALSE)
  }
  if (n == 2L && !.have1) {
    stop("'peripheral1' must be added before 'peripheral2'", call. = FALSE)
  }
  .sfx <- if (n == 1L) "" else "2"
  .periph <- paste0("peripheral", if (n == 1L) "1" else "2")
  .q <- paste0("q", .sfx)
  .vp <- paste0("vp", .sfx)
  .kIn <- if (n == 1L) "k12" else "k13"
  .kOut <- if (n == 1L) "k21" else "k31"
  # q/vp/kIn/kOut are this function's own variables: refuse when the
  # model already uses any of them, instead of silently overwriting
  # the user's lines and thetas
  for (.v in c(.q, .vp, .kIn, .kOut, paste0("l", .q), paste0("l", .vp))) {
    if (rxode2::testVariableExists(.ui, .v)) {
      stop("'", .v, "' is already present in the model", call. = FALSE)
    }
  }
  .modelLines <- .ui$lstExpr
  .w <- .whichDdt(.modelLines, central)
  .tmp <- .extractModelLinesAtW(.modelLines, .w)
  # extend the central ODE in place: - kIn*central + kOut*peripheralN
  .centralNew <- str2lang(paste0(
    deparse1(.tmp$w),
    " - ",
    .kIn,
    "*",
    central,
    " + ",
    .kOut,
    "*",
    .periph
  ))
  .periphLine <- str2lang(paste0(
    "d/dt(",
    .periph,
    ") <- ",
    .kIn,
    "*",
    central,
    " - ",
    .kOut,
    "*",
    .periph
  ))
  .rateLines <- list(
    str2lang(paste0(.q, " <- exp(l", .q, ")")),
    str2lang(paste0(.vp, " <- exp(l", .vp, ")")),
    str2lang(paste0(.kIn, " <- ", .q, "/vc")),
    str2lang(paste0(.kOut, " <- ", .q, "/", .vp))
  )
  # NOTE: the central volume is always `vc` (the seed convention);
  # models using another name (V, V1) fail loudly at solve time with
  # an undefined-variable error
  .modelLines <- c(.tmp$pre, .rateLines, list(.centralNew), list(.periphLine), .tmp$post)
  .ui <- rxode2::rxUiDecompress(.ui)
  if (exists("description", envir = .ui$meta)) {
    rm("description", envir = .ui$meta)
  }
  rxode2::model(.ui) <- .modelLines
  .ui <- .iniAddTheta(.ui, paste0("l", .q), est = 0.1, label = paste0("Intercompartmental clearance (", .q, ")"))
  .ui <- .iniAddTheta(.ui, paste0("l", .vp), est = 5, label = paste0("Peripheral volume of distribution (", .vp, ")"))
  rxode2::rxUiCompress(rxode2::as.rxUi(.ui))
}

#' Remove a peripheral compartment from a model
#'
#' Removes `peripheral2` (or `peripheral1` when it is the only one),
#' rewiring the central compartment and dropping the associated
#' parameters.  `peripheral1` cannot be removed while `peripheral2` is
#' still present.
#'
#' @inheritParams addDepot
#' @param n which peripheral to remove: 1 or 2; defaults to the
#'   highest-numbered one present
#' @return a model with the peripheral compartment removed
#' @family distribution
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' readModelDb("PK_2cmt_des") |> removePeriph()
#'
removePeriph <- function(ui, n = NULL, central = "central", model) {
  .useModelAsUi()
  .ui <- rxode2::assertRxUi(ui)
  .cmt <- as.character(substitute(central))
  central <- try(force(central), silent = TRUE)
  if (inherits(central, "try-error")) {
    central <- .cmt
  }
  central <- rxode2::assertCompartmentExists(.ui, central)
  .have1 <- rxode2::testCompartmentExists(.ui, "peripheral1")
  .have2 <- rxode2::testCompartmentExists(.ui, "peripheral2")
  if (is.null(n)) {
    n <- if (.have2) {
      2L
    } else if (.have1) {
      1L
    } else {
      stop("no peripheral compartment is present", call. = FALSE)
    }
  }
  checkmate::assertIntegerish(n, lower = 1L, upper = 2L, len = 1L, any.missing = FALSE)
  n <- as.integer(n)
  if (n == 1L && !.have1) {
    stop("'peripheral1' is not present", call. = FALSE)
  }
  if (n == 2L && !.have2) {
    stop("'peripheral2' is not present", call. = FALSE)
  }
  if (n == 1L && .have2) {
    stop("'peripheral1' cannot be removed while 'peripheral2' is present", call. = FALSE)
  }
  .sfx <- if (n == 1L) "" else "2"
  .periph <- paste0("peripheral", if (n == 1L) "1" else "2")
  .q <- paste0("q", .sfx)
  .vp <- paste0("vp", .sfx)
  .lqs <- paste0("l", .q)
  .lvp <- paste0("l", .vp)
  .modelLines <- .rmDdt(.ui$lstExpr, .periph)
  .modelLines <- .rmCmtPropLines(.modelLines, .periph)
  # strip the peripheral's terms from the central ODE; when nothing
  # but peripheral terms remain the ODE is gone too, and only an
  # explicit central input (e.g. an IV model) can follow — refuse
  # instead of emitting an empty d/dt() line
  .w <- .whichDdt(.modelLines, central)
  .tmp <- .extractModelLinesAtW(.modelLines, .w)
  .kIn <- if (n == 1L) "k12" else "k13"
  .kOut <- if (n == 1L) "k21" else "k31"
  .centralNew <- .dropDotAddExpr(.replaceMult(
    .replaceMult(.tmp$w, .kIn, central, "."),
    .kOut,
    .periph,
    "."
  ))[[1L]]
  # after the peripheral terms drop out, the right-hand side must
  # still mention the compartment or another state; a bare "." means
  # the ODE would be empty (the rxode2 error for `d/dt(central) <- .`
  # is cryptic, so refuse with a message instead)
  .rhs <- deparse1(.centralNew[[3L]])
  .mentions <- grepl(central, .rhs, fixed = TRUE) ||
    any(vapply(rxode2::rxModelVars(.ui)$state, function(s) grepl(s, .rhs, fixed = TRUE), logical(1), USE.NAMES = FALSE))
  if (!.mentions) {
    stop("removing '", .periph, "' would leave '", central, "' with no input", call. = FALSE)
  }
  .modelLines <- c(.tmp$pre, .centralNew, .tmp$post)
  .tmp2 <- .getEtaTheta(.ui)
  .theta <- .tmp2$theta
  .eta <- .tmp2$eta
  # drop the peripheral's own rate lines and q/vp assignments by
  # template match (removeLinesAndInis cascades each dropped variable
  # to its theta); the shared kel/vc lines survive because nothing
  # else references the removed names.  NOTE: .dropLines() cannot be
  # used here — its .removeLines() template leaves k12/k21 in place
  # and its cascade strips vc <-, which kel <- cl/vc still needs.
  .rm <- c(.kIn, .kOut, .q, .vp)
  .exprs <- unlist(
    lapply(.rm, function(v) {
      list(str2lang(paste0(v, "<- .")), str2lang(paste0(v, "= .")))
    }),
    recursive = FALSE
  )
  .w <- which(vapply(
    seq_along(.modelLines),
    function(i) {
      .cur <- .modelLines[[i]]
      any(vapply(.exprs, function(e) rxode2::.matchesLangTemplate(.cur, e), logical(1), USE.NAMES = FALSE))
    },
    logical(1),
    USE.NAMES = FALSE
  ))
  if (length(.w) > 0L) {
    .modelLines <- .modelLines[-.w]
  }
  .theta <- .dropTheta(.theta, c(paste0("l", .q), paste0("l", .vp)))
  .eta <- .dropEta(.eta, c(paste0("l", .q), paste0("l", .vp)))
  .ui <- rxode2::rxUiDecompress(.ui)
  .ui$iniDf <- rbind(.theta, .eta)
  if (exists("description", envir = .ui$meta)) {
    rm("description", envir = .ui$meta)
  }
  rxode2::model(.ui) <- .modelLines
  rxode2::rxUiCompress(rxode2::as.rxUi(.ui))
}
