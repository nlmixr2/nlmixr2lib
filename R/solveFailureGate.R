# Render gate for subjects rxode2 could not solve.
#
# When rxode2 cannot solve a subject (the ODE solver fails, or, from rxode2
# PR 1542 on, a state becomes non-finite), it sets every state of that subject
# to NA at every time (`badSolveExit()` in src/par_solve.cpp), raises the
# solve-wide warning "some ID(s) could not solve the ODEs correctly" and
# returns. The failure flag (`op->badSolve`) is internal and shared by every
# subject, so the result carries no per-subject record: a failed subject is
# recognised by all of its state columns being NA in all of its rows. An
# article that then filters NA rows or summarises with `na.rm = TRUE` renders
# as if nothing had failed.
#
# The article render scripts (.github/scripts/pkgdown-parallel-build.R and the
# consolidation merge gate) call .installSolveFailureGate() in the rendering
# process, which traces the exit of rxode2:::rxSolve.default and stops at the
# first solve with a failed subject. An article whose failures are expected,
# counted and asserted sets `options(nlmixr2lib.expectFailedSubjects = TRUE)`
# around that solve.

.solveFailureGateState <- new.env(parent = emptyenv())
.solveFailureGateState$label <- NA_character_
.solveFailureGateState$records <- list()

#' Find the subjects a simulation could not solve
#'
#' @param res The value of `rxSolve()`: a data frame, data.table, tibble,
#'   `rxSolve` object or matrix with column names.
#' @param states The model's state names (`rxode2::rxState()`).
#' @param lhs The model's calculated variables (`rxode2::rxLhs()`), used only
#'   when no state column is in the output.
#' @return Labels of the failed subjects, e.g. `"id=3"` or `"sim.id=2,id=3"`;
#'   `character(0)` when every subject solved or when the output has no state
#'   or calculated column to judge by.
#' @noRd
.solveFailedSubjects <- function(res, states, lhs = character(0)) {
  if (is.matrix(res) && !is.null(colnames(res))) {
    res <- as.data.frame(res)
  }
  if (!is.data.frame(res) || nrow(res) == 0L) {
    return(character(0))
  }
  cols <- intersect(states, names(res))
  if (length(cols) == 0L) {
    cols <- intersect(lhs, names(res))
  }
  if (length(cols) == 0L) {
    return(character(0))
  }
  allNa <- Reduce(`&`, lapply(cols, function(col) is.na(res[[col]])))
  key <- .solveSubjectKey(res)
  failed <- tapply(allNa, factor(key, levels = unique(key)), all)
  names(failed)[failed]
}

# One label per output row naming its subject: `sim.id` and `id`, whichever are
# in the output (a single-subject solve has neither).
.solveSubjectKey <- function(res) {
  keyCols <- intersect(c("sim.id", "id"), names(res))
  if (length(keyCols) == 0L) {
    return(rep("subject 1", nrow(res)))
  }
  do.call(
    paste,
    c(lapply(keyCols, function(k) paste0(k, "=", res[[k]])), sep = ",")
  )
}

#' Decide whether a solve trips the render gate
#'
#' @inheritParams .solveFailedSubjects
#' @param expectFailed `TRUE` when the article declared, with
#'   `options(nlmixr2lib.expectFailedSubjects = TRUE)`, that failed subjects
#'   are expected and accounted for.
#' @return A list: `failed`, the failed subject labels; `nSubjects`, the number
#'   of subjects in the output; and `fire`, `TRUE` when the render must stop.
#' @noRd
.solveFailureGateDecision <- function(
  res,
  states,
  lhs = character(0),
  expectFailed = isTRUE(getOption("nlmixr2lib.expectFailedSubjects", FALSE))
) {
  failed <- .solveFailedSubjects(res, states = states, lhs = lhs)
  nSubjects <- 0L
  if (length(failed) > 0L) {
    nSubjects <- length(unique(.solveSubjectKey(as.data.frame(res))))
  }
  list(
    failed = failed,
    nSubjects = nSubjects,
    fire = length(failed) > 0L && !isTRUE(expectFailed)
  )
}

#' The message the gate stops with
#'
#' @param label The article being rendered.
#' @param decision The value of `.solveFailureGateDecision()`.
#' @return A single string.
#' @noRd
.solveFailureGateMessage <- function(label, decision) {
  shown <- utils::head(decision$failed, 10L)
  more <- length(decision$failed) - length(shown)
  sprintf(
    paste0(
      "%s: rxSolve() could not solve %d of %d subject(s) (%s%s); rxode2 ",
      "replaced their output with NA. Fix the model or the article. If the ",
      "failures are expected, count and assert them in the article and set ",
      "options(nlmixr2lib.expectFailedSubjects = TRUE) around that solve."
    ),
    label,
    length(decision$failed),
    decision$nSubjects,
    paste(shown, collapse = " "),
    if (more > 0L) sprintf(" and %d more", more) else ""
  )
}

#' Exit hook traced into rxode2:::rxSolve.default
#'
#' Records every solve with a failed subject and stops unless the article opted
#' in.
#'
#' @param res The solve's return value.
#' @param object The `object` argument of `rxSolve.default()`, the compiled
#'   model.
#' @return `NULL`, invisibly; stops when the gate fires.
#' @noRd
.solveFailureGateExit <- function(res, object) {
  states <- tryCatch(rxode2::rxState(object), error = function(e) character(0))
  lhs <- tryCatch(rxode2::rxLhs(object), error = function(e) character(0))
  decision <- .solveFailureGateDecision(res, states = states, lhs = lhs)
  if (length(decision$failed) == 0L) {
    return(invisible(NULL))
  }
  label <- .solveFailureGateState$label
  msg <- .solveFailureGateMessage(label, decision)
  .solveFailureGateState$records[[length(.solveFailureGateState$records) + 1L]] <-
    list(label = label, accepted = !decision$fire, message = msg)
  if (decision$fire) {
    stop(structure(
      class = c("nlmixr2libSolveFailure", "error", "condition"),
      list(message = msg, call = NULL)
    ))
  }
  invisible(NULL)
}

#' Install the render gate for failed subjects
#'
#' Traces the exit of `rxode2:::rxSolve.default` in the calling process so that
#' a solve with a subject rxode2 could not solve stops with a message naming
#' `label` and the failed subjects. Used by the article render scripts; not for
#' interactive use.
#'
#' @param label The name of the article being rendered, used in the message.
#' @return `NULL`, invisibly.
#' @noRd
.installSolveFailureGate <- function(label) {
  checkmate::assertString(label, min.chars = 1L)
  .solveFailureGateState$label <- label
  .solveFailureGateState$records <- list()
  # The hook is embedded as a function object, so the traced body needs no
  # nlmixr2lib::: lookup; `object` is rxSolve.default()'s first argument.
  exitCall <- as.call(list(
    .solveFailureGateExit,
    quote(returnValue()),
    quote(object)
  ))
  suppressMessages(trace(
    "rxSolve.default",
    exit = exitCall,
    where = asNamespace("rxode2"),
    print = FALSE
  ))
  invisible(NULL)
}

#' Remove the render gate
#'
#' @return `NULL`, invisibly.
#' @noRd
.removeSolveFailureGate <- function() {
  suppressMessages(untrace("rxSolve.default", where = asNamespace("rxode2")))
  .solveFailureGateState$label <- NA_character_
  .solveFailureGateState$records <- list()
  invisible(NULL)
}

#' Failed solves the gate saw that the article did not opt in to
#'
#' A render can finish although the gate stopped a solve, when the article
#' caught the error (`try()`, `tryCatch()`, a chunk with `error = TRUE`); the
#' render script checks this after the render.
#'
#' @return The messages of the solves that fired the gate, a character vector.
#' @noRd
.solveFailureGateFired <- function() {
  recs <- .solveFailureGateState$records
  fired <- recs[!vapply(recs, `[[`, logical(1), "accepted")]
  vapply(fired, `[[`, character(1), "message")
}
