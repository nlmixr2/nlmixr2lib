# Run-time assertion that an `eta ~ variance | G` random effect was drawn per
# subject and level of G in a simulation.
#
# rxode2 decides from the simulation data whether G sits below or above the
# subject:
#
#   * some subject has two or more levels of G  -> one draw per subject and level
#     (what an occasion or a catheter needs);
#   * every subject has exactly one level       -> G is taken as a level ABOVE
#     the subject, and every subject on a level shares one draw, with no error
#     or warning;
#   * the whole table has a single level        -> rxSolve() fails with
#     "parameter(s) are required for solving: THETA[1]".
#
# The second case removes the between-occasion spread without any message, so an
# article that simulates a model with a below-subject `| G` eta calls this
# function on the simulation. checkModelConventions() requires the call (rule
# `nesting_assertion`).

#' Assert that a nested random effect carries one draw per subject and level
#'
#' Errors unless every (subject, level) pair of a simulation carries its own
#' draw of each named eta. With `n` given, errors unless each level carries
#' exactly `n` distinct draws.
#'
#' @param sim The simulation, a data frame (or `rxSolve()` result) with one row
#'   per record.
#' @param eta Character vector of the eta columns to check.
#' @param by The nesting column of the `eta ~ variance | by` declaration.
#' @param n `NULL` (the default) to require one distinct draw per (subject,
#'   level) pair; a positive whole number to require exactly that many distinct
#'   draws on every level instead.
#' @param id The subject column. The default is `sim.id` and `id`, whichever are
#'   present (both when they are).
#' @return `sim`, invisibly.
#' @noRd
assertDrawsPerLevel <- function(sim, eta, by, n = NULL, id = NULL) {
  checkmate::assertDataFrame(sim, min.rows = 1L)
  checkmate::assertCharacter(eta, any.missing = FALSE, min.len = 1L, min.chars = 1L, unique = TRUE)
  checkmate::assertString(by, min.chars = 1L)
  checkmate::assertCount(n, positive = TRUE, null.ok = TRUE)
  checkmate::assertCharacter(id, any.missing = FALSE, min.len = 1L, min.chars = 1L, null.ok = TRUE)
  df <- as.data.frame(sim)
  id <- .drawsSubjectColumns(df, id)
  absent <- setdiff(c(eta, by, id), names(df))
  if (length(absent) > 0L) {
    .drawsStop("the simulation has no column ", paste0("`", absent, "`", collapse = ", "))
  }
  level <- df[[by]]
  subject <- do.call(paste, c(unname(as.list(df[id])), sep = "\r"))
  if (anyNA(level) || anyNA(df[id])) {
    .drawsStop("the subject or the nesting column `", by, "` has NA values")
  }
  for (e in eta) {
    .drawsAssertOne(df[[e]], e, level, subject, by, n)
  }
  invisible(sim)
}

# The subject columns of a simulation: `id` when given, else `sim.id` and `id`
# as far as they exist.
.drawsSubjectColumns <- function(df, id) {
  if (!is.null(id)) {
    return(id)
  }
  cols <- intersect(c("sim.id", "id"), names(df))
  if (length(cols) == 0L) {
    .drawsStop("the simulation has no `id` or `sim.id` column; give `id`")
  }
  cols
}

# One eta: a single value per (subject, level) pair, and then distinct values
# across the pairs (`n` NULL) or exactly `n` distinct values on each level.
.drawsAssertOne <- function(x, eta, level, subject, by, n) {
  if (!is.numeric(x)) {
    .drawsStop("`", eta, "` is not numeric")
  }
  if (anyNA(x)) {
    .drawsStop("`", eta, "` has ", sum(is.na(x)), " NA value(s)")
  }
  # Compared exactly: the subjects of a level that share a draw carry the same
  # double, and a tiny-variance eta still differs between independent draws.
  draw <- as.numeric(x)
  pair <- paste(subject, level, sep = "\n")
  per_pair <- tapply(draw, pair, function(v) length(unique(v)))
  if (any(per_pair > 1L)) {
    .drawsStop(
      "`",
      eta,
      "` takes more than one value within ",
      sum(per_pair > 1L),
      " (subject, `",
      by,
      "`) pair(s)"
    )
  }
  value <- draw[!duplicated(pair)]
  # Only the levels that occur: a factor `by` may carry unused levels, which
  # tapply() would report as NA.
  lev <- factor(as.character(level[!duplicated(pair)]), levels = as.character(sort(unique(level))))
  distinct <- tapply(value, lev, function(v) length(unique(v)))
  subjects <- tapply(seq_along(lev), lev, length)
  if (is.null(n)) {
    if (anyDuplicated(value) > 0L) {
      .drawsStop(
        "`",
        eta,
        "` is shared: ",
        length(value),
        " (subject, `",
        by,
        "`) pairs carry ",
        length(unique(value)),
        " distinct draws; ",
        .drawsPerLevel(distinct, subjects)
      )
    }
    return(invisible(NULL))
  }
  if (any(distinct != n)) {
    .drawsStop(
      "`",
      eta,
      "` should have ",
      n,
      " distinct draws on each `",
      by,
      "` level: ",
      .drawsPerLevel(distinct, subjects)
    )
  }
  invisible(NULL)
}

.drawsPerLevel <- function(distinct, subjects) {
  paste0(
    "level ",
    names(distinct),
    ": ",
    distinct,
    " distinct of ",
    subjects,
    " subject(s)",
    collapse = "; "
  )
}

.drawsStop <- function(...) {
  stop(
    errorCondition(
      paste0("assertDrawsPerLevel(): ", ...),
      class = "nlmixr2lib_drawsPerLevel"
    )
  )
}
