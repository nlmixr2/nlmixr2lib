# Shard 4 of 8 of the whole-database convention sweep.
#
# Every model in the library is convention-checked, but split across files so
# testthat's workers can share the work: it parallelises across files, never
# within a test_that() block. See helper-conventionShard.R for how the shards
# are formed and why they interleave rather than follow modeldb's directories.

test_that("shard 4 of the model database satisfies the naming conventions", {
  # Building every model in the shard is far too slow for CRAN, but this is
  # the enumerating check that a convention violation cannot reach a release,
  # so it must run in CI. skip_on_cran() does exactly that: it skips unless
  # NOT_CRAN=true, which the check workflow sets.
  skip_on_cran()
  res <- .conventionShardCheck(4L)
  expect_s3_class(res, "data.frame")
  expect_named(
    res,
    c("model", "category", "severity", "name", "message", "suggestion"),
    ignore.order = TRUE
  )
  expect_true(all(res$severity %in% c("error", "warning", "info")))
  # No model ships with an error-severity violation. buildModelDb() gates on
  # this, but its gate runs at build time only; asserting it here means a
  # model edited after the last rebuild cannot slip through.
  #
  # Compared as text rather than as a count: a count tells you a shard of 364
  # models has one violation somewhere, which is not enough to act on. The
  # diff on failure names the model, the rule and the parameter.
  .err <- res[res$severity == "error", , drop = FALSE]
  # paste0() recycles zero-length arguments against the literals, so on a
  # clean shard paste0(.err$model, " [", ...) is " [] : " -- length 1, not
  # length 0. Build the vector only when there is something to describe.
  .errText <- if (nrow(.err) == 0L) {
    character(0)
  } else {
    paste0(.err$model, " [", .err$category, "] ", .err$name, ": ", .err$message)
  }
  expect_identical(.errText, character(0))
})
