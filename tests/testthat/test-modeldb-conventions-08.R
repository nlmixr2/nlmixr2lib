# Shard 8 of 8 of the whole-database convention sweep.
#
# Every model in the library is convention-checked, but split across files so
# testthat's workers can share the work: it parallelises across files, never
# within a test_that() block. See helper-conventionShard.R for how the shards
# are formed and why they interleave rather than follow modeldb's directories.

test_that("shard 8 of the model database satisfies the naming conventions", {
  res <- .conventionShardCheck(8L)
  expect_s3_class(res, "data.frame")
  expect_named(
    res,
    c("model", "category", "severity", "name", "message", "suggestion"),
    ignore.order = TRUE
  )
  expect_true(all(res$severity %in% c("error", "warning", "info")))
})
