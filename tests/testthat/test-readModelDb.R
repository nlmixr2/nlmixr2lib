test_that("readModelDb", {
  expect_error(
    readModelDb("foo"),
    regexp = "'name' not in database"
  )
  expect_true(is.function(readModelDb("PK_1cmt")))
})
