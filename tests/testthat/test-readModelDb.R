test_that("readModelDb", {
  expect_error(
    readModelDb("foo"),
    regexp = "'name' not in database"
  )
  expect_true(is.function(readModelDb("PK_1cmt")))
  expect_message(
    readModelDb("oncology_xenograft_simeoni_2004"),
    "You can modify the number of damaged cell compartments"
  )
})
