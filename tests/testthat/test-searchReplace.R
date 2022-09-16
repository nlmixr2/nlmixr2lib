test_that("searchReplace", {
  model <- readModelDb("PK_1cmt")
  expect_equal(
    functionBody(
      searchReplace(model, find = "lka", replace = "lka + etalka")
    )[[4]][[2]][[2]],
    str2lang("ka <- exp(lka + etalka)")
  )
})
