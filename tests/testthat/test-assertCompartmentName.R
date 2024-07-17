test_that("assertCompartmentName", {
  expect_equal(assertCompartmentName("x"), "x")
  expect_equal(assertCompartmentName("x.y"), "x.y")
  # This is a valid R variable name, but it probably doesn't work with nlmixr2
  # expect_equal(assertCompartmentName("."), ".")
  expect_error(assertCompartmentName("9"))
  expect_error(assertCompartmentName(9))
  expect_error(assertCompartmentName(NULL))
  expect_error(assertCompartmentName(c("A", "B")))
})
