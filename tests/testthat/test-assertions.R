test_that("assertCompartmentName", {
  expect_equal(assertCompartmentName("x"), "x")
  expect_equal(assertCompartmentName("x.y"), "x.y")
  expect_equal(assertCompartmentName("x_y"), "x_y")
  # This is a valid R variable name, but it probably doesn't work with nlmixr2
  # expect_equal(assertCompartmentName("."), ".")
  expect_error(assertCompartmentName("9"))
  expect_error(assertCompartmentName(9))
  expect_error(assertCompartmentName(NULL))
  expect_error(assertCompartmentName(c("A", "B")))
})

test_that("assertVariableName", {
  expect_equal(assertVariableName("x"), "x")
  expect_equal(assertVariableName("x.y"), "x.y")
  expect_equal(assertVariableName("x_y"), "x_y")
  # This is a valid R variable name, but it probably doesn't work with nlmixr2
  # expect_equal(assertVariableName("."), ".")
  expect_error(assertVariableName("9"))
  expect_error(assertVariableName(9))
  expect_error(assertVariableName(NULL))
  expect_error(assertVariableName(c("A", "B")))
})

test_that("assertParameterValue", {
  expect_equal(assertParameterValue(1), 1)
  expect_equal(assertParameterValue(-9), -9)
  expect_equal(assertParameterValue(0), 0)
  expect_error(assertParameterValue(Inf))
  expect_equal(assertParameterValue(NA))
})
