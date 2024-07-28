
# Define the tests
test_that("defaultCombine works correctly", {

  # By default the function works with camel case, not snake case

  # Test when no arguments are provided
  expect_error(defaultCombine(), "no arguments provided")

  # Test when a single character string of length 1 is provided
  expect_equal(defaultCombine("a"), "a")
  expect_equal(camelCombine("a"), "a")
  expect_equal(snakeCombine("a"), "a")
  expect_equal(dotCombine("a"), "a")
  expect_equal(blankCombine("a"), "a")


  # Test when a single character vector of length 2 or more is provided
  expect_equal(defaultCombine(c("a", "b")), "aB")
  expect_equal(camelCombine(c("a", "b")), "aB")
  expect_equal(snakeCombine(c("a", "b")), "a_b")
  expect_equal(dotCombine(c("a", "b")), "a.b")
  expect_equal(blankCombine(c("a", "b")), "ab")

  # Test when a single list is provided
  expect_equal(defaultCombine(list("a", "b")), "aB")
  expect_equal(camelCombine(list("a", "b")), "aB")
  expect_equal(dotCombine(list("a", "b")), "a.b")
  expect_equal(snakeCombine(list("a", "b")), "a_b")
  expect_equal(blankCombine(list("a", "b")), "ab")

  # Test when multiple arguments are provided
  expect_equal(defaultCombine("a", "b"), "aB")
  expect_equal(camelCombine("a", "b"), "aB")
  expect_equal(snakeCombine("a", "b"), "a_b")
  expect_equal(dotCombine("a", "b"), "a.b")
  expect_equal(blankCombine("a", "b"), "ab")

  # Test invalid arguments
  expect_error(defaultCombine(1), "invalid argument")
  expect_error(snakeCombine(1), "invalid argument")
  expect_error(dotCombine(1), "invalid argument")
  expect_error(blankCombine(1), "invalid argument")
  expect_error(camelCombine(1), "invalid argument")

  # now test changing the default combine changes the method
  setCombineType("snake")
  expect_equal(defaultCombine("a", "b"), "a_b")
  setCombineType("dot")
  expect_equal(defaultCombine("a", "b"), "a.b")
  setCombineType("blank")
  expect_equal(defaultCombine("a", "b"), "ab")
  # camel case should be called last to restore the default for this
  # package
  setCombineType("camel")
  expect_equal(defaultCombine("a", "b"), "aB")

})
