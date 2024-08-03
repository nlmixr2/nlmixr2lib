test_that("whichDdt", {

  f <- function() {
    model({
      d/dt(central) <- kel * central - k12 * central + k21 *
        peripheral1
      d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    })
  }
  f <- rxode2::rxode2(f)


  expect_equal(.whichDdt(f$lstExpr, "central", ddt=TRUE), 1L)
  expect_error(.whichDdt(f$lstExpr, "matt", ddt=TRUE))
  expect_error(.whichDdt(f$lstExpr, "central", ddt=FALSE))


  f <- function() {
    model({
      matt <- 3
      d/dt(central) <- kel * central - k12 * central + k21 *
        peripheral1
      d/dt(central) <- d/dt(central) + 1
      d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    })
  }
  f <- rxode2::rxode2(f)

  expect_error(.whichDdt(f$lstExpr, "central", ddt=TRUE))
  expect_error(.whichDdt(f$lstExpr, "peripheral1", ddt=FALSE))
  expect_equal(.whichDdt(f$lstExpr, "matt", ddt=FALSE), 1L)

})
