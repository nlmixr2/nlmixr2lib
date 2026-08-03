# buildModelDb() must FAIL on error-severity convention issues, not just
# print a count. It printed for a long time: a consolidation branch carried
# 579 error-severity violations through a green build because nothing
# consumed the number.

test_that("the convention-error accumulator starts empty and records errors", {
  nlmixr2lib:::.conventionErrorsReset()
  expect_silent(nlmixr2lib:::.conventionErrorsStopIfAny())

  issues <- data.frame(
    model = "Fake_2026_drug", category = "fixed_label_redundant",
    severity = "error", name = "e_wt_cl",
    message = "m", suggestion = "s", stringsAsFactors = FALSE
  )
  nlmixr2lib:::.conventionErrorsAdd("Fake_2026_drug", issues)
  expect_error(nlmixr2lib:::.conventionErrorsStopIfAny(),
               "1 convention error\\(s\\) in 1 model\\(s\\)")
  expect_error(nlmixr2lib:::.conventionErrorsStopIfAny(), "Fake_2026_drug")
  expect_error(nlmixr2lib:::.conventionErrorsStopIfAny(), "e_wt_cl")
  nlmixr2lib:::.conventionErrorsReset()
})

test_that("the accumulator aggregates across models rather than stopping at the first", {
  # Failing at the end means one run reports every offending model, instead
  # of surfacing one per rebuild.
  nlmixr2lib:::.conventionErrorsReset()
  mk <- function(nm, param) data.frame(
    model = nm, category = "deprecated_names", severity = "error", name = param,
    message = "m", suggestion = "s", stringsAsFactors = FALSE)
  nlmixr2lib:::.conventionErrorsAdd("A_2026_one", mk("A_2026_one", "allo_cl"))
  nlmixr2lib:::.conventionErrorsAdd("B_2026_two", mk("B_2026_two", "Km"))
  expect_error(nlmixr2lib:::.conventionErrorsStopIfAny(),
               "2 convention error\\(s\\) in 2 model\\(s\\)")
  expect_error(nlmixr2lib:::.conventionErrorsStopIfAny(), "A_2026_one")
  expect_error(nlmixr2lib:::.conventionErrorsStopIfAny(), "B_2026_two")
  nlmixr2lib:::.conventionErrorsReset()
})

test_that("reset clears a recorded error", {
  nlmixr2lib:::.conventionErrorsReset()
  nlmixr2lib:::.conventionErrorsAdd("X_2026_y", data.frame(
    model = "X_2026_y", category = "c", severity = "error", name = "n",
    message = "m", suggestion = "s", stringsAsFactors = FALSE))
  expect_error(nlmixr2lib:::.conventionErrorsStopIfAny())
  nlmixr2lib:::.conventionErrorsReset()
  expect_silent(nlmixr2lib:::.conventionErrorsStopIfAny())
})
