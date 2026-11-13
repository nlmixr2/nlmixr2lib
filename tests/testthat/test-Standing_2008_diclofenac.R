# Standing_2008_diclofenac, late dose record.
#
# rxode2 cannot solve the two transit() inputs of this model when they are dosed
# at a late absolute time (a dose record from several hundred hours after time
# zero for the faster children; the same dose at time zero solves). Depending on
# the rxode2 build the solve returns NaN with no message, or errors (one subject)
# or warns and replaces the values with NA (a population). The model's note and
# the article's two-admission check (admissions 100 h apart) rely on that.
#
# The second test goes red when the late dose returns a finite concentration:
# rxode2 fixed it, so remove the "Simulation limit" sentence from the model's
# population notes, the late-dose bullet and the 100-h spacing comment in
# vignettes/articles/Standing_2008_diclofenac.Rmd, and this file.
#
# The reproducer is one subject with every eta zeroed and ka2 doubled, which
# moves the failure from about 1800 h to 1000 h; the model's own typical ka2
# fails from a dose time between 1000 h and 2000 h.

# The concentration 1 h after a dose into both depots at `dose_time`: NA_real_
# when the solve errors, plus the warnings other than the model parse's.
.standingSolve <- function(dose_time, lka2_shift = log(2)) {
  covs <- list(WT = 70, CHILD = 0L, OCC = 1L)
  doses <- do.call(
    rbind,
    lapply(c("depot1", "depot2"), function(cm) {
      data.frame(id = 1L, time = dose_time, evid = 1L, amt = 1e6, cmt = cm, covs)
    })
  )
  obs <- data.frame(id = 1L, time = dose_time + 1, evid = 0L, amt = 0, cmt = "Cc", covs)
  events <- rbind(doses, obs)
  events <- events[order(events$time, -events$evid), ]
  lka2 <- log(log(2) / 1.06) + lka2_shift
  seen <- new.env(parent = emptyenv())
  seen$warnings <- character()
  sim <- withCallingHandlers(
    tryCatch(
      {
        mod <- rxode2::rxode2(readModelDb("Standing_2008_diclofenac"))
        rxode2::rxSolve(mod, events, omega = NA, params = c(lka2 = lka2))
      },
      error = function(e) NULL
    ),
    warning = function(w) {
      seen$warnings <- c(seen$warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(
    cc = if (is.null(sim)) NA_real_ else as.data.frame(sim)$Cc,
    # Parsing the model warns about its occasion etas' mu-referencing; that is not a solve warning.
    warnings = seen$warnings[!grepl("non-mu referenced", seen$warnings, fixed = TRUE)]
  )
}

test_that("a Standing 2008 dose record at time zero or 500 h solves, identically", {
  skip_on_cran()
  at_0 <- .standingSolve(0)
  at_500 <- .standingSolve(500)
  expect_length(at_0$cc, 1L)
  expect_true(is.finite(at_0$cc))
  expect_equal(at_0$cc, 8894.934, tolerance = 1e-6)
  expect_equal(at_500$cc, at_0$cc, tolerance = 1e-8)
  expect_identical(at_0$warnings, character())
  expect_identical(at_500$warnings, character())
})

test_that("a Standing 2008 dose record at 1000 h still has no finite concentration", {
  skip_on_cran()
  late <- .standingSolve(1000)
  expect_length(late$cc, 1L)
  expect_false(
    is.finite(late$cc),
    info = paste(
      "rxode2 now solves a Standing 2008 dose at 1000 h; see the header of",
      "tests/testthat/test-Standing_2008_diclofenac.R for what to remove"
    )
  )
})
