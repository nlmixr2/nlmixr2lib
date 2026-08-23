test_that("Eq 37 multiplier matches the printed formula", {
  n <- 1:6
  expected <- 1 / (5 * (1 - (n + 1)^2 / (1.3^2 + (n + 1)^2)))
  expect_equal(Poels_2025_elranatamab_qsp_cauc_multiplier(n), expected)
  # Spot values computed by hand from Supplementary Eq 37:
  # divisor at N = 1 is 5*(1 - 4/(1.69+4)) = 8.45/5.69
  expect_equal(Poels_2025_elranatamab_qsp_cauc_multiplier(1), 5.69 / 8.45)
  expect_equal(Poels_2025_elranatamab_qsp_cauc_multiplier(2), 10.69 / 8.45)
  expect_equal(Poels_2025_elranatamab_qsp_cauc_multiplier(3), 17.69 / 8.45)
  # The divisor shrinks monotonically, so the multiplier grows monotonically:
  # cumulative cytokine exposure is progressively inflated, which drives the
  # inhibition term IH toward Imax and attenuates release dose over dose.
  expect_true(all(diff(Poels_2025_elranatamab_qsp_cauc_multiplier(1:20)) > 0))
  expect_error(Poels_2025_elranatamab_qsp_cauc_multiplier(0))
})

test_that("event table has the dose, reset and observation records", {
  ev <- Poels_2025_elranatamab_qsp_events(
    dose_time = c(0, 72, 168),
    dose_mg = c(12, 32, 76),
    obs_time = c(0, 24, 200)
  )
  expect_named(ev, c("id", "time", "amt", "evid", "cmt"))
  expect_equal(nrow(ev), 3L + 3L + 3L)
  expect_false(is.unsorted(ev$time))

  doses <- ev[ev$evid == 1L, ]
  expect_equal(doses$cmt, rep("depot", 3L))
  # 76 mg / 148500 g/mol = 5.11785e-7 mol = 5.11785e5 pmol
  expect_equal(doses$amt, c(12, 32, 76) / 148500 * 1e9)
  expect_equal(doses$amt[3], 511784.5, tolerance = 1e-5)

  resets <- ev[ev$evid == 6L, ]
  expect_equal(resets$cmt, rep("cauc", 3L))
  expect_equal(resets$amt, Poels_2025_elranatamab_qsp_cauc_multiplier(1:3))
  # the reset must follow its dose at the same time
  expect_true(all(which(ev$evid == 6L) > which(ev$evid == 1L)[seq_len(3L)]))

  obs <- ev[ev$evid == 0L, ]
  expect_equal(obs$cmt, rep("central", 3L))
  expect_true(all(is.na(obs$amt)))
})

test_that("cytokine_reset = FALSE drops only the reset records", {
  args <- list(dose_time = c(0, 168), dose_mg = 76, obs_time = c(0, 100))
  with_reset <- do.call(Poels_2025_elranatamab_qsp_events, args)
  without <- do.call(
    Poels_2025_elranatamab_qsp_events, c(args, list(cytokine_reset = FALSE))
  )
  expect_equal(nrow(with_reset) - nrow(without), 2L)
  expect_false(any(without$evid == 6L))
  expect_equal(
    without[, c("time", "evid", "cmt")],
    with_reset[with_reset$evid != 6L, c("time", "evid", "cmt")],
    ignore_attr = TRUE
  )
})

test_that("dose_mg is recycled and inputs are validated", {
  ev <- Poels_2025_elranatamab_qsp_events(
    dose_time = c(0, 168, 336), dose_mg = 76, obs_time = 0
  )
  expect_equal(ev$amt[ev$evid == 1L], rep(76 / 148500 * 1e9, 3L))

  expect_error(
    Poels_2025_elranatamab_qsp_events(c(0, 168), c(12, 32, 76), 0),
    "same length"
  )
  expect_error(
    Poels_2025_elranatamab_qsp_events(c(168, 0), 76, 0),
    "strictly increasing"
  )
  expect_error(Poels_2025_elranatamab_qsp_events(0, -1, 0))
})

test_that("the emitted evid = 6 records really rescale cauc", {
  # Guards the helper against a change in rxode2's multiply-event semantics,
  # using a minimal model so the test stays fast.
  skip_if_not_installed("rxode2")
  toy <- rxode2::rxode2({
    d / dt(depot) <- -0.001 * depot
    d / dt(cauc) <- 1
  })
  ev <- Poels_2025_elranatamab_qsp_events(
    dose_time = c(0, 10, 20), dose_mg = 76, obs_time = c(9.999, 19.999, 29.999),
    obs_cmt = "depot"
  )
  out <- as.data.frame(
    rxode2::rxSolve(toy, ev, params = c(dummy = 1), returnType = "data.frame")
  )
  mult <- Poels_2025_elranatamab_qsp_cauc_multiplier(1:3)
  # Closed form: reset at dose k multiplies, then 10 h of unit accumulation.
  expected <- Reduce(function(prev, k) prev * mult[k] + 10, 1:3, accumulate = TRUE, init = 0)[-1]
  expect_equal(out$cauc, expected, tolerance = 1e-3)
})
