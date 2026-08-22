#' Build an event table for the Poels 2025 elranatamab QSP model
#'
#' `Poels_2025_elranatamab_qsp` reproduces the dose-to-dose attenuation of
#' cytokine release described by Poels et al. (2025) Supplementary Eq 37, which
#' rescales the cumulative cytokine exposure state `cauc` at the start
#' of every dosing interval:
#'
#' \deqn{C_{auc,N}(0) = \frac{C_{auc,N-1}(\tau)}{5\left(1 - \frac{(N+1)^2}{1.3^2 + (N+1)^2}\right)}}
#'
#' where \eqn{N} is the number of doses given. That is a discrete state reset,
#' which cannot be written inside an rxode2 ODE. This helper therefore emits it
#' as an `evid = 6` (multiply) record on `cauc` at each dose time,
#' carrying the reciprocal of the Eq 37 divisor as its `amt`, alongside the
#' subcutaneous dose records and the requested observation records.
#'
#' Doses are supplied in mg and converted to pmol using the elranatamab
#' molecular weight, because the model carries `depot` as a pmol amount
#' (Poels 2025 Supplementary Eq 1 divides the absorption flux by `Vc`).
#'
#' @param dose_time Numeric vector of dose times, in hours.
#' @param dose_mg Numeric vector of subcutaneous doses, in mg. Either length 1
#'   (recycled) or the same length as `dose_time`.
#' @param obs_time Numeric vector of observation times, in hours.
#' @param mw Elranatamab molecular weight in g/mol. Defaults to 148500, the
#'   "approximately 148.5 kDa" of the ELREXFIO US prescribing information
#'   (Description section).
#' @param id Subject identifier written to the `id` column.
#' @param cytokine_reset Logical; when `TRUE` (the default) the Eq 37 reset
#'   records are included. Set to `FALSE` to simulate the model without the
#'   dose-to-dose attenuation of cytokine release.
#' @param obs_cmt Compartment used for observation records. Must be an ODE state
#'   of the model; every algebraic observable (`Cc`, `mProtein`, `tumorBurden`,
#'   ...) is returned as a column at those rows regardless of which state is
#'   named.
#'
#' @return A `data.frame` with columns `id`, `time`, `amt`, `evid` and `cmt`,
#'   ordered by time, suitable for [rxode2::rxSolve()].
#'
#' @examples
#' # MagnetisMM-3 Cohort A: 12/32/76 mg on C1D1/C1D4/C1D8, then 76 mg weekly
#' dose_time <- c(0, 72, 168, seq(336, by = 168, length.out = 4))
#' dose_mg <- c(12, 32, rep(76, 5))
#' ev <- Poels_2025_elranatamab_qsp_events(
#'   dose_time = dose_time,
#'   dose_mg = dose_mg,
#'   obs_time = seq(0, 1344, by = 24)
#' )
#' head(ev)
#'
#' @references
#' Poels KE, Elmeliegy M, Hibma J, Wang D, Musante CJ, Shtylla B. Leveraging
#' quantitative systems pharmacology modeling for elranatamab regimen
#' optimization in relapsed or refractory multiple myeloma. npj Syst Biol Appl.
#' 2025;11:102. \doi{10.1038/s41540-025-00585-z}
#'
#' @export
Poels_2025_elranatamab_qsp_events <- function(dose_time,
                                              dose_mg,
                                              obs_time,
                                              mw = 148500,
                                              id = 1L,
                                              cytokine_reset = TRUE,
                                              obs_cmt = "central") {
  stopifnot(
    is.numeric(dose_time), length(dose_time) > 0L, !anyNA(dose_time),
    is.numeric(dose_mg), !anyNA(dose_mg), all(dose_mg > 0),
    is.numeric(obs_time), length(obs_time) > 0L, !anyNA(obs_time),
    is.numeric(mw), length(mw) == 1L, mw > 0,
    is.logical(cytokine_reset), length(cytokine_reset) == 1L,
    is.character(obs_cmt), length(obs_cmt) == 1L
  )
  if (length(dose_mg) == 1L) {
    dose_mg <- rep(dose_mg, length(dose_time))
  }
  if (length(dose_mg) != length(dose_time)) {
    stop("`dose_mg` must be length 1 or the same length as `dose_time`.",
         call. = FALSE)
  }
  if (is.unsorted(dose_time, strictly = TRUE)) {
    stop("`dose_time` must be strictly increasing.", call. = FALSE)
  }

  # mg -> pmol: (mg / 1000) g / (mw g/mol) * 1e12 pmol/mol
  dose_pmol <- dose_mg / mw * 1e9

  doses <- data.frame(
    id = id,
    time = dose_time,
    amt = dose_pmol,
    evid = 1L,
    cmt = "depot",
    .ord = 1L,
    stringsAsFactors = FALSE
  )

  obs <- data.frame(
    id = id,
    time = obs_time,
    amt = NA_real_,
    evid = 0L,
    cmt = obs_cmt,
    .ord = 3L,
    stringsAsFactors = FALSE
  )

  out <- if (cytokine_reset) {
    resets <- data.frame(
      id = id,
      time = dose_time,
      amt = Poels_2025_elranatamab_qsp_cauc_multiplier(seq_along(dose_time)),
      evid = 6L,
      cmt = "cauc",
      .ord = 2L,
      stringsAsFactors = FALSE
    )
    rbind(doses, resets, obs)
  } else {
    rbind(doses, obs)
  }

  out <- out[order(out$time, out$.ord), , drop = FALSE]
  out$.ord <- NULL
  rownames(out) <- NULL
  out
}

#' Eq 37 cumulative-cytokine-exposure multiplier
#'
#' Reciprocal of the Poels 2025 Supplementary Eq 37 divisor,
#' `5 * (1 - (n + 1)^2 / (1.3^2 + (n + 1)^2))`, evaluated at dose index `n`.
#' Exported so the reset applied by
#' [Poels_2025_elranatamab_qsp_events()] can be inspected and tested directly.
#'
#' @param n Integer vector of dose indices (1 for the first dose).
#'
#' @return Numeric vector of multipliers applied to the `cauc` state.
#'
#' @examples
#' Poels_2025_elranatamab_qsp_cauc_multiplier(1:5)
#'
#' @export
Poels_2025_elranatamab_qsp_cauc_multiplier <- function(n) {
  stopifnot(is.numeric(n), !anyNA(n), all(n >= 1))
  # Constants 5 and 1.3 are printed inside Supplementary Eq 37 itself; they
  # have no Supplementary Table 3 row.
  scale <- 5
  shape <- 1.3
  1 / (scale * (1 - (n + 1)^2 / (shape^2 + (n + 1)^2)))
}
