# Mann 2022 translational opioid-overdose model: package-level helpers
# required by inst/modeldb/specificDrugs/Mann_2022_*.R and the Laffont
# 2025 vignette inline composed chain.
#
# Contents:
#   .onLoad                  - registers the lagState C user function at
#                              package load (idempotent)
#   register_lagState        - rxFun() registration wrapper
#   .lagState_c_code         - C source for the transport-delay function
#   Mann2022Equilibrate      - FDA pre-dose-equilibration helper

.onLoad <- function(libname, pkgname) {
  try(register_lagState(), silent = TRUE)
}

# ============================================================================
# lagState C user function
# ============================================================================

#' Register the lagState C user function with rxode2
#'
#' lagState() implements a true transport delay on an algebraic or state
#' expression evaluated during ODE integration. It is required by the
#' Mann_2022_respiratory_physiology model to encode the FDA delaymymod.c
#' peripheral (~7 s) and central (~11 s) chemoreflex transport delays
#' exactly, without resorting to first-order-lag approximations (which
#' damp amplitude as well as delaying, distorting the overdose
#' trajectory).
#'
#' The function signature is:
#'   `lagState(t, val, lag, channel, init_val)`
#'
#' Arguments inside an rxode2 `model()` block:
#'   * `t`        - current integration time (use the model symbol `t`)
#'   * `val`      - the algebraic / state expression to be delayed
#'   * `lag`      - delay duration, in the same time units as the model
#'   * `channel`  - integer channel index in `[0, 8)` so multiple delays
#'                  may coexist in one model without buffer collision
#'   * `init_val` - value returned when `t - lag < 0` (the pre-history
#'                  window) - typically the steady-state of `val`
#'
#' Behaviour:
#'   * Per-thread, per-channel ring buffer of recent (time, value) pairs.
#'   * Detects new-subject start (t returning to ~0) and resets the buffer.
#'   * Handles solver rollback (entries with t greater than the current
#'     integration time are discarded).
#'   * Linear interpolation between bracketing buffer entries.
#'   * If the buffer overflows (the requested lagged time has fallen off
#'     the oldest still-kept entry), `Rf_error` is raised with a clear
#'     diagnostic message - never silently return a wrong value.
#'
#' The C code uses GCC thread-local storage (`__thread`), so each
#' rxode2 worker thread maintains independent buffers. This is safe
#' because rxode2 parallelises across subjects (one subject per thread
#' at a time).
#'
#' This function is called automatically at package load via `.onLoad`,
#' so end users do not need to invoke it. Registration is idempotent.
#'
#' @keywords internal
#' @export
register_lagState <- function() {
  if (isTRUE(.lagState_registered$flag)) {
    return(invisible(NULL))
  }
  rxode2::rxFun(
    name = "lagState",
    args = c("t", "val", "lag", "channel", "init_val"),
    cCode = .lagState_c_code
  )
  .lagState_registered$flag <- TRUE
  invisible(NULL)
}

.lagState_registered <- new.env(parent = emptyenv())
.lagState_registered$flag <- FALSE

.lagState_c_code <- "
#include <R.h>
#include <Rinternals.h>
#include <string.h>

#define LAGSTATE_NUM_CHANNELS 8
#define LAGSTATE_BUFFER_SIZE  10000

static __thread double ls_buf_t[LAGSTATE_NUM_CHANNELS][LAGSTATE_BUFFER_SIZE];
static __thread double ls_buf_v[LAGSTATE_NUM_CHANNELS][LAGSTATE_BUFFER_SIZE];
static __thread int    ls_buf_n[LAGSTATE_NUM_CHANNELS] = {0};
static __thread double ls_last_t[LAGSTATE_NUM_CHANNELS] = {0};

double lagState(double t, double val, double lag, double channel, double init_val) {
  int ch = (int) channel;
  if (ch < 0 || ch >= LAGSTATE_NUM_CHANNELS) {
    Rf_error(\"lagState: channel %d out of range [0, %d)\", ch, LAGSTATE_NUM_CHANNELS);
  }

  /* New-subject detection: t returns to near zero after non-trivial progress */
  if (t < 1e-9 && ls_last_t[ch] > 1e-6) {
    ls_buf_n[ch] = 0;
  }

  /* Solver rollback: discard entries with time strictly greater than t */
  while (ls_buf_n[ch] > 0 && ls_buf_t[ch][ls_buf_n[ch]-1] > t + 1e-12) {
    ls_buf_n[ch]--;
  }

  /* Append or update */
  if (ls_buf_n[ch] > 0 && ls_buf_t[ch][ls_buf_n[ch]-1] >= t - 1e-12) {
    ls_buf_v[ch][ls_buf_n[ch]-1] = val;
  } else {
    if (ls_buf_n[ch] >= LAGSTATE_BUFFER_SIZE) {
      /* Drop oldest half; if interpolation later needs the dropped history we error */
      int keep = LAGSTATE_BUFFER_SIZE / 2;
      memmove(&ls_buf_t[ch][0], &ls_buf_t[ch][LAGSTATE_BUFFER_SIZE - keep], keep * sizeof(double));
      memmove(&ls_buf_v[ch][0], &ls_buf_v[ch][LAGSTATE_BUFFER_SIZE - keep], keep * sizeof(double));
      ls_buf_n[ch] = keep;
    }
    ls_buf_t[ch][ls_buf_n[ch]] = t;
    ls_buf_v[ch][ls_buf_n[ch]] = val;
    ls_buf_n[ch]++;
  }
  ls_last_t[ch] = t;

  /* Compute lagged value */
  double target_t = t - lag;
  if (target_t < 0.0) return init_val;
  if (ls_buf_n[ch] == 0) return init_val;
  double t_oldest = ls_buf_t[ch][0];
  if (target_t < t_oldest - 1e-9) {
    Rf_error(\"lagState: channel %d requested time %.6f, oldest buffer entry is %.6f - buffer overflow lost history\",
             ch, target_t, t_oldest);
  }
  if (target_t <= t_oldest) return ls_buf_v[ch][0];
  if (target_t >= ls_buf_t[ch][ls_buf_n[ch]-1]) return val;

  int lo = 0, hi = ls_buf_n[ch] - 1;
  while (hi - lo > 1) {
    int mid = (lo + hi) / 2;
    if (ls_buf_t[ch][mid] <= target_t) lo = mid;
    else hi = mid;
  }
  double t_lo = ls_buf_t[ch][lo];
  double t_hi = ls_buf_t[ch][hi];
  double v_lo = ls_buf_v[ch][lo];
  double v_hi = ls_buf_v[ch][hi];
  if (t_hi - t_lo < 1e-12) return v_lo;
  return v_lo + (v_hi - v_lo) * (target_t - t_lo) / (t_hi - t_lo);
}
"

# ============================================================================
# Mann2022Equilibrate: FDA pre-dose equilibration helper
# ============================================================================

#' Pre-equilibrate the Mann 2022 respiratory-physiology states
#'
#' The Mann 2022 / FDA delaymymod.c physiology layer is sensitive to
#' initial conditions because the chemoreflex feedback sees a 48.6 mm Hg
#' error in `p_b_co2 - p_b_co2_0` at t = 0 when the FDA delaystates.R
#' nominal values (`palv_co2 = 40.28`, `cb_co2 = 0.645`) are used
#' verbatim. The transient over-ventilation that follows leaves the
#' system in a state that subsequent opioid doses cannot push into
#' sustained PaO2 < 15 mm Hg, severely under-estimating cardiac arrest
#' incidence.
#'
#' The FDA `simulateToGetOD_IM.R` reference script handles this by
#' running a no-drug pre-simulation via `fundede` per subject and using
#' the FINAL state as the dose-time initial state. This function does
#' the equivalent: it runs the supplied `model` for `duration_min`
#' minutes with zero opioid dose and zero antagonist input, then returns
#' the final-time values of the standard Mann 2022 physiology states as
#' a named list suitable for passing as `inits` to `rxode2::rxSolve`.
#'
#' The equilibrium does NOT depend on opioid PK parameters (no drug is
#' present during pre-equilibration). The standalone
#' `Mann_2022_respiratory_physiology` model and the Laffont 2025
#' vignette inline chain ship with the FDA nominal initial conditions
#' as `ini()` defaults (matching `delaystates.R`), so any code that
#' uses these models for downstream overdose simulation MUST call this
#' function and pass the result via the `inits` argument of
#' `rxode2::rxSolve`; without it the FDA-published cardiac-arrest rates
#' cannot be reproduced.
#'
#' @param model An `rxode2` model containing at least the Mann 2022
#'   respiratory physiology states (`palv_co2`, `palv_o2`, `cb_co2`,
#'   `cb_o2`, `ct_co2`, `ct_o2`, `yco2`, `yo2`, `dp_state`, `dc_state`,
#'   `alpha_h`). Works with the standalone
#'   `Mann_2022_respiratory_physiology` model, with composed chains
#'   (PK + binding + physiology), or with any user-extended variant.
#' @param params Named vector of model parameters. For a composed
#'   chain, supply the full PK + binding + physiology parameter set
#'   (the opioid PK states stay at zero throughout the pre-equilibration
#'   because no dose is given).
#' @param duration_min Pre-equilibration duration in minutes. The
#'   default `90` is sufficient for the Mann 2022 physiology to settle
#'   to four-decimal precision; the slowest state is `alpha_h` with a
#'   5-minute time constant.
#' @return Named list of physiology equilibrium state values, suitable
#'   for passing as the `inits` argument to `rxode2::rxSolve` on
#'   subsequent dose-bearing simulations.
#' @importFrom stats setNames
#' @export
Mann2022Equilibrate <- function(model, params = NULL, duration_min = 90) {
  ev <- data.frame(
    time = seq(0, duration_min, length.out = 50L),
    evid = 0L,
    amt  = 0,
    cmt  = NA_character_,
    stringsAsFactors = FALSE
  )
  # Mann / Laffont composed chains expect L_ANT_pM and patient_type
  # covariates - default to zero antagonist, chronic patient (the
  # patient-type choice does not affect the pre-equilibration since
  # av1 = 1 with no opioid, but the column must exist).
  if (!"L_ANT_pM"     %in% colnames(ev)) ev$L_ANT_pM     <- 0
  if (!"patient_type" %in% colnames(ev)) ev$patient_type <- 1L
  if (!"CAR_OPIOID"   %in% colnames(ev)) ev$CAR_OPIOID   <- 0
  if (!"OPIOID_PATIENT_TYPE" %in% colnames(ev)) ev$OPIOID_PATIENT_TYPE <- 1L
  if (!"Q_TOTAL_LPM"  %in% colnames(ev)) ev$Q_TOTAL_LPM  <- 4.87

  sim <- as.data.frame(rxode2::rxSolve(model, params = params, events = ev))
  final <- sim[nrow(sim), , drop = FALSE]

  state_names <- c("palv_co2", "palv_o2", "cb_co2", "cb_o2",
                   "ct_co2", "ct_o2", "yco2", "yo2",
                   "dp_state", "dc_state", "alpha_h")
  available <- intersect(state_names, names(final))
  setNames(lapply(available, function(s) as.numeric(final[[s]])), available)
}
