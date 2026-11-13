# The circadian indirect-response templates make the drug-free response follow
# R(t) = rm + ra*cos(w*(t - tacro)), with w = 0.2616 rad/h (a 24-h period).
# The *_kin_t templates hold kout fixed and give production the rate
# kin_t = kout*R - R'; the *_kout_t templates hold kin fixed and give loss the
# rate kout_t = (kin - R')/R. Either way, a drug-free solve started at R(0)
# stays on R(t), and a dose moves the response in the direction the name says.

.circW <- 0.2616

# Values that keep both time-varying rates positive: rm > ra, kin > w*ra and
# kout*(rm - ra) > w*ra.
.circParams <- c(
  lka = log(1),
  lkel = log(0.2),
  lvc = log(10),
  lrm = log(10),
  lra = log(3),
  ltacro = log(8),
  lic50 = log(1),
  lec50 = log(1),
  limax = log(0.5),
  lemax = log(2),
  lkin = log(5),
  lkout = log(0.5),
  lfdepot = log(1),
  propSd = 0
)

.circTemplates <- paste0(
  "indirect_circ_1cpt_",
  c(
    "inhi_kin_kin_t",
    "inhi_kin_kout_t",
    "inhi_kout_kin_t",
    "inhi_kout_kout_t",
    "stim_kin_kin_t",
    "stim_kin_kout_t",
    "stim_kout_kin_t",
    "stim_kout_kout_t"
  )
)

# Inhibiting production or stimulating loss lowers the response; inhibiting
# loss or stimulating production raises it.
.circDirection <- c(inhi_kin = -1, inhi_kout = 1, stim_kin = 1, stim_kout = -1)

# R(t) for the log-scale parameters `p`.
.circBaseline <- function(t, p) {
  exp(p[["lrm"]]) + exp(p[["lra"]]) * cos(.circW * (t - exp(p[["ltacro"]])))
}

# The template's ini() values, with `override` replacing those it shares.
.circThetas <- function(name, override) {
  theta <- rxode2::as.rxUi(nlmixr2lib::readModelDb(name))$theta
  shared <- intersect(names(override), names(theta))
  theta[shared] <- override[shared]
  theta
}

# Solve one template at the log-scale parameters `p` from R(0), with or
# without a 100-unit dose into depot.
.circSolve <- function(name, p, dose) {
  times <- seq(0, 72, by = 0.25)
  ev <- if (dose) rxode2::et(rxode2::et(amt = 100, cmt = "depot"), times) else rxode2::et(times)
  rxode2::rxSolve(
    rxode2::as.rxUi(nlmixr2lib::readModelDb(name)),
    events = ev,
    params = p,
    inits = c(effect = .circBaseline(0, p)),
    atol = 1e-10,
    rtol = 1e-10,
    returnType = "data.frame",
    addDosing = FALSE
  )
}

test_that("each circadian template's drug-free response stays on its circadian baseline", {
  skip_on_cran()
  for (nm in .circTemplates) {
    p <- .circThetas(nm, .circParams)
    free <- .circSolve(nm, p, dose = FALSE)
    deviation <- max(abs(free$effect - .circBaseline(free$time, p))) / exp(p[["lrm"]])
    expect_lt(deviation, 1e-7, label = paste(nm, "relative deviation"))
  }
})

test_that("a dose moves each circadian template's response the way its name says", {
  skip_on_cran()
  for (nm in .circTemplates) {
    key <- sub("^indirect_circ_1cpt_(inhi|stim)_(kin|kout)_.*$", "\\1_\\2", nm)
    p <- .circThetas(nm, .circParams)
    free <- .circSolve(nm, p, dose = FALSE)
    dosed <- .circSolve(nm, p, dose = TRUE)
    expect_equal(dosed$time, free$time)
    window <- dosed$time >= 1 & dosed$time <= 48
    drugEffect <- dosed$effect[window] - free$effect[window]
    expect_true(all(sign(drugEffect) == .circDirection[[key]]), label = paste(nm, "direction"))
    # At least a tenth of the mesor, so the sign is not solver noise.
    expect_gt(max(abs(drugEffect)), 1, label = paste(nm, "largest drug effect"))
  }
})

test_that("each circadian template's own ini() values give a positive baseline it stays on", {
  # With ra >= rm the baseline touches zero once a day, where a *_kout_t
  # template's loss rate divides by zero.
  skip_on_cran()
  for (nm in .circTemplates) {
    p <- .circThetas(nm, numeric())
    expect_lt(exp(p[["lra"]]), exp(p[["lrm"]]), label = paste(nm, "amplitude"))
    free <- .circSolve(nm, p, dose = FALSE)
    deviation <- max(abs(free$effect - .circBaseline(free$time, p))) / exp(p[["lrm"]])
    expect_lt(deviation, 1e-7, label = paste(nm, "relative deviation at the template's values"))
  }
})
