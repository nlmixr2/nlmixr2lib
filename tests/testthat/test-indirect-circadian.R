# The circadian indirect-response templates make the drug-free response follow
# R(t) = rm + ra*cos(w*(t - tacro)), with w = 0.2616 rad/h (a 24-h period).
# The *_kin_t templates hold kout fixed and give production the rate
# kin_t = kout*R - R'; the *_kout_t templates hold kin fixed and give loss the
# rate kout_t = (kin - R')/R. Either way, a drug-free solve started at R(0)
# stays on R(t), and a dose moves the response in the direction the name says.

.circW <- 0.2616

# Values that keep both time-varying rates positive: rm > ra, kin > w*ra and
# kout*(rm - ra) > w*ra. The templates' own ini() values set rm = ra, where the
# baseline touches zero once a day.
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

.circBaseline <- function(t) {
  exp(.circParams[["lrm"]]) + exp(.circParams[["lra"]]) * cos(.circW * (t - exp(.circParams[["ltacro"]])))
}

# Solve one template from R(0), with or without a 100-unit dose into depot.
# Most of the templates take vc from outside (their `depends`), so it is
# supplied when the model has no lvc of its own.
.circSolve <- function(name, dose) {
  ui <- rxode2::as.rxUi(nlmixr2lib::readModelDb(name))
  pars <- .circParams[intersect(names(.circParams), ui$iniDf$name)]
  if ("vc" %in% ui$allCovs) {
    pars <- c(pars, vc = 10)
  }
  times <- seq(0, 72, by = 0.25)
  ev <- if (dose) rxode2::et(rxode2::et(amt = 100, cmt = "depot"), times) else rxode2::et(times)
  rxode2::rxSolve(
    ui,
    events = ev,
    params = pars,
    inits = c(effect = .circBaseline(0)),
    atol = 1e-10,
    rtol = 1e-10,
    returnType = "data.frame",
    addDosing = FALSE
  )
}

test_that("each circadian template's drug-free response stays on its circadian baseline", {
  skip_on_cran()
  for (nm in .circTemplates) {
    free <- .circSolve(nm, dose = FALSE)
    deviation <- max(abs(free$effect - .circBaseline(free$time))) / exp(.circParams[["lrm"]])
    expect_lt(deviation, 1e-7, label = paste(nm, "relative deviation"))
  }
})

test_that("a dose moves each circadian template's response the way its name says", {
  skip_on_cran()
  for (nm in .circTemplates) {
    key <- sub("^indirect_circ_1cpt_(inhi|stim)_(kin|kout)_.*$", "\\1_\\2", nm)
    free <- .circSolve(nm, dose = FALSE)
    dosed <- .circSolve(nm, dose = TRUE)
    expect_equal(dosed$time, free$time)
    window <- dosed$time >= 1 & dosed$time <= 48
    drugEffect <- dosed$effect[window] - free$effect[window]
    expect_true(all(sign(drugEffect) == .circDirection[[key]]), label = paste(nm, "direction"))
    # At least a tenth of the mesor, so the sign is not solver noise.
    expect_gt(max(abs(drugEffect)), 1, label = paste(nm, "largest drug effect"))
  }
})
