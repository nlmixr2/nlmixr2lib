# The precursor-pool templates indirect_prec_1cpt_{inhi,stim}_r0 are
#   dP/dt = kpin - (kin + kpout)*D*P,   dR/dt = kin*D*P - kout*R,
# with D = 1 - Imax*C/(C + IC50) or 1 + Emax*C/(C + EC50), and they derive
# kpin from the drug-free baseline R0. Holding the concentration constant
# gives closed forms for the baseline, the steady state and the initial
# response, checked here at each template's own ini() values.

# -1: the drug lowers the response at first; 1: it raises it.
.precTemplates <- c(indirect_prec_1cpt_inhi_r0 = -1, indirect_prec_1cpt_stim_r0 = 1)

# The drug-free quantities implied by the template's log-scale parameters.
.precBaseline <- function(theta) {
  kin <- exp(theta[["lkin"]])
  kout <- exp(theta[["lkout"]])
  kpout <- exp(theta[["lkpout"]])
  r0 <- exp(theta[["lrbase"]])
  kpin <- r0 * kout * (kin + kpout) / kin
  list(r0 = r0, p0 = kpin / (kin + kpout), kin = kin, kpout = kpout, kpin = kpin)
}

# Half the drug's maximum effect, signed: the drug term D is 1 + .precHalf()
# at a concentration equal to the template's IC50 or EC50.
.precHalf <- function(theta) {
  if ("limax" %in% names(theta)) -exp(theta[["limax"]]) / 2 else exp(theta[["lemax"]]) / 2
}

.precPotency <- function(theta) {
  exp(theta[[if ("lic50" %in% names(theta)) "lic50" else "lec50"]])
}

# Solve template `name` at its own ini() values from its drug-free state, with
# the plasma concentration held at `conc`: central starts at conc*vc and
# clearance is negligible.
.precSolveAt <- function(name, conc, times) {
  ui <- rxode2::as.rxUi(nlmixr2lib::readModelDb(name))
  theta <- ui$theta
  base <- .precBaseline(theta)
  theta[["lcl"]] <- log(1e-12)
  rxode2::rxSolve(
    ui,
    events = rxode2::et(times),
    params = theta,
    inits = c(central = conc * exp(theta[["lvc"]]), precursor1 = base$p0, effect = base$r0),
    atol = 1e-12,
    rtol = 1e-10,
    returnType = "data.frame",
    addDosing = FALSE
  )
}

test_that("each precursor template rests at R0 and P0 without drug", {
  skip_on_cran()
  for (nm in names(.precTemplates)) {
    base <- .precBaseline(rxode2::as.rxUi(nlmixr2lib::readModelDb(nm))$theta)
    expect_gt(base$kpin, 0, label = paste(nm, "derived kpin"))
    sol <- .precSolveAt(nm, 0, seq(0, 20, by = 0.5))
    expect_lt(max(abs(sol$effect - base$r0)) / base$r0, 1e-9, label = paste(nm, "relative drift from R0"))
    expect_lt(max(abs(sol$precursor1 - base$p0)) / base$p0, 1e-9, label = paste(nm, "relative drift from P0"))
  }
})

test_that("at a constant concentration each precursor template returns to R0 as its pool adapts", {
  skip_on_cran()
  for (nm in names(.precTemplates)) {
    theta <- rxode2::as.rxUi(nlmixr2lib::readModelDb(nm))$theta
    base <- .precBaseline(theta)
    d <- 1 + .precHalf(theta)
    sol <- .precSolveAt(nm, .precPotency(theta), c(0, 60))
    final <- sol[sol$time == 60, ]
    expect_lt(abs(final$effect - base$r0) / base$r0, 1e-6, label = paste(nm, "relative distance from R0"))
    poolSteady <- base$kpin / ((base$kin + base$kpout) * d)
    expect_lt(
      abs(final$precursor1 - poolSteady) / poolSteady,
      1e-6,
      label = paste(nm, "relative distance from the adapted pool")
    )
  }
})

test_that("each precursor template first moves the way its name says, at its closed-form rate", {
  # From the drug-free state dR/dt(0) = kin*D*P0 - kout*R0 = kin*(D - 1)*P0.
  skip_on_cran()
  for (nm in names(.precTemplates)) {
    theta <- rxode2::as.rxUi(nlmixr2lib::readModelDb(nm))$theta
    base <- .precBaseline(theta)
    expected <- base$kin * .precHalf(theta) * base$p0
    expect_identical(sign(expected), unname(.precTemplates[[nm]]), label = paste(nm, "direction"))
    sol <- .precSolveAt(nm, .precPotency(theta), c(0, 1e-5))
    slope <- diff(sol$effect) / 1e-5
    expect_lt(abs(slope - expected) / abs(expected), 1e-3, label = paste(nm, "relative error of the initial slope"))
  }
})
