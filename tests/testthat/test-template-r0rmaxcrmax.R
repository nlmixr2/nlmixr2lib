# The R0/Rmax/CRmax indirect-response templates are parameterized by the
# drug-free baseline R0, the extreme response Rmax and the concentration CRmax
# at which it is reached, the magnitude |S1| of the initial slope at high
# concentration, and Imax or Emax. Holding the concentration constant turns
# each template into a one-state system whose steady state has a closed form,
# so each defining relation can be checked directly, at the template's own
# ini() values.

.r0Templates <- c(
  indirect_1cpt_inhi_kin_r0rmaxcrmax = "inhi_kin",
  indirect_1cpt_inhi_kout_r0rmaxcrmax = "inhi_kout",
  indirect_1cpt_stim_kin_r0rmaxcrmax = "stim_kin",
  indirect_1cpt_stim_kout_r0rmaxcrmax = "stim_kout"
)

# Inhibiting production or stimulating loss lowers the response; inhibiting
# loss or stimulating production raises it.
.r0Direction <- c(inhi_kin = -1, inhi_kout = 1, stim_kin = 1, stim_kout = -1)

# The steady-state response as the concentration grows without bound, when
# the drug term reaches its maximum `maxEffect` (Imax or Emax).
.r0Asymptote <- function(type, r0, maxEffect) {
  switch(
    type,
    inhi_kin = r0 * (1 - maxEffect),
    inhi_kout = r0 / (1 - maxEffect),
    stim_kin = r0 * (1 + maxEffect),
    stim_kout = r0 / (1 + maxEffect)
  )
}

.r0Theta <- function(name) {
  rxode2::as.rxUi(nlmixr2lib::readModelDb(name))$theta
}

# Solve template `name` at its own ini() values from R0, with the plasma
# concentration held at `conc`: central starts at conc*vc and clearance is
# negligible.
.r0SolveAt <- function(name, conc, times) {
  theta <- .r0Theta(name)
  theta[["lcl"]] <- log(1e-12)
  rxode2::rxSolve(
    rxode2::as.rxUi(nlmixr2lib::readModelDb(name)),
    events = rxode2::et(times),
    params = theta,
    inits = c(central = conc * exp(theta[["lvc"]]), effect = exp(theta[["lrbase"]])),
    atol = 1e-12,
    rtol = 1e-10,
    returnType = "data.frame",
    addDosing = FALSE
  )
}

test_that("each R0/Rmax/CRmax template rests at R0 without drug", {
  skip_on_cran()
  for (nm in names(.r0Templates)) {
    r0 <- exp(.r0Theta(nm)[["lrbase"]])
    sol <- .r0SolveAt(nm, 0, seq(0, 20, by = 0.5))
    expect_lt(max(abs(sol$effect - r0)) / r0, 1e-9, label = paste(nm, "relative drift from R0"))
  }
})

test_that("at concentration CRmax each R0/Rmax/CRmax template settles at Rmax", {
  skip_on_cran()
  for (nm in names(.r0Templates)) {
    theta <- .r0Theta(nm)
    rmax <- exp(theta[["lrmax"]])
    sol <- .r0SolveAt(nm, exp(theta[["lcrmax"]]), c(0, 60))
    expect_lt(abs(sol$effect[sol$time == 60] - rmax) / rmax, 1e-6, label = paste(nm, "relative distance from Rmax"))
  }
})

test_that("at a high concentration each R0/Rmax/CRmax template moves past Rmax to its asymptote", {
  skip_on_cran()
  for (nm in names(.r0Templates)) {
    type <- .r0Templates[[nm]]
    theta <- .r0Theta(nm)
    r0 <- exp(theta[["lrbase"]])
    rmax <- exp(theta[["lrmax"]])
    maxEffect <- exp(theta[[if (startsWith(type, "inhi")) "limax" else "lemax"]])
    asymptote <- .r0Asymptote(type, r0, maxEffect)
    # Rmax lies between R0 and the asymptote, in the direction the name says.
    expect_identical(sign(rmax - r0), .r0Direction[[type]], label = paste(nm, "Rmax side of R0"))
    expect_identical(sign(asymptote - rmax), .r0Direction[[type]], label = paste(nm, "asymptote side of Rmax"))
    sol <- .r0SolveAt(nm, 1e6 * exp(theta[["lcrmax"]]), c(0, 60))
    expect_lt(
      abs(sol$effect[sol$time == 60] - asymptote) / asymptote,
      1e-4,
      label = paste(nm, "relative distance from the asymptote")
    )
  }
})

test_that("at a high concentration each R0/Rmax/CRmax template starts moving at |S1|", {
  skip_on_cran()
  for (nm in names(.r0Templates)) {
    type <- .r0Templates[[nm]]
    theta <- .r0Theta(nm)
    sol <- .r0SolveAt(nm, 1e6 * exp(theta[["lcrmax"]]), c(0, 1e-5))
    slope <- diff(sol$effect) / 1e-5
    expected <- .r0Direction[[type]] * exp(theta[["ls1"]])
    expect_lt(abs(slope - expected) / abs(expected), 1e-3, label = paste(nm, "relative error of the initial slope"))
  }
})
