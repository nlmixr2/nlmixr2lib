# Solves of templates whose model() once read a symbol it never defined,
# checked against the closed form each template implies, at the template's own
# ini() values. Each was unsolvable as shipped until the name was defined, so a
# solve that matches its closed form pins the fix.

.baselineTimes <- seq(0, 48, by = 0.25)

# Solve template `name` at its own ini() values, overridden by `params`, from
# `inits`. A template with more than one endpoint needs its observations placed
# on a named state, `obsCmt`.
.baselineSolve <- function(name, inits, params = numeric(), obsCmt = NULL) {
  ui <- rxode2::as.rxUi(nlmixr2lib::readModelDb(name))
  theta <- ui$theta
  theta[names(params)] <- params
  ev <- if (is.null(obsCmt)) rxode2::et(.baselineTimes) else rxode2::et(.baselineTimes, cmt = obsCmt)
  rxode2::rxSolve(
    ui,
    events = ev,
    params = theta,
    inits = inits,
    atol = 1e-12,
    rtol = 1e-10,
    returnType = "data.frame",
    addDosing = FALSE
  )
}

.baselineTheta <- function(name) {
  rxode2::as.rxUi(nlmixr2lib::readModelDb(name))$theta
}

test_that("the rate-constant indirect templates rest at kin/kout and settle at their closed-form drug effect", {
  skip_on_cran()
  # With the concentration held at ec50 the drug term is emax/2, so
  # stimulating production settles at R0*(1 + emax/2) and stimulating loss at
  # R0/(1 + emax/2), where R0 = kin/kout. Elimination is made negligible to
  # hold the concentration: the rate-constant templates estimate kel, the
  # _CLV one cl.
  stimulates <- c(
    indirect_1cpt_stim_kin = "kin",
    indirect_1cpt_stim_kin_CLV = "kin",
    indirect_1cpt_stim_kout = "kout"
  )
  for (nm in names(stimulates)) {
    theta <- .baselineTheta(nm)
    r0 <- exp(theta[["lkin"]]) / exp(theta[["lkout"]])
    free <- .baselineSolve(nm, c(effect = r0))
    expect_lt(max(abs(free$effect - r0)) / r0, 1e-7, label = paste(nm, "drug-free relative deviation"))

    elimination <- if ("lcl" %in% names(theta)) c(lcl = log(1e-12)) else c(lkel = log(1e-12))
    ec50 <- exp(theta[["lec50"]])
    h <- exp(theta[["lemax"]]) / 2
    expected <- if (stimulates[[nm]] == "kin") r0 * (1 + h) else r0 / (1 + h)
    held <- .baselineSolve(nm, c(central = ec50 * exp(theta[["lvc"]]), effect = r0), params = elimination)
    settled <- held$effect[held$time == max(.baselineTimes)]
    expect_lt(abs(settled - expected) / expected, 1e-7, label = paste(nm, "relative distance from the steady state"))
  }
})

test_that("the saturating logistic TGI templates follow their closed-form growth curves", {
  # dTS/dt = kge*TS*(1 - (TS/TSmax)^g) has the solution
  # TS(t) = TSmax / (1 + ((TSmax/TS0)^g - 1)*exp(-g*kge*t))^(1/g),
  # with g = 1 for the logistic and g = gamma for the generalized logistic.
  skip_on_cran()
  for (nm in c("tgi_sat_logistic", "tgi_sat_genLogistic")) {
    theta <- .baselineTheta(nm)
    ts0 <- exp(theta[["lrbase"]])
    tsmax <- exp(theta[["ltsmax"]])
    kge <- exp(theta[["lkge"]])
    g <- if ("lgamma" %in% names(theta)) exp(theta[["lgamma"]]) else 1
    sol <- .baselineSolve(nm, numeric(), obsCmt = "tumor_size")
    expected <- tsmax / (1 + ((tsmax / ts0)^g - 1) * exp(-g * kge * sol$time))^(1 / g)
    expect_lt(max(abs(sol$tumor_size - expected)) / tsmax, 1e-7, label = paste(nm, "relative deviation"))
    # The curve must actually move, or matching it would prove nothing.
    expect_gt(max(sol$tumor_size) - ts0, 0.5 * (tsmax - ts0), label = paste(nm, "growth"))
  }
})
