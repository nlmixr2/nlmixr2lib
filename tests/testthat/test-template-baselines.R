# Drug-free solves of templates whose model() once read a symbol it never
# defined, checked against the closed form each template implies, at the
# template's own ini() values. Each was unsolvable as shipped until the name
# was defined, so a solve that matches its closed form pins the fix.

.baselineTimes <- seq(0, 48, by = 0.25)

# Solve template `name` at its own ini() values, with `inits`, with or without
# a 100-unit dose into depot. A template with more than one endpoint needs its
# observations placed on a named state, `obsCmt`.
.baselineSolve <- function(name, inits, dose = FALSE, obsCmt = NULL) {
  ev <- if (dose) rxode2::et(amt = 100, cmt = "depot") else rxode2::et()
  ev <- if (is.null(obsCmt)) rxode2::et(ev, .baselineTimes) else rxode2::et(ev, .baselineTimes, cmt = obsCmt)
  rxode2::rxSolve(
    rxode2::as.rxUi(nlmixr2lib::readModelDb(name)),
    events = ev,
    inits = inits,
    atol = 1e-10,
    rtol = 1e-10,
    returnType = "data.frame",
    addDosing = FALSE
  )
}

.baselineTheta <- function(name) {
  rxode2::as.rxUi(nlmixr2lib::readModelDb(name))$theta
}

test_that("the rate-constant indirect templates hold kin/kout drug-free and respond to a dose", {
  skip_on_cran()
  # Stimulating production raises the response; stimulating loss lowers it.
  direction <- c(
    indirect_1cpt_stim_kin = 1,
    indirect_1cpt_stim_kin_CLV = 1,
    indirect_1cpt_stim_kout = -1
  )
  for (nm in names(direction)) {
    theta <- .baselineTheta(nm)
    r0 <- exp(theta[["lkin"]]) / exp(theta[["lkout"]])
    free <- .baselineSolve(nm, c(effect = r0))
    expect_lt(max(abs(free$effect - r0)) / r0, 1e-7, label = paste(nm, "drug-free relative deviation"))
    dosed <- .baselineSolve(nm, c(effect = r0), dose = TRUE)
    drugEffect <- dosed$effect[dosed$time >= 1] - r0
    expect_true(all(sign(drugEffect) == direction[[nm]]), label = paste(nm, "direction"))
    expect_gt(max(abs(drugEffect)), 0.1 * r0, label = paste(nm, "largest drug effect"))
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
