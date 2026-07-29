indirect_0cpt_transitEx <- function() {
  description <- "Two compartment PK model with Michealis-Menten clearance using differential equations"
  reference <- "nlmixr2lib template"
  units <- list(time = "time_unit", dosing = "dose_unit", concentration = "conc_unit/vol_unit")
  ini({
    lka <- 0.45 ; label("Absorption rate (Ka)")
    lktr1 <- c(0, 0.05);  label("Transit rate constant 1 (log-scale)")
    lktr2 <- c(0, 0.05);  label("Transit rate constant 2 (log-scale)")
    lktr3 <- c(0, 0.05);  label("Transit rate constant 3 (log-scale)")
    lktr4 <- c(0, 0.05);  label("Transit rate constant 4 (log-scale)")
    lktr5 <- c(0, 0.05);  label("Transit rate constant 5 (log-scale)")
    lktr6 <- c(0, 0.05);  label("Transit rate constant 6 (log-scale)")
    lktr7 <- c(0, 0.05);  label("Transit rate constant 7 (log-scale)")
    lktr8 <- c(0, 0.05);  label("Transit rate constant 8 (log-scale)")
    lktr9 <- c(0, 0.05);  label("Transit rate constant 9 (log-scale)")
    lktr10 <- c(0, 0.05); label("Transit rate constant 10 (log-scale)")
    lvmax <- 0.04
    label("Maximum target-mediated rate of elimination Vmax (mg/L/d)")
    lkm <- 0.01
    label("Michaelis-Menten constant (mg/L)")
    lvc <- 3
    label("central volume of distribution (Vc)")
    lvp <- 5
    label("Peripheral volume of distribution (Vp)")
    lq <- 0.1
    label("Intercompartmental clearance (Q)")
    propSd <- c(0, 0.5)
    label("Proportional residual error (fraction)")

  })
  model({
    ka <- exp(lka)
    ktr1 <- exp(lktr1)
    ktr2 <- exp(lktr2)
    ktr3 <- exp(lktr3)
    ktr4 <- exp(lktr4)
    ktr5 <- exp(lktr5)
    ktr6 <- exp(lktr6)
    ktr7 <- exp(lktr7)
    ktr8 <- exp(lktr8)
    ktr9 <- exp(lktr9)
    ktr10 <- exp(lktr10)
    km   <- exp(lkm)
    vmax <- exp(lvmax)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q <- exp(lq)
    k12 <- q/vc
    k21 <- q/vp
    d/dt(depot) <- -ka*depot
    d/dt(transit1) <- ka * depot - ktr1 * transit1
    d/dt(transit2) <- ktr1 * transit1 - ktr2 * transit2
    d/dt(transit3) <- ktr2 * transit2 - ktr3 * transit3
    d/dt(transit4) <- ktr3 * transit3 - ktr4 * transit4
    d/dt(transit5) <- ktr4 * transit4 - ktr5 * transit5
    d/dt(transit6) <- ktr5 * transit5 - ktr6 * transit6
    d/dt(transit7) <- ktr6 * transit6 - ktr7 * transit7
    d/dt(transit8) <- ktr7 * transit7 - ktr8 * transit8
    d/dt(transit9) <- ktr8 * transit8 - ktr9 * transit9
    d/dt(transit10) <- ktr9 * transit9 - ktr10 * transit10
    d/dt(central) <- ktr10 * transit10 - (vmax/(km + central/vc)) *
      central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    Cc <- central/vc
    Cc ~ prop(propSd)
  })
}