indirect_0cpt_transitEx <- function() {
  description <- "Two compartment PK model with Michealis-Menten clearance using differential equations"
  ini({
    lka <- 0.45 ; label("Absorption rate (Ka)")
    lktr1 <- c(0, 0.05)
    lktr2 <- c(0, 0.05)
    lktr3 <- c(0, 0.05)
    lktr4 <- c(0, 0.05)
    lktr5 <- c(0, 0.05)
    lktr6 <- c(0, 0.05)
    lktr7 <- c(0, 0.05)
    lktr8 <- c(0, 0.05)
    lktr9 <- c(0, 0.05)
    lktr10 <- c(0, 0.05)
    lvm <- 0.04
    label("maximum target-mediated rate of elimination (mg/L/d)")
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
    km <- exp(lkm)
    vm <- exp(lvm)
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
    d/dt(central) <- ktr10 * transit10 - (vm/(km + central/vc)) * 
      central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    Cc <- central/vc
    Cc ~ prop(propSd)
  })
}