ivsc_2cmt_mm <- function() {
  description <- "Two compartment TMDD model with Michaelis-Menten approximation.Absorption of the ligand is first-order (ka) "
  ini({
    lka  <- 0.45 ; label("Absorption rate (Ka)")
    lvc  <- 3.45 ; label("Central volume of distribution (Vc)")
    lvm  <- 0.04; label("maximum target-mediated rate of elimination (mg/L/d)")
    lkm  <- 0.01; label("Michaelis-Menten constant (mg/L)")
    lkel <- 0.534; label("elimination rate (1/d)")
    lk12 <- 0.48; label("central-to-peripheral rate (1/d)")
    lk21 <- 0.34; label("peripheral-to-central rate (1/d)")
    lfdepot <- 0.4; label("Bioavailability (F)")
    propSd <- 0.5 ; label("Proportional residual error (fraction)")
  })
  model({
    ka  <- exp(lka)
    vm  <- exp(lvm)
    km  <- exp(lkm)
    vc  <- exp(lvc)
    kel <- exp(lkel)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    fdepot   <- exp(lfdepot)
    
    d/dt(depot)      <- -ka*depot
    f(depot)         <- fdepot
    d/dt(central)    <- ka*depot -(vm/(km + central/vc))*central- k12*central + k21*peripheral1 - kel*central
    d/dt(peripheral1)<- k12*central - k21*peripheral1
    
    Cc <-  central/vc
    Cc ~ prop(propSd)
  })
}