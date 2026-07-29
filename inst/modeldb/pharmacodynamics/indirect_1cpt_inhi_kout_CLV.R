indirect_1cpt_inhi_kout_CLV <- function() {
  description <- "One compartment indirect response model with inhibition of kout."
  reference <- "nlmixr2lib template"
  units <- list(time = "time_unit", dosing = "dose_unit", concentration = "conc_unit/vol_unit")
  ini({
    lka  <- 0.45 ; label("Absorption rate (Ka)")
    lvc  <- 3.45 ; label("Central volume of distribution (Vc)")
    lcl <- 1 ; label("Clearance (Cl)")
    lic50 <- 0.67; label("Drug concentration producing 50% of maximum inhibition at effect site (ic50)")
    lkin <- 0.48; label("Zero-order rate constant for production of drug response(1/d)")
    lkout <- 0.34; label("First-order rate constant for loss of drug response")
    lfdepot <- 0.4; label("Bioavailability (F)")
    propSd <- 0.5 ; label("Proportional residual error (fraction)")
  })
  model({
    ka  <- exp(lka)
    vc  <- exp(lvc)
    cl  <- exp(lcl)
    ic50 <- exp(lic50)
    kin <- exp(lkin)
    kout <- exp(lkout)
    fdepot   <- exp(lfdepot)
    
    kel <- cl/vc
    
    d/dt(depot)      <- -ka*depot
    f(depot)         <- fdepot
    d/dt(central)    <- ka*depot -(kel)*central
    Cc <-  central/vc
    
    d/dt(effect) <- kin - kout*(1-Cc/(Cc + ic50))*effect
    
    Cc ~ prop(propSd)
  })
}
