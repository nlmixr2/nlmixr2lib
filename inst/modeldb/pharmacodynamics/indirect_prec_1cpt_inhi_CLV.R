indirect_prec_1cpt_inhi_CLV  <- function() {
  description <- "One compartment precursor-dependent indirect response model with inhibition of drug response. Parameterized with clearance and volume. (effect)."
  reference <- "nlmixr2lib template"
  units <- list(time = "time_unit", dosing = "dose_unit", concentration = "conc_unit/vol_unit")
  ini({
    lka  <- 0.45 ; label("Absorption rate (Ka)")
    lvc <- 3.45  ; label("Central volume of distribution (Vc)")
    lcl   <- 0.04 ; label("Clearance (CL)")
    limax <- 0.56 ; label("Maximum inhibitory factor attributed to drug (imax)")
    lic50 <- 0.67; label("Drug concentration producing 50% of maximum inhibition at effect site (ic50)")
    lkin <- 0.48; label("First-order rate constant for production of drug response(1/d)")
    lkout <- 0.34; label("First-order rate constant for loss of drug response")
    lkpin <- 0.45 ; label("Zero order rate constant for production of precursor (kpin)")
    lkpout <- 0.45 ; label("First order rate constant for loss of precursor (kpout)")
    lfdepot <- 0.4; label("Bioavailability (F)")
    propSd <- 0.5 ; label("Proportional residual error (fraction)")
  })
  model({
    ka  <- exp(lka)
    vc  <- exp(lvc)
    cl <- exp(lcl)
    ic50 <- exp(lic50)
    imax<- exp(limax)
    kin <- exp(lkin)
    kout <- exp(lkout)
    kpin <- exp(lkpin)
    kpout <- exp(lkpout)
    fdepot   <- exp(lfdepot)
    
    kel <- cl/vc
    Cc <-  central/vc
    
    d/dt(depot)      <- -ka*depot
    f(depot)         <- fdepot
    d/dt(central)    <- ka*depot-kel*central
    d/dt(precursor1)  <- kpin -(kin + kpout)*(1-imax*Cc/(Cc + ic50))*precursor1
    d/dt(effect)   <- kin*(1-imax*Cc/(Cc + ic50))*precursor1-kout*effect
    
    Cc ~ prop(propSd)
  })
}       
