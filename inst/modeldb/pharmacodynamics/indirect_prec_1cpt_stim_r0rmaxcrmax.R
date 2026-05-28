indirect_prec_1cpt_stim_r0rmaxcrmax <- function() {
  description <- "One compartment precursor-dependent indirect response model with inhibition of drug response (effect). Parameterized with clearance and volume"
  depends <- c("kin", "crmax")
  reference <- "nlmixr2lib template"
  units <- list(time = "time_unit", dosing = "dose_unit", concentration = "conc_unit/vol_unit")
  ini({
    lka  <- 0.45 ; label("Absorption rate (Ka)")
    lvc <- 3.45  ; label("Central volume of distribution (Vc)")
    lcl <- 0.534; label("clearance (CL)")
    lrbase     <- 0.2  ; label("Baseline response prior to drug administration (R0)")
    lrmax   <- 0.9  ; label("Maximal response (CRmax)")
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
    rbase<- exp(lrbase)
    rmax<- exp(lrmax)
    kout <- exp(lkout)
    kpin <- exp(lkpin)
    kpout <- exp(lkpout)
    fdepot   <- exp(lfdepot)
    
    kel <- cl/vc
    emax <- (rmax-rbase)/rmax
    kin <- (kout*(kin+kpout)*rbase)/kpin
    EC50 <- crmax*(rbase*(1+emax)-rmax)/(rmax-rbase)
    Cc <-  central/vc
    
    d/dt(depot)      <- -ka*depot
    f(depot)         <- fdepot
    d/dt(central)    <- ka*depot-kel*central
    d/dt(precursor1)  <- kpin -(kin + kpout)*(1+emax*Cc/(Cc + EC50))*precursor1
    d/dt(effect)   <- kin*(1+emax*Cc/(Cc + EC50))*precursor1-kout*effect
    
    Cc ~ prop(propSd)
  })
}
