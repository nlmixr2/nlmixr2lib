indirect_prec_1cpt_stim_r0rmaxcrmax <- function() {
  description <- "One compartment precursor-dependent indirect response model with inhibition of drug response (effect). Parameterized with clearance and volume"
  ini({
    lka  <- 0.45 ; label("Absorption rate (Ka)")
    lvc <- 3.45  ; label("Central volume of distribution (Vc)")
    lcl <- 0.534; label("clearance (CL)")
    lr0     <- 0.2  ; label("Baseline response prior to drug administration (R0)")
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
    r0<- exp(lr0)
    rmax<- exp(lrmax)
    kout <- exp(lkout)
    kpin <- exp(lkpin)
    kpout <- exp(lkpout)
    fdepot   <- exp(lfdepot)
    
    kel <- cl/vc
    emax <- (rmax-r0)/rmax
    kin <- (kout*(kin+kpout)*r0)/kpin
    EC50 <- crmax*(r0*(1+emax)-rmax)/(rmax-r0)
    
    d/dt(depot)      <- -ka*depot
    f(depot)         <- fdepot
    d/dt(central)    <- ka*depot-kel*central
    d/dt(precursor)  <- kpin -(kin + kpout)*(1+emax*Cc/(Cc + EC50))*precursor
    d/dt(effect)   <- kin*(1+emax*Cc/(Cc + EC50))*precursor-kout*effect
    
    Cc <-  central/vc
    Cc ~ prop(propSd)
  })
}
