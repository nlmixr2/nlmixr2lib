indirect_prec_1cpt_inhi_r0rmaxcrmax  <- function() {
  description <- "One compartment precursor-dependent indirect response model with inhibition of drug response (effect)."
  ini({
    lka  <- 0.45 ; label("Absorption rate (Ka)")
    lvc <- 3.45  ; label("Central volume of distribution (Vc)")
    lcl   <- 0.04 ; label("Clearance (CL)")
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
    imax <- (r0-rmax)/r0
    kin <- (kout*(kin+kpout)*r0)/kpin
    IC50 <- crmax*(rmax-(1-imax)*r0)/(r0-rmax)
    
    d/dt(depot)      <- -ka*depot
    f(depot)         <- fdepot
    d/dt(central)    <- ka*depot-kel*central
    d/dt(precursor)  <- kpin -(kin + kpout)*(1-imax*Cc/(Cc + IC50))*precursor
    d/dt(effect)   <- kin*(1-imax*Cc/(Cc + IC50))*precursor-kout*effect
    
    Cc <-  central/vc
    Cc ~ prop(propSd)
  })
}
