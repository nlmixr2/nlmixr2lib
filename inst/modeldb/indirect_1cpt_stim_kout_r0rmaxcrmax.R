indirect_1cpt_stim_kout_r0rmaxcrmax <- function() {
  description <- "One compartment indirect response model with stimulation of kout."
  ini({
    lka     <- 0.45 ; label("Absorption rate (Ka)")
    lvc     <- 3.45 ; label("Central volume of distribution (Vc)")
    lcl     <- 0.85 ; label("Clearance (Cl)")
    lr0     <- 0.2  ; label("Baseline response prior to drug administration (R0)")
    lrmax   <- 0.9  ; label("Maximal response (CRmax)")
    ls1     <- 1.0  ; label("Initial slope of the response versus time curve (S1)")
    lemax   <- 0.56 ; label("Maximum inhibitory factor attributed to drug (Imax)")
    lcrmax  <- 0.67 ; label("Plasma concentration of drug at the time of maximal response (CRmax)")
    lfdepot <- 0.4  ; label("Bioavailability (F)")
    propSd  <- 0.5  ; label("Proportional residual error (fraction)")
  })
  model({
    ka     <- exp(lka)
    vc     <- exp(lvc)
    cl     <- exp(lcl)
    r0     <- exp(lr0)
    rmax   <- exp(lrmax)
    s1     <- exp(ls1)
    emax   <- exp(lemax)
    crmax  <- exp(lcrmax)
    fdepot <- exp(lfdepot)
    
    kel <- cl/vc
    emax <- (r0-rmax)/rmax
    kin <- -s1/emax
    kout <- kin/r0
    IC50 <- crmax*(rmax*(1+emax)-r0)/(r0-rmax)
    
    d/dt(depot)      <- -ka*depot
    f(depot)         <- fdepot
    d/dt(central)    <- ka*depot -(kel)*central
    Cc <-  central/vc
    
  
    d/dt(effect) <- kin - kout*(1+Emax*Cc/(Cc + EC50))*effect
    
    
    Cc ~ prop(propSd)
  })
}
