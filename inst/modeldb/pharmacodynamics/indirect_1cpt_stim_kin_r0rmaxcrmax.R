indirect_1cpt_stim_kin_r0rmaxcrmax <- function() {
  description <- "One compartment indirect response model with stimulation of kin."
  depends <- c("Emax")
  reference <- "nlmixr2lib template"
  units <- list(time = "time_unit", dosing = "dose_unit", concentration = "conc_unit/vol_unit")
  ini({
    lka     <- 0.45 ; label("Absorption rate (Ka)")
    lvc     <- 3.45 ; label("Central volume of distribution (Vc)")
    lcl     <- 0.85 ; label("Clearance (Cl)")
    lrbase     <- 0.2  ; label("Baseline response prior to drug administration (R0)")
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
    rbase     <- exp(lrbase)
    rmax   <- exp(lrmax)
    s1     <- exp(ls1)
    emax   <- exp(lemax)
    crmax  <- exp(lcrmax)
    fdepot <- exp(lfdepot)
    
    kel <- cl/vc
    emax <- (rmax-rbase)/rmax
    kin <- s1/emax
    kout <- kin/rbase
    IC50 <- crmax*(rbase*(1+emax)-rmax)/(rmax-rbase)
    
    d/dt(depot)      <- -ka*depot
    f(depot)         <- fdepot
    d/dt(central)    <- ka*depot -(kel)*central
    
    Cc <-  central/vc
    
    d/dt(effect) <- kin*(1+Emax*Cc/(Cc + IC50)) - kout*effect
    
    Cc ~ prop(propSd)
  })
}
