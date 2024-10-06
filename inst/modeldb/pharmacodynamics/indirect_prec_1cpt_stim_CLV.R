indirect_prec_1cpt_stim_CLV  <- function() {
  description <- "One compartment precursor-dependent indirect response model with inhibition of drug response (effect). Parameterized with clearance and volume"
  ini({
    lka  <- 0.45 ; label("Absorption rate (Ka)")
    lvc <- 3.45  ; label("Central volume of distribution (Vc)")
    lcl <- 0.534; label("clearance (CL)")
    lemax <- 0.56 ; label("Maximum inhibitory factor attributed to drug (Imax)")
    lEC50 <- 0.67; label("Drug concentration producing 50% of maximum inhibition at effect site (IC50)")
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
    EC50 <- exp(lEC50)
    emax<- exp(lemax)
    kin <- exp(lkin)
    kout <- exp(lkout)
    kpin <- exp(lkpin)
    kpout <- exp(lkpout)
    fdepot   <- exp(lfdepot)
    
    kel <- cl/vc
    
    d/dt(depot)      <- -ka*depot
    f(depot)         <- fdepot
    d/dt(central)    <- ka*depot-kel*central
    d/dt(precursor)  <- kpin -(kin + kpout)*(1+emax*Cc/(Cc + EC50))*precursor
    d/dt(effect)   <- kin*(1+emax*Cc/(Cc + EC50))*precursor-kout*effect
    
    Cc <-  central/vc
    Cc ~ prop(propSd)
  })
}
