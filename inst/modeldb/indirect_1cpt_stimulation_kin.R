indirect_1cpt_stimulation_kin <- function() {
  description <- "One compartment indirect response model with stimulation of kin.Parameterized using rate cosntants"
  ini({
    lka  <- 0.45 ; label("Absorption rate (Ka)")
    lkel <- 0.534; label("elimination rate (1/d)")
    lEC50 <- 0.67; label("Drug concentration producing 50% of maximum stimulation at effect site (EC50)")
    lEmax <- 0.85; label("Maximum effect attributed to drug (Emax)")
    lkin <- 0.48; label("Zero-order rate constant for production of drug response(1/d)")
    lkout <- 0.34; label("First-order rate constant for loss of drug response")
    lfdepot <- 0.4; label("Bioavailability (F)")
    propSd <- 0.5 ; label("Proportional residual error (fraction)")
  })
  model({
    ka  <- exp(lka)
    kel <- exp(lkel)
    EC50 <- exp(lEC50)
    Emax <- exp(lEmax)
    kin <- exp(lkin)
    kout <- exp(lkout)
    fdepot   <- exp(lfdepot)
    
    d/dt(depot)      <- -ka*depot
    f(depot)         <- fdepot
    d/dt(central)    <- ka*depot -kel*central
    d/dt(effect) <- kin*(1+Emax*Cc/(Cc + IC50)) - kout*effect
    
    Cc <-  central/vc
    Cc ~ prop(propSd)
  })
}
