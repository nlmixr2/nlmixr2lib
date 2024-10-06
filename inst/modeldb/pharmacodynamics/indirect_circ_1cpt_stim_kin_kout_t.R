indirect_circ_1cpt_stim_kin_kout_t <- function() {
  description <- "One compartment indirect response model with stimulation of kin and circadian kout_t.Parameterized using rate constants"
  ini({
    lka  <- 0.45 ; label("Absorption rate (Ka)")
    lkel <- 0.534; label("elimination rate (1/d)")
    lrm  <- 0.62; label ("Mean Baseline for drug response (mesor)(Rm)")
    lra  <- 0.62; label ("Amplitude of drug response (Ra)")
    ltz  <- 0.62; label ("peak time (Acrophase) (Tz)")
    lEC50 <- 0.67; label("Drug concentration producing 50% of maximum stimulation at effect site (EC50)")
    lEmax <- 0.85; label("Maximum effect attributed to drug (Emax)")
    lkin <- 0.34; label("Zero-order rate constant for production of drug response")
    lfdepot <- 0.4; label("Bioavailability (F)")
    propSd <- 0.5 ; label("Proportional residual error (fraction)")
  })
  model({
    ka  <- exp(lka)
    kel <- exp(lkel)
    rm <- exp(lrm)
    ra <- exp(lra)
    tz <- exp(ltz)
    EC50 <- exp(lEC50)
    Emax <- exp(lEmax)
    kin <- exp(lkin)
    fdepot   <- exp(lfdepot)
    
    kout_t <- kin + 0.2616*ra*sin(0.2616*(t-tz))/(rm+ra*cos*(0.2616*(t-tz)))
    
    d/dt(depot)      <- -ka*depot
    f(depot)         <- fdepot
    d/dt(central)    <- ka*depot -kel*central
    d/dt(effect) <- kin*(1+Emax*Cc/(Cc + IC50)) - kout_t*effect
    
    Cc <-  central/vc
    Cc ~ prop(propSd)
  })
}
