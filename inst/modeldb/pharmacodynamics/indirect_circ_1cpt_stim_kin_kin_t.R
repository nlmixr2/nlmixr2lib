indirect_circ_1cpt_stim_kin_kin_t <- function() {
  description <- "One compartment indirect response model with stimulation of kin and circadian kin_t.Parameterized using rate cosntants"
  depends <- c("vc", "ic50")
  reference <- "nlmixr2lib template"
  units <- list(time = "time_unit", dosing = "dose_unit", concentration = "conc_unit/vol_unit")
  ini({
    lka  <- 0.45 ; label("Absorption rate (Ka)")
    lkel <- 0.534; label("elimination rate (1/d)")
    lrm  <- 0.62; label ("Mean Baseline for drug response (mesor)(Rm)")
    lra  <- 0.62; label ("Amplitude of drug response (Ra)")
    ltacro  <- 0.62; label ("peak time (Acrophase) (Tz)")
    lec50 <- 0.67; label("Drug concentration producing 50% of maximum stimulation at effect site (ec50)")
    lemax <- 0.85; label("Maximum effect attributed to drug (emax)")
    lkout <- 0.34; label("First-order rate constant for loss of drug response")
    lfdepot <- 0.4; label("Bioavailability (F)")
    propSd <- 0.5 ; label("Proportional residual error (fraction)")
  })
  model({
    ka  <- exp(lka)
    kel <- exp(lkel)
    rm <- exp(lrm)
    ra <- exp(lra)
    tacro <- exp(ltacro)
    ec50 <- exp(lec50)
    emax <- exp(lemax)
    kout <- exp(lkout)
    fdepot   <- exp(lfdepot)
    
    kin_t <- kout*rm+kout*ra*cos(0.2616*(t-tacro))-0.2616*ra*sin(0.2616*(t-tacro))
    Cc <-  central/vc
    
    d/dt(depot)      <- -ka*depot
    f(depot)         <- fdepot
    d/dt(central)    <- ka*depot -kel*central
    d/dt(effect) <- kin_t*(1+emax*Cc/(Cc + ic50)) - kout*effect
    
    Cc ~ prop(propSd)
  })
}
