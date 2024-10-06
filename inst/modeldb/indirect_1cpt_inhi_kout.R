indirect_1cpt_inhi_kout  <- function() {
  description <- "One compartment indirect response model with inhibition of kout."
  ini({
    lka  <- 0.45 ; label("Absorption rate (Ka)")
    lvc  <- log(90) ; label("Central volume of distribution (Vc)")
    lkel <- 0.534; label("elimination rate (1/d)")
    lIC50 <- 0.67; label("Drug concentration producing 50% of maximum inhibition at effect site (IC50)")
    lkin <- 0.48; label("Zero-order rate constant for production of drug response(1/d)")
    lkout <- 0.34; label("First-order rate constant for loss of drug response")
    lfdepot <- 0.4; label("Bioavailability (F)")
    propSd <- 0.5 ; label("Proportional residual error (fraction)")
  })
  model({
    ka  <- exp(lka)
    vc  <- exp(lvc)
    kel <- exp(lkel)
    IC50 <- exp(lIC50)
    kin <- exp(lkin)
    kout <- exp(lkout)
    fdepot   <- exp(lfdepot)

    
    d/dt(depot)      <- -ka*depot
    f(depot)         <- fdepot
    d/dt(central)    <- ka*depot -(cl/vc)*central
    
    Cc <-  central/vc
    
    d/dt(effect) <- kin - kout*(1-Cc/(Cc + IC50))*effect
    
    
    Cc ~ prop(propSd)
  })
}
