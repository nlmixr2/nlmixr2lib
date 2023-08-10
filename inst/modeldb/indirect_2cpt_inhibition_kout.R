indirect_2cpt_inhibition_kout <- function() {
  description <- "Two compartment indirect response model with inhibition of kout.Parameterized using rate constants"
  ini({
    lka  <- 0.45 ; label("Absorption rate (Ka)")
    lvc  <- 3.45 ; label("Central volume of distribution (Vc)")
    lkel <- 0.534; label("elimination rate (1/d)")
    lk12 <- 0.48; label("central-to-peripheral rate (1/d)")
    lk21 <- 0.34; label("peripheral-to-central rate (1/d)")
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
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    IC50 <- exp(lIC50)
    kin <- exp(lkin)
    kout <- exp(lkout)
    fdepot   <- exp(lfdepot)

    
    d/dt(depot)      <- -ka*depot
    f(depot)         <- fdepot
    d/dt(central)    <- ka*depot -(kel)*central- k12*central + k21*peripheral1
    d/dt(peripheral1) <- k12*central - k21*peripheral1
    d/dt(effect) <- kin - kout*(1-Cc/(Cc + IC50))*effect
    
    Cc <-  central/vc
    Cc ~ prop(propSd)
  })
}
