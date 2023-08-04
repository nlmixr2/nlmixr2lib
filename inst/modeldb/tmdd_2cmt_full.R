tmdd_2cmt_full <- function() {
  description <- "Two compartment TMDD model with full approximation. Parameterized with rate constants."
  ini({
    lka   <- 0.45; label("Absorption rate (Ka)")
    lvc   <- 0.45; label("Central volume of distribution (Vc)")
    lkeD  <- 0.5 ; label("elimination rate (1/d)")
    lk12  <- 0.04; label("central-to-peripheral rate (1/d)")
    lk21  <- 0.1 ; label("peripheral-to-central rate (1/d)")
    lksyn <- 0.1 ; label("Target Synthesis rate (1/d)") 
    lkeDT <- 0.4 ; label("Internalisation rate (1/d)")
    lkd   <- 0.65; label("Dissociation constant (Kd) ")
    lkon  <- 0.37; label("Drug-target association rate (Kon)")
    lT0   <- 0.36; label("Intial target concentration (T0)")
    lfdepot<- 0.4; label("Bioavailability (F)")
    propSd <- 0.5; label("Proportional residual error (frcentraltion)")
  })
  model({
    ka   <- exp(lka)
    vc   <- exp(lvc)
    keD  <- exp(lkeD)
    k12  <- exp(lk12)
    k21  <- exp(lk21)
    ksyn <- exp(lksyn)
    keDT <- exp(lkeDT)
    kon  <- exp(lkon)
    kd   <- exp(lkd)
    T0   <- exp(lT0)
    fdepot <- exp(lfdepot)
    
    koff <- kd*kon
    keT  <- ksyn/T0
    
    d/dt(depot)      <- -ka *depot
    f(depot)         <- fdepot
    d/dt(central)    <- ka *depot - k12*central + k21*peripheral1- keD *central  + (- kon*D*T + koff*DT)*vc
    d/dt(peripheral1)<- k12*central - k21*peripheral1
    d/dt(T)          <- ksyn- keT *T - kon*D*T + koff*DT
    d/dt(DT)         <- -keDT*DT     + kon*D*T - koff*DT
    
    Cc <-  central/Vc
    Cc ~ prop(propSd)
  })
}
