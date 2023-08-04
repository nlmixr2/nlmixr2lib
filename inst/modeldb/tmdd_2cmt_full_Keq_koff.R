tmdd_2cmt_full_Keq_koff <- function() {
  description <- "Two compartment TMDD model with full approximation. Parameterized with rate constants "
  ini({
    lka <- 0.45  ; label("Absorption rate (Ka)")
    lvc <- 0.45  ; label("Central volume of distribution (Vc)")
    lkel  <- 0.5 ; label("elimination rate (1/d)")
    lk12 <- 0.04 ; label("central-to-peripheral rate (1/d)")
    lk21  <- 0.1 ; label("peripheral-to-central rate (1/d)")
    lksyn <- 0.1 ; label("Target Synthesis rate (1/d)") 
    lkeT  <- 0.1 ; label("Target elimination rate (1/d)")
    lkeDT <- 0.4 ; label("Internalisation rate (1/d)")
    lkoff <- 0.35 ; label("Drug-target dissociation rate (Koff)")
    lkeq <- 0.37 ; label("")
    lfdepot <- 0.4; label("Bioavailability (F)")
    propSd <- 0.5 ; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    vc <- exp(lvc)
    kel <- exp(lkel)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    ksyn  <- exp(lksyn)
    keT <- exp(lkeT)
    keDT <- exp(lkeDT)
    koff <- exp(lkoff)
    keq <- exp(lkeq)
    fdepot <- exp(lfdepot)
    
    kon <- koff/keq
    
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
