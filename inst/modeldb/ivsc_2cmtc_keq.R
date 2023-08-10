ivsc_2cmtc_keq <- function() {
  description <- "Two compartment TMDD model with QSS approximation. Parameterized with rate constants "
  ini({
    lka <- 0.45  ; label("Absorption rate (Ka)")
    lvc <- 0.45  ; label("Central volume of distribution (Vc)")
    lkel  <- 0.5 ; label("elimination rate (1/d)")
    lk12 <- 0.04 ; label("central-to-peripheral rate (1/d)")
    lk21  <- 0.1 ; label("peripheral-to-central rate (1/d)")
    lksyn <- 0.1 ; label("Target Synthesis rate (1/d)") 
    lkeT  <- 0.1 ; label("Target elimination rate (1/d)")
    lkeDT <- 0.4 ; label("Internalisation rate (1/d)")
    lkeq <- 0.35 ; label("")
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
    keq <- exp(lkeq)
    fdepot <- exp(lfdepot)
    
    B <- Cc  - Ttot - keq
    D <- 0.5*( B + (B^2 + 4*keq*Cc)^0.5 )
    Ac<- D*vc
    DT<- Ttot*D/(keq+D) 
    T <- Ttot - DT    
    
    d/dt(depot)     <- -ka *depot
    f(depot)        <- fdepot
    d/dt(central)  <- ka *Ad - k12*Ac + k21*Ap - kel*Ac - keDT*DT*vc
    d/dt(Ap)     <- k12*Ac - k21*Ap
    d/dt(Ttot)   <- ksyn - keT *T  - keDT*DT
    
    Cc <- central/vc
    Cc ~ prop(propSd)
  })
}