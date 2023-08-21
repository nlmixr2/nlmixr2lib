tmdd_1cmt_mm_CLV_Keq <- function() {
  description <- "One compartment TMDD model with Michaelis-Menten approximation. Parametrized using clearances and volumes"
  ini({  
    lka    <- 0.48; label("Absorption rate (Ka)")
    lvc    <- 0.45; label("Central volume of distribution (Vc)")
    lvm    <- 0.5 ; label("maximum target-mediated rate of elimination (mg/L/d)") 
    lkeq   <- 0.1 ; label ("Michaelis-Menten constant (mg/L)") 
    lcl    <- 0.1 ; label("Clearance (CL)") 
    lfdepot<- 0.4 ; label("Bioavailability (F)")
    propSd <- 0.5 ; label("Proportional residual error (fraction)")
  })
  model({
    ka  <- exp(lka)
    vc  <- exp(lvc)
    vm  <- exp(lvm)
    keq <- exp(lkeq)
    cl  <- exp(lcl)
    fdepot <- exp(lfdepot)
    
    kel <- cl/vc
    k12 <- 0
    k21 <- 0
    km <- keq
    
    d/dt(depot)      <- -ka*depot
    f(depot)         <- fdepot
    d/dt(central)    <- ka*depot -(vm/(km + central/vc))*central- k12*central + k21*peripheral1 - kel*central
    d/dt(peripheral1)<- k12*central - k21*peripheral1
    
    Cc <-  central/vc
    Cc ~ prop(propSd)
  })
}
