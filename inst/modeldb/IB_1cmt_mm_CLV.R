IB_1cmt_mm_CLV <- function() {
  description <- "One compartment TMDD model with Michaelis-Menten approximation. Parametrized using clearances and volumes"
  ini({  
    lvc    <- 0.45; label("Central volume of distribution (Vc)")
    lvm    <- 0.5 ; label("maximum target-mediated rate of elimination (mg/L/d)") 
    lkm  <- 0.01; label("Michaelis-Menten constant (mg/L)")
    lcl    <- 0.1 ; label("Clearance (CL)") 
    lfdepot<- 0.4 ; label("Bioavailability (F)")
    propSd <- 0.5 ; label("Proportional residual error (fraction)")
  })
  model({
    vc  <- exp(lvc)
    vm  <- exp(lvm)
    km <- exp(lkm)
    cl  <- exp(lcl)
    fdepot <- exp(lfdepot)
    
    kel <- cl/vc
    k12 <- 0
    k21 <- 0
    ka <- 0
    
    d/dt(depot)      <- -ka*depot
    f(depot)         <- fdepot
    d/dt(central)    <- ka*depot -(vm/(km + central/vc))*central- k12*central + k21*peripheral1 - kel*central
    d/dt(peripheral1)<- k12*central - k21*peripheral1
    
    Cc <-  central/vc
    Cc ~ prop(propSd)
  })
}