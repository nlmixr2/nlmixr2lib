ivsc_2cmt_mm_CLV_Keq <- function() {
  description <- "Two compartment TMDD model with Michaelis-Menten approximation. Parametrized using clearances and volumes"
  ini({
    lka <- 0.45 ; label("Absorption rate (Ka)")
    lvc <- 0.45 ; label("Central volume of distribution (Vc)")
    lvp  <- 0.5 ; label("Peripheral volume of distribution (Vp)")
    lvm <- 0.04; label("maximum target-mediated rate of elimination (mg/L/d)")
    lkeq  <- 0.1;label("Michaelis-Menten constant (mg/L)")
    lcl <- 0.1; label("Clearance (CL)")
    lq  <- 0.1 ;label("Intercompartmental clearance (Q)")
    lfdepot   <- 0.4; label("Bioavailability (F)")
    propSd <- 0.5 ; label("Proportional residual error (fraction)")
  })
  model({
    ka  <- exp(lka)
    vc  <- exp(lvc)
    vp  <- exp(lvp)
    vm  <- exp(lvm)
    keq <- exp(lkeq)
    cl  <- exp(lcl)
    q   <- exp(lq)
    fdepot <- exp(lfdepot)
    
    kel <- cl/vc
    k12 <- q/vc
    k21 <- q/vp
    km <- keq
    
    d/dt(depot)      <- -ka*depot
    f(depot)         <- fdepot
    d/dt(central)    <- ka*depot -(vm/(km + central/vc))*central- k12*central + k21*peripheral1 - kel*central
    d/dt(peripheral1)<- k12*central - k21*peripheral1
    
    
    Cc <-  central/vc
    Cc ~ prop(propSd)
  })
}