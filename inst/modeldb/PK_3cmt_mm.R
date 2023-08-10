PK_3cmt_mm <- function() {
  description <- "Three compartment PK model with Michealis-Menten clearance using differential equations"
  ini({
    lka <- 0.45 ; label("Absorption rate (Ka)")
    lvm <- 0.04; label("maximum target-mediated rate of elimination (mg/L/d)")
    lkm <- 0.01; label("Michaelis-Menten constant (mg/L)")
    lvc  <- 3 ; label("Central volume of distribution (V)")
    lvp  <- 5 ; label("Peripheral volume of distribution (Vp)")
    lvp2  <- 8 ; label("Second peripheral volume of distribution (Vp2)")
    lq  <- 0.1 ; label("Intercompartmental clearance (Q)")
    lq2  <- 0.5 ; label("Second intercompartmental clearance (Q2)")
    propSd <- 0.5 ; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    vm <- exp(lvm)
    km <- exp(lkm)
    vc <- exp(lvc)
    vp <- exp(lvp)
    vp2 <- exp(lvp2)
    q  <- exp(lq)
    q2  <- exp(lq2)
    
    k12 <- q/vc
    k21 <- q/vp
    k13 <- q2/vc
    k31 <- q2/vp2
    
    d/dt(depot) <- -ka*depot
    d/dt(central) <-  ka*depot - (vm/(km + central/vc))*central - k12*central + k21*peripheral1 - k13*central + k31*peripheral2
    d/dt(peripheral1) <- k12*central - k21*peripheral1
    d/dt(peripheral2) <- k13*central - k31*peripheral2
    Cc <- central / vc
    
    Cc ~ prop(propSd)
  })
}
