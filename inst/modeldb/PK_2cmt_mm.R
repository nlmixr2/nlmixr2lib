PK_2cmt_mm <- function() {
  description <- "Two compartment PK model with Michealis-Menten clearance using differential equations"
  ini({
    lka <- 0.45 ; label("Absorption rate (Ka)")
    lvm <- 0.04; label("maximum target-mediated rate of elimination (mg/L/d)")
    lkm <- 0.01; label("Michaelis-Menten constant (mg/L)")
    lvc  <- 3 ; label("Central volume of distribution (V)")
    lvp  <- 5 ; label("Peripheral volume of distribution (Vp)")
    lq  <- 0.1 ; label("Intercompartmental clearance (Q)")
    propSd <- 0.5 ; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    km <- exp(lkm)
    vm <- exp(lvm)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q  <- exp(lq)
    
    k12 <- q/vc
    k21 <- q/vp
    
    d/dt(depot) <- -ka*depot
    d/dt(central) <-  ka*depot - (vm/(km + central/vc))*central - k12*central + k21*peripheral1
    d/dt(peripheral1) <- k12*central - k21*peripheral1
    Cc <- central / vc
    
    Cc ~ prop(propSd)
  })
}
