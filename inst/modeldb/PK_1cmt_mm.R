PK_1cmt_mm <- function() {
  description <- "One compartment PK model with Michaelis-Menten clearance using differential equations"
  ini({
    lka <- 0.45 ; label("Absorption rate (Ka)")
    lvm <- 0.04; label("maximum target-mediated rate of elimination (mg/L/d)")
    Km <- 0.01; label("Michaelis-Menten constant (mg/L)")
    lvc  <- 3.45 ; label("Central volume of distribution (V)")
    propSd <- 0.5 ; label("Proportional residual error (fraction)")
  })
  model({
    ka <- exp(lka)
    Vm <- exp(lvm)
    vc  <- exp(lvc)
    
    d/dt(depot) <- -ka*depot
    d/dt(central) <- ka*depot-central*(Vm/(Km + central/vc))
    
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}