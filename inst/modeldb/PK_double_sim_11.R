PK_double_sim_11 <- function() {
  description <- "PK double absorption model with simultaneous first order absorptions"
  ini({
    lka1 <- 0.45 ; label("First order Absorption rate from first site (Ka)")
    lka2 <- 0.45 ; label("First order Absorption rate from second site (Ka)")
    lcl <- 1 ; label("Clearance (CL)")
    lvc  <- 3 ; label("Central volume of distribution (V)")
    propSd <- 0.5 ; label("Proportional residual error (fraction)")
  })
  model({
    ka1 <- exp(lka1)
    ka2 <- exp(lka2)
    cl <- exp(lcl)
    vc <- exp(lvc)
    
    kel <- cl/vc
    
    d/dt(depot1) <- -ka1*depot1
    d/dt(depot2) <- -ka2*depot2
    d/dt(central) <-  ka1*depot1+ka2*depot2 - kel*central 
    
    Cc <- central / vc
    
    Cc ~ prop(propSd)
  })
}
