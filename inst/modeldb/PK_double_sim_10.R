PK_double_sim_10 <- function() {
  description <-  "PK double absorption model with simultaneous first order and zero order absorptions"
  ini({
    lka1 <- 0.45 ; label("First order Absorption rate (Ka)")
    tk02 <- 0.4 ; label("Zero order absorption rate from second site (K02)")
    lcl <- 1 ; label("Clearance (CL)")
    lvc  <- 3 ; label("Central volume of distribution (V)")
    propSd <- 0.5 ; label("Proportional residual error (fraction)")
  })
  model({
    ka1 <- exp(lka1)
    k02 <- exp(tk02)
    cl <- exp(lcl)
    vc <- exp(lvc)
    
    kel <- cl/vc
    
    d/dt(depot1) <- -ka1*depot1
    d/dt(depot2) <- -k02
    d/dt(central) <-  ka1*depot1+ k02- kel*central 
    
    Cc <- central / vc
    
    Cc ~ prop(propSd)
  })
}
