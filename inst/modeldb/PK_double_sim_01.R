PK_double_sim_01 <- function() {
  description <- "PK double absorption model with simultaneous zero order and first order absorptions"
  ini({
    tk01 <- 0.4 ; label("Zero order absorption rate from first site (K01)")
    lka2 <- 0.45 ; label("First order Absorption rate (Ka)")
    lcl <- 1 ; label("Clearance (CL)")
    lvc  <- 3 ; label("Central volume of distribution (V)")
    propSd <- 0.5 ; label("Proportional residual error (fraction)")
  })
  model({
    k01 <- exp(tk01)
    ka2 <- exp(lka2)
    cl <- exp(lcl)
    vc <- exp(lvc)
    
    kel <- cl/vc
    
    d/dt(depot1) <- -k01
    d/dt(depot2) <- -ka2*depot2
    d/dt(central) <-  k01+ka2*depot2-kel*central 
    
    Cc <- central / vc
    
    Cc ~ prop(propSd)
  })
}
