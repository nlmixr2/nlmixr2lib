PK_double_seq_01 <- function() {
  description <- "PK double absorption model with sequential zero order absorption followed by first order absorption"
  ini({
    tk01 <- 0.4 ; label("Zero order absorption rate from first site (K01)")
    ka2 <- 0.45 ; label("First order Absorption rate (Ka)")
    lcl <- 1 ; label("Clearance (CL)")
    lvc  <- 3 ; label("Central volume of distribution (V)")
    propSd <- 0.5 ; label("Proportional residual error (fraction)")
  })
  model({
    k01 <- exp(tk01)
    k02 <- exp(tk02)
    cl <- exp(lcl)
    vc <- exp(lvc)
    
    kel <- cl/vc
    lag2 <- k01
    
    d/dt(depot1) <- -k01
    d/dt(depot2) <- k01-ka2*depot2
    alag(depot2)  <- lag2
    d/dt(central) <-  ka2*depot2 - kel*central 
    
    Cc <- central / vc
    
    Cc ~ prop(propSd)
  })
}
