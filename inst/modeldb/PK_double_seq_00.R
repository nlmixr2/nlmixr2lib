PK_double_seq_00 <- function() {
  description <- "PK double absorption model with sequential zero order absorptions"
  ini({
    tk01 <- 0.4 ; label("Zero order absorption rate (K01)")
    tk02 <- 0.45 ; label("Zero order absorption rate (K02)")
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
    d/dt(depot2) <- -k02
    alag(depot2)  <- lag2
    d/dt(central) <-  k01+k02 - kel*central 
    
    Cc <- central / vc
    
    Cc ~ prop(propSd)
  })
}
