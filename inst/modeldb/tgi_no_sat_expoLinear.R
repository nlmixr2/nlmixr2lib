tgi_no_sat_expoLinear <- function() {
  description <- "One compartment TGI model with with exponential tumor growth, without saturation."
  ini({
    lts0 <- 0.8; label("Initial tumor size (TS0)")
    lka <- 0.45 ; label("Absorption rate (Ka)")
    lcl <- 1 ; label("Clearance (CL)")
    lvc  <- 3.45 ; label("Central volume of distribution (V)")
    lkgl <- 0.7; label("Zero-order linear growth rate")
    lkge <- 0.7; label("First-order exponential growth rate")
    propSd <- 0.5 ; label("Proportional residual error (fraction)")
    addSd <- 30 ; label("Additive residual error (tumor volume)")
  })
  model({
    ts0 <- exp(lts0)
    ka <- exp(lka)
    cl <- exp(lcl)
    vc  <- exp(lvc)
    kge <- exp(lkge)
    kgl <- exp(lkgl)
    
    kel <- cl / vc
    tumorSize(0) <- ts0
    tau <- 
    
    d/dt(depot) <- -ka*depot
    d/dt(central) <- ka*depot-kel*central
    d/dt(tumorSize) <- kge*tumorSize
    
    Cc <- central / vc
    Cc ~ prop(propSd)
    tumorSize ~ prop(propSd) + add(addSd)
  })
}
