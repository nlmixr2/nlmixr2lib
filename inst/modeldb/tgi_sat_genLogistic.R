tgi_sat_genLogistic <- function() {
  description <- "One compartment TGI model with tumor growth proportional to tumor size through a generalized logistic function, with saturation."
  ini({
    lts0 <- 0.3; label("Initial tumor size (TS0)") 
    ltsmax <- 0.9; label("Maximum tumor size at saturation (TSmax)")
    lka <- 0.45 ; label("Absorption rate (Ka)")
    lcl <- 1 ; label("Clearance (CL)")
    lvc  <- 3.45 ; label("Central volume of distribution (V)")
    lkgl <- 0.7; label("Zero-order linear growth rate")
    lgamma <- 0.95; label("proliferative cells as a fraction of the full tumor volume (gamma)")
    propSd <- 0.5 ; label("Proportional residual error (fraction)")
    addSd <- 30 ; label("Additive residual error (tumor volume)")
  })
  model({
    ts0 <- exp(lts0)
    tsmax <- exp(ltsmax)
    ka <- exp(lka)
    cl <- exp(lcl)
    vc  <- exp(lvc)
    kgl <- exp(lkgl)
    gamma <- exp(lgamma)
    
    kel <- cl / vc
    tumorSize(0) <- ts0
    
    
    d/dt(depot) <- -ka*depot
    d/dt(central) <- ka*depot-kel*central
    d/dt(tumorSize) <- kge*tumorSize*(1-(tumorSize/tsmax)^gamma)
    
    Cc <- central / vc
    Cc ~ prop(propSd)
    tumorSize ~ prop(propSd) + add(addSd)
  })
}
