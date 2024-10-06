tgi_sat_genLogistic <- function() {
  description <- "One compartment TGI model where tumor growth is limited by a loss term, with saturation."
  ini({
    lts0 <- 0.3; label("Initial tumor size (TS0)") 
    ltsmax <- 0.9; label("Maximum tumor size at saturation (TSmax)")
    lka <- 0.45 ; label("Absorption rate (Ka)")
    lcl <- 1 ; label("Clearance (CL)")
    lvc  <- 3.45 ; label("Central volume of distribution (V)")
    lkg <- 0.7; label("Birth rate")
    lkd <- 0.7; label ("Death rate")
    CcpropSd <- 0.5 ; label("PK proportional residual error (fraction)")
    tumorSizepropSd <- 0.5 ; label("Tumor size proportional residual error (fraction)")
    tumorSizeaddSd <- 30 ; label("Tumor size additive residual error (tumor volume)")
  })
  model({
    ts0 <- exp(lts0)
    tsmax <- exp(ltsmax)
    ka <- exp(lka)
    cl <- exp(lcl)
    vc  <- exp(lvc)
    kg <- exp(lkg)
    kd <-exp(lkd)
    
    kel <- cl / vc
    tumorSize(0) <- ts0
    
    
    d/dt(depot) <- -ka*depot
    d/dt(central) <- ka*depot-kel*central
    d/dt(tumorSize) <- kg*tumorSize^(2/3)-kd*tumorSize
    
    Cc <- central / vc
    Cc ~ prop(CcpropSd)
    tumorSize ~ prop(tumorSizepropSd) + add(tumorSizeaddSd)
  })
}
