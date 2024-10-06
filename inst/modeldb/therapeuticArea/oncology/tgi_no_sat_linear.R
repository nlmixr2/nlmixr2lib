tgi_no_sat_linear <- function() {
  description <- "One compartment TGI model with with linear tumor growth, without saturation."
  ini({
    lts0 <- 0.8; label("Initial tumor size (TS0)") 
    lka <- 0.45 ; label("Absorption rate (Ka)")
    lcl <- 1 ; label("Clearance (CL)")
    lvc  <- 3.45 ; label("Central volume of distribution (V)")
    lkgl <- 0.7; label("Zero-order linear growth rate")
    CcpropSd <- 0.5 ; label("PK proportional residual error (fraction)")
    tumorSizepropSd <- 0.5 ; label("Tumor size proportional residual error (fraction)")
    tumorSizeaddSd <- 30 ; label("Tumor size additive residual error (tumor volume)")
  })
  model({
    ts0 <- exp(lts0)
    ka <- exp(lka)
    cl <- exp(lcl)
    vc  <- exp(lvc)
    kgl <- exp(lkgl)
    
    kel <- cl / vc
    tumorSize(0) <- ts0
    
    d/dt(depot) <- -ka*depot
    d/dt(central) <- ka*depot-kel*central
    d/dt(tumorSize) <- kgl
    
    Cc <- central / vc
    Cc ~ prop(CcpropSd)
    tumorSize ~ prop(tumorSizepropSd) + add(tumorSizeaddSd)
  })
}
