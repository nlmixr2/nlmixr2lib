tgi_no_sat_powerLaw <- function() {
  description <- "One compartment TGI model with with exponential tumor growth, without saturation."
  ini({
    lts0 <- 0.8; label("Initial tumor size (TS0)")
    lka <- 0.45 ; label("Absorption rate (Ka)")
    lcl <- 1 ; label("Clearance (CL)")
    lvc  <- 3.45 ; label("Central volume of distribution (V)")
    lgamma <- 0.95; label("proliferative cells as a fraction of the full tumor volume (gamma)")
    lkgl <- 0.7; label("Zero-order linear growth rate")
    lkge <- 0.7; label("First-order exponential growth rate")
    CcpropSd <- 0.5 ; label("PK proportional residual error (fraction)")
    tumorSizepropSd <- 0.5 ; label("Tumor size proportional residual error (fraction)")
    tumorSizeaddSd <- 30 ; label("Tumor size additive residual error (tumor volume)")
  })
  model({
    ts0 <- exp(lts0)
    ka <- exp(lka)
    cl <- exp(lcl)
    vc  <- exp(lvc)
    gamma<- exp(lgamma)
    kge <- exp(lkge)
    kgl <- exp(lkgl)
    
    kel <- cl / vc
    tumorSize(0) <- ts0
    
    d/dt(depot) <- -ka*depot
    d/dt(central) <- ka*depot-kel*central
    d/dt(tumorSize) <- kge*tumorSize^gamma
    
    Cc <- central / vc
    Cc ~ prop(CcpropSd)
    tumorSize ~ prop(tumorSizepropSd) + add(tumorSizeaddSd)
  })
}
