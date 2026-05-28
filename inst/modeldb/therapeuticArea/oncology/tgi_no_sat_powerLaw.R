tgi_no_sat_powerLaw <- function() {
  description <- "One compartment TGI model with with exponential tumor growth, without saturation."
  reference <- "nlmixr2lib template"
  units <- list(time = "time_unit", dosing = "dose_unit", concentration = "conc_unit/vol_unit")
  ini({
    lrbase <- 0.8; label("Initial tumor size (TS0)")
    lka <- 0.45 ; label("Absorption rate (Ka)")
    lcl <- 1 ; label("Clearance (CL)")
    lvc  <- 3.45 ; label("Central volume of distribution (V)")
    lgamma <- 0.95; label("proliferative cells as a fraction of the full tumor volume (gamma)")
    lkgl <- 0.7; label("Zero-order linear growth rate")
    lkge <- 0.7; label("First-order exponential growth rate")
    propSd <- 0.5 ; label("PK proportional residual error (fraction)")
    propSd_tumorSize <- 0.5 ; label("Tumor size proportional residual error (fraction)")
    addSd_tumorSize <- 30 ; label("Tumor size additive residual error (tumor volume)")
  })
  model({
    rbase <- exp(lrbase)
    ka <- exp(lka)
    cl <- exp(lcl)
    vc  <- exp(lvc)
    gamma<- exp(lgamma)
    kge <- exp(lkge)
    kgl <- exp(lkgl)
    
    kel <- cl / vc
    tumorSize(0) <- rbase
    
    d/dt(depot) <- -ka*depot
    d/dt(central) <- ka*depot-kel*central
    d/dt(tumorSize) <- kge*tumorSize^gamma
    
    Cc <- central / vc
    Cc ~ prop(propSd)
    tumorSize ~ prop(propSd_tumorSize) + add(addSd_tumorSize)
  })
}
