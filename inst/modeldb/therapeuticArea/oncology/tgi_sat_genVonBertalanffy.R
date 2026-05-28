tgi_sat_genVonBertalanffy <- function() {
  description <- "One compartment TGI model where tumor growth is limited by a loss term, with saturation."
  reference <- "nlmixr2lib template"
  units <- list(time = "time_unit", dosing = "dose_unit", concentration = "conc_unit/vol_unit")
  ini({
    lrbase <- 0.3; label("Initial tumor size (TS0)") 
    ltsmax <- 0.9; label("Maximum tumor size at saturation (TSmax)")
    lka <- 0.45 ; label("Absorption rate (Ka)")
    lcl <- 1 ; label("Clearance (CL)")
    lvc  <- 3.45 ; label("Central volume of distribution (V)")
    lkg <- 0.7; label("Birth rate")
    lkd <- 0.7; label ("Death rate")
    lgamma <- 0.95; label("proliferative cells as a fraction of the full tumor volume (gamma)")
    propSd <- 0.5 ; label("PK proportional residual error (fraction)")
    propSd_tumorSize <- 0.5 ; label("Tumor size proportional residual error (fraction)")
    addSd_tumorSize <- 30 ; label("Tumor size additive residual error (tumor volume)")
  })
  model({
    rbase <- exp(lrbase)
    tsmax <- exp(ltsmax)
    ka <- exp(lka)
    cl <- exp(lcl)
    vc  <- exp(lvc)
    kg <- exp(lkg)
    kd <-exp(lkd)
    gamma <- exp(lgamma)
    
    kel <- cl / vc
    tumorSize(0) <- rbase
    
    
    d/dt(depot) <- -ka*depot
    d/dt(central) <- ka*depot-kel*central
    d/dt(tumorSize) <- kg*tumorSize^(gamma)-kd*tumorSize
    
    Cc <- central / vc
    Cc ~ prop(propSd)
    tumorSize ~ prop(propSd_tumorSize) + add(addSd_tumorSize)
  })
}
