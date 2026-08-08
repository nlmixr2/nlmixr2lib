tgi_sat_Gompertz<- function() {
  description <- "One compartment TGI model with tumor growth proportional to tumor size through a generalized logistic function, with saturation."
  reference <- "nlmixr2lib template"
  units <- list(time = "time_unit", dosing = "dose_unit", concentration = "conc_unit/vol_unit")
  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    depot      = list(analyte = "drug", units = NA_character_, specimen = "administration site", verified = FALSE),
    central    = list(analyte = "drug", units = NA_character_, specimen = "plasma", verified = FALSE),
    tumor_size = list(analyte = "tumour cells", units = NA_character_, specimen = "not applicable", verified = FALSE)
  )

  ini({
    lrbase <- 0.3; label("Initial tumor size (TS0)") 
    ltsmax <- 0.9; label("Maximum tumor size at saturation (TSmax)")
    lka <- 0.45 ; label("Absorption rate (Ka, 1/time_unit)")
    lcl <- 1 ; label("Clearance (CL)")
    lvc  <- 3.45 ; label("Central volume of distribution (V)")
    lkgl <- 0.7; label("Zero-order linear growth rate (tumor size/time_unit)")
    lalpha <- 0.6; label("parameter one")
    lbeta <- 0.8; label("parameter two")
    propSd <- 0.5 ; label("PK proportional residual error (fraction)")
    propSd_tumor_size <- 0.5 ; label("Tumor size proportional residual error (fraction)")
    addSd_tumor_size <- 30 ; label("Tumor size additive residual error (tumor volume)")
  })
  model({
    rbase <- exp(lrbase)
    tsmax <- exp(ltsmax)
    ka <- exp(lka)
    cl <- exp(lcl)
    vc  <- exp(lvc)
    kgl <- exp(lkgl)
    alpha <- exp(lalpha)
    beta <- exp(lbeta)
    
    kel <- cl / vc
    tumor_size(0) <- rbase

    d/dt(depot) <- -ka*depot
    d/dt(central) <- ka*depot-kel*central
    d/dt(tumor_size) <- (alpha-beta*log(tumor_size))*tumor_size
    
    Cc <- central / vc
    Cc ~ prop(propSd)
    tumor_size ~ prop(propSd_tumor_size) + add(addSd_tumor_size)
  })
}
