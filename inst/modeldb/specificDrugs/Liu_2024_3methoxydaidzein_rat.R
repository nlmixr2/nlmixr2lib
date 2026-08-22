Liu_2024_3methoxydaidzein_rat <- function() {
  description <- paste(
    "Preclinical (rat). Two-compartment oral pharmacokinetic model with",
    "first-order absorption and no lag time for 3'-methoxydaidzein, one of four",
    "active constituents of Spatholobi Caulis (SPC) quantified in rat plasma",
    "after a single gavage dose of SPC aqueous extract at 60 g crude drug/kg",
    "body weight, equivalent to 0.124 mg/kg of 3'-methoxydaidzein (Liu 2024).",
    "Parameters come from the Supplementary Materials Table S3, where Drug and",
    "Statistics (DAS) v2.0 fitted each of five rats individually and the table",
    "reports the arithmetic mean of the per-rat estimates. Clearance and central",
    "volume are apparent (CL/F, V1/F) because the route is oral and no",
    "bioavailability was determined. The paper's headline model is a whole-body",
    "GastroPlus PBPK model; that layer is NOT reproduced here because the",
    "tissue volumes, blood flows, ACAT absorption parameters and enterohepatic",
    "recirculation settings are GastroPlus built-ins that the paper never",
    "prints. No between-subject variability or residual error was reported;",
    "every parameter is fixed at the published mean and the residual SD is",
    "fixed at zero."
  )
  reference <- paste(
    "Liu X, Du R, Zhang T, Li Y, Li L, Yang Z, Zhang Y, Wang Q (2024).",
    "Predicting pharmacokinetics of active constituents in Spatholobi caulis by",
    "using physiologically based pharmacokinetic models.",
    "Pharmaceuticals (Basel) 17(12):1621. doi:10.3390/ph17121621"
  )
  vignette <- "Liu_2024_spatholobiCaulis"
  units <- list(time = "h", dosing = "mg/kg", concentration = "ug/mL")

  covariateData <- list()

  compartmentData <- list(
    depot       = list(analyte = "3'-methoxydaidzein", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "3'-methoxydaidzein", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "3'-methoxydaidzein", units = "mg/kg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species       = "rat (strain not reported)",
    n_subjects    = 5,
    n_studies     = 1,
    disease_state = "healthy",
    dose_range    = "single oral gavage of 1.5 mL Spatholobi Caulis suspension at 60 g crude drug/kg body weight, delivering 0.124 mg/kg of 3'-methoxydaidzein",
    regions       = "China",
    notes         = paste(
      "Liu 2024 Section 4.2.5 and Supplementary Materials Table S2: five rats",
      "('Parallel rats: 5') dosed by gavage, with 300 uL orbital venous blood",
      "sampled at 0.083, 0.25, 0.5, 0.75, 1.0, 1.5, 2, 3, 4, 6, 8 and 12 h and",
      "assayed by UFLC-MS/MS. The dosing experiment followed the authors'",
      "earlier study (Liu 2021, J Pharm Biomed Anal 204:114267), which is where",
      "the strain, sex, age and body weight are recorded; none of them are",
      "restated in Liu 2024. Ethics approval: Experimental Animal Protection and",
      "Use Committee of Peking University, No. LA2021275."
    )
  )

  ini({
    # Two-compartment oral disposition fitted by DAS 2.0 to the 60 g crude
    # drug/kg validation experiment (Liu 2024 Supplementary Table S3, mean of
    # five individually fitted rats). Ka is taken directly as published and K12 / K21
    # are converted below; CL/F and V1/F are the table's macro parameters, which for this
    # constituent are coherent with the micro-constants (K10 * V1/F = 6.077 vs
    # the published CL/F = 6.06, a 0.3% difference).
    lka  <- fixed(log(4.239));  label("Absorption rate constant Ka (1/h)")  # Table S3, Ka = 4.239 1/h (SD 1.338)
    lcl  <- fixed(log(6.06));   label("Apparent clearance CL/F (L/h/kg)")  # Table S3, CL/F = 6.06 L/h/kg (SD 0.162)
    lvc  <- fixed(log(12.821)); label("Apparent central volume V1/F (L/kg)")  # Table S3, V1/F = 12.821 L/kg (SD 0.71)
    # DAS 2.0 reports the peripheral transfer as micro-constants. They are
    # converted to the canonical q / vp form here, with the arithmetic left
    # visible: Q/F = K12 * V1/F and V2/F = Q/F / K21. This is not cosmetic --
    # a two-compartment model written with k12 / k21 is silently reduced to
    # ONE compartment by rxSolve()'s default useLinCmt = TRUE, which drops
    # peripheral1 and shortens the terminal half-life several-fold.
    lq   <- fixed(log(0.207 * 12.821)); label("Apparent intercompartmental clearance Q/F (L/h/kg)")  # Table S3, K12 = 0.207 1/h (SD 0.113) times V1/F = 12.821 L/kg
    lvp  <- fixed(log(0.207 * 12.821 / 0.135)); label("Apparent peripheral volume V2/F (L/kg)")  # Table S3, Q/F divided by K21 = 0.135 1/h (SD 0.064)

    # Liu 2024 Supplementary Table S3 reports Tlag = 0 (SD 0) for every rat, so
    # no absorption lag is encoded.

    addSd <- fixed(0); label("Additive residual SD on 3'-methoxydaidzein plasma concentration (ug/mL; not reported)")  # Liu 2024 reports no residual-error model
  })

  model({
    ka  <- exp(lka)
    cl  <- exp(lcl)
    vc  <- exp(lvc)
    q   <- exp(lq)
    vp  <- exp(lvp)

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - (cl / vc) * central - (q / vc) * central + (q / vp) * peripheral1
    d/dt(peripheral1) <-  (q / vc) * central - (q / vp) * peripheral1

    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
