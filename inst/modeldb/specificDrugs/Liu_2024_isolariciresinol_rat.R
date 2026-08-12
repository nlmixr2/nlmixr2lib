Liu_2024_isolariciresinol_rat <- function() {
  description <- paste(
    "Preclinical (rat). Two-compartment oral pharmacokinetic model with",
    "first-order absorption and no lag time for isolariciresinol, one of four",
    "active constituents of Spatholobi Caulis (SPC) quantified in rat plasma",
    "after a single gavage dose of SPC aqueous extract at 60 g crude drug/kg",
    "body weight, equivalent to 2.62 mg/kg of isolariciresinol (Liu 2024).",
    "Parameters come from the Supplementary Materials Table S9, where Drug and",
    "Statistics (DAS) v2.0 fitted each of five rats individually and the table",
    "reports the arithmetic mean of the per-rat estimates. CAVEAT: the",
    "published apparent clearance CL/F of 99.561 L/h/kg carries an SD of",
    "130.836 - larger than the mean - and is irreconcilable with the same",
    "table's own AUC(0-inf); dose/CL/F gives 26.3 ug/L*h against a reported",
    "52.8 (compartmental) and 59.3 (statistical moment), so simulated exposure",
    "is ~2.3-fold low. The parameters are nevertheless encoded exactly as",
    "published; see the vignette Assumptions section. Clearance and central",
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
    depot       = list(analyte = "isolariciresinol", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "isolariciresinol", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "isolariciresinol", units = "mg/kg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species       = "rat (strain not reported)",
    n_subjects    = 5,
    n_studies     = 1,
    disease_state = "healthy",
    dose_range    = "single oral gavage of 1.5 mL Spatholobi Caulis suspension at 60 g crude drug/kg body weight, delivering 2.62 mg/kg of isolariciresinol",
    regions       = "China",
    notes         = paste(
      "Liu 2024 Section 4.2.5 and Supplementary Materials Table S2: five rats",
      "('Parallel rats: 5') dosed by gavage, with 300 uL orbital venous blood",
      "sampled at 0.083, 0.25, 0.5, 0.75, 1.0, 1.5, 2, 3, 4, 6 and 8 h and",
      "assayed by UFLC-MS/MS. The dosing experiment followed the authors'",
      "earlier study (Liu 2021, J Pharm Biomed Anal 204:114267), which is where",
      "the strain, sex, age and body weight are recorded; none of them are",
      "restated in Liu 2024. Ethics approval: Experimental Animal Protection and",
      "Use Committee of Peking University, No. LA2021275."
    )
  )

  ini({
    # Two-compartment oral disposition fitted by DAS 2.0 to the 60 g crude
    # drug/kg validation experiment (Liu 2024 Supplementary Table S9, mean of
    # five individually fitted rats). Ka is taken directly as published and K12 / K21
    # are converted below; CL/F and V1/F are the table's macro parameters. For this
    # constituent the macro and micro sets are NOT coherent - K10 * V1/F =
    # 222.15 against the published CL/F = 99.561, a 2.2-fold disagreement - and
    # the published SDs show why: K10 is 7.512 (SD 14.088) and CL/F is 99.561
    # (SD 130.836), so the arithmetic mean of five individual fits is dominated
    # by one extreme rat. The macro set is used because it is the canonical
    # nlmixr2lib parameterisation and because it scores better against the
    # observed data (simulated AUC(0-inf) fold-error 2.25 vs 5.03 for the micro
    # set). The same table's statistical-moment clearance CLz/F = 44.282
    # L/h/kg (SD 2.557) reproduces the observed AUC exactly, but it is an NCA
    # quantity rather than part of the fitted compartmental model, so it is
    # reported in the vignette rather than substituted here.
    lka  <- fixed(log(5.032));  label("Absorption rate constant Ka (1/h)")  # Table S9, Ka = 5.032 1/h (SD 1.943)
    lcl  <- fixed(log(99.561)); label("Apparent clearance CL/F (L/h/kg)")  # Table S9, CL/F = 99.561 L/h/kg (SD 130.836)
    lvc  <- fixed(log(29.573)); label("Apparent central volume V1/F (L/kg)")  # Table S9, V1/F = 29.573 L/kg (SD 11.601)
    # DAS 2.0 reports the peripheral transfer as micro-constants. They are
    # converted to the canonical q / vp form here, with the arithmetic left
    # visible: Q/F = K12 * V1/F and V2/F = Q/F / K21. This is not cosmetic --
    # a two-compartment model written with k12 / k21 is silently reduced to
    # ONE compartment by rxSolve()'s default useLinCmt = TRUE, which drops
    # peripheral1 and shortens the terminal half-life several-fold.
    lq   <- fixed(log(0.448 * 29.573)); label("Apparent intercompartmental clearance Q/F (L/h/kg)")  # Table S9, K12 = 0.448 1/h (SD 0.297) times V1/F = 29.573 L/kg
    lvp  <- fixed(log(0.448 * 29.573 / 0.247)); label("Apparent peripheral volume V2/F (L/kg)")  # Table S9, Q/F divided by K21 = 0.247 1/h (SD 0.111)

    # Liu 2024 Supplementary Table S9 reports Tlag = 0 (SD 0) for every rat, so
    # no absorption lag is encoded.

    addSd <- fixed(0); label("Additive residual SD on isolariciresinol plasma concentration (ug/mL; not reported)")  # Liu 2024 reports no residual-error model
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
