Liu_2024_8Omethylretusin_rat <- function() {
  description <- paste(
    "Preclinical (rat). Two-compartment oral pharmacokinetic model with",
    "first-order absorption and no lag time for 8-O-methylretusin, one of four",
    "active constituents of Spatholobi Caulis (SPC) quantified in rat plasma",
    "after a single gavage dose of SPC aqueous extract at 60 g crude drug/kg",
    "body weight, equivalent to 0.22 mg/kg of 8-O-methylretusin (Liu 2024).",
    "Parameters come from the Supplementary Materials Table S5, where Drug and",
    "Statistics (DAS) v2.0 fitted each of five rats individually and the table",
    "reports the arithmetic mean of the per-rat estimates. CAVEAT: this is the",
    "one constituent of the four whose published parameters do not reproduce",
    "the authors' own observed data - the apparent central volume V1/F of",
    "12.043 L/kg caps the attainable concentration at 0.22/12.043 = 18.3 ug/L",
    "even under instantaneous complete absorption, yet the observed mean peak",
    "in Table S2 is 22.9 ug/L, and the simulated Cmax falls 3.1-fold short of",
    "the reported Cmax. The parameters are nevertheless encoded exactly as",
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
    depot       = list(analyte = "8-O-methylretusin", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "8-O-methylretusin", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "8-O-methylretusin", units = "mg/kg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species       = "rat (strain not reported)",
    n_subjects    = 5,
    n_studies     = 1,
    disease_state = "healthy",
    dose_range    = "single oral gavage of 1.5 mL Spatholobi Caulis suspension at 60 g crude drug/kg body weight, delivering 0.22 mg/kg of 8-O-methylretusin",
    regions       = "China",
    notes         = paste(
      "Liu 2024 Section 4.2.5 and Supplementary Materials Table S2: five rats",
      "('Parallel rats: 5') dosed by gavage, with 300 uL orbital venous blood",
      "sampled at 0.083, 0.25, 0.5, 0.75, 1.0, 1.5, 2, 3, 4 and 6 h and assayed",
      "by UFLC-MS/MS; 8-O-methylretusin was the only constituent whose plasma",
      "concentrations fell below quantifiable levels by 8 h. The dosing",
      "experiment followed the authors' earlier study (Liu 2021, J Pharm Biomed",
      "Anal 204:114267), which is where the strain, sex, age and body weight are",
      "recorded; none of them are restated in Liu 2024. Ethics approval:",
      "Experimental Animal Protection and Use Committee of Peking University,",
      "No. LA2021275."
    )
  )

  ini({
    # Two-compartment oral disposition fitted by DAS 2.0 to the 60 g crude
    # drug/kg validation experiment (Liu 2024 Supplementary Table S5, mean of
    # five individually fitted rats). Ka is taken directly as published and K12 / K21
    # are converted below; CL/F and V1/F are the table's macro parameters. For this
    # constituent the macro and micro sets are NOT coherent - K10 * V1/F =
    # 80.34 against the published CL/F = 12.552, a 6.4-fold disagreement - and
    # the published SDs show why: K10 is 6.671 (SD 10.483) and V1/F is 12.043
    # (SD 17.445), so the arithmetic mean of five unstable individual fits does
    # not describe any rat. The macro set is used because it is the canonical
    # nlmixr2lib parameterisation and because it scores better against the
    # observed data (simulated Cmax fold-error 3.13 vs 7.87 for the micro set).
    lka  <- fixed(log(2.056));  label("Absorption rate constant Ka (1/h)")  # Table S5, Ka = 2.056 1/h (SD 2.947)
    lcl  <- fixed(log(12.552)); label("Apparent clearance CL/F (L/h/kg)")  # Table S5, CL/F = 12.552 L/h/kg (SD 5.493)
    lvc  <- fixed(log(12.043)); label("Apparent central volume V1/F (L/kg)")  # Table S5, V1/F = 12.043 L/kg (SD 17.445)
    # DAS 2.0 reports the peripheral transfer as micro-constants. They are
    # converted to the canonical q / vp form here, with the arithmetic left
    # visible: Q/F = K12 * V1/F and V2/F = Q/F / K21. This is not cosmetic --
    # a two-compartment model written with k12 / k21 is silently reduced to
    # ONE compartment by rxSolve()'s default useLinCmt = TRUE, which drops
    # peripheral1 and shortens the terminal half-life several-fold.
    lq   <- fixed(log(0.394 * 12.043)); label("Apparent intercompartmental clearance Q/F (L/h/kg)")  # Table S5, K12 = 0.394 1/h (SD 0.665) times V1/F = 12.043 L/kg
    lvp  <- fixed(log(0.394 * 12.043 / 0.474)); label("Apparent peripheral volume V2/F (L/kg)")  # Table S5, Q/F divided by K21 = 0.474 1/h (SD 0.495)

    # Liu 2024 Supplementary Table S5 reports Tlag = 0 (SD 0) for every rat, so
    # no absorption lag is encoded.

    addSd <- fixed(0); label("Additive residual SD on 8-O-methylretusin plasma concentration (ug/mL; not reported)")  # Liu 2024 reports no residual-error model
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
