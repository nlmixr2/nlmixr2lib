Janssen_2017_midazolam <- function() {
  description <- paste(
    "One-compartment population PK model for intravenous midazolam used as",
    "a CYP3A metabolic-phenotyping probe in 10 men with metastatic",
    "castration-resistant prostate cancer (Janssen 2017 Table 2A). Each",
    "patient received a single 2.5 mg IV midazolam bolus 1-7 days before",
    "their scheduled cabazitaxel infusion; the empirical Bayes estimates",
    "of individual midazolam clearance from this fit are then used as",
    "covariate input to the companion Janssen 2017 cabazitaxel model",
    "(see modellib('Janssen_2017_cabazitaxel'))."
  )
  reference <- paste(
    "Janssen A, Verkleij CPM, van der Vlist A, Mathijssen RHJ,",
    "Bloemendal HJ, ter Heine R (2017).",
    "Towards better dose individualisation: metabolic phenotyping to",
    "predict cabazitaxel pharmacokinetics in men with prostate cancer.",
    "Br J Cancer 116(10):1312-1317.",
    "doi:10.1038/bjc.2017.91.",
    sep = " "
  )
  vignette <- "Janssen_2017_cabazitaxel_phenotyping"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 10L,
    n_studies      = 1L,
    age_range      = "65-77 years (median 67)",
    weight_range   = "Body surface area 1.76-2.34 m^2 (median 1.95-1.96 m^2); body weight not reported in the paper",
    sex_female_pct = 0,
    disease_state  = "Metastatic castration-resistant prostate cancer; all previously treated with docetaxel; scheduled to receive cabazitaxel as routine palliative care",
    dose_range     = "Single 2.5 mg intravenous bolus of midazolam (CYP3A phenotyping probe), administered 1-7 days before the patient's first cabazitaxel infusion",
    regions        = "Meander Medical Center, Amersfoort, The Netherlands (single-centre)",
    notes          = paste(
      "Six plasma samples per subject collected at 0, 30, 60, 120, 240,",
      "and 360 min post-injection (Methods). Bioanalysis by validated",
      "LC-MS/MS with LOQ = 0.4 ng/mL; inter- and intra-run precision <",
      "12.5%, accuracy within 12% of nominal (Bioanalytics). Hepatic",
      "function (bilirubin) not decreased in any patient. Estimation",
      "with NONMEM 7.3.0, first-order conditional estimation with",
      "interaction; goodness-of-fit and prediction-corrected VPC used",
      "for model evaluation (Methods)."
    )
  )

  ini({
    # Final-model parameter values from Janssen 2017 Table 2A. The paper
    # describes the midazolam disposition as a "first order one compartment
    # linear pharmacokinetic model" (Results, "Pharmacokinetic modeling")
    # fitted to single-IV-bolus data, so the structural model is one
    # compartment with bolus loading of central.
    lcl <- log(26.0)
    label("Midazolam clearance CLmdz (L/h)")            # Table 2A: CLmdz = 26.0 L/h, RSE 25.0%
    lvc <- log(75.4)
    label("Midazolam central volume of distribution V3 (L)")  # Table 2A: V3 = 75.4 L, RSE 24.8%

    # IIV is reported as %CV; convert to internal log-scale variance via
    # omega^2 = log(1 + CV^2).
    #   CL : 73.2% -> log(1 + 0.732^2) = 0.42907
    #   V3 : 77.3% -> log(1 + 0.773^2) = 0.46846
    etalcl ~ 0.42907   # Table 2A: IIV CLmdz 73.2% CV (RSE 53.7%)
    etalvc ~ 0.46846   # Table 2A: IIV V3   77.3% CV (RSE 80.4%)

    # Residual error (proportional). The paper does not specify the residual
    # form explicitly, but reports the residual variability as a single CV%
    # entry under Final model -- consistent with the standard proportional
    # error structure that NONMEM applies on the linear-concentration scale.
    propSd <- 0.354
    label("Proportional residual error for midazolam plasma concentration (fraction)")  # Table 2A: residual error 35.4% CV (RSE 40.6%)
  })

  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    # One-compartment IV-bolus disposition. Each patient received a single
    # 2.5 mg IV midazolam shot at time 0, followed by a 5 mL NaCl 0.9%
    # flush (Methods). Data: 6 plasma samples at 0, 30, 60, 120, 240, 360
    # min after the injection.
    Cc <- linCmt()
    Cc ~ prop(propSd)
  })
}
