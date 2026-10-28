Liu_2023_bintrafuspAlfa_tgi <- function() {
  description <- paste(
    "Tumour growth inhibition (TGI) model of the sum of longest diameters (SLD) of target lesions in",
    "adults with second-line non-small cell lung cancer treated with bintrafusp alfa 500 or 1200 mg every",
    "2 weeks (phase I trial NCT02517398, 2L NSCLC expansion cohort). Claret-type model with zero-order",
    "tumour growth and first-order drug-induced shrinkage that decays with time since treatment",
    "initiation (resistance); the time axis starts at the screening scan, and shrinkage switches on at",
    "the first dose. Neither dose nor exposure enters the model (no significant dose or clearance effect",
    "was found). IIV on baseline SLD, growth and shrinkage rates; proportional residual error whose",
    "magnitude carries its own between-subject random effect."
  )
  reference <- paste(
    "Liu H, Milenkovic-Grisic AM, Krishnan SM, Jonsson S, Friberg LE, Girard P, Venkatakrishnan K,",
    "Vugmeyster Y, Khandelwal A, Karlsson MO.",
    "A multistate modeling and simulation framework to learn dose-response of oncology drugs:",
    "Application to bintrafusp alfa in non-small cell lung cancer.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(11):1738-1750. doi:10.1002/psp4.12976.",
    "PMCID: PMC10681430. Parameter values from Supplementary Table S2; structure from Appendix S1",
    "'Population tumor growth inhibition model development'."
  )
  vignette <- "Liu_2023_bintrafuspAlfa"
  units <- list(
    time = "week",
    dosing = "n/a (no dosing events; drug effect switches on at T_SCAN_TO_DOSE)",
    concentration = "mm"
  )

  # etaRUV is the random effect on the magnitude of the residual error
  # (Supplementary Table S2 row 'sigma^2 EPS-TS'); see the model() comment.
  paper_specific_etas <- c("etaRUV")

  covariateData <- list(
    T_SCAN_TO_DOSE = list(
      description = "Time from the screening tumour scan (the model's time origin) to the first bintrafusp alfa dose.",
      units = "week",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed per subject, non-negative. Appendix S1: 't is time since screening, T is time since",
        "treatment initiation, which was set to 0 before treatment' and 'kD was fixed to 0 before treatment",
        "initiation'. The interval averaged 17.6 days (range 1.4-48.4 days), i.e. 2.51 weeks (0.2-6.9",
        "weeks) (Liu 2023 Results and Appendix S1). Enter it in WEEKS, the model's time unit. Use 0 when",
        "the simulation axis starts at the first dose and baseline is defined at that time."
      ),
      source_name = "T (time since treatment initiation) = t - T_SCAN_TO_DOSE"
    )
  )

  covariatesDataExcluded <- list(
    LMET = list(
      description = "Liver metastases at baseline (1 = yes).",
      units = "(binary)",
      type = "binary",
      notes = "Identified by the stepwise covariate search on KG (dOFV = -13.53) but driven by a single influential patient (individual dOFV = -13.28); not retained (Appendix S1)."
    ),
    PDL1_TUM = list(
      description = "PD-L1 expression on tumour cells at baseline (%).",
      units = "%",
      type = "continuous",
      notes = "Association with lambda (dOFV = -24.24) gave an unstable model (condition number > 80000); not retained (Appendix S1)."
    )
  )

  compartmentData <- list(
    tumor_size = list(
      analyte = "sum of longest diameters of target lesions (RECIST 1.1, independent review)",
      units = "mm",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 80L,
    n_studies = 1L,
    age_range = "38-85 years",
    weight_range = "41.9-109.0 kg",
    sex_female_pct = 28.75,
    race_ethnicity = c(White = 48.75, Black = 1.25, Asian = 50),
    disease_state = "Advanced non-small cell lung cancer, second line (2L NSCLC)",
    dose_range = "500 or 1200 mg IV every 2 weeks",
    regions = "USA, United Kingdom, Belgium, Italy, Spain, South Korea, Taiwan (25 sites)",
    notes = paste(
      "SLD per independent review was available at screening for 76 of 80 patients (median 57 mm, range",
      "13-262 mm); 73 patients had at least one post-baseline SLD; median tumour-size follow-up 12 weeks",
      "(range 2-177 weeks). Tumour assessments every 6 weeks in the first year and every 12 weeks",
      "thereafter (Liu 2023 Methods and Results; Table 1)."
    )
  )

  ini({
    lrbase <- log(56.3)
    label("Baseline SLD at the screening scan IBASE (mm)") # Supplementary Table S2: IBASE = 56.3 mm (RSE 2.6%)
    lkg <- log(1.43)
    label("Zero-order tumour growth rate KG (mm/week)") # Supplementary Table S2: KG = 1.43 mm/week (RSE 13%)
    lkd <- log(0.0306)
    label("First-order drug-induced tumour shrinkage rate constant kD (1/week)") # Supplementary Table S2: kD = 0.0306 /week (RSE 17%)
    llambda <- log(0.00559)
    label("Rate constant of disappearance of drug effect (resistance) lambda (1/week)") # Supplementary Table S2: lambda = 0.00559 /week (RSE 25%)

    etalrbase ~ 0.373 # Supplementary Table S2: omega^2 IBASE = 0.373
    etalkg ~ 1.00 # Supplementary Table S2: omega^2 KG = 1.00
    etalkd ~ 0.650 # Supplementary Table S2: omega^2 kD = 0.650

    etaRUV ~ 0.0288 # Supplementary Table S2: 'sigma^2 EPS-TS (Variance)' = 0.0288 (RSE 63%), random effect on the residual-error magnitude

    propSd <- 0.153
    label("Proportional residual error on SLD (fraction)") # Supplementary Table S2: TS-proportional residual error (CV) = 0.153
  })

  model({
    rbase <- exp(lrbase + etalrbase)
    kg <- exp(lkg + etalkg)
    kd <- exp(lkd + etalkd)
    lambda <- exp(llambda)

    # Appendix S1: dSLD/dt = KG - kD * exp(-lambda * T) * SLD, SLD(0) = IBASE,
    # with t (the integration time) from the screening scan and T = time since
    # treatment initiation, 0 before the first dose. kD is 0 before the first
    # dose, so the untreated tumour grows at KG over the screening interval.
    onTrt <- t >= T_SCAN_TO_DOSE
    tTrt <- (t - T_SCAN_TO_DOSE) * onTrt

    d/dt(tumor_size) <- kg - kd * onTrt * exp(-lambda * tTrt) * tumor_size
    tumor_size(0) <- rbase

    # Proportional error with a subject-level random effect on its magnitude.
    propSdi <- propSd * exp(etaRUV)
    tumor_size ~ prop(propSdi)
  })
}
