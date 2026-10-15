Tran_2026_pioglitazone_alt_mbma <- function() {
  description <- "MBMA. Weibull time-course model of study-arm mean serum alanine aminotransferase (ALT) under pioglitazone or placebo in adults with metabolic dysfunction-associated steatohepatitis (MASH), fitted in NONMEM to arm-mean ALT digitised from eight published trials. ALT(t) = rmin + wamp * exp(-wlambda * t^wgamma): the arm starts at rmin + wamp (73.0 U/L) and declines towards the plateau rmin (33.9 U/L, which the paper calls 'Baseline'). Pioglitazone acts only on the Weibull rate, raising wlambda from 0.0689 (placebo) to 0.534 (any dose; no dose-response was found). Between-study random effects (exponential) on the plateau and the amplitude; combined proportional + additive residual error, reported unweighted -- the source divides the residual by sqrt(N) for an arm of N patients, which downstream simulation of arm means must apply. Suitable for simulating study-arm mean ALT trajectories, not individual patients. Companion models from the same paper: modellib('Tran_2026_pioglitazone_ast_mbma') and the four histology models modellib('Tran_2026_pioglitazone_steatosis_mbma'), modellib('Tran_2026_pioglitazone_inflammation_mbma'), modellib('Tran_2026_pioglitazone_ballooning_mbma') and modellib('Tran_2026_pioglitazone_fibrosis_mbma')."

  reference <- paste(
    "Tran QT, Bui TT, Ngo LT, Yang BR, Baek IH, Nguyen VH, Lee KA, Yun HY,",
    "Chae JW, Lee S, Kim JH, Jung W. Model-Based Meta-Analysis of the",
    "Relationship Between Pioglitazone and Histological Outcomes in Metabolic",
    "Dysfunction-Associated Steatohepatitis Patients.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15:e70034.",
    "doi:10.1002/psp4.70034.",
    sep = " "
  )
  vignette <- "Tran_2026_pioglitazone_mash"

  units <- list(
    time = "month",
    dosing = "(none; pioglitazone treatment enters as the DOSE_PIO_MGD covariate in mg/day, not as dosing events)",
    concentration = "U/L (the observation alt is the study-arm mean serum ALT activity, not a drug concentration)"
  )

  covariateData <- list(
    DOSE_PIO_MGD = list(
      description = "Daily pioglitazone dose of the treatment arm.",
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = "Used only as an on-drug indicator: any dose > 0 switches the Weibull rate from the placebo to the drug value (Tran 2026 Discussion: no dose-response relationship was detected for the liver enzymes). 0 for placebo arms. At time 0 the prediction is rmin + wamp whatever the dose, so the value on a pre-treatment record does not matter.",
      source_name = "dose"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 540L,
    n_studies = 8L,
    age_range = "arm mean ages 40.4-55 years (Tran 2026 Table S2)",
    weight_range = "not reported",
    sex_female_pct = NA_real_,
    race_ethnicity = "not reported",
    disease_state = "biopsy-proven metabolic dysfunction-associated steatohepatitis (MASH, formerly NASH), with or without type 2 diabetes or prediabetes",
    dose_range = "pioglitazone 30 or 45 mg/day orally (Belfort 2006 and Cusi 2016 started at 30 mg/day and increased to 45 mg/day after 2 months), or placebo; follow-up to 36 months",
    regions = "not reported at study level",
    notes = "ALT model: all eight trials reported ALT (Tran 2026 Table S2); the arm-mean time courses are shown in Tran 2026 Figure S2. External validation used 27 MASLD patients from Jeonbuk National University Hospital on 15 or 30 mg/day (59 ALT values over up to 24 months; Tran 2026 Methods 2.3, Figure 4), where 40-50% of the observations fell inside the prediction interval."
  )

  ini({
    # Tran 2026 Methods 2.2.2 and Table 2 note (Weibull model):
    #   Effect = Baseline + A * exp(-lambda * time^gamma) + eps / sqrt(N)
    # The paper's 'Baseline' is the ASYMPTOTE the arm declines to: at time 0
    # exp(0) = 1, so the pre-treatment value is Baseline + A = 73.0 U/L.
    lrmin_alt <- log(33.9)
    label("Log asymptotic plateau of the arm-mean ALT (paper: Baseline; U/L)") # Tran 2026 Table 2, ALT Baseline = 33.9 (RSE 6.8%)
    lwamp_alt <- log(39.1)
    label("Log amplitude of the Weibull decline from time 0 to the plateau (paper: A; U/L)") # Tran 2026 Table 2, ALT A = 39.1 (RSE 14.1%)
    lwlambda_placebo <- log(0.0689)
    label("Log Weibull rate coefficient, placebo arms (paper: lambda1; month^-gamma)") # Tran 2026 Table 2, ALT lambda1 (Placebo group) = 0.0689 (RSE 50.5%)
    lwlambda_drug <- log(0.534)
    label("Log Weibull rate coefficient, pioglitazone arms (paper: lambda2; month^-gamma)") # Tran 2026 Table 2, ALT lambda2 (Drug group) = 0.534 (RSE 22.3%)
    lwgamma <- log(0.703)
    label("Log Weibull time exponent (paper: gamma; unitless)") # Tran 2026 Table 2, ALT gamma = 0.703 (RSE 23.6%)

    # Table 2 'Omega' rows. Observations are study-arm means, so these are
    # between-study random effects; taken as exponential (see vignette).
    eta_study_lrmin_alt ~ 0.0332 # Tran 2026 Table 2, ALT Omega Baseline = 0.0332 (RSE 34.9%)
    eta_study_lwamp_alt ~ 0.157 # Tran 2026 Table 2, ALT Omega A = 0.157 (RSE 39%)

    # Residual error, unweighted. The source variance is sigma^2 / N for an
    # arm of N patients (Methods 2.2.2); divide both SDs by sqrt(N) when
    # simulating arm means of a given size.
    propSd <- 0.292
    label("Proportional residual SD for N = 1 (fraction)") # Tran 2026 Table 2, ALT Proportion error = 0.292 (RSE 17%)
    addSd <- 10.5
    label("Additive residual SD for N = 1 (U/L)") # Tran 2026 Table 2, ALT Additive error = 10.5 (RSE 20.5%)
  })

  model({
    rmin <- exp(lrmin_alt + eta_study_lrmin_alt)
    wamp <- exp(lwamp_alt + eta_study_lwamp_alt)
    wgamma <- exp(lwgamma)

    # Pioglitazone at any dose switches the Weibull rate (Table 2,
    # lambda1 placebo group vs lambda2 drug group).
    wlambda <- exp(lwlambda_placebo)
    if (DOSE_PIO_MGD > 0) {
      wlambda <- exp(lwlambda_drug)
    }

    alt <- rmin + wamp * exp(-wlambda * t^wgamma)
    alt ~ prop(propSd) + add(addSd)
  })
}
