Tran_2026_pioglitazone_ast_mbma <- function() {
  description <- "MBMA. Weibull time-course model of study-arm mean serum aspartate aminotransferase (AST) under pioglitazone or placebo in adults with metabolic dysfunction-associated steatohepatitis (MASH), fitted in NONMEM to arm-mean AST digitised from seven published trials. AST(t) = rmin + wamp * exp(-wlambda * t^wgamma): the arm starts at rmin + wamp (45.7 U/L) and declines towards the plateau rmin (29.6 U/L, which the paper calls 'Baseline'). Pioglitazone acts only on the Weibull rate, raising wlambda from 0.0352 (placebo) to 0.694 (any dose; no dose-response was found). Between-study random effects (exponential) on the plateau and the amplitude; proportional residual error (the additive component was fixed to 0), reported unweighted -- the source divides the residual by sqrt(N) for an arm of N patients, which downstream simulation of arm means must apply. Suitable for simulating study-arm mean AST trajectories, not individual patients. Companion models from the same paper: modellib('Tran_2026_pioglitazone_alt_mbma') and the four histology models modellib('Tran_2026_pioglitazone_steatosis_mbma'), modellib('Tran_2026_pioglitazone_inflammation_mbma'), modellib('Tran_2026_pioglitazone_ballooning_mbma') and modellib('Tran_2026_pioglitazone_fibrosis_mbma')."

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
    concentration = "U/L (the observation ast is the study-arm mean serum AST activity, not a drug concentration)"
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
    notes = "AST model: seven of the eight trials reported AST (all but Aithal 2008; Tran 2026 Table S2); the arm-mean time courses are shown in Tran 2026 Figure S2. External validation used 27 MASLD patients from Jeonbuk National University Hospital on 15 or 30 mg/day (59 AST values over up to 24 months; Tran 2026 Methods 2.3, Figure 5), where 40-50% of the observations fell inside the prediction interval."
  )

  ini({
    # Tran 2026 Methods 2.2.2 and Table 2 note (Weibull model):
    #   Effect = Baseline + A * exp(-lambda * time^gamma) + eps / sqrt(N)
    # The paper's 'Baseline' is the ASYMPTOTE the arm declines to: at time 0
    # exp(0) = 1, so the pre-treatment value is Baseline + A = 45.7 U/L.
    lrmin_ast <- log(29.6)
    label("Log asymptotic plateau of the arm-mean AST (paper: Baseline; U/L)") # Tran 2026 Table 2, AST Baseline = 29.6 (RSE 3.2%)
    lwamp_ast <- log(16.1)
    label("Log amplitude of the Weibull decline from time 0 to the plateau (paper: A; U/L)") # Tran 2026 Table 2, AST A = 16.1 (RSE 12.9%)
    lwlambda_placebo <- log(0.0352)
    label("Log Weibull rate coefficient, placebo arms (paper: lambda1; month^-gamma)") # Tran 2026 Table 2, AST lambda1 (Placebo group) = 0.0352 (RSE 28.4%)
    lwlambda_drug <- log(0.694)
    label("Log Weibull rate coefficient, pioglitazone arms (paper: lambda2; month^-gamma)") # Tran 2026 Table 2, AST lambda2 (Drug group) = 0.694 (RSE 16.1%)
    lwgamma <- log(0.591)
    label("Log Weibull time exponent (paper: gamma; unitless)") # Tran 2026 Table 2, AST gamma = 0.591 (RSE 21.2%)

    # Table 2 'Omega' rows. Observations are study-arm means, so these are
    # between-study random effects; taken as exponential (see vignette).
    eta_study_lrmin_ast ~ 0.0072 # Tran 2026 Table 2, AST Omega Baseline = 0.0072 (RSE 46.3%)
    eta_study_lwamp_ast ~ 0.0782 # Tran 2026 Table 2, AST Omega A = 0.0782 (RSE 38.4%)

    # Residual error, unweighted. The source variance is sigma^2 / N for an
    # arm of N patients (Methods 2.2.2); divide the SD by sqrt(N) when
    # simulating arm means of a given size. Table 2 fixes the additive
    # component to 0 ('Additive error 0 FIX'), so only the proportional
    # term is encoded.
    propSd <- 0.352
    label("Proportional residual SD for N = 1 (fraction)") # Tran 2026 Table 2, AST Proportion error = 0.352 (RSE 9.8%); Additive error = 0 FIX
  })

  model({
    rmin <- exp(lrmin_ast + eta_study_lrmin_ast)
    wamp <- exp(lwamp_ast + eta_study_lwamp_ast)
    wgamma <- exp(lwgamma)

    # Pioglitazone at any dose switches the Weibull rate (Table 2,
    # lambda1 placebo group vs lambda2 drug group).
    wlambda <- exp(lwlambda_placebo)
    if (DOSE_PIO_MGD > 0) {
      wlambda <- exp(lwlambda_drug)
    }

    ast <- rmin + wamp * exp(-wlambda * t^wgamma)
    ast ~ prop(propSd)
  })
}
