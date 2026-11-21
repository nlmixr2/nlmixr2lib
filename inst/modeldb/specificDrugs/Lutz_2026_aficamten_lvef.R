Lutz_2026_aficamten_lvef <- function() {
  description <- paste(
    "Exposure-response (safety) model for left ventricular ejection fraction",
    "(LVEF) in adults with obstructive hypertrophic cardiomyopathy given",
    "aficamten (Lutz 2026). Log-linear direct-effect model: log LVEF is the",
    "log typical baseline, shifted by the observed baseline LVEF (power term",
    "centred at 75%) and by the reading source (clinical-site read vs central",
    "core-laboratory read), minus an exposure slope times Cavg,24 (CAV), the",
    "average aficamten plasma concentration over the 24 h before the",
    "echocardiogram, supplied from a population PK model such as",
    "Xu_2025_aficamten. Log-normal IIV on the slope and a log-normal residual",
    "error. The paired efficacy model is",
    "Lutz_2026_aficamten_lvot_gradient_valsalva."
  )
  reference <- paste(
    "Lutz JD, Wada R, Jacoby DL, Heitner SB, Kupfer S, German P. Towards a",
    "Quantitative Understanding of Aficamten Clinical Pharmacology:",
    "Pharmacokinetic-Cardiodynamic Modeling to Support Safety and Efficacy.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15(5):e70237. doi:10.1002/psp4.70237.",
    "PMID: 42101453. PMCID: PMC13155137."
  )
  vignette <- "Lutz_2026_aficamten"
  units <- list(
    time = "h",
    dosing = "mg (aficamten, given to the companion PK model that supplies CAV)",
    concentration = "ng/mL (CAV, the Cavg,24 exposure covariate)",
    lvef = "%"
  )

  covariateData <- list(
    CAV = list(
      description = "Average aficamten plasma concentration over the 24 h before the echocardiographic assessment (Cavg,24).",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Lutz 2026 Methods 'Analysis Endpoints and Exposure Metrics': the",
        "population PK model-estimated average concentration 24 h prior to the",
        "PD assessment, from the aficamten population PK model of Xu 2025",
        "(modellib('Xu_2025_aficamten')). Divided by 1000 inside the model so",
        "the slope is per ug/mL (Data S2 'CAVG1/1000'). 0 for placebo and at",
        "the pre-dose baseline visit. Time-varying: it follows dose titrations."
      ),
      source_name = "Cavg,24 (CAVG1 in Data S2)"
    ),
    LVEF_BL = list(
      description = "Observed left ventricular ejection fraction at baseline (week 0, pre-dose).",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power covariate on the typical LVEF, (BL / 75)^0.56, entered as",
        "(log(BASE) - log(MEDBASE)) * theta3 (Data S2). The paper does not",
        "report MEDBASE, the analysis-set median baseline; the maintainers",
        "use 75% (the rounded typical core-laboratory baseline, Table 1",
        "'InterceptTypical' 75.5%). The paper does not say which reading",
        "supplies BASE on site-read records; the maintainers supply the same",
        "(core-laboratory) baseline on both reading types, so the reading",
        "difference is carried by the READ_SITE effects alone."
      ),
      source_name = "BASE (Data S2); 'Base' (Table 1)"
    ),
    READ_SITE = list(
      description = "Echocardiogram reading source, 1 = read at the clinical site, 0 = read by the central core laboratory.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (central core-laboratory read)",
      notes = paste(
        "Lutz 2026 Table 1 'Type on log intercept' (theta4) and 'Type on log",
        "slope' (theta5); Data S2 codes it as TYPE == 'LOCAL'. Per",
        "observation: the trials read every echocardiogram both ways, and",
        "site reads guided dose titration while core-laboratory reads were",
        "the study endpoints. The type effect on the LVEF slope was not",
        "significant and was retained for consistency with the LVOT-G model",
        "(Results). The residual and between-subject variability are shared",
        "by both reading types (Table 1 note)."
      ),
      source_name = "TYPE == 'LOCAL' (Data S2); 'Type' (Table 1)"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Significant (p < 0.05) on the slope or intercept of the LVEF model but not retained, because it only minimally decreased residual variability and no subgroup simulations were planned (Lutz 2026 Results)."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Significant (p < 0.05) on the slope or intercept of the LVEF model but not retained (Lutz 2026 Results)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Significant (p < 0.05) on the slope or intercept of the LVEF model but not retained (Lutz 2026 Results)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 183L,
    n_studies = 2L,
    n_observations = 2876L,
    age_range = "Adults (not reported for the PKPD analysis set)",
    weight_range = "Not reported for the PKPD analysis set",
    disease_state = "Symptomatic obstructive hypertrophic cardiomyopathy",
    dose_range = "Aficamten 5 mg once daily titrated at weeks 2, 4 and 6 to a maximum of 15, 20 or 30 mg once daily by site-read echocardiography",
    regions = "Multinational",
    notes = paste(
      "All aficamten-treated participants with oHCM in phase 2 REDWOOD-HCM",
      "(41 participants, 10 weeks) and phase 3 SEQUOIA-HCM (NCT05186818; 142",
      "participants, 24 weeks) with PK parameters from the population PK",
      "analysis (Lutz 2026 Methods 'Study Populations'; Results). The 2876",
      "observations include both the clinical-site and the core-laboratory",
      "read of each echocardiogram, taken pre-dose on day 1 and 2-2.5 h",
      "after dosing at later visits. Race, region and study were also",
      "significant on the LVEF model and not retained. Fitted with nlme in",
      "R 4.2.2."
    )
  )

  ini({
    # Table 1 values are referenced to the core-laboratory read (Table 1 note).
    lrbase <- 4.32; label("Log typical core-laboratory LVEF at zero exposure for a baseline of 75% (log %)") # Table 1 'Log intercept' theta1 4.32 (SE < 0.01); back-transformed 75.5%
    lslope <- -1.16; label("Log typical exposure-response slope on Cavg,24, core-laboratory read (log 1/(ug/mL))") # Table 1 'Log slope' theta2 -1.16 (SE 0.08); back-transformed 0.31 1/(ug/mL)
    e_lvef_bl_rbase <- 0.56; label("Power exponent of observed baseline LVEF/75% on the typical LVEF (unitless)") # Table 1 'Base on intercept' theta3 0.56 (SE 0.02)
    e_read_site_rbase <- -0.08; label("Additive shift in log LVEF for a site read vs a core-laboratory read (log units)") # Table 1 'Type on log intercept' theta4 -0.08 (SE < 0.01); ratio 0.92
    e_read_site_slope <- -0.12; label("Additive shift in log slope for a site read vs a core-laboratory read (log units)") # Table 1 'Type on log slope' theta5 -0.12 (SE 0.05); ratio 0.88

    # Table 1 note: %CV = 100 x omega, so omega (the SD of eta) is 0.759.
    etalslope ~ 0.576081 # Table 1 'Between-subject CV on slope' 75.9% -> omega = 0.759, omega^2 = 0.576081

    # Within-subject variability, log-additive (Methods), %CV = 100 x sigma.
    expSd_lvef <- 0.076; label("Log-normal residual SD on LVEF (unitless)") # Table 1 'Within-subject CV' 7.6% -> sigma = 0.076
  })

  model({
    # Data S2: log(AVAL) = logInt + (log(BASE) - log(MEDBASE)) * CovBase +
    #   type_effect_int * (TYPE == 'LOCAL') -
    #   exp(logK + eta2 + type_effect_slope * (TYPE == 'LOCAL')) * CAVG1 / 1000
    slope <- exp(lslope + etalslope + e_read_site_slope * READ_SITE)
    log_lvef <- lrbase +
      e_lvef_bl_rbase * log(LVEF_BL / 75) +
      e_read_site_rbase * READ_SITE -
      slope * CAV / 1000

    lvef <- exp(log_lvef)

    lvef ~ lnorm(expSd_lvef)
  })
}
