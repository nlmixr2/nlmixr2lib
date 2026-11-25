Lutz_2026_aficamten_lvot_gradient_valsalva <- function() {
  description <- paste(
    "Exposure-response (efficacy) model for the post-Valsalva left",
    "ventricular outflow tract gradient (LVOT-G) in adults with obstructive",
    "hypertrophic cardiomyopathy given aficamten (Lutz 2026). Log-linear",
    "direct-effect model: log LVOT-G is the log typical baseline, shifted by",
    "the observed baseline gradient (power term centred at 90 mmHg) and by",
    "the reading source (clinical-site read vs central core-laboratory read),",
    "minus an exposure slope times Cavg,24 (CAV), the average aficamten",
    "plasma concentration over the 24 h before the echocardiogram, supplied",
    "from a population PK model such as Xu_2025_aficamten. Log-normal IIV on",
    "the slope and a log-normal residual error. The paired safety model is",
    "Lutz_2026_aficamten_lvef."
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
    lvot_gradient_valsalva = "mmHg"
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
    LVOT_GRADIENT_VALSALVA_BL = list(
      description = "Observed post-Valsalva LVOT gradient at baseline (week 0, pre-dose).",
      units = "mmHg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power covariate on the typical gradient, (BL / 90)^0.63, entered as",
        "(log(BASE) - log(MEDBASE)) * theta3 (Data S2). The paper does not",
        "report MEDBASE, the analysis-set median baseline; the maintainers",
        "use 90 mmHg (the rounded typical core-laboratory baseline, Table 1",
        "'InterceptTypical' 89.4 mmHg). The paper does not say which reading",
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
        "the study endpoints. The residual and between-subject variability",
        "are shared by both reading types (Table 1 note)."
      ),
      source_name = "TYPE == 'LOCAL' (Data S2); 'Type' (Table 1)"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened in the stepwise covariate search (Lutz 2026 Methods 'Covariates') and not retained in the LVOT-G model."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened in the stepwise covariate search and not retained in the LVOT-G model."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened in the stepwise covariate search and not retained in the LVOT-G model."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 183L,
    n_studies = 2L,
    n_observations = 2880L,
    age_range = "Adults (not reported for the PKPD analysis set)",
    weight_range = "Not reported for the PKPD analysis set",
    disease_state = "Symptomatic obstructive hypertrophic cardiomyopathy",
    dose_range = "Aficamten 5 mg once daily titrated at weeks 2, 4 and 6 to a maximum of 15, 20 or 30 mg once daily by site-read echocardiography",
    regions = "Multinational",
    notes = paste(
      "All aficamten-treated participants with oHCM in phase 2 REDWOOD-HCM",
      "(41 participants, 10 weeks) and phase 3 SEQUOIA-HCM (NCT05186818; 142",
      "participants, 24 weeks) with PK parameters from the population PK",
      "analysis (Lutz 2026 Methods 'Study Populations'; Results). The 2880",
      "observations include both the clinical-site and the core-laboratory",
      "read of each post-Valsalva echocardiogram, taken pre-dose on day 1 and",
      "2-2.5 h after dosing at later visits. Fitted with nlme in R 4.2.2."
    )
  )

  ini({
    # Table 1 values are referenced to the core-laboratory read (Table 1 note).
    lrbase <- 4.49; label("Log typical core-laboratory LVOT-G at zero exposure for a baseline of 90 mmHg (log mmHg)") # Table 1 'Log intercept' theta1 4.49 (SE 0.02); back-transformed 89.4 mmHg
    lslope <- 1.44; label("Log typical exposure-response slope on Cavg,24, core-laboratory read (log 1/(ug/mL))") # Table 1 'Log slope' theta2 1.44 (SE 0.05); back-transformed 4.23 1/(ug/mL)
    e_lvot_gradient_valsalva_bl_rbase <- 0.63; label("Power exponent of observed baseline LVOT-G/90 mmHg on the typical LVOT-G (unitless)") # Table 1 'Base on intercept' theta3 0.63 (SE 0.03)
    e_read_site_rbase <- 0.12; label("Additive shift in log LVOT-G for a site read vs a core-laboratory read (log units)") # Table 1 'Type on log intercept' theta4 0.12 (SE 0.03); ratio 1.13
    e_read_site_slope <- 0.07; label("Additive shift in log slope for a site read vs a core-laboratory read (log units)") # Table 1 'Type on log slope' theta5 0.07 (SE 0.03); ratio 1.07

    # Table 1 note: %CV = 100 x omega, so omega (the SD of eta) is 0.596.
    etalslope ~ 0.355216 # Table 1 'Between-subject CV on slope' 59.6% -> omega = 0.596, omega^2 = 0.355216

    # Within-subject variability, log-additive (Methods), %CV = 100 x sigma.
    expSd_lvot_gradient_valsalva <- 0.496; label("Log-normal residual SD on LVOT-G (unitless)") # Table 1 'Within-subject CV' 49.6% -> sigma = 0.496
  })

  model({
    # Data S2: log(AVAL) = logInt + (log(BASE) - log(MEDBASE)) * CovBase +
    #   type_effect_int * (TYPE == 'LOCAL') -
    #   exp(logK + eta2 + type_effect_slope * (TYPE == 'LOCAL')) * CAVG1 / 1000
    slope <- exp(lslope + etalslope + e_read_site_slope * READ_SITE)
    log_lvot <- lrbase +
      e_lvot_gradient_valsalva_bl_rbase * log(LVOT_GRADIENT_VALSALVA_BL / 90) +
      e_read_site_rbase * READ_SITE -
      slope * CAV / 1000

    lvot_gradient_valsalva <- exp(log_lvot)

    lvot_gradient_valsalva ~ lnorm(expSd_lvot_gradient_valsalva)
  })
}
