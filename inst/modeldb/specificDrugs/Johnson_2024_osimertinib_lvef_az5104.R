Johnson_2024_osimertinib_lvef_az5104 <- function() {
  description <- paste0(
    "Penalized (Firth) logistic-regression exposure-safety model for ",
    "left ventricular ejection fraction (LVEF) decrease events (a ",
    "post-baseline fall of at least 10 percentage points to a value ",
    "below 50%) in patients with EGFR-mutated non-small cell lung ",
    "cancer treated with osimertinib in the first-line or ",
    "second-or-later-line setting, driven by exposure to the active ",
    "metabolite AZ5104 (Johnson 2024, final model of Table 2). The ",
    "probability of an event is expit(a + b * log(AUC_AZ5104) + ",
    "c * LVEF_BL), where AUC_AZ5104 is the individual steady-state ",
    "AZ5104 AUC in nmol*h/L derived from the companion population PK ",
    "model Johnson_2025_osimertinib and LVEF_BL is the baseline LVEF ",
    "in percent. There is no PK layer and no ODE: exposure is supplied ",
    "as data. Osimertinib and AZ5104 exposures are highly correlated ",
    "and were never entered together; the parent-driven alternative is ",
    "Johnson_2024_osimertinib_lvef_osimertinib."
  )
  reference <- paste(
    "Johnson M, Lin YW, Schmidt H, Sunnaker M, Van Maanen E, Huang X,",
    "Rukazenkov Y, Tomkinson H, Vishwanathan K. Exposure-response",
    "modelling of osimertinib in patients with non-small cell lung",
    "cancer. Br J Clin Pharmacol. 2024;90(12):3263-3276.",
    "doi:10.1111/bcp.16199. PMCID PMC11602960.",
    sep = " "
  )
  vignette <- "Johnson_2024_osimertinib_exposure_response"

  units <- list(
    time = "n/a (static landmark exposure-response model; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the AUC_AZ5104 covariate column)",
    concentration = "prob_lvef_decrease (probability of an LVEF decrease event, 0-1)"
  )

  covariateData <- list(
    AUC_AZ5104 = list(
      description = "Individual steady-state AUC of the osimertinib active metabolite AZ5104 over the once-daily dosing interval, derived from the individual apparent AZ5104 clearance of the population PK model and supplied as data (this model has no PK layer).",
      units = "nmol*h/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Johnson 2024 Methods 'Estimation of plasma exposure': AUCss of osimertinib and/or AZ5104 derived from the individual apparent clearances of the population PK analysis (reference 23 = Johnson 2025, modellib('Johnson_2025_osimertinib')) and the first dose received. At steady state the molar AZ5104 AUC is fm * dose / (CLmetabolite/F)_i with the fixed conversion fraction fm = 0.25, i.e. AUC_AZ5104 = 0.25 * dose_mg / cl_az5104 * 1e6 / 499.62 using the cl_az5104 (L/h) of Johnson_2025_osimertinib. Enters on the natural-log scale, uncentred ('log-transformed osimertinib or AZ5104 AUCss', Johnson 2024 Results).",
      source_name = "AUCss (AZ5104)"
    ),
    LVEF_BL = list(
      description = "Left ventricular ejection fraction at baseline (before the first osimertinib dose).",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters the logit linearly and uncentred with slope -0.087 per percentage point (Johnson 2024 Table 2 and Table S11, 'Baseline LVEF (all patients)'). The source does not tabulate the baseline LVEF distribution. A restricted refit on patients with baseline LVEF 50 to < 70% (Table S11: intercept -2.76, slope 0.907, LVEF_BL -0.11) is a sensitivity analysis and is not encoded.",
      source_name = "Baseline LVEF"
    )
  )

  covariatesDataExcluded <- list(
    RACE_ASIAN = list(
      description = "Asian race indicator (1 = Asian, 0 = White; other races tested as a separate level).",
      units = "(binary)",
      type = "binary",
      notes = "Screened on the intercept of the AZ5104 base model and not retained at the P < .001 threshold: slope -0.882 (95% CI -1.57, -0.205), p = 0.010 (Johnson 2024 Table S11; Others vs White -1.61, p = 0.020). Age (p = 0.224), weight (p = 0.551), body mass index (p = 0.372), sex (p = 0.537), WHO performance status (p = 0.260), history of hypertension (p = 0.279), line of therapy (p = 0.433), nicotine use (p = 0.534), statin use (p = 0.044), and histories of diabetes, ischaemic heart disease, heart failure, coronary artery disease, hypothyroidism and hypoalbuminaemia (all p > 0.13) were screened the same way and also dropped."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 936L,
    n_studies = 5L,
    n_events = "36 LVEF decrease events: 0/4 at 40 mg, 35/909 at 80 mg, 1/21 at 160 mg and 0/2 at 240 mg (Johnson 2024 Results and Figure 3B)",
    age_range = "25-91 years (median 62; Johnson 2024 Table S3, first-line plus second-or-later-line osimertinib population)",
    weight_range = "29-122 kg (median 60-62; Johnson 2024 Table S3)",
    sex_female_pct = 64.7,
    race_ethnicity = c(
      Japanese = 18.6,
      Chinese = 20.3,
      `Asian non-Chinese non-Japanese` = 23.7,
      White = 27.5,
      Other = 9.9
    ),
    disease_state = "EGFR mutation-positive advanced non-small cell lung cancer (first-line, or T790M-positive after prior EGFR-TKI therapy)",
    dose_range = "osimertinib 40-240 mg once daily (first dose received)",
    regions = "Multinational; 62% Asian",
    notes = "Same analysis set as Johnson_2024_osimertinib_lvef_osimertinib. n_subjects is the sum of the per-dose denominators of the LVEF analysis in Johnson 2024 Results (4 + 909 + 21 + 2). Demographics are those of the full first-line plus second-or-later-line osimertinib population (Table S3, N = 1364). Baseline LVEF is not tabulated."
  )

  ini({
    # Johnson 2024 Table 2, AZ5104 final model (also Table S11 row
    # 'Baseline LVEF (all patients)'); C_i is the natural log of the AZ5104
    # AUCss in nM*h and the covariate is the baseline LVEF in percent.
    logit_ref                 <- -4.38;  label("Logit of the LVEF-decrease probability at AUC_AZ5104 = 1 nmol*h/L and LVEF_BL = 0% (unitless logit)")  # Johnson 2024 Table 2, AZ5104 final model: Intercept -4.38 (95% CI -9.88, 1.07)
    e_auc_lvef_decrease       <- 0.927;  label("Log-odds of an LVEF decrease per e-fold increase in AZ5104 AUCss (unitless logit)")                        # Johnson 2024 Table 2, AZ5104 final model: Slope 0.927 (95% CI 0.292, 1.56)
    e_lvef_bl_lvef_decrease   <- -0.087; label("Log-odds of an LVEF decrease per percentage point of baseline LVEF (unitless logit per %)")               # Johnson 2024 Table 2, AZ5104 final model: Baseline LVEF -0.087 (95% CI -0.134, -0.041)

    # Binary landmark endpoint with a penalized Bernoulli likelihood and no
    # random effect. The placeholder additive residual exists only so that
    # rxode2 accepts an observation declaration.
    addSd_prob_lvef_decrease <- fixed(0.001); label("Placeholder additive residual SD on the LVEF-decrease probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # Natural-log AUCss and baseline LVEF, both uncentred, so logit_ref is
    # an extrapolation far outside the data and not a typical-patient
    # value.
    logit_lvef_decrease <- logit_ref + e_auc_lvef_decrease * log(AUC_AZ5104) +
      e_lvef_bl_lvef_decrease * LVEF_BL
    prob_lvef_decrease <- expit(logit_lvef_decrease)

    prob_lvef_decrease ~ add(addSd_prob_lvef_decrease)
  })
}
