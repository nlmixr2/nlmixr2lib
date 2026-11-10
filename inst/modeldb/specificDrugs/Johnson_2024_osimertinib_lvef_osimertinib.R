Johnson_2024_osimertinib_lvef_osimertinib <- function() {
  description <- paste0(
    "Penalized (Firth) logistic-regression exposure-safety model for ",
    "left ventricular ejection fraction (LVEF) decrease events (a ",
    "post-baseline fall of at least 10 percentage points to a value ",
    "below 50%) in patients with EGFR-mutated non-small cell lung ",
    "cancer treated with osimertinib in the first-line or ",
    "second-or-later-line setting (Johnson 2024, final model of ",
    "Table 2). The probability of an event is expit(a + ",
    "b * log(AUC_OSIMERTINIB) + c * LVEF_BL), where AUC_OSIMERTINIB is ",
    "the individual steady-state osimertinib AUC in nmol*h/L (dose ",
    "divided by the individual apparent clearance of the companion ",
    "population PK model Johnson_2025_osimertinib) and LVEF_BL is the ",
    "baseline LVEF in percent. There is no PK layer and no ODE: ",
    "exposure is supplied as data. A sibling model driven by the ",
    "AZ5104 metabolite AUC is Johnson_2024_osimertinib_lvef_az5104."
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
    dosing = "n/a (no dose events; exposure enters as the AUC_OSIMERTINIB covariate column)",
    concentration = "prob_lvef_decrease (probability of an LVEF decrease event, 0-1)"
  )

  covariateData <- list(
    AUC_OSIMERTINIB = list(
      description = "Individual steady-state osimertinib AUC over the once-daily dosing interval, AUCss = dose / (CL/F)_i, supplied as data (this model has no PK layer).",
      units = "nmol*h/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Johnson 2024 Methods 'Estimation of plasma exposure': AUCss = dose/clearance from the individual apparent clearances of the population PK analysis (reference 23 = Johnson 2025, modellib('Johnson_2025_osimertinib')); for the safety analysis the FIRST dose a patient received was used. The unit is nM*h, the axis label of Johnson 2024 Figure 3A; Johnson_2025_osimertinib returns cl in L/h, so AUC_OSIMERTINIB = dose_mg / cl * 1e6 / 499.62. Enters on the natural-log scale, uncentred ('the model using log-transformed osimertinib or AZ5104 AUCss, with baseline LVEF as a covariate, was chosen as the final model', Johnson 2024 Results). Figure 3A quartile boundaries of the LVEF analysis set: 7012.1, 9584.5, 12092 and 17567.1 nM*h.",
      source_name = "AUCss (osimertinib)"
    ),
    LVEF_BL = list(
      description = "Left ventricular ejection fraction at baseline (before the first osimertinib dose).",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters the logit linearly and uncentred with slope -0.087 per percentage point (Johnson 2024 Table 2 and Table S10, 'Baseline LVEF (all patients)'): a higher baseline LVEF lowers the probability of a fall to below 50%, which is partly definitional because the event requires an end value below 50%. The source does not tabulate the baseline LVEF distribution. A restricted refit on patients with baseline LVEF 50 to < 70% (Table S10: intercept -5.34, slope 0.973, LVEF_BL -0.111) is a sensitivity analysis and is not encoded.",
      source_name = "Baseline LVEF"
    )
  )

  covariatesDataExcluded <- list(
    RACE_ASIAN = list(
      description = "Asian race indicator (1 = Asian, 0 = White; other races tested as a separate level).",
      units = "(binary)",
      type = "binary",
      notes = "Screened on the intercept of the osimertinib base model and not retained at the P < .001 threshold: slope -0.938 (95% CI -1.63, -0.263), p = 0.0067 (Johnson 2024 Table S10; Others vs White -1.66, p = 0.015). Age (p = 0.213), weight (p = 0.399), body mass index (p = 0.249), sex (p = 0.448), WHO performance status (p = 0.199), history of hypertension (p = 0.294), line of therapy (p = 0.396), nicotine use (p = 0.469), statin use (p = 0.052), and histories of diabetes, ischaemic heart disease, heart failure, coronary artery disease, hypothyroidism and hypoalbuminaemia (all p > 0.14) were screened the same way and also dropped."
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
    notes = "n_subjects is the sum of the per-dose denominators of the LVEF analysis in Johnson 2024 Results (4 + 909 + 21 + 2). Demographics are those of the full first-line plus second-or-later-line osimertinib population (Table S3, N = 1364), of which the LVEF set is the subset with LVEF assessments and a computable AUCss. Baseline LVEF is not tabulated. Comparator arms (chemotherapy 1/94, comparator EGFR-TKI 4/256 events) and the ADAURA adjuvant arms were not included in the fit."
  )

  ini({
    # Johnson 2024 Table 2, osimertinib final model (also Table S10 row
    # 'Baseline LVEF (all patients)'). Supplementary Methods:
    # log(p / (1 - p)) = alpha + beta_0 * C_i + beta_j * Z_j,i with C_i the
    # natural log of the AUCss in nM*h and Z the baseline LVEF in percent.
    logit_ref                 <- -6.36;  label("Logit of the LVEF-decrease probability at AUC_OSIMERTINIB = 1 nmol*h/L and LVEF_BL = 0% (unitless logit)")  # Johnson 2024 Table 2, osimertinib final model: Intercept -6.36 (95% CI -13.6, 1.01)
    e_auc_lvef_decrease       <- 0.921;  label("Log-odds of an LVEF decrease per e-fold increase in osimertinib AUCss (unitless logit)")                        # Johnson 2024 Table 2, osimertinib final model: Slope 0.921 (95% CI 0.209, 1.61)
    e_lvef_bl_lvef_decrease   <- -0.087; label("Log-odds of an LVEF decrease per percentage point of baseline LVEF (unitless logit per %)")                    # Johnson 2024 Table 2, osimertinib final model: Baseline LVEF -0.087 (95% CI -0.134, -0.042)

    # Binary landmark endpoint with a penalized Bernoulli likelihood and no
    # random effect. The placeholder additive residual exists only so that
    # rxode2 accepts an observation declaration.
    addSd_prob_lvef_decrease <- fixed(0.001); label("Placeholder additive residual SD on the LVEF-decrease probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # Natural-log AUCss and baseline LVEF, both uncentred, so logit_ref is
    # an extrapolation far outside the data and not a typical-patient
    # value.
    logit_lvef_decrease <- logit_ref + e_auc_lvef_decrease * log(AUC_OSIMERTINIB) +
      e_lvef_bl_lvef_decrease * LVEF_BL
    prob_lvef_decrease <- expit(logit_lvef_decrease)

    prob_lvef_decrease ~ add(addSd_prob_lvef_decrease)
  })
}
