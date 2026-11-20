Tang_2026_datopotamab_sae <- function() {
  description <- paste0(
    "Logistic-regression exposure-safety model for any serious ",
    "treatment-emergent adverse event (SAE) in adults with breast or ",
    "non-small cell lung cancer treated with datopotamab deruxtecan ",
    "(Dato-DXd) 0.27-10 mg/kg every 3 weeks (Tang 2026, n = 1081 pooled ",
    "from TROPION-PanTumor01, -Lung01, -Lung05 and -Breast01). The ",
    "probability is expit(1.362 + 1.383 * CAV - 2.944 * ALB / 39 - ",
    "0.861 * breast cancer + race shifts), where CAV is the individual ",
    "average plasma concentration of the released PAYLOAD DXd (ng/mL) ",
    "up to the end of the cycle of the event, and baseline albumin is ",
    "normalised by the cohort median of 39 g/L. Lung cancer and White ",
    "race are the reference levels. There is no PK layer and no ODE: ",
    "exposure is supplied as data and was derived in the source from ",
    "the companion population PK model packaged as ",
    "Hong_2025_datopotamab. No between-subject random effect and no ",
    "residual error are estimated. One of ten exposure-response models ",
    "in the Tang_2026_datopotamab_* family."
  )
  reference <- paste(
    "Tang Z, Lim K, Jiang Y, Denduluri N, Rokutanda N, Pan Y, Hong Y, Ren S,",
    "Vajjah P, Zhou D. Exposure-Response Analyses of Datopotamab Deruxtecan",
    "(Dato-DXd) in Patients with Hormone Receptor (HR) Positive, Human",
    "Epidermal Growth Factor Receptor 2 (HER2)-Negative Breast Cancer.",
    "J Clin Pharmacol. 2026;66(4):245-258. doi:10.1002/jcph.70188.",
    "Coefficients are transcribed from Supporting Information Table S4.",
    sep = " "
  )
  vignette <- "Tang_2026_datopotamab_exposure_response"
  units <- list(
    time = "n/a (static landmark exposure-response model; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "prob_sae (probability of any serious TEAE, 0-1)"
  )

  covariateData <- list(
    CAV = list(
      description = "Individual average plasma concentration of the released PAYLOAD DXd, cumulative AUC divided by the time from the first dose to the end of the cycle of the event (or, without an event, to the cycle of the last treatment or the data cut-off, whichever came first)",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "PAYLOAD (DXd), not the conjugate: ng/mL, whereas the companion Tang_2026_datopotamab_stomatitis / _ose / _ose_grade2 models carry an intact Dato-DXd Cavg in ug/mL under the same column name. Enters LINEARLY and UNCENTRED. Derived in the source from empirical Bayes estimates of the Hong 2025 population PK model (modellib('Hong_2025_datopotamab'); its Cc_dxd is in ug/mL, multiply by 1000). Tang 2026 Figure 2B shows the analysis range spanning roughly 0.3-1.3 ng/mL.",
      source_name = "DXd Cavg"
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters LINEARLY after division by the safety-set median of 39 g/L (Tang 2026 Methods: theta_j * x_j / x_median; Results: 'relative to the median albumin level (39 g/L)'). The published odds ratios fix the scaling: exp(-2.94387 * (46/39 - 1)) = 0.590 and exp(-2.94387 * (29/39 - 1)) = 2.13 against the printed 0.587 and 2.14 at the 95th and 5th percentiles. Table S2 range 19.0-75.4 g/L.",
      source_name = "Baseline albumin"
    ),
    TUMTP_BREAST = list(
      description = "Breast cancer indicator (1 = breast cancer, HR-positive/HER2-negative or triple-negative; 0 = lung cancer)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (lung cancer: non-small cell lung cancer, plus 1 large-cell neuroendocrine carcinoma)",
      notes = "Tang 2026 Table S4 term 'Breast Cancer'; Results: lower probability of SAE in breast versus lung cancer (odds ratio 0.425). Table S2: 437 of 1081 (40.4%) had breast cancer.",
      source_name = "Breast Cancer"
    ),
    RACE_BLACK = list(
      description = "Black or African American race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0; the reference level is White race (all RACE_* indicators 0)",
      notes = "One of four non-reference race levels in Table S4 (Black, Asian, Multiple or Other, Unknown or Unreported vs White). Not significant (p = 0.592); 19 of 1081 (1.8%).",
      source_name = "Race - Black"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0; the reference level is White race",
      notes = "Odds ratio 0.424 versus White (Tang 2026 Results). 432 of 1081 (40.0%).",
      source_name = "Race - Asian"
    ),
    RACE_OTHER = list(
      description = "Race recorded as 'Multiple or other'",
      units = "(binary)",
      type = "binary",
      reference_category = "0; the reference level is White race",
      notes = "The source level pools multiracial patients with those whose race was recorded as other; 79 of 1081 (7.3%). Table S2 also lists American Indian or Alaska Native (1) and Native Hawaiian or other Pacific Islander (1) patients and 8 with race missing; Tang 2026 does not say which model level these 10 patients were assigned to.",
      source_name = "Race - Multiple or Other"
    ),
    RACE_MISSING = list(
      description = "Race recorded as 'Unknown or not reported'",
      units = "(binary)",
      type = "binary",
      reference_category = "0; the reference level is White race",
      notes = "31 of 1081 (2.9%), all from TROPION-Breast01. Coefficient 0.039894, p = 0.924.",
      source_name = "Race - Unknown or Unreported"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1081L,
    n_studies = 4L,
    age_range = "26-86 years (median 60)",
    weight_range = "35.6-156 kg (median 64.2)",
    sex_female_pct = 67.2,
    race_ethnicity = "47.2% White, 40% Asian, 1.8% Black, 7.3% multiple or other, 2.9% unknown or not reported (Table S2)",
    disease_state = "advanced or metastatic non-small cell lung cancer (59.5%), HR-positive/HER2-negative breast cancer (36.4%) or triple-negative breast cancer (4.1%)",
    dose_range = "datopotamab deruxtecan 0.27-10 mg/kg IV every 3 weeks; 919 of 1081 (85%) at 6 mg/kg, 50 at 4 mg/kg and 82 at 8 mg/kg",
    regions = "USA 29.7%, Europe 27.8%, Japan 17.8%, mainland China 4.5%, rest of world 20.2%",
    notes = "Safety analysis set: patients with at least one Dato-DXd dose and at least one post-baseline PK record in TROPION-PanTumor01, TROPION-Lung01, TROPION-Lung05 and TROPION-Breast01 (data cut-off 29 April 2024). Observed SAE incidence 20%, 23% and 50% at 4, 6 and 8 mg/kg (Table 2)."
  )

  ini({
    # ==================================================================
    # Tang 2026 Supporting Information Table S4, 'Serious AE' rows, all
    # printed on the log-odds scale. Internal checks against the Results
    # odds ratios: exp(-0.85655) = 0.425 (Asian, printed 0.424) and
    # exp(-0.86092) = 0.423 (breast cancer, printed 0.425); the albumin
    # check is in the ALB covariateData note.
    # ==================================================================
    logit_ref <- 1.361896; label("Logit of the event probability at zero DXd exposure and zero albumin for a White patient with lung cancer (unitless logit)") # Tang 2026 Table S4, Serious AE (Intercept): 1.361896 [0.155; 2.57]
    e_cav_logit <- 1.382566; label("Log-odds of a serious TEAE per 1 ng/mL increase in DXd Cavg (unitless logit)") # Tang 2026 Table S4, Serious AE DXd Cavg: 1.382566 [0.974; 1.79]
    alb_ref <- fixed(39); label("Median baseline albumin used to normalise ALB (g/L)") # Tang 2026 Results: 'relative to the median albumin level (39 g/L)'
    e_alb_logit <- -2.94387; label("Log-odds shift per unit of median-normalised baseline albumin (unitless logit)") # Tang 2026 Table S4, Serious AE Baseline albumin: -2.94387 [-4.11; -1.78]
    e_tumtp_breast_logit <- -0.86092; label("Log-odds shift for breast cancer vs lung cancer (unitless logit)") # Tang 2026 Table S4, Serious AE Breast Cancer: -0.86092 [-1.21; -0.516]
    e_race_black_logit <- -0.2853; label("Log-odds shift for Black vs White race (unitless logit)") # Tang 2026 Table S4, Serious AE Race - Black: -0.2853 [-1.33; 0.757]
    e_race_asian_logit <- -0.85655; label("Log-odds shift for Asian vs White race (unitless logit)") # Tang 2026 Table S4, Serious AE Race - Asian: -0.85655 [-1.19; -0.521]
    e_race_other_logit <- -0.59128; label("Log-odds shift for multiple or other race vs White (unitless logit)") # Tang 2026 Table S4, Serious AE Race - Multiple or Other: -0.59128 [-1.16; -0.0191]
    e_race_missing_logit <- 0.039894; label("Log-odds shift for unknown or unreported race vs White (unitless logit)") # Tang 2026 Table S4, Serious AE Race - Unknown or Unreported: 0.039894 [-0.777; 0.857]

    addSd_prob_sae <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # CAV enters uncentred and albumin is divided by its median, so
    # logit_ref is the logit at zero exposure and zero albumin and is not
    # itself a patient; at median albumin the reference logit is
    # 1.361896 - 2.94387 = -1.582. The RACE_* indicators are mutually
    # exclusive; all zero selects White.
    logit_sae <- logit_ref +
      e_cav_logit * CAV +
      e_alb_logit * ALB / alb_ref +
      e_tumtp_breast_logit * TUMTP_BREAST +
      e_race_black_logit * RACE_BLACK +
      e_race_asian_logit * RACE_ASIAN +
      e_race_other_logit * RACE_OTHER +
      e_race_missing_logit * RACE_MISSING

    prob_sae <- expit(logit_sae)

    prob_sae ~ add(addSd_prob_sae)
  })
}
