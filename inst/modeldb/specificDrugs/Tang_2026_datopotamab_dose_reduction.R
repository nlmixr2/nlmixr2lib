Tang_2026_datopotamab_dose_reduction <- function() {
  description <- paste0(
    "Logistic-regression exposure-safety model for a treatment-emergent ",
    "adverse event leading to dose reduction in adults with breast or ",
    "non-small cell lung cancer treated with datopotamab deruxtecan ",
    "(Dato-DXd) 0.27-10 mg/kg every 3 weeks (Tang 2026, n = 1081 pooled ",
    "from TROPION-PanTumor01, -Lung01, -Lung05 and -Breast01). The ",
    "probability is expit(-2.706 + 0.002585 * AUC_DATODXD - 0.404 * ",
    "TUMSZ / 60), where AUC_DATODXD is the individual Cycle-1 intact ",
    "Dato-DXd AUC in ug*day/mL and baseline tumor size is normalised by ",
    "the safety-set median of 60 mm; patients with larger tumours were ",
    "less likely to have a dose reduction. There is no PK layer and no ",
    "ODE: exposure is supplied as data and was derived in the source ",
    "from the companion population PK model packaged as ",
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
    dosing = "n/a (no dose events; exposure enters as the AUC_DATODXD covariate column)",
    concentration = "prob_dose_reduction (probability of a TEAE leading to dose reduction, 0-1)"
  )

  covariateData <- list(
    AUC_DATODXD = list(
      description = "Model-predicted Cycle-1 (day 0 to day 21) area under the plasma concentration-time curve of intact datopotamab deruxtecan",
      units = "ug*day/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters LINEARLY and UNCENTRED. Derived in the source from empirical Bayes estimates of the Hong 2025 population PK model (modellib('Hong_2025_datopotamab')).",
      source_name = "Dato-DXd AUC1"
    ),
    TUMSZ = list(
      description = "Baseline tumor size (sum of target-lesion diameters)",
      units = "mm",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters LINEARLY after division by the safety-set median of 60 mm (Tang 2026 Methods: theta_j * x_j / x_median; Results: 'relative to the median tumor size baseline (60 mm)'). The published odds ratios fix the scaling: exp(-0.40425 * (166/60 - 1)) = 0.490 and exp(-0.40425 * (18/60 - 1)) = 1.327 against the printed 0.489 and 1.33 at the 95th and 5th percentiles. Table S2 range 10-341 mm. Source column name TUMSBL.",
      source_name = "TUMSBL"
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
    notes = "Safety analysis set (data cut-off 29 April 2024). Observed incidence 2% (1/50), 21% (193/918) and 28% (23/82) at 4, 6 and 8 mg/kg (Table 2)."
  )

  ini({
    # Tang 2026 Supporting Information Table S4, 'Dose Reduction' rows,
    # printed on the log-odds scale. The tumor-size check against the
    # printed odds ratios is in the TUMSZ covariateData note.
    logit_ref <- -2.70582; label("Logit of the event probability at zero Cycle-1 Dato-DXd AUC and zero tumor size (unitless logit)") # Tang 2026 Table S4, Dose Reduction (Intercept): -2.70582 [-3.32; -2.09]
    e_auc_datodxd_logit <- 0.002585; label("Log-odds of a TEAE leading to dose reduction per 1 ug*day/mL increase in Cycle-1 Dato-DXd AUC (unitless logit)") # Tang 2026 Table S4, Dose Reduction Dato-DXd AUC1: 0.002585 [0.00185; 0.00332]
    tumsz_ref <- fixed(60); label("Median baseline tumor size used to normalise TUMSZ (mm)") # Tang 2026 Results: 'relative to the median tumor size baseline (60 mm)'
    e_tumsz_logit <- -0.40425; label("Log-odds shift per unit of median-normalised baseline tumor size (unitless logit)") # Tang 2026 Table S4, Dose Reduction TUMSBL: -0.40425 [-0.635; -0.174]

    addSd_prob_dose_reduction <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # At the median tumor size the logit offset is -0.40425, so the
    # median-patient intercept is -3.110.
    logit_dose_reduction <- logit_ref +
      e_auc_datodxd_logit * AUC_DATODXD +
      e_tumsz_logit * TUMSZ / tumsz_ref
    prob_dose_reduction <- expit(logit_dose_reduction)

    prob_dose_reduction ~ add(addSd_prob_dose_reduction)
  })
}
