Sheng_2026_sugemalimab_os <- function() {
  description <- "Parametric time-to-event exposure-response model for OVERALL SURVIVAL in adults with stage IV non-small-cell lung cancer receiving sugemalimab 1200 mg IV every 3 weeks with platinum-based chemotherapy (Sheng 2026; GEMSTONE-302 sugemalimab arm, n = 320). Constant (exponential) hazard, h = BL * HR_C^(CTROUGH / 10), with BL 0.137 deaths per month and a hazard ratio of 0.798 per 10 ug/mL of the model-predicted Cycle-1 trough serum sugemalimab concentration CTROUGH; no other covariate. No PK layer: CTROUGH is a time-fixed data column, the individual empirical-Bayes prediction from the companion population PK model Sheng_2026_sugemalimab. Outputs `hazard` (deaths per month), `cumhaz`, `sur` (survival probability) and `medtte` (closed-form median survival, months). Sister models: Sheng_2026_sugemalimab_pfs_inv, Sheng_2026_sugemalimab_pfs_bicr."
  reference <- paste(
    "Sheng Y, Yue Z, Xu F, Wang K, Wang J, Shi Q.",
    "Model-Informed Dosing Regimen of Sugemalimab for European Patients",
    "With Non-Small Cell Lung Cancer: Bridging From Asian Clinical Data.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15(3):e70220.",
    "doi:10.1002/psp4.70220.",
    "Estimates from Supplementary Table S3; hazard structure from the",
    "Supplementary Data S2 NONMEM control stream (OS).",
    sep = " "
  )
  vignette <- "Sheng_2026_sugemalimab"
  units <- list(
    time = "month",
    dosing = "n/a (no drug-dosing events; the exposure driver is supplied as the CTROUGH data covariate, in ug/mL)",
    concentration = "probability (the model output `sur` is the overall survival probability, not a drug concentration)"
  )

  compartmentData <- list(
    cumhaz = list(
      analyte = "cumulative hazard of death",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    CTROUGH = list(
      description = "Model-predicted Cycle-1 trough serum sugemalimab concentration",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Required input, time-fixed per subject. TOTAL serum sugemalimab at the end of the first 21-day cycle of 1200 mg IV every 3 weeks (source column CMINC1, 'Cmin,C1' in Table S3), the individual empirical-Bayes prediction from the companion population PK model modellib('Sheng_2026_sugemalimab'). Enters LINEARLY, UNCENTRED and RESCALED per 10 ug/mL: exp(e_ctrough_haz * CTROUGH / 10).",
      source_name = "CMINC1"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 320L,
    n_studies = 1L,
    disease_state = "Stage IV (metastatic) non-small-cell lung cancer, first-line treatment with platinum-based chemotherapy (Sheng 2026 Methods; GEMSTONE-302, Zhou 2022).",
    dose_range = "Sugemalimab 1200 mg IV every 3 weeks with platinum-based chemotherapy (single dose level).",
    regions = "Asia; the pivotal study enrolled only Asian patients (Sheng 2026 Discussion).",
    weight_range = "41-96 kg, mean 62 kg (Sheng 2026 Discussion).",
    notes = "Sugemalimab arm of the phase III GEMSTONE-302 study (CS1001-302). GEMSTONE-302 observed OS 0.716 at 12 months and 0.522 at 24 months (Sheng 2026 Table 3 note). The authors caution that the OS model may be biased because 45 placebo patients crossed over to sugemalimab after progression."
  )

  ini({
    lhaz_base <- log(0.137); label("Baseline hazard of death at zero exposure (1/month)") # Table S3 OS Baseline hazard (BL) = 0.137 (95% CI 0.071-0.261)
    e_ctrough_haz <- log(0.798); label("Log hazard ratio per 10 ug/mL of Cycle-1 trough sugemalimab (HR 0.798)") # Table S3 OS Linear effect of Cmin,C1 per 10 ug/mL, HR = 0.798 (95% CI 0.728-0.875)
  })
  model({
    # Control stream $DES: DADT(1) = EXP(LBASE)*EXP(CMINC1/10*B1) for T > 0.
    hazard <- exp(lhaz_base + e_ctrough_haz * CTROUGH / 10)

    d/dt(cumhaz) <- hazard
    cumhaz(0) <- 0
    sur <- exp(-cumhaz)

    # Constant hazard: the median survival time solves exp(-hazard * t) = 0.5.
    medtte <- log(2) / hazard
  })
}
