Sheng_2026_sugemalimab_pfs_bicr <- function() {
  description <- "Parametric time-to-event exposure-response model for progression-free survival assessed by BLINDED INDEPENDENT CENTRAL REVIEW (BICR) in adults with stage IV non-small-cell lung cancer receiving sugemalimab 1200 mg IV every 3 weeks with platinum-based chemotherapy (Sheng 2026; GEMSTONE-302 sugemalimab arm, n = 320). Gompertz hazard with a piecewise-constant baseline (0-1, 1-4 and >4 months), h(t) = BL_k * exp(gomp * t) * HR_C^(CTROUGH / 10) * HR_PDL1^(1 - PDL1_TUM_POS), where CTROUGH is the model-predicted Cycle-1 trough serum sugemalimab concentration (ug/mL) and the hazard ratio per 10 ug/mL is 0.876; tumour PD-L1 < 1% raises the hazard by 1.485. Same structure as the investigator-assessed model, refitted to the BICR event times. No PK layer: CTROUGH is a time-fixed data column, the individual empirical-Bayes prediction from the companion population PK model Sheng_2026_sugemalimab. Outputs `hazard` (events per month), `cumhaz` and `sur`. Sister models: Sheng_2026_sugemalimab_pfs_inv, Sheng_2026_sugemalimab_os."
  reference <- paste(
    "Sheng Y, Yue Z, Xu F, Wang K, Wang J, Shi Q.",
    "Model-Informed Dosing Regimen of Sugemalimab for European Patients",
    "With Non-Small Cell Lung Cancer: Bridging From Asian Clinical Data.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15(3):e70220.",
    "doi:10.1002/psp4.70220.",
    "Estimates from Supplementary Table S3; hazard structure from the",
    "Supplementary Data S2 NONMEM control stream, whose header states that",
    "the BICR fit reuses the investigator control stream with PFSICR as",
    "TIME and PFSICRST as DV.",
    sep = " "
  )
  vignette <- "Sheng_2026_sugemalimab"
  units <- list(
    time = "month",
    dosing = "n/a (no drug-dosing events; the exposure driver is supplied as the CTROUGH data covariate, in ug/mL)",
    concentration = "probability (the model output `sur` is the progression-free survival probability, not a drug concentration)"
  )

  compartmentData <- list(
    cumhaz = list(
      analyte = "cumulative hazard of progression or death (BICR-assessed)",
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
      notes = "Required input, time-fixed per subject. Serum sugemalimab trough of cycle 1 (end of the first 3-week interval of 1200 mg IV every 3 weeks; the paper does not state total vs unbound) (source column CMINC1, 'Cmin,C1' in Table S3), the individual empirical-Bayes prediction from the companion population PK model modellib('Sheng_2026_sugemalimab'). Enters LINEARLY, UNCENTRED and RESCALED per 10 ug/mL: exp(e_ctrough_haz * CTROUGH / 10).",
      source_name = "CMINC1"
    ),
    PDL1_TUM_POS = list(
      description = "Tumour PD-L1 expression >= 1% (or not available)",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (PD-L1 >= 1% or missing)",
      notes = "The hazard is multiplied by 1.485 when PDL1_TUM_POS = 0 (tumour proportion score < 1%). Table S3 footnote: 'PDL1 of missing grouped with >1%', so a subject with no PD-L1 result must be coded PDL1_TUM_POS = 1.",
      source_name = "PDL1"
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
    notes = "Sugemalimab arm of the phase III GEMSTONE-302 study (study 302 in Table S1). GEMSTONE-302 observed BICR-assessed PFS 0.410 at 12 months and 0.230 at 24 months (Sheng 2026 Table 3 note)."
  )

  ini({
    # Sheng 2026 Table S3, 'PFS (BICR)' column. BL values are hazards per
    # month; hazard ratios carried as log hazard ratios.
    lhaz_base_p1 <- log(0.043); label("Baseline hazard of progression or death, 0-1 month (1/month)") # Table S3 PFS (BICR) Baseline hazard (0-1 Months) = 0.043 (95% CI 0.015-0.125)
    lhaz_base_p2 <- log(0.16); label("Baseline hazard of progression or death, 1-4 months (1/month)") # Table S3 PFS (BICR) Baseline hazard (1-4 Months) = 0.16 (95% CI 0.083-0.309)
    lhaz_base_p3 <- log(0.347); label("Baseline hazard of progression or death, after 4 months (1/month)") # Table S3 PFS (BICR) Baseline hazard (4+ Months) = 0.347 (95% CI 0.173-0.698)
    gomp_haz <- -0.0768; label("Gompertz shape: coefficient of time on the log hazard (1/month)") # Table S3 PFS (BICR) Shape parameter = -0.0768 (95% CI -0.107 to -0.047)
    e_ctrough_haz <- log(0.876); label("Log hazard ratio per 10 ug/mL of Cycle-1 trough sugemalimab (HR 0.876)") # Table S3 PFS (BICR) Linear effect of Cmin,C1 per 10 ug/mL, HR = 0.876 (95% CI 0.799-0.961)
    e_pdl1neg_haz <- log(1.485); label("Log hazard ratio for tumour PD-L1 < 1% versus >= 1% or missing (HR 1.485)") # Table S3 PFS (BICR) PDL1 (<1% vs >1%) HR = 1.485 (95% CI 1.111-1.986)
  })
  model({
    lhaz_base <- lhaz_base_p3
    if (t <= 1) {
      lhaz_base <- lhaz_base_p1
    } else if (t <= 4) {
      lhaz_base <- lhaz_base_p2
    }

    hazard <- exp(lhaz_base + gomp_haz * t +
      e_ctrough_haz * CTROUGH / 10 +
      e_pdl1neg_haz * (1 - PDL1_TUM_POS))

    d/dt(cumhaz) <- hazard
    cumhaz(0) <- 0
    sur <- exp(-cumhaz)
  })
}
