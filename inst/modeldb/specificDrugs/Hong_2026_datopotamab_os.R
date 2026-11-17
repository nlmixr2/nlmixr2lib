Hong_2026_datopotamab_os <- function() {
  description <- paste(
    "Parametric time-to-event exposure-efficacy model for OVERALL SURVIVAL",
    "in adults with previously treated advanced non-small cell lung cancer",
    "receiving datopotamab deruxtecan (Dato-DXd) 0.27-10 mg/kg IV every 3",
    "weeks (Hong 2026, n = 644 pooled from TROPION-PanTumor01, -Lung05 and",
    "-Lung01). Gompertz hazard h(t) = 0.412 * exp(0.434 * t) * exp(-0.0216",
    "* AUC1 + covariates), t in years, where AUC1 is the individual Cycle-1",
    "Dato-DXd AUC in mg*h/mL (supplied as AUC_DATODXD in ug*day/mL and",
    "converted inside the model). Retained covariates on the log hazard:",
    "baseline albumin (centred at 38 g/L), baseline tumour size (centred at",
    "67 mm), ECOG performance status > 0, squamous histology, liver",
    "metastasis and female sex. There is no PK layer: exposure is supplied",
    "as data and was derived in the source from the companion population",
    "PK model packaged as Hong_2025_datopotamab. No between-subject random",
    "effect and no residual error. Outputs `hazard` (1/year), `cumhaz` and",
    "`sur`. One of six exposure-response models in the",
    "Hong_2026_datopotamab_* family."
  )
  reference <- paste(
    "Hong Y, Mc Laughlin AM, Pan Y, Hietala SF, Peigne S, Chargualaf M,",
    "Tang Z, Zhou D, Tajima N, Abutarif M, Garimella T. Exposure-Response",
    "Analyses of Datopotamab Deruxtecan (Dato-DXd), a TROP2-Directed",
    "Antibody-Drug Conjugate, in Advanced Non-small Cell Lung Cancer.",
    "Clin Pharmacol Ther. 2026. doi:10.1002/cpt.70477. PMID 42733211.",
    "Estimates from Table 2; hazard and covariate structure from",
    "Supplementary Equation S1. Individual exposures derive from the",
    "population PK model of Hong 2025 (doi:10.1002/psp4.70118); see",
    "modellib('Hong_2025_datopotamab').",
    sep = " "
  )
  vignette <- "Hong_2026_datopotamab_exposure_response"
  units <- list(
    time = "year",
    dosing = "n/a (no dose events; exposure enters as the AUC_DATODXD data column)",
    concentration = "probability (the output `sur` is the probability of being alive, not a drug concentration)"
  )

  compartmentData <- list(
    cumhaz = list(
      analyte = "cumulative hazard of death from any cause",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    AUC_DATODXD = list(
      description = "Model-predicted Cycle-1 (day 0 to day 21) area under the plasma concentration-time curve of intact datopotamab deruxtecan",
      units = "ug*day/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Hong 2026 reports this exposure (AUC1) in mg/mL*h; the model",
        "converts with AUC_DATODXD * 24 / 1000, so 633 ug*day/mL = 15.2",
        "mg*h/mL, the Figure 2 reference value (median AUC1 in the 6 mg/kg",
        "group). Enters LINEARLY and UNCENTRED on the log hazard",
        "(Supplementary Equation S1, h = h0 * exp(f(Exposure) + ...)), so",
        "the baseline hazard is the value at zero exposure. Figure 2",
        "percentiles in the 6 mg/kg group: 7.38 / 12.1 / 18.7 / 25.3 mg*h/mL",
        "(5th / 25th / 75th / 95th); medians 10.2 and 21.4 mg*h/mL at 4 and",
        "8 mg/kg. Derived in the source from empirical Bayes estimates of the",
        "Hong 2025 population PK model, not observed."
      ),
      source_name = "Dato-DXd AUC1"
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed baseline value. Linear on the log hazard, centred at the analysis-set median 38 g/L (Table S3 overall median 38.0; Figure 2 reference patient 38 g/L).",
      source_name = "Albumin"
    ),
    TUMSZ = list(
      description = "Baseline tumour size (sum of diameters of target lesions)",
      units = "mm",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed baseline value. Linear on the log hazard, centred at the analysis-set median 67 mm (Table S3 overall median 67.0; Figure 2 reference patient 67 mm).",
      source_name = "Tumor size"
    ),
    ECOG_GE1 = list(
      description = "Baseline ECOG performance status > 0 (1 = ECOG 1 or 2, 0 = ECOG 0)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Hong 2026 Table 2 'ECOG > 0 on the baseline hazard'. ECOG 2 was recorded in a single patient (Table S3).",
      source_name = "ECOG"
    ),
    TUMTP_NSCLC_SQUAM = list(
      description = "Squamous non-small cell lung cancer histology (1 = squamous, 0 = nonsquamous or other)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Hong 2026 Table 2 'Squamous histology on the baseline hazard'; the reference patient has nonsquamous histology (Figure 2).",
      source_name = "Histology"
    ),
    LMET = list(
      description = "Liver metastasis at baseline (1 = present, 0 = absent)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Hong 2026 Table 2 'Liver metastasis on the baseline hazard'.",
      source_name = "Liver metastases"
    ),
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Hong 2026 Table 2 'Female sex on the baseline hazard'; males are the reference.",
      source_name = "Sex"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 644L,
    n_studies = 3L,
    age_range = "26-84 years (median 63)",
    weight_range = "37.0-156 kg (median 65.0)",
    sex_female_pct = 45.8,
    race_ethnicity = "44% White, 42% Asian, 1.9% Black, 11% other, 0.9% missing (Table S3)",
    disease_state = "Previously treated advanced or metastatic non-small cell lung cancer (528 nonsquamous, 98 squamous, 18 other histology; one large-cell neuroendocrine carcinoma); 36% with actionable genomic alterations",
    dose_range = "Datopotamab deruxtecan 0.27-10 mg/kg IV every 3 weeks; 75% at 6 mg/kg",
    regions = "USA 35%, Europe 26%, Japan 22%, rest of world 16%, China 0.8%",
    notes = "Exposure-efficacy dataset of TROPION-PanTumor01 (NSCLC only, n = 210), TROPION-Lung05 (n = 137) and TROPION-Lung01 (Dato-DXd arm, n = 297). Baseline albumin median 38.0 g/L (19.0-50.0), tumour size median 67.0 mm (13.0-314), ECOG > 0 in 71%, liver metastasis in 23% (Table S3)."
  )

  ini({
    # Hong 2026 Table 2 (run 20106, OFV 5818.07). Supplementary Equation S1:
    # h(t) = lambda * exp(gamma * t) * exp(f(Exposure) + sum(beta_i * X_i)).
    # The printed Table 1 / Figure 2 hazard ratios are an internal check on
    # the log-hazard reading: exp(0.569) = 1.77 (ECOG, printed 1.80),
    # exp(-0.0612 * (28 - 38)) = 1.84 (albumin, printed 1.84),
    # exp(0.00477 * (173 - 67)) = 1.66 (tumour size, printed 1.65) and
    # exp(-0.0216 * (7.38 - 15.2)) = 1.18 (AUC1, printed 1.17).
    lhaz_base <- log(0.412); label("Gompertz scale: baseline hazard of death at zero exposure and reference covariates (1/year)") # Table 2 'Baseline hazard' = 0.412 Year-1 (RSE 20.3%)
    gomp_haz <- 0.434; label("Gompertz shape: coefficient of time on the log hazard (1/year)") # Table 2 'Shape parameter' = 0.434 Year-1 (RSE 18.9%)
    e_auc_datodxd_haz <- -0.0216; label("Log-hazard slope per unit of Cycle-1 Dato-DXd AUC (1/(mg*h/mL))") # Table 2 'Slope for linear effect of Dato-DXd AUC1 on the baseline hazard' = -0.0216 (RSE 42.5%)
    e_alb_haz <- -0.0612; label("Log-hazard slope per unit of baseline albumin above 38 g/L (1/(g/L))") # Table 2 'Albumin at baseline on the baseline hazard' = -0.0612 (RSE 17.5%)
    e_ecog_ge1_haz <- 0.569; label("Log-hazard shift for ECOG performance status > 0 (unitless)") # Table 2 'ECOG > 0 on the baseline hazard' = 0.569 (RSE 20.5%)
    e_tumtp_nsclc_squam_haz <- 0.393; label("Log-hazard shift for squamous histology (unitless)") # Table 2 'Squamous histology on the baseline hazard' = 0.393 (RSE 37.0%)
    e_lmet_haz <- 0.411; label("Log-hazard shift for liver metastasis (unitless)") # Table 2 'Liver metastasis on the baseline hazard' = 0.411 (RSE 27.8%)
    e_sexf_haz <- -0.279; label("Log-hazard shift for female sex (unitless)") # Table 2 'Female sex on the baseline hazard' = -0.279 (RSE 37.1%)
    e_tumsz_haz <- 0.00477; label("Log-hazard slope per unit of baseline tumour size above 67 mm (1/mm)") # Table 2 'Tumor size at baseline on the baseline hazard' = 0.00477 mm-1 (RSE 22.1%)
  })
  model({
    # Cycle-1 AUC in the paper's unit: ug*day/mL * 24 h/day / 1000 ug/mg.
    auc1_mghml <- AUC_DATODXD * 24 / 1000

    hazard <- exp(lhaz_base + gomp_haz * t +
      e_auc_datodxd_haz * auc1_mghml +
      e_alb_haz * (ALB - 38) +
      e_tumsz_haz * (TUMSZ - 67) +
      e_ecog_ge1_haz * ECOG_GE1 +
      e_tumtp_nsclc_squam_haz * TUMTP_NSCLC_SQUAM +
      e_lmet_haz * LMET +
      e_sexf_haz * SEXF)

    d/dt(cumhaz) <- hazard
    cumhaz(0) <- 0
    sur <- exp(-cumhaz)
  })
}
