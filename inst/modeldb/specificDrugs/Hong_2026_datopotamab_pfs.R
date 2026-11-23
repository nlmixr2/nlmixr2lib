Hong_2026_datopotamab_pfs <- function() {
  description <- paste(
    "Parametric time-to-event exposure-efficacy model for PROGRESSION-FREE",
    "SURVIVAL in adults with previously treated advanced non-small cell",
    "lung cancer receiving datopotamab deruxtecan (Dato-DXd) 0.27-10 mg/kg",
    "IV every 3 weeks (Hong 2026, n = 644 pooled from TROPION-PanTumor01,",
    "-Lung05 and -Lung01). Exponential (constant) baseline hazard 5.85 per",
    "year with an Emax effect of the TIME-VARYING cycle-wise Dato-DXd AUC",
    "over the dosing interval on the log hazard, -3.45 * AUCtau / (17.6 +",
    "AUCtau) with AUCtau in mg*h/mL (supplied as AUC_DATODXD_TAU in",
    "ug*day/mL and converted inside the model). Retained covariates on the",
    "log hazard: age (centred at 63 years), baseline albumin (centred at",
    "38 g/L), squamous histology, last prior line of immuno-oncology",
    "therapy, liver metastasis and female sex. There is no PK layer:",
    "exposure is supplied as data and was derived in the source from the",
    "companion population PK model packaged as Hong_2025_datopotamab. No",
    "between-subject random effect and no residual error. Outputs `hazard`",
    "(1/year), `cumhaz` and `sur`. One of six exposure-response models in",
    "the Hong_2026_datopotamab_* family."
  )
  reference <- paste(
    "Hong Y, Mc Laughlin AM, Pan Y, Hietala SF, Peigne S, Chargualaf M,",
    "Tang Z, Zhou D, Tajima N, Abutarif M, Garimella T. Exposure-Response",
    "Analyses of Datopotamab Deruxtecan (Dato-DXd), a TROP2-Directed",
    "Antibody-Drug Conjugate, in Advanced Non-small Cell Lung Cancer.",
    "Clin Pharmacol Ther. 2026. doi:10.1002/cpt.70477. PMID 42733211.",
    "Estimates from Supplementary Table S4; hazard and covariate structure",
    "from Supplementary Equation S1. Individual exposures derive from the",
    "population PK model of Hong 2025 (doi:10.1002/psp4.70118); see",
    "modellib('Hong_2025_datopotamab').",
    sep = " "
  )
  vignette <- "Hong_2026_datopotamab_exposure_response"
  units <- list(
    time = "year",
    dosing = "n/a (no dose events; exposure enters as the time-varying AUC_DATODXD_TAU data column)",
    concentration = "probability (the output `sur` is the probability of being alive and progression-free, not a drug concentration)"
  )

  compartmentData <- list(
    cumhaz = list(
      analyte = "cumulative hazard of radiographic progression or death",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    AUC_DATODXD_TAU = list(
      description = "Model-predicted intact datopotamab deruxtecan plasma AUC over the current dosing interval (cycle-wise, time-varying)",
      units = "ug*day/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "TIME-VARYING: one value per treatment cycle, reflecting the",
        "individual dosing history (interruptions, reductions) and exposure",
        "changes during treatment (Hong 2026 Methods). Supply it as a",
        "time-varying column, carried forward within each cycle. Hong 2026",
        "reports it in mg/mL*h; the model converts with AUC_DATODXD_TAU *",
        "24 / 1000, so 579 ug*day/mL = 13.9 mg*h/mL, the Figure S2 reference",
        "value (median of the individual mean AUCtau in the 6 mg/kg group).",
        "Figure S2 percentiles: 2.52 / 9.89 / 18.3 / 26.9 mg*h/mL (5th /",
        "25th / 75th / 95th); medians 10.2 and 18.0 mg*h/mL at 4 and 8",
        "mg/kg. Enters through an Emax function that is zero at zero",
        "exposure, so the baseline hazard is the value without drug."
      ),
      source_name = "Dato-DXd AUCtau"
    ),
    AGE = list(
      description = "Age at baseline",
      units = "year",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear on the log hazard, centred at the analysis-set median 63 years (Table S3 overall median 63.0; Figure S2 reference patient 63 years).",
      source_name = "Age"
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed baseline value. Linear on the log hazard, centred at 38 g/L (Figure S2 reference patient). SIGN: Table S4 prints +0.0266, but the Results text, Table 1 (albumin 28 vs 38 g/L, HR 1.29; 45 vs 38 g/L, HR 0.836) and Figure S2 (HR 1.29 / 1.10 / 0.926 / 0.836 at 28.0 / 34.4 / 41.0 / 45.0 g/L) all show the hazard FALLING as albumin rises, with magnitudes that exp(-0.0266 * (ALB - 38)) reproduces. The model applies that reading; see the vignette Assumptions and deviations.",
      source_name = "Albumin"
    ),
    TUMTP_NSCLC_SQUAM = list(
      description = "Squamous non-small cell lung cancer histology (1 = squamous, 0 = nonsquamous or other)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Hong 2026 Table S4 'Squamous histology on the baseline hazard'; the reference patient has nonsquamous histology (Figure S2).",
      source_name = "Histology"
    ),
    LAST_PRIOR_CHECKPOINT_INHIBITOR = list(
      description = "Last prior line of therapy contained an immune checkpoint inhibitor (1 = yes, 0 = no; the source calls it immuno-oncology, LPIO)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Hong 2026 'last prior line therapy containing immuno-oncology (LPIO)'; the Discussion describes it as an immune checkpoint inhibitor (PD-1/PD-L1 blockade) in the last prior line. 54% of the analysis set (Table S3).",
      source_name = "LPIO"
    ),
    LMET = list(
      description = "Liver metastasis at baseline (1 = present, 0 = absent)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Hong 2026 Table S4 'Liver metastasis on the base hazard'.",
      source_name = "Liver metastases"
    ),
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Hong 2026 Table S4 'Female sex on the baseline hazard'; males are the reference.",
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
    notes = "Exposure-efficacy dataset of TROPION-PanTumor01 (NSCLC only, n = 210), TROPION-Lung05 (n = 137) and TROPION-Lung01 (Dato-DXd arm, n = 297). Progression assessed by blinded independent central review (RECIST v1.1). Last prior line containing immuno-oncology therapy in 54%, liver metastasis in 23% (Table S3)."
  )

  ini({
    # Hong 2026 Supplementary Table S4 (run 20211, OFV 5310.1).
    # Supplementary Equation S1 with an exponential baseline:
    # h = lambda * exp(Emax * AUCtau / (EC50 + AUCtau) + sum(beta_i * X_i)).
    # Table S4 prints the EC50 unit as '(mg/mL*h)-1'; an EC50 has the
    # unit of the exposure it divides, mg/mL*h.
    lhaz_base <- log(5.85); label("Baseline hazard of progression or death without drug at reference covariates (1/year)") # Table S4 'Baseline hazard' = 5.85 year-1 (RSE 12.8%)
    emax_haz <- -3.45; label("Maximal change in the log hazard from Dato-DXd AUCtau (unitless)") # Table S4 'Emax for AUCtau on the base hazard' = -3.45 (RSE 20.4%)
    lec50 <- log(17.6); label("Dato-DXd AUCtau at half the maximal log-hazard change (mg*h/mL)") # Table S4 'EC50 for AUCtau on the base hazard' = 17.6 (RSE 41.9%)
    e_age_haz <- -0.0156; label("Log-hazard slope per year of age above 63 years (1/year)") # Table S4 'Age on the base hazard' = -0.0156 year-1 (RSE 30.9%)
    e_alb_haz <- -0.0266; label("Log-hazard slope per unit of baseline albumin above 38 g/L (1/(g/L))") # Table S4 'Albumin at baseline on the base hazard' printed 0.0266 (RSE 37.1%); sign applied as negative per Table 1, Figure S2 and Results (see covariateData ALB notes)
    e_tumtp_nsclc_squam_haz <- 0.598; label("Log-hazard shift for squamous histology (unitless)") # Table S4 'Squamous histology on the baseline hazard' = 0.598 (RSE 22.6%)
    e_last_prior_checkpoint_inhibitor_haz <- -0.303; label("Log-hazard shift for an immune checkpoint inhibitor in the last prior line (unitless)") # Table S4 'Last prior line of IO on the base hazard' = -0.303 (RSE 34.3%)
    e_lmet_haz <- 0.537; label("Log-hazard shift for liver metastasis (unitless)") # Table S4 'Liver metastasis on the base hazard' = 0.537 (RSE 23.4%)
    e_sexf_haz <- -0.312; label("Log-hazard shift for female sex (unitless)") # Table S4 'Female sex on the baseline hazard' = -0.312 (RSE 34.0%)
  })
  model({
    # Cycle-wise AUC in the paper's unit: ug*day/mL * 24 h/day / 1000 ug/mg.
    auctau_mghml <- AUC_DATODXD_TAU * 24 / 1000
    ec50 <- exp(lec50)
    drug_lhaz <- emax_haz * auctau_mghml / (ec50 + auctau_mghml)

    hazard <- exp(lhaz_base + drug_lhaz +
      e_age_haz * (AGE - 63) +
      e_alb_haz * (ALB - 38) +
      e_tumtp_nsclc_squam_haz * TUMTP_NSCLC_SQUAM +
      e_last_prior_checkpoint_inhibitor_haz * LAST_PRIOR_CHECKPOINT_INHIBITOR +
      e_lmet_haz * LMET +
      e_sexf_haz * SEXF)

    d/dt(cumhaz) <- hazard
    cumhaz(0) <- 0
    sur <- exp(-cumhaz)
  })
}
