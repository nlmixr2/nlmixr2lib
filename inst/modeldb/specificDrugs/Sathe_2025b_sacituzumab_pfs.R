Sathe_2025b_sacituzumab_pfs <- function() {
  description <- paste0(
    "Cox proportional-hazards exposure-efficacy model for ",
    "PROGRESSION-FREE SURVIVAL in adults with metastatic triple-negative ",
    "breast cancer treated with sacituzumab govitecan 8 or 10 mg/kg on ",
    "days 1 and 8 of 21-day cycles (Sathe 2025b, n = 275, pooled ",
    "IMMU-132-01 and ASCENT). The log relative hazard is ",
    "-0.00730 * CAV + 0.495 (prior checkpoint inhibitor without platinum) ",
    "+ 0.736 (prior platinum without checkpoint inhibitor) + 1.01 (both), ",
    "relative to patients whose prior therapy included neither (the ",
    "'Others' category), where CAV is the individual average serum ",
    "concentration of TOTAL ANTIBODY (tAB) from the first dose to the ",
    "event or censoring, in ug/mL. Outputs: lhr, hr and prob_pfs_12mo, ",
    "the probability of being progression-free at 12 months, computed ",
    "with a 12-month baseline cumulative hazard BACK-SOLVED from the ",
    "eight model-predicted probabilities of Table S11 (all eight agree on ",
    "3.42-3.44); a Cox fit leaves the baseline hazard nonparametric and ",
    "the paper does not print it. There is no PK layer and no ODE; ",
    "exposure derives from the companion population PK model ",
    "Sathe_2024_sacituzumab. One of six exposure-response models in the ",
    "Sathe_2025b_sacituzumab_* family."
  )
  reference <- paste(
    "Sathe AG, Diderichsen PM, Fauchet F, Phan SC, Girish S, Othman AA.",
    "Exposure-Response Analyses of Sacituzumab Govitecan Efficacy and",
    "Safety in Patients With Metastatic Triple-Negative Breast Cancer.",
    "Clin Pharmacol Ther. 2025;117(2):570-578. doi:10.1002/cpt.3495.",
    "PMCID PMC11739744. Coefficients from Supplementary Table S10; the",
    "12-month baseline cumulative hazard is back-solved from Table S11.",
    sep = " "
  )
  vignette <- "Sathe_2025b_sacituzumab_exposure_response"

  units <- list(
    time = "months (prob_pfs_12mo is a single landmark at 12 months; lhr and hr are time-invariant under proportional hazards)",
    dosing = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "prob_pfs_12mo (probability of being alive and progression-free at 12 months, 0-1); also hr (relative hazard, unitless) and lhr (its log)"
  )

  covariateData <- list(
    CAV = list(
      description = "Individual average serum concentration of sacituzumab govitecan TOTAL ANTIBODY (tAB, hRS7 conjugated or unconjugated to SN-38) from the first dose to the progression or death event or to censoring. Supplied as data: this model has no PK layer.",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Sathe 2025b Eq. 1 with t_event from Table S1 (time of event, or 21 days after the last SG dose if no event was reported). CAVG_tAB was the most significant of nine exposure metrics for PFS (Table S8, -2LL drop 25.9). Reproduce the column from the total-antibody output (Cc_tab) of modellib('Sathe_2024_sacituzumab'). Enters LINEARLY and UNCENTRED. Table 2 PFS exposure quartiles: medians 95.1, 130, 158 and 205 ug/mL (overall range 50.9-450); Table S11 footnote: mean 122 ug/mL at 8 mg/kg and 152 ug/mL at 10 mg/kg.",
      source_name = "CAVG_tAB"
    ),
    PRIOR_PLATIN = list(
      description = "Prior platinum-based chemotherapy (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (platinum-naive)",
      notes = "Sathe 2025b encodes prior treatment as one four-level factor -- 'Others' (reference), 'CPI', 'platinum based' and 'CPI + platinum based' (Table S3 and Table S10). The factor is decomposed here into PRIOR_PLATIN and PRIOR_CHECKPOINT_INHIBITOR, and the three printed coefficients are applied to the three non-reference cells exactly as printed (they are not additive: 0.495 + 0.736 = 1.231, not 1.01). Analysis-set counts: platinum based only 133 (48.0%), CPI only 25 (9.0%), CPI + platinum based 53 (19.1%), Others 66 (23.8%).",
      source_name = "prior treatment"
    ),
    PRIOR_CHECKPOINT_INHIBITOR = list(
      description = "Prior immune checkpoint inhibitor therapy (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (checkpoint-inhibitor-naive)",
      notes = "The 'CPI' component of the Sathe 2025b four-level prior-treatment factor; see PRIOR_PLATIN for the decomposition. All CPI-pretreated patients came from ASCENT (Table S3). The paper does not list which checkpoint inhibitors were received.",
      source_name = "prior treatment"
    )
  )

  covariatesDataExcluded <- list(
    LDH = list(
      description = "Baseline serum lactate dehydrogenase",
      units = "IU/L",
      type = "continuous",
      notes = "Screened for PFS (Sathe 2025b Methods, Covariate selection) but not retained in the final PFS model; it was retained for OS (see Sathe_2025b_sacituzumab_os). The same applies to baseline body weight, age, race, region, number of prior lines of therapy, ECOG performance status, baseline Trop-2 H score and UGT1A1 genotype. No point estimate is reported for any of them."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 275L,
    n_studies = 2L,
    age_range = "27-82 years (median 54)",
    weight_range = "37.2-132 kg (median 66.7)",
    sex_female_pct = 99.3,
    race_ethnicity = c(White = 80.5, Black = 11.2, Asian = 4.0, Other = 4.3),
    disease_state = "metastatic triple-negative breast cancer, previously treated",
    dose_range = "sacituzumab govitecan 8, 10 or 12 mg/kg starting dose IV on days 1 and 8 of 21-day cycles (10 mg/kg in about 92%)",
    regions = "North America 67.9%, Europe 32.1%",
    notes = "The 277-patient exposure-response set less two ASCENT patients without PFS data. 204 PFS events (74.2%; Table S4). PFS was by blinded independent central review in ASCENT and by investigator assessment in IMMU-132-01, where central review was missing for most patients. The proportional-hazards assumption was formally violated for CAVG_tAB, but Schoenfeld and Martingale residuals showed limited bias and the model was accepted (Supplementary Section 1)."
  )

  ini({
    # Sathe 2025b Supplementary Table S10, PFS rows (log-hazard scale;
    # positive = higher risk of progression or death).
    e_cav_pfs_haz <- -0.00730; label("Log hazard ratio for progression or death per 1 ug/mL increase in CAVG_tAB (per ug/mL)") # Sathe 2025b Table S10, PFS Slope for CAVG_tAB (per ug/mL): -0.00730 [-0.0104 to -0.00423], p < 0.001
    e_cpi_pfs_haz <- 0.495; label("Log hazard ratio for progression or death, prior checkpoint inhibitor without platinum versus Others (log scale)") # Sathe 2025b Table S10, PFS Coefficient for prior treatment CPI vs. others: 0.495 [-0.0898 to 1.08], p = 0.0971
    e_platin_pfs_haz <- 0.736; label("Log hazard ratio for progression or death, prior platinum without checkpoint inhibitor versus Others (log scale)") # Sathe 2025b Table S10, PFS Coefficient for prior treatment platinum based vs. others: 0.736 [0.357-1.12], p < 0.001
    e_cpi_platin_pfs_haz <- 1.01; label("Log hazard ratio for progression or death, prior checkpoint inhibitor and platinum versus Others (log scale)") # Sathe 2025b Table S10, PFS Coefficient for prior treatment CPI + platinum based vs. others: 1.01 [0.559-1.45], p < 0.001

    # Baseline (all covariates zero) cumulative hazard at 12 months, NOT
    # printed. Back-solved from the eight Table S11 predictions as
    # -log(S12) / exp(lhr) at the footnoted mean CAVG_tAB (122 ug/mL at
    # 8 mg/kg, 152 ug/mL at 10 mg/kg) for each prior-treatment category:
    # 3.433, 3.435 (CPI), 3.433, 3.439 (platinum), 3.420, 3.426
    # (CPI + platinum), 3.437, 3.437 (Others); mean 3.43.
    ch0_pfs_12mo <- 3.43; label("Baseline cumulative hazard of progression or death at 12 months for the all-zero covariate reference (unitless)") # back-solved from Sathe 2025b Table S11 (8 rows, 3.420-3.439); no printed value exists
  })

  model({
    # Prior-treatment categories: Others (both indicators 0, reference),
    # CPI only, platinum only, and both. Each non-reference cell carries
    # its own printed coefficient.
    lhr <- e_cav_pfs_haz * CAV +
      e_cpi_pfs_haz * PRIOR_CHECKPOINT_INHIBITOR * (1 - PRIOR_PLATIN) +
      e_platin_pfs_haz * PRIOR_PLATIN * (1 - PRIOR_CHECKPOINT_INHIBITOR) +
      e_cpi_platin_pfs_haz * PRIOR_CHECKPOINT_INHIBITOR * PRIOR_PLATIN
    hr <- exp(lhr)

    # Probability of being progression-free at the 12-month landmark.
    # Valid at 12 months only: the baseline hazard at other times is not
    # recoverable from the paper.
    prob_pfs_12mo <- exp(-ch0_pfs_12mo * hr)
  })
}
