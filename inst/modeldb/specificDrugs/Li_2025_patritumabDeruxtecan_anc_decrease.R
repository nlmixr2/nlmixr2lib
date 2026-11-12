Li_2025_patritumabDeruxtecan_anc_decrease <- function() {
  description <- paste(
    "Bayesian logistic-regression exposure-safety model for grade >= 3",
    "neutrophil count decrease in adults with breast cancer, non-small-cell",
    "lung cancer or colorectal cancer given patritumab deruxtecan (HER3-DXd)",
    "in the J101, U102, U201 and U202 studies (Li 2025, n = 698;",
    "1.6-8 mg/kg IV Q3W, up-titration or an alternative Q2W/Q3W regimen).",
    "logit(P) is log-linear in the Cycle 1 DXd (unconjugated payload) Cmax,",
    "linear in baseline body weight, baseline neutrophil count and baseline",
    "hemoglobin, with additive shifts for country of origin (Japan and",
    "non-Japanese Asian country vs non-Asian country). Li 2025 prints no",
    "coefficient table; every value is back-solved from the typical-value",
    "forest plot (Figure S3B). Covariate-by-exposure interaction terms of the",
    "source model are not recoverable and are omitted. There is no PK layer:",
    "exposure is a data column from the companion population PK model (see",
    "Lee_2023_patritumab). One of nine Li_2025_patritumabDeruxtecan",
    "exposure-response models."
  )
  reference <- paste(
    "Li L, Lee M, Joshi R, Xu Y, Garcia R, Baron K, Anderson E, Waterhouse T,",
    "Sternberg D, Shrestha P, Vaddady P, Abutarif M, Garimella T.",
    "Exposure-Response Relationships in Patients with Non-Small-Cell Lung",
    "Cancer and Other Solid Tumors Treated with Patritumab Deruxtecan",
    "(HER3-DXd). Clin Pharmacol Ther. 2025;118(2):408-417.",
    "doi:10.1002/cpt.3674."
  )
  vignette <- "Li_2025_patritumabDeruxtecan_exposure_response"
  units <- list(
    time = "n/a (static landmark exposure-safety model; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CMAX covariate column)",
    concentration = "prob_anc_decrease (probability of a grade >= 3 neutrophil count decrease, 0-1)"
  )

  covariateData <- list(
    CMAX = list(
      description = "Cycle 1 maximum plasma concentration of DXd, the unconjugated topoisomerase I inhibitor payload, per subject.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Unconjugated payload (DXd), not the conjugate. Cycle 1 (first-dose) landmark; median onset of the event was day 16 (Li 2025 Table S4). Empirical-Bayes predictions from the population PK model 'previously developed' (uncited in Li 2025; the same group's model is Lee_2023_patritumab). Centred at 23.8 ng/mL inside model(), the reference exposure implied by Figure S3B; Table S5 reports a median of 24.2 ng/mL at 5.6 mg/kg Q3W. Figure S3B 5th/25th/75th/95th percentiles: 8.0 / 15.5 / 33.9 / 58.5 ng/mL.",
      source_name = "Cycle 1 DXd Cmax (ng/mL)"
    ),
    WT = list(
      description = "Baseline body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Influential but not significant (Li 2025 Table 2). Centred at the 60 kg safety-population median inside model() (Results; the Figure S3B rows imply 60.4 kg). Figure S3B 5th/25th/75th/95th percentiles: 42 / 52 / 72 / 93 kg.",
      source_name = "Baseline weight (kg)"
    ),
    NEUT = list(
      description = "Baseline absolute neutrophil count.",
      units = "10^9/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline only, in 10^9 cells/L (not cells/mm^3). Significant covariate (Li 2025 Table 2). Centred at 4.1 x 10^9/L inside model(), the reference value implied by Figure S3B; Table S1 mean 4.88 (SD 3.36). Figure S3B 5th/25th/75th/95th percentiles: 2.1 / 3.0 / 5.8 / 9.7.",
      source_name = "Baseline neutrophils (10^9/L)"
    ),
    HGB = list(
      description = "Baseline hemoglobin concentration.",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline only. Influential but not significant (Li 2025 Table 2). Figure S3B and Table S1 label the unit g/dL, but the printed values (5th-95th percentile 95-144; mean 118, SD 15.1) are in g/L. Supply g/L. Centred at 118 g/L inside model(), the reference value implied by Figure S3B and the Table S1 mean.",
      source_name = "Baseline hemoglobin (g/dL as labelled; values in g/L)"
    ),
    REGION_JAPAN = list(
      description = "Country of origin Japan indicator (1 = Japan, 0 = otherwise).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian country of origin when REGION_ASIA_NONJAPAN is also 0).",
      notes = "Significant covariate. Figure S3B 'Country of origin: Japanese' 49.1% vs reference 18.4%; observed rates 47.6% in Japan, 21.3% in other Asian countries and 17.7% in non-Asian countries (Results). 34.4% of the safety population (Table S1).",
      source_name = "Country of origin: Japanese"
    ),
    REGION_ASIA_NONJAPAN = list(
      description = "Country of origin in Asia other than Japan indicator (1 = non-Japanese Asian country, 0 = otherwise).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian country of origin when REGION_JAPAN is also 0).",
      notes = "Figure S3B 'Country of origin: Non-Japanese Asian' 20.6% vs reference 18.4%. 95 of 715 patients (13.3%), all from U102 and U201 (Table S1).",
      source_name = "Country of origin: Non-Japanese Asian"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 698L,
    n_studies = 4L,
    age_range = "29-83 years (median 61)",
    weight_range = "32.4-137 kg (median 60.0)",
    sex_female_pct = 69.2,
    race_ethnicity = c(Asian = 54.7, White = 37.5, Other = 5.7, Black = 2.1),
    disease_state = "HER3-expressing advanced breast cancer (J101), metastatic or unresectable NSCLC (U102, U201; EGFR-mutated in 446 of 493) and advanced colorectal cancer (U202).",
    dose_range = "Patritumab deruxtecan 1.6-8 mg/kg IV Q3W, up-titration (3.2 / 4.8 / 6.4 mg/kg) or 4.2 mg/kg Q2W for Cycles 1-3 then 6.4 mg/kg Q3W.",
    regions = "Global; country of origin Japan 34.4%, other Asian 13.3%, non-Asian 52.3%.",
    biomarkers = "Baseline neutrophils mean 4.88 (SD 3.36) x 10^9/L; hemoglobin mean 118 (SD 15.1) g/L (Table S1).",
    notes = "Exposure-safety population of 715 patients (J101 n = 182, U102 n = 219, U201 n = 274, U202 n = 40; Li 2025 Table S1); 698 were evaluable for this endpoint and 196 (28.1%) had an event (Table S4)."
  )

  ini({
    # Li 2025 prints no coefficient table. Values below are reconstructed
    # from the typical-value forest plot of Figure S3B under the
    # logit-log-linear exposure form of Supplementary methods (residuals at
    # most 0.010 logit across the four exposure rows; the linear form leaves
    # 0.40 logit). Continuous covariates are linear on the logit (residuals
    # at most 0.013 logit). Categorical shifts are read directly from the
    # printed probabilities.
    logit_ref <- log(0.184 / (1 - 0.184)); label("Logit of grade >= 3 neutrophil decrease probability for the reference patient (unitless logit)") # Figure S3B 'Reference patient' 18.4 [11.5, 25.7]
    e_cmax_logit <- 1.7082; label("Change in the logit per unit natural log of Cycle 1 DXd Cmax (unitless logit)") # back-solved from Figure S3B exposure rows 8.0/15.5/33.9/58.5 ng/mL -> 3.4/9.8/29.1/51.4 (reference 23.8 ng/mL implied)
    e_wt_logit <- 0.00878; label("Change in the logit per kg of baseline body weight (unitless logit)") # back-solved from Figure S3B rows 42/52/72/93 kg -> 16.1/17.2/20.2/23.0 (reference 60 kg)
    e_neut_logit <- -0.15536; label("Change in the logit per 10^9/L of baseline neutrophil count (unitless logit)") # back-solved from Figure S3B rows 2.1/3.0/5.8/9.7 -> 23.5/20.9/14.8/8.6 (reference 4.1 implied)
    e_hgb_logit <- -0.01835; label("Change in the logit per g/L of baseline hemoglobin (unitless logit)") # back-solved from Figure S3B rows 95/107/128/144 -> 25.6/21.8/15.9/12.3 (reference 118 implied)
    e_region_japan_logit <- log(0.491 / 0.509) - log(0.184 / 0.816); label("Shift in the logit for Japan vs non-Asian country of origin (unitless logit)") # Figure S3B 'Country of origin: Japanese' 49.1 [24.8, 62.7]
    e_region_asia_nonjapan_logit <- log(0.206 / 0.794) - log(0.184 / 0.816); label("Shift in the logit for non-Japanese Asian vs non-Asian country of origin (unitless logit)") # Figure S3B 'Country of origin: Non-Japanese Asian' 20.6 [12.2, 36.0]

    addSd_prob_anc_decrease <- fixed(0.001); label("Placeholder additive residual SD on the event probability (fraction)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # Reference values are those implied by the Figure S3B rows (back-solved).
    logit_anc_decrease <- logit_ref +
      e_cmax_logit * log(CMAX / 23.8) +
      e_wt_logit * (WT - 60) +
      e_neut_logit * (NEUT - 4.1) +
      e_hgb_logit * (HGB - 118) +
      e_region_japan_logit * REGION_JAPAN +
      e_region_asia_nonjapan_logit * REGION_ASIA_NONJAPAN
    prob_anc_decrease <- expit(logit_anc_decrease)
    prob_anc_decrease ~ add(addSd_prob_anc_decrease)
  })
}
