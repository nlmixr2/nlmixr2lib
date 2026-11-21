Li_2025_patritumabDeruxtecan_plt_decrease <- function() {
  description <- paste(
    "Bayesian logistic-regression exposure-safety model for grade >= 3",
    "platelet count decrease in adults with breast cancer, non-small-cell",
    "lung cancer or colorectal cancer given patritumab deruxtecan (HER3-DXd)",
    "in the J101, U102, U201 and U202 studies (Li 2025, n = 714;",
    "1.6-8 mg/kg IV Q3W, up-titration or an alternative Q2W/Q3W regimen).",
    "logit(P) is log-linear in the Cycle 1 DXd (unconjugated payload) Cmax,",
    "linear in baseline platelet count and baseline hemoglobin, with",
    "additive shifts for race (Asian, Black, Other vs White) and tumor type",
    "(breast cancer, colorectal cancer, EGFR-wild-type NSCLC vs",
    "EGFR-mutated NSCLC). Li 2025 prints no coefficient table; every value",
    "is back-solved from the typical-value forest plot (Figure S3A).",
    "Covariate-by-exposure interaction terms of the source model are not",
    "recoverable and are omitted. There is no PK layer: exposure is a data",
    "column from the companion population PK model (see",
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
    concentration = "prob_plt_decrease (probability of a grade >= 3 platelet count decrease, 0-1)"
  )

  covariateData <- list(
    CMAX = list(
      description = "Cycle 1 maximum plasma concentration of DXd, the unconjugated topoisomerase I inhibitor payload, per subject.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Unconjugated payload (DXd), not the conjugate. Cycle 1 (first-dose) landmark, because grade >= 3 platelet decreases occurred early (median onset day 8, Li 2025 Table S4). Empirical-Bayes predictions from the population PK model 'previously developed' (uncited in Li 2025; the same group's model is Lee_2023_patritumab). Centred at 23.5 ng/mL inside model(), the reference exposure implied by Figure S3A; Table S5 reports a median of 24.2 ng/mL at 5.6 mg/kg Q3W. Figure S3A 5th/25th/75th/95th percentiles: 8.0 / 15.5 / 33.9 / 58.5 ng/mL.",
      source_name = "Cycle 1 DXd Cmax (ng/mL)"
    ),
    PLT = list(
      description = "Baseline platelet count.",
      units = "10^9/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline only. Significant covariate (Li 2025 Table 2). Centred at 243 x 10^9/L inside model(), the reference value implied by Figure S3A; Table S1 mean 259 (SD 92.3). Figure S3A 5th/25th/75th/95th percentiles: 138 / 193 / 308 / 427.",
      source_name = "Baseline platelets (10^9/L)"
    ),
    HGB = list(
      description = "Baseline hemoglobin concentration.",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline only. Significant covariate (Li 2025 Table 2). Figure S3A and Table S1 label the unit g/dL, but the printed values (5th-95th percentile 95-144; mean 118, SD 15.1) are in g/L. Supply g/L. Centred at 118 g/L inside model(), the reference value implied by Figure S3A and the Table S1 mean.",
      source_name = "Baseline hemoglobin (g/dL as labelled; values in g/L)"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator (1 = Asian, 0 = otherwise).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White when RACE_BLACK and RACE_OTHER are also 0).",
      notes = "Significant covariate. Figure S3A 'Race: Asian' 29.5% vs reference 11.3%. 54.7% of the safety population was Asian (Table S1).",
      source_name = "Race: Asian"
    ),
    RACE_BLACK = list(
      description = "Black or African American race indicator (1 = Black, 0 = otherwise).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White when RACE_ASIAN and RACE_OTHER are also 0).",
      notes = "Figure S3A 'Race: Black' 11.4% vs reference 11.3%; 2.1% of the safety population (Table S1).",
      source_name = "Race: Black"
    ),
    RACE_OTHER = list(
      description = "Other race indicator (1 = race other than White, Asian or Black, 0 = otherwise).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White when RACE_ASIAN and RACE_BLACK are also 0).",
      notes = "Figure S3A 'Race: Other' 11.2% vs reference 11.3%; 5.7% of the safety population (Table S1).",
      source_name = "Race: Other"
    ),
    TUMTP_BREAST = list(
      description = "Breast cancer tumor-type indicator (1 = breast cancer, 0 = otherwise).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (EGFR-mutated NSCLC, the Figure S3A reference tumor type).",
      notes = "Influential but not significant. Figure S3A 'Tumor type: Breast cancer' 18.9% vs reference 11.3%. All 182 J101 patients (25.5% of the safety population).",
      source_name = "Tumor type: Breast cancer"
    ),
    TUMTP_CRC = list(
      description = "Colorectal cancer tumor-type indicator (1 = colorectal cancer, 0 = otherwise).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (EGFR-mutated NSCLC).",
      notes = "Figure S3A 'Tumor type: CRC' 10.8% vs reference 11.3%. All 40 U202 patients (5.6%).",
      source_name = "Tumor type: CRC"
    ),
    TUMTP_NSCLC = list(
      description = "Non-small-cell lung cancer tumor-type indicator (1 = NSCLC, 0 = otherwise).",
      units = "(binary)",
      type = "binary",
      reference_category = "1 with TUM_EGFR_MUT = 1 is the reference tumor type.",
      notes = "Used only through the EGFR-wild-type NSCLC contrast TUMTP_NSCLC * (1 - TUM_EGFR_MUT). Set TUMTP_NSCLC = 1 for both NSCLC strata.",
      source_name = "Tumor type: NSCLC"
    ),
    TUM_EGFR_MUT = list(
      description = "Tumour EGFR mutation indicator (1 = EGFR-mutated, 0 = EGFR wild type).",
      units = "(binary)",
      type = "binary",
      reference_category = "1 in THIS model (EGFR-mutated NSCLC is the reference tumor type); the contrast is applied as TUMTP_NSCLC * (1 - TUM_EGFR_MUT).",
      notes = "Figure S3A 'Tumor type: NSCLC non-EGFRm' 10.9% vs reference 11.3%. 47 of 219 U102 patients had EGFR-wild-type NSCLC (6.6% of the safety population; Table S1). The value is irrelevant for non-NSCLC tumor types.",
      source_name = "NSCLC non-EGFRm"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 714L,
    n_studies = 4L,
    age_range = "29-83 years (median 61)",
    weight_range = "32.4-137 kg (median 60.0)",
    sex_female_pct = 69.2,
    race_ethnicity = c(Asian = 54.7, White = 37.5, Other = 5.7, Black = 2.1),
    disease_state = "HER3-expressing advanced breast cancer (J101), metastatic or unresectable NSCLC (U102, U201; EGFR-mutated in 446 of 493) and advanced colorectal cancer (U202).",
    dose_range = "Patritumab deruxtecan 1.6-8 mg/kg IV Q3W, up-titration (3.2 / 4.8 / 6.4 mg/kg) or 4.2 mg/kg Q2W for Cycles 1-3 then 6.4 mg/kg Q3W.",
    regions = "Global; country of origin Japan 34.4%, other Asian 13.3%, non-Asian 52.3%.",
    biomarkers = "Baseline platelets mean 259 (SD 92.3) x 10^9/L; hemoglobin mean 118 (SD 15.1) g/L (Table S1).",
    notes = "Exposure-safety population of 715 patients (J101 n = 182, U102 n = 219, U201 n = 274, U202 n = 40; Li 2025 Table S1); 714 were evaluable for this endpoint and 162 (22.7%) had an event (Table S4)."
  )

  ini({
    # Li 2025 prints no coefficient table. Values below are reconstructed
    # from the typical-value forest plot of Figure S3A under the
    # logit-log-linear exposure form of Supplementary methods (residuals at
    # most 0.019 logit across the four exposure rows; the linear form leaves
    # 0.33 logit). Continuous covariates are linear on the logit (residuals
    # at most 0.012 logit). Categorical shifts are read directly from the
    # printed probabilities.
    logit_ref <- log(0.113 / (1 - 0.113)); label("Logit of grade >= 3 platelet decrease probability for the reference patient (unitless logit)") # Figure S3A 'Reference patient' 11.3 [6.7, 16.9]
    e_cmax_logit <- 1.4613; label("Change in the logit per unit natural log of Cycle 1 DXd Cmax (unitless logit)") # back-solved from Figure S3A exposure rows 8.0/15.5/33.9/58.5 ng/mL -> 2.6/6.4/17.6/32.9 (reference 23.5 ng/mL implied)
    e_plt_logit <- -0.00751; label("Change in the logit per 10^9/L of baseline platelet count (unitless logit)") # back-solved from Figure S3A rows 138/193/308/427 -> 21.9/15.7/7.3/3.1 (reference 243 implied)
    e_hgb_logit <- -0.01728; label("Change in the logit per g/L of baseline hemoglobin (unitless logit)") # back-solved from Figure S3A rows 95/107/128/144 -> 16.1/13.4/9.7/7.6 (reference 118 implied)
    e_race_asian_logit <- log(0.295 / 0.705) - log(0.113 / 0.887); label("Shift in the logit for Asian vs White race (unitless logit)") # Figure S3A 'Race: Asian' 29.5 [18.2, 42.7]
    e_race_black_logit <- log(0.114 / 0.886) - log(0.113 / 0.887); label("Shift in the logit for Black vs White race (unitless logit)") # Figure S3A 'Race: Black' 11.4 [5.8, 19.6]
    e_race_other_logit <- log(0.112 / 0.888) - log(0.113 / 0.887); label("Shift in the logit for Other vs White race (unitless logit)") # Figure S3A 'Race: Other' 11.2 [6.2, 17.4]
    e_tumtp_breast_logit <- log(0.189 / 0.811) - log(0.113 / 0.887); label("Shift in the logit for breast cancer vs EGFR-mutated NSCLC (unitless logit)") # Figure S3A 'Tumor type: Breast cancer' 18.9 [10.3, 33.6]
    e_tumtp_crc_logit <- log(0.108 / 0.892) - log(0.113 / 0.887); label("Shift in the logit for colorectal cancer vs EGFR-mutated NSCLC (unitless logit)") # Figure S3A 'Tumor type: CRC' 10.8 [5.6, 16.8]
    e_nsclc_egfrwt_logit <- log(0.109 / 0.891) - log(0.113 / 0.887); label("Shift in the logit for EGFR-wild-type vs EGFR-mutated NSCLC (unitless logit)") # Figure S3A 'Tumor type: NSCLC non-EGFRm' 10.9 [5.6, 17.1]

    addSd_prob_plt_decrease <- fixed(0.001); label("Placeholder additive residual SD on the event probability (fraction)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # Reference values are those implied by the Figure S3A rows (back-solved).
    logit_plt_decrease <- logit_ref +
      e_cmax_logit * log(CMAX / 23.5) +
      e_plt_logit * (PLT - 243) +
      e_hgb_logit * (HGB - 118) +
      e_race_asian_logit * RACE_ASIAN +
      e_race_black_logit * RACE_BLACK +
      e_race_other_logit * RACE_OTHER +
      e_tumtp_breast_logit * TUMTP_BREAST +
      e_tumtp_crc_logit * TUMTP_CRC +
      e_nsclc_egfrwt_logit * TUMTP_NSCLC * (1 - TUM_EGFR_MUT)
    prob_plt_decrease <- expit(logit_plt_decrease)
    prob_plt_decrease ~ add(addSd_prob_plt_decrease)
  })
}
