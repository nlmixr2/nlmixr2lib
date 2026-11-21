Li_2025_mosunetuzumab_crs_cov <- function() {
  description <- paste(
    "Reduced exploratory covariate logistic-regression model for the",
    "occurrence of at least one GRADE >= 2 CYTOKINE RELEASE SYNDROME event",
    "versus the maximum model-predicted CD20 receptor occupancy over the",
    "first two cycles (ROMAX, percent), baseline age > 64 years (AGE_GT64)",
    "and 'other' NHL histology (TUMTP_OTHER_NHL: mantle cell lymphoma,",
    "Richter's syndrome and rarer types, versus DLBCL / transformed FL /",
    "FL) in adults with relapsed/refractory B-cell NHL on intravenous",
    "mosunetuzumab with Cycle 1 step-up dosing (Li 2025, study GO29781",
    "Group B, n = 387; Figure S4). The odds ratios are printed in the",
    "Figure S4 forest plot; the intercept is not printed and was digitised",
    "from the plateaus of the Figure S10 predicted-probability curves."
  )
  reference <- paste(
    "Li C-C, Bender B, Wilkins J, Li F, Turner DC, Wang B, Deng R,",
    "Vadhavkar S, Li Z, Kwan A, Huang H, Peng K, Penuel E, Huw L-Y, Chanu P,",
    "Li C, Yin S, Wei MC. A Novel Step-Up Dosage Regimen for Enhancing the",
    "Benefit-to-Risk Ratio of Mosunetuzumab in Relapsed or Refractory",
    "Follicular Lymphoma. Clin Pharmacol Ther. 2025;117(2):465-474.",
    "doi:10.1002/cpt.3445. Intercept digitised from Figure S10."
  )
  vignette <- "Li_2025_mosunetuzumab"
  units <- list(
    time = "n/a (static two-cycle landmark regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the ROMAX data column)",
    concentration = "prob_crs_grade2 (probability of grade >= 2 cytokine release syndrome, 0-1)"
  )

  covariateData <- list(
    ROMAX = list(
      description = "Maximum model-predicted mosunetuzumab CD20 receptor occupancy over days 0-42",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Maximum model-predicted CD20 receptor occupancy (percent, 0-100)",
        "over days 0-42 (the first two 21-day cycles); the `RO` output of",
        "Bender_2024_mosunetuzumab maximised over days 0-42. Figure S4",
        "reports the odds ratio 'per 10%' of ROmax; model() converts it",
        "to the per-percentage-point slope. Enters linearly and uncentred."
      ),
      source_name = "ROmax (day 42)"
    ),
    AGE_GT64 = list(
      description = "Baseline age above 64 years indicator (1 = older than 64 years, 0 = 64 years or younger)",
      units = "(binary)",
      type = "binary",
      reference_category = "64 years or younger (AGE_GT64 = 0)",
      notes = paste(
        "Li 2025 Figure S4 row 'Baseline age (>64 y)'; Table S11 and",
        "Figure S10 split the cohort as '<=64' versus '>64' years (45.1",
        "percent of the safety population older than 64). Derive as",
        "as.integer(AGE > 64)."
      ),
      source_name = "Baseline age (>64 y)"
    ),
    TUMTP_OTHER_NHL = list(
      description = "'Other' NHL histology indicator (1 = neither DLBCL, transformed FL nor FL; 0 = DLBCL, transformed FL or FL)",
      units = "(binary)",
      type = "binary",
      reference_category = "DLBCL, transformed FL or FL (TUMTP_OTHER_NHL = 0)",
      notes = paste(
        "Li 2025 Figure S4 row 'Histology: other'. Per the Figure S4 and",
        "Table S9 footnotes the 'other' group (n = 51) is MCL 34,",
        "Richter's syndrome 10, MZL 2, and 1 each of DLBCL/MCL, FL grade",
        "3B, melanoma, SLL and transformed MZL; the reference pools FL",
        "(148), DLBCL (142) and transformed FL (46). Note that transformed",
        "FL is in the REFERENCE group here."
      ),
      source_name = "Histology: other"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 387,
    n_studies = 1,
    age_range = "19-96 years (median 63; whole safety population, Table S11)",
    weight_range = "37.1-163 kg (median 77.9)",
    sex_female_pct = 35.3,
    race_ethnicity = c(
      White = 76.2,
      Asian = 17.6,
      Black = 2.75,
      `American Indian or Alaska Native` = 0.458,
      Unknown = 2.97
    ),
    disease_state = "Relapsed/refractory B-cell non-Hodgkin lymphoma (FL 148, DLBCL 142, transformed FL 46, other histologies 51)",
    dose_range = "Cycle 1 step-up dosing from 0.4/1.0/2.8 to 1/2/60 and 1/2/60/30 mg IV q3w (Group B)",
    regions = "Multinational",
    notes = paste(
      "Cohort B (step-up dosing) of the phase I/II GO29781 study",
      "(NCT02500407), ROmax-evaluable; the histology counts of the Figure",
      "S4 footnote sum to 387. Candidate covariates (Table S9) were Ann",
      "Arbor stage, log LDH, sqrt(SPD), histology, age > 64 y, white",
      "blood cell count, bone-marrow involvement and cellularity, and",
      "prior CAR-T therapy; stepwise forward addition (dOFV >= 3.84) and",
      "backward deletion (dOFV <= 6.63) retained ROmax0-42, age and",
      "histology. Demographics are those of the whole 439-patient safety",
      "population (Table S11)."
    )
  )

  ini({
    # ======================================================================
    # Li 2025 prints the three odds ratios (Figure S4 forest plot) but not
    # the intercept. The intercept was recovered from Figure S10, which
    # plots the model-predicted probability against log10(ROmax0-42) for the
    # four age x histology strata. The curves were pixel-traced on the
    # calibrated probability axis; each is flat at low ROmax (the ROmax term
    # is < 0.002 logit at the leftmost 0.03 percent), giving four plateau
    # probabilities:
    #   DLBCL/trFL/FL, age <= 64: 0.0595 -> logit -2.760 (= beta0)
    #   DLBCL/trFL/FL, age  > 64: 0.0194 -> logit -3.923 (- ln 0.323 = -2.793)
    #   other,         age <= 64: 0.285  -> logit -0.920 (- ln 6.193 = -2.743)
    #   other,         age  > 64: 0.110  -> logit -2.096 (-> -2.789)
    # The four intercept estimates agree within 0.05 logit; their mean is
    # -2.77. The plateau offsets also reproduce the printed odds ratios
    # (-1.163 vs ln 0.323 = -1.130; 1.840 vs ln 6.193 = 1.823), confirming
    # the reading.
    # ======================================================================
    logite0 <- -2.77; label("Logit of the probability of grade >= 2 CRS at ROMAX = 0 percent, age <= 64 years, DLBCL/trFL/FL histology (unitless logit)") # digitised from the four Figure S10 curve plateaus (mean of -2.760, -2.793, -2.743, -2.789)
    e_romax_logit <- log(1.289) / 10; label("Log-odds of grade >= 2 CRS per 1 percentage-point increase in ROMAX (unitless logit)") # Li 2025 Figure S4, 'ROmax (day 42, per 10%)' OR = 1.289 (p < 0.01), divided by 10 for the per-1% slope
    e_age_gt64_logit <- log(0.323); label("Log-odds of grade >= 2 CRS for age > 64 years versus <= 64 years (unitless logit)") # Li 2025 Figure S4, 'Baseline age (>64 y)' OR = 0.323 (p < 0.01)
    e_tumtp_other_nhl_logit <- log(6.193); label("Log-odds of grade >= 2 CRS for 'other' histology versus DLBCL/trFL/FL (unitless logit)") # Li 2025 Figure S4, 'Histology: other' OR = 6.193 (p < 0.001)

    addSd_prob_crs_grade2 <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # Li 2025 Methods, linear logistic regression with drug-independent
    # covariates:
    #   logit(p) = beta0 + beta_ROmax * ROmax + beta_1 * x_1 + beta_2 * x_2
    logit_crs <- logite0 + e_romax_logit * ROMAX +
      e_age_gt64_logit * AGE_GT64 + e_tumtp_other_nhl_logit * TUMTP_OTHER_NHL

    prob_crs_grade2 <- expit(logit_crs)

    prob_crs_grade2 ~ add(addSd_prob_crs_grade2)
  })
}
