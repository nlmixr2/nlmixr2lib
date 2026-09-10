Chan_2025_atezolizumab_isr <- function() {
  description <- paste0(
    "Binomial logistic-regression exposure-safety model for INJECTION ",
    "SITE REACTIONS in adults with locally advanced or metastatic ",
    "non-small cell lung cancer receiving SUBCUTANEOUS atezolizumab ",
    "1875 mg every 3 weeks (Chan 2025, n = 246, cohort 5 of the phase ",
    "III portion of IMscin001, NCT03735121, clinical cut-off 26 April ",
    "2022). The probability of the event is ",
    "expit(-4.50 + 0.581 * (AUC_ATEZO / 1000) - 2.55 * RACE_HISPANIC + ",
    "2.39 * RACE_HISPANIC_MISSING) where AUC_ATEZO is the ",
    "model-predicted Cycle-1 AUC from day 0 to day 21 in ug*day/mL and ",
    "the two indicators are the non-reference levels of a three-level ",
    "ethnicity factor whose reference is not-Hispanic-or-Latino. THE ",
    "EXPOSURE TERM IS NOT STATISTICALLY SIGNIFICANT (odds ratio 1.79 per ",
    "1000 ug*day/mL, 95 percent CI 0.933-3.40, p = 0.0793), though it is ",
    "the closest any exposure term in this paper comes to significance. ",
    "This endpoint is UNIQUE TO THE SUBCUTANEOUS ROUTE and has no ",
    "intravenous counterpart; conversely no patient in cohort 5 ",
    "experienced an infusion-related reaction, so the paper reports no ",
    "IRR model at all. Chan 2025 explicitly cautions that the ",
    "significant protective Hispanic effect rests on a small subgroup ",
    "(N = 61) and 'the physiological cause of this covariate effect ",
    "remains unclear and may be due to potential confounding'. There is ",
    "no PK layer and no ODE: the exposure metric is supplied as a data ",
    "column, derived in the source analysis from the individual ",
    "empirical-Bayes predictions of the companion population PK model ",
    "packaged as Chan_2025_atezolizumab. No between-subject random ",
    "effect and no residual error are estimated (Bernoulli likelihood). ",
    "Five companion exposure-response models in the ",
    "Chan_2025_atezolizumab_* family."
  )
  reference <- paste(
    "Chan P, Liu SN, Gosselin N, Sauve Z, Marchand M, Lin A,",
    "Herraez-Baranda L, Zanghi J, Shearer-Kang E, Liu X, Wu B, Chanu P.",
    "Population pharmacokinetics and exposure-response of subcutaneous",
    "atezolizumab in patients with non-small cell lung cancer.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(4):726-737.",
    "doi:10.1002/psp4.13310.",
    "Coefficients are transcribed from Supporting Information Table S10G.",
    "Individual Cycle-1 AUC0-21d values derive from the companion",
    "population PK model reported in the same paper; see",
    "modellib('Chan_2025_atezolizumab').",
    sep = " "
  )
  vignette <- "Chan_2025_atezolizumab_sc_nsclc"
  units <- list(
    time          = "n/a (static Cycle-1 landmark safety regression; no time dimension)",
    dosing        = "n/a (no dose events; exposure enters as the AUC_ATEZO data column)",
    concentration = "prob_isr (probability of an injection site reaction, 0-1; also logit_isr)"
  )

  covariateData <- list(
    AUC_ATEZO = list(
      description        = "Model-predicted Cycle-1 atezolizumab AUC from day 0 to day 21",
      units              = "ug*day/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters LINEARLY and is SCALED BY 1000, not centred: Chan 2025",
        "Table S10G labels the row 'AUC0-21d (Increase of 1000",
        "ug*day/mL)', so the coefficient 0.581 is the log-odds per 1000",
        "ug*day/mL and the intercept -4.50 is the logit at",
        "AUC_ATEZO = 0, outside the observed range. This is the largest",
        "exposure coefficient of the four safety endpoints and the only",
        "one with p < 0.10, which is biologically unsurprising for an",
        "injection-site endpoint whose driver is the subcutaneously",
        "administered amount; it nonetheless fails the paper's p < 0.05",
        "threshold. Cohort 5 distribution (Chan 2025 Table S5B):",
        "geometric mean 2907 ug*day/mL (geoCV 35.9 percent), median 2974,",
        "range 666-6572."
      ),
      source_name        = "AUC0-21d"
    ),
    RACE_HISPANIC = list(
      description        = "Hispanic or Latino ethnicity indicator (1 = Hispanic or Latino, 0 = otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not Hispanic or Latino, the model's reference level)",
      notes              = paste(
        "The non-reference level of a THREE-level ethnicity factor whose",
        "levels are not-Hispanic-or-Latino (reference), Hispanic or",
        "Latino, and Unknown. Chan 2025 Table S10G prints rows for only",
        "the latter two, which is how a treatment-contrast-coded factor",
        "is normally reported, so the omitted level is the reference.",
        "Strongly PROTECTIVE and statistically significant here: odds",
        "ratio 0.0782 (95 percent CI 0.000586-0.689, p = 0.0161). The",
        "Chan 2025 Discussion flags this as a small-subgroup finding",
        "(N = 61 Hispanic patients) whose 'physiological cause remains",
        "unclear and may be due to potential confounding with other",
        "statistically significant covariates'. This paper reports",
        "ethnicity as its own dimension separate from race, so the",
        "per-model interpretation of the canonical RACE_HISPANIC column",
        "here is strictly ETHNICITY, matching the precedent set by",
        "Overgaard_2016_liraglutide.R. Note that ethnicity was also",
        "screened on bioavailability in the companion population PK model",
        "and was NOT retained there (control-stream block F1ETHN, $THETA",
        "0 FIX)."
      ),
      source_name        = "ethnic: Hispanic or Latino"
    ),
    RACE_HISPANIC_MISSING = list(
      description        = "Ethnicity recorded as Unknown indicator (1 = ethnicity unknown / not reported, 0 = ethnicity recorded)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ethnicity recorded; mutually exclusive with RACE_HISPANIC)",
      notes              = paste(
        "The third level of the ethnicity factor, carried as its own",
        "indicator rather than collapsed onto the not-Hispanic reference,",
        "exactly as ADA_MISSING is carried alongside ADA_POS in",
        "Suri_2018_brentuximab.R. The point estimate is large and",
        "positive (2.39, odds ratio 10.9) but wholly uninformative: the",
        "95 percent CI spans 0.0684 to 258 and p = 0.2559, so this row",
        "records a missing-data stratum the regression could not resolve,",
        "not a real risk factor. It is retained because dropping it would",
        "silently change the reference group of RACE_HISPANIC from",
        "'recorded as not Hispanic' to 'not recorded as Hispanic",
        "(including unknown)' and thereby bias the Hispanic coefficient a",
        "downstream user reproduces. Interpret with caution: the",
        "missingness mechanism is not stated."
      ),
      source_name        = "ethnic: Unknown"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 246,
    n_studies      = 1,
    disease_state  = "locally advanced or metastatic non-small cell lung cancer",
    race_ethnicity = "61 of 246 patients (24.8 percent) are Hispanic or Latino (Chan 2025 Discussion); the remainder are not-Hispanic-or-Latino or Unknown",
    dose_range     = "subcutaneous atezolizumab 1875 mg every 3 weeks in the thigh, as a ready-to-use co-formulation with recombinant human hyaluronidase PH20",
    notes          = paste(
      "Cohort 5 of the phase III (part 2) portion of IMscin001",
      "(NCT03735121). This endpoint exists only in the subcutaneous arm.",
      "Chan 2025 Results records that 'No patient who received",
      "atezolizumab SC in the Phase II portion of IMscin001 (Cohort 5)",
      "experienced an IRR', so the paper's infusion-related-reaction",
      "endpoint yielded no fitted model and none is packaged here.",
      "Baseline medians (Chan 2025 Table 1): body weight 67.8 kg, tumor",
      "burden 79.5 mm, albumin 40.0 g/L, hemoglobin 123 g/L; 29.3 percent",
      "female. Cycle-1 exposure metrics were used deliberately 'to",
      "minimize the potential effect of response-dependent time-varying",
      "clearance'."
    )
  )

  ini({
    # ==================================================================
    # Chan 2025 Supporting Information Table S10G, "Logistic Regression
    # Model for ISR and Atezolizumab Exposure with Relevant Covariates in
    # Cohort 5 (Atezolizumab SC 1875 mg Q3W)".
    #
    #   Intercept                             -4.50               p < 0.001
    #   AUC0-21d (Increase of 1000 ug*day/mL)  0.581  OR 1.79      p = 0.0793
    #                                                 (0.933, 3.40)
    #   ethnic: Hispanic or Latino            -2.55   OR 0.0782    p = 0.0161
    #                                                 (0.000586, 0.689)
    #   ethnic: Unknown                        2.39   OR 10.9      p = 0.2559
    #                                                 (0.0684, 258)
    #
    # All three slopes are printed as Estimates on the log-odds scale, so
    # nothing here is back-solved or digitised. The printed odds ratios
    # are an internal check: exp(0.581) = 1.7878 rounds to 1.79,
    # exp(-2.55) = 0.078082 rounds to 0.0782, and exp(2.39) = 10.914
    # rounds to 10.9, all matching Table S10G.
    #
    # The two ethnicity rows are treatment contrasts against an omitted
    # not-Hispanic-or-Latino reference; see the covariateData notes.
    # ==================================================================
    logit_ref <- -4.50; label("Logit of the probability of an injection site reaction at AUC_ATEZO = 0 ug*day/mL in a not-Hispanic-or-Latino patient (unitless logit)")  # Chan 2025 Table S10G, Intercept row: Estimate -4.50, p < 0.001
    e_auc_atezo_logit <- 0.581; label("Log-odds of an injection site reaction per 1000 ug*day/mL increase in Cycle-1 AUC0-21d (unitless logit)")                        # Chan 2025 Table S10G, AUC0-21d row: Estimate 0.581, OR 1.79 (95% CI 0.933-3.40), p = 0.0793 (not significant)
    e_hispanic_logit <- -2.55; label("Log-odds of an injection site reaction for Hispanic or Latino ethnicity versus not Hispanic or Latino (unitless logit)")          # Chan 2025 Table S10G, "ethnic: Hispanic or Latino" row: Estimate -2.55, OR 0.0782 (95% CI 0.000586-0.689), p = 0.0161
    e_hispanic_missing_logit <- 2.39; label("Log-odds of an injection site reaction for unknown ethnicity versus not Hispanic or Latino (unitless logit)")              # Chan 2025 Table S10G, "ethnic: Unknown" row: Estimate 2.39, OR 10.9 (95% CI 0.0684-258), p = 0.2559

    # No between-subject variability and no residual error: the source is
    # a Bernoulli-likelihood logistic regression fitted in R 4.1.1, and
    # Table S10G reports no variance components. The placeholder additive
    # term exists only because rxode2 requires an observation
    # declaration; see the vignette Assumptions and deviations.
    addSd_prob_isr <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # ------------------------------------------------------------------
    # Linear predictor. The exposure is not centred, so logit_ref is the
    # logit at AUC_ATEZO = 0 in a not-Hispanic-or-Latino patient, which
    # is outside the observed data; only DIFFERENCES in the linear
    # predictor between two covariate settings are interpretable.
    #
    # The divisor 1000 converts the raw canonical data column
    # (ug*day/mL) into the per-1000 units in which Table S10G reports its
    # coefficient. RACE_HISPANIC and RACE_HISPANIC_MISSING are mutually
    # exclusive; setting both to 0 selects the reference level.
    # ------------------------------------------------------------------
    logit_isr <- logit_ref +
      e_auc_atezo_logit * (AUC_ATEZO / 1000) +
      e_hispanic_logit * RACE_HISPANIC +
      e_hispanic_missing_logit * RACE_HISPANIC_MISSING

    prob_isr <- expit(logit_isr)

    # ------------------------------------------------------------------
    # Observation. Deterministic probability of an injection site
    # reaction. Downstream callers can sample binary outcomes with
    # rbinom(n, 1, prob_isr) on the rxSolve output.
    # ------------------------------------------------------------------
    prob_isr ~ add(addSd_prob_isr)
  })
}
