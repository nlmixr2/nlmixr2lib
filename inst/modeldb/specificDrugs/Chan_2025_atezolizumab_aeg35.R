Chan_2025_atezolizumab_aeg35 <- function() {
  description <- paste0(
    "Binomial logistic-regression exposure-safety model for GRADE 3 TO 5 ",
    "ADVERSE EVENTS in adults with locally advanced or metastatic ",
    "non-small cell lung cancer receiving SUBCUTANEOUS atezolizumab ",
    "1875 mg every 3 weeks (Chan 2025, n = 246, cohort 5 of the phase ",
    "III portion of IMscin001, NCT03735121, clinical cut-off 26 April ",
    "2022). The probability of the event is ",
    "expit(-2.16 + 0.0879 * (AUC_ATEZO / 1000) + 0.166 * (CRP / 10)) ",
    "where AUC_ATEZO is the model-predicted Cycle-1 AUC from day 0 to ",
    "day 21 in ug*day/mL and CRP is baseline C-reactive protein in mg/L. ",
    "THE EXPOSURE TERM IS NOT STATISTICALLY SIGNIFICANT (odds ratio 1.09 ",
    "per 1000 ug*day/mL, 95 percent CI 0.785-1.52, p = 0.602) and that ",
    "flat exposure-response is the paper's headline result, not a ",
    "transcription gap. Baseline C-reactive protein, by contrast, IS ",
    "significant (odds ratio 1.18 per 10 mg/L, 95 percent CI 1.10-1.27, ",
    "p = 0.001). Note that this endpoint retained AUC0-21d whereas the ",
    "companion serious-adverse-event model retained Cmax: Chan 2025 ",
    "selected, per endpoint, whichever of the two safety exposure ",
    "metrics gave the lowest p-value. There is no PK layer and no ODE: ",
    "the exposure metric is supplied as a data column, derived in the ",
    "source analysis from the individual empirical-Bayes predictions of ",
    "the companion population PK model packaged as ",
    "Chan_2025_atezolizumab. No between-subject random effect and no ",
    "residual error are estimated (Bernoulli likelihood). Five companion ",
    "exposure-response models in the Chan_2025_atezolizumab_* family."
  )
  reference <- paste(
    "Chan P, Liu SN, Gosselin N, Sauve Z, Marchand M, Lin A,",
    "Herraez-Baranda L, Zanghi J, Shearer-Kang E, Liu X, Wu B, Chanu P.",
    "Population pharmacokinetics and exposure-response of subcutaneous",
    "atezolizumab in patients with non-small cell lung cancer.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(4):726-737.",
    "doi:10.1002/psp4.13310.",
    "Coefficients are transcribed from Supporting Information Table S10F.",
    "Individual Cycle-1 AUC0-21d values derive from the companion",
    "population PK model reported in the same paper; see",
    "modellib('Chan_2025_atezolizumab').",
    sep = " "
  )
  vignette <- "Chan_2025_atezolizumab_sc_nsclc"
  units <- list(
    time          = "n/a (static Cycle-1 landmark safety regression; no time dimension)",
    dosing        = "n/a (no dose events; exposure enters as the AUC_ATEZO data column)",
    concentration = "prob_aeg35 (probability of a grade 3-5 adverse event, 0-1; also logit_aeg35)"
  )

  covariateData <- list(
    AUC_ATEZO = list(
      description        = "Model-predicted Cycle-1 atezolizumab AUC from day 0 to day 21",
      units              = "ug*day/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters LINEARLY and is SCALED BY 1000, not centred: Chan 2025",
        "Table S10F labels the row 'AUC0-21d (Increase by 1000",
        "ug*day/mL)', so the coefficient 0.0879 is the log-odds per 1000",
        "ug*day/mL and the intercept -2.16 is the logit at",
        "AUC_ATEZO = 0, outside the observed range. Derived from the",
        "companion population PK model, NOT observed: this is the",
        "co-primary endpoint of the registrational study, and Chan 2025",
        "states it is 'the first popPK analysis to derive a measure for",
        "the primary analysis of a pivotal, Phase III PK non-inferiority",
        "study'. Cohort 5 distribution (Chan 2025 Table S5B): geometric",
        "mean 2907 ug*day/mL (geoCV 35.9 percent), median 2974, range",
        "666-6572."
      ),
      source_name        = "AUC0-21d"
    ),
    CRP = list(
      description        = "Baseline C-reactive protein",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters LINEARLY and is SCALED BY 10, not centred: Chan 2025",
        "Table S10F labels the row 'CRP (Increase by 10 mg/L)', so the",
        "coefficient 0.166 is the log-odds per 10 mg/L (odds ratio 1.18).",
        "This is the only significant term in the model (p = 0.001) and",
        "the Chan 2025 Discussion singles CRP out across endpoints. It",
        "also appears in the companion SAE, PFS and OS models. Reported",
        "in mg/L (SI)."
      ),
      source_name        = "CRP"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 246,
    n_studies      = 1,
    disease_state  = "locally advanced or metastatic non-small cell lung cancer",
    dose_range     = "subcutaneous atezolizumab 1875 mg every 3 weeks in the thigh, as a ready-to-use co-formulation with recombinant human hyaluronidase PH20",
    notes          = paste(
      "Cohort 5 of the phase III (part 2) portion of IMscin001",
      "(NCT03735121). Exposure-response was assessed on the subcutaneous",
      "arm ONLY -- Chan 2025 Methods notes that intravenous",
      "exposure-response has been evaluated extensively elsewhere and was",
      "therefore excluded. 'Data from all PK-evaluable patients were",
      "included in the ER analysis, and no additional exclusion was",
      "applied.' Baseline medians (Chan 2025 Table 1): body weight",
      "67.8 kg, tumor burden 79.5 mm, albumin 40.0 g/L, hemoglobin",
      "123 g/L; 29.3 percent female. Cycle-1 exposure metrics were used",
      "deliberately 'to minimize the potential effect of",
      "response-dependent time-varying clearance'."
    )
  )

  ini({
    # ==================================================================
    # Chan 2025 Supporting Information Table S10F, "Logistic Regression
    # Model for AEG35 and Atezolizumab Exposure with Relevant Covariates
    # in Cohort 5 (Atezolizumab SC 1875 mg Q3W)".
    #
    #   Intercept                             -2.16              p = 0.001
    #   AUC0-21d (Increase by 1000 ug*day/mL)  0.0879  OR 1.09    p = 0.602
    #                                                  (0.785, 1.52)
    #   CRP (Increase by 10 mg/L)              0.166   OR 1.18    p = 0.001
    #                                                  (1.10, 1.27)
    #
    # Both slopes are printed as Estimates on the log-odds scale, so
    # nothing here is back-solved or digitised. The printed odds ratios
    # are an internal check: exp(0.0879) = 1.0919 rounds to 1.09 and
    # exp(0.166) = 1.1806 rounds to 1.18, both matching Table S10F.
    #
    # NEITHER covariate is centred; both are RESCALED (per 1000
    # ug*day/mL and per 10 mg/L). The model() block divides the raw data
    # columns by those constants so the coefficients below are used
    # exactly as printed.
    # ==================================================================
    logit_ref <- -2.16; label("Logit of the probability of a grade 3-5 adverse event at AUC_ATEZO = 0 ug*day/mL and CRP = 0 mg/L (unitless logit)")  # Chan 2025 Table S10F, Intercept row: Estimate -2.16, p = 0.001
    e_auc_atezo_logit <- 0.0879; label("Log-odds of a grade 3-5 adverse event per 1000 ug*day/mL increase in Cycle-1 AUC0-21d (unitless logit)")     # Chan 2025 Table S10F, AUC0-21d row: Estimate 0.0879, OR 1.09 (95% CI 0.785-1.52), p = 0.602 (not significant)
    e_crp_logit <- 0.166; label("Log-odds of a grade 3-5 adverse event per 10 mg/L increase in baseline C-reactive protein (unitless logit)")        # Chan 2025 Table S10F, CRP row: Estimate 0.166, OR 1.18 (95% CI 1.10-1.27), p = 0.001

    # No between-subject variability and no residual error: the source is
    # a Bernoulli-likelihood logistic regression fitted in R 4.1.1, and
    # Table S10F reports no variance components. The placeholder additive
    # term exists only because rxode2 requires an observation
    # declaration; see the vignette Assumptions and deviations.
    addSd_prob_aeg35 <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # ------------------------------------------------------------------
    # Linear predictor. Neither covariate is centred, so logit_ref is the
    # logit at AUC_ATEZO = 0 and CRP = 0, which is outside the observed
    # data; only DIFFERENCES in the linear predictor between two
    # covariate settings are interpretable.
    #
    # The divisors 1000 and 10 convert the raw canonical data columns
    # (ug*day/mL and mg/L) into the per-1000 and per-10 units in which
    # Table S10F reports its coefficients.
    # ------------------------------------------------------------------
    logit_aeg35 <- logit_ref +
      e_auc_atezo_logit * (AUC_ATEZO / 1000) +
      e_crp_logit * (CRP / 10)

    prob_aeg35 <- expit(logit_aeg35)

    # ------------------------------------------------------------------
    # Observation. Deterministic probability of a grade 3-5 adverse
    # event. Downstream callers can sample binary outcomes with
    # rbinom(n, 1, prob_aeg35) on the rxSolve output.
    # ------------------------------------------------------------------
    prob_aeg35 ~ add(addSd_prob_aeg35)
  })
}
