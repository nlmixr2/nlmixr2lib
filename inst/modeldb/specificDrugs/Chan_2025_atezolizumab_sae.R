Chan_2025_atezolizumab_sae <- function() {
  description <- paste0(
    "Binomial logistic-regression exposure-safety model for SERIOUS ",
    "ADVERSE EVENTS in adults with locally advanced or metastatic ",
    "non-small cell lung cancer receiving SUBCUTANEOUS atezolizumab ",
    "1875 mg every 3 weeks (Chan 2025, n = 246, cohort 5 of the phase ",
    "III portion of IMscin001, NCT03735121, clinical cut-off 26 April ",
    "2022). The probability of the event is ",
    "expit(-2.69 + 0.210 * (CMAX / 100) + 0.130 * (CRP / 10)) where CMAX ",
    "is the model-predicted Cycle-1 maximum serum atezolizumab ",
    "concentration in ug/mL and CRP is baseline C-reactive protein in ",
    "mg/L. THE EXPOSURE TERM IS NOT STATISTICALLY SIGNIFICANT (odds ",
    "ratio 1.23 per 100 ug/mL, 95 percent CI 0.723-2.11, p = 0.441) and ",
    "that flat exposure-response is the paper's headline result, not a ",
    "transcription gap: it is one leg of the evidence that subcutaneous ",
    "1875 mg every 3 weeks carries a benefit-risk profile comparable to ",
    "intravenous 1200 mg every 3 weeks. Baseline C-reactive protein, by ",
    "contrast, IS significant (odds ratio 1.14 per 10 mg/L, 95 percent ",
    "CI 1.06-1.22, p < 0.001). There is no PK layer and no ODE: the ",
    "exposure metric is supplied as a data column, derived in the source ",
    "analysis from the individual empirical-Bayes predictions of the ",
    "companion population PK model packaged as Chan_2025_atezolizumab. ",
    "No between-subject random effect and no residual error are ",
    "estimated (Bernoulli likelihood). Five companion exposure-response ",
    "models in the Chan_2025_atezolizumab_* family."
  )
  reference <- paste(
    "Chan P, Liu SN, Gosselin N, Sauve Z, Marchand M, Lin A,",
    "Herraez-Baranda L, Zanghi J, Shearer-Kang E, Liu X, Wu B, Chanu P.",
    "Population pharmacokinetics and exposure-response of subcutaneous",
    "atezolizumab in patients with non-small cell lung cancer.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(4):726-737.",
    "doi:10.1002/psp4.13310.",
    "Coefficients are transcribed from Supporting Information Table S10D.",
    "Individual Cycle-1 Cmax values derive from the companion population",
    "PK model reported in the same paper; see",
    "modellib('Chan_2025_atezolizumab').",
    sep = " "
  )
  vignette <- "Chan_2025_atezolizumab_sc_nsclc"
  units <- list(
    time          = "n/a (static Cycle-1 landmark safety regression; no time dimension)",
    dosing        = "n/a (no dose events; exposure enters as the CMAX data column)",
    concentration = "prob_sae (probability of a serious adverse event, 0-1; also logit_sae)"
  )

  covariateData <- list(
    CMAX = list(
      description        = "Model-predicted Cycle-1 maximum serum atezolizumab concentration",
      units              = "ug/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters LINEARLY and is SCALED BY 100, not centred: Chan 2025",
        "Table S10D labels the row 'C_max (Increase by 100 ug /mL)', so",
        "the coefficient 0.210 is the log-odds per 100 ug/mL and the",
        "intercept -2.69 is the logit at CMAX = 0, outside the observed",
        "range. Derived from the companion population PK model, NOT",
        "observed: Chan 2025 Methods derives Cycle-1 exposure metrics",
        "'using individual empirical Bayes estimates'. Cohort 5",
        "distribution (Chan 2025 Table S5B): geometric mean 189 ug/mL",
        "(geoCV 36.2 percent), median 196, range 44.8-514. Chan 2025",
        "Methods selected 'the exposure metric with the lowest p-value'",
        "per endpoint from AUC0-21d and Cmax for the safety endpoints,",
        "and Cmax won for this one."
      ),
      source_name        = "Cmax"
    ),
    CRP = list(
      description        = "Baseline C-reactive protein",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters LINEARLY and is SCALED BY 10, not centred: Chan 2025",
        "Table S10D labels the row 'CRP (Increase by 10 mg/L)', so the",
        "coefficient 0.130 is the log-odds per 10 mg/L (odds ratio 1.14).",
        "This is the strongest predictor in the model (p < 0.001) and the",
        "Chan 2025 Discussion singles it out: 'baseline CRP level has a",
        "statistically significant impact on most of the efficacy",
        "endpoints'. It also appears in the companion AEG35, PFS and OS",
        "models. Reported in mg/L (SI)."
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
      "(NCT03735121); patients were randomized 2:1 to this subcutaneous",
      "arm versus intravenous 1200 mg every 3 weeks. Exposure-response",
      "was assessed on the subcutaneous arm ONLY -- Chan 2025 Methods",
      "notes that intravenous exposure-response has been evaluated",
      "extensively elsewhere and was therefore excluded. 'Data from all",
      "PK-evaluable patients were included in the ER analysis, and no",
      "additional exclusion was applied.' Baseline medians (Chan 2025",
      "Table 1): body weight 67.8 kg, tumor burden 79.5 mm, albumin",
      "40.0 g/L, hemoglobin 123 g/L; 29.3 percent female. Cycle-1",
      "exposure metrics were used deliberately 'to minimize the potential",
      "effect of response-dependent time-varying clearance'."
    )
  )

  ini({
    # ==================================================================
    # Chan 2025 Supporting Information Table S10D, "Logistic Regression
    # Model for SAE and Atezolizumab Exposure with Relevant Covariates in
    # Cohort 5 (Atezolizumab SC 1875 mg Q3W)".
    #
    #   Intercept                        -2.69                p < 0.01
    #   Cmax (Increase by 100 ug/mL)      0.210   OR 1.23      p = 0.441
    #                                             (0.723, 2.11)
    #   CRP (Increase by 10 mg/L)         0.130   OR 1.14      p < 0.001
    #                                             (1.06, 1.22)
    #
    # Both slopes are printed as Estimates on the log-odds scale, so
    # nothing here is back-solved or digitised. The printed odds ratios
    # are an internal check: exp(0.210) = 1.2337 rounds to 1.23 and
    # exp(0.130) = 1.1388 rounds to 1.14, both matching Table S10D.
    #
    # NEITHER covariate is centred; both are RESCALED (per 100 ug/mL and
    # per 10 mg/L). The model() block divides the raw data columns by
    # those constants so the coefficients below are used exactly as
    # printed.
    # ==================================================================
    logit_ref <- -2.69; label("Logit of the probability of a serious adverse event at CMAX = 0 ug/mL and CRP = 0 mg/L (unitless logit)")  # Chan 2025 Table S10D, Intercept row: Estimate -2.69, p < 0.01
    e_cmax_logit <- 0.210; label("Log-odds of a serious adverse event per 100 ug/mL increase in Cycle-1 Cmax (unitless logit)")           # Chan 2025 Table S10D, Cmax row: Estimate 0.210, OR 1.23 (95% CI 0.723-2.11), p = 0.441 (not significant)
    e_crp_logit <- 0.130; label("Log-odds of a serious adverse event per 10 mg/L increase in baseline C-reactive protein (unitless logit)")  # Chan 2025 Table S10D, CRP row: Estimate 0.130, OR 1.14 (95% CI 1.06-1.22), p < 0.001

    # No between-subject variability and no residual error: the source is
    # a Bernoulli-likelihood logistic regression fitted in R 4.1.1, and
    # Table S10D reports no variance components. The placeholder additive
    # term exists only because rxode2 requires an observation
    # declaration; see the vignette Assumptions and deviations.
    addSd_prob_sae <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # ------------------------------------------------------------------
    # Linear predictor. Neither covariate is centred, so logit_ref is the
    # logit at CMAX = 0 and CRP = 0, which is outside the observed data;
    # only DIFFERENCES in the linear predictor between two covariate
    # settings are interpretable.
    #
    # The divisors 100 and 10 convert the raw canonical data columns
    # (ug/mL and mg/L) into the per-100 and per-10 units in which
    # Table S10D reports its coefficients.
    # ------------------------------------------------------------------
    logit_sae <- logit_ref +
      e_cmax_logit * (CMAX / 100) +
      e_crp_logit * (CRP / 10)

    prob_sae <- expit(logit_sae)

    # ------------------------------------------------------------------
    # Observation. Deterministic probability of a serious adverse event.
    # Downstream callers can sample binary outcomes with
    # rbinom(n, 1, prob_sae) on the rxSolve output.
    # ------------------------------------------------------------------
    prob_sae ~ add(addSd_prob_sae)
  })
}
