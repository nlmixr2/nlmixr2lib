Liu_2024_saf189s_orr <- function() {
  description <- paste0(
    "Binomial logistic-regression exposure-efficacy model for the ",
    "independent-review-committee-assessed overall response rate (ORR, ",
    "complete or partial response by RECIST 1.1) in Chinese adults ",
    "with ALK-positive or ROS1-positive advanced non-small cell lung ",
    "cancer treated with the second-generation ALK/ROS1 tyrosine kinase ",
    "inhibitor SAF-189s (Liu 2024, n = 244, the phase II analysis set ",
    "of study SAF001, NCT04237805, 80-210 mg orally once daily). The ",
    "probability of response is expit(3.871 - 0.587 * log(CTROUGH)) ",
    "where CTROUGH is the individual SAF-189s steady-state trough ",
    "concentration in ng/mL at the patient's most prevalent dose and ",
    "the log is a natural log, so the printed odds ratio 0.556 is per ",
    "e-fold. The relationship is NOT statistically significant (95% CI ",
    "0.285-1.083, p = 0.084) and the paper's conclusion is that ",
    "SAF-189s exposures over 80-210 mg sit on the PLATEAU of the ",
    "exposure-response curve for efficacy; the model is packaged so ",
    "that this null result is reproducible rather than merely asserted. ",
    "There is no PK layer and no ODE: the exposure metric is supplied ",
    "as a data column, derived in the source analysis from the ",
    "individual post hoc parameters of the companion population PK ",
    "model, packaged as Liu_2024_saf189s. No between-subject random ",
    "effect and no residual error are estimated (Bernoulli ",
    "likelihood). The intercept is not printed by the paper and was ",
    "recovered by digitising the fitted curve of Figure 6; see the ",
    "in-file note on logit_ref and the vignette Errata. Six companion ",
    "exposure-response models in the Liu_2024_saf189s_* family."
  )
  reference <- paste(
    "Liu Y, Tan Y, Hu L, Li J, Yang J, Diao L, Yang J.",
    "Population pharmacokinetics and exposure-response analyses of",
    "SAF-189s in Chinese patients with ALK+/ROS1+ non-small cell lung",
    "cancer.",
    "Front Pharmacol. 2024;15:1418549.",
    "doi:10.3389/fphar.2024.1418549.",
    "Individual SAF-189s trough concentrations derive from the companion",
    "population pharmacokinetic model reported in the same paper; see",
    "modellib('Liu_2024_saf189s').",
    sep = " "
  )
  vignette <- "Liu_2024_saf189s"
  units <- list(
    time          = "n/a (static landmark exposure-efficacy regression; no time dimension)",
    dosing        = "n/a (no dose events; exposure enters as the CTROUGH covariate column)",
    concentration = "prob_orr_central (probability of an IRC-assessed complete or partial response, 0-1; also logit_orr_central)"
  )

  covariateData <- list(
    CTROUGH = list(
      description        = "Individual SAF-189s steady-state trough plasma concentration (Cmin,ss), per subject. Supplied as data: this model has no PK layer, and the source analysis used the individual post hoc parameters of the companion population PK model together with the patient's most prevalent dose level (Liu 2024 Methods, E-R analysis).",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "TOTAL (not unbound) plasma concentration, at STEADY STATE, and",
        "computed at the MOST PREVALENT dose the patient received --",
        "not at the first dose. That distinction is load-bearing in this",
        "paper: the companion exposure-SAFETY models use AUCss computed",
        "at the FIRST dose instead, so the two exposure columns are not",
        "interchangeable even after a unit conversion.",
        "Enters on the NATURAL LOG scale. Liu 2024 Figure 6 plots the",
        "fitted curve against 'C_min,ss (ng/mL) in Logarithmic scale',",
        "and digitising that curve returns a slope of -0.5883 per",
        "natural-log unit against the printed odds ratio 0.556",
        "(exp(-0.5883) = 0.5553), with the logit linear in log(CTROUGH)",
        "to R^2 = 0.99995 -- so the log is natural and the odds ratio is",
        "per e-fold, not per ng/mL.",
        "Units are load-bearing: the intercept absorbs the unit choice,",
        "so supplying ug/mL instead of ng/mL would shift the logit by",
        "-0.587 * log(1000) = -4.06.",
        "For calibration, Liu 2024 Table 3 gives the analysis set's",
        "Cmin,ss quartile boundaries as 18, 60.7, 79, 108 and",
        "182 ng/mL (61 patients per quartile)."
      ),
      source_name        = "C min,ss (steady-state minimum concentration) estimated using the most prevalent dose"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 244L,
    n_studies      = 1L,
    n_observations = "244 binary response records (one per patient; landmark analysis, no repeated measures)",
    age_range      = "SAF001 phase II: median 54.1 years, range 20.0-84.0 (Liu 2024 Supplementary Table 1)",
    weight_range   = "SAF001 phase II pooled with the rest of the PK-evaluable cohort: median 63.2 kg, range 37.3-92.5",
    sex_female_pct = 53.6,
    race_ethnicity = c(Asian = 100),
    disease_state  = "ALK-positive or ROS1-positive advanced non-small cell lung cancer; response assessed by an independent review committee under RECIST 1.1, ORR = complete response or partial response",
    dose_range     = "SAF-189s 80, 120, 160 or 210 mg orally once daily in 21-day cycles (phase IIa) or 160 mg once daily (phase IIb)",
    regions        = "China",
    notes          = paste0(
      "The exposure-efficacy analysis set is the 244 phase II patients ",
      "with IRC assessment, a subset of the 317 PK-evaluable subjects ",
      "used to fit the companion population PK model (Liu 2024 Results, ",
      "Exposure-response analysis). Observed outcome: overall ORR ",
      "78.69% (95% CI 73.92%-82.93%), i.e. roughly 192 responders and ",
      "52 non-responders. Best overall response by Cmin,ss quartile is ",
      "given in Liu 2024 Table 3 and is essentially flat: 88.5%, ",
      "72.1%, 77.0% and 77.0% complete-or-partial response from the ",
      "lowest to the highest quartile."
    )
  )

  ini({
    # ==================================================================
    # Liu 2024 reports the SLOPE of this regression but not its
    # intercept: Results, Exposure-response analysis gives only
    # "logistic regression analyses did not show a significant
    # relationship between the probability of achieving ORR and Cmin,ss
    # in all patients (odds ratio (OR): 0.556; 95% CI: 0.285-1.083;
    # p = 0.084)". There is no coefficient table for this fit anywhere
    # in the paper or in either supplementary document.
    #
    # The intercept was therefore recovered by DIGITISING the fitted
    # blue curve of Figure 6, which the paper draws over the full
    # observed Cmin,ss range on a log x-axis. Method and cross-checks:
    #
    #   - Axis calibration used the two unambiguous outer tick labels,
    #     18.01681 and 182.24820 ng/mL.
    #   - The curve is linear in logit versus log(Cmin,ss) to
    #     R^2 = 0.99995 over 1,395 digitised pixel columns, which is
    #     what a logistic fit on a log-transformed regressor must look
    #     like and confirms the functional form.
    #   - The FREE-fit slope is -0.5883 per natural-log unit, i.e. an
    #     odds ratio of exp(-0.5883) = 0.5553 against the printed
    #     0.556 -- a 0.1% agreement. Recovering the published slope to
    #     three significant figures from the figure alone is what
    #     licenses reading the intercept off the same curve.
    #   - Holding the slope at the printed log(0.556) and refitting the
    #     intercept alone gives 3.871.
    #   - INDEPENDENT CORROBORATION, from a source that is not the
    #     figure: at the analysis set's median Cmin,ss of 79 ng/mL
    #     (Liu 2024 Table 3 quartile boundary) this pair predicts
    #     P(response) = 0.787, against the paper's separately reported
    #     observed ORR of 78.69% for the same n = 244 set. A logistic
    #     MLE reproduces the observed event proportion at the centre of
    #     the exposure distribution, so this is the check the intercept
    #     has to pass, and it passes to a tenth of a percentage point.
    #
    # The slope below is the PRINTED value; only the intercept is
    # figure-derived. See the vignette Errata.
    # ==================================================================

    # ----- Logit intercept -----
    # digitised from Liu 2024 Figure 6 (fitted curve); NOT printed in the paper.
    # See the block comment above for the method and the two cross-checks.
    logit_ref <- 3.871 ; label("Logit of the probability of an IRC-assessed response at CTROUGH = 1 ng/mL (unitless logit)")  # figure-derived: Liu 2024 Figure 6 fitted curve, intercept recovered with the slope held at the printed log(0.556); cross-checks P(resp) = 0.787 at the median 79 ng/mL vs the observed ORR 78.69%

    # ----- Exposure slope on the logit -----
    e_ctrough_logit <- log(0.556) ; label("Log-odds of an IRC-assessed response per e-fold increase in SAF-189s steady-state trough concentration (unitless logit)")  # Liu 2024 Results, Exposure-response analysis: OR 0.556 (95% CI 0.285-1.083), p = 0.084; free-fit digitisation of Figure 6 independently returns OR 0.5553

    # ----- No between-subject variability, no residual error -----
    # The source model is a Bernoulli-likelihood logistic regression:
    # the exposure fully determines each subject's response probability,
    # and no omega or sigma is estimated. rxode2 requires an observation
    # declaration, so the deterministic probability is emitted with a
    # tiny placeholder additive residual, mirroring the
    # Chen_2021_lorlatinib_teae_grade3.R and
    # Fukae_2024_valemetostat_orr_central.R pattern. This does not
    # perturb the predicted probability; see the vignette Assumptions
    # and deviations section.
    addSd_prob_orr_central <- fixed(0.001) ; label("Placeholder additive residual SD on the typical-value response probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # ------------------------------------------------------------------
    # Linear predictor. The exposure is NOT centred: Liu 2024 fits the
    # raw log-transformed trough, so logit_ref is the logit at
    # CTROUGH = 1 ng/mL rather than at any clinically meaningful
    # patient. The negative slope means the fitted response probability
    # DECREASES with exposure, which is biologically implausible as a
    # causal statement and is exactly why the paper reports the result
    # as null (the confidence interval spans 1); the most likely driver
    # is confounding by dose reduction and by the treatment-resistant
    # subgroup, not a real inverse exposure-response.
    # ------------------------------------------------------------------
    logit_orr_central <- logit_ref + e_ctrough_logit * log(CTROUGH)

    prob_orr_central <- expit(logit_orr_central)

    # ------------------------------------------------------------------
    # Observation. Deterministic probability of an IRC-assessed complete
    # or partial response. Downstream callers can sample binary outcomes
    # with rbinom(n, 1, prob_orr_central) on the rxSolve output.
    # ------------------------------------------------------------------
    prob_orr_central ~ add(addSd_prob_orr_central)
  })
}
