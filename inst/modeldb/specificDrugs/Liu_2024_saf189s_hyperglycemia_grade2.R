Liu_2024_saf189s_hyperglycemia_grade2 <- function() {
  description <- paste0(
    "Binomial logistic-regression exposure-safety model for CTCAE ",
    "GRADE >= 2 hyperglycemia in Chinese adults with ALK-positive or ",
    "ROS1-positive advanced non-small cell lung cancer treated with the ",
    "second-generation ALK/ROS1 tyrosine kinase inhibitor SAF-189s ",
    "(Liu 2024, n = 296, the safety analysis set of study SAF001, ",
    "NCT04237805, 20-210 mg orally once daily). The probability of the ",
    "event is expit(-16.812 + 2.036 * log(AUC_SAF189S)) where ",
    "AUC_SAF189S is the individual SAF-189s steady-state daily AUC in ",
    "ng*h/mL computed at the patient's FIRST dose level and the log is ",
    "a natural log, so the printed odds ratio 7.662 is per e-fold. The ",
    "exposure dependence is roughly twice as steep as for the ",
    "any-grade companion model (2.036 versus 1.259 on the logit), so ",
    "exposure discriminates SEVERITY of hyperglycemia more sharply ",
    "than it discriminates occurrence; this is the quantitative basis ",
    "for the paper's conclusion that the 210 mg dose group was less ",
    "tolerated than the lower-dose groups while 160 mg once daily was ",
    "well tolerated. There is no PK layer and no ODE: the exposure ",
    "metric is supplied as a data column, derived in the source ",
    "analysis from the individual post hoc parameters of the companion ",
    "population PK model, packaged as Liu_2024_saf189s. No ",
    "between-subject random effect and no residual error are estimated ",
    "(Bernoulli likelihood). The intercept is not printed by the paper ",
    "and was recovered by digitising the fitted curve of Figure 10B; ",
    "see the in-file note on logit_ref and the vignette Errata. ",
    "Companion any-grade model in Liu_2024_saf189s_hyperglycemia; six ",
    "companion models in the Liu_2024_saf189s_* family."
  )
  reference <- paste(
    "Liu Y, Tan Y, Hu L, Li J, Yang J, Diao L, Yang J.",
    "Population pharmacokinetics and exposure-response analyses of",
    "SAF-189s in Chinese patients with ALK+/ROS1+ non-small cell lung",
    "cancer.",
    "Front Pharmacol. 2024;15:1418549.",
    "doi:10.3389/fphar.2024.1418549.",
    "Individual SAF-189s steady-state AUC values derive from the companion",
    "population pharmacokinetic model reported in the same paper; see",
    "modellib('Liu_2024_saf189s').",
    sep = " "
  )
  vignette <- "Liu_2024_saf189s"
  units <- list(
    time          = "n/a (static landmark exposure-safety regression; no time dimension)",
    dosing        = "n/a (no dose events; exposure enters as the AUC_SAF189S covariate column)",
    concentration = "prob_hyperglycemia_grade2 (probability of grade >= 2 hyperglycemia, 0-1; also logit_hyperglycemia_grade2)"
  )

  covariateData <- list(
    AUC_SAF189S = list(
      description        = "Individual SAF-189s area under the plasma concentration-time curve over the 24 h dosing interval at steady state (AUCss), per subject. Supplied as data: this model has no PK layer, and the source analysis used the individual post hoc parameters of the companion population PK model together with the patient's FIRST dose level (Liu 2024 Methods, E-R analysis).",
      units              = "ng*h/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "TOTAL (not unbound) exposure, at STEADY STATE, computed at the",
        "FIRST dose the patient received -- stated by the Figure 10 axis",
        "label ('AUCss (ng*h/mL) ... by first dose'). The companion",
        "exposure-EFFICACY model Liu_2024_saf189s_orr uses a different",
        "metric on a different basis (Cmin,ss at the most prevalent",
        "dose), so the two exposure columns are not interchangeable.",
        "Enters on the NATURAL LOG scale: digitising the Figure 10B",
        "fitted curve returns a slope of 2.0376 per natural-log unit",
        "against the printed odds ratio 7.662 (exp(2.0376) = 7.672),",
        "with the logit linear in log(AUCss) to R^2 = 0.99997. The",
        "paper's statement that the odds ratio corresponds to 'an",
        "increase in SAF-189s exposure of 1 ng*h*mL-1*d-1' is an",
        "erratum; the unit is one natural-log unit.",
        "Units are load-bearing: the intercept absorbs the unit choice,",
        "so supplying ug*h/mL instead of ng*h/mL would shift the logit",
        "by 2.036 * log(1000) = 14.07 -- more than twice the shift in",
        "the any-grade companion model, because the slope is steeper.",
        "For calibration, the analysis set's median AUCss is",
        "2,233 ng*h/mL and the Figure 10 x-axis spans",
        "174.0-8,338.1 ng*h/mL."
      ),
      source_name        = "AUCss estimated by the first dose"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 296L,
    n_studies      = 1L,
    n_observations = "296 binary event records (one per patient; first-occurrence landmark analysis, no repeated measures)",
    age_range      = "SAF001: phase I median 51 years (28.0-68.0), phase II median 54.1 years (20.0-84.0) (Liu 2024 Supplementary Table 1)",
    weight_range   = "median 63.2 kg, range 37.3-92.5 across the PK-evaluable cohort",
    sex_female_pct = 51.6,
    race_ethnicity = c(Asian = 100),
    disease_state  = "ALK-positive or ROS1-positive advanced non-small cell lung cancer; hyperglycemia graded by CTCAE",
    dose_range     = "SAF-189s 20, 40, 80, 120, 160 or 210 mg orally once daily in 21-day cycles",
    regions        = "China",
    notes          = paste0(
      "Same 296-patient safety analysis set as the any-grade companion ",
      "model. Observed outcome: 82 of 296 (27.70%) had grade >= 2 ",
      "hyperglycemia against 165 (55.74%) at any grade (Liu 2024 ",
      "Table 5), so roughly half of all hyperglycemia events reached ",
      "grade 2 or worse."
    )
  )

  ini({
    # ==================================================================
    # As for the any-grade companion model, Liu 2024 reports the SLOPE
    # but not the intercept: Results, Exposure-safety analysis gives
    # only "the ORs ... for ... grade >= 2 hyperglycemia ... were ...
    # 7.662 (95% CI: 3.523-16.661, p < 0.001)". No coefficient table
    # exists for any of the three exposure-safety fits.
    #
    # The intercept was recovered by DIGITISING the fitted blue curve
    # of Figure 10B. Method and cross-checks:
    #
    #   - Axis calibration used the two unambiguous outer tick labels,
    #     174.0126 and 8338.1020 ng*h/mL, shared by all three panels of
    #     Figure 10 (the panels' tick spacings agree to half a pixel).
    #     Validation: the figure's median-quartile dashed line lands at
    #     2,231 ng*h/mL against the 2,233 ng*h/mL printed in the
    #     Discussion.
    #   - The curve is linear in logit versus log(AUCss) to
    #     R^2 = 0.99997 over 325 digitised pixel columns.
    #   - The FREE-fit slope is 2.0376 per natural-log unit, i.e. an
    #     odds ratio of exp(2.0376) = 7.672 against the printed 7.662 --
    #     a 0.13% agreement.
    #   - Holding the slope at the printed log(7.662) and refitting the
    #     intercept alone gives -16.812.
    #   - INDEPENDENT CORROBORATION: at the analysis set's median AUCss
    #     of 2,233 ng*h/mL this pair predicts P(event) = 0.248, against
    #     the observed grade >= 2 incidence of 82/296 = 27.70% counted
    #     in Table 5. The 2.9-percentage-point gap is the largest of the
    #     three exposure-safety models and reflects that the fitted
    #     probability at the median of a right-skewed exposure
    #     distribution sits slightly below the mean fitted probability,
    #     which is the quantity a logistic MLE matches to the observed
    #     proportion.
    #
    # The slope below is the PRINTED value; only the intercept is
    # figure-derived. See the vignette Errata.
    # ==================================================================

    # ----- Logit intercept -----
    # digitised from Liu 2024 Figure 10B (fitted curve); NOT printed in the paper.
    # See the block comment above for the method and the three cross-checks.
    logit_ref <- -16.812 ; label("Logit of the probability of grade >= 2 hyperglycemia at AUC_SAF189S = 1 ng*h/mL (unitless logit)")  # figure-derived: Liu 2024 Figure 10B fitted curve, intercept recovered with the slope held at the printed log(7.662); cross-checks P(event) = 0.248 at the median 2,233 ng*h/mL vs the observed 27.70% in Table 5

    # ----- Exposure slope on the logit -----
    e_auc_saf189s_logit <- log(7.662) ; label("Log-odds of grade >= 2 hyperglycemia per e-fold increase in SAF-189s steady-state daily AUC (unitless logit)")  # Liu 2024 Results, Exposure-safety analysis: OR 7.662 (95% CI 3.523-16.661), p < 0.001; free-fit digitisation of Figure 10B independently returns OR 7.672

    # ----- No between-subject variability, no residual error -----
    # Bernoulli-likelihood logistic regression; no omega or sigma is
    # estimated. The tiny placeholder additive residual exists only
    # because rxode2 requires an observation declaration; see the
    # vignette Assumptions and deviations section.
    addSd_prob_hyperglycemia_grade2 <- fixed(0.001) ; label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # ------------------------------------------------------------------
    # Linear predictor. The exposure is NOT centred, so logit_ref is the
    # logit at AUC_SAF189S = 1 ng*h/mL.
    # ------------------------------------------------------------------
    logit_hyperglycemia_grade2 <- logit_ref + e_auc_saf189s_logit * log(AUC_SAF189S)

    prob_hyperglycemia_grade2 <- expit(logit_hyperglycemia_grade2)

    # ------------------------------------------------------------------
    # Observation. Deterministic probability of grade >= 2
    # hyperglycemia. Downstream callers can sample binary outcomes with
    # rbinom(n, 1, prob_hyperglycemia_grade2) on the rxSolve output.
    # ------------------------------------------------------------------
    prob_hyperglycemia_grade2 ~ add(addSd_prob_hyperglycemia_grade2)
  })
}
