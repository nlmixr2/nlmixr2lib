Liu_2024_saf189s_proteinuria <- function() {
  description <- paste0(
    "Binomial logistic-regression exposure-safety model for ANY-GRADE ",
    "proteinuria in Chinese adults with ALK-positive or ROS1-positive ",
    "advanced non-small cell lung cancer treated with the ",
    "second-generation ALK/ROS1 tyrosine kinase inhibitor SAF-189s ",
    "(Liu 2024, n = 296, the safety analysis set of study SAF001, ",
    "NCT04237805, 20-210 mg orally once daily). The probability of the ",
    "event is expit(-6.297 + 0.709 * log(AUC_SAF189S)) where ",
    "AUC_SAF189S is the individual SAF-189s steady-state daily AUC in ",
    "ng*h/mL computed at the patient's FIRST dose level and the log is ",
    "a natural log, so the printed odds ratio 2.031 is per e-fold. ",
    "This is the shallowest and least certain of the paper's three ",
    "significant exposure-safety relationships (95% CI 1.112-3.711, ",
    "p = 0.021, against p < 0.001 for both hyperglycemia endpoints). ",
    "Liu 2024 found no exposure dependence for the remaining adverse ",
    "events of interest -- hypercholesterolemia, nausea, vomiting and ",
    "diarrhea -- and reports no coefficients for them, so no companion ",
    "models exist for those endpoints. There is no PK layer and no ",
    "ODE: the exposure metric is supplied as a data column, derived in ",
    "the source analysis from the individual post hoc parameters of the ",
    "companion population PK model, packaged as Liu_2024_saf189s. No ",
    "between-subject random effect and no residual error are estimated ",
    "(Bernoulli likelihood). The intercept is not printed by the paper ",
    "and was recovered by digitising the fitted curve of Figure 10C; ",
    "see the in-file note on logit_ref and the vignette Errata. Six ",
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
    concentration = "prob_proteinuria (probability of any-grade proteinuria, 0-1; also logit_proteinuria)"
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
        "Enters on the NATURAL LOG scale: digitising the Figure 10C",
        "fitted curve returns a slope of 0.7088 per natural-log unit",
        "against the printed odds ratio 2.031 (exp(0.7088) = 2.032),",
        "with the logit linear in log(AUCss) to R^2 = 0.99997. The",
        "paper's statement that the odds ratio corresponds to 'an",
        "increase in SAF-189s exposure of 1 ng*h*mL-1*d-1' is an",
        "erratum; the unit is one natural-log unit.",
        "Units are load-bearing: the intercept absorbs the unit choice,",
        "so supplying ug*h/mL instead of ng*h/mL would shift the logit",
        "by 0.709 * log(1000) = 4.90.",
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
    disease_state  = "ALK-positive or ROS1-positive advanced non-small cell lung cancer; proteinuria graded by CTCAE",
    dose_range     = "SAF-189s 20, 40, 80, 120, 160 or 210 mg orally once daily in 21-day cycles",
    regions        = "China",
    notes          = paste0(
      "Same 296-patient safety analysis set as the two hyperglycemia ",
      "companion models. Observed outcome: 90 of 296 (30.41%) had ",
      "any-grade proteinuria and 28 (9.46%) had grade >= 2 (Liu 2024 ",
      "Table 5). Only the any-grade endpoint was modelled against ",
      "exposure; no grade >= 2 proteinuria coefficient is reported."
    )
  )

  ini({
    # ==================================================================
    # As for the two hyperglycemia companion models, Liu 2024 reports
    # the SLOPE but not the intercept: Results, Exposure-safety analysis
    # gives only "the ORs ... for ... proteinuria ... were ... 2.031
    # (95% CI: 1.112-3.711, p = 0.021)". No coefficient table exists
    # for any of the three exposure-safety fits.
    #
    # The intercept was recovered by DIGITISING the fitted blue curve
    # of Figure 10C. Method and cross-checks:
    #
    #   - Axis calibration used the two unambiguous outer tick labels,
    #     174.0126 and 8338.1020 ng*h/mL, shared by all three panels of
    #     Figure 10. Validation: the figure's median-quartile dashed
    #     line lands at 2,231 ng*h/mL against the 2,233 ng*h/mL printed
    #     in the Discussion.
    #   - The curve is linear in logit versus log(AUCss) to
    #     R^2 = 0.99997 over 469 digitised pixel columns.
    #   - The FREE-fit slope is 0.7088 per natural-log unit, i.e. an
    #     odds ratio of exp(0.7088) = 2.032 against the printed 2.031 --
    #     a 0.05% agreement.
    #   - Holding the slope at the printed log(2.031) and refitting the
    #     intercept alone gives -6.297.
    #   - INDEPENDENT CORROBORATION, and the tightest of the three: at
    #     the analysis set's median AUCss of 2,233 ng*h/mL this pair
    #     predicts P(event) = 0.303, against the observed any-grade
    #     proteinuria incidence of 90/296 = 30.41% counted in Table 5.
    #     The two agree to a tenth of a percentage point.
    #
    # The slope below is the PRINTED value; only the intercept is
    # figure-derived. See the vignette Errata.
    # ==================================================================

    # ----- Logit intercept -----
    # digitised from Liu 2024 Figure 10C (fitted curve); NOT printed in the paper.
    # See the block comment above for the method and the three cross-checks.
    logit_ref <- -6.297 ; label("Logit of the probability of any-grade proteinuria at AUC_SAF189S = 1 ng*h/mL (unitless logit)")  # figure-derived: Liu 2024 Figure 10C fitted curve, intercept recovered with the slope held at the printed log(2.031); cross-checks P(event) = 0.303 at the median 2,233 ng*h/mL vs the observed 30.41% in Table 5

    # ----- Exposure slope on the logit -----
    e_auc_saf189s_logit <- log(2.031) ; label("Log-odds of any-grade proteinuria per e-fold increase in SAF-189s steady-state daily AUC (unitless logit)")  # Liu 2024 Results, Exposure-safety analysis: OR 2.031 (95% CI 1.112-3.711), p = 0.021; free-fit digitisation of Figure 10C independently returns OR 2.032

    # ----- No between-subject variability, no residual error -----
    # Bernoulli-likelihood logistic regression; no omega or sigma is
    # estimated. The tiny placeholder additive residual exists only
    # because rxode2 requires an observation declaration; see the
    # vignette Assumptions and deviations section.
    addSd_prob_proteinuria <- fixed(0.001) ; label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # ------------------------------------------------------------------
    # Linear predictor. The exposure is NOT centred, so logit_ref is the
    # logit at AUC_SAF189S = 1 ng*h/mL.
    # ------------------------------------------------------------------
    logit_proteinuria <- logit_ref + e_auc_saf189s_logit * log(AUC_SAF189S)

    prob_proteinuria <- expit(logit_proteinuria)

    # ------------------------------------------------------------------
    # Observation. Deterministic probability of any-grade proteinuria.
    # Downstream callers can sample binary outcomes with
    # rbinom(n, 1, prob_proteinuria) on the rxSolve output.
    # ------------------------------------------------------------------
    prob_proteinuria ~ add(addSd_prob_proteinuria)
  })
}
