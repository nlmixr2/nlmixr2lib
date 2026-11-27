Liu_2026_venetoclax_orr <- function() {
  description <- paste0(
    "Univariate binomial logistic-regression exposure-efficacy model for the ",
    "probability of overall response (complete or partial remission on the ",
    "bone-marrow assessment after one 28-day venetoclax cycle) in adult ",
    "Chinese patients with acute myeloid leukemia receiving venetoclax 100-400 ",
    "mg once daily with chemotherapy in routine practice (Liu 2026, n = 63, ",
    "37 responders). The probability is expit of a linear predictor in the ",
    "individual steady-state average venetoclax concentration Cave, with a ",
    "positive slope (p = 0.0239, AIC 81.81). The paper publishes this ",
    "relationship only as a fitted curve (Figure 4B); the intercept and slope ",
    "were recovered by digitising that curve. No PK layer and no ODE: Cave is ",
    "supplied as a data column, derived in the source analysis from the ",
    "post hoc parameters of the companion population PK model ",
    "Liu_2026_venetoclax. No random effect and no residual error are ",
    "estimated (Bernoulli likelihood)."
  )
  reference <- paste0(
    "Liu Z, Liu X, Gong Q, Qin S, Zhu X, Kuan IHS, Mak WY, Xiang X, Jia C, ",
    "Wang Q, Cheng L, Li L. Real-World Pharmacokinetic and Exposure-Response ",
    "Characterization of Venetoclax in Chinese Patients with Hematological ",
    "Malignancies. Drug Des Devel Ther. 2026;20:563629. ",
    "doi:10.2147/DDDT.S563629. Coefficients digitised from the fitted curve ",
    "of Figure 4B."
  )
  vignette <- "Liu_2026_venetoclax"
  units <- list(
    time = "n/a (static landmark exposure-response regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "prob_orr_investigator (probability of CR or PR after one cycle, 0-1; also logit_orr)"
  )

  covariateData <- list(
    CAV = list(
      description = "Individual steady-state average venetoclax plasma concentration (Cave). Supplied as data: this model has no PK layer.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Liu 2026 (Methods, Exposure-Efficacy Relationship Evaluation)",
        "derived Cave from the individual post hoc PK parameters of the",
        "companion population PK model (Liu_2026_venetoclax) and the actual",
        "dose, i.e. Cave = daily dose / (24 h x individual CL/F) x 1000 for",
        "once-daily dosing. Enters linearly and uncentred; the slope is per",
        "1 ng/mL. The analysis-set Cave spanned about 330-7330 ng/mL",
        "(Figure 4A/4B); extrapolation outside that range is not supported",
        "by the source."
      ),
      source_name = "Cave"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 63L,
    n_studies = 1L,
    n_observations = "63 binary response records, one per AML patient (37 responders, 26 non-responders; Results and Figure 4)",
    race_ethnicity = c(Asian = 100),
    disease_state = paste0(
      "Acute myeloid leukemia, the AML subset of the 76-patient PK cohort; ",
      "every patient received venetoclax with one or two chemotherapy agents ",
      "(azacitidine alone in 41 of 63; Table S1)."
    ),
    dose_range = "Venetoclax 100 mg (n = 7), 200 mg (n = 13), 300 mg (n = 2) or 400 mg (n = 41) once daily (Table S1)",
    regions = "Single centre, Southwest Hospital, Chongqing, China (July 2021 to March 2024)",
    response_definition = paste0(
      "Overall response = complete remission (CR) or partial remission (PR) ",
      "versus no response or relapse, on a bone-marrow biopsy after one ",
      "28-day cycle (typically day 28-45), judged by the treating physician ",
      "against Chinese adult-AML guideline criteria (Methods). CR with ",
      "incomplete count recovery is not a separate category."
    ),
    notes = paste0(
      "Logistic regression with Cave as the single predictor (Cmin was also ",
      "evaluated; Cave was the metric the authors carried forward). ",
      "A sensitivity analysis in all 76 patients (Figures S8 and S9) and the ",
      "Cmin analysis in the AML subset (Figure S7) showed the same trend and ",
      "are not encoded here."
    )
  )

  ini({
    # ================================================================
    # Liu 2026 Figure 4B: the fitted ('Predicted') probability-of-efficacy
    # curve versus Cave in the 63 AML patients, annotated AIC = 81.81 and
    # p = 0.0239. Neither coefficient is printed anywhere in the paper or
    # its supplement.
    #
    #   logit(p) = logit_ref + e_cav_logit * CAV
    #
    # Digitised by the maintainers from the publisher's native-resolution
    # Figure 4 raster: x and y calibrated from the axis ticks, the orange
    # curve tracked column by column (2283 points over Cave 338-7323 ng/mL),
    # and logit(p) regressed on untransformed Cave (R^2 = 0.99992 against
    # 0.891 for log Cave, so the exposure enters linearly). Cross-checked
    # against the 63 individual Cave values digitised from Figure 4A, whose
    # maximum-likelihood refit reproduces AIC 81.82 and Wald p = 0.0240;
    # see the vignette.
    # ================================================================

    logit_ref <- -0.9055; label("Logit of the overall-response probability at zero venetoclax Cave (unitless logit)")  # Figure 4B fitted curve, digitised; intercept -0.905465 (not printed)

    e_cav_logit <- 0.0006443; label("Log-odds of overall response per 1 ng/mL increase in venetoclax Cave (unitless logit per ng/mL)")  # Figure 4B fitted curve, digitised; slope 0.000644271 per ng/mL (not printed; p = 0.0239 printed in the panel)

    # Bernoulli likelihood in the source; the fixed additive residual is a
    # placeholder so rxode2 has an error model, not a published quantity.
    addSd_prob_orr_investigator <- fixed(0.001); label("Placeholder additive residual SD on the response probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_orr <- logit_ref + e_cav_logit * CAV

    prob_orr_investigator <- expit(logit_orr)

    prob_orr_investigator ~ add(addSd_prob_orr_investigator)
  })
}
