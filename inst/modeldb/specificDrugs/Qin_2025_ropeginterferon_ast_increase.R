Qin_2025_ropeginterferon_ast_increase <- function() {
  description <- paste0(
    "Binomial logistic-regression exposure-safety model for a ",
    "treatment-related ASPARTATE AMINOTRANSFERASE INCREASE during the ",
    "DOSE-TITRATION PHASE of subcutaneous ropeginterferon alfa-2b ",
    "(ropeg) in 78 Chinese and Japanese patients with polycythaemia ",
    "vera (Qin 2025, phase II studies A19-201 and A20-202). The ",
    "probability of the event is expit(-4.75194 + 0.143564 * CAV), ",
    "where CAV is the individual average total serum ropeg ",
    "concentration over the TITRATION PHASE in ng/mL -- weeks 0-16 in ",
    "A19-201 and weeks 0-6 in A20-202, a study-dependent window that ",
    "is NOT the same column as the weeks 0-24 or 0-52 exposure used by ",
    "the efficacy models. The exposure term is significant ",
    "(p = 0.001546) and its slope is the steepest of the three ",
    "significant exposure-safety relationships. No covariate survived ",
    "screening. There is no PK layer and no ODE. No between-subject ",
    "random effect and no residual error are estimated (Bernoulli ",
    "likelihood). This is the transaminase sibling of ",
    "Qin_2025_ropeginterferon_alt_increase; companion models in the ",
    "Qin_2025_ropeginterferon_* family."
  )
  reference <- paste(
    "Qin A, Shimoda K, Suo S, Fu R, Kirito K, Wu D, Liao J, Chen H, Wu L,",
    "Su X, Gao Y, Sato T, Li Y, Zhang J, Shen W, Wang W, Zhang L, Jin J,",
    "Komatsu N.",
    "Population pharmacokinetics-pharmacodynamics and exposure-response of",
    "ropeginterferon alfa-2b in Chinese and Japanese patients with",
    "polycythemia vera.",
    "Pharmacol Res Perspect. 2025;13(3):e70109.",
    "doi:10.1002/prp2.70109.",
    sep = " "
  )
  vignette <- "Qin_2025_ropeginterferon"
  units <- list(
    time          = "n/a (static titration-phase landmark regression; no time dimension)",
    dosing        = "n/a (no dose events; the dosing history enters only through the CAV exposure column)",
    concentration = "prob_ast_increase (probability of a treatment-related AST increase during the titration phase, 0-1; also logit_ast_increase)"
  )

  covariateData <- list(
    CAV = list(
      description        = "Individual average total serum ropeginterferon alfa-2b concentration over the dose-titration phase.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "TOTAL (free plus target-bound) serum ropeg. THE AVERAGING",
        "WINDOW IS STUDY-DEPENDENT and follows each protocol's titration",
        "schedule (Qin 2025 Methods 2.4.7.1): weeks 0-16 for A19-201 and",
        "weeks 0-6 for A20-202. A data assembler must therefore compute",
        "this column per study, not with one common window. It is the",
        "SAME column as the one used by the companion ALT model, but is",
        "NOT interchangeable with the weeks 0-24 or weeks 0-52 exposure",
        "columns of the efficacy models, nor with the maintenance-phase",
        "column of the companion anemia model. Qin 2025 also matched the",
        "OUTCOME to the window: 'The incidence of safety indicators in",
        "each phase matched exposure in the corresponding phase', so an",
        "AST increase counts for this model only if it occurred during",
        "titration. Derived by simulation from the actual dosing records",
        "and the individual post hoc (empirical Bayes) PK parameters of",
        "modellib('Qin_2025_ropeginterferon'). NOT centred and NOT",
        "scaled, so the intercept is the logit at CAV = 0, an",
        "extrapolated anchor."
      ),
      source_name        = "Cavg, titration phase (average concentration within the titration phase)"
    )
  )

  covariatesDataExcluded <- list(
    AST = list(
      description = "Baseline aspartate aminotransferase activity.",
      units       = "U/L",
      type        = "continuous",
      notes       = paste(
        "Screened by forward inclusion at p < 0.05 and NOT retained:",
        "Qin 2025 Results 3.4 states 'No significant covariates were",
        "included in either model' for the ALT and AST exposure-safety",
        "regressions. Pooled baseline median 20.0 U/L, range 11.4-36.0",
        "(Table 1, Overall)."
      )
    ),
    WT = list(
      description = "Baseline body weight.",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened by forward inclusion at p < 0.05 and not retained for this endpoint."
    ),
    BMI = list(
      description = "Baseline body mass index.",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened by forward inclusion at p < 0.05 and not retained for this endpoint."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 78L,
    n_studies      = 2L,
    n_observations = "78 evaluable binary titration-phase AST-increase records, one per patient (Qin 2025 Results 3.4: 'Seventy-eight PV patients from A19-201 and A20-202 were included in the exposure-safety population')",
    age_range      = "median 54.0 years, range 26.0-72.0 (A19-201) and median 56.0 years, range 29.0-70.0 (A20-202) (Qin 2025 Table 1)",
    weight_range   = "median 56.0 kg, range 43.6-76.5 (A19-201) and median 67.9 kg, range 44.0-91.0 (A20-202) (Qin 2025 Table 1)",
    sex_female_pct = 43.6,
    race_ethnicity = "Japanese (A19-201, n = 29) and Chinese (A20-202, n = 49)",
    disease_state  = "polycythaemia vera; A20-202 enrolled patients resistant to or intolerant of hydroxyurea. Baseline AST median 22.0 U/L, range 13.0-36.0 (A19-201) and median 21.9 U/L, range 11.4-34.2 (A20-202) (Qin 2025 Table 1)",
    dose_range     = "A19-201 (slow titration): 100 ug every 2 weeks, or 50 ug on prior cytoreductive therapy, titrated in 50 ug steps to a 500 ug maximum. A20-202 (fast titration): 250 ug at week 0, 350 ug at week 2, 500 ug from week 4",
    regions        = "Japan (A19-201) and China (A20-202)",
    endpoint_definition = "Treatment-related adverse event of increased aspartate aminotransferase occurring during the study's dose-titration phase, graded by CTCAE version 5.0 (Qin 2025 Methods 2.3 and 2.4.7.1). AST increase is one of the ten most frequent treatment-related adverse events carried into the exposure-safety screen",
    notes          = paste0(
      "Model selection was evidence-led. Qin 2025 Results 3.4 reports ",
      "the exploratory analysis (Figure S7) found a significant ",
      "positive trend between titration-phase exposure and AST-increase ",
      "incidence (p = 0.009), and that 'the incidence of all three AEs ",
      "mentioned above increased linearly with increasing exposure', ",
      "which is why a LINEAR logistic form was fitted. The AST slope ",
      "(0.143564) is steeper than the ALT slope (0.1156) while the AST ",
      "intercept is lower, so the two curves cross: AST increase is the ",
      "rarer event at low exposure and the more exposure-sensitive one ",
      "across the observed range. The Discussion concludes that ",
      "transaminase elevation under rapid titration is manageable and ",
      "that 'elevated transaminase levels and anemia could serve as PD ",
      "markers that may not be associated with major adverse clinical ",
      "outcomes of ropeg treatment'."
    )
  )

  ini({
    # ==================================================================
    # Qin 2025 Results 3.4, reported in RUNNING TEXT rather than in a
    # table: "for AST increase: intercept beta [SE]: -4.75194 [1.123],
    # p < 0.0001; exposure beta [SE]: 0.143564 [0.0453], p = 0.001546."
    #
    # Note the printed precision differs between the two transaminase
    # models -- ALT is given to four decimal places and AST to five --
    # and both are transcribed here exactly as printed rather than
    # rounded to a common width.
    #
    # Fitted in R 4.2.2 (Methods 2.4.8). Methods 2.4.7.2 states
    # "Exposure-safety modeling procedures were consistent with those in
    # the development of the exposure-efficacy model", so the model form
    # is Qin 2025 Equation (3):
    #
    #   logit(P_i) = beta0 + beta1*Exposure_i + beta^T * X_i
    #
    # with the beta^T * X_i covariate block empty because no covariate
    # survived forward inclusion at p < 0.05.
    #
    # Only an estimate, a standard error and a p-value are printed --
    # no odds ratio and no confidence interval -- so there is no
    # redundant column to cross-check the coefficients against.
    #
    # The exposure regressor is NOT centred and NOT scaled.
    # ==================================================================

    # ----- Logit intercept -----
    # expit(-4.75194) = 0.00860, the extrapolated event probability at
    # zero exposure -- about half the corresponding ALT value.
    logit_ref <- -4.75194 ; label("Logit of the probability of a titration-phase AST increase at CAV = 0 ng/mL (unitless logit)")  # Qin 2025 Results 3.4 running text: intercept beta = -4.75194, standard error 1.123, p < 0.0001

    # ----- Exposure effect on the logit -----
    # exp(0.143564) = 1.154, so the odds of a titration-phase AST
    # increase rise by 15.4% per ng/mL, i.e. by
    # exp(0.143564*10) = 4.20-fold per 10 ng/mL.
    e_cav_logit <- 0.143564 ; label("Log-odds of a titration-phase AST increase per 1 ng/mL increase in average total serum ropeg concentration over the titration phase (unitless logit)")  # Qin 2025 Results 3.4 running text: exposure beta = 0.143564, standard error 0.0453, p = 0.001546 -- significant

    # ----- No between-subject variability, no residual error -----
    # The source is a binomial logistic regression; a Bernoulli
    # likelihood has no sigma, and the analysis estimates no random
    # effects. The tiny fixed additive residual below exists only so
    # rxode2 has an error model to attach to the typical-value
    # probability; it is NOT a published quantity.
    addSd_prob_ast_increase <- fixed(0.001) ; label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # ----- Linear predictor, Qin 2025 Equation (3) -----
    logit_ast_increase <- logit_ref + e_cav_logit * CAV

    prob_ast_increase <- expit(logit_ast_increase)

    # ----- Observation -----
    prob_ast_increase ~ add(addSd_prob_ast_increase)
  })
}
