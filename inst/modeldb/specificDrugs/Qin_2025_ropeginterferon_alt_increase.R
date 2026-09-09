Qin_2025_ropeginterferon_alt_increase <- function() {
  description <- paste0(
    "Binomial logistic-regression exposure-safety model for a ",
    "treatment-related ALANINE AMINOTRANSFERASE INCREASE during the ",
    "DOSE-TITRATION PHASE of subcutaneous ropeginterferon alfa-2b ",
    "(ropeg) in 78 Chinese and Japanese patients with polycythaemia ",
    "vera (Qin 2025, phase II studies A19-201 and A20-202). The ",
    "probability of the event is expit(-4.099 + 0.1156 * CAV), where ",
    "CAV is the individual average total serum ropeg concentration ",
    "over the TITRATION PHASE in ng/mL -- weeks 0-16 in A19-201 and ",
    "weeks 0-6 in A20-202, a study-dependent window that is NOT the ",
    "same column as the weeks 0-24 or 0-52 exposure used by the ",
    "efficacy models. The exposure term is significant ",
    "(p = 0.000378). No covariate survived screening. Qin 2025 judges ",
    "the resulting risk manageable: only one grade 3 ALT increase ",
    "occurred, under fast titration, and all others were mild or ",
    "moderate and reversible. There is no PK layer and no ODE. No ",
    "between-subject random effect and no residual error are estimated ",
    "(Bernoulli likelihood). Companion models in the ",
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
    concentration = "prob_alt_increase (probability of a treatment-related ALT increase during the titration phase, 0-1; also logit_alt_increase)"
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
        "schedule (Qin 2025 Methods 2.4.7.1): 'For A19-201, the titration",
        "phase was defined as the start of treatment to Week 16, and the",
        "maintenance phase ranged from Week 16 to the end of treatment.",
        "For A20-202, the titration phase was from the start of treatment",
        "to Week 6, and the maintenance phase was from Week 6 to the end",
        "of treatment.' A data assembler must therefore compute this",
        "column per study, not with one common window. It is NOT",
        "interchangeable with the weeks 0-24 or weeks 0-52 exposure",
        "columns of the efficacy models, nor with the maintenance-phase",
        "column of the companion anemia model. Qin 2025 also matched the",
        "OUTCOME to the window: 'The incidence of safety indicators in",
        "each phase matched exposure in the corresponding phase', so an",
        "ALT increase counts for this model only if it occurred during",
        "titration. Derived by simulation from the actual dosing records",
        "and the individual post hoc (empirical Bayes) PK parameters of",
        "modellib('Qin_2025_ropeginterferon'). NOT centred and NOT",
        "scaled, so the intercept is the logit at CAV = 0, an",
        "extrapolated anchor. Qin 2025 Discussion notes the titration",
        "phase exposure ranges of the two regimens overlap: 'the exposure",
        "range for the fast-dose titration regimen was similar to that",
        "for slow titration'."
      ),
      source_name        = "Cavg, titration phase (average concentration within the titration phase)"
    )
  )

  covariatesDataExcluded <- list(
    ALT = list(
      description = "Baseline alanine aminotransferase activity.",
      units       = "U/L",
      type        = "continuous",
      notes       = paste(
        "Screened by forward inclusion at p < 0.05 and NOT retained:",
        "Qin 2025 Results 3.4 states 'No significant covariates were",
        "included in either model' for the ALT and AST exposure-safety",
        "regressions. Pooled baseline median 17.1 U/L, range 7.00-51.0",
        "(Table 1, Overall). Note the contrast with the companion anemia",
        "model, which DID retain baseline JAK2 V617F allele burden, so",
        "the covariate screen was capable of finding an effect here and",
        "did not."
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
    n_observations = "78 evaluable binary titration-phase ALT-increase records, one per patient (Qin 2025 Results 3.4: 'Seventy-eight PV patients from A19-201 and A20-202 were included in the exposure-safety population')",
    age_range      = "median 54.0 years, range 26.0-72.0 (A19-201) and median 56.0 years, range 29.0-70.0 (A20-202) (Qin 2025 Table 1)",
    weight_range   = "median 56.0 kg, range 43.6-76.5 (A19-201) and median 67.9 kg, range 44.0-91.0 (A20-202) (Qin 2025 Table 1)",
    sex_female_pct = 43.6,
    race_ethnicity = "Japanese (A19-201, n = 29) and Chinese (A20-202, n = 49)",
    disease_state  = "polycythaemia vera; A20-202 enrolled patients resistant to or intolerant of hydroxyurea. Baseline ALT median 15.0 U/L, range 8.00-46.0 (A19-201) and median 20.5 U/L, range 7.00-51.0 (A20-202) (Qin 2025 Table 1)",
    dose_range     = "A19-201 (slow titration): 100 ug every 2 weeks, or 50 ug on prior cytoreductive therapy, titrated in 50 ug steps to a 500 ug maximum. A20-202 (fast titration): 250 ug at week 0, 350 ug at week 2, 500 ug from week 4",
    regions        = "Japan (A19-201) and China (A20-202)",
    endpoint_definition = "Treatment-related adverse event of increased alanine aminotransferase occurring during the study's dose-titration phase, graded by CTCAE version 5.0 (Qin 2025 Methods 2.3 and 2.4.7.1). ALT increase is one of the ten most frequent treatment-related adverse events carried into the exposure-safety screen",
    notes          = paste0(
      "Model selection was evidence-led. Qin 2025 Results 3.4 reports ",
      "the exploratory analysis (Figure S7) found 'a significant ",
      "positive trend between exposure in the titration phase and ",
      "incidence (p = 0.021 and p = 0.009, respectively)' for ALT and ",
      "AST, and that 'the incidence of all three AEs mentioned above ",
      "increased linearly with increasing exposure', which is why a ",
      "LINEAR logistic form was fitted. Only three of the ten screened ",
      "adverse events reached a significant exposure-response ",
      "relationship: ALT increase and AST increase in the titration ",
      "phase, and anemia in the maintenance phase. The Discussion ",
      "concludes 'the risk of transaminase elevation caused by rapid ",
      "titration is manageable'."
    )
  )

  ini({
    # ==================================================================
    # Qin 2025 Results 3.4, reported in RUNNING TEXT rather than in a
    # table: "Logistic regression models for Cavg, titration phase and
    # increases in ALT and AST indicated a statistically significant E-R
    # relationship: for ALT increase: intercept beta [Standard error
    # (SE)]: -4.099 [0.9843], p < 0.0001; exposure beta [SE]: 0.1156
    # [0.0399], p = 0.000378."
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
    # redundant column to cross-check the coefficients against. Each
    # value below is the printed estimate with its standard error in the
    # comment.
    #
    # The exposure regressor is NOT centred and NOT scaled.
    # ==================================================================

    # ----- Logit intercept -----
    # expit(-4.099) = 0.0163, the extrapolated event probability at zero
    # exposure.
    logit_ref <- -4.099 ; label("Logit of the probability of a titration-phase ALT increase at CAV = 0 ng/mL (unitless logit)")  # Qin 2025 Results 3.4 running text: intercept beta = -4.099, standard error 0.9843, p < 0.0001

    # ----- Exposure effect on the logit -----
    # exp(0.1156) = 1.123, so the odds of a titration-phase ALT increase
    # rise by 12.3% per ng/mL, i.e. by exp(0.1156*10) = 3.18-fold per
    # 10 ng/mL. This is the steepest exposure-safety slope of the three
    # significant safety endpoints after AST.
    e_cav_logit <- 0.1156 ; label("Log-odds of a titration-phase ALT increase per 1 ng/mL increase in average total serum ropeg concentration over the titration phase (unitless logit)")  # Qin 2025 Results 3.4 running text: exposure beta = 0.1156, standard error 0.0399, p = 0.000378 -- significant

    # ----- No between-subject variability, no residual error -----
    # The source is a binomial logistic regression; a Bernoulli
    # likelihood has no sigma, and the analysis estimates no random
    # effects. The tiny fixed additive residual below exists only so
    # rxode2 has an error model to attach to the typical-value
    # probability; it is NOT a published quantity.
    addSd_prob_alt_increase <- fixed(0.001) ; label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # ----- Linear predictor, Qin 2025 Equation (3) -----
    logit_alt_increase <- logit_ref + e_cav_logit * CAV

    prob_alt_increase <- expit(logit_alt_increase)

    # ----- Observation -----
    prob_alt_increase ~ add(addSd_prob_alt_increase)
  })
}
