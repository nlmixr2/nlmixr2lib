Qin_2025_ropeginterferon_anemia <- function() {
  description <- paste0(
    "Binomial logistic-regression exposure-safety model for ANEMIA ",
    "during the DOSE-MAINTENANCE PHASE of subcutaneous ropeginterferon ",
    "alfa-2b (ropeg) in 78 Chinese and Japanese patients with ",
    "polycythaemia vera (Qin 2025, phase II studies A19-201 and ",
    "A20-202). The probability of the event is expit(-6.95416 + ",
    "0.056548 * CAV + 0.041351 * TUM_JAK2_V617F_VAF), where CAV is the ",
    "individual average total serum ropeg concentration over the ",
    "MAINTENANCE PHASE in ng/mL -- weeks 16 to end of treatment in ",
    "A19-201 and weeks 6 to end of treatment in A20-202, a ",
    "study-dependent window that is NOT the same column as the ",
    "titration-phase exposure used by the transaminase safety models ",
    "nor the weeks 0-24 or 0-52 exposure used by the efficacy models. ",
    "Both the exposure term (p = 0.0284) and the baseline JAK2 V617F ",
    "allele-burden covariate (p = 0.0183) are significant; this is the ",
    "only one of the three significant exposure-safety endpoints in ",
    "which a covariate survived forward inclusion. Anemia risk rises ",
    "with exposure and with the size of the driver-mutant clone at ",
    "baseline, which is mechanistically coherent: ropeg acts by ",
    "depleting that clone. There is no PK layer and no ODE. No ",
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
    time          = "n/a (static maintenance-phase landmark regression; no time dimension)",
    dosing        = "n/a (no dose events; the dosing history enters only through the CAV exposure column)",
    concentration = "prob_anemia (probability of a treatment-related anemia event during the maintenance phase, 0-1; also logit_anemia)"
  )

  covariateData <- list(
    CAV = list(
      description        = "Individual average total serum ropeginterferon alfa-2b concentration over the dose-maintenance phase.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "TOTAL (free plus target-bound) serum ropeg. THE AVERAGING",
        "WINDOW IS STUDY-DEPENDENT and is the COMPLEMENT of the",
        "titration-phase window used by the companion transaminase",
        "models (Qin 2025 Methods 2.4.7.1): 'For A19-201, the titration",
        "phase was defined as the start of treatment to Week 16, and the",
        "maintenance phase ranged from Week 16 to the end of treatment.",
        "For A20-202, the titration phase was from the start of treatment",
        "to Week 6, and the maintenance phase was from Week 6 to the end",
        "of treatment.' A data assembler must therefore compute this",
        "column per study, not with one common window. Maintenance-phase",
        "exposure is systematically HIGHER than titration-phase exposure",
        "because patients have reached their titrated dose, so this",
        "column is not interchangeable with the titration-phase column",
        "of Qin_2025_ropeginterferon_alt_increase even though both are",
        "named CAV, and neither is interchangeable with the weeks 0-24",
        "or weeks 0-52 columns of the efficacy models. Qin 2025 matched",
        "the OUTCOME to the window: 'The incidence of safety indicators",
        "in each phase matched exposure in the corresponding phase', so",
        "an anemia event counts for this model only if it occurred",
        "during maintenance. Derived by simulation from the actual",
        "dosing records and the individual post hoc (empirical Bayes) PK",
        "parameters of modellib('Qin_2025_ropeginterferon'). NOT centred",
        "and NOT scaled, so the intercept is the logit at CAV = 0 and",
        "zero allele burden, a doubly extrapolated anchor."
      ),
      source_name        = "Cavg, maintenance phase (average concentration within the maintenance phase)"
    ),
    TUM_JAK2_V617F_VAF = list(
      description        = "Baseline JAK2 V617F somatic variant allele frequency (allele burden) in peripheral blood, before the first ropeg dose.",
      units              = "percent (0-100)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "THE SCALE IS PER CENT, NOT A FRACTION. Qin 2025 Table 1 prints",
        "the row as 'JAK2 V617F (%)' with medians of 77.8 (range",
        "0.0210-97.5) in A19-201 and 61.2 (range 4.70-96.4) in A20-202;",
        "the 48 phase I healthy volunteers are 0 by construction, which",
        "is why the Overall column reads 31.9 (range 0-97.5) and must",
        "NOT be used as the patient-population summary. Enters the logit",
        "LINEARLY and UNCENTRED, so the coefficient below is per",
        "PERCENTAGE POINT and the intercept is the logit at a burden of",
        "zero -- an extrapolated anchor, since all A20-202 patients and",
        "all but two A19-201 patients carried the mutation (Qin 2025",
        "Methods 2.1). A model fitted on the 0-1 fraction scale would",
        "need the coefficient multiplied by 100. This is a SOMATIC clone",
        "fraction, not a germline genotype: it is the proportion of",
        "hematopoiesis carried by the driver-mutant clone, it moves",
        "under treatment, and the companion",
        "Qin_2025_ropeginterferon_jak2_week24 and _week52 models predict",
        "that movement. Only the BASELINE value is used here.",
        "Quantitated as described in Qin 2025 references 39, 42 and 44",
        "(Methods 2.2)."
      ),
      source_name        = "JAK2 V617F (%) (baseline JAK2 V617F allele burden)"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Baseline body weight.",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened by forward inclusion at p < 0.05 and not retained for",
        "this endpoint, although it WAS retained by the companion",
        "week-24 JAK2 V617F efficacy regression",
        "(Qin_2025_ropeginterferon_jak2_week24)."
      )
    ),
    BMI = list(
      description = "Baseline body mass index.",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened by forward inclusion at p < 0.05 and not retained for this endpoint."
    ),
    HCT = list(
      description = "Baseline hematocrit.",
      units       = "percent (volume fraction times 100), per the HCT register entry. Note that Qin 2025 Table 1 mixes the two scales inside a single cell, printing an A20-202 median of 45.7 against a range of [0.421, 64.1]",
      type        = "continuous",
      notes       = paste(
        "Part of the Table 1 covariate set screened for the",
        "exposure-response models (Qin 2025 Methods 2.4.6.2: 'Covariates",
        "tested were the same as those in the PopPK analysis'). Not",
        "retained. Of direct interest for an anemia endpoint -- a lower",
        "baseline red-cell mass would plausibly predispose to anemia --",
        "so its absence from the final model is informative rather than",
        "incidental."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 78L,
    n_studies      = 2L,
    n_observations = "78 evaluable binary maintenance-phase anemia records, one per patient (Qin 2025 Results 3.4: 'Seventy-eight PV patients from A19-201 and A20-202 were included in the exposure-safety population')",
    age_range      = "median 54.0 years, range 26.0-72.0 (A19-201) and median 56.0 years, range 29.0-70.0 (A20-202) (Qin 2025 Table 1)",
    weight_range   = "median 56.0 kg, range 43.6-76.5 (A19-201) and median 67.9 kg, range 44.0-91.0 (A20-202) (Qin 2025 Table 1)",
    sex_female_pct = 43.6,
    race_ethnicity = "Japanese (A19-201, n = 29) and Chinese (A20-202, n = 49)",
    disease_state  = "polycythaemia vera; A20-202 enrolled patients resistant to or intolerant of hydroxyurea. Baseline JAK2 V617F allele burden median 77.8%, range 0.0210-97.5 (A19-201) and median 61.2%, range 4.70-96.4 (A20-202) (Qin 2025 Table 1). All A20-202 patients and all but two A19-201 patients carried the mutation",
    dose_range     = "A19-201 (slow titration): 100 ug every 2 weeks, or 50 ug on prior cytoreductive therapy, titrated in 50 ug steps to a 500 ug maximum. A20-202 (fast titration): 250 ug at week 0, 350 ug at week 2, 500 ug from week 4",
    regions        = "Japan (A19-201) and China (A20-202)",
    endpoint_definition = "Treatment-related adverse event of anemia occurring during the study's dose-maintenance phase, graded by CTCAE version 5.0 (Qin 2025 Methods 2.3 and 2.4.7.1). Anemia is one of the ten most frequent treatment-related adverse events carried into the exposure-safety screen",
    notes          = paste0(
      "Model selection was evidence-led. Qin 2025 Results 3.4 reports ",
      "that 'In the maintenance phase, anemia occurrence was associated ",
      "with average exposure in this phase (p = 0.022)' in the ",
      "exploratory analysis (Figure S7), and that 'the incidence of all ",
      "three AEs mentioned above increased linearly with increasing ",
      "exposure', which is why a LINEAR logistic form was fitted. Only ",
      "three of the ten screened adverse events reached a significant ",
      "exposure-response relationship: ALT increase and AST increase in ",
      "the titration phase, and anemia in the maintenance phase. This ",
      "is the only one of the three that retained a covariate. The ",
      "Discussion puts the finding in clinical context: 'In the ",
      "dose-maintenance phase, the incidence of anemia increased with ",
      "increasing exposure, suggesting a risk of anemia in participants ",
      "with higher exposure under long-term administration. However, no ",
      "grade 3 or clinically significant anemia was observed in our ",
      "clinical Phase II studies.'"
    )
  )

  ini({
    # ==================================================================
    # Qin 2025 Results 3.4, reported in RUNNING TEXT rather than in a
    # table: "Logistic regression model for Cavg, maintenance phase and
    # anemia indicated a statistically significant E-R relationship for
    # this AE (intercept beta [SE]: -6.95416 [1.9606], p = 0.00039;
    # exposure beta [SE]: 0.056548 [0.0258], p = 0.0284). Inclusion of
    # baseline JAK2 V617F identified as a significant covariate effect
    # (JAK2 V617F beta [SE]: 0.041351 [0.0175], p = 0.0183)."
    #
    # Fitted in R 4.2.2 (Methods 2.4.8). Methods 2.4.7.2 states
    # "Exposure-safety modeling procedures were consistent with those in
    # the development of the exposure-efficacy model", so the model form
    # is Qin 2025 Equation (3):
    #
    #   logit(P_i) = beta0 + beta1*Exposure_i + beta^T * X_i
    #
    # with a single retained covariate in the beta^T * X_i block.
    # Methods 2.4.6.2 governs how it got there: "if exposure was a
    # significant predictor ... a covariate model was evaluated ...
    # using forward inclusion, with the p-value set at 0.05."
    #
    # Only an estimate, a standard error and a p-value are printed --
    # no odds ratio and no confidence interval -- so there is no
    # redundant column to cross-check the coefficients against. Each
    # value below is the printed estimate with its standard error in the
    # comment. Note that the intercept and the covariate coefficient are
    # printed to SIX significant figures while their standard errors
    # carry three; the extra digits are reproduced verbatim rather than
    # rounded, since rounding a logit intercept is not neutral.
    #
    # BOTH regressors are UNCENTRED and UNSCALED.
    # ==================================================================

    # ----- Logit intercept -----
    # expit(-6.95416) = 0.000953, the extrapolated event probability for
    # a patient with zero maintenance-phase exposure AND zero baseline
    # allele burden. Neither condition occurs in the analysis
    # population, so this is an anchor and not a predicted placebo rate.
    logit_ref <- -6.95416 ; label("Logit of the probability of a maintenance-phase anemia event at CAV = 0 ng/mL and a baseline JAK2 V617F allele burden of 0 percent (unitless logit)")  # Qin 2025 Results 3.4 running text: intercept beta = -6.95416, standard error 1.9606, p = 0.00039

    # ----- Exposure effect on the logit -----
    # exp(0.056548) = 1.058, so the odds of a maintenance-phase anemia
    # event rise by 5.8% per ng/mL, i.e. by exp(0.056548*10) = 1.76-fold
    # per 10 ng/mL. This is the shallowest of the three significant
    # exposure-safety slopes in the paper (compare 0.1156 for ALT and
    # 0.143564 for AST), consistent with anemia being the endpoint the
    # Discussion describes as never reaching grade 3. The exposure
    # metric is not the same column as the transaminase models', so the
    # slopes are comparable only as per-ng/mL sensitivities.
    e_cav_logit <- 0.056548 ; label("Log-odds of a maintenance-phase anemia event per 1 ng/mL increase in average total serum ropeg concentration over the maintenance phase (unitless logit)")  # Qin 2025 Results 3.4 running text: exposure beta = 0.056548, standard error 0.0258, p = 0.0284 -- significant

    # ----- Baseline allele-burden effect on the logit -----
    # exp(0.041351) = 1.0422, so the odds rise by 4.22% per PERCENTAGE
    # POINT of baseline allele burden, i.e. by exp(0.041351*10) = 1.51
    # per 10 points. Across the observed patient span this is a large
    # effect: the A19-201 median of 77.8% carries exp(0.041351*77.8) =
    # 25-fold higher odds than a burden of zero, and the A19-201-versus-
    # A20-202 median gap of 16.6 points alone is a 1.99-fold difference.
    # Qin 2025 Figure 5D draws the fitted surface at the first quartile,
    # median and third quartile of baseline JAK2 V617F.
    e_jak2vaf_logit <- 0.041351 ; label("Log-odds of a maintenance-phase anemia event per 1 percentage point increase in baseline JAK2 V617F allele burden (unitless logit)")  # Qin 2025 Results 3.4 running text: JAK2 V617F beta = 0.041351, standard error 0.0175, p = 0.0183 -- significant

    # ----- No between-subject variability, no residual error -----
    # The source is a binomial logistic regression; a Bernoulli
    # likelihood has no sigma, and the analysis estimates no random
    # effects. The tiny fixed additive residual below exists only so
    # rxode2 has an error model to attach to the typical-value
    # probability; it is NOT a published quantity.
    addSd_prob_anemia <- fixed(0.001) ; label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # ----- Linear predictor, Qin 2025 Equation (3) -----
    # TUM_JAK2_V617F_VAF is in PER CENT (0-100), matching the printed
    # coefficient. Supplying a 0-1 fraction here understates the term
    # 100-fold.
    logit_anemia <-
      logit_ref +
      e_cav_logit * CAV +
      e_jak2vaf_logit * TUM_JAK2_V617F_VAF

    prob_anemia <- expit(logit_anemia)

    # ----- Observation -----
    prob_anemia ~ add(addSd_prob_anemia)
  })
}
