Qin_2025_ropeginterferon_jak2_week52 <- function() {
  description <- paste0(
    "Linear-regression exposure-efficacy model for the change from ",
    "baseline in JAK2 V617F allele burden at WEEK 52 of subcutaneous ",
    "ropeginterferon alfa-2b (ropeg) in 72 Chinese and Japanese ",
    "patients with polycythaemia vera (Qin 2025, phase II studies ",
    "A19-201 and A20-202). The predicted change is ",
    "-8.528 - 0.43 * CAV percentage points, where CAV is the ",
    "individual average total serum ropeg concentration over weeks ",
    "0-52 in ng/mL. The exposure term remains SIGNIFICANT ",
    "(p = 0.0295) with a slope only slightly shallower than at week 24 ",
    "(-0.43 versus -0.51), which is the paper's evidence that ropeg's ",
    "disease-modifying effect on the driver clone does NOT plateau -- ",
    "in pointed contrast to the complete-hematologic-response ",
    "endpoint, which had gone flat by this landmark. NO covariate was ",
    "retained here, so body weight, which is significant at week 24, ",
    "is absent. The outcome is a CHANGE FROM BASELINE in percentage ",
    "points, so negative values are the desired direction. There is no ",
    "PK layer and no ODE, and no between-subject random effect is ",
    "estimated. Companion models in the Qin_2025_ropeginterferon_* ",
    "family."
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
    time          = "n/a (static week-52 landmark regression; no time dimension)",
    dosing        = "n/a (no dose events; the dosing history enters only through the CAV exposure column)",
    concentration = "djak2v617f (change from baseline in JAK2 V617F allele burden at week 52, percentage points; negative = reduction)"
  )

  covariateData <- list(
    CAV = list(
      description        = "Individual average total serum ropeginterferon alfa-2b concentration over weeks 0 to 52.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "TOTAL (free plus target-bound) serum ropeg. Averaging window is",
        "weeks 0-52, from the FIRST DOSE to the landmark, so it includes",
        "the titration period and is NOT a steady-state interval average",
        "(Qin 2025 Methods 2.4.6.1). This is a DIFFERENT column from the",
        "week-24 companion model's CAV and the two must not be",
        "interchanged. Derived by simulation from the actual dosing",
        "records and the individual post hoc (empirical Bayes) PK",
        "parameters of modellib('Qin_2025_ropeginterferon'). NOT centred",
        "and NOT scaled, so the intercept is the predicted change at",
        "CAV = 0. Observed range about 10-70 ng/mL (Qin 2025 Figure 3D),",
        "with study medians near 28 ng/mL (A19-201) and 42 ng/mL",
        "(A20-202) -- a narrower separation than at week 24, because the",
        "longer averaging window dilutes the titration difference."
      ),
      source_name        = "Cavg,0-52W (average concentration of participants from 0 to 52 weeks)"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Baseline body weight.",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened by forward inclusion at p < 0.05 and NOT retained at",
        "this landmark: Qin 2025 Results 3.3 states 'No significant",
        "covariates were identified in the model at Week 52.' Weight IS",
        "significant at week 24 (p = 0.0119) in the companion",
        "Qin_2025_ropeginterferon_jak2_week24 model, so the omission",
        "here is the paper's own result and the correct encoding is the",
        "absence of the term rather than a fixed(0) coefficient."
      )
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
    n_subjects     = 72L,
    n_studies      = 2L,
    n_observations = "72 evaluable JAK2 V617F allele-burden change records, one per patient (Qin 2025 Results 3.3: 'For JAK2 V617F allele burden, 74 and 72 patients underwent measurements at Weeks 24 and 52, respectively')",
    age_range      = "median 54.0 years, range 26.0-72.0 (A19-201) and median 56.0 years, range 29.0-70.0 (A20-202) (Qin 2025 Table 1)",
    weight_range   = "median 56.0 kg, range 43.6-76.5 (A19-201) and median 67.9 kg, range 44.0-91.0 (A20-202) (Qin 2025 Table 1)",
    sex_female_pct = 43.6,
    race_ethnicity = "Japanese (A19-201, n = 29) and Chinese (A20-202, n = 49)",
    disease_state  = "polycythaemia vera. Baseline JAK2 V617F allele burden median 77.8%, range 0.0210-97.5 (A19-201) and median 61.2%, range 4.70-96.4 (A20-202) (Qin 2025 Table 1)",
    dose_range     = "A19-201 (slow titration): 100 ug every 2 weeks, or 50 ug on prior cytoreductive therapy, titrated in 50 ug steps to a 500 ug maximum. A20-202 (fast titration): 250 ug at week 0, 350 ug at week 2, 500 ug from week 4",
    regions        = "Japan (A19-201) and China (A20-202)",
    endpoint_definition = "Change from baseline in JAK2 V617F allele burden (percentage points) at week 52; the molecular response was the secondary efficacy endpoint of both phase II studies (Qin 2025 Methods 2.3)",
    notes          = paste0(
      "Qin 2025 Results 3.3 reports the exploratory scatterplot ",
      "(Figure S6) showed a significant trend at this landmark too ",
      "(p = 0.0295 by F-test) with no plateau, which is why a LINEAR ",
      "model was fitted. The Discussion is explicit that this ",
      "persistence is the point: 'we found that the average ",
      "concentration remained significantly correlated with the ",
      "reduction in JAK2 V617F allele burden from baseline at Week 52, ",
      "despite the slight flattening of the exposure slope ... The ",
      "modeling results reflect the need to maintain sufficient PK ",
      "exposure of ropeg for its disease-modifying action.' Note the ",
      "intercept is NEGATIVE here (-8.528) against +30.728 at week 24, ",
      "but the two are not comparable because the week-24 model also ",
      "carries an uncentred body-weight term that shifts its intercept ",
      "by about -27 percentage points at a typical weight."
    )
  )

  ini({
    # ==================================================================
    # Qin 2025 Table 4, row block "Linear regression of JAK2 V617F at
    # Week 52". Fitted in R 4.2.2 (Methods 2.4.8).
    #
    # Model form is Qin 2025 Equation (5):
    #
    #   JAK2_reduction_from_baseline = beta0 + beta1*Exposure_i
    #                                        + beta^T * X_i
    #
    # with the beta^T * X_i covariate block EMPTY -- no covariate
    # survived forward inclusion at p < 0.05 at this landmark.
    #
    # CROSS-CHECK AGAINST FIGURE 3D. The figure plots the fitted line
    # over the observed exposure range. At CAV = 0 the model predicts
    # -8.528 and the figure's line begins at about -8; at CAV = 70 it
    # predicts -8.528 - 0.43*70 = -38.6 and the figure reads about -38.
    # The uncentred reading reproduces the published line.
    #
    # Table 4 prints an Estimate, a standard error and a p-value per row
    # and nothing else, so each value below is the Estimate with the
    # standard error in the comment.
    # ==================================================================

    # ----- Regression intercept -----
    # The predicted week-52 change at CAV = 0 ng/mL. Not significant on
    # its own (p = 0.2787), and an extrapolation below the observed
    # exposure range of roughly 10-70 ng/mL, so it is meaningful only
    # in combination with the slope.
    intercept <- -8.528 ; label("Predicted change from baseline in JAK2 V617F allele burden at week 52 at CAV = 0 ng/mL (percentage points; an extrapolated anchor)")  # Qin 2025 Table 4: intercept beta0 = -8.528, standard error 7.812, p = 0.2787

    # ----- Exposure effect -----
    # Slightly shallower than the week-24 slope of -0.51, which is the
    # "slight flattening of the exposure slope" the Discussion notes,
    # but still significant.
    e_cav_jak2 <- -0.43 ; label("Change in the week-52 JAK2 V617F allele-burden response per 1 ng/mL increase in average total serum ropeg concentration over weeks 0-52 (percentage points per ng/mL)")  # Qin 2025 Table 4: exposure effect beta1 = -0.43, standard error 0.194, p = 0.0295 -- significant

    # ----- Residual error -----
    # Qin 2025 does not print a residual standard deviation, a residual
    # sum of squares or an R^2 for either JAK2 regression, so the
    # residual magnitude is NOT recoverable from the paper. The small
    # fixed additive residual below exists only so rxode2 has an error
    # model to attach to the typical-value prediction; it is NOT a
    # published quantity and must not be used to characterise
    # prediction uncertainty.
    addSd_djak2v617f <- fixed(0.001) ; label("Placeholder additive residual SD on the typical-value predicted allele-burden change; not published by the source")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # ----- Linear predictor, Qin 2025 Equation (5) -----
    djak2v617f <- intercept + e_cav_jak2 * CAV

    # ----- Observation -----
    djak2v617f ~ add(addSd_djak2v617f)
  })
}
