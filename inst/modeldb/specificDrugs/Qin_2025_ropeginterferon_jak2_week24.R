Qin_2025_ropeginterferon_jak2_week24 <- function() {
  description <- paste0(
    "Linear-regression exposure-efficacy model for the change from ",
    "baseline in JAK2 V617F allele burden at WEEK 24 of subcutaneous ",
    "ropeginterferon alfa-2b (ropeg) in 74 Chinese and Japanese ",
    "patients with polycythaemia vera (Qin 2025, phase II studies ",
    "A19-201 and A20-202). The predicted change is ",
    "30.728 - 0.51 * CAV - 0.435 * WT percentage points, where CAV is ",
    "the individual average total serum ropeg concentration over weeks ",
    "0-24 in ng/mL and WT is baseline body weight in kg. BOTH ",
    "regressors are significant (exposure p = 0.0006, weight ",
    "p = 0.0119), so greater exposure AND greater body weight each ",
    "predict a larger reduction in the driver-mutation allele burden. ",
    "This is the paper's disease-modification result and the only one ",
    "of its four exposure-efficacy models to retain a covariate. The ",
    "outcome is a CHANGE FROM BASELINE in percentage points, so ",
    "negative values are the desired direction. There is no PK layer ",
    "and no ODE, and no between-subject random effect is estimated. ",
    "Companion models in the Qin_2025_ropeginterferon_* family."
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
    time          = "n/a (static week-24 landmark regression; no time dimension)",
    dosing        = "n/a (no dose events; the dosing history enters only through the CAV exposure column)",
    concentration = "djak2v617f (change from baseline in JAK2 V617F allele burden at week 24, percentage points; negative = reduction)"
  )

  covariateData <- list(
    CAV = list(
      description        = "Individual average total serum ropeginterferon alfa-2b concentration over weeks 0 to 24.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "TOTAL (free plus target-bound) serum ropeg. Averaging window is",
        "weeks 0-24, from the FIRST DOSE to the landmark, not a",
        "steady-state dosing interval (Qin 2025 Methods 2.4.6.1).",
        "Identical column to the one used by the week-24 CHR companion",
        "model. Derived by simulation from the actual dosing records and",
        "the individual post hoc (empirical Bayes) PK parameters of",
        "modellib('Qin_2025_ropeginterferon'). NOT centred and NOT",
        "scaled. Observed range about 8-64 ng/mL (Qin 2025 Figure 3B),",
        "with study medians near 20 ng/mL (A19-201, slow titration) and",
        "37.5 ng/mL (A20-202, fast titration)."
      ),
      source_name        = "Cavg,0-24W (average concentration of participants from 0 to 24 weeks)"
    ),
    WT = list(
      description        = "Baseline body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters the regression LINEARLY and UNCENTRED, as a plain",
        "additive term -0.435 * WT, not as an allometric or",
        "median-normalised covariate -- Equation (5) writes the covariate",
        "block as beta^T * X_i with no centring, and the printed",
        "intercept of 30.728 is only reconcilable with the Figure 3C",
        "curves on the uncentred reading (see the ini() block for that",
        "arithmetic check). The coefficient is therefore per KILOGRAM,",
        "and the intercept is the predicted change at WT = 0 kg, an",
        "extrapolated anchor with no clinical meaning on its own.",
        "Figure 3C draws the fitted surface at three weights labelled",
        "56.05, 62.5 and 72.08 kg, which are the first quartile, median",
        "and third quartile of body weight in the 74-patient week-24",
        "analysis set. Note the sign: HEAVIER patients are predicted to",
        "have a LARGER allele-burden reduction at any given exposure,",
        "which cannot be a disguised exposure effect because exposure is",
        "already in the model and heavier patients clear ropeg faster.",
        "Body weight was screened but NOT retained in the companion",
        "population PK model (which kept BMI on clearance) nor in the",
        "week-52 JAK2 companion model."
      ),
      source_name        = "baseline body weight"
    )
  )

  covariatesDataExcluded <- list(
    BMI = list(
      description = "Baseline body mass index.",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened by forward inclusion at p < 0.05 and not retained for this endpoint, even though it is the covariate retained on clearance in the companion population PK model. Body weight entered instead."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 74L,
    n_studies      = 2L,
    n_observations = "74 evaluable JAK2 V617F allele-burden change records, one per patient (Qin 2025 Results 3.3: 'For JAK2 V617F allele burden, 74 and 72 patients underwent measurements at Weeks 24 and 52, respectively')",
    age_range      = "median 54.0 years, range 26.0-72.0 (A19-201) and median 56.0 years, range 29.0-70.0 (A20-202) (Qin 2025 Table 1)",
    weight_range   = "median 62.5 kg with first and third quartiles 56.05 and 72.08 kg in this analysis set, read from the Qin 2025 Figure 3C legend; the pooled study-level medians are 56.0 kg (A19-201) and 67.9 kg (A20-202) (Table 1)",
    sex_female_pct = 43.6,
    race_ethnicity = "Japanese (A19-201, n = 29) and Chinese (A20-202, n = 49)",
    disease_state  = "polycythaemia vera. Baseline JAK2 V617F allele burden median 77.8%, range 0.0210-97.5 (A19-201) and median 61.2%, range 4.70-96.4 (A20-202) (Qin 2025 Table 1). All A20-202 patients and all but two A19-201 patients carried the mutation",
    dose_range     = "A19-201 (slow titration): 100 ug every 2 weeks, or 50 ug on prior cytoreductive therapy, titrated in 50 ug steps to a 500 ug maximum. A20-202 (fast titration): 250 ug at week 0, 350 ug at week 2, 500 ug from week 4",
    regions        = "Japan (A19-201) and China (A20-202)",
    endpoint_definition = "Change from baseline in JAK2 V617F allele burden (percentage points) at week 24; the molecular response was the secondary efficacy endpoint of both phase II studies (Qin 2025 Methods 2.3)",
    notes          = paste0(
      "Model selection was evidence-led. Qin 2025 Results 3.3 reports ",
      "the exploratory scatterplot (Figure S6) showed a significant ",
      "positive linear trend between exposure and allele-burden ",
      "reduction (p = 0.0004 by F-test) with 'No efficacy plateau ",
      "observed', which is why a LINEAR rather than an Emax model was ",
      "fitted. The Discussion draws the clinical inference: 'higher ",
      "ropeg exposure in patients is efficient in eradicating ",
      "neoplastic clones carrying JAK2 V617F, and persistent ",
      "efficacious exposure may help deplete neoplastic clones that ",
      "drive the PV phenotype and progression.' Unlike the CHR ",
      "endpoint, the exposure relationship for allele burden persists ",
      "at week 52 (see the companion ",
      "Qin_2025_ropeginterferon_jak2_week52)."
    )
  )

  ini({
    # ==================================================================
    # Qin 2025 Table 4, row block "Linear regression of JAK2 V617F at
    # Week 24". Fitted in R 4.2.2 (Methods 2.4.8).
    #
    # Model form is Qin 2025 Equation (5):
    #
    #   JAK2_reduction_from_baseline = beta0 + beta1*Exposure_i
    #                                        + beta^T * X_i
    #
    # with the covariate block carrying baseline body weight.
    #
    # ARITHMETIC CHECK OF THE UNCENTRED READING. Equation (5) does not
    # print a centring constant, and Figure 3C lets that be verified
    # rather than assumed. The figure draws the fitted line at three
    # body weights, 56.05, 62.5 and 72.08 kg. On the UNCENTRED reading
    # the predicted change at CAV = 0 is
    #     56.05 kg: 30.728 - 0.435*56.05 = +6.35 percentage points
    #     62.50 kg: 30.728 - 0.435*62.50 = +3.54
    #     72.08 kg: 30.728 - 0.435*72.08 = -0.63
    # and Figure 3C's three intercepts at x = 0 read approximately +6,
    # +3 and -1 in that order. At CAV = 60 ng/mL the 56.05 kg line
    # predicts 6.35 - 0.51*60 = -24.3 against a figure value of about
    # -25. The uncentred reading reproduces the published figure; any
    # centred reading would shift the intercept by 0.435*COV_med, i.e.
    # by about 27 percentage points, and cannot.
    #
    # Table 4 prints an Estimate, a standard error and a p-value per row
    # and nothing else, so each value below is the Estimate with the
    # standard error in the comment.
    # ==================================================================

    # ----- Regression intercept -----
    # The predicted week-24 change at CAV = 0 ng/mL AND WT = 0 kg. Both
    # regressors are uncentred, so this is a double extrapolation and is
    # not interpretable on its own; it is meaningful only in combination
    # with the two slopes below.
    intercept <- 30.728 ; label("Predicted change from baseline in JAK2 V617F allele burden at week 24 at CAV = 0 ng/mL and WT = 0 kg (percentage points; an extrapolated anchor)")  # Qin 2025 Table 4: intercept beta0 = 30.728, standard error 11.601, p = 0.01

    # ----- Exposure effect -----
    # Across the roughly 20 ng/mL gap between the slow- and
    # fast-titration study medians this predicts an additional
    # 0.51*20 = 10.2 percentage points of allele-burden reduction.
    e_cav_jak2 <- -0.51 ; label("Change in the week-24 JAK2 V617F allele-burden response per 1 ng/mL increase in average total serum ropeg concentration over weeks 0-24 (percentage points per ng/mL)")  # Qin 2025 Table 4: exposure effect beta1 = -0.51, standard error 0.141, p = 0.0006 -- significant

    # ----- Body-weight effect -----
    e_wt_jak2 <- -0.435 ; label("Change in the week-24 JAK2 V617F allele-burden response per 1 kg increase in baseline body weight (percentage points per kg)")  # Qin 2025 Table 4: weight effect beta2 = -0.435, standard error 0.169, p = 0.0119 -- significant

    # ----- Residual error -----
    # Qin 2025 does not print a residual standard deviation, a residual
    # sum of squares or an R^2 for either JAK2 regression, so the
    # residual magnitude is NOT recoverable from the paper. The small
    # fixed additive residual below exists only so rxode2 has an error
    # model to attach to the typical-value prediction; it is NOT a
    # published quantity and must not be used to characterise
    # prediction uncertainty. See the vignette's Assumptions and
    # deviations.
    addSd_djak2v617f <- fixed(0.001) ; label("Placeholder additive residual SD on the typical-value predicted allele-burden change; not published by the source")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # ----- Linear predictor, Qin 2025 Equation (5) -----
    # Both regressors enter uncentred and unscaled; see the arithmetic
    # check against Figure 3C in the ini() block above.
    djak2v617f <- intercept + e_cav_jak2 * CAV + e_wt_jak2 * WT

    # ----- Observation -----
    djak2v617f ~ add(addSd_djak2v617f)
  })
}
