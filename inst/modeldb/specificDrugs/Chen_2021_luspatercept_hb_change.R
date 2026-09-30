Chen_2021_luspatercept_hb_change <- function() {
  description <- paste0(
    "Linear exposure-response model relating the steady-state luspatercept AUC of the starting dose to the average change from baseline ",
    "in hemoglobin over the first 3-week dosing interval in adults with beta-thalassemia and a low transfusion burden ",
    "(< 12 RBC units/24 weeks) (Chen 2021, Figure 2A; n = 34 phase 2 patients with weekly hemoglobin during weeks 1-3). ",
    "d_hb = d_hb_ref + e_auc_d_hb * AUCss with the printed slope 0.011 g/dL per ug*day/mL (R = 0.72) and an intercept ",
    "digitized from the fitted line in Figure 2A. There is no PK layer and no ODE: exposure is supplied as the AUC_LUSP ",
    "data column, which the source analysis derived as starting dose / individual CL/F from the companion population PK ",
    "model packaged as Chen_2021_luspatercept. One of four exposure-response models in the Chen_2021_luspatercept_* family."
  )
  reference <- paste(
    "Chen N, Kassir N, Laadem A, Giuseppi AC, Shetty J, Maxwell SE, Sriraman P, Ritland S, Linde PG, Budda B, Reynolds JG, Zhou S, Palmisano M.",
    "Population Pharmacokinetics and Exposure-Response Relationship of Luspatercept, an Erythroid Maturation Agent, in Anemic Patients With beta-Thalassemia.",
    "J Clin Pharmacol. 2021;61(1):52-63. doi:10.1002/jcph.1696. PMCID: PMC7754485.",
    "Linear regression in Figure 2A (panel annotation 'Slope = 0.011; R = 0.72; P < 0.0001');",
    "Results, 'Exposure-Response for Hb'.",
    sep = " "
  )
  vignette <- "Chen_2021_luspatercept"
  units <- list(
    time = "n/a (static landmark exposure-response regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the AUC_LUSP covariate column)",
    concentration = "d_hb (average change from baseline in hemoglobin over weeks 1-3, g/dL)"
  )

  covariateData <- list(
    AUC_LUSP = list(
      description = "Individual luspatercept steady-state area under the serum concentration-time curve over the 21-day dosing interval for the starting dose (AUCss). Supplied as data: this model has no PK layer.",
      units = "ug*day/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Chen 2021 Methods: AUCss = starting dose / individual CL/F, with CL/F the empirical Bayes estimate from the final population PK model (Chen_2021_luspatercept). Entered untransformed and uncentred, so the intercept is an extrapolation to zero exposure (placebo patients were not part of this fit). Observed range in Figure 2A about 15-212 ug*day/mL; Results report a mean AUCss of 129 ug*day/mL after 1 mg/kg.",
      source_name = "AUCss"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 34L,
    n_studies = 1L,
    n_observations = "34 per-patient averages of weekly hemoglobin change from baseline during weeks 1-3",
    age_range = "Not reported for the exposure-response subset (PK population 18-66 years)",
    weight_range = "Not reported for the exposure-response subset (PK population 34.1-97.0 kg)",
    sex_female_pct = NA_real_,
    disease_state = "Adults with beta-thalassemia and a baseline transfusion burden < 12 RBC units/24 weeks; only hemoglobin values > 14 days after a transfusion were included.",
    dose_range = "Single dose levels 0.2, 0.4, 0.6, 0.8, 1 and 1.25 mg/kg subcutaneously (first dose of q3w dosing).",
    regions = "Phase 2 study A536-04 (NCT01749540), Greece and Italy.",
    notes = "The response is the average of the weekly hemoglobin changes from baseline during the first dosing interval (weeks 1-3), chosen because the greatest hemoglobin increase within a dosing interval followed the first dose and was not confounded by dose modification (Results)."
  )

  ini({
    # Chen 2021 Figure 2A: ordinary linear regression of the average weekly
    # hemoglobin change from baseline in cycle 1 on AUCss,
    #
    #   d_hb = d_hb_ref + e_auc_d_hb * AUCss
    #
    # The slope is printed in the panel. The intercept is not printed; it is
    # taken from the fitted line, which is drawn as a vector path in the PDF
    # and was read exactly against the axis tick marks. The intercept below
    # makes a line of the printed slope pass through the midpoint of the
    # drawn line (AUCss 113.8, d_hb 0.924 g/dL). The drawn line alone gives
    # slope 0.01057 and intercept -0.279, and an OLS refit of the 33
    # digitized scatter points gives 0.01069 and -0.303 with R = 0.722
    # (printed R = 0.72). See the vignette source trace.
    d_hb_ref <- -0.327; label("Average change from baseline in hemoglobin over weeks 1-3 extrapolated to zero luspatercept exposure (g/dL)")  # digitized from Figure 2A fitted line (not printed); see vignette
    e_auc_d_hb <- 0.011; label("Slope of the average hemoglobin change over weeks 1-3 on luspatercept AUCss (g/dL per ug*day/mL)")  # Figure 2A panel annotation: 'Slope = 0.011; R = 0.72; P < 0.0001'

    # Residual SD of the regression is not printed. The value below is the
    # residual standard error of an OLS refit to the 33 scatter points
    # digitized from Figure 2A (the refit reproduces the printed slope and R).
    addSd_d_hb <- 0.511; label("Residual SD of the linear regression (g/dL)")  # not printed; residual SE of an OLS refit to the digitized Figure 2A points; see vignette
  })

  model({
    d_hb <- d_hb_ref + e_auc_d_hb * AUC_LUSP
    d_hb ~ add(addSd_d_hb)
  })
}
