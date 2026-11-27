Qin_2026_ropeginterferon_jak2_month6 <- function() {
  description <- paste0(
    "Linear-regression exposure-efficacy model for the change from ",
    "baseline in JAK2 V617F allele burden at month 6 (SURPASS-ET) or ",
    "month 7 (EXCEED-ET) of subcutaneous ropeginterferon alfa-2b ",
    "(ropeg) in 87 patients with essential thrombocythemia. The ",
    "predicted change is 2.45 - 0.127 * CAV - 0.079 * ",
    "TUM_JAK2_V617F_VAF percentage points, where CAV is the individual ",
    "average total serum ropeg concentration from the first dose to ",
    "month 6 or 7 (ng/mL) and TUM_JAK2_V617F_VAF is the baseline allele ",
    "burden (per cent). Negative values are reductions. No PK layer, no ",
    "ODE and no between-subject random effect. Companion models in the ",
    "Qin_2026_ropeginterferon_* family."
  )
  reference <- paste(
    "Qin A, Masarova L, Chen H, Xiao Z, Jin J, Gill H, Reeves BN,",
    "El Chaer F, Kirito K, Komatsu N, Abu-Zeinah G, Shimoda K,",
    "Zagrijtschuk O, Sato T, Wu D, Wu L, Su X, Gao Y, Cui J, Zhang L,",
    "Mesa RA.",
    "Exposure-Response Analyses of Ropeginterferon Alfa-2b in",
    "Essential Thrombocythemia.",
    "J Clin Pharmacol. 2026;66(10):e70258. doi:10.1002/jcph.70258.",
    sep = " "
  )
  vignette <- "Qin_2026_ropeginterferon"
  units <- list(
    time = "n/a (static month-6/7 landmark regression; no time dimension)",
    dosing = "n/a (no dose events; the dosing history enters only through the CAV exposure column)",
    concentration = "djak2v617f (change from baseline in JAK2 V617F allele burden at month 6 or 7, percentage points; negative = reduction)"
  )

  covariateData <- list(
    CAV = list(
      description = "Individual average total serum ropeginterferon alfa-2b concentration from the first dose to month 6 (SURPASS-ET) or month 7 (EXCEED-ET).",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "TOTAL serum ropeg (the PK model's Cc). Derived by simulation",
        "from the actual dosing records and individual post hoc PK",
        "parameters of modellib('Qin_2026_ropeginterferon'); for",
        "patients treated for less than 6 or 7 months the average runs",
        "from the first dose to 2 weeks after the last dose (Qin 2026",
        "Methods, Exposure-Efficacy Analysis). Not centred or scaled.",
        "Observed range about 0-87 ng/mL (Qin 2026 Figure 3a)."
      ),
      source_name = "Cavg,0 to 6 or 7 months"
    ),
    TUM_JAK2_V617F_VAF = list(
      description = "Baseline JAK2 V617F variant allele frequency (allele burden).",
      units = "percent (0-100)",
      type = "continuous",
      reference_category = NULL,
      notes = "Uncentred linear covariate, beta2 = -0.079 per percentage point (Qin 2026 Table 4). Figure 3c draws the fitted lines at 19.55, 31.3 and 46.5%, the quartiles in the regression population.",
      source_name = "baseline JAK2V617F allele burden"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 87L,
    n_studies = 2L,
    n_observations = "87 patients with a JAK2 V617F measurement at month 6 or 7 (Qin 2026 Results, Exposure-Efficacy Analysis)",
    age_range = "ET patients: median 57.0 (EXCEED-ET) and 61.0 (SURPASS-ET) years, range 21.0-84.0 (Qin 2026 Table 1; the regression subset is not tabulated separately)",
    race_ethnicity = "SURPASS-ET 94.5% Asian; EXCEED-ET 77.3% Caucasian (Qin 2026 Table 1)",
    disease_state = "Essential thrombocythemia, JAK2 V617F-positive patients of SURPASS-ET and EXCEED-ET",
    dose_range = "HIDAT regimen: 250 ug at week 0, 350 ug at week 2, 500 ug every 2 weeks from week 4, adjusted for tolerability",
    regions = "China and other Asian sites (SURPASS-ET), United States and Canada (EXCEED-ET)",
    notes = "Ordinary least-squares linear regression in R 4.2.2; covariate forward selection at P = .01. Exposure effect P = .021 (Table 4); exploratory F-test P = .008."
  )

  ini({
    # Qin 2026 Table 4, rows 'Linear regression of JAK2V617F at Month 6
    # or 7'; model form Equation (3):
    #   reduction = beta0 + beta1*Exposure_i + beta2*JAK2_baseline,i
    # Both regressors uncentred: Figure 3c's line for a 19.55% baseline
    # burden starts at 2.45 - 0.079*19.55 = +0.91 at CAV = 0, and the
    # figure reads about +1.
    intercept <- 2.45; label("Predicted month 6/7 change in JAK2 V617F allele burden at CAV = 0 ng/mL and a baseline burden of 0% (percentage points)") # Qin 2026 Table 4: Intercept beta0 2.45 (SE 2.41), P = .310545
    e_cav_jak2 <- -0.127; label("Change in the month 6/7 allele-burden response per 1 ng/mL of CAV (percentage points per ng/mL)") # Qin 2026 Table 4: Exposure effect beta1 -0.127 (SE 0.054), P = .021201
    e_jak2vaf_jak2 <- -0.079; label("Change in the month 6/7 allele-burden response per percentage point of baseline JAK2 V617F burden (unitless)") # Qin 2026 Table 4: Baseline JAK2 effect beta2 -0.079 (SE 0.038), P = .043368

    # The residual standard deviation of the regression is not reported;
    # this fixed placeholder only gives rxode2 an error model for the
    # typical-value prediction.
    addSd_djak2v617f <- fixed(0.001); label("Placeholder additive residual SD on the typical-value predicted change; not published by the source") # not from source; see vignette Assumptions and deviations
  })

  model({
    djak2v617f <- intercept + e_cav_jak2 * CAV + e_jak2vaf_jak2 * TUM_JAK2_V617F_VAF
    djak2v617f ~ add(addSd_djak2v617f)
  })
}
