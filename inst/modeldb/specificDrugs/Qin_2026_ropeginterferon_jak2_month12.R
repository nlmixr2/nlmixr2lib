Qin_2026_ropeginterferon_jak2_month12 <- function() {
  description <- paste0(
    "Linear-regression exposure-efficacy model for the change from ",
    "baseline in JAK2 V617F allele burden at month 12 (SURPASS-ET) or ",
    "month 13 (EXCEED-ET) of subcutaneous ropeginterferon alfa-2b ",
    "(ropeg) in 80 patients with essential thrombocythemia. The ",
    "predicted change is 1.94 - 0.147 * CAV - 0.103 * ",
    "TUM_JAK2_V617F_VAF percentage points, where CAV is the individual ",
    "average total serum ropeg concentration from the first dose to ",
    "month 12 or 13 (ng/mL) and TUM_JAK2_V617F_VAF is the baseline allele ",
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
    time = "n/a (static month-12/13 landmark regression; no time dimension)",
    dosing = "n/a (no dose events; the dosing history enters only through the CAV exposure column)",
    concentration = "djak2v617f (change from baseline in JAK2 V617F allele burden at month 12 or 13, percentage points; negative = reduction)"
  )

  covariateData <- list(
    CAV = list(
      description = "Individual average total serum ropeginterferon alfa-2b concentration from the first dose to month 12 (SURPASS-ET) or month 13 (EXCEED-ET).",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "TOTAL serum ropeg (the PK model's Cc). Derived by simulation",
        "from the actual dosing records and individual post hoc PK",
        "parameters of modellib('Qin_2026_ropeginterferon'); for",
        "patients treated for less than 12 or 13 months the average runs",
        "from the first dose to 2 weeks after the last dose (Qin 2026",
        "Methods, Exposure-Efficacy Analysis). Not centred or scaled.",
        "Observed range about 5-92 ng/mL (Qin 2026 Figure 3b)."
      ),
      source_name = "Cavg,0 to 12 or 13 months"
    ),
    TUM_JAK2_V617F_VAF = list(
      description = "Baseline JAK2 V617F variant allele frequency (allele burden).",
      units = "percent (0-100)",
      type = "continuous",
      reference_category = NULL,
      notes = "Uncentred linear covariate, beta2 = -0.103 per percentage point (Qin 2026 Table 4). Figure 3d draws the fitted lines at 19.55, 31.3 and 46.5%, the quartiles in the regression population.",
      source_name = "baseline JAK2V617F allele burden"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 80L,
    n_studies = 2L,
    n_observations = "80 patients with a JAK2 V617F measurement at month 12 or 13 (Qin 2026 Results, Exposure-Efficacy Analysis)",
    age_range = "ET patients: median 57.0 (EXCEED-ET) and 61.0 (SURPASS-ET) years, range 21.0-84.0 (Qin 2026 Table 1; the regression subset is not tabulated separately)",
    race_ethnicity = "SURPASS-ET 94.5% Asian; EXCEED-ET 77.3% Caucasian (Qin 2026 Table 1)",
    disease_state = "Essential thrombocythemia, JAK2 V617F-positive patients of SURPASS-ET and EXCEED-ET",
    dose_range = "HIDAT regimen: 250 ug at week 0, 350 ug at week 2, 500 ug every 2 weeks from week 4, adjusted for tolerability",
    regions = "China and other Asian sites (SURPASS-ET), United States and Canada (EXCEED-ET)",
    notes = "Ordinary least-squares linear regression in R 4.2.2; covariate forward selection at P = .01. Exposure effect P = .077 (Table 4, not significant at .05 in the final model); exploratory F-test P = .035."
  )

  ini({
    # Qin 2026 Table 4, rows 'Linear regression of JAK2V617F at Month 12
    # or 13'; model form Equation (3):
    #   reduction = beta0 + beta1*Exposure_i + beta2*JAK2_baseline,i
    # Both regressors uncentred: Figure 3d's line for a 46.5% baseline
    # burden starts at 1.94 - 0.103*46.5 = -2.85 at CAV = 0, and the
    # figure reads about -3.
    intercept <- 1.94; label("Predicted month 12/13 change in JAK2 V617F allele burden at CAV = 0 ng/mL and a baseline burden of 0% (percentage points)") # Qin 2026 Table 4: Intercept beta0 1.94 (SE 3.807), P = .610909
    e_cav_jak2 <- -0.147; label("Change in the month 12/13 allele-burden response per 1 ng/mL of CAV (percentage points per ng/mL)") # Qin 2026 Table 4: Exposure effect beta1 -0.147 (SE 0.082), P = .077276
    e_jak2vaf_jak2 <- -0.103; label("Change in the month 12/13 allele-burden response per percentage point of baseline JAK2 V617F burden (unitless)") # Qin 2026 Table 4: Baseline JAK2 effect beta2 -0.103 (SE 0.058), P = .082337

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
