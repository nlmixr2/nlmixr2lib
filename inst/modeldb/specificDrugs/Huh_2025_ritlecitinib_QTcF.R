Huh_2025_ritlecitinib_QTcF <- function() {
  description <- paste(
    "Prespecified linear mixed-effects concentration-QTc model (Garnett",
    "2018 white-paper form) for the JAK3/TEC-family kinase inhibitor",
    "ritlecitinib in 79 healthy volunteers from the single-ascending-dose",
    "(5-800 mg, placebo-controlled) period of the first-in-human study",
    "B7981001. The endpoint is the change from baseline in the",
    "Fridericia-corrected QT interval (DeltaQTcF, ms):",
    "DeltaQTcF = (theta0 + eta0) + (theta1 + eta1) * C",
    "+ theta2 * ON_TREATMENT + theta3 * (QTC_BL - 400) + theta_k(NTIME) + eps,",
    "with intercept theta0 = -2.76 ms at the 0.5 h reference time, slope",
    "theta1 = 0.000213 ms per ng/mL (0.213 ms per 1000 ng/mL, 95% CI",
    "including zero), treatment-specific intercept theta2 = -0.354 ms,",
    "baseline-QTcF coefficient theta3 = -0.208 and seven nominal-time",
    "intercept shifts (1, 2, 4, 8, 12, 16 and 24 h). Correlated random",
    "effects on intercept (SD 4.37 ms) and slope (SD 3.75e-7 ms per",
    "ng/mL, correlation 5.84e-5) and an additive residual SD of 5.05 ms.",
    "The time-matched placebo-corrected effect is DeltaDeltaQTcF =",
    "theta2 + theta1 * C. PD-only model: the ritlecitinib plasma",
    "concentration is supplied as the time-varying covariate",
    "CP_RITLECITINIB_NGML (ng/mL; 0 on placebo records);",
    "Wojciechowski_2023_ritlecitinib_final can supply the concentration",
    "trajectory."
  )

  reference <- paste(
    "Huh Y, Qiu R, Prybylski J, Wojciechowski J, Wang Y, Purohit VS.",
    "Comprehensive Safety Exposure-Response Analysis to Support",
    "Ritlecitinib Dose Selection. CPT Pharmacometrics Syst Pharmacol.",
    "2025;14(6):1108-1118. doi:10.1002/psp4.70030.",
    sep = " "
  )

  vignette <- "Huh_2025_ritlecitinib_safety"

  units <- list(
    time = "h",
    dosing = "(none; PD-only model fed by an external ritlecitinib plasma-concentration covariate)",
    concentration = "(observation QTcF is the CHANGE FROM BASELINE in the Fridericia-corrected QT interval, DeltaQTcF, in ms; driving covariate CP_RITLECITINIB_NGML is in ng/mL)"
  )

  covariateData <- list(
    CP_RITLECITINIB_NGML = list(
      description = "Instantaneous ritlecitinib plasma concentration at the time of each ECG observation, supplied as a time-varying covariate from observed plasma samples or an upstream PK source.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying per event row. Drives the linear term (theta1 + eta1) * CP_RITLECITINIB_NGML; the slope is reported in ms per ng/mL (Huh 2025 Table 2), so no in-model unit rescaling is required.",
        "In Huh 2025 this was the observed ritlecitinib concentration from PK samples time-matched to the ECGs at 0, 0.5, 1, 2, 4, 8, 12, 16 and 24 h post-dose (Table 1, study B7981001).",
        "Placebo records carry 0 (Equation 1 legend: 'Cijk is 0 for placebo')."
      ),
      source_name = "Cijk / CON1"
    ),
    ON_TREATMENT = list(
      description = "Treatment indicator: 1 = ritlecitinib, 0 = placebo.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (placebo)",
      notes = "Huh 2025 Equation 1: 'TRTj is a categorical covariate taking the value of 0 for placebo and 1 for active treatment'; theta2 is the treatment effect on the intercept (Supplement 2 nlme fixed formula a ~ factor(NBSL) + factor(TRT) + DQTFB). The time-matched placebo-corrected effect is DeltaDeltaQTcF = theta2 + theta1 * C.",
      source_name = "TRT"
    ),
    QTC_BL = list(
      description = "Subject's pre-dose baseline Fridericia-corrected QT interval (QTcF_ij0). Enters the intercept as the centered term e_qtc_bl_e0 * (QTC_BL - 400); set QTC_BL = 400 ms for the typical subject.",
      units = "ms",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Day 1 pre-dose morning ECG of the single-ascending-dose period (Results, Concentration-QTc Model).",
        "Huh 2025 centers the baseline 'around the overall mean of baseline QTcFs' but does not print that mean; the model uses 400 ms, a rounded value typical of healthy-volunteer QTcF (the same value Wang_2022_abrocitinib_QTcF uses for a Pfizer healthy-volunteer cohort). A centering error of d ms shifts every prediction by -0.208 * d ms and cancels in DeltaDeltaQTcF."
      ),
      source_name = "QTcFij0 (DQTFB in the nlme call)"
    ),
    NTIME = list(
      description = "Nominal post-dose ECG time point in hours. Levels 0.5 (reference), 1, 2, 4, 8, 12, 16 and 24 h; each non-reference level carries its own intercept shift.",
      units = "h",
      type = "categorical",
      reference_category = "0.5 (the first post-dose nominal time point)",
      notes = paste(
        "Huh 2025 Equation 1 theta_k, 'the covariate effect associated with each sampling time k'; Table 2 lists time effects for 1, 2, 4, 8, 12, 16 and 24 h, so the first post-dose time 0.5 h (Table 1) is the reference carried by the intercept. The 0 h reading is the baseline and is not a modelled observation.",
        "model() bins the supplied value to the nearest scheduled level (cut points 0.75, 1.5, 3, 6, 10, 14 and 20 h).",
        "Supply the NOMINAL time, not rxode2's tad()."
      ),
      source_name = "NBSL (theta_k)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 79L,
    n_studies = 1L,
    n_observations = 946L,
    disease_state = "Healthy volunteers",
    dose_range = "Single oral doses of 5, 20, 50, 100, 200, 400 and 800 mg or placebo (single-ascending-dose period of B7981001)",
    regions = "Not reported",
    notes = "Huh 2025 Results, Concentration-QTc Model: 946 post-baseline observations from 79 subjects. Only the single-ascending-dose period was used because baseline ECGs were collected differently in the multiple-ascending-dose period (Day -1 afternoon). Demographics of the C-QTc subset are not reported. Fitted with nlme in R (Supplement 2). The upper bound of the 90% CI of DeltaDeltaQTcF stayed below 10 ms over the whole concentration range, 4.13 ms at the mean Cmax of 800 mg (Figure 1)."
  )

  ini({
    # Huh 2025 Equation 1 and Table 2 ('Concentration-QTc model' rows).
    e0 <- -2.76
    label("Intercept theta0 on DeltaQTcF: placebo, 0.5 h, mean baseline (ms)")
    # Table 2 'Intercept (ms)' = -2.76 (95% CI -4.94, -0.58)

    slope <- 0.000213
    label("Concentration-DeltaQTcF slope theta1 (ms per ng/mL)")
    # Table 2 'Slope (ms/ng/mL)' = 0.000213 (95% CI -0.000362, 0.000788)

    e_on_treatment_e0 <- -0.354
    label("Treatment-specific intercept theta2 for ritlecitinib (ms)")
    # Table 2 'Treatment effect (ms)' = -0.354 (95% CI -2.62, 1.91)

    e_qtc_bl_e0 <- -0.208
    label("Effect of centered baseline QTcF on the intercept, theta3 (ms per ms)")
    # Table 2 'Baseline effect (ms)' = -0.208 (95% CI -0.274, -0.143)

    # Nominal-time intercept shifts theta_k against the 0.5 h reference.
    e_ntime1_e0 <- 0.702
    label("Effect of the 1 h nominal time point on the intercept (ms)")
    # Table 2 'Time effect-1 h (ms)' = 0.702 (95% CI -0.779, 2.18)

    e_ntime2_e0 <- 1.91
    label("Effect of the 2 h nominal time point on the intercept (ms)")
    # Table 2 'Time effect-2 h (ms)' = 1.91 (95% CI 0.440, 3.39)

    e_ntime4_e0 <- 1.16
    label("Effect of the 4 h nominal time point on the intercept (ms)")
    # Table 2 'Time effect-4 h (ms)' = 1.16 (95% CI -0.319, 2.65)

    e_ntime8_e0 <- -2.73
    label("Effect of the 8 h nominal time point on the intercept (ms)")
    # Table 2 'Time effect-8 h (ms)' = -2.73 (95% CI -4.24, -1.22)

    e_ntime12_e0 <- 1.09
    label("Effect of the 12 h nominal time point on the intercept (ms)")
    # Table 2 'Time effect-12 h (ms)' = 1.09 (95% CI -0.431, 2.61)

    e_ntime16_e0 <- 1.96
    label("Effect of the 16 h nominal time point on the intercept (ms)")
    # Table 2 'Time effect-16 h (ms)' = 1.96 (95% CI 0.435, 3.48)

    e_ntime24_e0 <- 2.25
    label("Effect of the 24 h nominal time point on the intercept (ms)")
    # Table 2 'Time effect-24 h (ms)' = 2.25 (95% CI 0.727, 3.77)

    # Table 2 labels the random-effect rows omega^2 and 'Covariance', but the
    # rows are the nlme::intervals() output: SD(eta0) = 4.37 ms, SD(eta1) =
    # 3.75e-7 ms/(ng/mL) (CI 2.87e-74 to 4.91e+60, the signature of an SD
    # estimated at zero), correlation 5.84e-5 (CI -0.255 to 0.255, symmetric
    # and bounded as a correlation is) and residual sigma 5.05 ms.
    #   var(eta0) = 4.37^2 = 19.0969
    #   cov       = 5.84e-5 * 4.37 * 3.75e-7 = 9.5703e-11
    #   var(eta1) = (3.75e-7)^2 = 1.40625e-13
    etae0 + etaslope ~ c(
      19.0969,
      9.5703e-11, 1.40625e-13
    )

    addSd <- 5.05
    label("Additive residual error SD on DeltaQTcF (ms)")
    # Table 2 'Residual Variability' = 5.05 (95% CI 4.78, 5.33), the nlme sigma
  })

  model({
    # Nominal-time class effect: bin NTIME to the nearest scheduled level;
    # NTIME < 0.75 is the 0.5 h reference (all indicators 0).
    nt1 <- (NTIME >= 0.75) * (NTIME < 1.5)
    nt2 <- (NTIME >= 1.5) * (NTIME < 3)
    nt4 <- (NTIME >= 3) * (NTIME < 6)
    nt8 <- (NTIME >= 6) * (NTIME < 10)
    nt12 <- (NTIME >= 10) * (NTIME < 14)
    nt16 <- (NTIME >= 14) * (NTIME < 20)
    nt24 <- (NTIME >= 20)
    ntime_effect <-
      e_ntime1_e0 * nt1 + e_ntime2_e0 * nt2 + e_ntime4_e0 * nt4 +
      e_ntime8_e0 * nt8 + e_ntime12_e0 * nt12 + e_ntime16_e0 * nt16 +
      e_ntime24_e0 * nt24

    # Baseline QTcF centered on the (unreported) overall mean; 400 ms used.
    qtc_bl_ref <- 400
    qtc_bl_centered <- QTC_BL - qtc_bl_ref

    e0_i <- e0 + etae0 +
      e_on_treatment_e0 * ON_TREATMENT +
      e_qtc_bl_e0 * qtc_bl_centered +
      ntime_effect
    slope_i <- slope + etaslope

    # Change from baseline in QTcF (ms)
    QTcF <- e0_i + slope_i * CP_RITLECITINIB_NGML
    QTcF ~ add(addSd)
  })
}
