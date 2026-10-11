Wang_2022_abrocitinib_QTcF <- function() {
  description <- paste(
    "Prespecified linear mixed-effects concentration-QTc model for the",
    "Janus kinase 1 inhibitor abrocitinib in 36 healthy adult volunteers",
    "from a randomized, placebo- and moxifloxacin-controlled, 3-period",
    "crossover thorough-QT study (NCT03386279; single oral abrocitinib",
    "600 mg). The endpoint is the change from baseline in the",
    "Fridericia-corrected QT interval (DeltaQTcF, ms):",
    "DeltaQTcF = (theta0 + eta0) + (theta1 + eta1) * C",
    "+ theta2 * ON_TREATMENT + theta3 * (QTC_BL - 400) + theta_k(NTIME) + eps,",
    "with theta0 = -0.0062 ms, slope theta1 = 0.0026 ms/(ng/mL), a",
    "treatment-specific intercept theta2 = 0.3141 ms, a centered",
    "baseline-QTcF coefficient theta3 = -0.0577 ms/ms and seven",
    "nominal-time intercept shifts (0.5, 1, 2, 3, 6, 12 and 24 h against",
    "the 0.25 h reference). Additive, correlated random effects on the",
    "intercept (SD 3.315 ms) and slope (SD 0.0014 ms/(ng/mL), correlation",
    "-0.5294); additive residual SD 5.075 ms. The time-matched",
    "placebo-corrected effect is DeltaDeltaQTcF = theta2 + theta1 * C",
    "(6.0 ms at the 2156 ng/mL supratherapeutic concentration). PD-only",
    "model: abrocitinib plasma concentration is supplied as the",
    "time-varying covariate CP_ABROCITINIB_NGML (ng/mL); placebo records",
    "carry 0. The companion concentration-DeltaHR model is",
    "Wang_2022_abrocitinib_HR.R, and Wojciechowski_2022_abrocitinib.R is",
    "an abrocitinib population PK model that can supply the concentration",
    "trajectory.",
    sep = " "
  )

  reference <- paste(
    "Wang X, Gupta P, Malhotra BK, Farooqui SA, Le VH, Wojciechowski J,",
    "Mukherjee A, Nicholas T. Population Pharmacokinetic/Pharmacodynamic",
    "Modeling of the Effect of Abrocitinib on QT Intervals in Healthy",
    "Volunteers. Clin Pharmacol Drug Dev. 2022;11(9):1036-1045.",
    "doi:10.1002/cpdd.1111.",
    sep = " "
  )

  vignette <- "Wang_2022_abrocitinib_QTc"

  units <- list(
    time = "h",
    dosing = "(none; PD-only model fed by an external abrocitinib plasma-concentration covariate)",
    concentration = "(observation QTcF is the CHANGE FROM BASELINE in the Fridericia-corrected QT interval, DeltaQTcF, in ms; driving covariate CP_ABROCITINIB_NGML is in ng/mL)"
  )

  covariateData <- list(
    CP_ABROCITINIB_NGML = list(
      description = "Instantaneous total abrocitinib (parent) plasma concentration at the time of each ECG observation, supplied as a time-varying covariate from observed plasma samples or an upstream PK source.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying per event row. Drives the linear concentration-DeltaQTcF term (theta1 + eta1) * CP_ABROCITINIB_NGML; the slope is reported directly in ms per ng/mL (Wang 2022 Table 1), so no in-model unit rescaling is required.",
        "In Wang 2022 this was the observed abrocitinib concentration from blood samples drawn at the triplicate-ECG time points (0.25, 0.5, 1, 2, 3, 6, 12 and 24 h post-dose); concentrations below the 1 ng/mL quantitation limit were entered as 0 (Methods 'Pharmacokinetic Evaluations').",
        "Placebo records carry 0 ('Cijk is 0 for placebo', Methods Equation 1 legend).",
        "Reference values: maximum observed concentration > 4000 ng/mL after the 600 mg single dose; 1701 ng/mL quoted for healthy volunteers and 2156 ng/mL (the predicted 200 mg once-daily steady-state Cmax in atopic dermatitis, 1123 ng/mL, times the 1.92-fold fluconazole DDI increase) as the supratherapeutic concentration (Wang 2022 Methods and Table 3).",
        "Parent only: the active metabolites M1 and M2 were not measured or modelled (Discussion)."
      ),
      source_name = "Cijk (abrocitinib plasma concentration)"
    ),
    ON_TREATMENT = list(
      description = "Treatment-period indicator: 1 = abrocitinib 600 mg period, 0 = placebo period. Time-varying within subject in this crossover design (constant within each period).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (placebo period)",
      notes = paste(
        "Source paper writes this as TRT, 'a categorical covariate that takes the value of 0 with placebo and 1 for active treatment', whose coefficient theta2 is the 'treatment-specific intercept' included 'to capture any potential drug effect ... with minimal exposure' (Wang 2022 Methods, Equation 1 legend).",
        "Moxifloxacin periods were excluded from the exposure-response analysis (Methods 'Treatment Received'), so the indicator has only the placebo and abrocitinib levels.",
        "The time-matched placebo-corrected effect is theta2 * 1 + theta1 * C minus theta2 * 0, i.e. DeltaDeltaQTcF = theta2 + theta1 * C."
      ),
      source_name = "TRT"
    ),
    QTC_BL = list(
      description = "Subject's pre-dose baseline Fridericia-corrected QT interval for the treatment period (QTcF_ij0). Time-fixed within a period. Enters the intercept as the centered term e_qtc_bl_e0 * (QTC_BL - 400); set QTC_BL = 400 ms for the typical subject.",
      units = "ms",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Per subject AND per treatment period (subscript ij0 in Wang 2022 Equation 3): in the crossover each period has its own pre-dose baseline from triplicate ECGs at -1, -0.5 and 0 h (Methods 'Electrocardiography').",
        "Fridericia correction (beta = 0.333), selected over Bazett and over the estimated study-specific beta = 0.3375 (Results).",
        "Centering reference: the paper centers on QTcF0, 'the overall mean of all the baseline' values, but does not print that mean. It reports only that mean baseline QTcF ranged from 399.2 to 400.7 ms across the treatment groups (Results), so the model uses the rounded 400 ms. The overall mean lies inside that range, so the centering error is under 0.8 ms and shifts the prediction by under 0.05 ms (0.8 * 0.0577)."
      ),
      source_name = "QTcF_ij0 (baseline QTcF)"
    ),
    NTIME = list(
      description = "Nominal post-dose ECG time point in hours. Levels 0.25 (reference), 0.5, 1, 2, 3, 6, 12 and 24 h; each non-reference level carries its own intercept shift.",
      units = "h",
      type = "categorical",
      reference_category = "0.25 (the first post-dose nominal time point)",
      notes = paste(
        "Wang 2022 Equation 1/3 term theta_k: 'the covariate effect associated with each (apart from the first) nominal sampling time k across all treatment phases', included 'to account for diurnal variation'. The reference intercept theta0 refers to 'placebo, at the first postdose sampling time' (Methods).",
        "Nominal post-dose ECG times were 0.25, 0.5, 1, 2, 3, 6, 12 and 24 h (Methods 'Electrocardiography'). The reference level is 0.25 h, not 0 h: the pre-dose readings form the baseline and are not modelled observations.",
        "model() bins the supplied value to the nearest scheduled level (cut points 0.375, 0.75, 1.5, 2.5, 4.5, 9 and 18 h) so that a value off the nominal grid is assigned to a real level rather than silently falling through to the reference.",
        "Supply the NOMINAL time, not rxode2's tad()."
      ),
      source_name = "nominal sampling time k (theta_k)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 36L,
    n_studies = 1L,
    n_observations = 864L,
    age_range = "19-55 years (mean 33.2)",
    weight_range = "55.8-105.4 kg (mean 80.2)",
    sex_female_pct = NA_real_,
    race_ethnicity = c(White = 86.1, Black = 8.3, Other = 5.6),
    disease_state = "Healthy adult volunteers (BMI 17.5-30.5 kg/m^2, body weight > 50 kg, no clinically relevant ECG or medical abnormality).",
    dose_range = paste(
      "Single oral doses of abrocitinib 600 mg, placebo and moxifloxacin",
      "400 mg in a randomized 3-period crossover with >= 5-day washouts;",
      "moxifloxacin data (the positive control) were not used in the",
      "exposure-response analysis."
    ),
    regions = "Belgium (Pfizer Clinical Research Unit, Brussels)",
    notes = paste(
      "Thorough-QT study NCT03386279. Triplicate 12-lead ECGs at -1, -0.5 and 0 h pre-dose and 0.25, 0.5, 1, 2, 3, 6, 12 and 24 h post-dose, with time-matched PK samples (Methods).",
      "Mean baseline HR 54.1-54.7 bpm and mean baseline QTcF 399.2-400.7 ms across treatment groups; mean BMI 24.7 kg/m^2 (range 18-30.5) (Results).",
      "Sex distribution is not reported in the paper.",
      "Fitted with nlme::lme (method = 'ML') in R 3.4.1 (Methods).",
      "Model-predicted DeltaDeltaQTcF (90% CI) was 4.8 (3.63-5.98) ms at 1701 ng/mL and 6.0 (4.52-7.49) ms at the 2156 ng/mL supratherapeutic concentration (Table 3); the upper bound is below the 10 ms regulatory threshold."
    )
  )

  ini({
    # Wang 2022 Equation 3 (prespecified linear mixed-effects model):
    #   DeltaQTcF_ijk = (theta0 + eta0_i) + (theta1 + eta1_i) * C_ijk
    #                 + theta2 * TRT + theta3 * (QTcF_ij0 - QTcF0)
    #                 + theta_k + eps_ijk
    # All values from Wang 2022 Table 1 ('Parameter Estimates for the
    # Concentration-DeltaQTcF Modeling'). Every parameter is on the
    # linear ms scale and both random effects are additive, as in the
    # source nlme::lme fit.

    e0 <- -0.0062
    label("Intercept theta0 on DeltaQTcF: placebo, 0.25 h, mean baseline (ms)")
    # Table 1 'theta0: intercept (ms)' = -0.0062 (90% CI -1.41, 1.398)

    slope <- 0.0026
    label("Concentration-DeltaQTcF slope theta1 (ms per ng/mL)")
    # Table 1 'theta1: slope (ms/[ng/mL])' = 0.0026 (90% CI 0.0018, 0.0035)

    e_on_treatment_e0 <- 0.3141
    label("Treatment-specific intercept theta2 for the abrocitinib period (ms)")
    # Table 1 'theta2: treatment-specific intercept (ms)' = 0.3141 (90% CI -0.6115, 1.24)

    e_qtc_bl_e0 <- -0.0577
    label("Effect of centered baseline QTcF on the intercept, theta3 (ms per ms)")
    # Table 1 'theta3: baseline effect' = -0.0577 (90% CI -0.1043, -0.0112)

    # Nominal-time intercept shifts theta_k against the 0.25 h reference.
    e_ntime0p5_e0 <- -2.678
    label("Effect of the 0.5 h nominal time point on the intercept (ms)")
    # Table 1 'theta4: time effect (0.5 h) (ms)' = -2.678 (90% CI -4.088, -1.268)

    e_ntime1_e0 <- -1.851
    label("Effect of the 1 h nominal time point on the intercept (ms)")
    # Table 1 'theta5: time effect (1 h) (ms)' = -1.851 (90% CI -3.281, -0.4212)

    e_ntime2_e0 <- -3.152
    label("Effect of the 2 h nominal time point on the intercept (ms)")
    # Table 1 'theta6: time effect (2 h) (ms)' = -3.152 (90% CI -4.613, -1.6905)

    e_ntime3_e0 <- -2.666
    label("Effect of the 3 h nominal time point on the intercept (ms)")
    # Table 1 'theta7: time effect (3 h) (ms)' = -2.666 (90% CI -4.177, -1.155)

    e_ntime6_e0 <- -7.052
    label("Effect of the 6 h nominal time point on the intercept (ms)")
    # Table 1 'theta8: time effect (6 h) (ms)' = -7.052 (90% CI -8.472, -5.631)

    e_ntime12_e0 <- -3.402
    label("Effect of the 12 h nominal time point on the intercept (ms)")
    # Table 1 'theta9: time effect (12 h) (ms)' = -3.402 (90% CI -4.796, -2.007)

    e_ntime24_e0 <- -5.103
    label("Effect of the 24 h nominal time point on the intercept (ms)")
    # Table 1 'theta10: time effect (24 h) (ms)' = -5.103 (90% CI -6.499, -3.708)

    # Table 1 reports the random effects as standard deviations plus a
    # correlation ('IIV and RV terms are reported as standard deviations'):
    #   SD(eta0) = 3.315 ms, SD(eta1) = 0.0014 ms/(ng/mL), corr = -0.5294.
    # Converted to the variance-covariance block:
    #   var(eta0) = 3.315^2 = 10.989225
    #   cov       = -0.5294 * 3.315 * 0.0014 = -0.0024569454
    #   var(eta1) = 0.0014^2 = 1.96e-06
    etae0 + etaslope ~ c(
      10.989225,
      -0.0024569454, 1.96e-06
    )

    addSd <- 5.075
    label("Additive residual error SD on DeltaQTcF (ms)")
    # Table 1 'epsilon: RV' = 5.075 (90% CI 4.819, 5.345), an SD
  })

  model({
    # Nominal-time class effect: bin NTIME to the nearest scheduled level;
    # NTIME < 0.375 is the 0.25 h reference (all indicators 0).
    nt0p5 <- (NTIME >= 0.375) * (NTIME < 0.75)
    nt1 <- (NTIME >= 0.75) * (NTIME < 1.5)
    nt2 <- (NTIME >= 1.5) * (NTIME < 2.5)
    nt3 <- (NTIME >= 2.5) * (NTIME < 4.5)
    nt6 <- (NTIME >= 4.5) * (NTIME < 9)
    nt12 <- (NTIME >= 9) * (NTIME < 18)
    nt24 <- (NTIME >= 18)
    ntime_effect <-
      e_ntime0p5_e0 * nt0p5 + e_ntime1_e0 * nt1 + e_ntime2_e0 * nt2 +
      e_ntime3_e0 * nt3 + e_ntime6_e0 * nt6 + e_ntime12_e0 * nt12 +
      e_ntime24_e0 * nt24

    # Baseline QTcF centered on the overall mean baseline; 400 ms is the
    # rounded value inside the reported 399.2-400.7 ms group-mean range.
    qtc_bl_ref <- 400
    qtc_bl_centered <- QTC_BL - qtc_bl_ref

    e0_i <- e0 + etae0 +
      e_on_treatment_e0 * ON_TREATMENT +
      e_qtc_bl_e0 * qtc_bl_centered +
      ntime_effect
    slope_i <- slope + etaslope

    # Change from baseline in QTcF (ms)
    QTcF <- e0_i + slope_i * CP_ABROCITINIB_NGML
    QTcF ~ add(addSd)
  })
}
