Wang_2022_abrocitinib_HR <- function() {
  description <- paste(
    "Prespecified linear mixed-effects concentration-heart-rate model for",
    "the Janus kinase 1 inhibitor abrocitinib in 36 healthy adult",
    "volunteers from a randomized, placebo- and moxifloxacin-controlled,",
    "3-period crossover thorough-QT study (NCT03386279; single oral",
    "abrocitinib 600 mg). The endpoint is the change from baseline in",
    "heart rate (DeltaHR, bpm):",
    "DeltaHR = (theta0 + eta0) + (theta1 + eta1) * C",
    "+ theta2 * ON_TREATMENT + theta3 * (HR - 54.4) + theta_k(NTIME) + eps,",
    "with theta0 = -0.669 bpm, slope theta1 = 0.0031 bpm/(ng/mL), a",
    "treatment-specific intercept theta2 = -0.1979 bpm, a centered",
    "baseline-HR coefficient theta3 = -0.1544 bpm/bpm and seven",
    "nominal-time intercept shifts (0.5, 1, 2, 3, 6, 12 and 24 h against",
    "the 0.25 h reference). Additive, correlated random effects on the",
    "intercept (SD 1.594 bpm) and slope (SD 0.0016 bpm/(ng/mL),",
    "correlation -0.4389); additive residual SD 3.264 bpm. The",
    "time-matched placebo-corrected effect is DeltaDeltaHR = theta2 +",
    "theta1 * C (6.5 bpm at the 2156 ng/mL supratherapeutic",
    "concentration, below the 10 bpm level at which Fridericia correction",
    "is considered adequate). PD-only model: abrocitinib plasma",
    "concentration is supplied as the time-varying covariate",
    "CP_ABROCITINIB_NGML (ng/mL); placebo records carry 0. The companion",
    "concentration-DeltaQTcF model is Wang_2022_abrocitinib_QTcF.R.",
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
    concentration = "(observation d_hr is the change from baseline in heart rate, bpm; driving covariate CP_ABROCITINIB_NGML is in ng/mL)"
  )

  covariateData <- list(
    CP_ABROCITINIB_NGML = list(
      description = "Instantaneous total abrocitinib (parent) plasma concentration at the time of each ECG observation, supplied as a time-varying covariate from observed plasma samples or an upstream PK source.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying per event row. Drives the linear concentration-DeltaHR term (theta1 + eta1) * CP_ABROCITINIB_NGML; the slope is reported directly in bpm per ng/mL (Wang 2022 Table 2).",
        "Same time-matched PK-ECG dataset as the companion C-DeltaQTcF model: observed concentrations at 0.25-24 h post-dose, values below the 1 ng/mL quantitation limit entered as 0, placebo records 0 (Wang 2022 Methods).",
        "Reference values: 1701 ng/mL (healthy volunteers) and 2156 ng/mL (supratherapeutic concentration in atopic dermatitis) in Wang 2022 Table 3."
      ),
      source_name = "Cijk (abrocitinib plasma concentration)"
    ),
    ON_TREATMENT = list(
      description = "Treatment-period indicator: 1 = abrocitinib 600 mg period, 0 = placebo period. Time-varying within subject in this crossover design (constant within each period).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (placebo period)",
      notes = paste(
        "Source paper writes this as TRT (0 = placebo, 1 = active), whose coefficient theta2 is the 'treatment-specific intercept' included 'to capture any potential drug effect on HR with minimal exposure' (Wang 2022 Methods, Equation 1).",
        "Moxifloxacin periods were excluded from the exposure-response analysis.",
        "DeltaDeltaHR = theta2 + theta1 * C."
      ),
      source_name = "TRT"
    ),
    HR = list(
      description = "Subject's pre-dose baseline heart rate for the treatment period (HR_ij0). Time-fixed within a period. Enters the intercept as the centered term e_hr_bl_e0 * (HR - 54.4); set HR = 54.4 bpm for the typical subject.",
      units = "beats/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "This is the BASELINE reading of the canonical HR covariate, not an observation-time vital sign, following the C-DeltaHR convention of Mukker_2026_tuvusertib_HR.R; the model's observable is named d_hr so that it does not shadow this column.",
        "Per subject AND per treatment period (subscript ij0 in Wang 2022 Equation 1), from the triplicate pre-dose ECGs at -1, -0.5 and 0 h.",
        "Centering reference: the paper centers on HR0, 'the overall mean of all the baseline HR values', but does not print it. Mean baseline HR ranged from 54.1 to 54.7 bpm across the treatment groups (Results), so the model uses the midpoint 54.4 bpm. The overall mean lies inside that range, so the centering error is under 0.3 bpm and shifts the prediction by under 0.05 bpm (0.3 * 0.1544)."
      ),
      source_name = "HR_ij0 (baseline HR)"
    ),
    NTIME = list(
      description = "Nominal post-dose ECG time point in hours. Levels 0.25 (reference), 0.5, 1, 2, 3, 6, 12 and 24 h; each non-reference level carries its own intercept shift.",
      units = "h",
      type = "categorical",
      reference_category = "0.25 (the first post-dose nominal time point)",
      notes = paste(
        "Wang 2022 Equation 1 term theta_k, 'the covariate effect associated with each (apart from the first) nominal sampling time k', included 'to account for diurnal variation'; theta0 refers to 'placebo, at the first postdose sampling time' (Methods).",
        "model() bins the supplied value to the nearest scheduled level (cut points 0.375, 0.75, 1.5, 2.5, 4.5, 9 and 18 h).",
        "Supply the NOMINAL time, not rxode2's tad()."
      ),
      source_name = "nominal sampling time k (theta_k)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 36L,
    n_studies = 1L,
    age_range = "19-55 years (mean 33.2)",
    weight_range = "55.8-105.4 kg (mean 80.2)",
    sex_female_pct = NA_real_,
    race_ethnicity = c(White = 86.1, Black = 8.3, Other = 5.6),
    disease_state = "Healthy adult volunteers (BMI 17.5-30.5 kg/m^2, body weight > 50 kg, no clinically relevant ECG or medical abnormality).",
    dose_range = paste(
      "Single oral doses of abrocitinib 600 mg, placebo and moxifloxacin",
      "400 mg in a randomized 3-period crossover with >= 5-day washouts;",
      "moxifloxacin data were not used in the exposure-response analysis."
    ),
    regions = "Belgium (Pfizer Clinical Research Unit, Brussels)",
    notes = paste(
      "Thorough-QT study NCT03386279; same time-matched PK-ECG dataset as Wang_2022_abrocitinib_QTcF.R. Mean baseline HR 54.1-54.7 bpm across treatment groups (Results).",
      "Sex distribution is not reported in the paper.",
      "Fitted with nlme::lme (method = 'ML') in R 3.4.1 (Methods).",
      "Model-predicted DeltaDeltaHR (90% CI) was 5.1 (4.09-6.11) bpm at 1701 ng/mL and 6.51 (5.23-7.80) bpm at 2156 ng/mL (Table 3 and Results); being below 10 bpm, Fridericia correction was considered appropriate for the QTc analysis."
    )
  )

  ini({
    # Wang 2022 Equation 1 (prespecified linear mixed-effects model):
    #   DeltaHR_ijk = (theta0 + eta0_i) + (theta1 + eta1_i) * C_ijk
    #               + theta2 * TRT + theta3 * (HR_ij0 - HR0)
    #               + theta_k + eps_ijk
    # All values from Wang 2022 Table 2 ('Parameter Estimates for the
    # Concentration-DeltaHR Modeling'); linear bpm scale, additive random
    # effects. The Table 2 footnote prints the intercept eta as 'eta0,1';
    # the Methods text and Equation 3 make clear it is eta0,i.

    e0 <- -0.669
    label("Intercept theta0 on DeltaHR: placebo, 0.25 h, mean baseline (bpm)")
    # Table 2 'theta0: intercept (bpm)' = -0.669 (90% CI -1.485, 0.1466)

    slope <- 0.0031
    label("Concentration-DeltaHR slope theta1 (bpm per ng/mL)")
    # Table 2 'theta1: slope (bpm/[ng/mL])' = 0.0031 (90% CI 0.0024, 0.0038)

    e_on_treatment_e0 <- -0.1979
    label("Treatment-specific intercept theta2 for the abrocitinib period (bpm)")
    # Table 2 'theta2: treatment-specific intercept (bpm)' = -0.1979 (90% CI -0.7999, 0.4041)

    e_hr_bl_e0 <- -0.1544
    label("Effect of centered baseline HR on the intercept, theta3 (bpm per bpm)")
    # Table 2 'theta3: baseline effect' = -0.1544 (90% CI -0.2126, -0.0963)

    # Nominal-time intercept shifts theta_k against the 0.25 h reference.
    e_ntime0p5_e0 <- 0.0797
    label("Effect of the 0.5 h nominal time point on the intercept (bpm)")
    # Table 2 'theta4: time effect (0.5 h) (bpm)' = 0.0797 (90% CI -0.8296, 0.9891)

    e_ntime1_e0 <- -0.299
    label("Effect of the 1 h nominal time point on the intercept (bpm)")
    # Table 2 'theta5: time effect (1 h) (bpm)' = -0.299 (90% CI -1.222, 0.6243)

    e_ntime2_e0 <- -0.7794
    label("Effect of the 2 h nominal time point on the intercept (bpm)")
    # Table 2 'theta6: time effect (2 h) (bpm)' = -0.7794 (90% CI -1.724, 0.1652)

    e_ntime3_e0 <- -0.5218
    label("Effect of the 3 h nominal time point on the intercept (bpm)")
    # Table 2 'theta7: time effect (3 h) (bpm)' = -0.5218 (90% CI -1.500, 0.4569)

    e_ntime6_e0 <- 5.689
    label("Effect of the 6 h nominal time point on the intercept (bpm)")
    # Table 2 'theta8: time effect (6 h) (bpm)' = 5.689 (90% CI 4.774, 6.605)

    e_ntime12_e0 <- 4.274
    label("Effect of the 12 h nominal time point on the intercept (bpm)")
    # Table 2 'theta9: time effect (12 h) (bpm)' = 4.274 (90% CI 3.377, 5.172)

    e_ntime24_e0 <- 3.111
    label("Effect of the 24 h nominal time point on the intercept (bpm)")
    # Table 2 'theta10: time effect (24 h) (bpm)' = 3.111 (90% CI 2.213, 4.009)

    # Table 2 reports the random effects as standard deviations plus a
    # correlation: SD(eta0) = 1.594 bpm, SD(eta1) = 0.0016 bpm/(ng/mL),
    # corr = -0.4389. Converted to the variance-covariance block:
    #   var(eta0) = 1.594^2 = 2.540836
    #   cov       = -0.4389 * 1.594 * 0.0016 = -0.00111937056
    #   var(eta1) = 0.0016^2 = 2.56e-06
    etae0 + etaslope ~ c(
      2.540836,
      -0.00111937056, 2.56e-06
    )

    addSd <- 3.264
    label("Additive residual error SD on DeltaHR (bpm)")
    # Table 2 'epsilon: RV' = 3.264 (90% CI 3.096, 3.441), an SD
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

    # Baseline HR centered on the overall mean baseline; 54.4 bpm is the
    # midpoint of the reported 54.1-54.7 bpm group-mean range.
    hr_bl_ref <- 54.4
    hr_bl_centered <- HR - hr_bl_ref

    e0_i <- e0 + etae0 +
      e_on_treatment_e0 * ON_TREATMENT +
      e_hr_bl_e0 * hr_bl_centered +
      ntime_effect
    slope_i <- slope + etaslope

    # Change from baseline in heart rate (bpm)
    d_hr <- e0_i + slope_i * CP_ABROCITINIB_NGML
    d_hr ~ add(addSd)
  })
}
