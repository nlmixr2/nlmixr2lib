Sundell_2025_balcinrenone_QTcF <- function() {
  description <- paste(
    "Prespecified linear mixed-effects concentration-QTc model for the",
    "selective mineralocorticoid receptor modulator balcinrenone (AZD9977)",
    "in 62 healthy male participants of the phase I single ascending dose",
    "study D6400C00001 (NCT02484729, Part A; single oral suspension doses",
    "of 5-800 mg plus fractionated 800 and 1200 mg, or placebo, fasted),",
    "used as a thorough-QT study substitute. The endpoint is the change",
    "from baseline in the Fridericia-corrected QT interval (DeltaQTcF, ms):",
    "DeltaQTcF = (theta0 + eta0) + theta1 * ON_TREATMENT",
    "+ (theta2 + eta2) * C + theta3(NTIME) + theta4 * (QTC_BL - 400),",
    "with intercept theta0 = -0.267 ms, treatment-specific intercept",
    "theta1 = 0.359 ms, slope theta2 = -0.003 ms per umol/L (95% CI -0.181,",
    "0.176; not significant), twelve nominal-time intercept shifts",
    "(0.5-48 h against the pre-dose 0 h reference) and a centered",
    "baseline-QTcF coefficient theta4 = -0.052 ms/ms (primary analysis,",
    "all dose arms). The intercept and slope random effects (unstructured",
    "covariance) and the residual SD are not reported and are encoded as",
    "fixed(0), so the model returns typical-value predictions. The",
    "time-matched placebo-corrected effect is",
    "DeltaDeltaQTcF = theta1 + theta2 * C (0.35 ms at the 2.156 umol/L",
    "high clinical exposure scenario). PD-only model: balcinrenone plasma",
    "concentration is supplied as the time-varying covariate",
    "CP_BALCINRENONE_UM (umol/L); placebo records carry 0.",
    "Parkinson_2025_balcinrenone.R is a balcinrenone population PK model",
    "(in nmol/L) that can supply the concentration trajectory.",
    sep = " "
  )

  reference <- paste(
    "Sundell J, Rekic D, Melin J, Johansson S, Ebrahimi A, Dota C,",
    "Parkinson J. Concentration-QT modeling demonstrates that the selective",
    "mineralocorticoid receptor modulator, balcinrenone (AZD9977), does not",
    "prolong QT interval. CPT Pharmacometrics Syst Pharmacol.",
    "2025;14:292-301.",
    "doi:10.1002/psp4.13272. PMCID: PMC11812937.",
    "Fixed-effect estimates are from Supplementary Table S1 (primary",
    "analysis); the model structure is Equation 1 and the lme4 code in the",
    "supplementary materials.",
    sep = " "
  )

  vignette <- "Sundell_2025_balcinrenone"

  units <- list(
    time = "h",
    dosing = "(none; PD-only model fed by an external balcinrenone plasma-concentration covariate)",
    concentration = "(observation QTcF is the CHANGE FROM BASELINE in the Fridericia-corrected QT interval, DeltaQTcF, in ms; driving covariate CP_BALCINRENONE_UM is in umol/L)"
  )

  covariateData <- list(
    CP_BALCINRENONE_UM = list(
      description = "Instantaneous balcinrenone (parent) plasma concentration at the time of each ECG observation, supplied as a time-varying covariate from observed plasma samples or an upstream PK source.",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying per event row. Drives the linear concentration-DeltaQTcF term (theta2 + eta2) * CP_BALCINRENONE_UM; the slope is reported in ms per umol/L (Sundell 2025 Abstract and Table S1), so no in-model unit rescaling is required.",
        "In Sundell 2025 this was the observed (time-matched) balcinrenone concentration from the PK sample drawn immediately after each dECG recording (Methods 'ECG and PK sampling'); the assay quantitation limit was 10 nmol/L (0.01 umol/L).",
        "Placebo records carry 0 (Equation 2 evaluates the placebo arm at C = 0).",
        "Reference values: high clinical exposure scenario 2.156 umol/L (steady-state Cmax of 40 mg once daily, capsule, fed, with a strong CYP3A4 inhibitor) and highest observed concentration 16.7 umol/L (1200 mg fractionated) (Methods 'High clinical exposure scenario', Table 1).",
        "Convert from nmol/L (the units of Parkinson_2025_balcinrenone.R) by dividing by 1000."
      ),
      source_name = "C_ijl (observed balcinrenone concentration; 'DV' column in the supplementary lme4 code)"
    ),
    ON_TREATMENT = list(
      description = "Treatment-arm indicator: 1 = balcinrenone (any dose), 0 = placebo. Time-fixed within subject in this parallel-group design.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (placebo arm)",
      notes = paste(
        "Source paper writes this as TRT_j 'for participant i in treatment j (0 = placebo, 1 = active drug)', whose coefficient theta1 is 'the fixed effect associated with treatment' and 'the treatment-specific intercept' (Sundell 2025 Methods, Equation 1 legend).",
        "The supplementary lme4 code enters it as TRTID 'as factor'; Table S1 reports one coefficient, 'theta1 (active arm)' against the 'Placebo arm' reference, so the factor has two levels.",
        "The time-matched placebo-corrected effect is theta1 * 1 + theta2 * C minus theta1 * 0, i.e. DeltaDeltaQTcF = theta1 + theta2 * C (Equation 2)."
      ),
      source_name = "TRT_j / TRTID"
    ),
    QTC_BL = list(
      description = "Subject's pre-dose baseline Fridericia-corrected QT interval (QTcF_ijl=0). Time-fixed. Enters the intercept as the centered term e_qtc_bl_e0 * (QTC_BL - 400); set QTC_BL = 400 ms for the typical subject.",
      units = "ms",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline is the pre-dose dECG observation (Results 'Categorical evaluation': 'baseline was defined as the pre-dose observation').",
        "Fridericia correction, supported over Bazett by the QT-RR scatter plots of Figure S1B (Results).",
        "Centering reference: the paper centers on the overall mean of all baseline QTcF values ('BQTcF-mean BQTcF', Table S1; 'QTcF.mB' in the supplementary code) but does not print that mean or any baseline QTcF summary. The model uses the rounded 400 ms. A centering error of d ms shifts every prediction by 0.052 * d ms and cancels exactly in the placebo-corrected DeltaDeltaQTcF."
      ),
      source_name = "BQTCF (baseline QTcF)"
    ),
    NTIME = list(
      description = "Nominal post-dose dECG time point in hours. Levels 0 (pre-dose; reference), 0.5, 1, 2, 3, 3.5, 4, 6, 8, 12, 24, 36 and 48 h; each non-reference level carries its own intercept shift.",
      units = "h",
      type = "categorical",
      reference_category = "0 (pre-dose)",
      notes = paste(
        "Sundell 2025 Equation 1 term theta3 * TIME_l, 'the fixed effect associated with time', estimated at each nominal time 'regardless of its significance, to minimize the potential confounding effect of any diurnal change' (Methods 'C-QTc model'); entered as the factor NOMTIME in the supplementary lme4 code. Table S1 lists 'Time = 0' as the reference group.",
        "Nominal dECG/PK times were pre-dose and 0.5, 1, 2, 3, 4, 6, 8, 12, 24, 36 and 48 h, plus 3.5 h for the 800 and 1200 mg fractionated arms (Methods 'ECG and PK sampling').",
        "model() bins the supplied value to the nearest scheduled level (cut points 0.25, 0.75, 1.5, 2.5, 3.25, 3.75, 5, 7, 10, 18, 30 and 42 h).",
        "Supply the NOMINAL time, not rxode2's tad()."
      ),
      source_name = "TIME_l / NOMTIME"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 62L,
    n_studies = 1L,
    n_observations = 653L,
    age_range = "18-50 years (inclusion criterion)",
    weight_range = "50-100 kg (inclusion criterion)",
    sex_female_pct = 0,
    disease_state = "Healthy male participants (BMI 18-30 kg/m^2).",
    dose_range = paste(
      "Single oral suspension doses of balcinrenone 5, 25, 100, 200, 400",
      "and 800 mg, and fractionated doses of 800 and 1200 mg given as",
      "200 mg/h, or placebo, under fasted conditions (6 active and 2",
      "placebo per group; 5 active for 200 and 800 mg)."
    ),
    regions = NA_character_,
    notes = paste(
      "Phase I single ascending dose study D6400C00001 (NCT02484729), Part A: randomized, single-blind, placebo-controlled, sequential-group design, prospectively planned as a thorough-QT study substitute (Methods 'Study design').",
      "Twelve time-matched PK samples and 12-lead continuous dECG recordings per participant (pre-dose to 48 h); 653 time-matched PK/dECG pairs entered the model, and 61 post-dose samples after 5-200 mg were below the 10 nmol/L quantitation limit (Methods 'ECG and PK sampling').",
      "Fitted with lme4 1.1-7 / lmerTest in R 3.5.1 (Methods 'Analysis'; supplementary code).",
      "No baseline demographics or baseline QTcF/HR summaries are reported in the paper; age and weight are the inclusion ranges.",
      "Model-predicted DeltaDeltaQTcF (90% CI) was 0.35 (-1.29, 2.00) ms at 2.16 umol/L and 0.31 (-2.26, 2.88) ms at 16.70 umol/L (Table 1). A sensitivity analysis excluding the 1200 mg fractionated arm (Table S1, lower block) gave a slope of -0.140 ms per umol/L; it is not packaged as a separate model."
    )
  )

  ini({
    # Sundell 2025 Equation 1 (prespecified linear mixed-effects model):
    #   DeltaQTcF_ijl = theta0 + h0_i + theta1 * TRT_j + (theta2 + h2_i) * C_ijl
    #                 + theta3 * TIME_l + theta4 * (QTcF_ijl=0 - mean QTcF0)
    # lme4 code: DQTCF ~ NOMTIME + TRTID + I(BQTCF - QTcF.mB) + DV + (DV | SID).
    # All values are the primary-analysis rows of Supplementary Table S1
    # ('Fixed effects parameter estimates in the final model'). Every
    # parameter is on the linear ms scale.

    e0 <- -0.267
    label("Intercept theta0 on DeltaQTcF: placebo, pre-dose, mean baseline (ms)")
    # Table S1 'theta0 (intercept)' = -0.267 (95% CI -2.283, 1.750)

    e_on_treatment_e0 <- 0.359
    label("Treatment-specific intercept theta1 for the balcinrenone arm (ms)")
    # Table S1 'theta1 (active arm)' = 0.359 (95% CI -1.648, 2.366), reference 'Placebo arm'

    slope <- -0.003
    label("Concentration-DeltaQTcF slope theta2 (ms per umol/L)")
    # Table S1 'theta2 (slope concentration)' = -0.003 (95% CI -0.181, 0.176); Abstract 'ms/umol/L'

    # Nominal-time intercept shifts theta3 against the 'Time = 0' reference.
    e_ntime0p5_e0 <- -2.908
    label("Effect of the 0.5 h nominal time point on the intercept (ms)")
    # Table S1 'theta3 (time = 0.5)' = -2.908 (95% CI -4.579, -1.238)

    e_ntime1_e0 <- -0.822
    label("Effect of the 1 h nominal time point on the intercept (ms)")
    # Table S1 'theta3 (time = 1)' = -0.822 (95% CI -2.503, 0.858)

    e_ntime2_e0 <- -0.644
    label("Effect of the 2 h nominal time point on the intercept (ms)")
    # Table S1 'theta3 (time = 2)' = -0.644 (95% CI -2.315, 1.028)

    e_ntime3_e0 <- -0.671
    label("Effect of the 3 h nominal time point on the intercept (ms)")
    # Table S1 'theta3 (time = 3)' = -0.671 (95% CI -2.341, 1.000)

    e_ntime3p5_e0 <- 0.585
    label("Effect of the 3.5 h nominal time point on the intercept (ms)")
    # Table S1 'theta3 (time = 3.5)' = 0.585 (95% CI -2.315, 3.484); fractionated arms only

    e_ntime4_e0 <- 1.228
    label("Effect of the 4 h nominal time point on the intercept (ms)")
    # Table S1 'theta3 (time = 4)' = 1.228 (95% CI -0.454, 2.911)

    e_ntime6_e0 <- -0.658
    label("Effect of the 6 h nominal time point on the intercept (ms)")
    # Table S1 'theta3 (time = 6)' = -0.658 (95% CI -2.271, 0.954)

    e_ntime8_e0 <- -2.997
    label("Effect of the 8 h nominal time point on the intercept (ms)")
    # Table S1 'theta3 (time = 8)' = -2.997 (95% CI -4.611, -1.383)

    e_ntime12_e0 <- -1.475
    label("Effect of the 12 h nominal time point on the intercept (ms)")
    # Table S1 'theta3 (time = 12)' = -1.475 (95% CI -3.153, 0.203)

    e_ntime24_e0 <- -2.011
    label("Effect of the 24 h nominal time point on the intercept (ms)")
    # Table S1 'theta3 (time = 24)' = -2.011 (95% CI -3.730, -0.292)

    e_ntime36_e0 <- -1.774
    label("Effect of the 36 h nominal time point on the intercept (ms)")
    # Table S1 'theta3 (time = 36)' = -1.774 (95% CI -3.505, -0.044)

    e_ntime48_e0 <- -2.128
    label("Effect of the 48 h nominal time point on the intercept (ms)")
    # Table S1 'theta3 (time = 48)' = -2.128 (95% CI -3.920, -0.336)

    e_qtc_bl_e0 <- -0.052
    label("Effect of centered baseline QTcF on the intercept, theta4 (ms per ms)")
    # Table S1 'theta4 (BQTcF-mean BQTcF)' = -0.052 (95% CI -0.108, 0.003)

    # Random effects h0_i (intercept) and h2_i (slope) with an unstructured
    # covariance ('(DV | SID)' in the supplementary code; Methods Equation 1
    # legend). Their variances, covariance and the residual SD are not
    # reported anywhere in the paper or supplements (Table S1 lists fixed
    # effects only), so all three are fixed to zero and the model gives
    # typical-value predictions.
    etae0 ~ fixed(0)
    etaslope ~ fixed(0)

    addSd <- fixed(0)
    label("Additive residual error SD on DeltaQTcF (ms; zero because the source does not report it)")
  })

  model({
    # Nominal-time class effect: bin NTIME to the nearest scheduled level;
    # NTIME < 0.25 is the pre-dose 0 h reference (all indicators 0).
    nt0p5 <- (NTIME >= 0.25) * (NTIME < 0.75)
    nt1 <- (NTIME >= 0.75) * (NTIME < 1.5)
    nt2 <- (NTIME >= 1.5) * (NTIME < 2.5)
    nt3 <- (NTIME >= 2.5) * (NTIME < 3.25)
    nt3p5 <- (NTIME >= 3.25) * (NTIME < 3.75)
    nt4 <- (NTIME >= 3.75) * (NTIME < 5)
    nt6 <- (NTIME >= 5) * (NTIME < 7)
    nt8 <- (NTIME >= 7) * (NTIME < 10)
    nt12 <- (NTIME >= 10) * (NTIME < 18)
    nt24 <- (NTIME >= 18) * (NTIME < 30)
    nt36 <- (NTIME >= 30) * (NTIME < 42)
    nt48 <- (NTIME >= 42)
    ntime_effect <-
      e_ntime0p5_e0 * nt0p5 + e_ntime1_e0 * nt1 + e_ntime2_e0 * nt2 +
      e_ntime3_e0 * nt3 + e_ntime3p5_e0 * nt3p5 + e_ntime4_e0 * nt4 +
      e_ntime6_e0 * nt6 + e_ntime8_e0 * nt8 + e_ntime12_e0 * nt12 +
      e_ntime24_e0 * nt24 + e_ntime36_e0 * nt36 + e_ntime48_e0 * nt48

    # Baseline QTcF centered on the overall mean baseline, which the paper
    # does not print; 400 ms is the rounded standard (see covariateData).
    qtc_bl_ref <- 400
    qtc_bl_centered <- QTC_BL - qtc_bl_ref

    e0_i <- e0 + etae0 +
      e_on_treatment_e0 * ON_TREATMENT +
      e_qtc_bl_e0 * qtc_bl_centered +
      ntime_effect
    slope_i <- slope + etaslope

    # Change from baseline in QTcF (ms)
    QTcF <- e0_i + slope_i * CP_BALCINRENONE_UM
    QTcF ~ add(addSd)
  })
}
