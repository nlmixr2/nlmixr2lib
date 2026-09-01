Mukker_2026_tuvusertib_HR <- function() {
  description <- paste(
    "Garnett-type linear mixed-effects concentration-heart-rate model",
    "for the investigational ATR kinase inhibitor tuvusertib (M1774) in",
    "patients with advanced solid tumors enrolled in Part A1 of the",
    "phase I first-in-human study DDRiver Solid Tumors 301",
    "(NCT04170153; dose range 5-270 mg). The endpoint is the change",
    "from baseline in heart rate (DeltaHR, bpm). The model shares the",
    "structure of the companion C-DeltaQTcF model and adds a dosing-day",
    "term:",
    "DeltaHR = (theta0 + eta0) + (theta1 + eta1) * C",
    "+ theta2(nominal time) + theta_day8 * DAY8",
    "+ theta3 * (HR - 70) + eps,",
    "with theta0 = 3.00 bpm, theta1 = 0.00111 bpm/(ng/mL) whose 95% CI",
    "includes zero, four nominal post-dose-time intercept shifts",
    "(-1.19, -5.16, 0.694, 2.80 bpm at 0, 1, 2 and 3 h), a day-8-versus-",
    "day-1 shift of 2.84 bpm, and a centered-baseline-HR coefficient of",
    "-0.0214 bpm/bpm. Random effects sit additively on both the",
    "intercept (SD 4.43 bpm) and the concentration slope (SD 0.00173",
    "bpm/(ng/mL)) with a correlation of 0.400; residual SD is 5.40 bpm.",
    "This model exists to establish the FIRST of the four Garnett",
    "assumptions underlying the companion QTc analysis -- that the drug",
    "has no effect on heart rate. Its purpose is a negative result: the",
    "slope CI spans zero and the predicted DeltaHR at the maximum",
    "observed concentration is below the 10 bpm threshold of clinical",
    "relevance. PD-only model: tuvusertib plasma concentration is",
    "supplied as a time-varying covariate CP_TUVUSERTIB_NGML (ng/mL).",
    "Companion model files are Mukker_2026_tuvusertib_QTcF.R and",
    "Mukker_2026_tuvusertib_hERG.R.",
    sep = " "
  )

  reference <- paste(
    "Mukker JK, Yap TA, Tolcher AW, de Bono JS, Plummer R, Grosser G,",
    "van Amsterdam C, Schieferstein H, Witjes H, Diderichsen PM,",
    "Krebs-Brown A, Gao W, Strotmann R, Szucs Z, Gounaris I,",
    "Venkatakrishnan K. An Integrated Nonclinical and Clinical Risk",
    "Assessment of the Effects of Investigational ATRi Tuvusertib on QTc",
    "Interval in Patients With Solid Tumors.",
    "Clinical and Translational Science 2026;19(2):e70496.",
    "doi:10.1111/cts.70496.",
    "Parameter estimates for this model are in Supporting Information",
    "Table S1 (Data S1, cts70496-sup-0001-DataS1.docx).",
    sep = " "
  )

  vignette <- "Mukker_2026_tuvusertib_QTc"

  units <- list(
    time          = "h",
    dosing        = "(none; PD-only model fed by an external tuvusertib plasma-concentration covariate)",
    concentration = "(observation dHR is the change from baseline in heart rate, bpm; driving covariate CP_TUVUSERTIB_NGML is in ng/mL)"
  )

  covariateData <- list(
    CP_TUVUSERTIB_NGML = list(
      description        = "Instantaneous tuvusertib plasma concentration at the time of each ECG observation, supplied as a time-varying covariate from observed plasma samples or an upstream PK source.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying per event row. Drives the linear concentration-DeltaHR term (theta1 + eta1) * CP_TUVUSERTIB_NGML.",
        "Same time-matched PK-ECG dataset as the companion C-DeltaQTcF model (Mukker 2026 Methods 2.2.2); concentrations measured by validated LC/MS, LLOQ 0.5 ng/mL.",
        "The slope is reported directly in bpm per ng/mL (Mukker 2026 Table S1), so no in-model unit rescaling is needed.",
        "Reference values observed, used as the Table S2 validation grid: median 524 ng/mL, P90 1732 ng/mL, P95 2252 ng/mL, maximum 3290 ng/mL.",
        "Set to 0 outside the drug-exposure window."
      ),
      source_name        = "tuvusertib plasma concentration"
    ),
    HR = list(
      description        = "Subject's baseline (pre-dose, time-zero) heart rate, treated as a per-subject time-fixed covariate. Enters the linear-mixed-effects intercept as the centered term e_hr_bl_e0 * (HR - 70). Set HR = 70 bpm for the typical subject -- the centered term then collapses to 0.",
      units              = "beats/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. This is the BASELINE reading of the canonical HR covariate, not an observation-time vital sign; the HR register entry directs that baseline-versus-time-varying status be documented per-model, which is what this note does.",
        "Centering reference: the paper does NOT report the cohort mean baseline heart rate (Table 1 gives baseline QTcF but no baseline HR, and no HR summary appears in the supplement). Per the standing policy on undefined centering values, the model uses the rounded clinical standard hr_bl_ref = 70 bpm. This is an ASSUMPTION and is recorded in the vignette Errata.",
        "The assumption is low-impact: the coefficient is -0.0214 bpm/bpm with a 95% CI (-0.148, 0.105) that comfortably spans zero, so a mis-specified centering constant shifts the intercept by at most a fraction of a bpm over any plausible cohort mean (e.g. a true mean of 80 bpm would shift the typical-value prediction by 0.214 bpm). None of the paper's reported validation targets (Table S2) depend on this term, because they tabulate the drug effect corrected for intercept, baseline and time.",
        "Note that the model's observable is named dHR rather than HR precisely so that it does not shadow this covariate column."
      ),
      source_name        = "Baseline HR"
    ),
    T_LASTDOSE = list(
      description        = "Nominal (protocol-scheduled) time after the most recent tuvusertib dose at which the triplicate ECG was recorded, in hours. Selects one of the four estimated nominal-timepoint intercept shifts.",
      units              = "h",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying per event row. Mukker 2026 Table S1 reports four 'Nominal time after dose' estimates (0, 1, 2, 3 h), each footnoted as 'the estimated difference from the population mean intercept' -- mean-centered class effects, not contrasts against an omitted reference level.",
        "Because this is a PD-only model with no dosing records, rxode2's native tad() is unavailable and the post-dose clock is supplied as a covariate column; see the T_LASTDOSE register entry.",
        "Carries the NOMINAL scheduled hour. model() bins the supplied value to its nearest scheduled level (< 0.5 h -> 0 h; 0.5-1.5 h -> 1 h; 1.5-2.5 h -> 2 h; >= 2.5 h -> 3 h).",
        "The nominal-time effects are markedly larger here than in the companion QTcF model (a 7.96 bpm spread from the 1 h to the 3 h timepoint versus a 4.46 ms spread for QTcF), consistent with a real within-day heart-rate rhythm that is independent of drug concentration."
      ),
      source_name        = "nominal time after dose"
    ),
    DAY8 = list(
      description        = "Dosing-day landmark indicator: 1 = the ECG was recorded on day 8 of treatment, 0 = on day 1.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (day 1 of treatment)",
      notes              = paste(
        "Time-varying within subject. Mukker 2026 Table S1 row 'Dosing day (Day 8 versus Day 1), bpm' = 2.84 (95% CI 1.78, 3.91): heart rate is 2.84 bpm higher on day 8 than on day 1 at matched concentration and matched nominal post-dose time.",
        "Day 1 and day 8 are the two intensive PK-ECG sampling days of DDRiver Solid Tumors 301 Part A1 (Mukker 2026 Figure S1).",
        "This term appears ONLY in the C-DeltaHR model; the C-DeltaQTcF model in Table 2 has no dosing-day term, so the companion model file does not declare this covariate.",
        "The CI excludes zero, so the day effect is statistically resolved -- but it is a time-on-study effect, not a concentration effect, and therefore does not bear on the Garnett assumption being tested."
      ),
      source_name        = "Dosing day (Day 8 versus Day 1)"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 55L,
    n_studies        = 1L,
    age_range        = "mean (SD) 61.9 (10.9) years",
    weight_range     = "mean (SD) 78.6 (18.0) kg",
    sex_female_pct   = 58.2,
    race_ethnicity   = c(White = 76.4, Asian = 9.1, Black_or_African_American = 3.6, Other = 10.9),
    disease_state    = paste(
      "Patients with advanced solid tumors enrolled in Part A1 of DDRiver",
      "Solid Tumors 301 (NCT04170153), an open-label, first-in-human,",
      "multicenter phase I dose-escalation study of tuvusertib monotherapy.",
      "Same analysis population as the companion C-DeltaQTcF model",
      "(Mukker 2026 Table 1); the paper does not report a separate",
      "subject count for the C-DeltaHR analysis, which was run on the",
      "same time-matched PK-ECG dataset."
    ),
    dose_range       = paste(
      "Tuvusertib 5-270 mg orally: 5, 10, 20, 40, 80, 130, 180, 220 and",
      "270 mg once daily continuous; 180 and 220 mg QD 2 weeks on / 1",
      "week off; and 150 mg twice daily 4 days on / 3 days off."
    ),
    regions          = NA_character_,
    notes            = paste(
      "This model was developed to test the FIRST Garnett assumption -- 'there is no effect of tuvusertib on HR' -- after an exploratory analysis of the concentration-RR-interval relationship suggested a slightly positive trend in the C-DeltaHR relationship at higher concentrations (Mukker 2026 Supporting Information Results S1, Figure S2a).",
      "The conclusion is a negative result: the slope estimate 0.00111 bpm/(ng/mL) has a 95% CI (-0.000250, 0.00247) that includes zero, 'indicating the lack of a discernible tuvusertib plasma concentration effect on HR' (Results S1).",
      "Validation targets: predicted DeltaHR (90% CI) was 0.581 (-0.0168, 1.18) bpm at the median observed concentration of 524 ng/mL, 1.92 (-0.0554, 3.90) at P90 = 1732 ng/mL, 2.50 (-0.0720, 5.06) at P95 = 2252 ng/mL, and 3.65 (-0.105, 7.40) at the maximum observed 3290 ng/mL -- the upper bound at the maximum concentration being below the 10 bpm threshold of clinical relevance (Mukker 2026 Table S2).",
      "Parameter estimates were obtained by restricted maximum likelihood (REML), stated explicitly for this model in Results S1 (the main-paper C-DeltaQTcF model does not name its estimation method). Software was R 4.3.1.",
      "The second Garnett assumption was checked separately from baseline ECG data: uncorrected QT depended on heart rate, and Fridericia correction removed that dependence to non-significance (p = 0.11), so QTcF is independent of HR (Results S1, Figure S2b)."
    )
  )

  ini({
    # ==================================================================
    # Garnett linear mixed-effects concentration-heart-rate model.
    # Same structure as Mukker 2026 Equation 1 (the C-DeltaQTcF model),
    # with DeltaHR substituted for DeltaQTcF, baseline HR substituted
    # for baseline QTcF, and an extra additive dosing-day term:
    #
    #   DeltaHR_ij = (theta0 + eta0_i)
    #              + (theta1 + eta1_i) * C_ij
    #              + theta2 * Time_j
    #              + theta_day8 * DAY8_ij
    #              + theta3 * (HR_i0 - HR0)
    #              + eps_ij
    #
    # All values below are from Mukker 2026 Supporting Information
    # Table S1 ('Parameter estimates of the final C-DeltaHR model').
    # ==================================================================

    e0 <- 3.00
    label("Population mean intercept theta0 on DeltaHR (bpm)")
    # Mukker 2026 Table S1 'Mean intercept, bpm' = 3.00 (95% CI 1.49,
    # 4.50). Not log-transformed: a change-from-baseline intercept is
    # signed. (It happens to be positive here, unlike the QTcF model's
    # -0.125, but the encoding decision follows the endpoint's nature,
    # not the sign of one estimate.)

    slope <- 0.00111
    label("Concentration-DeltaHR slope theta1 (bpm per ng/mL)")
    # Mukker 2026 Table S1 'Tuvusertib plasma concentration slope,
    # bpm/(ng/mL)' = 0.00111 (95% CI -0.000250, 0.00247). The CI
    # INCLUDES zero -- this is the paper's evidence for the first
    # Garnett assumption (no drug effect on heart rate), quoted in the
    # main text at Results 3.2.3.
    #
    # Not log-transformed, for the same structural reason as the
    # companion QTcF model and more forcefully so: the additive slope
    # eta (SD 0.00173) is larger than the typical value (0.00111), and
    # the typical value itself is not distinguishable from zero. A
    # log-transform would be meaningless here.

    # ------------------------------------------------------------------
    # Nominal-post-dose-time intercept shifts. Mukker 2026 Table S1
    # rows 'Nominal time after dose (0h/1h/2h/3h)', each footnoted
    # (footnote a): 'The fixed-effect estimates for nominal time after
    # the dose represent the estimated difference from the population
    # mean intercept.'
    # ------------------------------------------------------------------

    e_tad0h_e0 <- -1.19
    label("Effect of the 0 h nominal post-dose timepoint on e0 (additive, bpm)")
    # Mukker 2026 Table S1 'Nominal time after dose (0 h), bpm' = -1.19
    # (95% CI -2.56, 0.181).

    e_tad1h_e0 <- -5.16
    label("Effect of the 1 h nominal post-dose timepoint on e0 (additive, bpm)")
    # Mukker 2026 Table S1 'Nominal time after dose (1 h), bpm' = -5.16
    # (95% CI -6.13, -4.19). CI excludes zero.

    e_tad2h_e0 <- 0.694
    label("Effect of the 2 h nominal post-dose timepoint on e0 (additive, bpm)")
    # Mukker 2026 Table S1 'Nominal time after dose (2 h), bpm' = 0.694
    # (95% CI -0.352, 1.74).

    e_tad3h_e0 <- 2.80
    label("Effect of the 3 h nominal post-dose timepoint on e0 (additive, bpm)")
    # Mukker 2026 Table S1 'Nominal time after dose (3 h), bpm' = 2.80
    # (95% CI 1.82, 3.79). CI excludes zero.

    e_day8_e0 <- 2.84
    label("Effect of dosing day 8 versus day 1 on e0 (additive, bpm)")
    # Mukker 2026 Table S1 'Dosing day (Day 8 versus Day 1), bpm' =
    # 2.84 (95% CI 1.78, 3.91). Present in the C-DeltaHR model only;
    # the C-DeltaQTcF model of Table 2 carries no dosing-day term.

    e_hr_bl_e0 <- -0.0214
    label("Effect of centered baseline HR on e0 (additive, bpm per bpm)")
    # Mukker 2026 Table S1 'Baseline HR, bpm' = -0.0214 (95% CI -0.148,
    # 0.105). Negative, matching the regression-to-the-mean sign of the
    # baseline-QTcF coefficient in the companion model, but with a CI
    # that spans zero.

    # ==================================================================
    # Random effects. Additive on both intercept and slope, matching
    # the companion C-DeltaQTcF model's structure (Mukker 2026 Methods
    # 2.2.3, applied to the DeltaHR endpoint per Results S1).
    #
    # Table S1 reports SDs plus a correlation:
    #   SD of intercept                              = 4.43    bpm
    #   SD of tuvusertib plasma concentration effect = 0.00173 bpm/(ng/mL)
    #   Correlation between the two                  = 0.400
    # Converted to the lower-triangular variance-covariance block:
    #   var(eta0)      = 4.43^2                = 19.6249
    #   cov(eta0,eta1) = 0.400 * 4.43 * 0.00173 = 3.065560e-03
    #   var(eta1)      = 0.00173^2             = 2.9929e-06
    # Positive-definiteness check: cov^2 = 9.398e-06 versus
    # var(eta0) * var(eta1) = 5.873e-05; their ratio is 0.16 = 0.400^2
    # exactly, as it must be, confirming the conversion.
    #
    # Note the contrast with the companion QTcF model, whose
    # intercept-slope correlation is -0.00273 (negligible). Here the
    # correlation is a substantial +0.400.
    # ==================================================================
    etae0 + etaslope ~ c(
      19.6249,
      3.065560e-03, 2.9929e-06
    )

    # ==================================================================
    # Residual unexplained variability, additive on the DeltaHR (bpm)
    # scale. Mukker 2026 Table S1 'SD of residual variability, bpm'
    # = 5.40.
    # ==================================================================
    addSd <- 5.40
    label("Additive residual error standard deviation on DeltaHR (bpm)")
  })

  model({
    # ==================================================================
    # 1. Per-individual intercept and slope; both etas additive.
    # ==================================================================
    e0_i    <- e0 + etae0
    slope_i <- slope + etaslope

    # ==================================================================
    # 2. Nominal-post-dose-time class effect, binned to the nearest
    #    scheduled level (see the companion QTcF model for the full
    #    rationale). Comparisons evaluate to 0/1 in rxode2 and the four
    #    indicators are mutually exclusive and exhaustive.
    # ==================================================================
    i_tad0h <- (T_LASTDOSE < 0.5)
    i_tad1h <- (T_LASTDOSE >= 0.5) * (T_LASTDOSE < 1.5)
    i_tad2h <- (T_LASTDOSE >= 1.5) * (T_LASTDOSE < 2.5)
    i_tad3h <- (T_LASTDOSE >= 2.5)

    tad_effect <-
      e_tad0h_e0 * i_tad0h +
      e_tad1h_e0 * i_tad1h +
      e_tad2h_e0 * i_tad2h +
      e_tad3h_e0 * i_tad3h

    # ==================================================================
    # 3. Centered baseline-HR term. The cohort mean baseline heart rate
    #    is NOT reported anywhere in Mukker 2026 or its supplement, so
    #    the rounded clinical standard 70 bpm is used as the centering
    #    reference per the standing policy on undefined centering
    #    values. Documented in the vignette Errata. Set HR = 70 to make
    #    the term vanish for a typical subject.
    # ==================================================================
    hr_bl_ref      <- 70
    hr_bl_centered <- HR - hr_bl_ref

    # ==================================================================
    # 4. Predicted change from baseline in heart rate (bpm). The
    #    observable is named dHR, not HR, so that it does not shadow
    #    the baseline-HR covariate column of the same canonical name.
    #
    #    As with the companion QTcF model, Mukker 2026 Table S2
    #    tabulates the DRUG effect only (slope * C); recover it by
    #    differencing this model's prediction against its own
    #    prediction at zero concentration.
    # ==================================================================
    dHR <- e0_i + slope_i * CP_TUVUSERTIB_NGML + tad_effect +
      e_day8_e0 * DAY8 + e_hr_bl_e0 * hr_bl_centered

    dHR ~ add(addSd)
  })
}
