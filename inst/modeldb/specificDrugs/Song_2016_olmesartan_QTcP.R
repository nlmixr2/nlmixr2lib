Song_2016_olmesartan_QTcP <- function() {
  description <- paste(
    "Linear mixed-effects concentration-QTc model for olmesartan (the active",
    "moiety of olmesartan medoxomil) in 55 healthy adults of a single-dose",
    "thorough QTc crossover study (olmesartan medoxomil 40 mg and 160 mg,",
    "placebo, moxifloxacin 400 mg). The endpoint is the baseline-adjusted,",
    "time-matched placebo-adjusted population-corrected QT interval",
    "(DeltaDeltaQTcP, ms; QTcP = QT / RR^(1/3.265)): DeltaDeltaQTcP = e0 * (1 + etae0) +",
    "e_rr_delta_e0 * (1 + etae_rr_delta_e0) * RR_DELTA + slope * (1 +",
    "etaslope) * CP_OLMESARTAN_NGML + eps, where RR_DELTA is the",
    "baseline-adjusted, placebo-adjusted RR interval (DeltaDeltaRR, ms) and",
    "CP_OLMESARTAN_NGML the time-matched olmesartan plasma concentration",
    "(ng/mL). Intercept -0.156 ms, olmesartan slope 0.000800 ms/(ng/mL),",
    "DeltaDeltaRR coefficient 0.00646 ms/ms. Proportional random effects on",
    "all three coefficients (SDs of the multiplier 4.63, 3.31 and 5.58) and an",
    "additive residual SD of 3.99 ms. PD-only model: the concentration and",
    "DeltaDeltaRR are supplied as time-varying covariates. Companion models",
    "from the same paper: Song_2016_olmesartan_QTc (uncorrected QT) and",
    "Song_2016_olmesartan_QTcF (Fridericia-corrected QT).",
    sep = " "
  )

  reference <- paste(
    "Song S, Matsushima N, Lee J, Mendell J.",
    "Linear mixed-effects model of QTc prolongation for olmesartan",
    "medoxomil. J Clin Pharmacol. 2016;56(1):96-100.",
    "doi:10.1002/jcph.572.",
    sep = " "
  )

  vignette <- "Song_2016_olmesartan_QTc"

  units <- list(
    time = "h",
    dosing = "(none; PD-only model fed by an external olmesartan plasma-concentration covariate)",
    concentration = "(observation QTcP is the baseline-adjusted, time-matched placebo-adjusted population-corrected QT interval, DeltaDeltaQTcP (QTcP = QT / RR^(1/3.265)), in ms; driving covariates CP_OLMESARTAN_NGML in ng/mL and RR_DELTA in ms)"
  )

  covariateData <- list(
    CP_OLMESARTAN_NGML = list(
      description = "Olmesartan plasma concentration at the time of each ECG observation, supplied as a time-varying covariate from observed plasma samples or an upstream PK model.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters as the linear term slope * (1 + etaslope) * CP_OLMESARTAN_NGML; Song 2016 Table 2 reports the slope in ms per ng/mL, so no in-model unit rescaling is required.",
        "In the source analysis the concentration was measured by LC-MS/MS in samples taken at the ECG time points (0.5 to 24 h postdose) after single oral doses of 40 mg and 160 mg olmesartan medoxomil.",
        "Mean (SD) Cmax was 863.3 (268.7) ng/mL after 40 mg and 2440 (678) ng/mL after 160 mg (companion thorough QTc paper, Mendell 2016, doi:10.1002/jcph.610)."
      ),
      source_name = "OMconc"
    ),
    RR_DELTA = list(
      description = "Baseline-adjusted, time-matched placebo-adjusted RR interval (DeltaDeltaRR): the change of RR from the mean of three predose values in the olmesartan period, minus the same change at the matching time in the placebo period of the same subject.",
      units = "ms",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "RR = 60000 / heart rate (ms), from digitally extracted Holter ECGs.",
        "Enters as the linear term e_rr_delta_e0 * (1 + etae_rr_delta_e0) * RR_DELTA, a regression-based heart-rate correction inside the model; RR_DELTA = 0 (no drug-induced change in RR) gives the heart-rate-neutral prediction.",
        "Observed range roughly -350 to +300 ms (Song 2016 Figure 1, panels 4-6). Positive values are RR lengthening (heart-rate slowing)."
      ),
      source_name = "DeltaDeltaRR_ij"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 55L,
    n_studies = 1L,
    age_range = "19-45 years (mean 31.9)",
    weight_range = "not reported; BMI 19.1-30.6 kg/m^2 (mean 26.2)",
    sex_female_pct = 50,
    race_ethnicity = c(White = 91.1, Black = 7.1, AmericanIndianAlaskaNative = 1.8),
    disease_state = "Healthy volunteers.",
    dose_range = paste(
      "Single oral doses of olmesartan medoxomil 40 mg (highest approved",
      "dose) and 160 mg (supratherapeutic), each in its own period of a",
      "4-period crossover with placebo and moxifloxacin 400 mg."
    ),
    regions = "not reported (single centre)",
    notes = paste(
      "Thorough QTc study per ICH E14; 56 subjects (28 men, 28 women) randomised, 51 completed; 55 in the cardiodynamic and pharmacokinetic analysis sets.",
      "Hispanic/Latino 85.7%. Demographics from Song 2016 Methods and the companion study report (Mendell 2016 Table 1, doi:10.1002/jcph.610).",
      "Triplicate ECGs extracted from Holter recordings predose and at 0.5, 1, 1.5, 2, 2.5, 3, 4, 6, 8, 10, 12 and 24 h postdose (Mendell 2016); the concentration-QTc model was fitted to the two olmesartan treatments in NONMEM 7.1.0 (FOCE with interaction)."
    )
  )

  ini({
    # Song 2016 Equation 1 (Results 'Exposure-Response Analysis'):
    #   ddQTc_ij = th1*(1 + w1) + th2*ddRR_ij*(1 + w2) + th3*OMconc_ij*(1 + w3)
    # with the additive residual of Equation 2. All values from Table 2,
    # QTcP columns. The table's 'omega2' rows are signed (here -5.58 for ddRR),
    # so they are the SDs of the proportional random effects (a NONMEM THETA
    # multiplying a unit-variance ETA), not variances; the variance of each
    # proportional eta below is the square of the printed value.

    e0 <- -0.156
    label("Intercept of DeltaDeltaQTcP at zero concentration and zero DeltaDeltaRR (ms)")
    # Table 2 'Intercept (ms)', QTcP = -0.156 (SE 0.134)

    slope <- 0.000800
    label("Olmesartan concentration-DeltaDeltaQTcP slope (ms per ng/mL)")
    # Table 2 'theta3 OMconc (ms/[ng/mL])', QTcP = 0.000800 (SE 0.000381); also Results text

    e_rr_delta_e0 <- 0.00646
    label("Coefficient of DeltaDeltaRR on DeltaDeltaQTcP (ms per ms)")
    # Table 2 'theta ddRR (ms)', QTcP = 0.00646 (SE 0.00538)

    etae0 ~ 21.4369
    # Table 2 'omega2 Intercept', QTcP = 4.63 (SE 4.02), proportional SD; 4.63^2 = 21.4369
    etaslope ~ 10.9561
    # Table 2 'omega2 beta OMconc', QTcP = 3.31 (SE 1.52), proportional SD; 3.31^2 = 10.9561
    etae_rr_delta_e0 ~ 31.1364
    # Table 2 'omega2 ddRR', QTcP = -5.58 (SE 4.74), proportional SD; 5.58^2 = 31.1364

    addSd <- 3.99
    label("Additive residual error SD on DeltaDeltaQTcP (ms)")
    # Table 2 'sigma2', QTcP = 3.99 (SE 0.209), read as an SD like the omega rows
  })

  model({
    # Individual coefficients, proportional random effects (Equation 1)
    e0_i <- e0 * (1 + etae0)
    slope_i <- slope * (1 + etaslope)
    rr_coef_i <- e_rr_delta_e0 * (1 + etae_rr_delta_e0)

    # Baseline-adjusted, placebo-adjusted QTcP (ms)
    QTcP <- e0_i + rr_coef_i * RR_DELTA + slope_i * CP_OLMESARTAN_NGML
    QTcP ~ add(addSd)
  })
}
