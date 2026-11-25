Minocha_2019_placebo_hr <- function() {
  description <- "Placebo (no-drug) heart-rate model for healthy adult volunteers in AbbVie phase I trials (Minocha 2019): subject baseline at midnight with an additive female-sex shift, plus three additive cosine circadian components (24-h, 12-h and 6-h periods) and a linear time drift, with additive between-subject variability on the baseline, the 24-h and 12-h amplitudes, the 24-h and 6-h acrophases and the drift slope, and a proportional residual error. Time is clock time in hours since the midnight before the subject's first recorded measurement. Intended as a placebo backbone for exposure-response analyses of heart rate."
  reference <- paste(
    "Minocha M, Li H, Chiu YL, Carter D, Othman AA.",
    "Models of Variability and Circadian Rhythm in Heart Rate, Blood Pressure,",
    "and QT Interval for Healthy Subjects Who Received Placebo in Phase I Trials.",
    "Clin Transl Sci. 2019;12(5):470-480. doi:10.1111/cts.12640.",
    "Corrigendum: Clin Transl Sci. 2020;13(3):638. doi:10.1111/cts.12773",
    "(Table S1, between-subject variability estimates missing from Table 2).",
    sep = " "
  )
  vignette <- "Minocha_2019_placebo_cardiovascular"
  units <- list(
    time = "h",
    dosing = "n/a (placebo baseline model with no drug input)",
    concentration = "beats/min (heart rate)"
  )

  covariateData <- list(
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Additive shift on the midnight baseline heart rate (Minocha 2019 Eq. 9 and Table 2). The authors' NONMEM SEX column enters the baseline as THETA(6)*SEX with males the most frequent (reference, value 0) category (Methods, Eq. 3: CATi = 0 for the most frequent category).",
      source_name = "SEX"
    )
  )

  covariatesDataExcluded <- list(
    RACE_BLACK = list(
      description = "Black race indicator (1 = Black, 0 = other)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White, Asian or other)",
      notes = "Screened as an additive baseline shift (Methods, Eq. 3) and not significant for heart rate (Table 2, 'NS').",
      source_name = "RACE"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator (1 = Asian, 0 = other)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White, Black or other)",
      notes = "Screened as an additive baseline shift (Methods, Eq. 3) and not significant for heart rate (Table 2, 'NS').",
      source_name = "RACE"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened as a power model normalised to 30 years (Methods, Eq. 3) and not significant for heart rate (Table 2, 'NS').",
      source_name = "AGE"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened as a power model normalised to 70 kg (Methods, Eq. 3) and not significant for heart rate (Table 2, 'NS').",
      source_name = "WTKG"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 996L,
    n_studies = 63L,
    n_observations = 16509L,
    age_range = "adults; mean 34 years (SD 10); range not reported",
    weight_range = "mean 75 kg (SD 13); range not reported",
    sex_female_pct = 24,
    race_ethnicity = c(White = 57, Black = 33, Asian = 8, Other = 2),
    disease_state = "Healthy volunteers randomised to placebo in single- and multiple-dose phase I trials",
    dose_range = "Placebo only (single or multiple administrations); no active drug",
    regions = "Not reported (AbbVie phase I units)",
    notes = "Minocha 2019 Results 'Data sources and baseline demographics' and Table 1 (heart rate and QT interval column). Heart rate was recorded as the mean of high-precision triplicate ECG measurements per time point, from 3 days before the first placebo dose to about 21 days after it. Studies were conducted 2007-2014."
  )

  ini({
    # Baseline at midnight for a male (Minocha 2019 Table 2 and Eq. 9).
    rbase <- 60.7; label("Typical midnight baseline heart rate for a male (beats/min)") # Table 2, 'Baseline for males' = 60.7 (95% CI 60.1, 61.2)
    e_sexf_rbase <- 6.8; label("Additive shift of the baseline for females (beats/min)") # Table 2, 'Additive factor for female baseline' = 6.8 (5.7, 7.9)

    # Circadian cosines AMP * cos(2 * 3.14 * (t - PS) / period) (Eq. 4).
    amp_24h <- 2.2; label("Amplitude of the 24-h cosine (beats/min)") # Table 2, '24 hours circadian amplitude' = 2.2 (2.0, 2.4)
    acrophase_24h <- 17.6; label("Phase shift (peak clock time) of the 24-h cosine (h after midnight)") # Table 2, '24 hours circadian phase shift' = 17.6 (17.1, 17.9)
    amp_12h <- 0.5; label("Amplitude of the 12-h cosine (beats/min)") # Table 2, '12 hours circadian amplitude' = 0.5 (0.4, 0.7)
    acrophase_12h <- fixed(0); label("Phase shift of the 12-h cosine (h after midnight)") # Table 2 'NS'; Data S1 control stream 'PHI12 = THETA(8)' with '0 FIX'; Results: not estimated because the estimate was close to 12 h (equivalent to 0 for a 12-h period)
    amp_6h <- 2.7; label("Amplitude of the 6-h cosine (beats/min)") # Table 2, '6 hours circadian amplitude' = 2.7 (2.6, 2.9)
    acrophase_6h <- 2.8; label("Phase shift (peak clock time) of the 6-h cosine (h after midnight)") # Table 2, '6 hours circadian phase shift' = 2.8 (2.8, 2.9)

    # Linear time drift TE = SLP * (TIME / 24) (Eq. 5; Data S1 'TE = TSL*(TIME/24)').
    slope_drift <- 0.15; label("Linear drift of heart rate with time (beats/min per day)") # Table 2, 'Time drift model' = 0.15 x (TIME); Results: 0.15 (+/- 0.45) bpm/day

    # Between-subject variability, all additive on the linear scale. The
    # baseline eta is the authors' ETA(1)*THETA(4) with OMEGA(1) held at 1,
    # re-expressed as a single variance = THETA(4)^2 = 7.1^2.
    etarbase ~ 50.41 # Table 2, 'Scaling factor for omega baseline' = 7.1 (6.7, 7.4) with omega^2 held at 1; variance = 7.1^2
    etaamp_24h ~ 3.52 # Corrigendum Table S1, 'Amplitude_24 hrs' variance = 3.52 (2.87, 4.34)
    etaacrophase_24h ~ 3.12 # Corrigendum Table S1, 'Phase shift_ 24 hrs' variance = 3.12 (2.05, 4.17)
    etaamp_12h ~ 1.27 # Corrigendum Table S1, 'Amplitude_12 hrs' variance = 1.27 (0.83, 1.72)
    etaacrophase_6h ~ 0.638 # Corrigendum Table S1, 'Phase shift_ 6 hrs' variance = 0.638 (0.55, 0.75)
    etaslope_drift ~ 0.204 # Corrigendum Table S1, 'Time Drift' variance = 0.204 (0.15, 0.26); Results quote the SD as 0.45 bpm/day

    # Proportional residual error (Table 2 footnote c; Data S1 Y = HR*(1+EPS(2)*THETA(12)), THETA(12) = 1 FIX, additive EPS(1) variance 0 FIX).
    propSd <- 0.0632; label("Proportional residual error (fraction)") # Table 2, 'sigma^2 (additive residual error)' for heart rate = 0.004 (0.004, 0.005), proportional per footnote c; SD = sqrt(0.004) = 0.0632
  })

  model({
    # Individual baseline (Eq. 9).
    rbase_i <- rbase + e_sexf_rbase * SEXF + etarbase

    # Individual circadian parameters (Data S1 control stream).
    amp_24h_i <- amp_24h + etaamp_24h
    acrophase_24h_i <- acrophase_24h + etaacrophase_24h
    amp_12h_i <- amp_12h + etaamp_12h
    acrophase_6h_i <- acrophase_6h + etaacrophase_6h
    slope_drift_i <- slope_drift + etaslope_drift

    # Cosines use 2 * 3.14 exactly as in Eq. 4 and the Data S1 control stream.
    # t is clock time (h) since the midnight before the first measurement.
    circ_24h <- amp_24h_i * cos(2 * 3.14 * (t - acrophase_24h_i) / 24)
    circ_12h <- amp_12h_i * cos(2 * 3.14 * (t - acrophase_12h) / 12)
    circ_6h <- amp_6h * cos(2 * 3.14 * (t - acrophase_6h_i) / 6)

    # Linear drift, slope per day.
    drift <- slope_drift_i * (t / 24)

    hr <- rbase_i + circ_24h + circ_12h + circ_6h + drift
    hr ~ prop(propSd)
  })
}
