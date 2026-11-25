Minocha_2019_placebo_QTcF <- function() {
  description <- "Placebo (no-drug) Fridericia-corrected QT interval (QTcF) model for healthy adult volunteers in AbbVie phase I trials (Minocha 2019): subject baseline at midnight with a power effect of age and additive female-sex and Black-race shifts, plus three additive cosine circadian components (24-h, 12-h and 6-h periods) and an exponential time drift that lowers QTcF towards a plateau, with additive between-subject variability on the baseline and on the three acrophases, and an additive residual error. Time is clock time in hours since the midnight before the subject's first recorded measurement. Intended as a placebo backbone for concentration-QTc analyses."
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
    concentration = "ms (QTcF interval)"
  )

  covariateData <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power model on the typical baseline, (AGE / 30)^0.031 (Minocha 2019 Eq. 10 and Table 2; reference age 30 years per Methods, Eq. 3).",
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Additive shift on the midnight baseline QTcF (Eq. 10 and Table 2); males are the reference (most frequent) category.",
      source_name = "SEX"
    ),
    RACE_BLACK = list(
      description = "Black race indicator (1 = Black, 0 = other)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White, Asian or other; White is the most frequent category)",
      notes = "Additive shift on the midnight baseline QTcF (Eq. 10 and Table 2). The Asian and other-race shifts were not significant, so those subjects share the White reference baseline.",
      source_name = "RACE"
    )
  )

  covariatesDataExcluded <- list(
    RACE_ASIAN = list(
      description = "Asian race indicator (1 = Asian, 0 = other)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White, Black or other)",
      notes = "Screened as an additive baseline shift (Methods, Eq. 3) and not significant for QTcF (Table 2, 'NS').",
      source_name = "RACE"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened as a power model normalised to 70 kg and not significant for QTcF (Table 2, 'NS'; Results 'QT interval').",
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
    notes = "Minocha 2019 Results 'Data sources and baseline demographics' and Table 1 (heart rate and QT interval column). QT was corrected for heart rate with Fridericia's formula, QTcF = QT / (RR / 1000)^(1/3) (Eq. 1), each observation the mean of high-precision triplicate ECG measurements, from 3 days before the first placebo dose to about 21 days after it. Studies were conducted 2007-2014."
  )

  ini({
    # Baseline at midnight for a 30-year-old White male (Eq. 10, Table 2).
    rbase <- 400; label("Typical midnight baseline QTcF for a 30-year-old White male (ms)") # Table 2, 'Baseline for males' = 400 (95% CI 398, 402)
    e_age_rbase <- 0.031; label("Power exponent of age (normalised to 30 years) on the typical baseline (unitless)") # Table 2, 'Exponent for age; power model (normalized to 30 years)' = 0.031 (0.025, 0.039)
    e_sexf_rbase <- 14.7; label("Additive shift of the baseline for females (ms)") # Table 2, 'Additive factor for female baseline' = 14.7 (12.8, 16.5)
    e_race_black_rbase <- -3.4; label("Additive shift of the baseline for Black subjects (ms)") # Table 2, 'Additive factor for black race' = -3.4 (-5.4, -1.5)

    # Circadian cosines AMP * cos(2 * 3.14 * (t - PS) / period) (Eq. 4).
    amp_24h <- 2.6; label("Amplitude of the 24-h cosine (ms)") # Table 2, '24 hours circadian amplitude' = 2.6 (2.2, 2.9)
    acrophase_24h <- 2.2; label("Phase shift (peak clock time) of the 24-h cosine (h after midnight)") # Table 2, '24 hours circadian phase shift' = 2.2 (1.9, 2.6)
    amp_12h <- 2.0; label("Amplitude of the 12-h cosine (ms)") # Table 2, '12 hours circadian amplitude' = 2.0 (1.8, 2.3)
    acrophase_12h <- 3.6; label("Phase shift (peak clock time) of the 12-h cosine (h after midnight)") # Table 2, '12 hours circadian phase shift' = 3.6 (3.4, 3.8)
    amp_6h <- 1.2; label("Amplitude of the 6-h cosine (ms)") # Table 2, '6 hours circadian amplitude' = 1.2 (1.1, 1.4)
    acrophase_6h <- 1.7; label("Phase shift (peak clock time) of the 6-h cosine (h after midnight)") # Table 2, '6 hours circadian phase shift' = 1.7 (1.5, 1.8)

    # Exponential drift TE = MAX * (1 - exp(-(TIME / 24) * SLP)) (Eq. 6), subtracted (Eq. 10).
    drift_max <- 5.0; label("Maximum reduction of QTcF by the time drift (ms)") # Table 2, 'Time drift model' = 5.0 x (1-exp(-1.4 x TIME)); Results: 'maximum reduction of 5 msec'
    lkdrift <- log(1.4); label("First-order rate constant of the time drift (1/day)") # Table 2, 'Time drift model' rate 1.4; Results: '1.4 day-1 (half-life of ~ 12 hours)'

    # Between-subject variability, additive on the linear scale.
    etarbase ~ 193.21 # Table 2, 'Scaling factor for omega baseline' = 13.9 (13.1, 14.7) with omega^2 held at 1; variance = 13.9^2
    etaacrophase_24h ~ 8.09 # Corrigendum Table S1, 'Phase shift_ 24 hrs' variance = 8.09 (5.65, 11.79)
    etaacrophase_12h ~ 2.24 # Corrigendum Table S1, 'Phase shift_ 12 hrs' variance = 2.24 (1.73, 2.76)
    etaacrophase_6h ~ 0.967 # Corrigendum Table S1, 'Phase shift_ 6 hrs' variance = 0.967 (0.65, 1.24)

    # Additive residual error, epsilon * THETA with sigma^2 = 1 (Eq. 8).
    addSd <- 5.7; label("Additive residual error (ms)") # Table 2, 'Scaling factor for sigma' = 5.7 (5.5, 5.8) with sigma^2 = 1
  })

  model({
    # Individual baseline (Eq. 10): the power term scales only the typical
    # value; the categorical shifts and eta are additive (Eq. 3).
    rbase_i <- rbase * (AGE / 30)^e_age_rbase + e_sexf_rbase * SEXF +
      e_race_black_rbase * RACE_BLACK + etarbase

    acrophase_24h_i <- acrophase_24h + etaacrophase_24h
    acrophase_12h_i <- acrophase_12h + etaacrophase_12h
    acrophase_6h_i <- acrophase_6h + etaacrophase_6h

    # Cosines use 2 * 3.14 exactly as in Eq. 4 and the Data S1 control stream.
    # t is clock time (h) since the midnight before the first measurement.
    circ_24h <- amp_24h * cos(2 * 3.14 * (t - acrophase_24h_i) / 24)
    circ_12h <- amp_12h * cos(2 * 3.14 * (t - acrophase_12h_i) / 12)
    circ_6h <- amp_6h * cos(2 * 3.14 * (t - acrophase_6h_i) / 6)

    kdrift <- exp(lkdrift)
    drift <- drift_max * (1 - exp(-kdrift * t / 24))

    QTcF <- rbase_i + circ_24h + circ_12h + circ_6h - drift
    QTcF ~ add(addSd)
  })
}
