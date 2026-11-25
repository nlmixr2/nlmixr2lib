Minocha_2019_placebo_dbp <- function() {
  description <- "Placebo (no-drug) diastolic blood pressure model for healthy adult volunteers in AbbVie phase I trials (Minocha 2019): subject baseline at midnight with a power effect of age and additive female-sex, Black-race and Asian-race shifts, plus a single 24-h cosine circadian component and an exponential time drift that lowers diastolic pressure towards a plateau, with additive between-subject variability on the baseline, the 24-h amplitude and the 24-h acrophase, and an additive residual error. Time is clock time in hours since the midnight before the subject's first recorded measurement. Intended as a placebo backbone for exposure-response analyses of blood pressure."
  reference <- paste(
    "Minocha M, Li H, Chiu YL, Carter D, Othman AA.",
    "Models of Variability and Circadian Rhythm in Heart Rate, Blood Pressure,",
    "and QT Interval for Healthy Subjects Who Received Placebo in Phase I Trials.",
    "Clin Transl Sci. 2019;12(5):470-480. doi:10.1111/cts.12640.",
    "Corrigendum: Clin Transl Sci. 2020;13(3):638. doi:10.1111/cts.12773",
    "(Table S1, between-subject variability estimates missing from Table 2;",
    "revised Figure 2).",
    sep = " "
  )
  vignette <- "Minocha_2019_placebo_cardiovascular"
  units <- list(
    time = "h",
    dosing = "n/a (placebo baseline model with no drug input)",
    concentration = "mmHg (diastolic blood pressure)"
  )

  covariateData <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power model on the typical baseline, (AGE / 30)^0.08 (Minocha 2019 Eq. 12 and Table 2; reference age 30 years per Methods, Eq. 3).",
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Additive shift on the midnight baseline diastolic pressure (Eq. 12 and Table 2); males are the reference (most frequent) category.",
      source_name = "SEX"
    ),
    RACE_BLACK = list(
      description = "Black race indicator (1 = Black, 0 = other)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White, Asian or other; White is the most frequent category)",
      notes = "Additive shift on the midnight baseline diastolic pressure (Eq. 12 and Table 2). The other-race shift was not significant, so those subjects share the White reference baseline.",
      source_name = "RACE"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator (1 = Asian, 0 = other)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White, Black or other; White is the most frequent category)",
      notes = "Additive shift on the midnight baseline diastolic pressure (Eq. 12 and Table 2).",
      source_name = "RACE"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened as a power model normalised to 70 kg and not significant for diastolic pressure (Table 2, 'NS'; Results 'Diastolic BP').",
      source_name = "WTKG"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1035L,
    n_studies = 65L,
    n_observations = 18807L,
    age_range = "adults; mean 34 years (SD 10); range not reported",
    weight_range = "mean 75 kg (SD 12); range not reported",
    sex_female_pct = 25,
    race_ethnicity = c(White = 57, Black = 33, Asian = 8, Other = 2),
    disease_state = "Healthy volunteers randomised to placebo in single- and multiple-dose phase I trials",
    dose_range = "Placebo only (single or multiple administrations); no active drug",
    regions = "Not reported (AbbVie phase I units)",
    notes = "Minocha 2019 Results 'Data sources and baseline demographics' and Table 1 (systolic and diastolic blood pressure column). Blood pressure was recorded from 3 days before the first placebo dose to 12 weeks after it. Studies were conducted 2007-2014."
  )

  ini({
    # Baseline at midnight for a 30-year-old White male (Eq. 12, Table 2).
    rbase <- 78.3; label("Typical midnight baseline diastolic pressure for a 30-year-old White male (mmHg)") # Table 2, 'Baseline for males' = 78.3 (95% CI 76.9, 80.2)
    e_age_rbase <- 0.08; label("Power exponent of age (normalised to 30 years) on the typical baseline (unitless)") # Table 2, 'Exponent for age; power model (normalized to 30 years)' = 0.08 (0.07, 1.0); the printed upper bound 1.0 is presumably 0.10
    e_sexf_rbase <- -4.7; label("Additive shift of the baseline for females (mmHg)") # Table 2, 'Additive factor for female baseline' = -4.7 (-5.7, -3.9)
    e_race_black_rbase <- 2.1; label("Additive shift of the baseline for Black subjects (mmHg)") # Table 2, 'Additive factor for black race' = 2.1 (1.5, 2.9)
    e_race_asian_rbase <- -1.5; label("Additive shift of the baseline for Asian subjects (mmHg)") # Table 2, 'Additive factor for Asian race' = -1.5 (-1.0, -2.1)

    # Circadian cosine AMP * cos(2 * 3.14 * (t - PS) / 24) (Eq. 4).
    amp_24h <- 1.6; label("Amplitude of the 24-h cosine (mmHg)") # Table 2, '24 hours circadian amplitude' = 1.6 (1.4, 1.9)
    acrophase_24h <- 2.3; label("Phase shift (peak clock time) of the 24-h cosine (h after midnight)") # Table 2, '24 hours circadian phase shift' = 2.3 (1.9, 2.8)

    # Exponential drift TE = MAX * (1 - exp(-(TIME / 24) * SLP)) (Eq. 6), subtracted (Eq. 12).
    drift_max <- 4.8; label("Maximum reduction of diastolic pressure by the time drift (mmHg)") # Table 2, 'Time drift model' = 4.8 x (1-exp(-1.6 x TIME)); Results: 'maximum reduction of 4.8 mmHg'
    lkdrift <- log(1.6); label("First-order rate constant of the time drift (1/day)") # Table 2, 'Time drift model' rate 1.6; Results: '1.6 day-1 (half-life of ~ 10 hours)'

    # Between-subject variability, additive on the linear scale.
    etarbase ~ 36 # Table 2, 'Scaling factor for omega baseline' = 6.0 (5.6, 6.3) with omega^2 held at 1; variance = 6.0^2
    etaamp_24h ~ 0.93 # Corrigendum Table S1, 'Amplitude_24 hrs' variance = 0.93 (0.30, 2.37)
    etaacrophase_24h ~ 6.12 # Corrigendum Table S1, 'Phase shift_ 24 hrs' variance = 6.12 (2.57, 8.92)

    # Additive residual error, epsilon * THETA with sigma^2 = 1 (Eq. 8).
    addSd <- 5.2; label("Additive residual error (mmHg)") # Table 2, 'Scaling factor for sigma' = 5.2 (5.1, 5.3) with sigma^2 = 1
  })

  model({
    # Individual baseline (Eq. 12): the power term scales only the typical
    # value; the categorical shifts and eta are additive (Eq. 3).
    rbase_i <- rbase * (AGE / 30)^e_age_rbase + e_sexf_rbase * SEXF +
      e_race_black_rbase * RACE_BLACK + e_race_asian_rbase * RACE_ASIAN +
      etarbase

    amp_24h_i <- amp_24h + etaamp_24h
    acrophase_24h_i <- acrophase_24h + etaacrophase_24h

    # Cosine uses 2 * 3.14 exactly as in Eq. 4 and the Data S1 control stream.
    # t is clock time (h) since the midnight before the first measurement.
    circ_24h <- amp_24h_i * cos(2 * 3.14 * (t - acrophase_24h_i) / 24)

    kdrift <- exp(lkdrift)
    drift <- drift_max * (1 - exp(-kdrift * t / 24))

    dbp <- rbase_i + circ_24h - drift
    dbp ~ add(addSd)
  })
}
