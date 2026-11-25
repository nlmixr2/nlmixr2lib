Minocha_2019_placebo_sbp <- function() {
  description <- "Placebo (no-drug) systolic blood pressure model for healthy adult volunteers in AbbVie phase I trials (Minocha 2019): subject baseline at midnight with a power effect of body weight and additive female-sex and Black-race shifts, plus a single 24-h cosine circadian component and an exponential time drift that lowers systolic pressure towards a plateau, with additive between-subject variability on the baseline and the 24-h amplitude, and an additive residual error. Time is clock time in hours since the midnight before the subject's first recorded measurement. Intended as a placebo backbone for exposure-response analyses of blood pressure."
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
    concentration = "mmHg (systolic blood pressure)"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power model on the typical baseline, (WT / 70)^0.09 (Minocha 2019 Eq. 11 and Table 2; reference weight 70 kg per Methods, Eq. 3).",
      source_name = "WTKG"
    ),
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Additive shift on the midnight baseline systolic pressure (Eq. 11 and Table 2); males are the reference (most frequent) category.",
      source_name = "SEX"
    ),
    RACE_BLACK = list(
      description = "Black race indicator (1 = Black, 0 = other)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White, Asian or other; White is the most frequent category)",
      notes = "Additive shift on the midnight baseline systolic pressure (Eq. 11 and Table 2). The Asian and other-race shifts were not significant, so those subjects share the White reference baseline.",
      source_name = "RACE"
    )
  )

  covariatesDataExcluded <- list(
    RACE_ASIAN = list(
      description = "Asian race indicator (1 = Asian, 0 = other)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White, Black or other)",
      notes = "Screened as an additive baseline shift (Methods, Eq. 3) and not significant for systolic pressure (Table 2, 'NS').",
      source_name = "RACE"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened as a power model normalised to 30 years and not significant for systolic pressure (Table 2, 'NS'; Results 'Systolic BP').",
      source_name = "AGE"
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
    # Baseline at midnight for a 70-kg White male (Eq. 11, Table 2).
    rbase <- 123; label("Typical midnight baseline systolic pressure for a 70-kg White male (mmHg)") # Table 2, 'Baseline for males' = 123 (95% CI 122, 125)
    e_wt_rbase <- 0.09; label("Power exponent of body weight (normalised to 70 kg) on the typical baseline (unitless)") # Table 2, 'Exponent for weight; power model (normalized to 70 kg)' = 0.09 (0.06, 0.12)
    e_sexf_rbase <- -8.3; label("Additive shift of the baseline for females (mmHg)") # Table 2, 'Additive factor for female baseline' = -8.3 (-9.6, -7.2)
    e_race_black_rbase <- 3.8; label("Additive shift of the baseline for Black subjects (mmHg)") # Table 2, 'Additive factor for black race' = 3.8 (2.8, 4.9)

    # Circadian cosine AMP * cos(2 * 3.14 * (t - PS) / 24) (Eq. 4).
    amp_24h <- 1.0; label("Amplitude of the 24-h cosine (mmHg)") # Table 2, '24 hours circadian amplitude' = 1.0 (0.7, 1.3)
    acrophase_24h <- 19.3; label("Phase shift (peak clock time) of the 24-h cosine (h after midnight)") # Table 2, '24 hours circadian phase shift' = 19.3 (18.4, 20.3)

    # Exponential drift TE = MAX * (1 - exp(-(TIME / 24) * SLP)) (Eq. 6), subtracted (Eq. 11).
    drift_max <- 7.6; label("Maximum reduction of systolic pressure by the time drift (mmHg)") # Table 2, 'Time drift model' = 7.6 x (1-exp(-1.2 x TIME)); Results: 'maximum reduction of 7.6 mmHg'
    lkdrift <- log(1.2); label("First-order rate constant of the time drift (1/day)") # Table 2, 'Time drift model' rate 1.2; Results: '1.2 day-1 (half-life of ~ 14 hours)'

    # Between-subject variability, additive on the linear scale.
    etarbase ~ 54.76 # Table 2, 'Scaling factor for omega baseline' = 7.4 (6.9, 7.8) with omega^2 held at 1; variance = 7.4^2
    etaamp_24h ~ 3.44 # Corrigendum Table S1, 'Amplitude_24 hrs' variance = 3.44 (2.3, 4.7)

    # Additive residual error, epsilon * THETA with sigma^2 = 1 (Eq. 8).
    addSd <- 7.5; label("Additive residual error (mmHg)") # Table 2, 'Scaling factor for sigma' = 7.5 (7.4, 7.7) with sigma^2 = 1
  })

  model({
    # Individual baseline (Eq. 11): the power term scales only the typical
    # value; the categorical shifts and eta are additive (Eq. 3).
    rbase_i <- rbase * (WT / 70)^e_wt_rbase + e_sexf_rbase * SEXF +
      e_race_black_rbase * RACE_BLACK + etarbase

    amp_24h_i <- amp_24h + etaamp_24h

    # Cosine uses 2 * 3.14 exactly as in Eq. 4 and the Data S1 control stream.
    # t is clock time (h) since the midnight before the first measurement.
    circ_24h <- amp_24h_i * cos(2 * 3.14 * (t - acrophase_24h) / 24)

    kdrift <- exp(lkdrift)
    drift <- drift_max * (1 - exp(-kdrift * t / 24))

    sbp <- rbase_i + circ_24h - drift
    sbp ~ add(addSd)
  })
}
