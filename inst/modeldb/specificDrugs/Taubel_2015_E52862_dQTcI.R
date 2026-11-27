Taubel_2015_E52862_dQTcI <- function() {
  description <- paste(
    "Concentration-DeltaQTcI linear mixed-effects PD model with nominal",
    "time as a factor, for the sigma-1 receptor antagonist E-52862 (S1RA)",
    "in 32 healthy White adults (17 men, 15 women) given single oral doses",
    "of 500, 600 and 800 mg and placebo in a four-period crossover. The",
    "endpoint is the change from the averaged day -1 baseline in the",
    "individually rate-corrected QT interval (DeltaQTcI, ms), with placebo",
    "periods entering at concentration 0. Typical-value structure:",
    "DeltaQTcI = e0 + theta(NTIME) + slope * CP_E52862_NGML, with slope =",
    "-0.00079 ms per ng/mL (best-fitting model 2 of the paper's Table 5:",
    "time, period, sequence and sex as fixed factors, random intercept and",
    "random slope per subject). The drug-free time course e0 + theta(NTIME)",
    "over 0-12 h, which carries the QTc shortening after the 5 h lunch, is",
    "printed only graphically (Figure 5) and was digitised as the mean of",
    "the six analyte-specific curves (parent and metabolites M1-M5, which",
    "the figure does not label); it is accurate to about +/-0.8 ms. PD-only",
    "model: E-52862 plasma concentration is supplied as a time-varying",
    "covariate; the paper characterises the PK by NCA only. Companion files",
    "Taubel_2015_E52862_dQTcF.R, Taubel_2015_E52862_ddQTcI.R and",
    "Taubel_2015_E52862_ddQTcF.R.",
    sep = " "
  )

  reference <- paste(
    "Taubel J, Ferber G, Lorch U, Wang D, Sust M, Camm AJ. (2015).",
    "Single Doses up to 800 mg of E-52862 Do Not Prolong the QTc Interval",
    "- A Retrospective Validation by Pharmacokinetic-Pharmacodynamic",
    "Modelling of Electrocardiography Data Utilising the Effects of a Meal",
    "on QTc to Demonstrate ECG Assay Sensitivity.",
    "PLoS One 10(8):e0136369. doi:10.1371/journal.pone.0136369.",
    sep = " "
  )

  vignette <- "Taubel_2015_E52862"

  units <- list(
    time = "h",
    dosing = "(none; PD-only model fed by an external E-52862 plasma-concentration covariate)",
    concentration = "(observation QTcI is DeltaQTcI, the change from averaged day -1 baseline in the individually corrected QT interval, ms; driving covariate CP_E52862_NGML is in ng/mL)"
  )

  covariateData <- list(
    CP_E52862_NGML = list(
      description = "Instantaneous E-52862 (parent) plasma concentration at the time of each ECG, supplied as a time-varying covariate from observed plasma samples or an external PK source.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives the slope * CP_E52862_NGML term. Taubel 2015 regressed on",
        "the observed LC-MS/MS plasma concentration drawn at the ECG time",
        "points; placebo data entered with concentration 0 (Methods",
        "'Statistical Methods'). Observed exposure (Table 3, mean (SD)):",
        "Cmax 3938.5 (791.6), 4357.9 (709.5) and 5580.1 (1177.2) ng/mL",
        "after 500, 600 and 800 mg; geometric-mean Cmax used for the",
        "Table 6 predictions 3859.6, 4302.4 and 5464.0 ng/mL."
      ),
      source_name = "plasma concentration of unchanged E-52862"
    ),
    NTIME = list(
      description = "Nominal post-dose ECG time point in hours. Levels 0 (pre-dose; reference), 0.5, 1, 1.5, 2, 3, 4, 6, 8 and 12 h; each non-reference level carries its own intercept shift.",
      units = "h",
      type = "categorical",
      reference_category = "0 (pre-dose)",
      notes = paste(
        "Taubel 2015 Methods: 'time can and should be introduced as factor",
        "in the model'. Figure 5 plots the estimated time course at the",
        "0, 0.5, 1, 1.5, 2, 3, 4, 6, 8 and 12 h levels; the 0.25, 0.75,",
        "1.25, 1.75 and 5 h ECG slots and the 24-96 h slots are not shown,",
        "so their time effects are unknown. model() bins the supplied value",
        "to the nearest plotted level (cut points 0.25, 0.75, 1.25, 1.75,",
        "2.5, 3.5, 5, 7 and 10 h); any NTIME >= 10 h takes the 12 h level,",
        "so the model is valid over 0-12 h only. Supply the nominal time,",
        "not rxode2's tad()."
      ),
      source_name = "time (factor)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 32L,
    n_studies = 1L,
    age_range = "18-35 years (inclusion); mean 24.22 (SD 4.668)",
    weight_range = "mean 67.39 kg (SD 9.654); BMI 18-25 kg/m^2 (inclusion), mean 22.41",
    sex_female_pct = round(100 * 15 / 32, 1),
    race_ethnicity = c(White = 100),
    disease_state = "Healthy adult volunteers (all White, non-Hispanic).",
    dose_range = paste(
      "Single oral doses of E-52862 500, 600 and 800 mg and matching",
      "placebo, randomised four-period crossover with escalating active",
      "doses and >= 7-day washouts; each period had a placebo baseline",
      "ECG day (day -1) time-matched to the dosing day (day 1)."
    ),
    regions = "United Kingdom (Richmond Pharmacology, St George's University of London)",
    notes = paste(
      "EudraCT 2010-020343-13. 31 of 32 subjects completed all four",
      "periods. Triplicate 12-lead ECGs and PK samples at pre-dose and",
      "0.25-96 h post-dose; standardised lunch at 5 h and dinner at 9 h",
      "post-dose on day -1 and day 1. QTcI from 24-h screening Holter",
      "(linear or log-log-linear per subject). Fitted with nlme in R 2.13.",
      "The paper also fitted the same model to each of five metabolites",
      "(M1-M5) but reports slopes for the parent only."
    )
  )

  ini({
    # Taubel 2015 Table 5, row 'QTcI, model 2' (best fit by AIC, 11535.7):
    # slope -0.00079 ms per ng/mL, 95% CI -0.00143 to -0.00016. Model 2 =
    # 'time, period, sequence and sex as factors, random intercept' plus a
    # random slope (Table 5 footnote 1). Table 6 rows 'Based on change from
    # average baseline, QTcI' equal slope * Cmax (-4.3, -3.4, -3.1 ms at
    # 5464.0, 4302.4, 3859.6 ng/mL).
    slope <- -0.00079
    label("Linear concentration-DeltaQTcI slope (ms per ng/mL)")

    # Drug-free time course, digitised from Taubel 2015 Figure 5 (axis
    # calibrated on the 5, 0, -5 and -10 ms gridlines and the 0-12 h
    # ticks). The figure overlays the six analyte-specific models without
    # a legend, so each value is the mean of the six curves at that time;
    # the half-spread of the six curves is about 0.6-0.9 ms. Increments are
    # relative to the pre-dose level. Digitised absolute levels (ms):
    # 0 h 0.9, 0.5 h -0.8, 1 h 0.9, 1.5 h 1.2, 2 h 1.2, 3 h -0.5,
    # 4 h 0.0, 6 h -6.9, 8 h -6.1, 12 h -3.7.
    e0 <- 0.9
    label("DeltaQTcI at the pre-dose nominal time, drug-free (ms)")
    e_ntime0p5_e0 <- -1.7
    label("Time effect at 0.5 h relative to pre-dose (ms)")
    e_ntime1_e0 <- 0.0
    label("Time effect at 1 h relative to pre-dose (ms)")
    e_ntime1p5_e0 <- 0.3
    label("Time effect at 1.5 h relative to pre-dose (ms)")
    e_ntime2_e0 <- 0.3
    label("Time effect at 2 h relative to pre-dose (ms)")
    e_ntime3_e0 <- -1.4
    label("Time effect at 3 h relative to pre-dose (ms)")
    e_ntime4_e0 <- -0.9
    label("Time effect at 4 h relative to pre-dose (ms)")
    e_ntime6_e0 <- -7.8
    label("Time effect at 6 h relative to pre-dose, 1 h after lunch (ms)")
    e_ntime8_e0 <- -7.0
    label("Time effect at 8 h relative to pre-dose, 3 h after lunch (ms)")
    e_ntime12_e0 <- -4.6
    label("Time effect at 12 h relative to pre-dose (ms)")

    # Random intercept and random slope per subject (model 2) are not
    # reported numerically, so no etas are declared (typical-value model).
    # Residual error is not reported either.
    addSd <- fixed(0)
    label("Additive residual error SD on DeltaQTcI (ms; not reported in source)")
  })

  model({
    # Bin NTIME to the nearest level plotted in Taubel 2015 Figure 5;
    # NTIME < 0.25 is the pre-dose reference (all indicators 0).
    nt0p5 <- (NTIME >= 0.25) * (NTIME < 0.75)
    nt1 <- (NTIME >= 0.75) * (NTIME < 1.25)
    nt1p5 <- (NTIME >= 1.25) * (NTIME < 1.75)
    nt2 <- (NTIME >= 1.75) * (NTIME < 2.5)
    nt3 <- (NTIME >= 2.5) * (NTIME < 3.5)
    nt4 <- (NTIME >= 3.5) * (NTIME < 5)
    nt6 <- (NTIME >= 5) * (NTIME < 7)
    nt8 <- (NTIME >= 7) * (NTIME < 10)
    nt12 <- (NTIME >= 10)
    ntime_effect <-
      e_ntime0p5_e0 * nt0p5 + e_ntime1_e0 * nt1 + e_ntime1p5_e0 * nt1p5 +
      e_ntime2_e0 * nt2 + e_ntime3_e0 * nt3 + e_ntime4_e0 * nt4 +
      e_ntime6_e0 * nt6 + e_ntime8_e0 * nt8 + e_ntime12_e0 * nt12

    # Taubel 2015 Methods 'Statistical Methods', second series of models:
    # change from averaged baseline with time as a factor and placebo at
    # concentration 0, so the slope is corrected for the spontaneous
    # (meal-driven) time course.
    QTcI <- e0 + ntime_effect + slope * CP_E52862_NGML
    QTcI ~ add(addSd)
  })
}
