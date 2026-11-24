Taubel_2015_E52862_dQTcF <- function() {
  description <- paste(
    "Concentration-DeltaQTcF linear mixed-effects PD model with nominal",
    "time as a factor, for the sigma-1 receptor antagonist E-52862 (S1RA)",
    "in 32 healthy White adults (17 men, 15 women) given single oral doses",
    "of 500, 600 and 800 mg and placebo in a four-period crossover. The",
    "paper fits the change from the averaged day -1 baseline in the",
    "Fridericia-corrected QT interval with time, period, sequence and sex",
    "as fixed factors and a random intercept and slope per subject (best",
    "fit, model 2 of Table 5), and reports only the slope (-0.00007 ms per",
    "ng/mL) and the meal-induced time effects at 6 h and 8 h relative to",
    "pre-dose (-8.1 and -7.2 ms, Table 7, the study's assay-sensitivity",
    "test). This file therefore returns the change in QTcF from the",
    "pre-dose nominal time: QTcF = theta(NTIME) + slope * CP_E52862_NGML,",
    "with theta = 0 at pre-dose, -8.1 ms at 6 h and -7.2 ms at 8 h. The",
    "time effects at all other nominal times are NOT reported, so at those",
    "times the output carries only the concentration term; use the file",
    "for the slope, and for the 6 h and 8 h meal contrasts. PD-only model:",
    "E-52862 plasma concentration is supplied as a time-varying covariate.",
    "Companion files Taubel_2015_E52862_dQTcI.R,",
    "Taubel_2015_E52862_ddQTcI.R and Taubel_2015_E52862_ddQTcF.R.",
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
    concentration = "(observation QTcF is the change in DeltaQTcF from the pre-dose nominal time in the Fridericia-corrected QT interval, ms; driving covariate CP_E52862_NGML is in ng/mL)"
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
        "'Statistical Methods'). The pre-dose concentration is 0, so the",
        "change from pre-dose carries the full concentration term.",
        "Observed exposure (Table 3, mean (SD)): Cmax 3938.5 (791.6),",
        "4357.9 (709.5) and 5580.1 (1177.2) ng/mL after 500, 600 and",
        "800 mg; geometric-mean Cmax used for the Table 6 predictions",
        "3859.6, 4302.4 and 5464.0 ng/mL."
      ),
      source_name = "plasma concentration of unchanged E-52862"
    ),
    NTIME = list(
      description = "Nominal post-dose ECG time point in hours. Only 0 (pre-dose; reference), 6 and 8 h carry a reported time effect.",
      units = "h",
      type = "categorical",
      reference_category = "0 (pre-dose)",
      notes = paste(
        "Taubel 2015 Table 7 prints the QTcF time effects only as the",
        "change from pre-dose at 6 h (1 h after the 5 h lunch) and 8 h",
        "(3 h after lunch). model() assigns 5 <= NTIME < 7 to the 6 h",
        "level and 7 <= NTIME < 10 to the 8 h level; every other NTIME",
        "gets a zero time effect because the paper does not report it.",
        "Supply the nominal time, not rxode2's tad()."
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
      "post-dose on day -1 and day 1. QTcF = QT / RR^0.33.",
      "Fitted with nlme in R 2.13. The paper also fitted the same model",
      "to each of five metabolites (M1-M5) but reports slopes for the",
      "parent only."
    )
  )

  ini({
    # Taubel 2015 Table 5, row 'QTcF, model 2' (best fit by AIC, 11266.3):
    # slope -0.00007 ms per ng/mL, 95% CI -0.00049 to 0.00035. Table 6
    # rows 'Based on change from average baseline, QTcF' equal
    # slope * Cmax (-0.4, -0.3, -0.3 ms at 5464.0, 4302.4, 3859.6 ng/mL).
    slope <- -0.00007
    label("Linear concentration-DeltaQTcF slope (ms per ng/mL)")

    # Taubel 2015 Table 7 'QTcF prediction of change from pre-dose':
    # 6 h -8.1 ms (95% CI -10.4 to -5.9); 8 h -7.2 ms (95% CI -9.4 to -5.0).
    e_ntime6_e0 <- -8.1
    label("Time effect at 6 h relative to pre-dose, 1 h after lunch (ms)")
    e_ntime8_e0 <- -7.2
    label("Time effect at 8 h relative to pre-dose, 3 h after lunch (ms)")

    # Random intercept and random slope per subject (model 2) are not
    # reported numerically, so no etas are declared (typical-value model).
    # Residual error is not reported either.
    addSd <- fixed(0)
    label("Additive residual error SD on the change in DeltaQTcF (ms; not reported in source)")
  })

  model({
    nt6 <- (NTIME >= 5) * (NTIME < 7)
    nt8 <- (NTIME >= 7) * (NTIME < 10)
    ntime_effect <- e_ntime6_e0 * nt6 + e_ntime8_e0 * nt8

    # Change from the pre-dose nominal time: the pre-dose time effect is
    # the reference and the pre-dose concentration is 0, so both cancel.
    QTcF <- ntime_effect + slope * CP_E52862_NGML
    QTcF ~ add(addSd)
  })
}
