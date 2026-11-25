Taubel_2015_E52862_ddQTcF <- function() {
  description <- paste(
    "Concentration-DeltaDeltaQTcF linear mixed-effects PD model for the",
    "sigma-1 receptor antagonist E-52862 (S1RA) in 32 healthy White adults",
    "(17 men, 15 women) given single oral doses of 500, 600 and 800 mg and",
    "placebo in a four-period crossover. The endpoint is the difference to",
    "time-matched placebo of the change from the averaged day -1 baseline",
    "in the Fridericia-corrected QT interval (DeltaDeltaQTcF, ms).",
    "Typical-value structure: DeltaDeltaQTcF = e0 + slope *",
    "CP_E52862_NGML, with slope = 0.00051 ms per ng/mL (best-fitting",
    "model 2 of the paper's Table 4: period, sequence and sex as fixed",
    "factors, random intercept and random slope per subject) and e0 =",
    "-1.58 ms, the population-averaged intercept (fixed intercept plus the",
    "period, sequence and sex effects, none printed individually)",
    "back-solved from the model predictions of Table 6. The slope is not",
    "significantly different from zero (90% CI -0.00011 to 0.00113). PD-only",
    "model: E-52862 plasma concentration is supplied as a time-varying",
    "covariate; the paper characterises the PK by NCA only. Companion files",
    "Taubel_2015_E52862_ddQTcI.R (individual correction),",
    "Taubel_2015_E52862_dQTcI.R and Taubel_2015_E52862_dQTcF.R (change",
    "from baseline with time as a factor).",
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
    concentration = "(observation QTcF is DeltaDeltaQTcF, the placebo-corrected change from averaged day -1 baseline in the Fridericia-corrected QT interval, ms; driving covariate CP_E52862_NGML is in ng/mL)"
  )

  covariateData <- list(
    CP_E52862_NGML = list(
      description = "Instantaneous E-52862 (parent) plasma concentration at the time of each ECG, supplied as a time-varying covariate from observed plasma samples or an external PK source.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives DeltaDeltaQTcF = e0 + slope * CP_E52862_NGML. Taubel 2015",
        "regressed on the observed LC-MS/MS plasma concentration drawn at",
        "the ECG time points (Methods 'Pharmacokinetic assessments').",
        "Observed exposure (Table 3, mean (SD)): Cmax 3938.5 (791.6),",
        "4357.9 (709.5) and 5580.1 (1177.2) ng/mL after 500, 600 and",
        "800 mg; geometric-mean Cmax used for the Table 6 predictions",
        "3859.6, 4302.4 and 5464.0 ng/mL. Set to 0 for placebo."
      ),
      source_name = "plasma concentration of unchanged E-52862"
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
      "to each of five metabolites (M1-M5) but reports results for the",
      "parent only."
    )
  )

  ini({
    # Taubel 2015 Table 4, row 'QTcF, model 2' (best fit by AIC, 8768.2):
    # slope 0.00051 ms per ng/mL, 90% CI -0.00011 to 0.00113. Model 2 =
    # 'period, sequence and sex as factors, random intercept' plus a
    # random slope (Table 4 footnote a).
    slope <- 0.00051
    label("Linear concentration-DeltaDeltaQTcF slope (ms per ng/mL)")

    # Not printed. Back-solved from Taubel 2015 Table 6 rows 'Based on
    # difference to placebo of change from average baseline, QTcF'
    # (predictions 1.2, 0.6 and 0.4 ms at 5464.0, 4302.4 and
    # 3859.6 ng/mL). With the printed slope every e0 in [-1.618, -1.544]
    # reproduces all three printed predictions to their 0.1 ms rounding;
    # -1.58 is the interval centre. It is the population-averaged sum of
    # the fixed intercept and the period, sequence and sex effects, which
    # the paper does not print individually.
    e0 <- -1.58
    label("Population-averaged intercept on DeltaDeltaQTcF (ms)")

    # Random intercept and random slope per subject (model 2) are not
    # reported numerically, so no etas are declared (typical-value model).
    # Residual error is not reported either.
    addSd <- fixed(0)
    label("Additive residual error SD on DeltaDeltaQTcF (ms; not reported in source)")
  })

  model({
    # Taubel 2015 Methods 'Statistical Methods', first series of models:
    # placebo-corrected change from averaged baseline regressed on the
    # concentration of one analyte, with period, sequence and sex as fixed
    # effects. No ODE states: purely algebraic in the concentration column.
    QTcF <- e0 + slope * CP_E52862_NGML
    QTcF ~ add(addSd)
  })
}
