Abbas_2020_sertraline_dQTcF <- function() {
  description <- paste(
    "Concentration-DeltaQTcF linear mixed-effects PD model for the selective",
    "serotonin reuptake inhibitor sertraline in 54 healthy adults from a",
    "randomized, double-blind, placebo- and moxifloxacin-controlled, 3-period",
    "crossover thorough-QT study (sertraline titrated to a supratherapeutic",
    "400 mg/day, 200 mg twice daily, at steady state). The endpoint is the",
    "change from baseline in the Fridericia-corrected QT interval",
    "(DeltaQTcF, ms), not placebo-corrected. Typical-value structure:",
    "DeltaQTcF = e0 + slope * CP_SERTRALINE_NGML, the population regression",
    "line the paper draws in Figure 3A, with the printed slope 0.036 ms per",
    "ng/mL and e0 = -5.9 ms digitised from the line's value at zero",
    "concentration. e0 is a population average over the treatment, nominal",
    "time and baseline terms of the paper's model, none of which is printed",
    "individually. Random effects and residual error are not reported, so the",
    "model is typical-value only. PD-only model: sertraline plasma",
    "concentration is supplied as a time-varying covariate. Companion files",
    "Abbas_2020_sertraline_ddQTcF.R (placebo-corrected) and",
    "Abbas_2020_desmethylsertraline_dQTcF.R (N-desmethylsertraline).",
    sep = " "
  )

  reference <- paste(
    "Abbas R, Riley S, LaBadie RR, Bachinsky M, Chappell PB, Crownover PH,",
    "Damle B. A Thorough QT Study to Evaluate the Effects of a",
    "Supratherapeutic Dose of Sertraline on Cardiac Repolarization in",
    "Healthy Subjects. Clin Pharmacol Drug Dev. 2020;9(3):307-320.",
    "doi:10.1002/cpdd.749.",
    sep = " "
  )

  vignette <- "Abbas_2020_sertraline_QTc"

  units <- list(
    time = "h",
    dosing = "(none; PD-only model fed by an external sertraline plasma-concentration covariate)",
    concentration = "(observation QTcF is DeltaQTcF, the change from the period baseline in the Fridericia-corrected QT interval, ms; driving covariate CP_SERTRALINE_NGML is in ng/mL)"
  )

  covariateData <- list(
    CP_SERTRALINE_NGML = list(
      description = "Instantaneous total sertraline (parent) plasma concentration at the time of each ECG, supplied as a time-varying covariate from observed plasma samples or an external PK source.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives DeltaQTcF = e0 + slope * CP_SERTRALINE_NGML. Abbas 2020",
        "regressed on the observed LC-MS/MS plasma concentration drawn at",
        "the ECG time points on days 1 and 14; concentrations below the",
        "0.500 ng/mL quantitation limit were set to 0. Figure 3A spans",
        "0-500 ng/mL on day 14; geometric-mean Cmax 234.2 ng/mL (Table 3)."
      ),
      source_name = "Cijk (sertraline plasma concentration)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 54L,
    n_studies = 1L,
    n_observations = 1066L,
    age_range = "20-54 years (mean 36.5, SD 9.7)",
    weight_range = "mean 74.4 kg (SD 10.3); total body weight > 50 kg and BMI 17.5-30.5 kg/m^2 (inclusion); observed BMI 18.1-30.5 kg/m^2",
    sex_female_pct = 20.4,
    race_ethnicity = c(White = 85.2, Black = 5.6, Asian = 1.9, Other = 7.4),
    disease_state = "Healthy adult volunteers (no clinically significant cardiovascular history, screening QT <= 450 ms and QRS <= 120 ms).",
    dose_range = paste(
      "Sertraline titrated from 50 mg on day 1 to 200 mg twice daily",
      "(400 mg/day) on days 7-13, with a single 200 mg morning dose on",
      "day 14; matching placebo; single open-label moxifloxacin 400 mg on",
      "day 14 (not used in the exposure-response analysis)."
    ),
    regions = "Belgium (Pfizer Clinical Research Unit, Erasme Hospital, Brussels)",
    notes = paste(
      "See Abbas_2020_sertraline_ddQTcF.R for the full study description.",
      "Table 2 demographics are for the 54 enrolled subjects (43 men,",
      "11 women); the PK/PD data set held 1066 sertraline-QTcF pairs."
    )
  )

  ini({
    # Abbas 2020 Results 'Pharmacokinetic/Pharmacodynamic Analysis': slope
    # of DeltaQTcF versus sertraline concentration 0.036 ms/ng/mL (90%CI,
    # 0.029-0.043; P < .0001). The Figure 3A line digitises to 0.0362.
    slope <- 0.036
    label("Linear concentration-DeltaQTcF slope (ms per ng/mL)")

    # Not printed. Digitised from the solid regression line of Figure 3A
    # ('Observed and predicted (regression line) change from baseline in
    # QTcF versus sertraline'), calibrated on the 100-500 ng/mL and -40 to
    # 20 ms gridlines: the line crosses zero concentration at -5.9 ms
    # (one pixel = 0.23 ms). Population average of mu + TRT + t_k +
    # theta0 * (BQTc - Bbar) over the plotted day-14 records.
    e0 <- -5.9
    label("Population-averaged intercept on DeltaQTcF (ms)")

    # Random effects and residual error are not reported numerically.
    addSd <- fixed(0)
    label("Additive residual error SD on DeltaQTcF (ms; not reported in source)")
  })

  model({
    # Abbas 2020 Methods equation, collapsed to the population regression
    # line of Figure 3A.
    QTcF <- e0 + slope * CP_SERTRALINE_NGML
    QTcF ~ add(addSd)
  })
}
