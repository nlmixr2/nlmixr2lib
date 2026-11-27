Abbas_2020_desmethylsertraline_dQTcF <- function() {
  description <- paste(
    "Concentration-DeltaQTcF linear mixed-effects PD model for",
    "N-desmethylsertraline, the principal metabolite of sertraline, in 54",
    "healthy adults from a randomized, double-blind, placebo- and",
    "moxifloxacin-controlled, 3-period crossover thorough-QT study",
    "(sertraline titrated to a supratherapeutic 400 mg/day, 200 mg twice",
    "daily, at steady state). The metabolite model was fitted separately",
    "from the parent model. The endpoint is the change from baseline in the",
    "Fridericia-corrected QT interval (DeltaQTcF, ms), not",
    "placebo-corrected. Typical-value structure: DeltaQTcF = e0 + slope *",
    "CP_DESMETHYLSERTRALINE_NGML, the population regression line of the",
    "paper's Figure 3B, with the printed slope 0.036 ms per ng/mL and",
    "e0 = -6.4 ms digitised from the line's value at zero concentration.",
    "Random effects and residual error are not reported, so the model is",
    "typical-value only. PD-only model: N-desmethylsertraline plasma",
    "concentration is supplied as a time-varying covariate. Companion files",
    "Abbas_2020_sertraline_dQTcF.R and Abbas_2020_sertraline_ddQTcF.R",
    "(parent sertraline).",
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
    dosing = "(none; PD-only model fed by an external N-desmethylsertraline plasma-concentration covariate)",
    concentration = "(observation QTcF is DeltaQTcF, the change from the period baseline in the Fridericia-corrected QT interval, ms; driving covariate CP_DESMETHYLSERTRALINE_NGML is in ng/mL)"
  )

  covariateData <- list(
    CP_DESMETHYLSERTRALINE_NGML = list(
      description = "Instantaneous total N-desmethylsertraline (the principal sertraline metabolite) plasma concentration at the time of each ECG, supplied as a time-varying covariate from observed plasma samples or an external PK source.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives DeltaQTcF = e0 + slope * CP_DESMETHYLSERTRALINE_NGML. Abbas",
        "2020 regressed on the observed LC-MS/MS plasma concentration drawn",
        "at the ECG time points on days 1 and 14; concentrations below the",
        "0.500 ng/mL quantitation limit were set to 0. Day 14 (400 mg/day",
        "sertraline) geometric-mean Cmax 249.6 ng/mL and Cmin 184.1 ng/mL",
        "(Table 3); Figure 3B spans 0-450 ng/mL."
      ),
      source_name = "Cijk (N-desmethylsertraline plasma concentration)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 54L,
    n_studies = 1L,
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
      "The paper does not print the number of metabolite-QTcF pairs."
    )
  )

  ini({
    # Abbas 2020 Results 'Pharmacokinetic/Pharmacodynamic Analysis': 'for
    # DeltaQTcF versus N-desmethylsertraline plasma concentration, the
    # estimated slope was 0.036 ms/ng/mL (90%CI, 0.029-0.042; P < .0001)'.
    # The Figure 3B line digitises to 0.0356.
    slope <- 0.036
    label("Linear concentration-DeltaQTcF slope (ms per ng/mL)")

    # Not printed. Digitised from the solid regression line of Figure 3B,
    # calibrated on the 100-400 ng/mL and -40 to 20 ms gridlines: the line
    # crosses zero concentration at -6.4 ms (one pixel = 0.23 ms).
    # Population average of mu + TRT + t_k + theta0 * (BQTc - Bbar) over
    # the plotted day-14 records.
    e0 <- -6.4
    label("Population-averaged intercept on DeltaQTcF (ms)")

    # Random effects and residual error are not reported numerically.
    addSd <- fixed(0)
    label("Additive residual error SD on DeltaQTcF (ms; not reported in source)")
  })

  model({
    # Abbas 2020 Methods equation, fitted with the metabolite concentration
    # and collapsed to the population regression line of Figure 3B.
    QTcF <- e0 + slope * CP_DESMETHYLSERTRALINE_NGML
    QTcF ~ add(addSd)
  })
}
