Abbas_2020_sertraline_ddQTcF <- function() {
  description <- paste(
    "Concentration-DeltaDeltaQTcF linear mixed-effects PD model for the",
    "selective serotonin reuptake inhibitor sertraline in 54 healthy adults",
    "from a randomized, double-blind, placebo- and moxifloxacin-controlled,",
    "3-period crossover thorough-QT study (sertraline titrated to a",
    "supratherapeutic 400 mg/day, 200 mg twice daily, at steady state).",
    "The endpoint is the placebo-corrected change from baseline in the",
    "Fridericia-corrected QT interval (DeltaDeltaQTcF, ms). Typical-value",
    "structure: DeltaDeltaQTcF = e0 + slope * CP_SERTRALINE_NGML, with the",
    "printed slope 0.036 ms per ng/mL and e0 = 0.49 ms, the treatment-specific",
    "intercept back-solved from the paper's three model-predicted",
    "DeltaDeltaQTcF values (3.57 ms at 86 ng/mL, 8.55 ms at 223.5 ng/mL and",
    "8.93 ms at 234 ng/mL). Random effects and residual error are not",
    "reported, so the model is typical-value only. PD-only model: sertraline",
    "plasma concentration is supplied as a time-varying covariate; the paper",
    "characterises the PK by NCA only. Companion files",
    "Abbas_2020_sertraline_dQTcF.R (change from baseline, sertraline) and",
    "Abbas_2020_desmethylsertraline_dQTcF.R (change from baseline,",
    "N-desmethylsertraline).",
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
    concentration = "(observation QTcF is DeltaDeltaQTcF, the placebo-corrected change from baseline in the Fridericia-corrected QT interval, ms; driving covariate CP_SERTRALINE_NGML is in ng/mL)"
  )

  covariateData <- list(
    CP_SERTRALINE_NGML = list(
      description = "Instantaneous total sertraline (parent) plasma concentration at the time of each ECG, supplied as a time-varying covariate from observed plasma samples or an external PK source.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives DeltaDeltaQTcF = e0 + slope * CP_SERTRALINE_NGML. Abbas 2020",
        "regressed on the observed LC-MS/MS plasma concentration drawn at",
        "the ECG time points on days 1 and 14 (Methods 'Pharmacokinetics' and",
        "'Pharmacokinetic/Pharmacodynamic Analysis'); concentrations below",
        "the 0.500 ng/mL quantitation limit were set to 0. Day 14",
        "(400 mg/day) geometric-mean Cmax 234.2 ng/mL and Cmin 135.3 ng/mL",
        "(Table 3); mean therapeutic Cmax at 200 mg/day 86 ng/mL (Results).",
        "The paper assayed sertraline hydrochloride; concentrations are",
        "reported as ng/mL of the analyte without a salt correction."
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
      "Sertraline titrated from 50 mg on day 1 through 50-200 mg twice",
      "daily on days 2-6 to 200 mg twice daily (400 mg/day) on days 7-13,",
      "with a single 200 mg morning dose on day 14; matching placebo; and a",
      "single open-label moxifloxacin 400 mg dose on day 14. Three-period",
      "crossover with >= 14-day washouts. Moxifloxacin periods were not",
      "used in the exposure-response analysis."
    ),
    regions = "Belgium (Pfizer Clinical Research Unit, Erasme Hospital, Brussels)",
    notes = paste(
      "Table 2 demographics are for the 54 enrolled subjects (43 men, 11",
      "women). The PK/PD data set held 54 subjects and 1066",
      "sertraline-QTcF pairs plus 1050 QTcF records on placebo; ECGs from",
      "days 1 and 14 were used. Triplicate 12-lead ECGs at -1, -0.5 and 0 h",
      "pre-dose (the period baseline) and 1-24 h (day 1) or 1-72 h (day 14)",
      "post-dose. Fitted with nlme::lme (REML, unstructured random-effects",
      "covariance) in R 3.2.2."
    )
  )

  ini({
    # Abbas 2020 Results 'Pharmacokinetic/Pharmacodynamic Analysis': 'The
    # estimated slope of the regression for DeltaQTcF versus sertraline
    # plasma concentration was 0.036 ms/ng/mL (90%CI, 0.029-0.043;
    # P < .0001)' (also quoted in the Abstract).
    slope <- 0.036
    label("Linear concentration-DeltaDeltaQTcF slope (ms per ng/mL)")

    # Not printed. Back-solved from the three model predictions of mean
    # DeltaDeltaQTcF in the same Results paragraph: 3.57 ms at 86 ng/mL,
    # 8.93 ms at 234 ng/mL and 8.55 ms at 223.5 ng/mL. At the printed slope
    # the three anchors imply e0 = 0.474, 0.506 and 0.504 ms; 0.49 is the
    # centre of that range and reproduces every prediction within
    # 0.016 ms. The residual is the rounding of the printed slope: the
    # unrounded line through the anchors (slope 0.03622 ms per ng/mL,
    # e0 0.455 ms) reproduces all three to 0.005 ms, and the Figure 3A
    # regression line digitises to the same 0.0362 slope. e0 is the
    # treatment-specific intercept TRT of the Methods equation (sertraline
    # minus placebo); the time effects and the baseline term cancel in the
    # placebo-corrected contrast.
    e0 <- 0.49
    label("Treatment-specific intercept on DeltaDeltaQTcF (ms)")

    # The random intercept and random slope (unstructured covariance) are
    # not reported numerically, so no etas are declared. Residual error is
    # not reported either.
    addSd <- fixed(0)
    label("Additive residual error SD on DeltaDeltaQTcF (ms; not reported in source)")
  })

  model({
    # Abbas 2020 Methods 'Pharmacokinetic/Pharmacodynamic Analysis':
    #   DeltaQTc_ijk = (mu + TRT_j + t_k + eta_mu,i) + theta0 * (BQTc_i - Bbar)
    #                  + (theta1 + eta_C,i) * C_ijk + eps_ijk
    # Placebo-corrected, the typical contrast is TRT + theta1 * C.
    QTcF <- e0 + slope * CP_SERTRALINE_NGML
    QTcF ~ add(addSd)
  })
}
