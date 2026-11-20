Wang_2023_polymyxinB_aki_cmid <- function() {
  description <- paste0(
    "Landmark univariate logistic-regression exposure-safety model for the ",
    "occurrence of acute kidney injury (KDIGO creatinine or urine-output ",
    "criteria) during intravenous polymyxin B therapy in critically ill ",
    "adults (Wang 2023; n = 55 intensive-care patients, 30 with acute ",
    "kidney injury). The probability of acute kidney injury is ",
    "1 / (1 + exp(2.742 - 0.845 * CMID_INTERVAL)), where CMID_INTERVAL is ",
    "the measured steady-state polymyxin B concentration 6 h after the ",
    "start of a q12h dose (the paper's mid-dosing-interval concentration ",
    "C1/2t) in ug/mL; the probability reaches 0.5 at 3.24 ug/mL. The paper ",
    "offers it as a single-sample alternative to its final two-predictor ",
    "model, Wang_2023_polymyxinB_aki, because C1/2t was the best ",
    "single-timepoint predictor (area under the ROC curve 0.747). There is ",
    "no PK layer and no ODE: exposure enters as a per-subject data column. ",
    "No between-subject random effect and no residual error are estimated ",
    "(Bernoulli likelihood)."
  )
  reference <- paste(
    "Wang J, Li Y, Huang S, Wang M, Jin L, Luo X, Cheng X, Yang N, Zhu H.",
    "Mid-dosing interval concentration is important for polymyxin B",
    "exposure and acute kidney injury in critically ill patients.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12:1911-1921.",
    "doi:10.1002/psp4.13040.",
    sep = " "
  )
  vignette <- "Wang_2023_polymyxinB"
  units <- list(
    time = "n/a (static landmark exposure-response regression evaluated once per subject over the polymyxin B course; no time dimension)",
    dosing = "n/a (no dose events; the dose enters only through the CMID_INTERVAL exposure column)",
    concentration = "prob_aki (probability of acute kidney injury during polymyxin B therapy, 0-1; also logit_aki)"
  )

  covariateData <- list(
    CMID_INTERVAL = list(
      description = "Measured steady-state total polymyxin B plasma concentration 6 h after the start of a q12h dose (the paper's mid-dosing-interval concentration C1/2t).",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters the logit LINEARLY and UNCENTRED, slope 0.845 per ug/mL",
        "(Wang 2023 Equation 5). Total polymyxin B, the sum of polymyxin",
        "B1 and B2 by UPLC-MS/MS. Sampled after at least 48 h of therapy",
        "(steady state) under a q12h regimen with infusions of 1 h (87.3",
        "percent of patients), 2 h or 2.5 h. Observed range 1.438-7.030",
        "ug/mL; mean +/- SD 3.6 +/- 1.3 ug/mL overall, 3.0 +/- 0.9 without",
        "and 4.1 +/- 1.3 with acute kidney injury (Table 1). The",
        "Youden-optimal ROC cutoff is 3.47 ug/mL (sensitivity 66.7",
        "percent, specificity 72 percent); the model probability crosses",
        "0.5 at 2.742 / 0.845 = 3.24 ug/mL, the value quoted in the",
        "Results."
      ),
      source_name = "C1/2t(B) (mid-dosing interval concentration)"
    )
  )

  covariatesDataExcluded <- list(
    CMAX = list(
      description = "Measured steady-state polymyxin B peak concentration at the end of the infusion.",
      units = "ug/mL",
      type = "continuous",
      notes = "Rival single-concentration predictor; area under the ROC curve 0.675 against 0.747 for C1/2t (Figure 3b). No single-concentration logistic coefficients are published for it."
    ),
    CTROUGH = list(
      description = "Measured steady-state polymyxin B trough concentration before the next dose.",
      units = "ug/mL",
      type = "continuous",
      notes = "Rival single-concentration predictor; area under the ROC curve 0.672 against 0.747 for C1/2t (Figure 3b). No single-concentration logistic coefficients are published for it."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 55L,
    n_studies = 1L,
    n_observations = "55 binary records, one per subject; 30 subjects (54.5 percent) developed acute kidney injury, 3-17 days after the start of polymyxin B (Wang 2023 Results, Figure 2).",
    age_range = "Adults over 18 years; mean +/- SD 61.0 +/- 16.7 years (Table 1).",
    weight_range = "Median (IQR) 62.0 (55.0, 73.0) kg (Table 1).",
    sex_female_pct = 36.4,
    race_ethnicity = "Chinese (single centre, Nanjing Drum Tower Hospital, China).",
    disease_state = paste(
      "Critically ill intensive-care adults with confirmed multidrug-resistant",
      "Gram-negative infections; APACHE II 20.6 +/- 5.8 and SOFA 6.3 +/- 3.0",
      "at baseline. Baseline serum creatinine above 1.5 mg/dL and renal",
      "replacement therapy before polymyxin B were exclusion criteria."
    ),
    dose_range = "Polymyxin B sulfate intravenously every 12 h, maintenance 100-200 mg per day in two divided doses (mostly 150 mg/day; Methods: 2.5-3 mg/kg/day), loading dose 75-150 mg in 34 patients (Table 1).",
    regions = "China",
    notes = "Same 55-patient cohort as Wang_2023_polymyxinB_aki. Univariate logistic regression (SPSS 25), plotted in Wang 2023 Figure 4."
  )

  ini({
    # Wang 2023 Equation 5 (univariate model on C1/2t):
    #   ProbabilityAKI = 1 / (1 + exp(2.742 - 0.845 * C1/2t))
    # so logit(p) = -2.742 + 0.845 * C1/2t. A logistic refit of the
    # deposited patient-level dataset (Data S1) returns -2.7420 and 0.8450.
    logit_ref <- -2.742; label("Logit of the acute-kidney-injury probability at zero mid-dosing-interval polymyxin B concentration (unitless logit)") # Wang 2023 Equation 5: constant 2.742 inside exp(), i.e. intercept -2.742
    e_cmid_logit <- 0.845; label("Log-odds of acute kidney injury per ug/mL of steady-state mid-dosing-interval polymyxin B concentration (unitless logit)") # Wang 2023 Equation 5: 0.845 x C1/2t

    # Bernoulli likelihood in the source; placeholder residual only.
    addSd_prob_aki <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_aki <- logit_ref + e_cmid_logit * CMID_INTERVAL

    prob_aki <- expit(logit_aki)

    prob_aki ~ add(addSd_prob_aki)
  })
}
