Wang_2023_polymyxinB_aki <- function() {
  description <- paste0(
    "Landmark multivariable logistic-regression exposure-safety model for ",
    "the occurrence of acute kidney injury (KDIGO creatinine or urine-output ",
    "criteria) during intravenous polymyxin B therapy in critically ill ",
    "adults (Wang 2023; n = 55 intensive-care patients, 30 with acute ",
    "kidney injury). The probability of acute kidney injury is ",
    "1 / (1 + exp(4.623 - 0.116 * WBC - 0.040 * AUC_POLYMYXINB)), where WBC ",
    "is the baseline white blood cell count in 10^9/L and AUC_POLYMYXINB is ",
    "the steady-state polymyxin B AUC0-24 in ug*h/mL estimated by a ",
    "one-compartment first-order equation from the measured peak and ",
    "mid-dosing-interval concentrations. This is the paper's final model, ",
    "selected by forward stepwise logistic regression. There is no PK layer ",
    "and no ODE: exposure enters as a per-subject data column. No ",
    "between-subject random effect and no residual error are estimated ",
    "(Bernoulli likelihood). The single-concentration alternative from the ",
    "same paper is Wang_2023_polymyxinB_aki_cmid."
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
    dosing = "n/a (no dose events; the dose enters only through the AUC_POLYMYXINB exposure column)",
    concentration = "prob_aki (probability of acute kidney injury during polymyxin B therapy, 0-1; also logit_aki)"
  )

  covariateData <- list(
    WBC = list(
      description = "White blood cell count at baseline, before the start of polymyxin B.",
      units = "10^9/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline only (Wang 2023 Methods: baseline routine laboratory",
        "examination data). Enters the logit LINEARLY and UNCENTRED, odds",
        "ratio 1.123 per 10^9/L (95% CI 1.000-1.260, p = 0.0499; Table 2).",
        "Observed median (IQR) 9.6 (7.3, 13.4) x 10^9/L overall, 8.5 (5.8,",
        "11.3) without and 10.8 (8.1, 17.2) with acute kidney injury",
        "(Table 1). On its own it discriminates poorly (area under the ROC",
        "curve 0.661)."
      ),
      source_name = "White blood cell (WBC)"
    ),
    AUC_POLYMYXINB = list(
      description = "Steady-state polymyxin B AUC over 24 h estimated by the one-compartment first-order (Pai) equation from the measured peak and mid-dosing-interval concentrations (the paper's AUC1/2t,max).",
      units = "ug*h/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters the logit LINEARLY and UNCENTRED, odds ratio 1.040 per",
        "ug*h/mL (95% CI 1.013-1.068, p = 0.003; Table 2), so the intercept",
        "is an extrapolated anchor at zero exposure. Total polymyxin B,",
        "measured as the sum of polymyxin B1 and B2 by UPLC-MS/MS. Derived",
        "per subject from two steady-state samples drawn after at least 48",
        "h of q12h therapy: Cmax at the end of the infusion and C1/2t 6 h",
        "after the start of the dose (Wang 2023 Equations 2 and 3, after Pai",
        "2014): Ke = ln(Cmax / C1/2t) / dt, the trough is extrapolated as",
        "Cmax * exp(-Ke * (tau - t_inf)), the interval AUC is",
        "t_inf * (Cmax + Ctrough) / 2 + (Cmax - Ctrough) / Ke, and the",
        "12-h value is doubled to a 24-h value. The values in the deposited",
        "dataset (Data S1) are reproduced to within 1.0 percent with",
        "dt = 5 h for every subject, including those infused over 2 or",
        "2.5 h. Observed mean +/- SD 91.9 +/- 30.1 ug*h/mL overall,",
        "77.7 +/- 21.8 without and 103.7 +/- 31.2 with acute kidney injury",
        "(Table 1). A popPK-derived steady-state AUC0-24 can be supplied in",
        "its place, but the regression was fitted to this two-sample",
        "estimate, which runs below the three-point trapezoidal AUC",
        "(99.1 ug*h/mL)."
      ),
      source_name = "AUC1/2t,max (AUC0-24, first-order pharmacokinetic equation, based on C1/2t(B) and Cmax(B))"
    )
  )

  covariatesDataExcluded <- list(
    CMID_INTERVAL = list(
      description = "Measured steady-state polymyxin B concentration 6 h after the start of a q12h dose (C1/2t).",
      units = "ug/mL",
      type = "continuous",
      notes = paste(
        "Significant in univariate analysis (p = 0.001) but not retained in",
        "the multivariable model once the AUC it nearly determines entered",
        "(Pearson r = 0.985 with AUC1/2t,max; Figure 5). It is the sole",
        "driver of the companion univariate model",
        "Wang_2023_polymyxinB_aki_cmid."
      )
    ),
    CMAX = list(
      description = "Measured steady-state polymyxin B peak concentration at the end of the infusion.",
      units = "ug/mL",
      type = "continuous",
      notes = "Significant in univariate analysis (p = 0.014; Table 1) and not retained by forward stepwise selection; single-concentration area under the ROC curve 0.675 (Figure 3b)."
    ),
    CTROUGH = list(
      description = "Measured steady-state polymyxin B trough concentration before the next dose.",
      units = "ug/mL",
      type = "continuous",
      notes = "Significant in univariate analysis (p = 0.019; Table 1) and not retained by forward stepwise selection; single-concentration area under the ROC curve 0.672 (Figure 3b)."
    ),
    NEUT = list(
      description = "Absolute neutrophil count at baseline.",
      units = "10^9/L",
      type = "continuous",
      notes = "Not significant in univariate analysis (p = 0.106; Table 1). The neutrophil PERCENTAGE was significant (p = 0.041) and, like the maintenance dose (p = 0.031), entered the stepwise selection without being retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 55L,
    n_studies = 1L,
    n_observations = "55 binary records, one per subject; 30 subjects (54.5 percent) developed acute kidney injury, 3-17 days after the start of polymyxin B (Wang 2023 Results, Figure 2). 165 polymyxin B concentrations (three per subject) were measured to derive the exposure metrics.",
    age_range = "Adults over 18 years; mean +/- SD 61.0 +/- 16.7 years (Table 1).",
    weight_range = "Median (IQR) 62.0 (55.0, 73.0) kg (Table 1).",
    sex_female_pct = 36.4,
    race_ethnicity = "Chinese (single centre, Nanjing Drum Tower Hospital, China).",
    disease_state = paste(
      "Critically ill intensive-care adults with confirmed multidrug-resistant",
      "Gram-negative infections (Acinetobacter baumannii 69.1 percent,",
      "Klebsiella pneumoniae 40.0 percent), mostly pulmonary (69.1",
      "percent); APACHE II 20.6 +/- 5.8 and SOFA 6.3 +/- 3.0 at baseline.",
      "Baseline serum creatinine above 1.5 mg/dL and renal replacement",
      "therapy before polymyxin B were exclusion criteria."
    ),
    dose_range = paste(
      "Polymyxin B sulfate intravenously every 12 h, maintenance 100 mg/day",
      "(27.3 percent), 150 mg/day (70.9 percent) or 200 mg/day (1.8",
      "percent) given as two divided doses (Methods: 2.5-3 mg/kg/day),",
      "infused over 1 h (87.3 percent), 2 h or 2.5 h; a loading dose of",
      "75-150 mg in 34 patients (Table 1). Treatment duration 4-42 days."
    ),
    regions = "China",
    renal_function = "Baseline creatinine clearance (Cockcroft-Gault) 129.6 +/- 62.8 mL/min; serum creatinine median 52.0 umol/L (Table 1).",
    notes = paste(
      "Single-centre real-world cohort, August 2020 to July 2022. Acute",
      "kidney injury follows KDIGO: serum creatinine rise >= 0.3 mg/dL within",
      "48 h, > 1.5-fold within 7 days, or urine output < 0.5 mL/kg/h for",
      "> 6 h. Analysed in SPSS 25 with univariate screening (p < 0.05) then",
      "forward stepwise multivariable logistic regression; Hosmer-Lemeshow",
      "p = 0.749. Area under the ROC curve 0.805 (95% CI 0.690-0.921); the",
      "Youden-optimal cutoff on the predicted probability is 0.57",
      "(sensitivity 70 percent, specificity 80 percent)."
    )
  )

  ini({
    # ==================================================================
    # Wang 2023 Equation 4 (final multivariable model):
    #   ProbabilityAKI = 1 / (1 + exp(4.623 - 0.116 * WBC - 0.040 * AUC1/2t,max))
    # so logit(p) = -4.623 + 0.116 * WBC + 0.040 * AUC. Table 2 prints the
    # same two slopes as odds ratios: 1.123 per 10^9/L WBC
    # (log 1.123 = 0.116) and 1.040 per ug*h/mL AUC (log 1.040 = 0.039).
    # A logistic refit of the deposited patient-level dataset (Data S1)
    # returns -4.6234, 0.1157 and 0.0396.
    # ==================================================================

    logit_ref <- -4.623; label("Logit of the acute-kidney-injury probability at zero WBC and zero polymyxin B AUC (unitless logit)") # Wang 2023 Equation 4: constant 4.623 inside exp(), i.e. intercept -4.623
    e_wbc_logit <- 0.116; label("Log-odds of acute kidney injury per 10^9/L of baseline white blood cell count (unitless logit)") # Wang 2023 Equation 4: 0.116 x WBC; Table 2 OR 1.123 (95% CI 1.000-1.260)
    e_auc_logit <- 0.040; label("Log-odds of acute kidney injury per ug*h/mL of steady-state polymyxin B AUC0-24 (unitless logit)") # Wang 2023 Equation 4: 0.040 x AUC1/2t,max; Table 2 OR 1.040 (95% CI 1.013-1.068)

    # No between-subject variability and no residual error: the source is a
    # logistic regression with a Bernoulli likelihood. The fixed additive
    # residual below exists only so that rxode2 has an error model to
    # attach to the probability; it is not a published quantity.
    addSd_prob_aki <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_aki <- logit_ref + e_wbc_logit * WBC + e_auc_logit * AUC_POLYMYXINB

    prob_aki <- expit(logit_aki)

    prob_aki ~ add(addSd_prob_aki)
  })
}
