Mukherjee_2018_tofacitinib_clinrem_cav <- function() {
  description <- paste(
    "Landmark logistic exposure-response model for the probability of",
    "clinical remission at week 8 (total Mayo score <= 2 with no",
    "individual subscore > 1) in adults with moderately to severely active",
    "ulcerative colitis given placebo or tofacitinib 0.5, 3, 10 or 15 mg",
    "twice daily in the phase 2 dose-ranging induction trial A3921063",
    "(Mukherjee 2018). The logit is a placebo intercept plus a drug effect",
    "that is an Emax function of steady-state average plasma concentration",
    "Cav,ss (ng/mL), Emax * X / (E50 + X), with E50 estimated on the log",
    "scale (0 for placebo), plus linear terms in baseline C-reactive",
    "protein and baseline total Mayo score. The paper prints the covariate",
    "coefficients without their scaling; the total Mayo score enters as",
    "SCORE_MAYO / 12 and C-reactive protein as CRP / 7 mg/L, as identified",
    "by the maintainers from the paper's placebo rates, subgroup",
    "predictions and confidence-interval widths (see the vignette). No",
    "between-subject variability and no residual error in the source",
    "(Bernoulli likelihood, NONMEM 7.3). AIC 190. One of nine",
    "exposure-response models (three endpoints by dose, Cav,ss and",
    "Ctrough,ss) in the Mukherjee_2018_tofacitinib_* family; the paper",
    "uses Cav,ss as its primary concentration metric."
  )
  reference <- paste(
    "Mukherjee A, Hazra A, Smith MK, Martin SW, Mould DR, Su C,",
    "Niezychowski W. Exposure-response characterization of tofacitinib",
    "efficacy in moderate to severe ulcerative colitis: Results from a",
    "dose-ranging phase 2 trial. Br J Clin Pharmacol.",
    "2018;84(6):1136-1145. doi:10.1111/bcp.13523. PMID: 29377257. PMCID:",
    "PMC5980553. Parameter values are the 'Clinical remission' block,",
    "'Cav,ss' column of Table 2."
  )
  vignette <- "Mukherjee_2018_tofacitinib"
  units <- list(
    time = "n/a (static landmark logistic regression at week 8; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "prob_clinrem (probability of clinical remission at week 8, 0-1; not a drug concentration)"
  )

  covariateData <- list(
    CAV = list(
      description = "Individual steady-state average tofacitinib plasma concentration over the 12 h twice-daily dosing interval (0 for placebo).",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Equation 1 of Mukherjee 2018: Cav,ss = Dose / (12 * CL/F), with the",
        "individual post hoc CL/F of modellib('Mukherjee_2018_tofacitinib')",
        "and the dose converted to ug so the result is in ng/mL (1000 *",
        "Dose[mg] / (12 * CL/F[L/h])). Set to 0 for placebo (Results, 'Data').",
        "Time-fixed per patient; CL/F did not change between baseline and week",
        "8."
      ),
      source_name = "Cav,ss"
    ),
    CRP = list(
      description = "Baseline C-reactive protein (standard assay).",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline-only value. Enters the logit as CRP / 7 mg/L, uncentered.",
        "The paper prints only the coefficient and does not state the scaling.",
        "Read per raw mg/L it predicts a placebo clinical-response rate of",
        "about 0.23 against the observed 0.43. A Fisher-information back-solve",
        "of the printed 95% CIs (Table 1 mean 13.9, SD 21 mg/L) puts the",
        "divisor at about 6-7.5 mg/L, close to the population median CRP. 7",
        "mg/L is the maintainers' rounded value. Using 10 (CRP in mg/dL)",
        "instead moves predicted response rates by less than 2 percentage",
        "points (see the vignette section on covariate scaling)."
      ),
      source_name = "baseline CRP"
    ),
    SCORE_MAYO = list(
      description = "Baseline total Mayo score for ulcerative colitis (0-12).",
      units = "(score, 0-12)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters the logit as SCORE_MAYO / 12 (the scale maximum), uncentered.",
        "The paper prints only the coefficient; the divisor is not stated and",
        "was identified by the maintainers (see the vignette section on",
        "covariate scaling): read per raw point the coefficient is",
        "implausible, and dividing by 12 reproduces the observed placebo rates",
        "of all three endpoints, the paper's subgroup prediction",
        "(placebo-adjusted remission at 10 mg twice daily of 39% versus 21%",
        "for baseline Mayo <= 8 versus > 8), and the printed 95% CI widths of",
        "the intercept and coefficient. Mean 8.2 (SD 1.6), range 3-12 in the",
        "analysed patients (Table 1)."
      ),
      source_name = "baseline total Mayo score"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 188L,
    n_studies = 1L,
    age_range = "19-77 years (Table 1)",
    age_mean = "42.1 years (SD 13.6) (Table 1)",
    weight_range = "41.4-126 kg (Table 1)",
    weight_mean = "74.8 kg (SD 15.5) (Table 1)",
    sex_female_pct = 45.2,
    disease_state = paste(
      "Moderately to severely active ulcerative colitis (total Mayo score",
      "6-12, endoscopic subscore >= 2, diagnosis >= 3 months); baseline total",
      "Mayo score mean 8.2 (range 3-12), CRP mean 13.9 mg/L, 43.1% on",
      "corticosteroids and 18.6% with prior TNF-inhibitor failure (Table 1)."
    ),
    dose_range = "Placebo (n = 47) or tofacitinib 0.5 (n = 31), 3 (n = 32), 10 (n = 30) or 15 mg (n = 48) orally twice daily for 8 weeks.",
    regions = "Multinational phase 2 trial A3921063 (NCT00787202).",
    notes = paste(
      "188 of 194 treated patients (five without plasma concentrations and",
      "one never treated were excluded). Patients missing a covariate in a",
      "model were excluded from that model; no imputation. Fitted by maximum",
      "likelihood in NONMEM 7.3."
    )
  )

  ini({
    # Mukherjee 2018 Methods: logit p = beta0 + betaD + sum(beta_m * X_m),
    # with betaD = Emax * C / (E50 + C) or betaD = slope * C. Table 2 prints
    # each value with its 95% CI; asymmetric CIs on E50 and on the slope show
    # those two were estimated on the log scale (Table 2 footnote: 'E50 ...
    # (log scale)'), so they are carried as log parameters here.
    logite0 <- 1.83
    label("Logit of clinical remission at zero exposure with SCORE_MAYO = 0 and CRP = 0 (uncentered intercept; unitless)")  # Table 2, Clinical remission, Cav,ss column, 'Placebo' = 1.83 (95% CI -0.33-3.99)
    emax_clinrem <- 2.35
    label("Maximum tofacitinib effect on the clinical remission logit (unitless)")  # Table 2, Clinical remission, Cav,ss column, 'E max' = 2.35 (95% CI 1.25-3.45)
    lec50_clinrem <- log(3.86)
    label("Log of the steady-state average concentration Cav,ss giving half the maximum effect (log ng/mL)")  # Table 2, Clinical remission, Cav,ss column, 'E 50 (mg or ng ml-1)' = 3.86 (95% CI 0.63-23.55), estimated on the log scale
    e_crp_clinrem <- -0.105
    label("Log-odds shift per unit of CRP / 7 mg/L (unitless)")  # Table 2, Clinical remission, Cav,ss column, 'Coefficient for baseline CRP' = -0.105 (95% CI -0.26-0.05)
    e_score_mayo_clinrem <- -6
    label("Log-odds shift per unit of SCORE_MAYO / 12 (unitless)")  # Table 2, Clinical remission, Cav,ss column, 'Coefficient for baseline total Mayo score' = -6 (95% CI -9.25 to -2.75)

    # The source maximises a Bernoulli likelihood on a 0/1 outcome; this
    # placeholder residual only lets nlmixr2 accept the model for simulation.
    addSd_prob_clinrem <- fixed(0.001)
    label("Placeholder additive residual SD on prob_clinrem; not from the source (unitless)")
  })

  model({
    # SCORE_MAYO is divided by 12 (the scale maximum) and CRP by 7 mg/L; the
    # paper prints the coefficients without the scaling. See the vignette
    # section on covariate scaling for how both divisors were identified.
    drug <- emax_clinrem * CAV / (exp(lec50_clinrem) + CAV)
    logit_clinrem <- logite0 +
      drug +
      e_crp_clinrem * CRP / 7 +
      e_score_mayo_clinrem * SCORE_MAYO / 12
    prob_clinrem <- expit(logit_clinrem)
    prob_clinrem ~ add(addSd_prob_clinrem)
  })
}
