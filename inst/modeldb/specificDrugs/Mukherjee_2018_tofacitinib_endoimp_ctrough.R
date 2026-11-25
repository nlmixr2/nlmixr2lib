Mukherjee_2018_tofacitinib_endoimp_ctrough <- function() {
  description <- paste(
    "Landmark logistic exposure-response model for the probability of",
    "mucosal healing at week 8 (Mayo endoscopic subscore of 0 or 1,",
    "derived post hoc) in adults with moderately to severely active",
    "ulcerative colitis given placebo or tofacitinib 0.5, 3, 10 or 15 mg",
    "twice daily in the phase 2 dose-ranging induction trial A3921063",
    "(Mukherjee 2018). The logit is a placebo intercept plus a drug effect",
    "that is an Emax function of steady-state trough plasma concentration",
    "Ctrough,ss (ng/mL), Emax * X / (E50 + X), with E50 estimated on the",
    "log scale (0 for placebo), plus linear terms in baseline total Mayo",
    "score and prior TNF-inhibitor failure. The paper prints the covariate",
    "coefficients without their scaling; the total Mayo score enters as",
    "SCORE_MAYO / 12, as identified by the maintainers from the paper's",
    "placebo rates, subgroup predictions and confidence-interval widths",
    "(see the vignette). No between-subject variability and no residual",
    "error in the source (Bernoulli likelihood, NONMEM 7.3). AIC 236. One",
    "of nine exposure-response models (three endpoints by dose, Cav,ss and",
    "Ctrough,ss) in the Mukherjee_2018_tofacitinib_* family; the paper",
    "uses Cav,ss as its primary concentration metric."
  )
  reference <- paste(
    "Mukherjee A, Hazra A, Smith MK, Martin SW, Mould DR, Su C,",
    "Niezychowski W. Exposure-response characterization of tofacitinib",
    "efficacy in moderate to severe ulcerative colitis: Results from a",
    "dose-ranging phase 2 trial. Br J Clin Pharmacol.",
    "2018;84(6):1136-1145. doi:10.1111/bcp.13523. PMID: 29377257. PMCID:",
    "PMC5980553. Parameter values are the 'Mucosal healing' block,",
    "'Ctrough,ss' column of Table 2."
  )
  vignette <- "Mukherjee_2018_tofacitinib"
  units <- list(
    time = "n/a (static landmark logistic regression at week 8; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CTROUGH covariate column)",
    concentration = "prob_endoimp (probability of mucosal healing at week 8, 0-1; not a drug concentration)"
  )

  covariateData <- list(
    CTROUGH = list(
      description = "Individual steady-state trough tofacitinib plasma concentration at the end of the 12 h twice-daily dosing interval (0 for placebo).",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "TOTAL plasma concentration at steady state, from Equation 2 of",
        "Mukherjee 2018 (the closed-form one-compartment",
        "first-order-absorption steady-state concentration at tau = 12 h,",
        "absorption lag ignored) with the individual post hoc Ka, CL/F and V/F",
        "of modellib('Mukherjee_2018_tofacitinib'). Set to 0 for placebo",
        "(Results, 'Data'). Time-fixed per patient."
      ),
      source_name = "Ctrough,ss"
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
    ),
    BIOLOGIC_IR = list(
      description = "Prior tumour necrosis factor inhibitor (TNFi) failure within the previous 12 months (1 = yes, 0 = no).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no prior TNFi failure).",
      notes = paste(
        "The source covariate restricts the failure to anti-TNF agents within",
        "the previous 12 months (Table 1: 35/188, 18.6%); BIOLOGIC_IR is",
        "reused with that restriction recorded here. Log-odds shift on the",
        "mucosal-healing logit; its 95% CI includes 0 (Table 2)."
      ),
      source_name = "prior TNFi failure"
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
    logite0 <- 3.12
    label("Logit of mucosal healing at zero exposure with SCORE_MAYO = 0 and BIOLOGIC_IR = 0 (uncentered intercept; unitless)")  # Table 2, Mucosal healing, Ctrough,ss column, 'Placebo' = 3.12 (95% CI 1.27-4.97)
    emax_endoimp <- 1.41
    label("Maximum tofacitinib effect on the mucosal healing logit (unitless)")  # Table 2, Mucosal healing, Ctrough,ss column, 'E max' = 1.41 (95% CI 0.41-2.41)
    lec50_endoimp <- log(1.71)
    label("Log of the steady-state trough concentration Ctrough,ss giving half the maximum effect (log ng/mL)")  # Table 2, Mucosal healing, Ctrough,ss column, 'ED 50' = 1.71 (95% CI 0.13-21.85), estimated on the log scale
    e_score_mayo_endoimp <- -5.35
    label("Log-odds shift per unit of SCORE_MAYO / 12 (unitless)")  # Table 2, Mucosal healing, Ctrough,ss column, 'Baseline total Mayo score' = -5.35 (95% CI -7.98 to -2.72)
    e_biologic_ir_endoimp <- -0.789
    label("Log-odds shift for prior TNF-inhibitor failure (unitless)")  # Table 2, Mucosal healing, Ctrough,ss column, 'Prior TNFi failure' = -0.789 (95% CI -1.62-0.04)

    # The source maximises a Bernoulli likelihood on a 0/1 outcome; this
    # placeholder residual only lets nlmixr2 accept the model for simulation.
    addSd_prob_endoimp <- fixed(0.001)
    label("Placeholder additive residual SD on prob_endoimp; not from the source (unitless)")
  })

  model({
    # SCORE_MAYO is divided by 12 (the scale maximum); the paper prints the
    # coefficient without the scaling. See the vignette section on
    # covariate scaling for how the divisor was identified.
    drug <- emax_endoimp * CTROUGH / (exp(lec50_endoimp) + CTROUGH)
    logit_endoimp <- logite0 +
      drug +
      e_score_mayo_endoimp * SCORE_MAYO / 12 +
      e_biologic_ir_endoimp * BIOLOGIC_IR
    prob_endoimp <- expit(logit_endoimp)
    prob_endoimp ~ add(addSd_prob_endoimp)
  })
}
