Yao_2025_acalabrutinib_teae_grade2 <- function() {
  description <- paste(
    "Binomial logistic-regression exposure-safety model for any grade >= 2",
    "adverse event as a function of the individual steady-state daily",
    "acalabrutinib AUC(0-24) in adults with B-cell malignancies treated with",
    "acalabrutinib 100-250 mg QD or 100-200 mg BID (Yao 2025, seven pooled",
    "global and China studies, n = 476). The probability is",
    "expit(1.851 + 0.0003197 * AUC_ACALABRUTINIB), with the exposure in",
    "ng*h/mL entering linearly and uncentred. The relationship is NOT",
    "statistically significant (P = 0.12 annotated on Figure 6). The paper",
    "prints no coefficient; intercept and slope were recovered from the",
    "fitted curve of Figure 6 (left panel). There is no PK layer and no",
    "ODE: the exposure is supplied as a data column, derived in the source",
    "from the post hoc estimates of the population PK model packaged as",
    "Yao_2025_acalabrutinib. Companion of Yao_2025_acalabrutinib_teae_grade3."
  )
  reference <- paste(
    "Yao T, Ding J, Sharma S, Li Y, Xu Q, Ma P. Assessment of ethnic differences",
    "in pharmacokinetics and clinical responses of acalabrutinib between",
    "Chinese and White patients with B-cell malignancies.",
    "Br J Clin Pharmacol. 2025;91(7):2008-2019. doi:10.1002/bcp.70018"
  )
  vignette <- "Yao_2025_acalabrutinib"
  units <- list(
    time = "n/a (static landmark exposure-response regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the AUC_ACALABRUTINIB covariate column)",
    concentration = "prob_teae_grade2 (probability of any grade >= 2 adverse event, 0-1; also logit_teae_grade2)"
  )

  covariateData <- list(
    AUC_ACALABRUTINIB = list(
      description = "Individual steady-state acalabrutinib AUC over 24 h (two BID dosing intervals or one QD interval), from the post hoc estimates of the companion population PK model",
      units = "ng*h/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Yao 2025 Methods 'E-R analyses': AUC24h,ss 'derived from post hoc",
        "estimates'. Compute from modellib('Yao_2025_acalabrutinib') as",
        "1000 * (total daily dose, mg) * F1 / (individual CL/F, L/h).",
        "Enters LINEARLY and uncentred. Figure 6 plots the fitted line over",
        "about 200-7600 ng*h/mL; the 100 mg BID group sits near 1700",
        "ng*h/mL."
      ),
      source_name = "Acalabrutinib AUC24h,ss"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 476L,
    n_studies = 7L,
    n_observations = "476 binary records, one per patient",
    age_range = "not reported for the exposure-safety cohort",
    weight_range = "not reported for the exposure-safety cohort",
    sex_female_pct = NA_real_,
    race_ethnicity = "Pooled global and Chinese patients; Supplementary Table S1 lists 105 China-study patients (Results text: 66 Chinese R/R CLL/SLL and 34 R/R MCL)",
    disease_state = "Adults with B-cell malignancies: CLL 390, MCL 57, DLBCL 16, FL 11, SLL 2 (Supplementary Table S1)",
    dose_range = "Acalabrutinib 100 mg BID (440), 100 mg QD (10), 175 mg QD (3), 200 mg BID (18), 200 mg QD (4), 250 mg QD (1) (Supplementary Table S1)",
    regions = "Global (ACE-CL-001, ACE-CL-003, ACE-CL-007, ACE-LY-002, ACE-LY-003, ACE-LY-004) and China (D8220C00007)",
    notes = "Yao 2025 Methods 'Data' and Results 'E-S analysis'. Safety endpoint: incidence of any grade >= 2 adverse event."
  )

  ini({
    # Figure 6 (left panel, 'P = 0.12'): logit(p) = logit_ref + e_auc_logit *
    # AUC_ACALABRUTINIB. No coefficient is printed; the fitted line was
    # digitised from the figure image (135 pixel columns, AUC 195-7650
    # ng*h/mL) and fitted as a straight line on the logit scale, which
    # reproduces every digitised point within 0.003 in probability.
    logit_ref <- 1.851; label("Logit of the probability of any grade >= 2 adverse event at zero acalabrutinib exposure (unitless logit)") # figure-derived: Figure 6 left panel fitted line; p = 0.864 at AUC 0
    e_auc_logit <- 0.0003197; label("Log-odds of any grade >= 2 adverse event per 1 ng*h/mL increase in AUC_ACALABRUTINIB (unitless logit)") # figure-derived: Figure 6 left panel fitted line (no coefficient printed); panel annotates P = 0.12

    # The source is a binomial logistic regression (no random effects, no
    # residual). The tiny additive residual exists only so rxode2 has an
    # error model to attach to the probability.
    addSd_prob_teae_grade2 <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # Exposure enters linearly and uncentred, so logit_ref is the logit at
    # zero exposure, an extrapolation below the observed range.
    logit_teae_grade2 <- logit_ref + e_auc_logit * AUC_ACALABRUTINIB
    prob_teae_grade2 <- expit(logit_teae_grade2)
    prob_teae_grade2 ~ add(addSd_prob_teae_grade2)
  })
}
