Biesdorf_2026_cofetuzumabPelidotin_orr <- function() {
  description <- paste0(
    "Logistic-regression exposure-efficacy model for the confirmed ",
    "objective response rate (ORR, RECIST 1.1 best overall response) in ",
    "adults with PTK7-expressing recurrent non-small cell lung cancer ",
    "treated with cofetuzumab pelidotin (Cofe-P, ABBV-647) 2.8 mg/kg ",
    "every 3 weeks (Biesdorf 2026, n = 60, phase 1b study M19-611). The ",
    "probability of response is expit(a + (b + b_target * ",
    "PTK7_HIGH_NONSQUAMOUS_EGFR_WT) * CAV), linear in CAV, the individual ",
    "average Cofe-P CONJUGATE serum concentration up to the time of best overall ",
    "response, in ug/mL. Membership of the target population ",
    "(PTK7-high, >= 90% of tumour cells at >= 2+ staining intensity; ",
    "non-squamous histology; EGFR wild-type) steepens the exposure ",
    "slope. There is no PK layer and no ODE: exposure is supplied as ",
    "data, and the source obtained it from post hoc estimates of the ",
    "companion population PK model Biesdorf_2026_cofetuzumabPelidotin. ",
    "Conjugate exposure predicted response better than payload exposure."
  )
  reference <- paste(
    "Biesdorf C, Rinas M, Engelhardt B, Saab R, Guo C, Ferlini C, Freise KJ,",
    "Menon RM, Polepally AR. Population Pharmacokinetics and Exposure-Response",
    "Analyses to Support Cofetuzumab Pelidotin (ABBV-647) Dose Optimization in",
    "Non-Small Cell Lung Cancer Patients. Clin Pharmacol Ther. 2026.",
    "doi:10.1002/cpt.70432. PMCID PMC13473755.",
    "Coefficients from Supplementary Table S4; checked against the Table 3",
    "predicted ORR.",
    sep = " "
  )
  vignette <- "Biesdorf_2026_cofetuzumabPelidotin"

  units <- list(
    time = "n/a (static landmark exposure-response model; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "prob_orr (probability of a confirmed objective response, 0-1)"
  )

  covariateData <- list(
    CAV = list(
      description = "Individual average serum concentration of the cofetuzumab pelidotin CONJUGATE up to the time of best overall response (or the last scan for non-responders). Supplied as data: this model has no PK layer.",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Biesdorf 2026 Methods: 'The Cavg values were estimated from the cumulative exposure up to the time to occurrence of event (time of best overall response [BOR], last scan, or start of TEAE) divided by the time to the occurrence of respective efficacy or safety event', using empirical Bayes estimates from the conjugate popPK model (modellib('Biesdorf_2026_cofetuzumabPelidotin')). For dose projections Biesdorf 2026 Table 3 substitutes the steady-state Cavg (AUCtau / tau), equal to cycle-1 Cavg because the conjugate does not accumulate: median 10.1, 11.6, 13.5 and 15.4 ug/mL at 2.1, 2.4, 2.8 and 3.2 mg/kg Q3W in the phase 1b target population. Figure 1a quartile boundaries 5.74, 8.73, 12, 15.6 and 39.3 ug/mL. Enters LINEARLY and uncentred.",
      source_name = "Cofe-P (ADC) Cavg up to BOR"
    ),
    PTK7_HIGH_NONSQUAMOUS_EGFR_WT = list(
      description = "1 = patient in the trial's predefined target population, NSCLC with PTK7 expression in >= 90% of tumour cells at >= 2+ immunohistochemistry staining intensity, non-squamous histology and EGFR wild-type; 0 = any other NSCLC patient in the analysis.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-target NSCLC)",
      notes = "Biesdorf 2026 Results: 'Target population eligibility (PTK7-high, non-squamous, EGFR WT) was a statistically significant covariate on the exposure-efficacy slope'. Table S4: 28 target and 32 non-target patients. Acts additively on the logit-scale Cavg slope (0.35 + 0.03 = 0.38 per ug/mL); reading 0.03 as the target slope itself would give a target ORR near 0.2% at 13.5 ug/mL, contradicting Table 3's 20.6%. The individual components (PTK7 cut-off, histology, EGFR status) were tested separately and were not significant.",
      source_name = "Target population eligibility"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex; 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      notes = "Listed in Biesdorf 2026 Table S2 among the exposure-efficacy covariates (with age, albumin, soluble PTK7, PTK7 staining cut-off, tumour histology and EGFR status, each tested separately) and not retained: 'None of the other tested covariates ... evaluated in the model were statistically significant'. No point estimate is reported."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 60L,
    n_studies = 1L,
    n_observations = "60 binary confirmed-response records (one per patient; landmark analysis)",
    disease_state = "Recurrent PTK7-expressing NSCLC from the phase 1b study M19-611 (NCT04189614) with PK and efficacy data; 28 in the target population (PTK7 >= 90% of cells at >= 2+ intensity, non-squamous, EGFR wild-type) and 32 outside it. Observed ORR 12/60 (20%).",
    dose_range = "Cofe-P 2.8 mg/kg every 3 weeks, 60-minute intravenous infusion",
    regions = "M19-611 international sites",
    notes = paste0(
      "NSCLC patients from the first-in-human study were excluded ",
      "because their PTK7 assay was not comparable and histology and ",
      "EGFR status were unavailable. Only confirmed responses counted. ",
      "Logistic regression fitted with the R package gnm (R 4.2.2)."
    )
  )

  ini({
    # ==================================================================
    # Biesdorf 2026 Table S4, 'Logistic Regression Analysis of
    # Cofe-P(ADC) Exposure-Response', row header 'Cofetuzumab pelidotin
    # (ADC) Cavg [Linear Model]'. The Results text says 'Log model as
    # link function was selected based on AIC' for the efficacy
    # endpoints, but Table S4 labels the ORR model Linear and the linear
    # form reproduces the Table 3 predictions (6.47 / 10.97 / 20.4 /
    # 34.7% for the target population at Cavg 10.1 / 11.6 / 13.5 / 15.4
    # ug/mL with the slope at the top of its rounding interval, against
    # printed 6.52 / 11.1 / 20.6 / 35.1%), so the table wins.
    # ==================================================================
    logit_ref <- -6.56; label("Logit of the response probability at zero conjugate Cavg (unitless logit)") # Table S4 'Intercept' -6.56, RSE 28.4, p 0.0004
    e_cav_logit <- 0.35; label("Log-odds of response per ug/mL of conjugate Cavg outside the target population (per ug/mL)") # Table S4 'Slope Non-Target' 0.35, RSE 35.1, p 0.004
    e_ptk7_high_nonsquamous_egfr_wt_cav_logit <- 0.03; label("Additional log-odds of response per ug/mL of conjugate Cavg in the target population (per ug/mL)") # Table S4 'Slope Target' 0.03, RSE 156, p 0.521

    # ==================================================================
    # The source likelihood is Bernoulli: no random effect and no
    # residual error. rxode2 requires an observation declaration, so the
    # deterministic probability carries a tiny placeholder additive
    # residual, as in the other landmark logistic models of the library.
    # ==================================================================
    addSd_prob_orr <- fixed(0.001); label("Placeholder additive residual SD on the response probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_orr <- logit_ref + (e_cav_logit + e_ptk7_high_nonsquamous_egfr_wt_cav_logit * PTK7_HIGH_NONSQUAMOUS_EGFR_WT) * CAV
    prob_orr <- expit(logit_orr)

    # Downstream callers can sample binary outcomes with
    # rbinom(n, 1, prob_orr) on the rxSolve output.
    prob_orr ~ add(addSd_prob_orr)
  })
}
