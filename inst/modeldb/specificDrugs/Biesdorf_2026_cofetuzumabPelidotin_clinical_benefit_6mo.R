Biesdorf_2026_cofetuzumabPelidotin_clinical_benefit_6mo <- function() {
  description <- paste0(
    "Logistic-regression exposure-efficacy model for clinical benefit ",
    "lasting at least 6 months (CBR6: complete response, partial ",
    "response or stable disease for >= 6 months, RECIST 1.1) in adults ",
    "with PTK7-expressing recurrent non-small cell lung cancer treated ",
    "with cofetuzumab pelidotin (Cofe-P, ABBV-647) 2.8 mg/kg every 3 ",
    "weeks (Biesdorf 2026, n = 60, phase 1b study M19-611). The ",
    "probability is expit(a + b * CAV), linear in CAV, the individual ",
    "average Cofe-P CONJUGATE serum concentration up to the time of best ",
    "overall response, in ug/mL. Target-population membership was not ",
    "significant for this endpoint. There is no PK layer and no ODE: ",
    "exposure is supplied as data from the companion population PK model ",
    "Biesdorf_2026_cofetuzumabPelidotin. NOTE: Biesdorf 2026 does not ",
    "print the CBR6 coefficients; the two values here were back-solved ",
    "from the four model-predicted CBR6 probabilities of Table 3, which a ",
    "linear-in-Cavg logit reproduces to rounding; see the vignette."
  )
  reference <- paste(
    "Biesdorf C, Rinas M, Engelhardt B, Saab R, Guo C, Ferlini C, Freise KJ,",
    "Menon RM, Polepally AR. Population Pharmacokinetics and Exposure-Response",
    "Analyses to Support Cofetuzumab Pelidotin (ABBV-647) Dose Optimization in",
    "Non-Small Cell Lung Cancer Patients. Clin Pharmacol Ther. 2026.",
    "doi:10.1002/cpt.70432. PMCID PMC13473755.",
    "The coefficients are not tabulated by the source; they were back-solved",
    "from the Table 3 predicted CBR6 probabilities and checked against the",
    "fitted curve of Figure 1c.",
    sep = " "
  )
  vignette <- "Biesdorf_2026_cofetuzumabPelidotin"

  units <- list(
    time = "n/a (static landmark exposure-response model; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "prob_clinical_benefit_6mo (probability of clinical benefit lasting at least 6 months, 0-1)"
  )

  covariateData <- list(
    CAV = list(
      description = "Individual average serum concentration of the cofetuzumab pelidotin CONJUGATE up to the time of best overall response (or the last scan). Supplied as data: this model has no PK layer.",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Derived from empirical Bayes estimates of the conjugate popPK model (modellib('Biesdorf_2026_cofetuzumabPelidotin')) as cumulative exposure divided by time up to best overall response. Table 3 dose projections substitute the steady-state Cavg: median 10.1, 11.6, 13.5 and 15.4 ug/mL at 2.1, 2.4, 2.8 and 3.2 mg/kg Q3W. Figure 1d quartile boundaries 5.44, 9.76, 12.1, 16.3 and 39.3 ug/mL. Enters LINEARLY and uncentred.",
      source_name = "Cofe-P (ADC) Cavg up to BOR"
    )
  )

  covariatesDataExcluded <- list(
    TARGET_POPULATION = list(
      description = "1 = PTK7-high (>= 90% of cells at >= 2+), non-squamous, EGFR wild-type NSCLC; 0 = other NSCLC.",
      units = "(binary)",
      type = "binary",
      notes = "Tested on the CBR6 model and not retained: 'The influence of the target population was not significant on CBR or CBR6' (Biesdorf 2026 Results). The other Table S2 efficacy covariates (sex, age, albumin, soluble PTK7, PTK7 cut-off, histology, EGFR status) were likewise not significant. Used by the companion ORR model."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 60L,
    n_studies = 1L,
    n_observations = "60 binary CBR6 records (one per patient; landmark analysis)",
    disease_state = "Recurrent PTK7-expressing NSCLC from the phase 1b study M19-611 (NCT04189614) with PK and efficacy data. Observed CBR6 24/60 (40%).",
    dose_range = "Cofe-P 2.8 mg/kg every 3 weeks, 60-minute intravenous infusion",
    regions = "M19-611 international sites",
    notes = "Same 60-patient analysis set as the ORR model. Plain clinical benefit rate (without the 6-month requirement) had no significant exposure relationship and has no model. Logistic regression fitted with the R package gnm (R 4.2.2)."
  )

  ini({
    # ==================================================================
    # Biesdorf 2026 Results report a significant CBR6 exposure-response
    # relationship (P < 0.001, Figure 1c) but print no coefficients.
    # Table 3 gives the model-predicted CBR6 at four printed Cavg values:
    #   Cavg 10.1 / 11.6 / 13.5 / 15.4 ug/mL -> 24.3 / 34.2 / 48.8 / 63.6%
    # A least-squares line through logit(p) against Cavg fits these four
    # points with a residual SD of 0.0014 logit units - pure rounding of
    # the printed percentages and Cavg values - and reproduces them as
    # 24.32 / 34.17 / 48.79 / 63.62%. The log(Cavg) and log(Cavg + 1)
    # forms used for the other endpoints leave residuals 35-40 times
    # larger, so the CBR6 model is linear in Cavg.
    # Cross-check: the back-solved curve gives 23.7% at 10 ug/mL and 80%
    # at 18 ug/mL, on the fitted line of Figure 1c.
    # ==================================================================
    logit_ref <- -4.363; label("Logit of the CBR6 probability at zero conjugate Cavg (unitless logit)") # back-solved from Table 3 CBR6 predictions; no printed value exists
    e_cav_logit <- 0.3196; label("Log-odds of CBR6 per ug/mL of conjugate Cavg (per ug/mL)") # back-solved from Table 3 CBR6 predictions; no printed value exists

    addSd_prob_clinical_benefit_6mo <- fixed(0.001); label("Placeholder additive residual SD on the CBR6 probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_clinical_benefit_6mo <- logit_ref + e_cav_logit * CAV
    prob_clinical_benefit_6mo <- expit(logit_clinical_benefit_6mo)

    prob_clinical_benefit_6mo ~ add(addSd_prob_clinical_benefit_6mo)
  })
}
