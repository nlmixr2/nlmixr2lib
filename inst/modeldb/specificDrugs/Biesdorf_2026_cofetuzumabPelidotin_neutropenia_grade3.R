Biesdorf_2026_cofetuzumabPelidotin_neutropenia_grade3 <- function() {
  description <- paste0(
    "Logistic-regression exposure-safety model for grade >= 3 neutropenia (treatment-emergent adverse event) in ",
    "adults with advanced solid tumours treated with cofetuzumab ",
    "pelidotin (Cofe-P, ABBV-647) (Biesdorf 2026, n = 198, pooled ",
    "first-in-human study B7661001 and phase 1b NSCLC study M19-611). ",
    "The probability is expit(a + b * log(CAV + 1)), where CAV is the individual average ",
    "Cofe-P CONJUGATE serum concentration up to the onset of the event ",
    "(or the last available time for patients without it), in ug/mL. ",
    "No covariate was retained. There is no PK layer and no ODE: ",
    "exposure is supplied as data from the companion population PK model ",
    "Biesdorf_2026_cofetuzumabPelidotin. Conjugate exposure predicted ",
    "the event better than unconjugated Aur0101 payload exposure."
  )
  reference <- paste(
    "Biesdorf C, Rinas M, Engelhardt B, Saab R, Guo C, Ferlini C, Freise KJ,",
    "Menon RM, Polepally AR. Population Pharmacokinetics and Exposure-Response",
    "Analyses to Support Cofetuzumab Pelidotin (ABBV-647) Dose Optimization in",
    "Non-Small Cell Lung Cancer Patients. Clin Pharmacol Ther. 2026.",
    "doi:10.1002/cpt.70432. PMCID PMC13473755.",
    "Coefficients from Supplementary Table S5; checked against the Table 3",
    "predicted incidences.",
    sep = " "
  )
  vignette <- "Biesdorf_2026_cofetuzumabPelidotin"

  units <- list(
    time = "n/a (static landmark exposure-response model; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "prob_neutropenia_grade3 (probability of grade 3 or worse neutropenia, 0-1)"
  )

  covariateData <- list(
    CAV = list(
      description = "Individual average serum concentration of the cofetuzumab pelidotin CONJUGATE up to the onset of the safety event (or the last available time without it). Supplied as data: this model has no PK layer.",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Biesdorf 2026 Methods: Cavg is the cumulative exposure up to the start of the TEAE divided by the time to it, from empirical Bayes estimates of the conjugate popPK model (modellib('Biesdorf_2026_cofetuzumabPelidotin')); patients without the event use the last available time. Table 3 dose projections substitute the steady-state Cavg: median 10.1, 11.6, 13.5 and 15.4 ug/mL at 2.1, 2.4, 2.8 and 3.2 mg/kg Q3W. Figure 2 quartile boundaries 0.287, 7.61, 12.8, 17.3 and 42.6 ug/mL. Enters as the natural log of (CAV + 1), uncentred.",
      source_name = "Cofe-P (ADC) Cavg up to the event"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex; 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      notes = "Listed in Biesdorf 2026 Table S2 among the exposure-safety covariates (with age, weight, study) and not retained: 'None of the tested covariates ... evaluated in the models were statistically significant and therefore were not included'. No point estimate is reported."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 198L,
    n_studies = 2L,
    n_observations = "198 binary event records (one per patient)",
    disease_state = "Advanced solid tumours from the first-in-human study B7661001 (NCT02222922) and the phase 1b NSCLC study M19-611 (NCT04189614) with PK and safety data. Observed incidence 60/198 (30.3%).",
    dose_range = "Cofe-P 0.2-3.7 mg/kg every 3 weeks and 2.1-3.2 mg/kg every 2 weeks (B7661001); 2.8 mg/kg every 3 weeks (M19-611); 60-minute intravenous infusion.",
    notes = "Adverse events graded with NCI CTCAE v4.03 (B7661001) and v5.0 (M19-611). Logistic regression fitted with the R package gnm (R 4.2.2)."
  )

  ini({
    # ==================================================================
    # Biesdorf 2026 Table S5, 'Logistic Regression Analyses of Cofe-P
    # (ADC) Exposure-Safety', Neutropenia rows. The row header 'Cofe-P(ADC) Cavg [log(Cavg +1)]' fixes the exposure transform as log(Cavg + 1); the natural logarithm is the one that reproduces Table 3.
    # Table 3 cross-check (Cavg 10.1 / 11.6 / 13.5 / 15.4 ug/mL):
    # model 16.75 / 22.72 / 30.91 / 39.27% against printed 16.8 / 22.8 / 31.0 / 39.3%.
    # ==================================================================
    logit_ref <- -8.80; label("Logit of the grade 3 or worse neutropenia probability at zero conjugate Cavg, where log(Cavg + 1) = 0 (unitless logit)") # Table S5 Neutropenia 'Intercept' -8.80, RSE 16.7, p < 0.001
    e_cav_logit <- 2.99; label("Log-odds of grade 3 or worse neutropenia per unit of log(conjugate Cavg in ug/mL + 1) (unitless)") # Table S5 Neutropenia 'Slope' 2.99, RSE 17.6, p < 0.001

    # The source likelihood is Bernoulli; rxode2 needs an observation
    # declaration, so the probability carries a placeholder residual.
    addSd_prob_neutropenia_grade3 <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_neutropenia_grade3 <- logit_ref + e_cav_logit * log(CAV + 1)
    prob_neutropenia_grade3 <- expit(logit_neutropenia_grade3)

    prob_neutropenia_grade3 ~ add(addSd_prob_neutropenia_grade3)
  })
}
