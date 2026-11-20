Gaurav_2026_favezelimab_orr <- function() {
  description <- paste(
    "Logistic-regression exposure-efficacy model for objective response",
    "rate (ORR) in gastric cancer patients randomised to favezelimab 200 or",
    "700 mg IV every 3 weeks with pembrolizumab 200 mg in the phase 1 P001",
    "study (Gaurav 2026, n = 80, 9 responders). The probability of response",
    "is expit(a + b * (AUC - AUCc)), linear in the individual average",
    "favezelimab AUC per 21-day dosing interval (ng*h/mL) derived from the",
    "companion population PK model Gaurav_2026_favezelimab using each",
    "patient's dosing history. No covariate (baseline tumor size, sex, age,",
    "race) was retained. The intercept and slope are printed in Table S7;",
    "the centring value AUCc is not printed and was recovered by digitising",
    "the fitted curve of Figure 3. The Grade >= 3 drug-related AE and",
    "drug-related AEOSI analyses found no exposure-response relationship",
    "and were not modelled further by the source."
  )
  reference <- paste(
    "Gaurav M, Barcomb H, Maxwell KF, Kandala B, Chatterjee MS.",
    "Population Pharmacokinetic, Exposure-Response Efficacy and Safety",
    "Analyses of Favezelimab in Patients With Solid Tumors.",
    "Clin Transl Sci. 2026;19(7):e70658. doi:10.1111/cts.70658.",
    "PMCID PMC13330136. Intercept and slope from Supplementary Table S7;",
    "the AUC centring value was digitised from the fitted line of Figure 3."
  )
  vignette <- "Gaurav_2026_favezelimab"

  units <- list(
    time = "n/a (static landmark exposure-response model; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the AUC_FAVEZELIMAB covariate column)",
    concentration = "prob_orr (probability of objective response, 0-1; also logit_orr)"
  )

  covariateData <- list(
    AUC_FAVEZELIMAB = list(
      description = "Individual average favezelimab serum AUC per 21-day dosing interval over the treatment period, computed from the population PK model with the patient's actual dosing history",
      units = "ng*h/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Gaurav 2026 Methods: 'average area under the concentration-time curve [AUCavg] considering individual-specific dosing history', derived from the final popPK model (modellib('Gaurav_2026_favezelimab')). Enters LINEARLY and centred on 1.54e7 ng*h/mL (digitised). Figure 3 shows the analysis-set range spanning about 0.6e7 to 7.1e7 ng*h/mL; a typical patient at steady state has about 0.87e7 at 200 mg and 3.0e7 at 700 mg Q3W (dose / CL). Table S7 labels the term 'Centered Favezelimab Average AUC (ng h/mL)'.",
      source_name = "AUCavg"
    )
  )

  covariatesDataExcluded <- list(
    TUMSZ = list(
      description = "Baseline tumor size (sum of target-lesion diameters)",
      units = "mm",
      type = "continuous",
      notes = "Tested by forward selection (p < 0.01) and backward elimination (p < 0.001) and not retained: 'None of the evaluated covariates (baseline tumor size [mm], sex, age [years], and race) were found to be statistically significant predictors of ORR'. No estimate is reported. The same applies to sex, age and race."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 80L,
    n_studies = 1L,
    n_observations = "80 binary objective-response records (one per patient)",
    disease_state = "Advanced gastric cancer (P001 randomised dose-finding cohort)",
    dose_range = "Favezelimab 200 mg (n = 40) or 700 mg (n = 40) IV every 3 weeks, each with pembrolizumab 200 mg Q3W",
    regions = "Multinational (MK-4280-001 / NCT02720068)",
    notes = "Observed ORR 7.5% (3/40) at 200 mg and 15.0% (6/40) at 700 mg; 11.3% (9/80) overall (Table S6). The response assessor is not stated. The phase 3 KEYFORM-007 colorectal data were deliberately excluded from the exposure-response analysis because only one dose was studied there."
  )

  ini({
    # Supplementary Table S7 (glm, binomial, R 4.3.2):
    #   (Intercept)                               -2.48928  SE 0.46387
    #   Centered Favezelimab Average AUC (ng h/mL) 4.36E-08 SE 2.07E-08, p 0.03511
    # The centring constant AUCc is not printed. It was recovered from
    # Figure 3 by reading the fitted line pixel by pixel (about 470
    # columns) and solving logit(p) = -2.48928 + 4.36e-8 * (AUC - AUCc):
    # median AUCc 1.54e7 ng*h/mL (interquartile 1.525e7 to 1.555e7). A
    # free two-parameter fit of the digitised curve returns a slope of
    # 4.373e-8, within 0.3% of the printed 4.36e-8, which confirms the
    # axis calibration.
    logit_ref <- -2.48928; label("Logit of the response probability at the centring AUC (unitless logit)") # Table S7 (Intercept) -2.48928
    e_auc_favezelimab_logit <- 4.36e-8; label("Log-odds of objective response per 1 ng*h/mL increase in average favezelimab AUC (unitless logit)") # Table S7 'Centered Favezelimab Average AUC (ng h/mL)' 4.36E-08
    auc_favezelimab_ref <- fixed(1.54e7); label("Centring value of the average favezelimab AUC (ng*h/mL)") # digitised from Gaurav 2026 Figure 3; no printed value exists

    # The source fits a binomial logistic regression with an exact
    # Bernoulli likelihood: no random effect and no residual error. This
    # placeholder residual exists only so rxode2 accepts an observation.
    addSd_prob_orr <- fixed(0.001); label("Placeholder additive residual SD on the response probability; the source likelihood is Bernoulli") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_orr <- logit_ref +
      e_auc_favezelimab_logit * (AUC_FAVEZELIMAB - auc_favezelimab_ref)
    prob_orr <- expit(logit_orr)

    prob_orr ~ add(addSd_prob_orr)
  })
}
