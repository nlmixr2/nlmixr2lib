Quintairos_2021_kidneyTransplantRejection <- function() {
  description <- paste(
    "Logistic regression of the risk of biopsy-proven acute rejection (AR) on",
    "urinary-pellet miR155-5p relative expression in adult de novo kidney",
    "transplant recipients during the first 6 months post-transplant",
    "(Quintairos 2021). logit(P) = -5.89 + 3.51 * MIR155_URINE, where",
    "MIR155_URINE is the 2^-dCq relative expression measured at a study visit and",
    "P is the probability that AR is diagnosed before the next visit (rejection",
    "events were attributed to the visit preceding their occurrence). No",
    "tacrolimus or mycophenolic acid exposure term is present: individual",
    "cumulative AUC and mean trough concentrations from the companion popPK",
    "models, and urinary CXCL-10, were tested and not retained. No",
    "between-subject variability (the control stream fixes the only omega at 0)",
    "and no residual error (Bernoulli likelihood)."
  )
  reference <- paste(
    "Quintairos L, Colom H, Millan O, Fortuna V, Espinosa C, Guirado L, Budde K,",
    "Sommerer C, Lizana A, Lopez-Pua Y, Brunet M. Early prognostic performance of",
    "miR155-5p monitoring for the risk of rejection: Logistic regression with a",
    "population pharmacokinetic approach in adult kidney transplant patients.",
    "PLoS ONE. 2021;16(1):e0245880. doi:10.1371/journal.pone.0245880.",
    "Parameter values from Table 5; model structure from the S1 Appendix NONMEM",
    "control stream ('miR155-5p Logistic regression').",
    sep = " "
  )
  vignette <- "Quintairos_2021_tacrolimus_mycophenolic_acid"
  units <- list(
    time = "n/a (static per-visit landmark regression; time only orders the visits)",
    dosing = "n/a (no dose events and no drug exposure covariate)",
    concentration = "prob_acute_rejection (probability of acute rejection before the next visit, 0-1; also logit_acute_rejection)"
  )

  covariateData <- list(
    MIR155_URINE = list(
      description = "Urinary-pellet miR155-5p relative expression measured by qPCR at the visit (2^-dCq against the reference control).",
      units = "(relative expression, 2^-dCq)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying: measured in first-morning urine at week 1 and months 1, 2,",
        "3 and 6. Enters the logit linearly and untransformed, so the intercept is",
        "the logit at zero expression. Table 2 labels the unit 'dCt', but the",
        "Methods define the value as 2^-dCq (higher value = higher expression) and",
        "the deposited S3 Table values are positive, 0.01-2.80 (median 0.08, mean",
        "0.32 over the 183 modelled visits), which is the 2^-dCq scale. Visits",
        "with a missing (9999) or zero value were excluded from the fit (control",
        "stream IGNORE records)."
      ),
      source_name = "M155"
    )
  )

  covariatesDataExcluded <- list(
    AUC = list(
      description = "Individual predicted cumulative tacrolimus or mycophenolic acid AUC from the companion popPK models.",
      units = "ng*h/mL (tacrolimus); mg*h/L (MPA)",
      type = "continuous",
      notes = paste(
        "Tested as an explanatory variable and not retained (Results, logistic",
        "regression model). Integrated in a dummy compartment in the source; no",
        "point estimate exists because the term did not enter the final model."
      )
    ),
    CTROUGH = list(
      description = "Individual predicted mean tacrolimus or mycophenolic acid trough concentration.",
      units = "ng/mL (tacrolimus); mg/L (MPA)",
      type = "continuous",
      notes = "Tested and not retained (Results, logistic regression model)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 58L,
    n_studies = 1L,
    n_observations = "183 visit-level binary AR records (8 AR events) after the control stream excluded visits with missing or zero miR155-5p",
    age_range = "Median 48 years (IQR 38-58); recipients older than 70 years were excluded",
    weight_range = "Median 73 kg (IQR 62.9-86.8)",
    sex_female_pct = 34.5,
    race_ethnicity = "Caucasian only (Discussion)",
    disease_state = paste(
      "Adult de novo kidney transplant recipients (living donor 30, deceased donor",
      "28) on tacrolimus, mycophenolate mofetil, methylprednisolone and basiliximab",
      "induction. Eight of 58 (14%) had biopsy-proven (Banff 2011) cellular acute",
      "rejection: four in week 1, three at the end of month 1 and one in month 6."
    ),
    dose_range = "No drug input in this model; tacrolimus and MMF doses were adjusted by therapeutic drug monitoring",
    regions = "Germany (Charite Berlin, Heidelberg) and Spain (Fundacio Puigvert, Barcelona); EudraCT 2013-001817-33",
    notes = "Visits: week 1 (day 1-11), months 1, 2, 3 and 6. Urinary miR155-5p global mean 0.39 (IQR 0.03-0.55) (Table 2)."
  )

  ini({
    # Quintairos 2021 Table 5; control stream '$PRED ... A1 = B0 + B1 + ETA(1)'
    # with B1 = THETA(2) * M155 and '$OMEGA 0 FIX', METHOD=1 LIKELIHOOD LAPLACE.
    # The stream $THETA records (-7, 4) are initial estimates only.
    logit_ref <- -5.89
    label("Logit of the probability of acute rejection before the next visit at zero urinary miR155-5p expression (unitless logit)") # Table 5, 'beta0' -5.89 (RSE 15%); bootstrap mean -5.94
    e_mir155_urine_logit <- 3.51
    label("Change in the logit of acute rejection per unit of urinary miR155-5p relative expression (unitless logit per 2^-dCq unit)") # Table 5, 'beta1' 3.51 (RSE 24%); bootstrap mean 3.45

    # No residual error in the source (Bernoulli likelihood). This placeholder
    # only gives rxode2 an error model to attach to the typical-value probability.
    addSd_prob_acute_rejection <- fixed(0.001)
    label("Placeholder additive residual SD on the typical-value acute-rejection probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # Quintairos 2021 Table 5: LOGIT(Pi) = beta0 + beta1 * miR155-5p
    logit_acute_rejection <- logit_ref + e_mir155_urine_logit * MIR155_URINE
    prob_acute_rejection <- expit(logit_acute_rejection)

    prob_acute_rejection ~ add(addSd_prob_acute_rejection)
  })
}
