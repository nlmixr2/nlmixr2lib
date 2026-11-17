Sinnollareddy_2026_epcoritamab_3step_crs_auc <- function() {
  description <- paste(
    "Logistic exposure-safety model for the probability of cytokine release",
    "syndrome of any grade (ASTCT grading) over Cycle 1 and the Cycle 2 Day",
    "1 dosing interval versus model-predicted plasma epcoritamab AUC over",
    "Cycle 1 of the 3-step regimen, including the Cycle 2 Day 1 dose",
    "interval (AUC_EPCORITAMAB, day*ug/mL) in adults with relapsed or",
    "refractory follicular lymphoma in the EPCORE NHL-1 optimization cohort",
    "treated with the 3-step step-up regimen (0.16/0.8/3 mg then 48 mg; Arm",
    "B 0.16/0.8/6/48 mg in 6 patients) with hydration and dexamethasone",
    "prophylaxis in Cycle 1 (Sinnollareddy 2026, supplementary Figure S8",
    "(top); 40 of 81 with event, p = 0.124). The probability is",
    "expit(intercept + slope * AUC_EPCORITAMAB), a linear logistic",
    "regression. The coefficients are not printed in the paper; they were",
    "digitised by the maintainers from the fitted regression curve of the",
    "cited figure. There is no PK layer and no ODE; AUC_EPCORITAMAB is",
    "supplied as a data column, in the source an individual prediction from",
    "the updated popPK model (modellib('Sinnollareddy_2026_epcoritamab'))."
  )
  reference <- paste(
    "Sinnollareddy M, Sanghavi K, Gibiansky L, Li T, Putnins M, Mohamed MEF,",
    "Patah P, Favaro E, Gupta M, Parikh A, Xu S. Exposure-Response Analyses to",
    "Inform the Optimal Epcoritamab Monotherapy Dosing Regimen in Relapsed or",
    "Refractory Follicular Lymphoma. Clin Pharmacokinet. 2026;65(9):1415-1428.",
    "doi:10.1007/s40262-026-01661-1."
  )
  vignette <- "Sinnollareddy_2026_epcoritamab"
  units <- list(
    time = "n/a (static landmark regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the AUC_EPCORITAMAB data column)",
    concentration = "prob_crs (probability of cytokine release syndrome of any grade, 0-1)"
  )

  covariateData <- list(
    AUC_EPCORITAMAB = list(
      description = "Model-predicted plasma epcoritamab AUC over Cycle 1 of the 3-step regimen, including the Cycle 2 Day 1 dose interval",
      units = "day*ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Individual model-predicted AUC over the priming dose, the two",
        "intermediate doses and the first two 48 mg full doses (the Cycle 2 Day",
        "1 dose interval included; supplementary Figure S8 legend), from the",
        "updated popPK model (modellib('Sinnollareddy_2026_epcoritamab')).",
        "Enters linearly and uncentred. Observed range about 0.7-148 day*ug/mL."
      ),
      source_name = "Cycle 1 AUC"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 81,
    n_studies = 1,
    disease_state = "Adults with relapsed or refractory follicular lymphoma in the EPCORE NHL-1 optimization cohort treated with the 3-step step-up regimen (0.16/0.8/3 mg then 48 mg; Arm B 0.16/0.8/6/48 mg in 6 patients) with hydration and dexamethasone prophylaxis in Cycle 1",
    n_events = 40,
    event_pct = 49.4,
    notes = paste(
      "EPCORE NHL-1 optimization cohort (supplementary Table S1, N = 81).",
      "EPCORE NHL-1 NCT03625037; EPCORE NHL-3 NCT04542824."
    )
  )

  ini({
    # supplementary Figure S8 (top) fitted logistic curve, digitised by the maintainers
    # (coefficients not printed): logit(p) = intercept + slope * AUC_EPCORITAMAB.
    # Rounded to four significant figures.
    logit_ref <- 0.7079; label("Logit of the event probability at zero Cycle 1 AUC (unitless logit)") # supplementary Figure S8 (top) fitted curve, digitised; intercept not printed
    e_auc_logit <- -0.01641; label("Log-odds of the event per 1 day*ug/mL increase in Cycle 1 AUC (unitless logit per day*ug/mL)") # supplementary Figure S8 (top) fitted curve, digitised; slope not printed (p = 0.124 printed in the panel)

    # Bernoulli likelihood in the source; the fixed additive residual is a
    # placeholder so rxode2 has an error model, not a published quantity.
    addSd_prob_crs <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_p <- logit_ref + e_auc_logit * AUC_EPCORITAMAB
    prob_crs <- expit(logit_p)

    prob_crs ~ add(addSd_prob_crs)
  })
}
