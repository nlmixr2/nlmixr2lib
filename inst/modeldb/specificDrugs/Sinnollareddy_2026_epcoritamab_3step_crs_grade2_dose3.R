Sinnollareddy_2026_epcoritamab_3step_crs_grade2_dose3 <- function() {
  description <- paste(
    "Logistic exposure-safety model for the probability of cytokine release",
    "syndrome of grade >= 2 within the dosing interval of dose 3 (3 mg",
    "second intermediate dose; Cycle 1 Day 15) versus maximum plasma",
    "epcoritamab concentration within the dosing interval of dose 3 (3 mg,",
    "Cycle 1 Day 15) of the 3-step regimen (CMAX, ug/mL) in adults with",
    "relapsed or refractory follicular lymphoma in the EPCORE NHL-1",
    "optimization cohort treated with the 3-step step-up regimen (0.16/0.8/3",
    "mg then 48 mg; Arm B 0.16/0.8/6/48 mg in 6 patients) with hydration and",
    "dexamethasone prophylaxis in Cycle 1 (Sinnollareddy 2026, Figure 5B",
    "(Dose 3); 1 of 81 with event, p = 0.599). The probability is",
    "expit(intercept + slope * CMAX), a linear logistic regression. The",
    "coefficients are not printed in the paper; they were digitised by the",
    "maintainers from the fitted regression curve of the cited figure. There",
    "is no PK layer and no ODE; CMAX is supplied as a data column, in the",
    "source an individual prediction from the updated popPK model",
    "(modellib('Sinnollareddy_2026_epcoritamab'))."
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
    dosing = "n/a (no dose events; exposure enters as the CMAX data column)",
    concentration = "prob_crs_stepup3_grade2 (probability of cytokine release syndrome of grade >= 2 within the dosing interval of dose 3, 0-1)"
  )

  covariateData <- list(
    CMAX = list(
      description = "Maximum plasma epcoritamab concentration within the dosing interval of dose 3 (3 mg, Cycle 1 Day 15) of the 3-step regimen",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Individual model-predicted Cmax of the Cycle 1 Day 15 second",
        "intermediate-dose interval, from the updated popPK model",
        "(modellib('Sinnollareddy_2026_epcoritamab')). Free (QSS) concentration,",
        "enters linearly and uncentred. Observed range about 0-0.305 ug/mL."
      ),
      source_name = "Dose 3 Cmax"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 81,
    n_studies = 1,
    disease_state = "Adults with relapsed or refractory follicular lymphoma in the EPCORE NHL-1 optimization cohort treated with the 3-step step-up regimen (0.16/0.8/3 mg then 48 mg; Arm B 0.16/0.8/6/48 mg in 6 patients) with hydration and dexamethasone prophylaxis in Cycle 1",
    n_events = 1,
    event_pct = 1.2,
    notes = paste(
      "EPCORE NHL-1 optimization cohort (supplementary Table S1, N = 81).",
      "EPCORE NHL-1 NCT03625037; EPCORE NHL-3 NCT04542824."
    )
  )

  ini({
    # Figure 5B (Dose 3) fitted logistic curve, digitised by the maintainers
    # (coefficients not printed): logit(p) = intercept + slope * CMAX.
    # Rounded to four significant figures.
    logit_ref <- -3.661; label("Logit of the event probability at zero Dose 3 Cmax (unitless logit)") # Figure 5B (Dose 3) fitted curve, digitised; intercept not printed
    e_cmax_logit <- -8.748; label("Log-odds of the event per 1 ug/mL increase in Dose 3 Cmax (unitless logit per ug/mL)") # Figure 5B (Dose 3) fitted curve, digitised; slope not printed (p = 0.599 printed in the panel)

    # Bernoulli likelihood in the source; the fixed additive residual is a
    # placeholder so rxode2 has an error model, not a published quantity.
    addSd_prob_crs_stepup3_grade2 <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_p <- logit_ref + e_cmax_logit * CMAX
    prob_crs_stepup3_grade2 <- expit(logit_p)

    prob_crs_stepup3_grade2 ~ add(addSd_prob_crs_stepup3_grade2)
  })
}
