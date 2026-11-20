Sinnollareddy_2026_epcoritamab_maintained_q4w <- function() {
  description <- paste(
    "Logistic exposure-efficacy model for the probability of maintained or",
    "improved response during the Q4W dosing period (cycles 10+) among",
    "patients who were responders at entry to or during that period versus",
    "steady-state trough plasma epcoritamab concentration during the Q4W",
    "dosing period (cycles 10+) (CTROUGH, ug/mL) in responders with relapsed",
    "or refractory follicular lymphoma in the expansion cohorts of EPCORE",
    "NHL-1 and EPCORE NHL-3 who entered or responded during the Q4W 48 mg",
    "dosing period (cycles 10+) (Sinnollareddy 2026, supplementary Figure S4",
    "(bottom left); 59 of 61 with event, p = 0.15). The probability is",
    "expit(intercept + slope * CTROUGH), a linear logistic regression. The",
    "coefficients are not printed in the paper; they were digitised by the",
    "maintainers from the fitted regression curve of the cited figure. There",
    "is no PK layer and no ODE; CTROUGH is supplied as a data column, in the",
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
    dosing = "n/a (no dose events; exposure enters as the CTROUGH data column)",
    concentration = "prob_response_maintained (probability of maintained or improved response during the Q4W dosing period, 0-1)"
  )

  covariateData <- list(
    CTROUGH = list(
      description = "Steady-state trough plasma epcoritamab concentration during the Q4W dosing period (cycles 10+)",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Individual model-predicted trough at steady state under 48 mg every 4",
        "weeks, from the updated popPK model",
        "(modellib('Sinnollareddy_2026_epcoritamab')). Free (QSS) concentration,",
        "enters linearly and uncentred. Observed range about 0.4-3.8 ug/mL",
        "(supplementary Figure S4)."
      ),
      source_name = "Q4W Ctr"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 61,
    n_studies = 2,
    disease_state = "Responders with relapsed or refractory follicular lymphoma in the expansion cohorts of EPCORE NHL-1 and EPCORE NHL-3 who entered or responded during the Q4W 48 mg dosing period (cycles 10+)",
    n_events = 59,
    event_pct = 96.7,
    notes = paste(
      "Supplementary Figure S4, all patients (N = 61). EPCORE NHL-1",
      "NCT03625037; EPCORE NHL-3 NCT04542824."
    )
  )

  ini({
    # supplementary Figure S4 (bottom left) fitted logistic curve, digitised by the maintainers
    # (coefficients not printed): logit(p) = intercept + slope * CTROUGH.
    # Rounded to four significant figures.
    logit_ref <- 5.896; label("Logit of the event probability at zero Q4W Ctr (unitless logit)") # supplementary Figure S4 (bottom left) fitted curve, digitised; intercept not printed
    e_ctrough_logit <- -1.17; label("Log-odds of the event per 1 ug/mL increase in Q4W Ctr (unitless logit per ug/mL)") # supplementary Figure S4 (bottom left) fitted curve, digitised; slope not printed (p = 0.15 printed in the panel)

    # Bernoulli likelihood in the source; the fixed additive residual is a
    # placeholder so rxode2 has an error model, not a published quantity.
    addSd_prob_response_maintained <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_p <- logit_ref + e_ctrough_logit * CTROUGH
    prob_response_maintained <- expit(logit_p)

    prob_response_maintained ~ add(addSd_prob_response_maintained)
  })
}
