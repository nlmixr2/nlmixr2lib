Sinnollareddy_2026_epcoritamab_maintained_q2w <- function() {
  description <- paste(
    "Logistic exposure-efficacy model for the probability of maintained or",
    "improved response during the Q2W dosing period (cycles 4-9) among",
    "patients who were responders at entry to or during that period versus",
    "steady-state trough plasma epcoritamab concentration during the Q2W",
    "dosing period (cycles 4-9) (CTROUGH, ug/mL) in responders with relapsed",
    "or refractory follicular lymphoma in the expansion cohorts of EPCORE",
    "NHL-1 and EPCORE NHL-3 who entered or responded during the Q2W 48 mg",
    "dosing period (cycles 4-9) (Sinnollareddy 2026, supplementary Figure S4",
    "(top left); 86 of 100 with event, p = 0.265). The probability is",
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
    concentration = "prob_response_maintained (probability of maintained or improved response during the Q2W dosing period, 0-1)"
  )

  covariateData <- list(
    CTROUGH = list(
      description = "Steady-state trough plasma epcoritamab concentration during the Q2W dosing period (cycles 4-9)",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Individual model-predicted trough at steady state under 48 mg every 2",
        "weeks, from the updated popPK model",
        "(modellib('Sinnollareddy_2026_epcoritamab')). Free (QSS) concentration,",
        "enters linearly and uncentred. Observed range about 0.7-9.9 ug/mL",
        "(supplementary Figure S4)."
      ),
      source_name = "Q2W Ctr"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 100,
    n_studies = 2,
    disease_state = "Responders with relapsed or refractory follicular lymphoma in the expansion cohorts of EPCORE NHL-1 and EPCORE NHL-3 who entered or responded during the Q2W 48 mg dosing period (cycles 4-9)",
    n_events = 86,
    event_pct = 86,
    notes = paste(
      "Supplementary Figure S4, all patients (N = 100). EPCORE NHL-1",
      "NCT03625037; EPCORE NHL-3 NCT04542824."
    )
  )

  ini({
    # supplementary Figure S4 (top left) fitted logistic curve, digitised by the maintainers
    # (coefficients not printed): logit(p) = intercept + slope * CTROUGH.
    # Rounded to four significant figures.
    logit_ref <- 1.018; label("Logit of the event probability at zero Q2W Ctr (unitless logit)") # supplementary Figure S4 (top left) fitted curve, digitised; intercept not printed
    e_ctrough_logit <- 0.1814; label("Log-odds of the event per 1 ug/mL increase in Q2W Ctr (unitless logit per ug/mL)") # supplementary Figure S4 (top left) fitted curve, digitised; slope not printed (p = 0.265 printed in the panel)

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
