Sinnollareddy_2026_epcoritamab_fl_crs_tocilizumab <- function() {
  description <- paste(
    "Logistic exposure-safety model for the probability of cytokine release",
    "syndrome treated with tocilizumab versus average plasma epcoritamab",
    "concentration over Cycles 1-3 (the QW dosing period) (CAV, ug/mL) in",
    "adults with relapsed or refractory follicular lymphoma (grade 1-3A)",
    "from the escalation and expansion cohorts of EPCORE NHL-1 and EPCORE",
    "NHL-3 treated with subcutaneous epcoritamab on the 2-step step-up",
    "regimen (0.16/0.8 mg then full doses of 0.12-48 mg; 149 of 160 on 48",
    "mg; QW cycles 1-3, Q2W cycles 4-9, Q4W from cycle 10) (Sinnollareddy",
    "2026, supplementary Figure S6 (CRS with tocilizumab); 39 of 160 with",
    "event, p = 0.601). The probability is expit(intercept + slope * CAV), a",
    "linear logistic regression. The coefficients are not printed in the",
    "paper; they were digitised by the maintainers from the fitted",
    "regression curve of the cited figure. There is no PK layer and no ODE;",
    "CAV is supplied as a data column, in the source an individual",
    "prediction from the updated popPK model",
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
    dosing = "n/a (no dose events; exposure enters as the CAV data column)",
    concentration = "prob_crs_tocilizumab (probability of cytokine release syndrome treated with tocilizumab, 0-1)"
  )

  covariateData <- list(
    CAV = list(
      description = "Average plasma epcoritamab concentration over Cycles 1-3 (the QW dosing period)",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Individual model-predicted average concentration from the start of",
        "Cycle 1 to the end of Cycle 3 (first Cycle 4 dose, or last Cycle 3 dose",
        "plus 7 days), computed from each patient's actual dosing history and",
        "the individual parameters of the updated popPK model",
        "(modellib('Sinnollareddy_2026_epcoritamab')). Free (QSS) concentration,",
        "enters linearly and uncentred. Observed range about 0-14.5 ug/mL."
      ),
      source_name = "Cycle 1-3 Cavg"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 160,
    n_studies = 2,
    disease_state = "Adults with relapsed or refractory follicular lymphoma (grade 1-3A) from the escalation and expansion cohorts of EPCORE NHL-1 and EPCORE NHL-3 treated with subcutaneous epcoritamab on the 2-step step-up regimen (0.16/0.8 mg then full doses of 0.12-48 mg; 149 of 160 on 48 mg; QW cycles 1-3, Q2W cycles 4-9, Q4W from cycle 10)",
    n_events = 39,
    event_pct = 24.4,
    notes = paste(
      "Exposure-efficacy and R/R FL exposure-safety population (main text",
      "Results; supplementary Table S5 quartile counts). EPCORE NHL-1",
      "NCT03625037; EPCORE NHL-3 NCT04542824."
    )
  )

  ini({
    # supplementary Figure S6 (CRS with tocilizumab) fitted logistic curve, digitised by the maintainers
    # (coefficients not printed): logit(p) = intercept + slope * CAV.
    # Rounded to four significant figures.
    logit_ref <- -0.9606; label("Logit of the event probability at zero Cycle 1-3 Cavg (unitless logit)") # supplementary Figure S6 (CRS with tocilizumab) fitted curve, digitised; intercept not printed
    e_cav_logit <- -0.03467; label("Log-odds of the event per 1 ug/mL increase in Cycle 1-3 Cavg (unitless logit per ug/mL)") # supplementary Figure S6 (CRS with tocilizumab) fitted curve, digitised; slope not printed (p = 0.601 printed in the panel)

    # Bernoulli likelihood in the source; the fixed additive residual is a
    # placeholder so rxode2 has an error model, not a published quantity.
    addSd_prob_crs_tocilizumab <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_p <- logit_ref + e_cav_logit * CAV
    prob_crs_tocilizumab <- expit(logit_p)

    prob_crs_tocilizumab ~ add(addSd_prob_crs_tocilizumab)
  })
}
