Li_2025_mosunetuzumab_orr <- function() {
  description <- paste(
    "Emax logistic-regression exposure-efficacy model for INVESTIGATOR-",
    "ASSESSED OBJECTIVE RESPONSE (complete or partial response) versus",
    "the model-predicted two-cycle",
    "mosunetuzumab AUC0-42 (AUC_MOSUNETUZUMAB, day*ug/mL) in adults with",
    "relapsed/refractory follicular lymphoma after >= 2 prior therapies",
    "treated with intravenous mosunetuzumab monotherapy (Li 2025, study",
    "GO29781, n = 159; Table S7). The probability is expit(beta0 + Emax *",
    "AUC / (E50 + AUC)). This is the base model (no covariates) that the",
    "paper uses for all dose-regimen predictions (Figure 4). There is no PK",
    "layer and no ODE; AUC0-42 is supplied as a data column, in the source",
    "an empirical-Bayes prediction from the Bender 2024 popPK model",
    "(packaged as Bender_2024_mosunetuzumab, whose `auc` state gives it",
    "directly)."
  )
  reference <- paste(
    "Li C-C, Bender B, Wilkins J, Li F, Turner DC, Wang B, Deng R,",
    "Vadhavkar S, Li Z, Kwan A, Huang H, Peng K, Penuel E, Huw L-Y, Chanu P,",
    "Li C, Yin S, Wei MC. A Novel Step-Up Dosage Regimen for Enhancing the",
    "Benefit-to-Risk Ratio of Mosunetuzumab in Relapsed or Refractory",
    "Follicular Lymphoma. Clin Pharmacol Ther. 2025;117(2):465-474.",
    "doi:10.1002/cpt.3445."
  )
  vignette <- "Li_2025_mosunetuzumab"
  units <- list(
    time = "n/a (static two-cycle landmark regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the AUC_MOSUNETUZUMAB data column)",
    concentration = "prob_orr_investigator (probability of objective response, CR or PR, investigator-assessed, 0-1)"
  )

  covariateData <- list(
    AUC_MOSUNETUZUMAB = list(
      description = "Model-predicted mosunetuzumab serum AUC from day 0 to day 42 (the first two 21-day cycles)",
      units = "day*ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Li 2025 Methods: 'model-estimated two-cycle area-under-the-",
        "concentration curve (AUC0-42)', the primary efficacy exposure",
        "metric, selected over three-cycle, four-cycle and steady-state AUC",
        "by AIC. Derived from the individual empirical-Bayes parameters of",
        "the Bender 2024 popPK model (Bender_2024_mosunetuzumab; its `auc`",
        "state read at day 42 is this quantity). Observed range in the",
        "efficacy population 0.4-251 day*ug/mL; median 145 in cohort B11",
        "(1/2/60/30 mg) and 85.9 overall (Table S10). Enters through an",
        "Emax term, untransformed."
      ),
      source_name = "AUC0-42"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 159,
    n_studies = 1,
    age_range = "27-90 years (median 60)",
    weight_range = "42-140 kg (median 80.4)",
    sex_female_pct = 38.4,
    race_ethnicity = c(
      White = 76.1,
      Asian = 15.1,
      Black = 5.03,
      `American Indian or Alaska Native` = 0.629,
      Unknown = 3.14
    ),
    disease_state = "Relapsed/refractory follicular lymphoma after >= 2 prior systemic therapies (median 3)",
    dose_range = "0.2-2.8 mg IV q3w fixed dosing (Group A, n = 7) and Cycle 1 step-up dosing from 0.4/1.0/2.8 to 1/2/60/30 mg IV q3w (Group B, n = 152; cohort B11 1/2/60/30 mg n = 90)",
    regions = "Multinational",
    notes = paste(
      "Efficacy exposure-response population of the phase I/II GO29781",
      "study (NCT02500407), clinical cutoff 27 August 2021; baseline",
      "characteristics from Li 2025 Table S10. Investigator-assessed",
      "(International Harmonization Project criteria) response; ORR 77.8",
      "percent in B11 and 63.0 percent in B7 (1/2/13.5 mg)."
    )
  )

  ini({
    logite0 <- -4.76; label("Logit of the probability of objective response at AUC0-42 = 0, beta0 (unitless logit)") # Li 2025 Table S7, OR row, beta0 = -4.76 (%RSE 124)
    emax <- 6.23; label("Maximal logit shift of the objective-response probability produced by AUC0-42, Emax (unitless logit)") # Li 2025 Table S7, OR row, Emax = 6.23 (%RSE 92.3)
    lec50 <- log(3.12); label("AUC0-42 producing half of Emax, E50 (day*ug/mL)") # Li 2025 Table S7, OR row, E50 = 3.12 (%RSE 156)

    # No between-subject variability and no residual error: the source is a
    # Bernoulli-likelihood logistic regression. The placeholder additive term
    # exists only because rxode2 requires an observation declaration.
    addSd_prob_orr_investigator <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # Li 2025 Methods, Emax logistic regression:
    #   logit(p) = beta0 + Emax * AUC0-42 / (E50 + AUC0-42)
    ec50 <- exp(lec50)
    logit_orr <- logite0 + emax * AUC_MOSUNETUZUMAB / (ec50 + AUC_MOSUNETUZUMAB)

    prob_orr_investigator <- expit(logit_orr)

    # Deterministic probability of the event. Downstream callers can sample
    # binary outcomes with rbinom(n, 1, prob_orr_investigator) on the rxSolve
    # output.
    prob_orr_investigator ~ add(addSd_prob_orr_investigator)
  })
}
