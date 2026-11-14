Tachet_2026_ruxolitinib_ae_auc24 <- function() {
  description <- paste0(
    "Univariable logistic-regression exposure-toxicity model for the ",
    "occurrence of any adverse event (CTCAE v5.0, any grade, reported at ",
    "a routine visit) as a function of the individual model-predicted ",
    "steady-state ruxolitinib AUC over 24 h, in adults with ",
    "myeloproliferative neoplasms, graft-versus-host disease and other ",
    "haematologic or immune-mediated disorders treated in routine care ",
    "(Tachet 2026). logit(p) = logit_ref + log(1.88) * ",
    "log2(AUC_RUXOLITINIB): the odds of an adverse event rise 1.88-fold ",
    "per doubling of AUC24 (95% CI 1.08-3.29, P = 0.026). The paper prints ",
    "the odds ratio but not the intercept; the intercept here was ",
    "reconstructed by the maintainers from the Figure 3 per-visit ",
    "exposure distributions with the printed slope held fixed (see the ",
    "vignette). No PK layer and no ODE: AUC_RUXOLITINIB is supplied as ",
    "data, derived in the source from individual post hoc parameters of ",
    "the companion population PK model (Tachet_2026_ruxolitinib). ",
    "Companion of Tachet_2026_ruxolitinib_ae_ctrough."
  )
  reference <- "Tachet J, Meyer SC, Grandoni F, Haefliger E, Nagy-Hulliger M, Halter JP, Passweg JR, Haefliger D, Bardinet C, Decosterd LA, Guidi M, Girardin FR. Ruxolitinib Pharmacokinetics and Exposure-Toxicity Relationship in Hematologic Malignancies and Immune-Mediated Diseases: A Prospective Observational Study. Clin Pharmacol Ther. 2026;120:753-761. doi:10.1002/cpt.70367"
  vignette <- "Tachet_2026_ruxolitinib"
  units <- list(
    time = "n/a (static landmark exposure-toxicity regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the AUC_RUXOLITINIB covariate column)",
    concentration = "prob_ae (probability of any adverse event at a visit, 0-1; also logit_ae)"
  )

  covariateData <- list(
    AUC_RUXOLITINIB = list(
      description = "Individual model-predicted steady-state ruxolitinib AUC over 24 h (daily AUC) at the time of the toxicity assessment.",
      units = "ng*h/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "TOTAL plasma AUC over a 24 h day at STEADY STATE for the regimen",
        "the patient was on at the visit (daily dose / individual CL/F for",
        "a regular regimen), predicted from the individual post hoc",
        "parameters of the companion population PK model",
        "(modellib('Tachet_2026_ruxolitinib')); Tachet 2026 Methods,",
        "'Exploratory pharmacokinetic/pharmacodynamic analysis'. Enters on",
        "a log scale; the paper expresses the effect per two-fold increase,",
        "so the model uses log2(). The intercept absorbs the unit:",
        "supplying ug*h/mL instead of ng*h/mL would shift the logit by",
        "log(1.88) * log2(1000) = 6.3. Figure 3 per-visit quartiles read by",
        "the maintainers: no adverse event 891 / 1464 / 2355 ng*h/mL, any",
        "adverse event 1181 / 2314 / 3032 ng*h/mL."
      ),
      source_name = "AUC24 (area under the concentration-time curve over 24 h)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 77L,
    n_studies = 1L,
    n_observations = "visit-level toxicity assessments after removal of duplicate records with identical exposure and unchanged status; about 97 visits are plotted in Figure 3 (about 31 with any adverse event), the exact analysis count is not printed",
    age_range = "median 64 years, range 17-89",
    weight_range = "median 72.0 kg, range 40.7-116.0",
    sex_female_pct = 40,
    disease_state = "Graft-versus-host disease 43%; myeloproliferative neoplasms (myelofibrosis 30%, polycythaemia vera 16%, essential thrombocythaemia 3%); other haematologic or immune-mediated disease (VEXAS, APECED, HLH)",
    dose_range = "Ruxolitinib 5-25 mg orally once or twice daily, asymmetric schedules, or every 48-72 h",
    regions = "Switzerland",
    notes = "Any adverse event was documented in 101 of 209 samples with a recorded assessment (Tachet 2026 Table 1). Cmax was also tested and was not significantly associated with toxicity (P = 0.25); no odds ratio is printed for it, so no Cmax model is packaged."
  )

  ini({
    # Tachet 2026 Results, 'Exploratory PK/PD analysis': a two-fold increase in
    # AUC24 was associated with 1.9-fold higher odds of toxicity
    # (OR 1.88, 95% CI 1.08-3.29, P = 0.026).
    logit_ref <- -7.53; label("Logit of the probability of any adverse event at AUC_RUXOLITINIB = 1 ng*h/mL (unitless logit)") # not printed in the source; reconstructed by the maintainers from Tachet 2026 Figure 3 (score equation with the printed slope held fixed; see the vignette)
    e_auc_ruxolitinib_logit <- log(1.88); label("Log-odds of any adverse event per doubling of the steady-state ruxolitinib AUC24 (unitless logit)") # Results, 'Exploratory PK/PD analysis': OR 1.88 (95% CI 1.08-3.29) per two-fold increase in AUC24

    # A Bernoulli likelihood has no residual error and the source estimates no
    # random effect; this fixed placeholder only gives rxode2 an error model.
    addSd_prob_ae <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_ae <- logit_ref + e_auc_ruxolitinib_logit * log(AUC_RUXOLITINIB) / log(2)
    prob_ae <- expit(logit_ae)
    prob_ae ~ add(addSd_prob_ae)
  })
}
