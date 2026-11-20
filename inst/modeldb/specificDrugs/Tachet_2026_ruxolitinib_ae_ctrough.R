Tachet_2026_ruxolitinib_ae_ctrough <- function() {
  description <- paste0(
    "Univariable logistic-regression exposure-toxicity model for the ",
    "occurrence of any adverse event (CTCAE v5.0, any grade, reported at ",
    "a routine visit) as a function of the individual model-predicted ",
    "steady-state ruxolitinib trough concentration, in adults with ",
    "myeloproliferative neoplasms, graft-versus-host disease and other ",
    "haematologic or immune-mediated disorders treated in routine care ",
    "(Tachet 2026). logit(p) = logit_ref + log(1.51) * log2(CTROUGH): the ",
    "odds of an adverse event rise 1.51-fold per doubling of Cmin (95% CI ",
    "1.07-2.13, P = 0.019). The paper prints the odds ratio but not the ",
    "intercept; the intercept here was reconstructed by the maintainers ",
    "from the Figure 3 per-visit exposure distributions with the printed ",
    "slope held fixed (see the vignette). No PK layer and no ODE: CTROUGH ",
    "is supplied as data, derived in the source from individual post hoc ",
    "parameters of the companion population PK model ",
    "(Tachet_2026_ruxolitinib). Companion of Tachet_2026_ruxolitinib_ae_auc24."
  )
  reference <- "Tachet J, Meyer SC, Grandoni F, Haefliger E, Nagy-Hulliger M, Halter JP, Passweg JR, Haefliger D, Bardinet C, Decosterd LA, Guidi M, Girardin FR. Ruxolitinib Pharmacokinetics and Exposure-Toxicity Relationship in Hematologic Malignancies and Immune-Mediated Diseases: A Prospective Observational Study. Clin Pharmacol Ther. 2026;120:753-761. doi:10.1002/cpt.70367"
  vignette <- "Tachet_2026_ruxolitinib"
  units <- list(
    time = "n/a (static landmark exposure-toxicity regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CTROUGH covariate column)",
    concentration = "prob_ae (probability of any adverse event at a visit, 0-1; also logit_ae)"
  )

  covariateData <- list(
    CTROUGH = list(
      description = "Individual model-predicted steady-state ruxolitinib trough (minimum) plasma concentration at the time of the toxicity assessment.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "TOTAL plasma concentration at STEADY STATE for the regimen the",
        "patient was on at the visit, predicted from the individual post",
        "hoc parameters of the companion population PK model",
        "(modellib('Tachet_2026_ruxolitinib')); Tachet 2026 Methods,",
        "'Exploratory pharmacokinetic/pharmacodynamic analysis'. Enters on",
        "a log scale; the paper expresses the effect per two-fold increase,",
        "so the model uses log2(). The intercept absorbs the unit: supplying",
        "ug/mL instead of ng/mL would shift the logit by",
        "log(1.51) * log2(1000) = 4.1. Figure 3 per-visit quartiles read by",
        "the maintainers: no adverse event 4.8 / 12.7 / 25.2 ng/mL, any",
        "adverse event 9.4 / 25.2 / 48.5 ng/mL."
      ),
      source_name = "Cmin (trough concentration)"
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
    notes = "Any adverse event was documented in 101 of 209 samples with a recorded assessment (Tachet 2026 Table 1); metabolism and nutrition disorders and blood and lymphatic system disorders were the most frequent categories (Supplementary Table S3)."
  )

  ini({
    # Tachet 2026 Results, 'Exploratory PK/PD analysis': a two-fold increase in
    # Cmin was associated with 1.5-fold higher odds of toxicity
    # (OR 1.51, 95% CI 1.07-2.13, P = 0.019).
    logit_ref <- -2.40; label("Logit of the probability of any adverse event at CTROUGH = 1 ng/mL (unitless logit)") # not printed in the source; reconstructed by the maintainers from Tachet 2026 Figure 3 (score equation with the printed slope held fixed; see the vignette)
    e_ctrough_logit <- log(1.51); label("Log-odds of any adverse event per doubling of the steady-state ruxolitinib trough concentration (unitless logit)") # Results, 'Exploratory PK/PD analysis': OR 1.51 (95% CI 1.07-2.13) per two-fold increase in Cmin

    # A Bernoulli likelihood has no residual error and the source estimates no
    # random effect; this fixed placeholder only gives rxode2 an error model.
    addSd_prob_ae <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_ae <- logit_ref + e_ctrough_logit * log(CTROUGH) / log(2)
    prob_ae <- expit(logit_ae)
    prob_ae ~ add(addSd_prob_ae)
  })
}
