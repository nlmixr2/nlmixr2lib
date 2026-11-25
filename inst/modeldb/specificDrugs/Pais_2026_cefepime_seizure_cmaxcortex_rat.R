Pais_2026_cefepime_seizure_cmaxcortex_rat <- function() {
  description <- paste(
    "Preclinical (rat). Logistic exposure-toxicity model for the probability of a",
    "convulsive seizure (modified Racine stage > 1) in male Sprague-Dawley rats with",
    "folic-acid-induced acute kidney injury given intravenous cefepime, driven by the",
    "first-24-h cefepime cerebral cortex Cmax0-24 estimated from the companion population PK model",
    "(Pais_2026_cefepime_rat). The logit is linear in the exposure and crosses zero at the",
    "published TC50. One of five independently fitted logistic regressions in the paper",
    "(Figure 5); no PK layer and no ODE -- the exposure is supplied as data.",
    sep = " "
  )
  reference <- paste(
    "Pais GM, Lesnicki E, Marianski S, Valdez K, Gibson Z, Christopher J, Lepard K,",
    "Gilchrist A, Mizuno T, Scheetz MH.",
    "Elevated hippocampal and cortical cefepime concentrations correlate with seizures",
    "in an acute kidney injury rat model.",
    "Antimicrob Agents Chemother. 2026;70(4):e01005-25. doi:10.1128/aac.01005-25.",
    sep = " "
  )
  vignette <- "Pais_2026_cefepime_rat"
  units <- list(
    time = "n/a (static exposure-toxicity model; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CMAX_CEFEPIME_BRAIN_CORTEX covariate column)",
    concentration = "prob_racine_gt1 (probability of modified Racine seizure stage > 1, 0-1)"
  )

  covariateData <- list(
    CMAX_CEFEPIME_BRAIN_CORTEX = list(
      description = "Per-animal cefepime cerebral cortex Cmax0-24 over the first 24 h of treatment (from the first dose), in cerebral cortex tissue.",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Pais 2026 computed it from the empirical Bayes concentration profile of each rat",
        "given its exact dosing (Simulx, 0.1-h grid), using the Pais_2026_cefepime_rat PK",
        "model. Observed range about 5-91 mg/L (Figure 5d). Total (not unbound) cefepime.",
        "Brain-tissue values are in ug/g tissue, treated as mg/L by the source."
      ),
      source_name = "Cmax_cortex (Figure 5d x-axis)"
    )
  )

  population <- list(
    species = "rat (Sprague-Dawley; male)",
    n_subjects = 18L,
    n_studies = 1L,
    weight_range = "mean 282 g (SD 17.0 g)",
    sex_female_pct = 0,
    disease_state = "Folic-acid-induced acute kidney injury; seizures scored on a modified Racine scale (0-7) and dichotomised as stage <= 1 (non-convulsive) versus > 1 (convulsive)",
    dose_range = "Cefepime 1,593 mg/kg once or 1,250 mg/kg once daily for up to 5 days, intravenously",
    regions = "United States (Midwestern University, Downers Grove, Illinois)",
    notes = "Seizure stage was the last observed stage before brain harvest. Logistic regression fitted in Stata BE 17.0 (one binary outcome per rat)."
  )

  ini({
    # TC50 is printed in Pais 2026 Results and the Figure 5 legend. The slope is
    # not printed: it was digitised from the Figure 5d predicted-probability
    # markers (red, Pr(seizurestage1)) by least squares on logit(p) with the
    # printed TC50 held. Fitted freely, the same markers give slope 0.0590 and
    # TC50 44.9, reproducing the printed TC50 to within 1%.
    tc50 <- 45.0; label("Exposure (Cmax0-24) at a 50% probability of seizure stage > 1, TC50 (mg/L)") # Results / Figure 5 legend TC50 = 45.0 mg/L
    logit_slope <- 0.0590; label("Logistic slope of seizure log-odds on Cmax0-24 (per mg/L)") # digitised from Figure 5d; not printed in the source

    # The source likelihood is Bernoulli; this placeholder lets the model be used
    # with the nlmixr2 machinery and is not from the source.
    addSd_prob_racine_gt1 <- fixed(0.001); label("Placeholder additive residual SD on the seizure probability (not from the source)") # not from source
  })

  model({
    logit_racine_gt1 <- logit_slope * (CMAX_CEFEPIME_BRAIN_CORTEX - tc50)
    prob_racine_gt1 <- expit(logit_racine_gt1)
    prob_racine_gt1 ~ add(addSd_prob_racine_gt1)
  })
}
