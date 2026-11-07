Menon_2025_tofacitinib_hb_decrease <- function() {
  description <- paste(
    "Logistic exposure-safety model for the probability of a hemoglobin",
    "decrease of more than 2 g/dL from baseline in adults with active",
    "psoriatic arthritis given placebo or tofacitinib 5 or 10 mg twice daily",
    "in OPAL Broaden and OPAL Beyond (Menon 2025). The logit is linear in the",
    "individual steady-state average tofacitinib concentration (CAV,",
    "supplied from a population PK model such as Xie_2019_tofacitinib).",
    "The intercept and slope are not printed in the paper; they are",
    "back-solved from the printed placebo-adjusted incidences at the 5 and",
    "10 mg twice-daily median Cavg."
  )
  reference <- paste(
    "Menon S, Shoji S, Tsuchiwata S, Fallon L, Kanik K. Exposure-Response",
    "Analysis of Tofacitinib in Active Psoriatic Arthritis: Results from Two",
    "Phase 3 Studies. J Clin Pharmacol. 2025;65(3):369-377.",
    "doi:10.1002/jcph.6147. PMID: 39453735. PMCID: PMC11867917."
  )
  vignette <- "Menon_2025_tofacitinib"
  units <- list(
    time = "day (landmark model; the probability does not depend on time)",
    dosing = "mg",
    concentration = "prob_hb_decrease_2g (probability of a > 2 g/dL hemoglobin decrease from baseline during the analysed treatment period, 0-1)"
  )

  covariateData <- list(
    CAV = list(
      description = "Individual average steady-state tofacitinib plasma concentration over the twice-daily dosing interval.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Patient-specific Cavg predicted from the tofacitinib psoriatic",
        "arthritis population PK model of Xie 2019 (Menon 2025 Methods;",
        "modellib('Xie_2019_tofacitinib')). 0 for placebo, so the intercept",
        "is the placebo incidence. Enters the logit linearly and uncentred.",
        "Median (range) 17.6 (6.67-47.2) and 36.1 (15.1-88.1) ng/mL at 5 and",
        "10 mg twice daily (Results). Time-fixed per patient."
      ),
      source_name = "Cavg"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 685L,
    n_studies = 2L,
    age_range = "adults >= 18 years (median 50 years)",
    weight_range = "median 83 kg",
    sex_female_pct = 54.9,
    disease_state = paste(
      "Active psoriatic arthritis (CASPAR criteria, >= 3 swollen and >= 3",
      "tender joints, active plaque psoriasis) with inadequate response to",
      "conventional synthetic DMARDs (OPAL Broaden) or to TNF inhibitors",
      "(OPAL Beyond), on a stable background csDMARD"
    ),
    dose_range = "Placebo (up to Month 3) or tofacitinib 5 or 10 mg orally twice daily (up to Month 6 or Month 12)",
    regions = "Multinational (OPAL Broaden NCT01877668, OPAL Beyond NCT01882439)",
    baseline_hemoglobin = "median 13.8 g/dL",
    notes = paste(
      "Menon 2025 Results: the 685-patient hemoglobin analysis set (309",
      "male, 376 female). The paper does not state the number of patients",
      "with a > 2 g/dL decrease or the observation window per patient beyond",
      "the analysed hemoglobin data (placebo to Month 3, tofacitinib to Month",
      "6 or 12). Fitted in NONMEM 7.3."
    )
  )

  ini({
    # logit(Pr(Y = 1)) = theta0 + theta_cov * Cavg (Menon 2025 Methods).
    # theta0 and theta_cov are not printed. Solving
    #   expit(theta0 + 17.6 theta_cov) - expit(theta0) = 0.0108
    #   expit(theta0 + 36.1 theta_cov) - expit(theta0) = 0.0301
    # for the printed placebo-adjusted incidences at the median Cavg of the
    # 5 and 10 mg twice-daily groups (Results) gives the values below. They
    # predict 7.82% at the printed 90th-percentile Cavg of 54.8 ng/mL, where
    # the paper prints 7.97%, and about 19% at 88 ng/mL, the right edge of
    # the fitted curve in Figure 2a.
    logit_ref <- -4.168; label("Logit of the probability of a > 2 g/dL hemoglobin decrease at zero tofacitinib exposure (logit)")  # back-solved from Results placebo-adjusted incidences 1.08% (Cavg 17.6) and 3.01% (Cavg 36.1); not printed
    e_cav_logit <- 0.03105; label("Log-odds of a > 2 g/dL hemoglobin decrease per 1 ng/mL tofacitinib Cavg (logit per ng/mL)")  # back-solved from Results placebo-adjusted incidences 1.08% (Cavg 17.6) and 3.01% (Cavg 36.1); not printed

    # Bernoulli likelihood in the source; the fixed additive residual is a
    # placeholder so rxode2 has an error model, not a published quantity.
    addSd_prob_hb_decrease_2g <- fixed(0.001); label("Placeholder additive residual SD on prob_hb_decrease_2g (not from source)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_hb_decrease_2g <- logit_ref + e_cav_logit * CAV

    prob_hb_decrease_2g <- expit(logit_hb_decrease_2g)

    prob_hb_decrease_2g ~ add(addSd_prob_hb_decrease_2g)
  })
}
