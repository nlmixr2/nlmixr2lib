Yang_2025_axatilimab_sae <- function() {
  description <- paste0(
    "Binomial logistic-regression exposure-safety model for any ",
    "serious treatment-emergent adverse event in patients with chronic ",
    "graft-versus-host disease treated with axatilimab, an anti-CSF-1R ",
    "IgG4 monoclonal antibody (Yang 2025, n = 278, SNDX-6352-0503 and ",
    "AGAVE-201). The probability is expit(-0.736 + 0.000546 * ",
    "AUC_AXATILIMAB_SS - 0.0964 * (ALB - 40) + 0.687 * ECOG_GE2) where ",
    "AUC_AXATILIMAB_SS is the individual axatilimab steady-state AUC ",
    "over a 4-week treatment cycle in ug*day/mL and continuous ",
    "covariates are centred at their analysis-population medians (Yang ",
    "2025 Table S7). There is no PK layer: the exposure metric is ",
    "supplied as a data column derived from the companion population ",
    "PK/PD model (Yang_2024_axatilimab). No random effect and no ",
    "residual error are estimated (Bernoulli likelihood). One of ",
    "twelve exposure-response models in the Yang_2025_axatilimab_* ",
    "family."
  )
  reference <- paste(
    "Yang Y, Volkova A, Sokolov V, Liu X, Leon C, Kosinsky Y, Sheng J, Chen X.",
    "Exposure-Response Relationships for Axatilimab in Patients with Chronic",
    "Graft-Versus-Host Disease. Clin Pharmacol Ther. 2025;118(4):853-864.",
    "doi:10.1002/cpt.3731.",
    "Individual axatilimab exposure metrics derive from the companion",
    "population PK/PD model: Yang YO, Sokolov V, Volkova A, et al.",
    "Semimechanistic Population PK/PD Modeling of Axatilimab in Healthy",
    "Participants and Patients With Solid Tumors or Chronic Graft-Versus-Host",
    "Disease. Clin Pharmacol Ther. 2025;117(3):704-715. doi:10.1002/cpt.3503",
    "(packaged as Yang_2024_axatilimab).",
    sep = " "
  )
  vignette <- "Yang_2025_axatilimab_exposure_response"
  units <- list(
    time = "n/a (static landmark exposure-response regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the AUC_AXATILIMAB_SS covariate column)",
    concentration = "prob_sae (probability of any serious treatment-emergent adverse event, 0-1; also logit_sae)"
  )

  covariateData <- list(
    AUC_AXATILIMAB_SS = list(
      description = "Individual axatilimab steady-state area under the plasma concentration-time curve over one 4-week treatment cycle (AUCSS; two Q2W doses or one Q4W dose).",
      units = "ug*day/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Supplied as data: this model has no PK layer. Yang 2025 simulated ",
        "it per patient and per endpoint from the empirical Bayes ",
        "estimates of the companion population PK/PD model ",
        "(doi:10.1002/cpt.3503, packaged as Yang_2024_axatilimab), using ",
        "the patient's cycle-average dose from day 1 up to the time of ",
        "first event incidence (or end of treatment when no event ",
        "occurred) so that dose modifications are reflected. Enters ",
        "UNCENTRED and per 1 ug*day/mL; Yang 2025 tabulates the odds ratio ",
        "per 50 ug*day/mL (OR = exp(50 * slope)). Safety population median ",
        "about 115-133 ug*day/mL depending on the endpoint (5th-95th ",
        "percentile about 18-1440; Yang 2025 Table 2)."
      ),
      source_name = "AUCSS"
    ),
    ALB = list(
      description = "Baseline serum albumin.",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "CENTRED at the safety-population median 40 g/L: the analysis ",
        "script (Yang 2025 Data S2) subtracts the median from every ",
        "continuous covariate before fitting ('dataset with covariates ",
        "normalized to the median value'). The centring value is ",
        "self-pinned by Table 2: exp(-0.0964 * (32.085 - 40)) = 2.14 and ",
        "exp(-0.0964 * (47 - 40)) = 0.509 reproduce the printed ",
        "5th/95th-percentile odds ratios 2.144 and 0.509. Table 2 labels ",
        "the row 'ALB, U/L'; Table S7 and the Discussion give g/L. Safety ",
        "population median 40.0 g/L, range 22.0-53.0 (Table S3)."
      ),
      source_name = "ALB (baseline albumin)"
    ),
    ECOG_GE2 = list(
      description = "Decreased baseline functional status ('unable to work'): baseline Karnofsky (adults) or Lansky (children) performance score <= 70; 1 = yes, 0 = score >= 80.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = paste0(
        "Yang 2025 codes this as 'Unable to work (baseline ",
        "Karnofsky/Lansky scores <= 70%)'. Under the Oken 1982 crosswalk ",
        "Karnofsky 60-70 corresponds to ECOG 2 ('unable to carry out any ",
        "work activities') and Karnofsky <= 50 to ECOG >= 3, so a ",
        "Karnofsky/Lansky score <= 70 is ECOG >= 2. Supply 1 when the ",
        "Karnofsky/Lansky score is <= 70 (or ECOG >= 2). Safety ",
        "population: 99 of 277 patients with a known score (35.7 percent; ",
        "Table S3)."
      ),
      source_name = "PERFSC (Karnofsky/Lansky performance score, dichotomised at <= 70)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 278L,
    n_studies = 2L,
    n_observations = "118 with event of 278 (42.4 percent; Yang 2025 Table S6)",
    age_range = "median 55.0 years, range 7.0-81.0 (Yang 2025 Table S3)",
    weight_range = "median 72.4 kg, range 18.1-151 (Yang 2025 Table S3)",
    sex_female_pct = 37.4,
    race_ethnicity = c(White = 83.1, Asian = 5.8, Black = 3.6, Other = 1.8, Missing = 5.8),
    disease_state = "chronic graft-versus-host disease; cGVHD severity at baseline severe 80.6 percent, moderate 2.2 percent, mild 17.3 percent; Karnofsky/Lansky score <= 70 in 35.7 percent; mild hepatic impairment 21.2 percent (Yang 2025 Table S3)",
    dose_range = "axatilimab IV 0.3 mg/kg Q2W (n = 79), 1 mg/kg Q2W (n = 107), 3 mg/kg Q4W (n = 85), 3 mg/kg Q2W (n = 6), or 0.15 / 0.5 mg/kg Q2W (n = 2)",
    regions = "Phase 1/2 study SNDX-6352-0503 (NCT03604692, n = 39) and multinational phase 2 study AGAVE-201 (NCT04710576, n = 239)",
    notes = "Exposure-safety analysis set: safety-evaluable patients with at least one axatilimab dose and one quantifiable post-dose concentration (278 of 279). Median treatment duration 7.1 months (SNDX-6352-0503) and 7.8 months (AGAVE-201). Adverse events are treatment-emergent regardless of causality."
  )

  ini({
    # Yang 2025 Supplementary Table S5/S7 prints the fitted glm() coefficients
    # (intercept and per-unit slopes); main-text Table 1/2 prints the same
    # models as odds ratios, per 50 ug*day/mL for the exposure term.
    # Endpoint: any serious treatment-emergent adverse event.

    logit_ref <- -0.736 ; label("Logit of the event probability at zero axatilimab exposure with covariates at their reference values (unitless logit)")  # Yang 2025 Table S7, Intercept -0.736 (95% CI -1.1 to -0.39), P < 0.001
    e_auc_axatilimab_ss_logit <- 0.000546 ; label("Log-odds per 1 ug*day/mL increase in axatilimab AUCSS (unitless logit)")  # Yang 2025 Table S7, AUCSS estimate 0.000546 (95% CI 0.0000154-0.00112), P 0.0484; Table 2 OR per 50 ug*day/mL 1.028 = exp(50 * 0.000546)
    e_alb_logit <- -0.0964 ; label("Log-odds per 1 g/L increase in baseline ALB above the reference 40 g/L (unitless logit)")  # Yang 2025 Table S7, estimate -0.0964 (95% CI -0.157 to -0.0393), P 0.00125, OR 0.909
    e_ecog_ge2_logit <- 0.687 ; label("Log-odds for ECOG_GE2 = 1 versus 0 (unitless logit)")  # Yang 2025 Table S7, estimate 0.687 (95% CI 0.175-1.2), P 0.00874, OR 2.002

    # glm(family = binomial) has a Bernoulli likelihood with no residual
    # parameter and no random effects. This small additive residual exists
    # only so rxode2 can attach an error model to the typical-value
    # probability; it is not a published quantity.
    addSd_prob_sae <- fixed(0.001) ; label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # Continuous covariates are centred at the analysis-population median,
    # as in the analysis script (Yang 2025 Data S2); the exposure term is not.
    logit_sae <- logit_ref +
      e_auc_axatilimab_ss_logit * AUC_AXATILIMAB_SS +
      e_alb_logit * (ALB - 40) +
      e_ecog_ge2_logit * ECOG_GE2

    prob_sae <- expit(logit_sae)

    prob_sae ~ add(addSd_prob_sae)
  })
}
