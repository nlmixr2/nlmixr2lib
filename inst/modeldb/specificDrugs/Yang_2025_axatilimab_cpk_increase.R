Yang_2025_axatilimab_cpk_increase <- function() {
  description <- paste0(
    "Binomial logistic-regression exposure-safety model for creatine ",
    "phosphokinase elevation (grouped adverse event terms, any grade) ",
    "in patients with chronic graft-versus-host disease treated with ",
    "axatilimab, an anti-CSF-1R IgG4 monoclonal antibody (Yang 2025, n ",
    "= 278, SNDX-6352-0503 and AGAVE-201). The probability is ",
    "expit(-2.14 + 0.0042 * AUC_AXATILIMAB_SS + 0.0102 * (CPK - 64)) ",
    "where AUC_AXATILIMAB_SS is the individual axatilimab steady-state ",
    "AUC over a 4-week treatment cycle in ug*day/mL and continuous ",
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
    concentration = "prob_cpk_increase (probability of creatine phosphokinase elevation (grouped adverse event terms, any grade), 0-1; also logit_cpk_increase)"
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
    CPK = list(
      description = "Baseline serum creatine phosphokinase (creatine kinase) activity.",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "CENTRED at 64 U/L, the reference value the Yang 2025 Discussion ",
        "states ('relative to a patient with a reference CPK level of 64 ",
        "U/L') and the Table 2 median predictor. Table S3 prints a median ",
        "of 63.0 U/L over the 272 patients with a non-missing value; the ",
        "analysis dataset imputed the 6 missing values. Both Table 2 ",
        "percentile odds ratios are consistent with 64 (implied slope ",
        "0.0103 from each), not with 63."
      ),
      source_name = "CPK / BCPK (baseline CPK)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 278L,
    n_studies = 2L,
    n_observations = "98 with event of 278 (35.3 percent; Yang 2025 Table S6)",
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
    # Endpoint: creatine phosphokinase elevation (grouped adverse event terms, any grade).

    logit_ref <- -2.14 ; label("Logit of the event probability at zero axatilimab exposure with covariates at their reference values (unitless logit)")  # Yang 2025 Table S7, Intercept -2.14 (95% CI -2.64 to -1.69), P < 0.001
    e_auc_axatilimab_ss_logit <- 0.0042 ; label("Log-odds per 1 ug*day/mL increase in axatilimab AUCSS (unitless logit)")  # Yang 2025 Table S7, AUCSS estimate 0.0042 (95% CI 0.00308-0.00552), P < 0.001; Table 2 OR per 50 ug*day/mL 1.235 = exp(50 * 0.0042)
    e_cpk_logit <- 0.0102 ; label("Log-odds per 1 U/L increase in baseline CPK above the reference 64 U/L (unitless logit)")  # Yang 2025 Table S7, estimate 0.0102 (95% CI 0.00534-0.0157), P < 0.001, OR 1.01

    # glm(family = binomial) has a Bernoulli likelihood with no residual
    # parameter and no random effects. This small additive residual exists
    # only so rxode2 can attach an error model to the typical-value
    # probability; it is not a published quantity.
    addSd_prob_cpk_increase <- fixed(0.001) ; label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # Continuous covariates are centred at the analysis-population median,
    # as in the analysis script (Yang 2025 Data S2); the exposure term is not.
    logit_cpk_increase <- logit_ref +
      e_auc_axatilimab_ss_logit * AUC_AXATILIMAB_SS +
      e_cpk_logit * (CPK - 64)

    prob_cpk_increase <- expit(logit_cpk_increase)

    prob_cpk_increase ~ add(addSd_prob_cpk_increase)
  })
}
