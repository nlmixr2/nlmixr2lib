Yang_2025_axatilimab_orr <- function() {
  description <- paste0(
    "Binomial logistic-regression exposure-efficacy model for overall ",
    "response (complete or partial response of chronic ",
    "graft-versus-host disease) by Cycle 7 Day 1 in patients with ",
    "chronic graft-versus-host disease treated with axatilimab, an ",
    "anti-CSF-1R IgG4 monoclonal antibody (Yang 2025, n = 239, ",
    "AGAVE-201). The probability is expit(0.899 - 0.00157 * ",
    "AUC_AXATILIMAB_SD) where AUC_AXATILIMAB_SD is the individual ",
    "axatilimab single-dose AUC from time 0 to infinity in ug*day/mL ",
    "(Yang 2025 Table S5). There is no PK layer: the exposure metric ",
    "is supplied as a data column derived from the companion ",
    "population PK/PD model (Yang_2024_axatilimab). No random effect ",
    "and no residual error are estimated (Bernoulli likelihood). One ",
    "of twelve exposure-response models in the Yang_2025_axatilimab_* ",
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
    dosing = "n/a (no dose events; exposure enters as the AUC_AXATILIMAB_SD covariate column)",
    concentration = "prob_orr (probability of overall response (complete or partial response of chronic graft-versus-host disease) by Cycle 7 Day 1, 0-1; also logit_orr)"
  )

  covariateData <- list(
    AUC_AXATILIMAB_SD = list(
      description = "Individual axatilimab area under the plasma concentration-time curve from time 0 to infinity after a single dose (AUCSDINF) at the patient's nominal starting dose and regimen.",
      units = "ug*day/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Supplied as data: this model has no PK layer. Yang 2025 derived ",
        "it per patient from the empirical Bayes estimates of the ",
        "companion population PK/PD model (doi:10.1002/cpt.3503, packaged ",
        "as Yang_2024_axatilimab) with body weight, baseline CSF-1 and ADA ",
        "status (ADA treated as time-invariant: 1 if ADA-positive at any ",
        "time point) as individual covariates. Enters UNCENTRED and per 1 ",
        "ug*day/mL; Yang 2025 tabulates the odds ratio per 50 ug*day/mL ",
        "(OR = exp(50 * slope)). The paper selected AUCSDINF because it ",
        "equals AUCtau at steady state under linear kinetics. Efficacy ",
        "population median 64.0 ug*day/mL (5th-95th percentile 8.99-714.5; ",
        "Yang 2025 Table 1)."
      ),
      source_name = "AUCSDINF"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 239L,
    n_studies = 1L,
    n_observations = "153 responders of 239 (64 percent; Yang 2025 Table S4)",
    age_range = "median 53.0 years, range 7.0-81.0 (Yang 2025 Table S2)",
    weight_range = "median 71.6 kg, range 18.1-151 (Yang 2025 Table S2)",
    sex_female_pct = 37.7,
    race_ethnicity = c(White = 82.8, Asian = 6.7, Black = 2.1, Other = 1.7, Missing = 6.7),
    disease_state = "chronic graft-versus-host disease after at least two prior lines of systemic therapy; cGVHD severity at baseline severe 79.9 percent, mild 20.1 percent; more than 4 organs involved 32.2 percent; prior ruxolitinib 74.5 percent, prior ibrutinib 31.0 percent, prior belumosudil 23.0 percent (Yang 2025 Table S2)",
    dose_range = "axatilimab IV 0.3 mg/kg Q2W (n = 79), 1 mg/kg Q2W (n = 81) or 3 mg/kg Q4W (n = 79)",
    regions = "Multinational pivotal phase 2 study AGAVE-201 (NCT04710576)",
    notes = "Exposure-efficacy analysis set: AGAVE-201 efficacy-evaluable patients with at least one axatilimab dose and one quantifiable post-dose concentration (239 of 241). ADA-positive 32.6 percent. Efficacy assessed through Cycle 7 Day 1 (or end of treatment / start of new anti-cGVHD therapy)."
  )

  ini({
    # Yang 2025 Supplementary Table S5/S7 prints the fitted glm() coefficients
    # (intercept and per-unit slopes); main-text Table 1/2 prints the same
    # models as odds ratios, per 50 ug*day/mL for the exposure term.
    # Endpoint: overall response (complete or partial response of chronic graft-versus-host disease) by Cycle 7 Day 1.

    logit_ref <- 0.899 ; label("Logit of the event probability at zero axatilimab exposure (unitless logit)")  # Yang 2025 Table S5, Intercept 0.899 (95% CI 0.558-1.25), P < 0.001
    e_auc_axatilimab_sd_logit <- -0.00157 ; label("Log-odds per 1 ug*day/mL increase in axatilimab AUCSDINF (unitless logit)")  # Yang 2025 Table S5, AUCSDINF estimate -0.00157 (95% CI -0.00262 to -0.000542), P 0.00299; Table 1 OR per 50 ug*day/mL 0.924 = exp(50 * -0.00157)

    # glm(family = binomial) has a Bernoulli likelihood with no residual
    # parameter and no random effects. This small additive residual exists
    # only so rxode2 can attach an error model to the typical-value
    # probability; it is not a published quantity.
    addSd_prob_orr <- fixed(0.001) ; label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli")  # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_orr <- logit_ref +
      e_auc_axatilimab_sd_logit * AUC_AXATILIMAB_SD

    prob_orr <- expit(logit_orr)

    prob_orr ~ add(addSd_prob_orr)
  })
}
