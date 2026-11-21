Yang_2025_axatilimab_mlss_improve7 <- function() {
  description <- paste0(
    "Binomial logistic-regression exposure-efficacy model for ",
    "clinically meaningful symptom improvement, a >= 7-point reduction ",
    "from baseline in the modified Lee Symptom Scale (mLSS) summary ",
    "score within the first 6 cycles in patients with chronic ",
    "graft-versus-host disease treated with axatilimab, an anti-CSF-1R ",
    "IgG4 monoclonal antibody (Yang 2025, n = 239, AGAVE-201). The ",
    "probability is expit(-0.114 - 0.00118 * AUC_AXATILIMAB_SD + 0.77 ",
    "* CGVHD_ORGANS_GT4) where AUC_AXATILIMAB_SD is the individual ",
    "axatilimab single-dose AUC from time 0 to infinity in ug*day/mL ",
    "and continuous covariates are centred at their ",
    "analysis-population medians (Yang 2025 Table S5). There is no PK ",
    "layer: the exposure metric is supplied as a data column derived ",
    "from the companion population PK/PD model (Yang_2024_axatilimab). ",
    "No random effect and no residual error are estimated (Bernoulli ",
    "likelihood). One of twelve exposure-response models in the ",
    "Yang_2025_axatilimab_* family."
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
    concentration = "prob_mlss_improve7 (probability of clinically meaningful symptom improvement, a >= 7-point reduction from baseline in the modified Lee Symptom Scale (mLSS) summary score within the first 6 cycles, 0-1; also logit_mlss_improve7)"
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
    ),
    CGVHD_ORGANS_GT4 = list(
      description = "More than four organs involved by chronic graft-versus-host disease at baseline; 1 = more than 4 organs, 0 = 4 or fewer.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = paste0(
        "Efficacy population 77 of 239 (32.2 percent; Yang 2025 Table S2)."
      ),
      source_name = "Number of organs involved at baseline (> 4 vs <= 4)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 239L,
    n_studies = 1L,
    n_observations = "114 responders of 239 (47.7 percent; Yang 2025 Table S4)",
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
    # Endpoint: clinically meaningful symptom improvement, a >= 7-point reduction from baseline in the modified Lee Symptom Scale (mLSS) summary score within the first 6 cycles.

    logit_ref <- -0.114 ; label("Logit of the event probability at zero axatilimab exposure with covariates at their reference values (unitless logit)")  # Yang 2025 Table S5, Intercept -0.114 (95% CI -0.491-0.26), P 0.55
    e_auc_axatilimab_sd_logit <- -0.00118 ; label("Log-odds per 1 ug*day/mL increase in axatilimab AUCSDINF (unitless logit)")  # Yang 2025 Table S5, AUCSDINF estimate -0.00118 (95% CI -0.00228 to -0.000133), P 0.0304; Table 1 OR per 50 ug*day/mL 0.942 = exp(50 * -0.00118)
    e_cgvhd_organs_gt4_logit <- 0.77 ; label("Log-odds for CGVHD_ORGANS_GT4 = 1 versus 0 (unitless logit)")  # Yang 2025 Table S5, estimate 0.77 (95% CI 0.214-1.34), P 0.0071, OR 2.16

    # glm(family = binomial) has a Bernoulli likelihood with no residual
    # parameter and no random effects. This small additive residual exists
    # only so rxode2 can attach an error model to the typical-value
    # probability; it is not a published quantity.
    addSd_prob_mlss_improve7 <- fixed(0.001) ; label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli")  # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_mlss_improve7 <- logit_ref +
      e_auc_axatilimab_sd_logit * AUC_AXATILIMAB_SD +
      e_cgvhd_organs_gt4_logit * CGVHD_ORGANS_GT4

    prob_mlss_improve7 <- expit(logit_mlss_improve7)

    prob_mlss_improve7 ~ add(addSd_prob_mlss_improve7)
  })
}
