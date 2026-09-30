Shulgin_2020_ici_trae34_mbma <- function() {
  description <- "MBMA. Study-level logit meta-regression of the grade 3/4 treatment-related adverse event (trAE) rate on normalized anti-CTLA-4 exposure for immune checkpoint inhibitor (ICI) therapy with the anti-PD-1 antibodies nivolumab / pembrolizumab and the anti-CTLA-4 antibodies ipilimumab / tremelimumab, alone or in combination (Shulgin 2020 'Model 4', CTLA-4 inhibitor-driven with PD-1 inhibitor-dependent modulation; 147 cohorts). Exposure is the steady-state average concentration of the anti-CTLA-4 antibody divided by its in vitro IC50; PD-1 inhibitor exposure has no effect, and an anti-PD-1 antibody in the regimen steepens the CTLA-4 exposure slope roughly eight-fold. This is a STUDY-LEVEL model: it predicts the expected grade 3/4 trAE proportion of a trial cohort, not an individual-patient risk, and it consumes no dose events (exposure arrives as covariate columns)."
  reference <- paste(
    "Shulgin B, Kosinsky Y, Omelchenko A, Chu L, Mugundu G, Aksenov S,",
    "Pimentel R, DeYulia G, Kim G, Peskov K, Helmlinger G.",
    "Dose dependence of treatment-related adverse events for immune checkpoint",
    "inhibitor therapies: a model-based meta-analysis.",
    "Oncoimmunology. 2020;9(1):1748982. doi:10.1080/2162402X.2020.1748982.",
    sep = " "
  )
  vignette <- "Shulgin_2020_ici_adverse_events_mbma"

  # Algebraic study-level MBMA: no rxode2 dose events are consumed and the
  # output is a probability, not a drug concentration (same device as
  # Chen_2025_hemoporfin_patient_rating).
  units <- list(
    time = "n/a (static study-level meta-regression; no time dimension)",
    dosing = "n/a (no dose events; anti-CTLA-4 exposure enters as the covariates CAV and IC50_CTLA4)",
    concentration = "prob_trae_grade34 (expected proportion of a cohort with a grade 3/4 treatment-related adverse event, 0-1)"
  )

  covariateData <- list(
    CAV = list(
      description = "Cohort-level steady-state average serum concentration of the anti-CTLA-4 antibody (ipilimumab or tremelimumab) at the cohort's dosing regimen; 0 when the regimen contains no anti-CTLA-4 antibody.",
      units = "nM",
      type = "continuous",
      reference_category = NULL,
      notes = "Shulgin 2020 Supplemental Methods: concentrations (nM) were simulated with published two-compartment linear population PK models (Supplemental Methods table: ipilimumab Vc 4.15 L, CL 0.0150 L/h, Vp 3.11 L, Q 0.0411 L/h; tremelimumab Vc 3.72 L, CL 0.0109 L/h, Vp 3.31 L, Q 0.0172 L/h) and 'averaged over time intervals corresponding to between 3 and 4 delivered doses, to assume steady-state conditions' (read here as the dosing interval between the 3rd and 4th doses). Divided by IC50_CTLA4 inside model() to give the paper's normalized exposure 'dose-CTLA-4' (C_norm,CTLA4). Anti-PD-1 exposure is NOT an input: it was not retained in Model 4. The vignette shows how to derive CAV from a dosing regimen.",
      source_name = "C_av (Shulgin 2020 Supplemental Methods); 'Cnorm,CTLA4' / 'dose-CTLA-4' after normalization"
    ),
    IC50_CTLA4 = list(
      description = "In vitro IC50 of the anti-CTLA-4 antibody for inhibition of ligand binding to human CTLA-4, used to normalize CAV.",
      units = "nM",
      type = "continuous",
      reference_category = NULL,
      notes = "Shulgin 2020 Supplemental Methods table prints ipilimumab 13.3 nM and tremelimumab 0.65 nM (human targets expressed on CHO or HEK cells). The coefficients were fitted on the exposure axis of main-text Figure 2, and on that axis every ipilimumab regimen sits 10-fold higher than 13.3 nM gives (ipilimumab 1 mg/kg Q3W, 3 mg/kg Q3W and 10 mg/kg Q3W plot at about 45, 140 and 430) while tremelimumab and both anti-PD-1 antibodies plot where the printed values put them. Supply 1.33 nM for ipilimumab to reproduce the paper's calibration; the vignette shows the derivation. Must be strictly positive in every row, including rows whose regimen has no anti-CTLA-4 antibody (where CAV = 0); a zero or missing IC50 propagates NaN into the prediction rather than defaulting silently.",
      source_name = "IC50 (Shulgin 2020 Supplemental Methods table)"
    ),
    TRT_ANTIPD1 = list(
      description = "Indicator that the cohort's regimen includes an anti-PD-1 antibody (nivolumab or pembrolizumab), as monotherapy or in combination.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no anti-PD-1 antibody in the regimen -- anti-CTLA-4 monotherapy)",
      notes = "The paper's 'FactorPD1' / 'factor-PD-1': 'set to 1 if a PD-1 inhibitor drug was given as monotherapy or in combination, and to 0 otherwise' (Shulgin 2020 Methods). In Model 4 it enters only through the interaction with normalized anti-CTLA-4 exposure, so for anti-PD-1 monotherapy (CAV = 0) its value does not change the prediction.",
      source_name = "FactorPD1 / factor-PD-1 (Shulgin 2020 Methods, Supplemental Tables 1-3)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 21305L,
    n_studies = 80L,
    n_cohorts = 147L,
    age_range = "Cohort median ages 42-74 years (mean of cohort medians 60.6, SD 4.76; 128 cohorts reporting; Shulgin 2020 Supplemental Table 4).",
    sex_female_pct = 38.8,
    disease_state = "Advanced solid tumours treated with PD-1 and/or CTLA-4 immune checkpoint inhibitors -- predominantly melanoma and non-small cell lung cancer, plus renal cell, small-cell lung, head and neck, gastric, urothelial, prostate, mesothelioma, sarcoma and other cancers (Shulgin 2020 Supplemental Table 5). Cancer type was screened and not retained.",
    dose_range = "Nivolumab 0.1-10 mg/kg Q2W/Q3W or 240 mg Q2W; pembrolizumab 2-10 mg/kg Q2W/Q3W or 200 mg Q3W; ipilimumab 0.3-10 mg/kg Q3W/Q4W as monotherapy and 1 or 3 mg/kg Q3W/Q6W/Q12W in combination with an anti-PD-1 antibody; tremelimumab 10 mg/kg Q4W or 15 mg/kg Q12W; alone, as PD-1 + CTLA-4 combination (27 cohorts), or with standard chemotherapy (Shulgin 2020 Supplemental Tables 4 and 5).",
    regions = "Multinational; clinical trials published 2005-2018 (PubMed-Medline, Citeline Trialtrove, ASCO and ESMO abstracts).",
    notes = "Study-level (aggregate) meta-analytic data: 153 treatment cohorts from 80 trials (21,305 patients) in total, of which 147 cohorts contributed grade 3/4 trAE rates to this model (Shulgin 2020 Supplemental Table 1). n_subjects / n_studies are the whole-analysis totals. sex_female_pct is 100 minus the mean cohort percentage of males (61.2%, 125 cohorts; Supplemental Table 4). Meta-regression fitted with the R package metafor on logit-transformed rates."
  )

  ini({
    # Shulgin 2020 Supplemental Table 1, block 'logit(trAE Grade 3/4)', row
    # 'beta0 + beta1 * dose-CTLA-4 + beta3 * dose-CTLA-4 * factor-PD-1'
    # (Model 4, N = 147 cohorts, AIC 300.2). Model 4 is the selected model
    # (main text Results, 'Meta-regression analyses').
    logit_ref <- -1.44
    label("Logit of the grade 3/4 trAE rate at zero anti-CTLA-4 exposure (unitless logit)") # Supplemental Table 1 Model 4, beta0 = -1.44 (SE 0.06, p<0.001)

    e_cnorm_ctla4_logit <- 0.0015
    label("Log-odds per unit normalized anti-CTLA-4 exposure (CAV/IC50), no anti-PD-1 antibody (unitless)") # Supplemental Table 1 Model 4, beta1 = 0.0015 (SE 0.0003, p<0.001); main text Results '0.0015 (95% CI, 0.001-0.002)'

    e_cnorm_ctla4_antipd1_logit <- 0.011
    label("Additional log-odds per unit normalized anti-CTLA-4 exposure when an anti-PD-1 antibody is co-administered (unitless)") # Supplemental Table 1 Model 4, beta3 = 0.011 (SE 0.0015, p<0.001); main text Results gives the combined slope beta1 + beta3 as '0.0124 (95% CI, 0.0095-0.0153)'

    # The source is a random-effects meta-regression of logit rates; the
    # between-study variance was not reported. The tiny fixed additive
    # residual below exists only so rxode2 has an error model to attach to
    # the typical-value probability; it is NOT a published quantity.
    addSd_prob_trae_grade34 <- fixed(0.001)
    label("Placeholder additive residual SD on the typical-value cohort trAE proportion; no source residual") # not from source; see vignette Assumptions and deviations
  })

  model({
    # Normalized anti-CTLA-4 exposure (Supplemental Methods: C_av / IC50)
    cnorm_ctla4 <- CAV / IC50_CTLA4

    # Shulgin 2020 Methods, Model 4:
    # logit(P_AE) = beta0 + beta1 * Cnorm,CTLA4 + beta3 * Cnorm,CTLA4 * FactorPD1
    logit_trae_grade34 <- logit_ref +
      e_cnorm_ctla4_logit * cnorm_ctla4 +
      e_cnorm_ctla4_antipd1_logit * cnorm_ctla4 * TRT_ANTIPD1

    prob_trae_grade34 <- expit(logit_trae_grade34)

    prob_trae_grade34 ~ add(addSd_prob_trae_grade34)
  })
}
