Shulgin_2020_ici_trae34_covariate_mbma <- function() {
  description <- "MBMA. Study-level logit meta-regression of the grade 3/4 treatment-related adverse event (trAE) rate on normalized anti-CTLA-4 exposure for immune checkpoint inhibitor (ICI) therapy with the anti-PD-1 antibodies nivolumab / pembrolizumab and the anti-CTLA-4 antibodies ipilimumab / tremelimumab, alone or in combination extended with the two trial-level covariates retained by the paper's stepwise search (Shulgin 2020 Supplemental Table 3; 122 cohorts with a single line of therapy): standard chemotherapy added to an anti-PD-1 antibody raises the logit by 1.1, and first-line therapy steepens the anti-CTLA-4 exposure slope. This is the model behind the prospective predictions of main-text Figure 3 (anti-PD-1 + chemotherapy; first- vs later-line anti-CTLA-4; the untested anti-PD-1 + anti-CTLA-4 + chemotherapy triple combination). Exposure is the steady-state average concentration of the anti-CTLA-4 antibody divided by its in vitro IC50; PD-1 inhibitor exposure has no effect. See Shulgin_2020_ici_trae34_mbma for the covariate-free Model 4 fit to all 147 cohorts. This is a STUDY-LEVEL model: it predicts the expected grade 3/4 trAE proportion of a trial cohort, not an individual-patient risk, and it consumes no dose events (exposure arrives as covariate columns)."
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
      notes = "The paper's 'FactorPD1' / 'factor-PD-1': 'set to 1 if a PD-1 inhibitor drug was given as monotherapy or in combination, and to 0 otherwise' (Shulgin 2020 Methods). Enters through the interaction with normalized anti-CTLA-4 exposure (beta3) and through the chemotherapy interaction chemo * factor-PD-1 (beta5).",
      source_name = "FactorPD1 / factor-PD-1 (Shulgin 2020 Methods, Supplemental Tables 1-3)"
    ),
    CONMED_CHEMO = list(
      description = "Indicator that standard chemotherapy was given together with the immune checkpoint inhibitor(s) in the cohort.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (immune checkpoint inhibitor(s) without chemotherapy)",
      notes = "The paper's 'chemo' covariate (Shulgin 2020 Supplemental Tables 2-4; 9 anti-PD-1 + chemotherapy and 7 anti-CTLA-4 + chemotherapy cohorts). In the retained model it enters ONLY as the product chemo * factor-PD-1 (beta5), so chemotherapy added to anti-CTLA-4 monotherapy has no effect on the prediction; chemotherapy with an anti-PD-1 antibody (with or without anti-CTLA-4) adds 1.1 to the logit.",
      source_name = "chemo (Shulgin 2020 Supplemental Tables 2-4)"
    ),
    LINE_1L = list(
      description = "Indicator that the cohort received the immune checkpoint inhibitor regimen as first-line therapy.",
      units = "(binary)",
      type = "binary",
      reference_category = "1 in this model's parameterization (first-line); the source covariate is the opposite-polarity 'Line2+' indicator, formed here as (1 - LINE_1L)",
      notes = "The paper codes line of therapy as 'Line2+' (1 = second-line or later, 0 = first-line); this model forms Line2+ = 1 - LINE_1L inside model(), so the beta6 coefficient keeps the paper's sign. Cohorts enrolling a mix of lines ('All' in Supplemental Table 5) were excluded from the covariate analysis (122 of 147 cohorts retained); the model is defined only for LINE_1L in {0, 1}. Supplemental Table 4: 51 first-line and 76 second-line-or-later cohorts.",
      source_name = "Line2+ (Shulgin 2020 Supplemental Tables 2 and 3); 'Line' in Supplemental Table 5"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 21305L,
    n_studies = 80L,
    n_cohorts = 122L,
    age_range = "Cohort median ages 42-74 years (mean of cohort medians 60.6, SD 4.76; 128 cohorts reporting; Shulgin 2020 Supplemental Table 4).",
    sex_female_pct = 38.8,
    disease_state = "Advanced solid tumours treated with PD-1 and/or CTLA-4 immune checkpoint inhibitors -- predominantly melanoma and non-small cell lung cancer, plus renal cell, small-cell lung, head and neck, gastric, urothelial, prostate, mesothelioma, sarcoma and other cancers (Shulgin 2020 Supplemental Table 5). Cancer type was screened and not retained.",
    dose_range = "Nivolumab 0.1-10 mg/kg Q2W/Q3W or 240 mg Q2W; pembrolizumab 2-10 mg/kg Q2W/Q3W or 200 mg Q3W; ipilimumab 0.3-10 mg/kg Q3W/Q4W as monotherapy and 1 or 3 mg/kg Q3W/Q6W/Q12W in combination with an anti-PD-1 antibody; tremelimumab 10 mg/kg Q4W or 15 mg/kg Q12W; alone, as PD-1 + CTLA-4 combination (27 cohorts), or with standard chemotherapy (Shulgin 2020 Supplemental Tables 4 and 5).",
    regions = "Multinational; clinical trials published 2005-2018 (PubMed-Medline, Citeline Trialtrove, ASCO and ESMO abstracts).",
    notes = "Study-level (aggregate) meta-analytic data: 153 treatment cohorts from 80 trials (21,305 patients) in total, of which 122 cohorts with a single reported line of therapy contributed grade 3/4 trAE rates to this model (Shulgin 2020 Supplemental Table 3). n_subjects / n_studies are the whole-analysis totals. sex_female_pct is 100 minus the mean cohort percentage of males (61.2%, 125 cohorts; Supplemental Table 4). Meta-regression fitted with the R package metafor on logit-transformed rates."
  )

  ini({
    # Shulgin 2020 Supplemental Table 3 ('Selection of optimal predictive
    # regression models for grade 3/4 trAE'), row
    # 'beta0 + beta1 * dose-CTLA-4 + beta3 * dose-CTLA-4 * factor-PD-1
    #  + beta5 * chemo * factor-PD-1 + beta6 * dose-CTLA-4 * (Line2+)'
    # (N = 122 cohorts, AIC 207.8). The table does not flag its selected row.
    # This row is the one with exactly the 'two additional covariates' of the
    # main text ('Incorporation of these two additional covariates into the
    # final model did not change the functional form'), and it is the only
    # Table 3 row that reproduces both the main-text 17.7% (anti-PD-1
    # monotherapy: expit(-1.53) = 17.8%) and 39.0% (anti-PD-1 + chemotherapy:
    # expit(-1.53 + 1.1) = 39.4%) predictions and the curves of Figure 3.
    logit_ref <- -1.53
    label("Logit of the grade 3/4 trAE rate at zero anti-CTLA-4 exposure without chemotherapy (unitless logit)") # Supplemental Table 3 row 7, beta0 = -1.53 (SE 0.05, p<0.001)

    e_cnorm_ctla4_logit <- 0.004
    label("Log-odds per unit normalized anti-CTLA-4 exposure (CAV/IC50), first-line, no anti-PD-1 antibody (unitless)") # Supplemental Table 3 row 7, beta1 = 0.004 (SE 0.0003, p<0.001)

    e_cnorm_ctla4_antipd1_logit <- 0.01
    label("Additional log-odds per unit normalized anti-CTLA-4 exposure when an anti-PD-1 antibody is co-administered (unitless)") # Supplemental Table 3 row 7, beta3 = 0.01 (SE 0.0014, p<0.001)

    e_chemo_antipd1_logit <- 1.1
    label("Additional logit when standard chemotherapy is given with an anti-PD-1 antibody (unitless)") # Supplemental Table 3 row 7, beta5 (chemo * factor-PD-1) = 1.1 (SE 0.22, p<0.001)

    e_cnorm_ctla4_line2_logit <- -0.003
    label("Change in log-odds per unit normalized anti-CTLA-4 exposure for second-line-or-later therapy (unitless)") # Supplemental Table 3 row 7, beta6 (dose-CTLA-4 * Line2+) = -0.003 (SE 0.0004, p<0.001)

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

    # Source covariate 'Line2+' is the complement of the canonical LINE_1L
    line2plus <- 1 - LINE_1L

    # Shulgin 2020 Supplemental Table 3 row 7:
    # logit(P) = beta0 + beta1 * C + beta3 * C * FactorPD1
    #            + beta5 * chemo * FactorPD1 + beta6 * C * Line2+
    logit_trae_grade34 <- logit_ref +
      e_cnorm_ctla4_logit * cnorm_ctla4 +
      e_cnorm_ctla4_antipd1_logit * cnorm_ctla4 * TRT_ANTIPD1 +
      e_chemo_antipd1_logit * CONMED_CHEMO * TRT_ANTIPD1 +
      e_cnorm_ctla4_line2_logit * cnorm_ctla4 * line2plus

    prob_trae_grade34 <- expit(logit_trae_grade34)

    prob_trae_grade34 ~ add(addSd_prob_trae_grade34)
  })
}
