Zhang_2022_ici_irae_mbma <- function() {
  description <- "MBMA. Study-level logit meta-regression of the any-grade immune-related adverse event (irAE) rate in non-small cell lung cancer (NSCLC) cohorts treated with immune checkpoint inhibitors (ICIs): the anti-PD-1 antibodies nivolumab / pembrolizumab, the anti-PD-L1 antibodies atezolizumab / durvalumab / avelumab and the anti-CTLA-4 antibodies ipilimumab / tremelimumab, alone, in ICI combinations, or with chemotherapy or targeted therapy (Zhang 2022 final model, Equation 4 and Table 3; 126 cohorts). Anti-PD-(L)1 exposure has no effect: an anti-PD-L1 antibody lowers the logit by 0.35 relative to the anti-PD-1 reference, the logit rises 0.0013 per unit of normalized anti-CTLA-4 exposure (steady-state average concentration divided by the in vitro IC50), second-line-or-later therapy lowers it by 0.48, and combination with chemotherapy or targeted therapy raises it by 0.91. See Zhang_2022_ici_irae_grade3_mbma for the grade >= 3 endpoint. This is a STUDY-LEVEL model: it predicts the expected any-grade irAE proportion of a trial cohort, not an individual-patient risk, and it consumes no dose events (exposure arrives as covariate columns)."
  reference <- paste(
    "Zhang R, Kong D, Chen R, Guo Y, Jian W, Han M, Zhou T.",
    "A model-based meta-analysis of immune-related adverse events during immune",
    "checkpoint inhibitors treatment for NSCLC.",
    "CPT Pharmacometrics Syst Pharmacol. 2022;11(8):1135-1146. doi:10.1002/psp4.12834.",
    sep = " "
  )
  vignette <- "Zhang_2022_ici_irae_nsclc_mbma"

  # Algebraic study-level MBMA: no rxode2 dose events are consumed and the
  # output is a probability, not a drug concentration (same device as
  # Shulgin_2020_ici_trae34_covariate_mbma).
  units <- list(
    time = "n/a (static study-level meta-regression; no time dimension)",
    dosing = "n/a (no dose events; anti-CTLA-4 exposure enters as the covariates CAV and IC50_CTLA4)",
    concentration = "prob_irae (expected proportion of a cohort with an any-grade immune-related adverse event, 0-1)"
  )

  covariateData <- list(
    CAV = list(
      description = "Cohort-level steady-state average plasma concentration of the anti-CTLA-4 antibody (ipilimumab or tremelimumab) at the cohort's dosing regimen; 0 when the regimen contains no anti-CTLA-4 antibody.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Zhang 2022 Methods ('ICI exposure normalization'): Cav was simulated in NONMEM 7.4 from the mean parameters of published two-compartment linear popPK models (Supplementary Table S3: ipilimumab CL 0.360 L/day, Vc 4.15 L, Q 0.986 L/day, Vp 3.11 L; tremelimumab CL 0.262 L/day, Vc 3.72 L, Q 0.413 L/day, Vp 3.31 L) and divided by IC50_CTLA4 inside model() to give the paper's normalized exposure C_CTLA-4. The paper does not state the body weight or the averaging window. The normalized exposures the coefficients were fitted on are listed per cohort in the authors' supplementary dataset (column CCTLA4): ipilimumab 1 mg/kg Q12W 8, 1 mg/kg Q6W 17.5, 1 mg/kg Q3W 36, 3 mg/kg Q3W 109, 10 mg/kg Q3W 375; tremelimumab 1 mg/kg Q4W 115.61, 3 mg/kg Q4W 346.83, 10 mg/kg Q4W 1156.1. For ipilimumab these are 0.69-0.81 times Dose / (CL * tau) at 70 kg, so supply CAV = CCTLA4 * IC50_CTLA4 for these regimens to stay on the fitted axis (see the vignette). Must be in the same unit as IC50_CTLA4.",
      source_name = "C_av (Zhang 2022 Methods); 'C_CTLA-4' (Equations 1-5) / 'CCTLA4' (supplementary dataset) after normalization"
    ),
    IC50_CTLA4 = list(
      description = "In vitro IC50 of the anti-CTLA-4 antibody, used to normalize CAV.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Zhang 2022 Supplementary Table S3 prints ipilimumab 200 ng/mL and tremelimumab 95.15 ng/mL. Supplied here in ng/mL (the register's unit for this column is nM; only the ratio CAV / IC50_CTLA4 enters the model, so both columns must share one unit). 200 ng/mL is about 1.35 nM for a 148 kDa IgG, the same value Shulgin 2020 fitted on. Must be strictly positive in every row, including rows whose regimen has no anti-CTLA-4 antibody (where CAV = 0); a zero or missing IC50 propagates NaN into the prediction.",
      source_name = "IC50 (Zhang 2022 Supplementary Table S3)"
    ),
    TRT_ANTIPDL1 = list(
      description = "Indicator that the cohort's regimen includes an anti-PD-L1 antibody (atezolizumab, durvalumab or avelumab), as monotherapy or in combination.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no anti-PD-L1 antibody in the regimen -- the anti-PD-1 reference, or anti-CTLA-4 without an anti-PD-L1 antibody)",
      notes = "The paper's 'Factor_PD-L1': 'set to 1 if a PD-L1 inhibitor was given as monotherapy or in combination with a CTLA-4 inhibitor, and to 0 when not given' (Zhang 2022 Methods; column FactorPDL1 of the supplementary dataset). Shifts the intercept by beta1 = -0.3484.",
      source_name = "Factor_PD-L1 (Zhang 2022 Equations 2-5); FactorPDL1 (supplementary dataset)"
    ),
    LINE_1L = list(
      description = "Indicator that the cohort received the immune checkpoint inhibitor regimen as first-line therapy.",
      units = "(binary)",
      type = "binary",
      reference_category = "1 in this model's parameterization (first-line); the source covariate is the opposite-polarity 'line2+' indicator, formed here as (1 - LINE_1L)",
      notes = "The paper codes line of therapy as 'line2+' (1 = second-line or later, 0 = first-line); model() forms line2+ = 1 - LINE_1L so beta3 keeps the paper's sign. Cohorts that mixed first- and later-line patients were 'rounded off to 0 or 1 based on the percentage of patients receiving first-line versus second-line or later therapy', and the 12.4% of cohorts with missing line were imputed (Zhang 2022 Results and Supplementary Table S6), so the model is defined only for LINE_1L in {0, 1}.",
      source_name = "line2+ (Zhang 2022 Equations 4-5); Line2 (supplementary dataset)"
    ),
    CONMED_CHEMO = list(
      description = "Indicator that chemotherapy was given together with the immune checkpoint inhibitor(s) in the cohort.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no chemotherapy with the immune checkpoint inhibitor(s))",
      notes = "One half of the paper's single 'chemo/target' covariate (34 chemotherapy cohorts and one chemotherapy + targeted therapy cohort; Supplementary Tables S4 and S6). The paper fitted ONE coefficient for combination with chemotherapy OR targeted therapy, formed in model() as max(CONMED_CHEMO, CONMED_NONCHEMO_OTHER), so setting both to 1 gives the same prediction as setting either.",
      source_name = "chemo/target (Zhang 2022 Equations 4-5); Comb (supplementary dataset); 'chemo' in the Combination column of Supplementary Table S4"
    ),
    CONMED_NONCHEMO_OTHER = list(
      description = "Indicator that a targeted (non-chemotherapy) anticancer therapy was given together with the immune checkpoint inhibitor(s) in the cohort.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no targeted therapy with the immune checkpoint inhibitor(s))",
      notes = "Membership for this model: the 12 cohorts that Zhang 2022 Supplementary Table S4 labels 'target' in its Combination column (Supplementary Table S6); the agents are not named individually. The other half of the paper's single 'chemo/target' covariate, combined with CONMED_CHEMO as max(CONMED_CHEMO, CONMED_NONCHEMO_OTHER). Another ICI is NOT a member: ICI combinations are carried by TRT_ANTIPDL1 and CAV.",
      source_name = "chemo/target (Zhang 2022 Equations 4-5); Comb (supplementary dataset); 'target' in the Combination column of Supplementary Table S4"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 19322L,
    n_studies = 81L,
    n_cohorts = 126L,
    age_range = "Cohort median ages 50-72 years (median of cohort medians 64.3; 124 cohorts reporting; Zhang 2022 Supplementary Table S6).",
    sex_female_pct = 38.2,
    disease_state = "Non-small cell lung cancer treated with PD-1, PD-L1 and/or CTLA-4 immune checkpoint inhibitors, as monotherapy, as ICI combinations, or with chemotherapy or targeted therapy. Median across cohorts: 56.3% PD-L1 positive, 85.0% current or former smokers, 81.8% stage IV, 25.8% squamous histology (Supplementary Table S6); none of these was retained as a covariate.",
    dose_range = "Nivolumab 1-10 mg/kg Q2W/Q3W, 240 mg Q2W or 360/480 mg Q3W/Q4W; pembrolizumab 2-10 mg/kg Q2W/Q3W or 200 mg Q3W; atezolizumab 1200 mg Q3W; durvalumab 3-20 mg/kg Q2W/Q4W; avelumab 10 mg/kg Q2W; ipilimumab 3 or 10 mg/kg Q3W with chemotherapy or targeted therapy, or 1 mg/kg Q3W/Q6W/Q12W or 3 mg/kg Q3W with an anti-PD-1 antibody; tremelimumab 10 mg/kg Q4W alone or 1-10 mg/kg Q4W with durvalumab (Supplementary Table S4).",
    regions = "Multinational; clinical trials and real-world studies published up to April 30, 2021 (PubMed, Embase, Cochrane Library, ClinicalTrials.gov).",
    notes = "Study-level (aggregate) data from 129 treatment cohorts in 81 studies (19,322 patients): 53 RCTs, 18 dose-escalation trials, 21 single-arm trials, 20 nonrandomized trials and 17 real-world studies. 126 cohorts reported an any-grade irAE rate (73 anti-PD-1, 33 anti-PD-L1, 5 anti-CTLA-4, 10 anti-CTLA-4 + anti-PD-1, 5 anti-CTLA-4 + anti-PD-L1 monotherapy-class cohorts in Table 1). sex_female_pct is 100 minus the median cohort percentage of males (61.8%). Logit-transformed generalized linear mixed meta-regression with a random cohort intercept (R packages meta and metafor)."
  )

  ini({
    # Zhang 2022 Table 3 ('Parameter estimates of the final model'), any-grade
    # irAE block; Equation 4.
    logit_ref <- -1.2696
    label("Logit of the any-grade irAE rate for first-line anti-PD-1 monotherapy (unitless logit)") # Table 3, any grade, beta0 (intercept) = -1.2696 (SE 0.1255)

    e_antipdl1_logit <- -0.3484
    label("Additive logit shift for an anti-PD-L1 antibody in the regimen (unitless)") # Table 3, any grade, beta1 (on Factor_PD-L1) = -0.3484 (SE 0.1508; OR 0.7058)

    e_cnorm_ctla4_logit <- 0.0013
    label("Log-odds per unit normalized anti-CTLA-4 exposure CAV/IC50 (unitless)") # Table 3, any grade, beta2 (on C_CTLA-4) = 0.0013 (SE 0.0005; OR 1.0013)

    e_line2_logit <- -0.4757
    label("Additive logit shift for second-line-or-later therapy (unitless)") # Table 3, any grade, beta3 (on line2+) = -0.4757 (SE 0.1444; OR 0.6215)

    e_chemo_target_logit <- 0.9093
    label("Additive logit shift for combination with chemotherapy or targeted therapy (unitless)") # Table 3, any grade, beta4 (on chemo/target) = 0.9093 (OR 2.4826; printed SE 1.2124 is a misprint, see vignette Errata)

    # The source is a random-effects meta-regression; the between-cohort
    # variance was not reported. The tiny fixed additive residual below exists
    # only so rxode2 has an error model to attach to the typical-value
    # probability; it is NOT a published quantity.
    addSd_prob_irae <- fixed(0.001)
    label("Placeholder additive residual SD on the typical-value cohort irAE proportion; no source residual") # not from source; see vignette Assumptions and deviations
  })

  model({
    # Normalized anti-CTLA-4 exposure (Methods: C_av / IC50)
    cnorm_ctla4 <- CAV / IC50_CTLA4

    # Source covariate 'line2+' is the complement of the canonical LINE_1L
    line2plus <- 1 - LINE_1L

    # Source covariate 'chemo/target' is one indicator for chemotherapy OR
    # targeted therapy
    chemo_target <- max(CONMED_CHEMO, CONMED_NONCHEMO_OTHER)

    # Zhang 2022 Equation 4:
    # logit(P) = beta0 + beta1 * Factor_PD-L1 + beta2 * C_CTLA-4
    #            + beta3 * line2+ + beta4 * chemo/target
    logit_irae <- logit_ref +
      e_antipdl1_logit * TRT_ANTIPDL1 +
      e_cnorm_ctla4_logit * cnorm_ctla4 +
      e_line2_logit * line2plus +
      e_chemo_target_logit * chemo_target

    prob_irae <- expit(logit_irae)

    prob_irae ~ add(addSd_prob_irae)
  })
}
