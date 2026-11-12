Sathe_2025b_sacituzumab_orr <- function() {
  description <- paste0(
    "Logistic-regression exposure-efficacy model for OBJECTIVE RESPONSE ",
    "(ORR: complete or partial response, by blinded independent central review) in adults with ",
    "metastatic triple-negative breast cancer treated with sacituzumab ",
    "govitecan 8 or 10 mg/kg on days 1 and 8 of 21-day cycles (Sathe ",
    "2025b, n = 277, pooled IMMU-132-01 and ASCENT). The probability of ",
    "objective response is expit(-6.96 + 0.343 * CAV), where CAV is the individual ",
    "average serum concentration of the intact antibody-drug conjugate ",
    "(SG) from the first dose to the time of the response assessment, ",
    "in ug/mL. CAVG of SG was the most significant of 12 exposure ",
    "metrics; no covariate was retained. There is no PK layer and no ",
    "ODE: exposure is supplied as data, and the source derived it from ",
    "the companion population PK model packaged as ",
    "Sathe_2024_sacituzumab. One of six exposure-response models in the ",
    "Sathe_2025b_sacituzumab_* family."
  )
  reference <- paste(
    "Sathe AG, Diderichsen PM, Fauchet F, Phan SC, Girish S, Othman AA.",
    "Exposure-Response Analyses of Sacituzumab Govitecan Efficacy and",
    "Safety in Patients With Metastatic Triple-Negative Breast Cancer.",
    "Clin Pharmacol Ther. 2025;117(2):570-578. doi:10.1002/cpt.3495.",
    "PMCID PMC11739744. Coefficients from Supplementary Table S7 (ORR rows).",
    sep = " "
  )
  vignette <- "Sathe_2025b_sacituzumab_exposure_response"

  units <- list(
    time = "n/a (static landmark exposure-response model; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "prob_orr_central (probability of objective response by blinded independent central review, 0-1)"
  )

  covariateData <- list(
    CAV = list(
      description = "Individual average serum concentration of INTACT sacituzumab govitecan (the antibody-drug conjugate, SG) from the first dose to the time of the response assessment. Supplied as data: this model has no PK layer.",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Sathe 2025b Eq. 1: CAVG_SG = integral from 0 to t_event of C_SG / t_event, with C_SG the population-PK-predicted SG concentration accounting for dose reductions and delays. Per Table S1, t_event for CR and ORR is the time of the first response, or 21 days after the last SG dose if the first response (or no response) was reported later than that. Reproduce the column from the SG output (Cc) of modellib('Sathe_2024_sacituzumab'). Enters LINEARLY and UNCENTRED. Figure 1 caption: mean CAVG_SG for the ORR endpoint was 12.9 ug/mL in the 8 mg/kg and 17.5 ug/mL in the 10 mg/kg starting-dose group.",
      source_name = "CAVG_SG"
    )
  )

  covariatesDataExcluded <- list(
    LDH = list(
      description = "Baseline serum lactate dehydrogenase",
      units = "IU/L",
      type = "continuous",
      notes = "Screened for ORR (Sathe 2025b Methods, Covariate selection) but not retained: 'After the inclusion of CAVG_SG in the model for CR and ORR, none of the evaluated covariates were found to be significant.' The same applies to baseline body weight, age, race, region, prior treatment, number of prior lines of therapy, ECOG performance status, baseline Trop-2 H score and UGT1A1 genotype. No point estimate is reported for any of them."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 277L,
    n_studies = 2L,
    age_range = "27-82 years (median 54)",
    weight_range = "37.2-132 kg (median 66.7)",
    sex_female_pct = 99.3,
    race_ethnicity = c(White = 80.5, Black = 11.2, Asian = 4.0, Other = 4.3),
    disease_state = "metastatic triple-negative breast cancer, previously treated",
    dose_range = "sacituzumab govitecan 8 mg/kg (n = 20), 10 mg/kg (n = 254) or 12 mg/kg (n = 3) starting dose IV on days 1 and 8 of 21-day cycles",
    regions = "North America 67.9%, Europe 32.1%",
    notes = "IMMU-132-01 (NCT01631552) mTNBC patients (n = 24) and ASCENT (NCT02574455) SG-arm patients (n = 253) with population-PK exposure estimates and at least one efficacy or safety outcome (Sathe 2025b Results and Table S3). 86 objective responders (31.0%; Table S4). Data cut-offs 1 March 2019 (IMMU-132-01) and 11 March 2020 (ASCENT)."
  )

  ini({
    # Sathe 2025b Supplementary Table S7, ORR rows. Fitted by maximum
    # likelihood in R; binomial logit link, no random effects.
    orr_intercept <- -6.96; label("Logit of the probability of objective response at CAVG_SG = 0 (unitless logit)") # Sathe 2025b Table S7, ORR Intercept: -6.96 [-8.72 to -5.21], p < 0.001
    e_cav_orr <- 0.343; label("Log-odds of objective response per 1 ug/mL increase in CAVG_SG (per ug/mL)") # Sathe 2025b Table S7, ORR Slope for CAVG_SG (per ug/mL): 0.343 [0.249-0.437], p < 0.001

    # The source likelihood is Bernoulli with no residual error; rxode2
    # needs an observation declaration, so the deterministic probability
    # carries a tiny placeholder additive residual (Oniki_2018_nafld_risk
    # convention).
    addSd_prob_orr_central <- fixed(0.001); label("Placeholder additive residual SD on the response probability; the source likelihood is Bernoulli") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_orr <- orr_intercept + e_cav_orr * CAV
    prob_orr_central <- expit(logit_orr)
    prob_orr_central ~ add(addSd_prob_orr_central)
  })
}
