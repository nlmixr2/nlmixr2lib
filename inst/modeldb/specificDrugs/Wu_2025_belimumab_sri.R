Wu_2025_belimumab_sri <- function() {
  description <- "Landmark logistic Emax exposure-response model for the probability of an SLE Responder Index (SRI) response at Week 52 in adult and paediatric patients with systemic lupus erythematosus given placebo or belimumab 1 or 10 mg/kg IV every 4 weeks or 200 mg SC weekly (Wu 2025, five studies, n = 2897). The logit is a placebo effect plus an Emax function of the individual steady-state average belimumab concentration. The placebo effect shifts with study BEL112341 (the SC study), Black African ancestry, Hawaiian or Native American ancestry and a baseline SELENA-SLEDAI score of 10 or more; the SELENA-SLEDAI threshold also adds to Emax. No PK layer and no ODE: the exposure enters as the CAV covariate column, derived in the source from the companion population PK model Wu_2025_belimumab. No random effect and no residual error are estimated (Bernoulli likelihood)."
  reference <- paste(
    "Wu J, Zhou X, Dimelow R, Marshall S. Inter-regional pharmacokinetics",
    "and exposure-response analyses of belimumab in patients with system",
    "lupus erythematosus. Br J Clin Pharmacol. 2025;91(2):374-385.",
    "doi:10.1111/bcp.16263.",
    "Final exposure-response estimates: Table 3.",
    sep = " "
  )
  vignette <- "Wu_2025_belimumab"
  units <- list(
    time = "n/a (static landmark exposure-response regression at Week 52; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "prob_sri (probability of an SRI response at Week 52, 0-1; also logit_sri)"
  )

  covariateData <- list(
    CAV = list(
      description = "Individual steady-state average serum belimumab concentration (Cavg_ss). Supplied as data; 0 for placebo.",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Source: Cavg_ss = F * Dose / (CL * tau) for SC and Dose[mg/kg] * BW / (CL * tau) for IV, with individual empirical-Bayes CL from the companion PK model (Supporting Information equations). Derive from modellib('Wu_2025_belimumab'). Table 2 prints the Cavg_ss unit as ng/mL, but the text, Table 3 (EC50 in ug/mL) and Table 4 use ug/mL, which matches the PK model's scale.",
      source_name = "Cavg_ss"
    ),
    STUDY_BEL112341 = list(
      description = "Participant of study BEL112341 (1) versus any other study in the exposure-response set (0)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (studies BEL110751, BEL110752, BEL113750 and BEL114055, all IV)",
      notes = "Additive shift of +0.404 on the placebo-effect logit (Table 3, footnote a). BEL112341 is the only SC (200 mg weekly) study in the set, so the effect cannot be separated from route or dosing frequency (Discussion). Set to 0 when simulating a new patient unless reproducing that trial.",
      source_name = "BEL112341 study effect"
    ),
    RACE_BLACK = list(
      description = "Black African ancestry (African-American) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Asian or Caucasian; Hawaiian or Native American ancestry has its own indicator)",
      notes = "Additive shift of -0.375 on the placebo-effect logit (Table 3, footnote b). Source race stratification RAC1 = B (Figure S3).",
      source_name = "RAC1 = B (Race_African-American)"
    ),
    RACE_NATIVE_AMERICAN_HAWAIIAN = list(
      description = "Hawaiian or Native American ancestry indicator (composite fitted as one level)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Asian or Caucasian; Black African ancestry has its own indicator)",
      notes = "Additive shift of +0.394 on the placebo-effect logit (Table 3, footnote c). Source race stratification RAC1 = O, 'Native American/Hawaiian' (Figure S3).",
      source_name = "RAC1 = O (Race_Hawaiian or Native American)"
    ),
    SCORE_SLEDAI = list(
      description = "Baseline SELENA-SLEDAI score",
      units = "(score, 0-105)",
      type = "continuous",
      reference_category = NULL,
      notes = "SELENA-SLEDAI variant. Enters only as the indicator SCORE_SLEDAI >= 10 (the paper's SE10): +0.533 on the placebo-effect logit and +0.503 on Emax (Table 3, footnote d). A threshold of 8 fitted equally well and 12 worse (Supporting Information sensitivity analysis).",
      source_name = "SSBL (SE10 = SSBL >= 10)"
    )
  )

  covariatesDataExcluded <- list(
    REGION_CHINA = list(
      description = "Mainland China participant indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (all other regions)",
      notes = "Tested on the placebo effect and Emax but not significant after accounting for exposure (Results)."
    ),
    AGE = list(
      description = "Age at baseline (also tested as a paediatric indicator, age <= 17 years)",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Not retained (Results)."
    ),
    WT = list(
      description = "Baseline body weight (also tested as a >= 74.4 kg quartile indicator)",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Not significant as a categorical or continuous covariate (Results)."
    ),
    IGG = list(
      description = "Baseline serum IgG",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the stepwise covariate search; not retained (Methods, Table 3)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 2897L,
    n_studies = 5L,
    n_observations = "2897 binary Week-52 SRI response records, one per patient (1800 belimumab, 1097 placebo); withdrawals before Week 52 counted as non-responders",
    age_range = "6-77 years (median 36.0); 92 paediatric patients aged 6-18 years",
    weight_range = "17.0-169.5 kg (median 62.0)",
    sex_female_pct = 94.1,
    race_ethnicity = "Table 2 strata: adult Caucasian 1298, Mainland Chinese 195, other North East Asian (Hong Kong, Taiwan, Japan, Korea) 354, others 958; paediatric 92 (48 Caucasian, 14 Asian, 4 Black, 26 Native American/Hawaiian).",
    disease_state = "Active systemic lupus erythematosus; baseline SELENA-SLEDAI median 10 (range 0-30); 92.6% on baseline steroid or immunosuppressant.",
    dose_range = "Placebo, belimumab 1 mg/kg or 10 mg/kg IV every 4 weeks, or 200 mg SC weekly, for 52 weeks.",
    regions = "Global phase 2/3 studies BEL110751, BEL110752, BEL112341 and BEL114055 (paediatric) plus the North East Asian study BEL113750.",
    notes = "Exposure-response analysis set (Table 2). Cavg_ss 2.5th-97.5th percentiles: 4.4-19.3 ug/mL at 1 mg/kg IV and 50.0-205.3 ug/mL at 10 mg/kg IV or 200 mg SC."
  )

  ini({
    # Table 3 final population exposure-response estimates, with its
    # footnote: Epbo = -0.789 + 0.404 * BEL112341 - 0.375 * African-American
    # + 0.394 * Hawaiian or Native American + 0.533 * SE10;
    # Emax = exp(-1.32) + 0.503 * SE10; EC50 = exp(2.07) = 7.92 ug/mL.
    e0 <- -0.789; label("Placebo-effect logit of the Week-52 SRI response for the reference patient (unitless logit)") # Table 3 'Placebo effect' = -0.789 (RSE 12.1%)
    e_study_bel112341_e0 <- 0.404; label("Additive shift on the placebo-effect logit for study BEL112341 (unitless logit)") # Table 3 'BEL112341 on placebo effect' = 0.404 (RSE 23.2%)
    e_race_black_e0 <- -0.375; label("Additive shift on the placebo-effect logit for Black African ancestry (unitless logit)") # Table 3 'Race_African-American on placebo effect' = -0.375 (RSE 40.0%)
    e_race_native_american_hawaiian_e0 <- 0.394; label("Additive shift on the placebo-effect logit for Hawaiian or Native American ancestry (unitless logit)") # Table 3 'Race_Hawaiian or Native American on placebo effect' = 0.394 (RSE 26.9%)
    e_sledai10_e0 <- 0.533; label("Additive shift on the placebo-effect logit for baseline SELENA-SLEDAI of 10 or more (unitless logit)") # Table 3 'SELENA-SLEDAI >=10 on placebo effect' = 0.533 (RSE 24.6%)
    lemax <- -1.32; label("Log of the maximum belimumab effect on the SRI logit for the reference patient (log unitless logit)") # Table 3 'Ln(Emax)' = -1.32 (RSE 45.8%)
    e_sledai10_emax <- 0.503; label("Additive shift on the maximum belimumab logit effect for baseline SELENA-SLEDAI of 10 or more (unitless logit)") # Table 3 'SELENA-SLEDAI >=10 on Emax' = 0.503 (RSE 37.2%)
    lec50 <- 2.07; label("Log of the steady-state average concentration giving half the maximum effect, EC50 (log ug/mL)") # Table 3 'Ln(EC50)' = 2.07 (RSE 48.0%); EC50 = 7.92 ug/mL

    # Bernoulli likelihood in the source; the fixed additive residual is a
    # placeholder so rxode2 has an error model, not a published quantity.
    addSd_prob_sri <- fixed(0.001); label("Placeholder additive residual SD on the SRI response probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    sledai10 <- 0
    if (SCORE_SLEDAI >= 10) sledai10 <- 1

    epbo <- e0 + e_study_bel112341_e0 * STUDY_BEL112341 +
      e_race_black_e0 * RACE_BLACK +
      e_race_native_american_hawaiian_e0 * RACE_NATIVE_AMERICAN_HAWAIIAN +
      e_sledai10_e0 * sledai10
    # Table 3 footnote: Emax = exp(lnEmax) + 0.503 * SE10 (additive on the
    # natural scale, not inside the exponential)
    emax <- exp(lemax) + e_sledai10_emax * sledai10
    ec50 <- exp(lec50)

    logit_sri <- epbo + emax * CAV / (CAV + ec50)
    prob_sri <- expit(logit_sri)

    prob_sri ~ add(addSd_prob_sri)
  })
}
