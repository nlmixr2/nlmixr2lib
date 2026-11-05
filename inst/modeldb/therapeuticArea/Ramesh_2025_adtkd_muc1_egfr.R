Ramesh_2025_adtkd_muc1_egfr <- function() {
  description <- "Natural-history disease-progression model of estimated glomerular filtration rate (eGFR, CKD-EPI creatinine) versus age in adults with autosomal dominant tubulointerstitial kidney disease caused by MUC1 mutations (ADTKD-MUC1; Ramesh 2025). eGFR declines from an extrapolated early-adulthood value along an asymmetric sigmoid Imax function of age whose Hill coefficient differs before (shallow) and after (steep) the inflection age; the maximum fractional decrease is fixed. eGFR at the first clinic visit acts as a power covariate on the inflection age, both Hill coefficients and the early-adulthood eGFR, and age at the first clinic visit on the inflection age and the post-inflection Hill coefficient. No drug input; the model time is the subject's age in years."
  reference <- paste(
    "Ramesh SS, Rogge M, Kidd KO, Williams AH, Yoon DY, Roignot J, Blakeslee K,",
    "Bleyer AJ, Kim S. (2025). Quantifying clinical and genetic factors influencing",
    "rate and severity of autosomal dominant tubulointerstitial kidney disease",
    "progression. J Pharmacokinet Pharmacodyn 52(4):41.",
    "doi:10.1007/s10928-025-09989-0.",
    "Parameter estimates: Table 2A (ADTKD-MUC1 column); structure and covariate",
    "centring values: Supplementary Material S2 (Monolix project and model file,",
    "replaced by the publisher on 2026-04-10).",
    sep = " "
  )
  vignette <- "Ramesh_2025_adtkd"
  units <- list(
    time = "year (the subject's age; the model time IS age, not time since a dose)",
    dosing = "n/a (natural-history model, no drug input)",
    concentration = "mL/min/1.73 m^2 (eGFR)"
  )

  covariateData <- list(
    CRCL_BASE = list(
      description = "eGFR measured at the subject's first clinic visit (eGFR_FCV), time-fixed",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "CKD-EPI creatinine equation (Levey 2009) from serum creatinine, sex and race (Methods 'Dependent and time variables'). Enters as the power term (CRCL_BASE / 43.4044)^beta, 43.4044 being the observation-weighted mean printed in Supplementary S2 ('logtBL_EGFR = log(BL_EGFR/43.4044)'). Effects on the inflection age t50max (beta 0.04), the pre-inflection Hill coefficient (beta -0.67), the post-inflection Hill coefficient (beta -1.13) and the early-adulthood eGFR e0 (beta 0.96). Supplementary S2 prints the same centring value for the UMOD and MUC1 models. The time-varying eGFR is the model OUTPUT, not a covariate.",
      source_name = "BL_EGFR (eGFR_FCV in the paper)"
    ),
    AGE_BASE = list(
      description = "Age at the subject's first clinic visit (AFCV), time-fixed",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters as the power term (AGE_BASE / 35.2138)^beta, 35.2138 being the observation-weighted mean printed in Supplementary S2 ('logtBL_VA = log(BL_VA/35.2138)'). Effects on the inflection age t50max (beta 0.97) and the post-inflection Hill coefficient (beta 0.28). Distinct from the subject's current age, which is the model time.",
      source_name = "BL_VA (AFCV in the paper)"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Declared as a covariate in the Supplementary S2 Monolix project and screened by stepwise covariate modelling, but not retained in the final eGFR model (Results; Table 2A)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 233L,
    n_studies = 1L,
    n_observations = "UMOD and MUC1 together: 6447 visit records from 604 patients (Results 'Data summary'); models fitted to a random 4:1 training split, the remaining fifth held out for validation.",
    age_range = "18-80 years at assessment (whole data set); MUC1 median 45.65 years (95% interval 24.03-69.11)",
    age_median = "45.65 years at assessment; 37.68 years at the first clinic visit (95% interval 20.56-64.99)",
    sex_female_pct = 61,
    race_ethnicity = "92% 'American White' and 89% not Hispanic or Latinx across the whole data set; race and ethnicity were therefore not tested as covariates.",
    disease_state = "Autosomal dominant tubulointerstitial kidney disease with a pathogenic MUC1 (mucin-1) variant. eGFR (CKD-EPI) median 38.18 mL/min/1.73 m^2 over all visits; eGFR at the first clinic visit median 58.63 (95% interval 25.4-101.74). 145 of 233 reached kidney failure (eGFR <= 10, dialysis).",
    dose_range = "Not applicable -- natural-history model, no drug.",
    regions = "United States (Wake Forest University School of Medicine inherited-kidney-disease registry, referrals since 1996)",
    notes = "Table 1 ADTKD-MUC1 column: 142 female / 91 male; promoter SNP rs4072037 phase C 136, phase T 35, unknown 62. Covariates screened: baseline eGFR, baseline age, rs4072037 phase and sex; age at gout onset was also not significant (Results; Discussion). eGFR values reported below 10 mL/min/1.73 m^2 were handled as left-censored (M4 method)."
  )

  ini({
    # Final ADTKD-MUC1 eGFR model, Ramesh 2025 Table 2A ('Final estimates
    # (%RSE)' column). Monolix 2024R1, Supplementary S2 Model Code 2. The
    # Supplementary S2 <PARAMETER> block holds the INITIAL values of that run
    # (e.g. IC50_pop 38.41 vs Table 2A 38.36); only its FIXED entry is a final
    # value.
    le0 <- log(82.74); label("eGFR extrapolated to age 0, read as the early-adulthood (about 18 years) eGFR, DP0 (mL/min/1.73 m^2)") # Table 2A TV_DP0 = 82.74 (RSE 1.66%)
    logitimax <- fixed(log(0.999 / (1 - 0.999))); label("Maximum fractional decrease in eGFR, DPmax, on the logit scale (unitless)") # Supplementary S2 'Imax_pop = {value=0.999, method=FIXED}' (logitNormal, no variability); Table 2A and Results print 0.99
    lt50max <- log(38.36); label("Age at half of the maximum eGFR decrease, DPT50 (years)") # Table 2A TV_DPT50 = 38.36 (RSE 0.9%)
    lhill <- log(0.9); label("Hill coefficient of the age sigmoid before DPT50, gamma1 (unitless)") # Table 2A TV_gamma1 = 0.9 (RSE 22.6%)
    lhill_late <- log(10.23); label("Hill coefficient of the age sigmoid after DPT50, gamma2 (unitless)") # Table 2A TV_gamma2 = 10.23 (RSE 9.01%)

    # Power covariate effects (Supplementary S1 Equation S2), centred on the
    # observation-weighted means printed in Supplementary S2.
    e_crcl_base_t50max <- 0.04; label("Power exponent of eGFR at first clinic visit on DPT50 (unitless)") # Table 2A beta_eGFR_FCV_DPT50 = 0.04 (RSE 53.4%)
    e_age_base_t50max <- 0.97; label("Power exponent of age at first clinic visit on DPT50 (unitless)") # Table 2A beta_AFCV_DPT50 = 0.97 (RSE 2.41%)
    e_crcl_base_hill <- -0.67; label("Power exponent of eGFR at first clinic visit on gamma1 (unitless)") # Table 2A beta_eGFR_FCV_gamma1 = -0.67 (RSE 26.8%)
    e_crcl_base_hill_late <- -1.13; label("Power exponent of eGFR at first clinic visit on gamma2 (unitless)") # Table 2A beta_eGFR_FCV_gamma2 = -1.13 (RSE 16.4%)
    e_age_base_hill_late <- 0.28; label("Power exponent of age at first clinic visit on gamma2 (unitless)") # Table 2A beta_AFCV_gamma2 = 0.28 (RSE 83.1%)
    e_crcl_base_e0 <- 0.96; label("Power exponent of eGFR at first clinic visit on DP0 (unitless)") # Table 2A beta_eGFR_FCV_DP0 = 0.96 (RSE 4.34%)

    # Log-normal IIV. Table 2A prints omega as the SD of the random effect
    # with the %CV in brackets; the variances are log(1 + CV^2), which
    # carries more significant digits than the printed SD.
    etalt50max ~ 0.0058048 # Table 2A omega_DPT50 = 0.07 [7.63% CV]
    etalhill_late ~ 0.43671 # Table 2A omega_gamma2 = 0.66 [74% CV]
    etalhill + etale0 ~ c(0.045600, -0.012752, 0.0089845) # Table 2A omega_gamma1 = 0.21 [21.6% CV], omega_DP0 = 0.09 [9.5% CV], correlation gamma1-DP0 = -0.63 (RSE 105%)

    # Monolix proportional error model: DV = E + b * E * e.
    propSd <- 0.12; label("Proportional residual SD (fraction)") # Table 2A 'proportional' = 0.12 (RSE 1.70%)
  })
  model({
    # The model time is the subject's age in years (Monolix regressor
    # TIME = age at the visit, Supplementary S2).
    age <- t

    e0 <- exp(le0 + etale0) * (CRCL_BASE / 43.4044)^e_crcl_base_e0
    imax <- expit(logitimax)
    t50max <- exp(lt50max + etalt50max) *
      (CRCL_BASE / 43.4044)^e_crcl_base_t50max *
      (AGE_BASE / 35.2138)^e_age_base_t50max
    hill <- exp(lhill + etalhill) * (CRCL_BASE / 43.4044)^e_crcl_base_hill
    hill_late <- exp(lhill_late + etalhill_late) *
      (CRCL_BASE / 43.4044)^e_crcl_base_hill_late *
      (AGE_BASE / 35.2138)^e_age_base_hill_late

    # Asymmetric sigmoid of Supplementary S2 'Asymmetrical_Model_Code': the
    # Hill coefficient is gamma1 up to DPT50 and gamma2 beyond it. The two
    # branches meet at 0.5 when age = DPT50, so the curve is continuous.
    if (age <= t50max) {
      ratio <- (age / t50max)^hill
      frac <- ratio / (1 + ratio)
    } else {
      ratio <- (t50max / age)^hill_late
      frac <- 1 / (1 + ratio)
    }

    # Ramesh 2025 Equation 1: eGFR = DP0 * (1 - DPmax * f(age)).
    egfr <- e0 * (1 - imax * frac)
    egfr ~ prop(propSd)
  })
}
