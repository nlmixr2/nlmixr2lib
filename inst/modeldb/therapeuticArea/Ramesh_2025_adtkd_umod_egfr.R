Ramesh_2025_adtkd_umod_egfr <- function() {
  description <- "Natural-history disease-progression model of estimated glomerular filtration rate (eGFR, CKD-EPI creatinine) versus age in adults with autosomal dominant tubulointerstitial kidney disease caused by UMOD mutations (ADTKD-UMOD; Ramesh 2025). eGFR declines from an extrapolated early-adulthood value along an asymmetric sigmoid Imax function of age whose Hill coefficient differs before (shallow) and after (steep) the inflection age; the maximum fractional decrease is fixed. eGFR and age at the first clinic visit act as power covariates on the inflection age, the post-inflection Hill coefficient and the early-adulthood eGFR. No drug input; the model time is the subject's age in years."
  reference <- paste(
    "Ramesh SS, Rogge M, Kidd KO, Williams AH, Yoon DY, Roignot J, Blakeslee K,",
    "Bleyer AJ, Kim S. (2025). Quantifying clinical and genetic factors influencing",
    "rate and severity of autosomal dominant tubulointerstitial kidney disease",
    "progression. J Pharmacokinet Pharmacodyn 52(4):41.",
    "doi:10.1007/s10928-025-09989-0.",
    "Parameter estimates: Table 2A (ADTKD-UMOD column); structure and covariate",
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
      notes = "CKD-EPI creatinine equation (Levey 2009) from serum creatinine, sex and race (Methods 'Dependent and time variables'). Enters as the power term (CRCL_BASE / 43.4044)^beta, 43.4044 being the observation-weighted mean printed in Supplementary S2 ('logtBL_EGFR = log(BL_EGFR/43.4044)'). Effects on the post-inflection Hill coefficient (beta -0.68) and on the early-adulthood eGFR e0 (beta 1). The time-varying eGFR is the model OUTPUT, not a covariate.",
      source_name = "BL_EGFR (eGFR_FCV in the paper)"
    ),
    AGE_BASE = list(
      description = "Age at the subject's first clinic visit (AFCV), time-fixed",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters as the power term (AGE_BASE / 35.2138)^beta, 35.2138 being the observation-weighted mean printed in Supplementary S2 ('logtBL_VA = log(BL_VA/35.2138)'). Effects on the inflection age t50max (beta 0.93) and the post-inflection Hill coefficient (beta 1.14). Distinct from the subject's current age, which is the model time.",
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
    n_subjects = 371L,
    n_studies = 1L,
    n_observations = "UMOD and MUC1 together: 6447 visit records from 604 patients (Results 'Data summary'); models fitted to a random 4:1 training split, the remaining fifth held out for validation.",
    age_range = "18-80 years at assessment (whole data set); UMOD median 44.26 years (95% interval 22.48-63.10)",
    age_median = "44.26 years at assessment; 36.77 years at the first clinic visit (95% interval 18.34-58.93)",
    sex_female_pct = 58,
    race_ethnicity = "92% 'American White' and 89% not Hispanic or Latinx across the whole data set; race and ethnicity were therefore not tested as covariates.",
    disease_state = "Autosomal dominant tubulointerstitial kidney disease with a pathogenic UMOD (uromodulin) variant. eGFR (CKD-EPI) median 31.86 mL/min/1.73 m^2 over all visits; eGFR at the first clinic visit median 45.35 (95% interval 14.56-94.83). 105 of 371 reached kidney failure (eGFR <= 10, dialysis).",
    dose_range = "Not applicable -- natural-history model, no drug.",
    regions = "United States (Wake Forest University School of Medicine inherited-kidney-disease registry, referrals since 1996)",
    notes = "Table 1 ADTKD-UMOD column: 215 female / 156 male; promoter SNP rs4293393 phase A 359, phase G 12; mutation families p.H177-R185del 77, p.V93-G97delinsAASC 34, p.R178P 21, p.C106F 50, other 189. Covariates screened: baseline eGFR, baseline age, rs4293393 phase, sex and the regrouped common mutations (Results). eGFR values reported below 10 mL/min/1.73 m^2 were handled as left-censored (M4 method)."
  )

  ini({
    # Final ADTKD-UMOD eGFR model, Ramesh 2025 Table 2A ('Final estimates
    # (%RSE)' column). Monolix 2024R1, Supplementary S2 Model Code 1. The
    # Supplementary S2 <PARAMETER> block holds the INITIAL values of that run
    # (e.g. IC50_pop 38.95 vs Table 2A 39.76); only its FIXED entry is a final
    # value.
    le0 <- log(76.69); label("eGFR extrapolated to age 0, read as the early-adulthood (about 18 years) eGFR, DP0 (mL/min/1.73 m^2)") # Table 2A TV_DP0 = 76.69 (RSE 1.48%)
    logitimax <- fixed(log(0.999 / (1 - 0.999))); label("Maximum fractional decrease in eGFR, DPmax, on the logit scale (unitless)") # Supplementary S2 'Imax_pop = {value=0.999, method=FIXED}' (logitNormal, no variability); Table 2A and Results print 0.99
    lt50max <- log(39.76); label("Age at half of the maximum eGFR decrease, DPT50 (years)") # Table 2A TV_DPT50 = 39.76 (RSE 2.06%)
    lhill <- log(1.12); label("Hill coefficient of the age sigmoid before DPT50, gamma1 (unitless)") # Table 2A TV_gamma1 = 1.12 (RSE 15.2%)
    lhill_late <- log(6.34); label("Hill coefficient of the age sigmoid after DPT50, gamma2 (unitless)") # Table 2A TV_gamma2 = 6.34 (RSE 7.76%)

    # Power covariate effects (Supplementary S1 Equation S2), centred on the
    # observation-weighted means printed in Supplementary S2.
    e_age_base_t50max <- 0.93; label("Power exponent of age at first clinic visit on DPT50 (unitless)") # Table 2A beta_AFCV_DPT50 = 0.93 (RSE 2.05%)
    e_crcl_base_hill_late <- -0.68; label("Power exponent of eGFR at first clinic visit on gamma2 (unitless)") # Table 2A beta_eGFR_FCV_gamma2 = -0.68 (RSE 28.9%)
    e_age_base_hill_late <- 1.14; label("Power exponent of age at first clinic visit on gamma2 (unitless)") # Table 2A beta_AFCV_gamma2 = 1.14 (RSE 20.29%)
    e_crcl_base_e0 <- 1; label("Power exponent of eGFR at first clinic visit on DP0 (unitless)") # Table 2A beta_eGFR_FCV_DP0 = 1 (RSE 1.38%)

    # Log-normal IIV. Table 2A prints omega as the SD of the random effect
    # with the %CV in brackets; the variances are log(1 + CV^2), which
    # carries more significant digits than the printed SD (0.8 for gamma2,
    # whose 98.93% CV implies 0.826).
    etalt50max ~ 0.0057745 # Table 2A omega_DPT50 = 0.07 [7.61% CV]
    etalhill_late ~ 0.68245 # Table 2A omega_gamma2 = 0.8 [98.93% CV]
    etalhill + etale0 ~ c(1.5610, -0.14201, 0.016682) # Table 2A omega_gamma1 = 1.25 [194% CV], omega_DP0 = 0.13 [12.97% CV], correlation gamma1-DP0 = -0.88 (RSE 8.37%)

    # Monolix combined1 error model: DV = E + (a + b * E) * e.
    addSd <- 0.33; label("Additive residual SD (mL/min/1.73 m^2)") # Table 2A 'additive' = 0.33 (RSE 27.5%)
    propSd <- 0.11; label("Proportional residual SD (fraction)") # Table 2A 'proportional' = 0.11 (RSE 3.19%)
  })
  model({
    # The model time is the subject's age in years (Monolix regressor
    # TIME = age at the visit, Supplementary S2).
    age <- t

    e0 <- exp(le0 + etale0) * (CRCL_BASE / 43.4044)^e_crcl_base_e0
    imax <- expit(logitimax)
    t50max <- exp(lt50max + etalt50max) * (AGE_BASE / 35.2138)^e_age_base_t50max
    hill <- exp(lhill + etalhill)
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
    egfr ~ add(addSd) + prop(propSd) + combined1()
  })
}
