Sheng_2026_sugemalimab <- function() {
  description <- "Two-compartment population PK model for intravenous sugemalimab (anti-PD-L1 IgG4) in adults with advanced solid tumours or lymphomas from six Phase I-III trials (n = 1002, 97.6% Asian), with total clearance the sum of a time-independent component and a mono-exponentially decaying time-dependent component, CL(t) = CL_inf + CL_T * exp(-kdes * t). Covariates: body weight, albumin and sex on CL_inf and Vc; sex and time-varying anti-drug-antibody titre on CL_T; lymphoma on kdes and Vc; albumin and NSCLC on Vp (Sheng 2026). This is the popPK model the paper uses to bridge the 1200 mg Q3W dose to heavier European NSCLC patients; it is a different fit from Wang 2024. Sister exposure-response models: Sheng_2026_sugemalimab_pfs_inv, Sheng_2026_sugemalimab_pfs_bicr, Sheng_2026_sugemalimab_os."
  reference <- paste(
    "Sheng Y, Yue Z, Xu F, Wang K, Wang J, Shi Q.",
    "Model-Informed Dosing Regimen of Sugemalimab for European Patients",
    "With Non-Small Cell Lung Cancer: Bridging From Asian Clinical Data.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15(3):e70220.",
    "doi:10.1002/psp4.70220.",
    "Parameter estimates from Supplementary Table S2; structure and",
    "covariate coding from the Supplementary Data S2 NONMEM control stream.",
    sep = " "
  )
  vignette <- "Sheng_2026_sugemalimab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    central = list(analyte = "sugemalimab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "sugemalimab", units = "mg", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL_inf (exponent 0.770) and Vc (exponent 0.446), normalised to 61.5 kg (control stream WTCL0 = (WT/61.5)**THETA(9), WTVC = (WT/61.5)**THETA(11)). Median 61 kg (range 36-124) in the modelling dataset (Table S1).",
      source_name = "WT"
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL_inf (-0.980), Vc (-0.332) and Vp (-1.40), normalised to 41.9 g/L (control stream ALBCL0 / ALBVC / ALBVP). Median 41.7 g/L (range 26.1-96) in the modelling dataset (Table S1).",
      source_name = "ALB"
    ),
    SEXF = list(
      description = "Biological sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Multiplicative factor THETA**SEX on CL_inf (0.904), CL_T (0.713) and Vc (0.858). The main-text equations apply each factor 'if Female', so the source SEX column is 1 = female, identical to canonical SEXF.",
      source_name = "SEX"
    ),
    ADA_TITER = list(
      description = "Time-varying anti-drug-antibody titre, last observation carried forward",
      units = "titre (dilution factor); 0.5 for ADA-negative",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters CL_T linearly on the natural-log scale: CL_T * (1 + 0.092 * (log(ADA_TITER) - log(0.5))). Zero-encoding convention: an ADA-negative record carries 0.5 (half the lowest dilution factor), NOT 0 and NOT 1, so the effect term is exactly 1 for ADA-negative subjects (control stream comment 'Baseline (negative ADA) was set to 0.5 (half dilution factor) in dataset'). Supplying 0 gives log(0) = -Inf and an invalid clearance. Source column TITERLCF (titre with LOCF). 112 of 1002 subjects (11.2%) were ADA-positive at some time (Table S1). The paper's dosing simulations omit ADA (all subjects negative), which it argues is conservative.",
      source_name = "TITERLCF"
    ),
    TUMTP_LYMPH = list(
      description = "Lymphoma tumour-type indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (solid tumour: NSCLC or other solid tumours)",
      notes = "Replaces the typical kdes (0.0465 vs 0.0188 1/day) and multiplies Vc by 0.896. Source DSSTAT == 1. The lymphoma cohort (n = 169) pools extranodal NK/T-cell lymphoma (study CS1001-201) and relapsed or refractory classical Hodgkin lymphoma (study CS1001-202).",
      source_name = "DSSTAT"
    ),
    TUMTP_NSCLC = list(
      description = "Non-small-cell lung cancer (stage III or IV) tumour-type indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (lymphoma or other solid tumour)",
      notes = "Multiplies Vp by 3.21. Source DSSTAT >= 2 (stage III NSCLC, n = 259, and stage IV NSCLC, n = 354). Other solid tumours (n = 220) share the reference with lymphoma for Vp and with NSCLC for Vc and kdes. TUMTP_LYMPH and TUMTP_NSCLC must not both be 1.",
      source_name = "DSSTAT"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1002L,
    n_studies = 6L,
    age_range = "18-78 years",
    age_median = "59 years",
    weight_range = "36-124 kg",
    weight_median = "61 kg",
    sex_female_pct = 22.55,
    race_ethnicity = c(Asian = 97.6, White = 2.0, Other = 0.4),
    disease_state = "Adults with advanced solid tumours or lymphomas: lymphoma 169 (16.9%), stage III NSCLC 259 (25.9%), stage IV NSCLC 354 (35.3%), other solid tumours 220 (22.0%).",
    dose_range = "Sugemalimab IV 3, 10, 20 or 40 mg/kg (25 subjects) or 1200 mg (977 subjects, 97.5%) every 3 weeks.",
    regions = "Predominantly China (97.6% Asian).",
    ada_status = "ADA-negative 890 (88.8%), ADA-positive 112 (11.2%).",
    notes = "Baseline demographics from Sheng 2026 Table S1. Two Phase I (CS1001-101 parts a and b, CS1001-102), two Phase II (CS1001-201 ENKTL, CS1001-202 cHL) and two Phase III (CS1001-301 / GEMSTONE-301 stage III NSCLC consolidation; CS1001-302 / GEMSTONE-302 stage IV NSCLC with platinum chemotherapy, n = 320) studies. 75 subjects (7.5%) weighed 80 kg or more. Fitted in NONMEM 7.5 by FOCE-I; time in the data set is in hours, parameters are reported per day."
  )

  ini({
    # Sheng 2026 Table S2. Typical values are for a male, 61.5 kg, albumin
    # 41.9 g/L, ADA-negative patient with a non-lymphoma tumour; Vp is for a
    # non-NSCLC tumour (the NSCLC factor is applied in model()).
    lcl_exp_inf <- log(0.149); label("Time-independent clearance CL0 (L/day)") # Table S2 CL0 = 0.149 L/day (RSE 1.4%)
    lcl_exp_component <- log(0.108); label("Initial time-dependent clearance CLT at t = 0 (L/day)") # Table S2 CLT = 0.108 L/day (RSE 4.4%)
    lcl_exp_kdes <- log(0.0188); label("Decay rate constant of time-dependent clearance kdes, non-lymphoma (1/day)") # Table S2 Kdes = 0.0188 1/day (RSE 6.8%)
    lcl_exp_kdes_lymph <- log(0.0465); label("Decay rate constant of time-dependent clearance kdes in lymphoma (1/day; replaces the non-lymphoma value)") # Table S2 'Effect of Lymphoma Disease on Kdes' = 0.0465 (RSE 13.7%); control stream IF(DSSTAT.EQ.1) TVKDES = THETA(13)
    lvc <- log(3.42); label("Central volume of distribution Vc (L)") # Table S2 Vc = 3.42 L (RSE 1.1%)
    lq <- log(0.373); label("Intercompartmental clearance Q (L/day)") # Table S2 Q = 0.373 L/day (RSE 8.0%)
    lvp <- log(0.573); label("Peripheral volume of distribution Vp, non-NSCLC (L)") # Table S2 Vp = 0.573 L (RSE 22.2%)

    # Continuous covariates: power exponents on (WT / 61.5) and (ALB / 41.9).
    e_wt_cl_exp_inf <- 0.770; label("Power exponent of body weight on CL0 (unitless)") # Table S2 Effect of Weight on CL0 = 0.770 (RSE 7.2%)
    e_alb_cl_exp_inf <- -0.980; label("Power exponent of albumin on CL0 (unitless)") # Table S2 Effect of Albumin on CL0 = -0.980 (RSE 9.0%)
    e_wt_vc <- 0.446; label("Power exponent of body weight on Vc (unitless)") # Table S2 Effect of Weight on Vc = 0.446 (RSE 9.7%)
    e_alb_vc <- -0.332; label("Power exponent of albumin on Vc (unitless)") # Table S2 Effect of Albumin on Vc = -0.332 (RSE 21.5%)
    e_alb_vp <- -1.40; label("Power exponent of albumin on Vp (unitless)") # Table S2 Effect of Albumin on Vp = -1.40 (RSE 25.2%)

    # Categorical covariates: the control stream applies THETA**SEX (or a
    # THETA multiplier for a disease group); carried here as the log of that
    # multiplier so that exp(e * indicator) reproduces it exactly.
    e_sexf_cl_exp_inf <- log(0.904); label("Log multiplier of CL0 for female sex (unitless)") # Table S2 Effect of Sex on CL0 = 0.904 (RSE 2.6%)
    e_sexf_cl_exp_component <- log(0.713); label("Log multiplier of CLT for female sex (unitless)") # Table S2 Effect of Sex on CLT = 0.713 (RSE 6.5%)
    e_sexf_vc <- log(0.858); label("Log multiplier of Vc for female sex (unitless)") # Table S2 Effect of Sex on Vc = 0.858 (RSE 2.1%)
    e_lymph_vc <- log(0.896); label("Log multiplier of Vc for lymphoma (unitless)") # Table S2 Effect of Lymphoma on Vc = 0.896 (RSE 2.1%)
    e_nsclc_vp <- log(3.21); label("Log multiplier of Vp for NSCLC (unitless)") # Table S2 Effect of NSCLC Disease on Vp = 3.21 (RSE 23.8%)

    # Time-varying ADA titre on CLT: 1 + theta * (log(TITER) - log(0.5)).
    e_adatiter_cl_exp_component <- 0.092; label("Slope of CLT on log ADA titre above the ADA-negative value 0.5 (per log unit)") # Table S2 Effect of Time-varying ADA Titer on CLT = 0.092 (RSE 15.6%)

    # Between-subject variability. The control stream estimates each SD and
    # each correlation directly as a THETA (Cholesky reparameterisation over
    # an OMEGA BLOCK(5) fixed to the identity), with ETA_i entering as
    # EXP(ETA_i). Table S2 reports those THETAs: the SD x 100 and the
    # correlation itself (the residual-error row, printed 16.9, is likewise
    # THETA(7) x 100). Variance = (SD / 100)^2 and covariance =
    # correlation * SD_i * SD_j. Correlations fixed to 0 in the control
    # stream ($THETA (0 FIXED) for CL0-CLT, CL0-Kdes, Vc-Kdes and all Vp
    # pairs) are carried as 0.
    #   CL0  SD 0.249 -> 0.062001       Vc   SD 0.186 -> 0.034596
    #   CLT  SD 0.463 -> 0.214369       Kdes SD 1.042 -> 1.085764
    #   cov(CL0, Vc)   = 0.365 * 0.249 * 0.186 = 0.016904610
    #   cov(Vc, CLT)   = 0.610 * 0.186 * 0.463 = 0.052531980
    #   cov(CLT, Kdes) = 0.448 * 0.463 * 1.042 = 0.216135808
    etalcl_exp_inf + etalvc + etalcl_exp_component + etalcl_exp_kdes ~ c(
      0.062001,
      0.016904610, 0.034596,
      0, 0.052531980, 0.214369,
      0, 0, 0.216135808, 1.085764
    ) # Table S2 BSV CL0 24.9, Vc 18.6, CLT 46.3, Kdes 104.2; correlations CL0-Vc 0.365, CLT-Vc 0.610, CLT-Kdes 0.448
    etalvp ~ 0.379456 # Table S2 BSV Vp 61.6 -> 0.616^2

    propSd <- 0.169; label("Proportional residual error (fraction)") # Table S2 Proportional Residual Unexplained Variability 16.9% (RSE 0.9%); control stream W = IPRED*THETA(7)
  })
  model({
    # Covariate-adjusted individual parameters (control stream $PK).
    cl_exp_inf <- exp(lcl_exp_inf + etalcl_exp_inf) *
      (WT / 61.5)^e_wt_cl_exp_inf *
      (ALB / 41.9)^e_alb_cl_exp_inf *
      exp(e_sexf_cl_exp_inf * SEXF)

    adatiter_eff <- 1 + e_adatiter_cl_exp_component * (log(ADA_TITER) - log(0.5))
    cl_exp_component <- exp(lcl_exp_component + etalcl_exp_component) *
      exp(e_sexf_cl_exp_component * SEXF) *
      adatiter_eff

    cl_exp_kdes <- exp(
      lcl_exp_kdes * (1 - TUMTP_LYMPH) +
        lcl_exp_kdes_lymph * TUMTP_LYMPH +
        etalcl_exp_kdes
    )

    vc <- exp(lvc + etalvc) *
      (WT / 61.5)^e_wt_vc *
      (ALB / 41.9)^e_alb_vc *
      exp(e_sexf_vc * SEXF) *
      exp(e_lymph_vc * TUMTP_LYMPH)

    q <- exp(lq)
    vp <- exp(lvp + etalvp) *
      (ALB / 41.9)^e_alb_vp *
      exp(e_nsclc_vp * TUMTP_NSCLC)

    # Control stream $DES: CL = CLT*EXP(-KDES*T/24) + CL0 with T in hours;
    # here t is days since the first dose.
    cl <- cl_exp_inf + cl_exp_component * exp(-cl_exp_kdes * t)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, volume in L: mg/L = ug/mL (control stream S1 = VC).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
