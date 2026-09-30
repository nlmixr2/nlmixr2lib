Leil_2021_rheumatoidArthritis_das28_mbma <- function() {
  description <- paste0(
    "MBMA. Longitudinal model-based meta-analysis of the change from ",
    "baseline in the 28-joint Disease Activity Score (DAS28, on the DAS28-CRP ",
    "scale) for seven approved rheumatoid arthritis drugs (abatacept, ",
    "adalimumab, certolizumab, etanercept, rituximab, tocilizumab, ",
    "tofacitinib) on a background of conventional synthetic DMARDs, fitted ",
    "in NONMEM 7.3 to 994 study-arm-mean records from 130 randomized trials ",
    "(197 arms, 27,355 patients). The arm-mean change from baseline is the ",
    "sum of a sigmoid Hill-in-time placebo (background-therapy) response, a ",
    "hyperbolic drug-specific Emax-in-time treatment effect and a linear ",
    "disease-progression term. Placebo Emax depends on baseline DAS28 and ",
    "disease duration, placebo ET50 on baseline DAS28, the placebo Hill ",
    "coefficient on a low-male-proportion indicator (male < 18.5 percent) ",
    "and the progression slope on the trial year. Between-trial and ",
    "between-arm random effects and the residual error are all scaled by ",
    "sqrt(100 / N_ARM). Drug arms are selected with per-drug arm indicators ",
    "(all zero = background-therapy placebo arm). Suitable simulation scope is ",
    "study-arm-mean DAS28 change-from-baseline trajectories, NOT individual ",
    "patients."
  )

  reference <- paste(
    "Leil TA, Lu Y, Bouillon-Pichault M, Wong R, Nowak M.",
    "Model-Based Meta-Analysis Compares DAS28 Rheumatoid Arthritis Treatment",
    "Effects and Suggests an Expedited Trial Design for Early Clinical",
    "Development. Clin Pharmacol Ther. 2021;109(2):517-527.",
    "doi:10.1002/cpt.2023.",
    "Parameter estimates are in Supplementary Table S3 (CPT-109-517-s001.docx);",
    "the structural model is main-text Equations 6-10 and 12.",
    sep = " "
  )

  vignette <- "Leil_2021_rheumatoidArthritis_das28_mbma"

  units <- list(
    time = "week",
    dosing = "none (no rxode2 dose events; the treatment arm enters through the per-drug arm-indicator covariates)",
    concentration = "DAS28/arm (das28cfb is the study-arm-mean change from baseline in DAS28-CRP units, negative = improvement; NOT a drug concentration. The slash satisfies checkModelConventions unit parsing.)"
  )

  covariateData <- list(
    ABATACEPT = list(
      description = "Study-arm treatment indicator: 1 = the arm received abatacept (10 mg/kg IV q4w or 125 mg SC qw) on background csDMARD, 0 = it did not.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (arm did not receive abatacept)",
      notes = "MBMA arm indicator; the seven drug indicators are mutually exclusive and all-zero marks a background-therapy (placebo) arm. Only the approved maintenance regimens were included (Leil 2021 Methods and Table 3), so the model carries no dose term.",
      source_name = "drug = abatacept (Leil 2021 Table 2)"
    ),
    ADALIMUMAB = list(
      description = "Study-arm treatment indicator: 1 = the arm received adalimumab (40 mg SC q2w) on background csDMARD, 0 = it did not.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (arm did not receive adalimumab)",
      notes = "MBMA arm indicator; mutually exclusive with the other six drug indicators.",
      source_name = "drug = adalimumab (Leil 2021 Table 2)"
    ),
    CERTOLIZUMAB = list(
      description = "Study-arm treatment indicator: 1 = the arm received certolizumab pegol (200 mg SC q2w or 400 mg SC q4w) on background csDMARD, 0 = it did not.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (arm did not receive certolizumab)",
      notes = "MBMA arm indicator; mutually exclusive with the other six drug indicators.",
      source_name = "drug = certolizumab (Leil 2021 Table 2)"
    ),
    ETANERCEPT = list(
      description = "Study-arm treatment indicator: 1 = the arm received etanercept (50 mg SC qw) on background csDMARD, 0 = it did not.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (arm did not receive etanercept)",
      notes = "MBMA arm indicator; mutually exclusive with the other six drug indicators.",
      source_name = "drug = etanercept (Leil 2021 Table 2)"
    ),
    RITUXIMAB = list(
      description = "Study-arm treatment indicator: 1 = the arm received rituximab (1,000 mg IV infusion x 2 per year) on background csDMARD, 0 = it did not.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (arm did not receive rituximab)",
      notes = "MBMA arm indicator; mutually exclusive with the other six drug indicators.",
      source_name = "drug = rituximab (Leil 2021 Table 2)"
    ),
    TOCILIZUMAB = list(
      description = "Study-arm treatment indicator: 1 = the arm received tocilizumab (8 mg/kg IV q4w) on background csDMARD, 0 = it did not.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (arm did not receive tocilizumab)",
      notes = "MBMA arm indicator; mutually exclusive with the other six drug indicators. No tocilizumab trial reported both DAS28-ESR and DAS28-CRP; its DAS28-ESR records were converted with Eq. 1 after the CRP-vs-ESR regression showed no MoA dependence (Leil 2021 Methods and Results).",
      source_name = "drug = tocilizumab (Leil 2021 Table 2)"
    ),
    TOFACITINIB = list(
      description = "Study-arm treatment indicator: 1 = the arm received tofacitinib (5 mg PO bid) on background csDMARD, 0 = it did not.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (arm did not receive tofacitinib)",
      notes = "MBMA arm indicator; mutually exclusive with the other six drug indicators.",
      source_name = "drug = tofacitinib (Leil 2021 Table 2)"
    ),
    N_ARM = list(
      description = "Number of patients in the study arm (the arm-mean sample size).",
      units = "participants",
      type = "count",
      reference_category = NULL,
      reference_value = 100,
      notes = "Every random effect and the residual error is multiplied by (N / 100)^-0.5 = sqrt(100 / N_ARM) (Leil 2021 Equations 7, 8 and 12). The paper scales the between-TRIAL random effects by the trial's mean arm size N_i and the between-arm and residual terms by the arm's own N_ij; this file uses N_ARM for all of them, which is exact when the arms of a trial are equal-sized. N_ARM = 100 recovers the unscaled variances.",
      source_name = "N_i,j (Leil 2021 Equations 5, 6 and 12); N_i (Equations 7 and 8)"
    ),
    SCORE_DAS28CRP = list(
      description = "Study-arm mean baseline DAS28 on the DAS28-CRP scale.",
      units = "(score)",
      type = "continuous",
      reference_category = NULL,
      reference_value = 6.2,
      notes = "Arm-level MEAN of the per-subject score. Enters placebo Emax additively as 1.60 * log(SCORE_DAS28CRP / 6.2) and placebo ET50 exponentially as exp(1.30 * log(SCORE_DAS28CRP / 6.2)) (Leil 2021 Equation 9, Table S3). Arms reporting DAS28-ESR were converted with DAS28-CRP = 0.899 * DAS28-ESR - 0.194 (Eq. 1, Results). The data-set median used for centring is not printed; 6.2 is the typical-trial baseline of Table 3 footnote b / Figure 3 (Methods: the simulated population was 'based on the median values across the trials').",
      source_name = "baseline DAS28 (Leil 2021 Table 2, Table S3)"
    ),
    T_DIAG_RA = list(
      description = "Study-arm mean time since rheumatoid arthritis diagnosis (disease duration).",
      units = "year",
      type = "continuous",
      reference_category = NULL,
      reference_value = 8.2,
      notes = "Enters placebo Emax additively as -0.133 * log(T_DIAG_RA / 8.2) (Leil 2021 Equation 9, Table S3). Centred on the 8.2-year typical-trial value of Table 3 footnote b; the data-set median is not printed.",
      source_name = "disease duration (years) (Leil 2021 Table 2, Table S3)"
    ),
    SEXF_PCT = list(
      description = "Percentage of female patients in the study arm (0-100).",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      reference_value = 81,
      notes = "The paper's covariate is the percentage of MALE participants, dichotomised at the 18.5 percent data-set median; the retained effect acts on the placebo Hill coefficient for arms with male < 18.5 percent. Derived in model() as lowMale = (100 - SEXF_PCT) < 18.5. The Table S3 label ('gamma_placebo ~ male participants < 18.5%') and the Table 3 typical trial (19 percent male, reproduced with the unmodified gamma_placebo) fix the indicator to 1 for the LOW-male category; the Methods parenthetical coding ('< 18.5% = 0; >= 18.5% = 1') is the reverse and is not used. Reference value 81 = 100 - 19 percent male (Table 3 footnote b).",
      source_name = "% male (Leil 2021 Table 2, Methods, Table S3); SEXF_PCT = 100 - % male"
    ),
    YEAR_PUB = list(
      description = "Calendar year of the trial; in this model the value is the year the trial was CONDUCTED, not the publication year.",
      units = "year",
      type = "continuous",
      reference_category = NULL,
      reference_value = 2013,
      notes = "Leil 2021 Methods describe the covariate as 'the year of conduct of the trial', so supply the trial-conduct year; if only the publication year is known, it lags the conduct year by a few years and the steep exponent makes that lag matter. Enters the disease-progression slope as exp(-397 * log(YEAR_PUB / 2013)) (Leil 2021 Equation 9, Table S3). The centring year is NOT printed; 2013 is the median publication year of the 130 trials listed in Supplementary Table S1 (range 2002-2017). The effect is extremely steep (a trial 5 years before the centre has a ~2.7-fold faster progression) and does not reproduce the Discussion's '~0.31 DAS28 units/year in 2000 vs ~0.0025 in 2016' under any centring year; see the vignette Errata.",
      source_name = "year of trial conduct (Leil 2021 Methods, Table S3 'Slope of DAS28 progression ~ trial year')"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Study-arm mean age. Screened but NOT retained.",
      units = "year",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on Emax, ET50 and gamma for drug and placebo arms (Leil 2021 Methods); not in the final model (Table S3). Table 3 typical trial: 53 years.",
      source_name = "Age (years) (Leil 2021 Table 2)"
    ),
    MTX_IR_PCT = list(
      description = "Percentage of patients in the arm who had failed methotrexate. Screened but NOT retained.",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "Dichotomised at < 100 percent and tested (Leil 2021 Methods); not in the final model (Table S3).",
      source_name = "percentage of participants who failed methotrexate (Leil 2021 Methods)"
    ),
    CONMED_CSDMARD_PCT = list(
      description = "Percentage of patients in the arm on background csDMARDs. Screened but NOT retained.",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "Dichotomised at < 50 percent and tested (Leil 2021 Methods); not in the final model (Table S3).",
      source_name = "percentage of participants on background csDMARDs (Leil 2021 Methods)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 27355L,
    n_studies = 130L,
    age_range = "trial-mean age 45-61 years (weighted mean of arm means 51-53 by drug; Leil 2021 Table 2)",
    weight_range = "not reported",
    sex_female_pct = 81,
    race_ethnicity = "not reported",
    disease_state = "adults with active moderate-to-severe rheumatoid arthritis on stable background csDMARDs (predominantly methotrexate); most arms were methotrexate inadequate responders; prior biologic exposure allowed",
    dose_range = "approved maintenance regimens only: abatacept 10 mg/kg IV q4w or 125 mg SC qw; adalimumab 40 mg SC q2w; certolizumab 200 mg SC q2w or 400 mg SC q4w; etanercept 50 mg SC qw; rituximab 1,000 mg IV x 2 per year; tocilizumab 8 mg/kg IV q4w; tofacitinib 5 mg PO bid",
    regions = "multinational (published trials 1994-2018, Quantify RA Clinical Outcomes Database v03/08/2018)",
    baseline = "trial-mean baseline DAS28-CRP 3.7-6.5 (weighted mean 4.9-5.9 by drug); disease duration 0.24-14 years; 0-45 percent male (Leil 2021 Table 2)",
    notes = "Study-arm-level MBMA: 994 arm-mean DAS28 change-from-baseline records from 197 arms (91 active, 106 background-therapy placebo) of 130 randomized controlled trials; per-drug trials/patients: abatacept 16/3478, adalimumab 15/2802, certolizumab 10/2654, etanercept 8/1747, rituximab 14/1336, tocilizumab 18/3724, tofacitinib 10/1611, placebo 106/10003 (Table 2). DAS28-ESR records were converted to DAS28-CRP with Eq. 1. The model simulates arm-mean trajectories, NOT individual patients; the paper's separate within-arm SD model (Eq. 11) is not part of this file."
  )

  ini({
    # ============================================================
    # Final model, Leil 2021 Equation 12 (Hill function, Eq. 6, plus
    # a linear progression term):
    #
    #   dDAS28_ij(t) = -[ Emax,pbo_i * t^g_i / (ET50,pbo_i^g_i + t^g_i)
    #                   + Emax,drug_ij * t^g_drug / (ET50,drug^g_drug + t^g_drug) ]
    #                 + PROG_i * t + (N_ij / 100)^-0.5 * eps_ij(t)
    #
    # All estimates: Leil 2021 Supplementary Table S3 ('Estimated
    # parameters and 95% confidence intervals from the model-based
    # meta-analysis'), Estimate column. Emax values are in DAS28
    # units (see the drug Emax note below), ET50 in weeks.
    # ============================================================

    # ---- Placebo (background-therapy) response -------------------
    emax_pbo <- 2.4
    label("Maximum placebo (background-therapy) decrease from baseline in DAS28 for the typical trial (DAS28 units). Normal (additive) inter-trial variability, so not log-transformed.") # Table S3 'Emax Placebo' = 2.4 (95% CI 2.12 to 2.68)
    let50_pbo <- log(26.7)
    label("Log time to half-maximal placebo response (log week)") # Table S3 'ET50 Placebo' = 26.7 weeks (95% CI 20.5 to 34.4)
    lhill_pbo <- log(0.568)
    label("Log Hill coefficient of the placebo response (unitless)") # Table S3 'gamma_placebo' = 0.568 (95% CI 0.485 to 0.655)

    # ---- Disease progression ---------------------------------------
    lslope <- log(0.103)
    label("Log slope of the linear DAS28 progression (worsening) over time (log DAS28 units/year)") # Table S3 'Slope of DAS28 progression [DAS28 units/year]' = 0.103 (95% CI 0.0599 to 0.155); the Results text rounds it as 0.105

    # ---- Drug-specific Emax (DAS28 units) --------------------------
    # Table S3 footnote a describes Emax as 'maximum reduction in DAS28
    # as a proportion of the baseline value', but the tabulated values
    # (placebo 2.4, tocilizumab 2.34) cannot be fractions of a ~6-unit
    # baseline, and the absolute-unit reading reproduces Table 3 exactly
    # (tocilizumab 2.34 * t / (4.24 + t) = 1.14, 1.73, 1.99, 2.15 at
    # 4, 12, 24, 48 weeks). Encoded in absolute DAS28 units.
    lemax_abatacept <- log(1.3)
    label("Log maximum abatacept treatment effect (placebo-corrected decrease in DAS28) (log DAS28 units)") # Table S3 'Emax Abatacept' = 1.3 (95% CI 1.11 to 1.52)
    lemax_adalimumab <- log(0.946)
    label("Log maximum adalimumab treatment effect (log DAS28 units)") # Table S3 'Emax Adalimumab' = 0.946 (95% CI 0.765 to 1.14)
    lemax_certolizumab <- log(1.24)
    label("Log maximum certolizumab treatment effect (log DAS28 units)") # Table S3 'Emax Certolizumab' = 1.24 (95% CI 1.01 to 1.49)
    lemax_etanercept <- log(1.3)
    label("Log maximum etanercept treatment effect (log DAS28 units)") # Table S3 'Emax Etanercept' = 1.3 (95% CI 1.02 to 1.6)
    lemax_rituximab <- log(1.57)
    label("Log maximum rituximab treatment effect (log DAS28 units)") # Table S3 'Emax Rituximab' = 1.57 (95% CI 1.21 to 1.97)
    lemax_tocilizumab <- log(2.34)
    label("Log maximum tocilizumab treatment effect (log DAS28 units)") # Table S3 'Emax Tocilizumab' = 2.34 (95% CI 2.14 to 2.55)
    lemax_tofacitinib <- log(1.09)
    label("Log maximum tofacitinib treatment effect (log DAS28 units)") # Table S3 'Emax Tofacitinib' = 1.09 (95% CI 0.86 to 1.31)

    # ---- Drug-specific ET50 (weeks) ---------------------------------
    let50_abatacept <- log(3.42)
    label("Log time to half-maximal abatacept effect (log week)") # Table S3 'ET50 Abatacept' = 3.42 (95% CI 2.19 to 5.16)
    let50_adalimumab <- log(2.4)
    label("Log time to half-maximal adalimumab effect (log week)") # Table S3 'ET50 Adalimumab' = 2.4 (95% CI 1.32 to 3.75)
    let50_certolizumab <- log(2.2)
    label("Log time to half-maximal certolizumab effect (log week)") # Table S3 'ET50 Certolizumab' = 2.2 (95% CI 1.33 to 3.31)
    let50_etanercept <- log(2.9)
    label("Log time to half-maximal etanercept effect (log week)") # Table S3 'ET50 Etanercept' = 2.9 (95% CI 0.857 to 5.85)
    let50_rituximab <- log(13.5)
    label("Log time to half-maximal rituximab effect (log week)") # Table S3 'ET50 Rituximab' = 13.5 (95% CI 9.96 to 18.5)
    let50_tocilizumab <- log(4.24)
    label("Log time to half-maximal tocilizumab effect (log week)") # Table S3 'ET50 Tocilizumab' = 4.24 (95% CI 3.39 to 5.11)
    let50_tofacitinib <- log(1.77)
    label("Log time to half-maximal tofacitinib effect (log week)") # Table S3 'ET50 Tofacitinib' = 1.77 (95% CI 1.12 to 2.43)

    # Drug Hill coefficient. Equations 6 and 12 carry a gamma_drug and
    # Table S3 reports its inter-trial SD, but no typical value is
    # tabulated among the fixed effects. gamma_drug = 1 reproduces the
    # Table 3 tocilizumab medians to all three printed digits at every
    # time point, so it is encoded as a typical value fixed at 1.
    lhill_drug <- fixed(log(1))
    label("Log Hill coefficient of the drug treatment effect (unitless); typical value 1, not tabulated") # not in Table S3 fixed effects; value 1 inferred from Table 3 (tocilizumab row reproduced exactly)

    # ---- Covariate effects (Leil 2021 Equations 9 and 10, Table S3) ---
    # Continuous covariates: COVEFF = theta * log(COV / COVmed); binary:
    # COVEFF = theta * IND. 'Incorporated in the model additively' --
    # on the parameter scale for the normally distributed placebo Emax,
    # in the exponent for the log-normally distributed ET50, gamma and
    # progression slope.
    e_das28_emax_pbo <- 1.60
    label("Additive effect of log(baseline DAS28 / 6.2) on placebo Emax (DAS28 units)") # Table S3 'Emax,placebo ~ baseline DAS28' = 1.60 (95% CI 1.07 to 2.15)
    e_tdiag_emax_pbo <- -0.133
    label("Additive effect of log(disease duration / 8.2 years) on placebo Emax (DAS28 units)") # Table S3 'Emax,placebo ~ disease duration' = -0.133 (95% CI -0.199 to -0.0679)
    e_das28_et50_pbo <- 1.30
    label("Exponent of (baseline DAS28 / 6.2) on placebo ET50 (unitless)") # Table S3 'ET50,placebo ~ baseline DAS28' = 1.30 (95% CI 0.867 to 1.61)
    e_lowmale_hill_pbo <- -0.176
    label("Log-scale shift in the placebo Hill coefficient for arms with fewer than 18.5 percent male patients (unitless)") # Table S3 'gamma_placebo ~ male participants < 18.5%' = -0.176 (95% CI -0.256 to -0.0941)
    e_year_slope <- -397
    label("Exponent of (trial year / 2013) on the progression slope (unitless)") # Table S3 'Slope of DAS28 progression ~ trial year' = -397 (95% CI -696 to -80.6)

    # ============================================================
    # Random effects. Table S3 reports SDs; ini() holds variances.
    # Each eta is multiplied by sqrt(100 / N_ARM) in model() (Eqs. 7,
    # 8). Footnote c marks placebo Emax as ADDITIVE (normal, Eq. 8);
    # footnote b marks the rest as PROPORTIONAL (log-normal, Eq. 7).
    # These are MBMA between-TRIAL (study-level) and between-ARM
    # effects, not between-subject variability.
    # ============================================================
    eta_study_emax_pbo ~ 0.685584 # Table S3 between-study 'Emax,placebo' SD = 0.828 (additive, footnote c); 0.828^2
    eta_study_hill_pbo ~ 0.380689 # Table S3 between-study 'gamma_placebo' SD = 0.617 (proportional, footnote b); 0.617^2
    eta_study_emax_drug ~ 0.073984 # Table S3 between-study 'Emax,drug' SD = 0.272 (proportional); 0.272^2
    eta_study_hill_drug ~ 0.308025 # Table S3 between-study 'gamma_drug' SD = 0.555 (proportional); 0.555^2
    eta_study_slope ~ 1 # Table S3 between-study 'Slope of DAS28 progression' SD = 1.00 (proportional); 1.00^2
    eta_arm_emax_drug ~ 0.053361 # Table S3 'Between-arm random effects in Emax,drug' SD = 0.231 (proportional); 0.231^2

    # Residual within-arm error, additive, for an arm of 100 patients.
    addSd <- 0.101
    label("Additive residual SD of the arm-mean DAS28 change for an arm of 100 patients; per-record SD is addSd * sqrt(100 / N_ARM) (DAS28 units)") # Table S3 'Residual within-arm random effect' SD = 0.101 (95% CI 0.0959 to 0.106; additive, footnote c)
  })

  model({
    # Per-row arm inputs: the seven drug indicators (all 0 = placebo /
    # background-therapy arm), N_ARM, and the trial-level covariates
    # SCORE_DAS28CRP, T_DIAG_RA, SEXF_PCT and YEAR_PUB. The model is
    # algebraic in time, so all arms of one trial can share one ID
    # (and therefore one set of between-trial etas) as separate rows.
    wN <- sqrt(100 / N_ARM)

    # ---- Placebo response (Eqs. 8-10) ----
    lowMale <- 0
    if ((100 - SEXF_PCT) < 18.5) lowMale <- 1
    emaxPbo <- emax_pbo + e_das28_emax_pbo * log(SCORE_DAS28CRP / 6.2) +
      e_tdiag_emax_pbo * log(T_DIAG_RA / 8.2) + eta_study_emax_pbo * wN
    et50Pbo <- exp(let50_pbo + e_das28_et50_pbo * log(SCORE_DAS28CRP / 6.2))
    hillPbo <- exp(lhill_pbo + e_lowmale_hill_pbo * lowMale + eta_study_hill_pbo * wN)

    # ---- Drug treatment effect (Eqs. 6, 7) ----
    onDrug <- ABATACEPT + ADALIMUMAB + CERTOLIZUMAB + ETANERCEPT +
      RITUXIMAB + TOCILIZUMAB + TOFACITINIB
    emaxDrugTv <- exp(lemax_abatacept) * ABATACEPT +
      exp(lemax_adalimumab) * ADALIMUMAB +
      exp(lemax_certolizumab) * CERTOLIZUMAB +
      exp(lemax_etanercept) * ETANERCEPT +
      exp(lemax_rituximab) * RITUXIMAB +
      exp(lemax_tocilizumab) * TOCILIZUMAB +
      exp(lemax_tofacitinib) * TOFACITINIB
    emaxDrug <- emaxDrugTv * exp((eta_study_emax_drug + eta_arm_emax_drug) * wN)
    # ET50 of 1 week on a placebo arm is a placeholder that keeps the
    # Hill fraction finite; it is multiplied by emaxDrug = 0.
    et50Drug <- exp(let50_abatacept) * ABATACEPT +
      exp(let50_adalimumab) * ADALIMUMAB +
      exp(let50_certolizumab) * CERTOLIZUMAB +
      exp(let50_etanercept) * ETANERCEPT +
      exp(let50_rituximab) * RITUXIMAB +
      exp(let50_tocilizumab) * TOCILIZUMAB +
      exp(let50_tofacitinib) * TOFACITINIB + (1 - onDrug)
    hillDrug <- exp(lhill_drug + eta_study_hill_drug * wN)

    # ---- Progression (Eq. 12); PROG is per year, time in weeks ----
    slopeProg <- exp(lslope + e_year_slope * log(YEAR_PUB / 2013) + eta_study_slope * wN)

    effPbo <- emaxPbo * time^hillPbo / (et50Pbo^hillPbo + time^hillPbo)
    effDrug <- emaxDrug * time^hillDrug / (et50Drug^hillDrug + time^hillDrug)

    # Arm-mean change from baseline in DAS28-CRP (negative = improvement).
    das28cfb <- -(effPbo + effDrug) + slopeProg * time / (365.25 / 7)

    addSdArm <- addSd * wN
    das28cfb ~ add(addSdArm)
  })
}
