Morales_2023_type1Diabetes <- function() {
  description <- paste(
    "Type 1 diabetes disease-progression joint model (Morales 2023) for",
    "individuals with two or more islet autoantibodies, fitted to the TrialNet",
    "Pathway to Prevention and TEDDY natural-history cohorts. The 2-h oral",
    "glucose tolerance test plasma glucose (GLU120, mg/dL) rises from a",
    "baseline DP0 as a sigmoid Emax function of time since the derived",
    "baseline (maximum rise DPmax, time to half-maximum DP50, Hill gamma), and",
    "the time to clinical type 1 diabetes diagnosis follows a Weibull hazard",
    "whose scale is shortened by the individual predicted GLU120 relative to",
    "112 mg/dL. Baseline covariates: GADA positivity and HbA1c on DPmax, IA-2A",
    "positivity on DP50, fasting glucose and the 2-h / fasting glucose ratio on",
    "DP0. No drug input.",
    sep = " "
  )
  reference <- paste(
    "Morales JF, Muse R, Podichetty JT, Burton J, David S, Lang P, Schmidt S,",
    "Romero K, O'Doherty I, Martin F, Campbell-Thompson M, Haller MJ,",
    "Atkinson MA, Kim S.",
    "Disease progression joint model predicts time to type 1 diabetes onset:",
    "Optimizing future type 1 diabetes prevention studies.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(7):1016-1028.",
    "doi:10.1002/psp4.12973.",
    "Parameter values from Table 2; the inter-individual correlations and the",
    "covariate centring values from the Supplementary Material S2 Monolix",
    "project file.",
    sep = " "
  )
  vignette <- "Morales_2023_type1Diabetes"
  units <- list(
    time = "year",
    dosing = "n/a (disease-progression model; no drug input)",
    concentration = "mg/dL (2-h OGTT plasma glucose glucose_2h); sur is a survival probability"
  )

  compartmentData <- list(
    cumhaz = list(
      analyte = "type 1 diabetes diagnosis hazard",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    GADA_POS = list(
      description = "Glutamic acid decarboxylase 65 autoantibody (GADA) present at the derived baseline (1 = positive, 0 = negative)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (GADA negative)",
      notes = "Time-fixed at the derived baseline. Categorical power effect on DPmax: DPmax = DPmax_pop * exp(e_gada_pos_emax * GADA_POS) (Supplementary Equation S2). 89.6% of the analysis cohort was GADA positive (Table 1).",
      source_name = "BSCORE_GADA"
    ),
    IA2A_POS = list(
      description = "Islet antigen-2 autoantibody (IA-2A) present at the derived baseline (1 = positive, 0 = negative)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (IA-2A negative)",
      notes = "Time-fixed at the derived baseline. Categorical power effect on DP50: DP50 = DP50_pop * exp(e_ia2a_pos_t50 * IA2A_POS) (Supplementary Equation S2). 61.9% of the analysis cohort was IA-2A positive (Table 1).",
      source_name = "BSCORE_IA2A"
    ),
    HBA1C = list(
      description = "Glycated hemoglobin at the derived baseline",
      units = "% (NGSP)",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed at the derived baseline. Power effect on DPmax centred on the weighted mean 5.10576% printed in the Supplementary Material S2 [COVARIATE] block (logtBSCORE_HbA1c = log(BSCORE_HbA1c / 5.10576)). Cohort mean (SD) 5.1 (0.3)% (Table 1); subjects with HbA1c >= 6.5% at the derived baseline were excluded.",
      source_name = "BSCORE_HbA1c"
    ),
    FPG = list(
      description = "Fasting plasma glucose (0-min OGTT glucose, GLU0) at the derived baseline",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed at the derived baseline. Power effect on DP0 centred on the weighted mean 88.0112 mg/dL printed in the Supplementary Material S2 [COVARIATE] block (logtBSCORE_GLU0 = log(BSCORE_GLU0 / 88.0112)). Cohort mean (SD) 88.8 (9.2) mg/dL (Table 1); subjects with fasting glucose >= 126 mg/dL were excluded.",
      source_name = "BSCORE_GLU0"
    ),
    OGTT_GLU_RATIO_BL = list(
      description = "Ratio of the 2-h to the fasting (0-min) OGTT plasma glucose at the derived baseline (GLU120 / GLU0)",
      units = "(unitless)",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed at the derived baseline. Power effect on DP0 centred on the weighted mean 1.25818 printed in the Supplementary Material S2 [COVARIATE] block (logtBSCORE_GLURATIO = log(BSCORE_GLURATIO / 1.25818)). Cohort mean (SD) 1.3 (0.4) (Table 1).",
      source_name = "BSCORE_GLURATIO"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at the derived baseline",
      units = "year",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested in the stepwise covariate search (Table 1; centring value 10.3268 years in Supplementary S2) but not retained in the final model.",
      source_name = "BAGE"
    ),
    BMI = list(
      description = "Body mass index at the derived baseline",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested in the stepwise covariate search (Table 1; centring value 19.2068 kg/m^2 in Supplementary S2) but not retained; the final Monolix project is named 'without_BMI'.",
      source_name = "BSCORE_BMI"
    ),
    SEXF = list(
      description = "Female sex",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Tested in the stepwise covariate search but not retained. The source column SEX was coded 0 = female, 1 = male (Table 1), so SEXF = 1 - SEX.",
      source_name = "SEX"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1265L,
    n_studies = 2L,
    age_range = "1-56 years at the derived baseline",
    age_median = "mean (SD) 12.7 (9.7) years at the derived baseline",
    weight_range = "not reported (BMI mean (SD) 20.1 (8.4) kg/m^2)",
    sex_female_pct = 43.9,
    race_ethnicity = "not reported",
    disease_state = paste(
      "Individuals at risk of type 1 diabetes with two or more islet",
      "autoantibodies (GADA, IA-2A, IAA, ZnT8) at the derived baseline and",
      "without diabetes at that visit (2-h OGTT glucose < 200 mg/dL, fasting",
      "glucose < 126 mg/dL, HbA1c < 6.5%)."
    ),
    dose_range = "n/a (natural-history cohorts; no treatment)",
    regions = "TrialNet Pathway to Prevention (United States, Canada, Finland, United Kingdom, Italy, Germany, Sweden, Australia, New Zealand) and TEDDY (United States, Finland, Sweden, Germany)",
    biomarkers = "Longitudinal 2-h OGTT plasma glucose (8341 observations); time to clinical type 1 diabetes diagnosis",
    notes = paste(
      "Table 1 summarises the curated dataset: 259 TEDDY and 1006 TrialNet",
      "participants; 56.1% male; 82.7% with a first-degree relative with type 1",
      "diabetes. Time zero is the derived baseline, the earliest visit with",
      ">= 2 autoantibodies and complete covariates; follow-up was truncated at",
      "6 years. Model development used a random 80% training split (n = 1016);",
      "the remaining 249 subjects formed the hold-out test set."
    )
  )

  ini({
    # Sigmoid Emax-in-time model for GLU120 (Eq. 1). Table 2 'Final parameter
    # estimates of the joint model'; reference subject GADA-, IA-2A- and at the
    # Supplementary S2 covariate centring values.
    le0 <- log(107.31); label("Typical baseline 2-h OGTT glucose DP0 (mg/dL)") # Table 2 'Typical value of DP0' = 107.31 mg/dL (RSE 0.414%)
    lemax <- log(208.72); label("Typical maximum rise in 2-h OGTT glucose DPmax (mg/dL)") # Table 2 'Typical value of DPmax' = 208.72 mg/dL (RSE 23.9%)
    lt50 <- log(6.09); label("Typical time from the derived baseline to half of DPmax, DP50 (year)") # Table 2 'Typical value of DP50' = 6.09 years (RSE 12.2%)
    lhill <- log(1.75); label("Typical Hill coefficient gamma of the GLU120 time course (unitless)") # Table 2 'Typical value of gamma' = 1.75 (RSE 21.5%)

    # Covariate effects (Supplementary Equations S1 / S2: exponent on the
    # centred log covariate, or log-scale shift for a 0/1 indicator).
    e_gada_pos_emax <- -0.73; label("Log-scale shift in DPmax for GADA positivity at the derived baseline (unitless)") # Table 2 'beta_DPmax_BSCORE_GADA_1' = -0.73 (RSE 29.5%)
    e_hba1c_emax <- 8.92; label("Power exponent of baseline HbA1c / 5.10576% on DPmax (unitless)") # Table 2 'beta_DPmax_logtBSCORE_HbA1c' = 8.92 (RSE 13.2%)
    e_ia2a_pos_t50 <- -0.71; label("Log-scale shift in DP50 for IA-2A positivity at the derived baseline (unitless)") # Table 2 'beta_DP50_BSCORE_IA2A_1' = -0.71 (RSE 16.5%)
    e_fpg_e0 <- 0.88; label("Power exponent of baseline fasting glucose / 88.0112 mg/dL on DP0 (unitless)") # Table 2 'beta_DP0_logtBSCORE_GLU0' = 0.88 (RSE 4.39%)
    e_ogtt_glu_ratio_bl_e0 <- 0.79; label("Power exponent of baseline GLU120/GLU0 ratio / 1.25818 on DP0 (unitless)") # Table 2 'beta_DP0_logtBSCORE_GLURATIO' = 0.79 (RSE 1.71%)

    # Weibull time-to-diagnosis model (Eqs. 2 and 4). The paper's scale Te is
    # a time; llam_haz carries its reciprocal so lam_haz = 1 / Te.
    llam_haz <- log(1 / 91.81); label("Reciprocal of the Weibull scale Te at GLU120 = 112 mg/dL (1/year)") # Table 2 'Te' = 91.81 years (RSE 4.33%)
    lalfa_haz <- log(1.20); label("Weibull shape p (unitless)") # Table 2 'p' = 1.20 (RSE 0.683%)
    e_glucose_2h_haz <- 6.13; label("Power exponent of predicted GLU120 / 112 mg/dL on the Weibull rate 1/Te (unitless)") # Table 2 'beta_Te' = 6.13 (RSE 2.29%)

    # Between-subject variability: Table 2 SDs, squared to variances. The
    # three correlations come from the Supplementary S2 [INDIVIDUAL] block
    # 'correlation = {level=id, r(S0, Emax), r(gamma, Emax), r(gamma, S0)}'
    # with the <PARAMETER> values corr_S0_Emax = 0.8487756,
    # corr_gamma_Emax = 0.7281359, corr_gamma_S0 = 0.3707610; Table 2 does
    # not report them. Order: DPmax, DP0, gamma.
    etalemax + etale0 + etalhill ~ c(
      1.69^2,
      0.8487756 * 1.69 * 0.043, 0.043^2,
      0.7281359 * 1.69 * 1.76, 0.3707610 * 0.043 * 1.76, 1.76^2
    ) # Table 2 'BSV DPmax' SD 1.69, 'BSV DP0' SD 0.043, 'BSV gamma' SD 1.76; correlations from Supplementary S2
    etalt50 ~ 0.72^2 # Table 2 'BSV DP50' SD 0.72 (CV 82%, shrinkage 77.4%); uncorrelated in Supplementary S2

    propSd_glucose_2h <- 0.16; label("Proportional residual error on 2-h OGTT glucose (fraction)") # Table 2 'RV proportional error' = 0.16 (RSE 1.23%)
  })

  model({
    # Individual GLU120 parameters (Supplementary Equations S1 / S2 and the
    # S2 [INDIVIDUAL] DEFINITION block; all log-normal).
    e0 <- exp(le0 + etale0 +
      e_fpg_e0 * log(FPG / 88.0112) +
      e_ogtt_glu_ratio_bl_e0 * log(OGTT_GLU_RATIO_BL / 1.25818))
    emax <- exp(lemax + etalemax +
      e_hba1c_emax * log(HBA1C / 5.10576) +
      e_gada_pos_emax * GADA_POS)
    t50 <- exp(lt50 + etalt50 + e_ia2a_pos_t50 * IA2A_POS)
    hill <- exp(lhill + etalhill)

    # Eq. 1: GLU120(t) = DP0 + DPmax * t^gamma / (DP50^gamma + t^gamma), with
    # t in years since the derived baseline (Supplementary S3 E = A + S).
    # The sigmoid fraction is evaluated in the two branches of Supplementary
    # S3 so that a large individual gamma (SD 1.76 on the log scale) never
    # forms Inf / Inf.
    if (t <= t50) {
      rsig <- (t / t50)^hill
      fsig <- rsig / (1 + rsig)
    } else {
      rsig <- (t50 / t)^hill
      fsig <- 1 / (1 + rsig)
    }
    glucose_2h <- e0 + emax * fsig

    # Eq. 4: Te' = exp(ln Te - beta_Te * ln(GLU120(t) / 112)), so the Weibull
    # rate 1/Te' = lam_haz * (GLU120 / 112)^beta_Te. The hazard uses the
    # individual prediction without residual error (Supplementary S3 'E').
    lam <- exp(llam_haz) * (glucose_2h / 112)^e_glucose_2h_haz
    alfa <- exp(lalfa_haz)
    hazard <- alfa * lam * (lam * t)^(alfa - 1)

    d/dt(cumhaz) <- hazard
    cumhaz(0) <- 0
    sur <- exp(-cumhaz)

    glucose_2h ~ prop(propSd_glucose_2h)
  })
}
