Jeong_2022_torsemide_pbpk <- function() {
  description <- paste0(
    "PBPK-PD (minimal, perfusion-limited, coded in Berkeley Madonna). ",
    "Torsemide after oral or intravenous dosing in healthy adults, CYP2C9 ",
    "extensive / intermediate / poor metabolizers, and patients with liver ",
    "cirrhosis (Child-Pugh A / B / C), chronic kidney disease (mild / ",
    "moderate / severe) or congestive heart failure (NYHA II / III / IV) ",
    "(Jeong, Jang and Lee 2022, Pharmaceutics). A single blood pool feeds ",
    "lung, kidney, liver, GI-tract and rest-of-body tissues in parallel; ",
    "every blood-tissue flow is scaled by the unbound fraction times the ",
    "plasma-to-blood-cell partition ratio. Oral dose is absorbed first-order ",
    "into the GI-tract tissue, which drains to the liver and loses drug by a ",
    "first-order GI elimination; Michaelis-Menten (CYP2C9) hepatic ",
    "metabolism and first-order urinary excretion from the kidney. The ",
    "urinary sodium excretion rate is a direct baseline sigmoid-Emax ",
    "function of the urinary torsemide excretion rate, with separate ",
    "parameter sets for healthy (also CKD and CYP2C9), cirrhosis and CHF. ",
    "Disease and phenotype scenarios are selected by binary covariates. ",
    "Deterministic: the paper reports no inter-individual or residual ",
    "variability, so the model is for typical-value simulation."
  )
  reference <- paste0(
    "Jeong S-H, Jang J-H, Lee Y-B. Torsemide Pharmacometrics in Healthy ",
    "Adult Populations Including CYP2C9 Genetic Polymorphisms and Various ",
    "Patient Groups through Physiologically Based Pharmacokinetic-",
    "Pharmacodynamic Modeling. Pharmaceutics. 2022;14(12):2720. ",
    "doi:10.3390/pharmaceutics14122720"
  )
  vignette <- "Jeong_2022_torsemide_pbpk"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    CYP2C9_EM = list(
      description = "CYP2C9 extensive-metabolizer phenotype indicator (1 = EM, e.g. *1/*1)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CYP2C9 phenotype not specified: the Table S2 population-average metabolism, when CYP2C9_EM, CYP2C9_IM and CYP2C9_PM are all 0)",
      notes = paste(
        "Section 2.6.1: 'extensive-metabolizers (such as *1*1; EM)'.",
        "CYP2C9_EM, CYP2C9_IM and CYP2C9_PM are mutually exclusive; all",
        "three 0 selects the Table S2 base Vmax 1.11 ug/mL/h and Km 3.90",
        "ug/mL, which every disease scenario of the paper scales.",
        "CYP2C9_EM = 1 replaces them with Table 1 Vmax,EM 2.216 ug/mL/h",
        "and Km,EM 1.951 ug/mL."
      ),
      source_name = "EM"
    ),
    CYP2C9_IM = list(
      description = "CYP2C9 intermediate-metabolizer phenotype indicator (1 = IM, e.g. *1/*3)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CYP2C9 phenotype not specified when CYP2C9_EM, CYP2C9_IM and CYP2C9_PM are all 0)",
      notes = paste(
        "Section 2.6.1: 'intermediate-metabolizers (such as *1*3; IM)'.",
        "CYP2C9_IM = 1 selects Table 1 Vmax,IM 0.886 ug/mL/h and Km,IM",
        "5.073 ug/mL. Mutually exclusive with CYP2C9_EM and CYP2C9_PM."
      ),
      source_name = "IM"
    ),
    CYP2C9_PM = list(
      description = "CYP2C9 poor-metabolizer phenotype indicator (1 = PM, e.g. *3/*3)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CYP2C9 phenotype not specified when CYP2C9_EM, CYP2C9_IM and CYP2C9_PM are all 0)",
      notes = paste(
        "Section 2.6.1: 'poor-metabolizers (such as *3*3; PM)'.",
        "CYP2C9_PM = 1 selects Table 1 Vmax,PM 0.222 ug/mL/h and Km,PM",
        "8.242 ug/mL. Mutually exclusive with CYP2C9_EM and CYP2C9_IM."
      ),
      source_name = "PM"
    ),
    HEPIMP_MILD = list(
      description = "Liver cirrhosis, Child-Pugh class A indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no cirrhosis when HEPIMP_MILD, HEPIMP_MOD and HEPIMP_SEV are all 0)",
      notes = paste(
        "Child-Pugh classification (Section 2.6.2). Selects the CP-A column",
        "of Table 2 (physiology) and Figure S10 (partition coefficients),",
        "and the cirrhosis PD parameters of Table 5 (fitted to CP-B data",
        "and applied to every Child-Pugh class, Section 3.4)."
      ),
      source_name = "CP-A"
    ),
    HEPIMP_MOD = list(
      description = "Liver cirrhosis, Child-Pugh class B indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no cirrhosis when HEPIMP_MILD, HEPIMP_MOD and HEPIMP_SEV are all 0)",
      notes = "Child-Pugh classification (Section 2.6.2). Selects the CP-B column of Table 2 and Figure S10 and the cirrhosis PD parameters of Table 5.",
      source_name = "CP-B"
    ),
    HEPIMP_SEV = list(
      description = "Liver cirrhosis, Child-Pugh class C indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no cirrhosis when HEPIMP_MILD, HEPIMP_MOD and HEPIMP_SEV are all 0)",
      notes = "Child-Pugh classification (Section 2.6.2). Selects the CP-C column of Table 2 and Figure S10 and the cirrhosis PD parameters of Table 5.",
      source_name = "CP-C"
    ),
    RENALIMP_MILD = list(
      description = "Mild chronic kidney disease indicator (GFR 60-89 mL/min/1.73 m^2)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no CKD when RENALIMP_MILD, RENALIMP_MOD and RENALIMP_SEV are all 0)",
      notes = paste(
        "Table 3 footnote: 'Patients with GFR values of 60-89, 30-59, and",
        "<30 mL/min/1.73 m2 were considered mild, moderate, and severe'.",
        "Selects the Table 3 Mild column. The paper has no CKD PD model;",
        "the healthy Table 5 PD parameters are used, as in its Figure 4."
      ),
      source_name = "CKD mild"
    ),
    RENALIMP_MOD = list(
      description = "Moderate chronic kidney disease indicator (GFR 30-59 mL/min/1.73 m^2)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no CKD when RENALIMP_MILD, RENALIMP_MOD and RENALIMP_SEV are all 0)",
      notes = "Selects the Table 3 Moderate column; healthy PD parameters (Section 3.9, Figure 4).",
      source_name = "CKD moderate"
    ),
    RENALIMP_SEV = list(
      description = "Severe chronic kidney disease indicator (GFR < 30 mL/min/1.73 m^2)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no CKD when RENALIMP_MILD, RENALIMP_MOD and RENALIMP_SEV are all 0)",
      notes = "Selects the Table 3 Severe column; healthy PD parameters (Section 3.9, Figure 4).",
      source_name = "CKD severe"
    ),
    DIS_CHF_NYHA2 = list(
      description = "Mild congestive heart failure, NYHA class II indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no heart failure when DIS_CHF_NYHA2, DIS_CHF_NYHA3 and DIS_CHF_NYHA4 are all 0)",
      notes = paste(
        "Table 4 footnote: 'NYHA classes II, III, and IV were considered",
        "mild, moderate, and severe'. Selects the Table 4 Mild column and",
        "the CHF PD parameters of Table 5 (fitted to mild-to-moderate data",
        "and applied to every NYHA class, Section 3.6)."
      ),
      source_name = "CHF mild"
    ),
    DIS_CHF_NYHA3 = list(
      description = "Moderate congestive heart failure, NYHA class III indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no heart failure when DIS_CHF_NYHA2, DIS_CHF_NYHA3 and DIS_CHF_NYHA4 are all 0)",
      notes = "Selects the Table 4 Moderate column and the CHF PD parameters of Table 5.",
      source_name = "CHF moderate"
    ),
    DIS_CHF_NYHA4 = list(
      description = "Severe congestive heart failure, NYHA class IV indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no heart failure when DIS_CHF_NYHA2, DIS_CHF_NYHA3 and DIS_CHF_NYHA4 are all 0)",
      notes = "Selects the Table 4 Severe column and the CHF PD parameters of Table 5.",
      source_name = "CHF severe"
    )
  )

  # Every state is an AMOUNT of torsemide in mg; model() divides by the
  # Table S2 volume (converted to L) to obtain the concentration in mg/L,
  # numerically equal to ug/mL.
  compartmentData <- list(
    depot = list(analyte = "torsemide", units = "mg", specimen = "administration site", verified = TRUE),
    a_blood = list(analyte = "torsemide", units = "mg", specimen = "whole blood", verified = TRUE),
    a_lung = list(analyte = "torsemide", units = "mg", specimen = "tissue", verified = TRUE),
    a_gut = list(analyte = "torsemide", units = "mg", specimen = "tissue", verified = TRUE),
    a_liver = list(analyte = "torsemide", units = "mg", specimen = "tissue", verified = TRUE),
    a_kidney = list(analyte = "torsemide", units = "mg", specimen = "tissue", verified = TRUE),
    a_urine = list(analyte = "torsemide", units = "mg", specimen = "urine", verified = TRUE),
    a_remainder = list(analyte = "torsemide", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 17L,
    age_range = "19-83 years across the pooled healthy-adult studies (Table S1)",
    weight_range = NULL,
    sex_female_pct = NA_real_,
    race_ethnicity = "Mixed: Korean bioequivalence and pop-PK cohorts plus European and US clinical-pharmacology studies (Table S1).",
    disease_state = "Healthy adults (including CYP2C9 EM / IM / PM) for model building; liver cirrhosis, chronic kidney disease and congestive heart failure patients by parameter substitution.",
    dose_range = "Healthy: 5, 10 and 20 mg single oral or IV, 10 mg once daily oral; CKD: 100 mg oral or IV; cirrhosis: 10 mg oral or IV; CHF: 5-200 mg oral or IV (Table S1).",
    regions = "Republic of Korea, Europe and USA (literature studies of Table S1).",
    notes = paste(
      "No individual-level data set. Mean (or median) plasma / serum",
      "concentration, cumulative urinary excretion and urinary sodium",
      "excretion profiles were digitised from 11 healthy-adult studies",
      "(Jeong 2022, Vormfelde 2004, Barr 1990, Kramer 1994, Vargo 1995,",
      "Schwartz 1993, Kang 2013, Barbanoj 2009, Shah 2017, Cho 2005,",
      "Spahn 1990), one CKD study (Gehr 1994, n = 24), one cirrhosis study",
      "(Schwartz 1993, n = 12) and six CHF studies (Vargo 1995, Vargo 1994,",
      "Ballester 2015, Bleske 1998, Kramer 1996, Hariman 1994), together",
      "with individual bioequivalence data (Table S1; Section 2.1). Most",
      "data came from male adults. Organ volumes and flows describe a",
      "reference adult (Table S2)."
    )
  )

  ini({
    # ----------------------------------------------------------------
    # Organ volumes -- Table S2, printed in mL and stored here in L so that
    # amount (mg) / volume (L) is mg/L = ug/mL, the unit of Vmax and Km.
    # ----------------------------------------------------------------
    v_blood <- fixed(5.4)
    label("Volume of the blood pool Vbl (L)") # Table S2 'Vbl' = 5400 mL (Literature)
    v_lung <- fixed(1.2)
    label("Volume of lung Vlu (L)") # Table S2 'Vlu' = 1200 mL (Literature)
    v_kidney <- fixed(0.3)
    label("Volume of kidney Vki (L)") # Table S2 'Vki' = 300 mL (Literature)
    v_liver <- fixed(1.5)
    label("Volume of liver Vli (L)") # Table S2 'Vli' = 1500 mL (Literature)
    v_gut <- fixed(2.4)
    label("Volume of GI-tract Vgi (L)") # Table S2 'Vgi' = 2400 mL (Literature)
    lv_remainder <- log(62.43)
    label("Log volume of rest of body Vrb (log L)") # Table S2 'Vrb' = 6.243 x 10^4 mL (Estimated)

    # ----------------------------------------------------------------
    # Blood flows -- Table S2, printed in mL/h and stored here in L/h.
    # Qbl (3.798 x 10^5 mL/h) equals Qlu and does not enter the ODEs of
    # Appendix A; the hepatic-vein flow Qlive (1.53 x 10^5 mL/h) is formed
    # in model() as Qliar + Qgi, which reproduces the printed value exactly
    # and keeps the blood pool in mass balance.
    # ----------------------------------------------------------------
    q_lung <- fixed(379.8)
    label("Blood flow to lung Qlu = cardiac output Qbl (L/h)") # Table S2 'Qlu' = 'Qbl' = 3.798 x 10^5 mL/h (Literature)
    q_kidney <- fixed(74.4)
    label("Blood flow to kidney Qki (L/h)") # Table S2 'Qki' = 7.44 x 10^4 mL/h (Literature)
    q_gut <- fixed(66)
    label("Blood flow to GI-tract Qgi (L/h)") # Table S2 'Qgi' = 6.6 x 10^4 mL/h (Literature)
    q_hepatic_artery <- fixed(87)
    label("Blood flow in the hepatic artery Qliar (L/h)") # Table S2 'Qliar' = 8.7 x 10^4 mL/h (Literature)
    lq_remainder <- log(152.4)
    label("Log blood flow to rest of body Qrb (log L/h)") # Table S2 'Qrb' = 1.524 x 10^5 mL/h (Estimated)

    # ----------------------------------------------------------------
    # Blood-to-tissue partition coefficients -- Table S2, predicted with
    # the Simcyp Simulator v21 from the physicochemical data of Table S3.
    # ----------------------------------------------------------------
    lkp_lung <- fixed(log(0.23))
    label("Log lung partition coefficient Klu (unitless)") # Table S2 'Klu' = 0.23 (Simulated)
    lkp_kidney <- fixed(log(0.15))
    label("Log kidney partition coefficient Kki (unitless)") # Table S2 'Kki' = 0.15 (Simulated)
    lkp_liver <- fixed(log(0.11))
    label("Log liver partition coefficient Kli (unitless)") # Table S2 'Kli' = 0.11 (Simulated)
    lkp_gut <- fixed(log(0.18))
    label("Log GI-tract partition coefficient Kgi (unitless)") # Table S2 'Kgi' = 0.18 (Simulated)
    lkp_remainder <- fixed(log(1.03))
    label("Log rest-of-body partition coefficient Krb (unitless)") # Table S2 'Krb' = 1.03 (Simulated)

    # ----------------------------------------------------------------
    # Binding and blood-cell distribution -- Table S2.
    # ----------------------------------------------------------------
    fu <- fixed(0.05)
    label("Unbound fraction in plasma or serum Fr (unitless)") # Table S2 'Fr' = 0.05 (Literature; human serum albumin equilibrium dialysis)
    krbc <- fixed(0.28)
    label("Plasma-to-blood-cell partition ratio Pr (unitless)") # Table S2 'Pr' = 0.28 (Literature; measured in rat blood)

    # ----------------------------------------------------------------
    # Hepatic Michaelis-Menten metabolism -- Table S2 base values (human
    # liver microsomes) and Table 1 CYP2C9-phenotype values, the latter as
    # log shifts so that exp(lvmax + e_cyp2c9em_lvmax) = 2.216 exactly.
    # ----------------------------------------------------------------
    lvmax <- fixed(log(1.11))
    label("Log maximum hepatic metabolic rate Vmax (log ug/mL/h)") # Table S2 'Vmax' = 1.11 ug/mL/h (Literature)
    lkm <- fixed(log(3.90))
    label("Log Michaelis-Menten constant Km (log ug/mL)") # Table S2 'Km' = 3.90 ug/mL (Literature)
    e_cyp2c9em_lvmax <- log(2.216 / 1.11)
    label("CYP2C9 EM log shift on Vmax (unitless)") # Table 1 'Vmax,EM' = 2.216 ug/mL/h (Estimated) vs Table S2 1.11
    e_cyp2c9em_lkm <- log(1.951 / 3.90)
    label("CYP2C9 EM log shift on Km (unitless)") # Table 1 'Km,EM' = 1.951 ug/mL (Estimated) vs Table S2 3.90
    e_cyp2c9im_lvmax <- log(0.886 / 1.11)
    label("CYP2C9 IM log shift on Vmax (unitless)") # Table 1 'Vmax,IM' = 0.886 ug/mL/h (Estimated) vs Table S2 1.11
    e_cyp2c9im_lkm <- log(5.073 / 3.90)
    label("CYP2C9 IM log shift on Km (unitless)") # Table 1 'Km,IM' = 5.073 ug/mL (Estimated) vs Table S2 3.90
    e_cyp2c9pm_lvmax <- log(0.222 / 1.11)
    label("CYP2C9 PM log shift on Vmax (unitless)") # Table 1 'Vmax,PM' = 0.222 ug/mL/h (Estimated) vs Table S2 1.11
    e_cyp2c9pm_lkm <- log(8.242 / 3.90)
    label("CYP2C9 PM log shift on Km (unitless)") # Table 1 'Km,PM' = 8.242 ug/mL (Estimated) vs Table S2 3.90

    # ----------------------------------------------------------------
    # First-order rate constants -- Table S2, fitted to the healthy-adult
    # profiles.
    # ----------------------------------------------------------------
    lka <- log(5.50)
    label("Log oral absorption rate constant Ka (log 1/h)") # Table S2 'Ka' = 5.50 1/h (Estimated)
    lkurine <- log(18.50)
    label("Log urinary excretion rate constant from kidney Ku (log 1/h)") # Table S2 'Ku' = 18.50 1/h (Estimated)
    lkgutex <- log(1.00)
    label("Log elimination rate constant from the GI-tract Ke (log 1/h)") # Table S2 'Ke' = 1.00 1/h (Estimated)

    # ----------------------------------------------------------------
    # PD: baseline sigmoid Emax of the urinary sodium excretion rate on
    # the urinary torsemide excretion rate -- Table 5. The printed slope
    # P is negative (-2.13); hill = -P (see model()). Cirrhosis and CHF
    # parameter sets enter as shifts from the healthy set.
    # ----------------------------------------------------------------
    le0 <- log(4.39)
    label("Log basal urinary sodium excretion rate E0, healthy (log mEq/h)") # Table 5 Healthy 'E0' = 4.39 mEq/h
    lemax <- log(105.1)
    label("Log maximal urinary sodium excretion rate Emax, healthy (log mEq/h)") # Table 5 Healthy 'Emax' = 105.1 mEq/h
    lec50 <- log(345.5)
    label("Log urinary torsemide excretion rate at half Emax EC50, healthy (log ug/h)") # Table 5 Healthy 'EC50' = 345.5 ug/h
    hill <- 2.13
    label("Sigmoid slope (-P), healthy (unitless)") # Table 5 Healthy 'P' = -2.13
    e_hepimp_le0 <- log(4.33 / 4.39)
    label("Cirrhosis log shift on E0 (unitless)") # Table 5 Cirrhosis 'E0' = 4.33 mEq/h
    e_hepimp_lemax <- log(65.24 / 105.1)
    label("Cirrhosis log shift on Emax (unitless)") # Table 5 Cirrhosis 'Emax' = 65.24 mEq/h
    e_hepimp_lec50 <- log(388.8 / 345.5)
    label("Cirrhosis log shift on EC50 (unitless)") # Table 5 Cirrhosis 'EC50' = 388.8 ug/h
    e_hepimp_hill <- 2.82 - 2.13
    label("Cirrhosis additive shift on the sigmoid slope (unitless)") # Table 5 Cirrhosis 'P' = -2.82
    e_chf_le0 <- log(3.25 / 4.39)
    label("CHF log shift on E0 (unitless)") # Table 5 CHF 'E0' = 3.25 mEq/h
    e_chf_lemax <- log(83.45 / 105.1)
    label("CHF log shift on Emax (unitless)") # Table 5 CHF 'Emax' = 83.45 mEq/h
    e_chf_lec50 <- log(3985 / 345.5)
    label("CHF log shift on EC50 (unitless)") # Table 5 CHF 'EC50' = 3985 ug/h
    e_chf_hill <- 1.16 - 2.13
    label("CHF additive shift on the sigmoid slope (unitless)") # Table 5 CHF 'P' = -1.16

    # ----------------------------------------------------------------
    # Residual error. The paper fitted mean profiles in Berkeley Madonna
    # and reports no residual-error model; encoded as zero rather than
    # invented.
    # ----------------------------------------------------------------
    propSd <- fixed(0)
    label("Proportional residual error (fraction; ZERO - not reported in source)")
  })

  model({
    # ---------------------------------------------------------------
    # Scenario selectors. Within each disease the indicators are mutually
    # exclusive; the paper simulates one disease group at a time, so at
    # most one disease group should be active.
    # ---------------------------------------------------------------
    hep_any <- HEPIMP_MILD + HEPIMP_MOD + HEPIMP_SEV
    ckd_any <- RENALIMP_MILD + RENALIMP_MOD + RENALIMP_SEV
    chf_any <- DIS_CHF_NYHA2 + DIS_CHF_NYHA3 + DIS_CHF_NYHA4

    # Table 2 (cirrhosis) fractions of the control value, CP-A / B / C.
    f_liver_vol_hep <- 1 - hep_any + 0.81 * HEPIMP_MILD + 0.65 * HEPIMP_MOD + 0.53 * HEPIMP_SEV
    f_co_hep <- 1 - hep_any + 1.16 * HEPIMP_MILD + 1.32 * HEPIMP_MOD + 1.41 * HEPIMP_SEV
    f_qki_hep <- 1 - hep_any + 0.70 * HEPIMP_MILD + 0.58 * HEPIMP_MOD + 0.55 * HEPIMP_SEV
    f_portal_hep <- 1 - hep_any + 0.91 * HEPIMP_MILD + 0.63 * HEPIMP_MOD + 0.55 * HEPIMP_SEV
    f_qliar_hep <- 1 - hep_any + 1.41 * HEPIMP_MILD + 1.62 * HEPIMP_MOD + 1.91 * HEPIMP_SEV
    f_villous_hep <- 1 - hep_any + 1.29 * HEPIMP_MILD + 1.52 * HEPIMP_MOD + 1.99 * HEPIMP_SEV
    # Table 2 hematocrit (%), albumin + alpha1-acid glycoprotein (g/L) and
    # hepatic CYP2C9 abundance (pmol/mg); the control column is 43 / 50.34 / 73.
    hct <- 43 * (1 - hep_any) + 38.46 * HEPIMP_MILD + 34.51 * HEPIMP_MOD + 33.51 * HEPIMP_SEV
    protein <- 50.34 * (1 - hep_any) + 46.33 * HEPIMP_MILD + 38.18 * HEPIMP_MOD + 29.69 * HEPIMP_SEV
    cyp2c9_abund <- 73 * (1 - hep_any) + 50.4 * HEPIMP_MILD + 38.0 * HEPIMP_MOD + 24.1 * HEPIMP_SEV

    # Table 3 (CKD) fractions, mild / moderate / severe.
    f_kidney_vol_ckd <- 1 - ckd_any + 0.99 * RENALIMP_MILD + 0.93 * RENALIMP_MOD + 0.87 * RENALIMP_SEV
    f_qki_ckd <- 1 - ckd_any + 0.99 * RENALIMP_MILD + 0.86 * RENALIMP_MOD + 0.72 * RENALIMP_SEV
    f_clint_ckd <- 1 - ckd_any + 1.00 * RENALIMP_MILD + 0.69 * RENALIMP_MOD + 0.52 * RENALIMP_SEV
    f_renal_ckd <- 1 - ckd_any + 0.99 * RENALIMP_MILD + 0.68 * RENALIMP_MOD + 0.37 * RENALIMP_SEV

    # Table 4 (CHF) fractions, mild / moderate / severe.
    f_co_chf <- 1 - chf_any + 0.69 * DIS_CHF_NYHA2 + 0.50 * DIS_CHF_NYHA3 + 0.42 * DIS_CHF_NYHA4
    f_qki_chf <- 1 - chf_any + 0.78 * DIS_CHF_NYHA2 + 0.55 * DIS_CHF_NYHA3 + 0.63 * DIS_CHF_NYHA4
    f_portal_chf <- 1 - chf_any + 0.76 * DIS_CHF_NYHA2 + 0.54 * DIS_CHF_NYHA3 + 0.46 * DIS_CHF_NYHA4
    f_qliar_chf <- 1 - chf_any + 0.76 * DIS_CHF_NYHA2 + 0.54 * DIS_CHF_NYHA3 + 0.46 * DIS_CHF_NYHA4
    f_other_chf <- 1 - chf_any + 0.57 * DIS_CHF_NYHA2 + 0.44 * DIS_CHF_NYHA3 + 0.28 * DIS_CHF_NYHA4
    f_clint_chf <- 1 - chf_any + 0.76 * DIS_CHF_NYHA2 + 0.54 * DIS_CHF_NYHA3 + 0.46 * DIS_CHF_NYHA4

    # ---------------------------------------------------------------
    # Physiology. The disease rules are those of Sections 2.6.2-2.6.4;
    # each reproduces the bars of Figures S10, S14 and S16.
    # ---------------------------------------------------------------
    vli <- v_liver * f_liver_vol_hep
    vki <- v_kidney * f_kidney_vol_ckd
    vrb <- exp(lv_remainder)

    qlu <- q_lung * f_co_hep * f_co_chf
    qki <- q_kidney * f_qki_hep * f_qki_ckd * f_qki_chf
    qgi <- q_gut * f_portal_hep * f_portal_chf
    qliar <- q_hepatic_artery * f_qliar_hep * f_qliar_chf
    # Hepatic-vein outflow Qlive = Qliar + Qgi (Table S2: 1.53e5 = 8.7e4 +
    # 6.6e4 mL/h; Section 2.6.4 'Qlive reflected results of Qliar and Qgi').
    qlive <- qliar + qgi
    # Rest-of-body flow. Cirrhosis: cardiac output minus kidney and hepatic
    # flows (Figure S10 Qrb). CKD: the kidney flow lost is added to the
    # rest of the body (Section 2.6.3). CHF: Table 4 'other tissues'
    # fraction. In the healthy state the estimated Qrb also equals
    # Qbl - Qki - Qlive (379.8 - 74.4 - 153 = 152.4 L/h).
    qrb <- hep_any * (qlu - qki - qlive) +
      (1 - hep_any) * (exp(lq_remainder) * f_other_chf + ckd_any * (q_kidney - qki))

    # Unbound fraction in cirrhosis, Section 2.6.2:
    #   Fr,cirr = 1 / (1 + (1 - Fr) * PTcirr / (PThealthy * Fr));
    # identity at the control protein concentration.
    fr <- 1 / (1 + (1 - fu) * (protein / 50.34) / fu)
    # Plasma-to-blood-cell ratio scaled by the healthy / cirrhosis
    # hematocrit ratio, Section 2.6.2.
    pr <- krbc * 43 / hct
    fp <- fr * pr

    # Partition coefficients. Cirrhosis values are the Simcyp predictions
    # with the changed Fr and Pr, shown only graphically in Figure S10 and
    # digitised from it (CP-A / B / C).
    kp_lung <- (1 - hep_any) * exp(lkp_lung) + 0.233 * HEPIMP_MILD + 0.237 * HEPIMP_MOD + 0.243 * HEPIMP_SEV
    kp_kidney <- (1 - hep_any) * exp(lkp_kidney) + 0.157 * HEPIMP_MILD + 0.160 * HEPIMP_MOD + 0.170 * HEPIMP_SEV
    kp_liver <- (1 - hep_any) * exp(lkp_liver) + 0.113 * HEPIMP_MILD + 0.117 * HEPIMP_MOD + 0.127 * HEPIMP_SEV
    kp_gut <- (1 - hep_any) * exp(lkp_gut) + 0.180 * HEPIMP_MILD + 0.187 * HEPIMP_MOD + 0.193 * HEPIMP_SEV
    kp_remainder <- (1 - hep_any) * exp(lkp_remainder) + 1.041 * HEPIMP_MILD + 1.078 * HEPIMP_MOD + 1.139 * HEPIMP_SEV

    # Hepatic metabolism. CYP2C9 phenotype sets Vmax and Km (Table 1);
    # cirrhosis scales Vmax by the CYP2C9 abundance ratio (Section 2.6.2),
    # CKD and CHF by the liver intrinsic-clearance fraction (Sections
    # 2.6.3 and 2.6.4).
    vmax <- exp(lvmax + e_cyp2c9em_lvmax * CYP2C9_EM + e_cyp2c9im_lvmax * CYP2C9_IM + e_cyp2c9pm_lvmax * CYP2C9_PM) *
      (cyp2c9_abund / 73) * f_clint_ckd * f_clint_chf
    km <- exp(lkm + e_cyp2c9em_lkm * CYP2C9_EM + e_cyp2c9im_lkm * CYP2C9_IM + e_cyp2c9pm_lkm * CYP2C9_PM)

    ka <- exp(lka)
    # Urinary excretion scaled by the CKD renal-elimination fraction.
    kurine <- exp(lkurine) * f_renal_ckd
    # GI elimination divided by the villous-flow fraction in cirrhosis
    # (Figure S10 Ke = 0.78 / 0.66 / 0.50 1/h = 1.00 / 1.29, 1.52, 1.99).
    kgutex <- exp(lkgutex) / f_villous_hep

    # ---------------------------------------------------------------
    # Mass balances, Appendix A Equations (A1)-(A10), written for amounts
    # (mg). Every blood-tissue flow is multiplied by Fr * Pr as printed.
    # ---------------------------------------------------------------
    c_blood <- a_blood / v_blood
    c_lung <- a_lung / v_lung
    c_gut <- a_gut / v_gut
    c_liver <- a_liver / vli
    c_kidney <- a_kidney / vki
    c_remainder <- a_remainder / vrb

    # (A3)-(A4): oral dose to the oral compartment, absorbed first-order.
    d/dt(depot) <- -ka * depot
    # (A1); an IV bolus (A2) is a dose into a_blood.
    d/dt(a_blood) <- fp * (qlu * c_lung / kp_lung + qki * c_kidney / kp_kidney +
      qlive * c_liver / kp_liver + qrb * c_remainder / kp_remainder) -
      fp * (qlu + qki + qliar + qgi + qrb) * c_blood
    # (A9) lung, in parallel with the other tissues.
    d/dt(a_lung) <- fp * qlu * (c_blood - c_lung / kp_lung)
    # (A5) GI-tract tissue.
    d/dt(a_gut) <- fp * qgi * (c_blood - c_gut / kp_gut) + ka * depot - kgutex * c_gut * v_gut
    # (A6) liver.
    d/dt(a_liver) <- fp * (qliar * c_blood - qlive * c_liver / kp_liver + qgi * c_gut / kp_gut) -
      vmax * c_liver / (km + c_liver) * vli
    # (A7) kidney.
    d/dt(a_kidney) <- fp * qki * (c_blood - c_kidney / kp_kidney) - kurine * c_kidney * vki
    # (A8) cumulative urinary excretion.
    d/dt(a_urine) <- kurine * c_kidney * vki
    # (A10) rest of body.
    d/dt(a_remainder) <- fp * qrb * (c_blood - c_remainder / kp_remainder)

    # ---------------------------------------------------------------
    # PD, Table 5. x is the urinary torsemide excretion rate in ug/h
    # (mg/h * 1000). The table prints E = E0 + Emax / (1 + x/EC50)^P, which
    # gives E = E0 + Emax at x = 0 and grows without bound for P < 0; the
    # fitted curves of Figure S1 and the return of E to E0 in Figures S6
    # and S12 require E = E0 + Emax / (1 + (x/EC50)^P), which for P < 0 is
    # the sigmoid E0 + Emax * x^hill / (EC50^hill + x^hill), hill = -P.
    # ---------------------------------------------------------------
    e0 <- exp(le0 + e_hepimp_le0 * hep_any + e_chf_le0 * chf_any)
    emax <- exp(lemax + e_hepimp_lemax * hep_any + e_chf_lemax * chf_any)
    ec50 <- exp(lec50 + e_hepimp_lec50 * hep_any + e_chf_lec50 * chf_any)
    hill_pd <- hill + e_hepimp_hill * hep_any + e_chf_hill * chf_any
    urine_excr_rate <- kurine * c_kidney * vki * 1000
    sodium_excr_rate <- e0 + emax * urine_excr_rate^hill_pd / (ec50^hill_pd + urine_excr_rate^hill_pd)

    # Observation: Appendix A has no separate plasma compartment; the
    # blood-pool concentration Cbl is the simulated 'plasma (or serum)'
    # concentration of Figures S3-S4 (mg/L = ug/mL).
    Cc <- c_blood
    Cc ~ prop(propSd)
  })
}
