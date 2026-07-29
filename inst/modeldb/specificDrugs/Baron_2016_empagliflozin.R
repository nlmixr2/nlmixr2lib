Baron_2016_empagliflozin <- function() {
  description <- paste(
    "Two-compartment population PK with lagged first-order absorption for",
    "empagliflozin in patients with type 2 diabetes (T2DM), coupled with",
    "two indirect-response PK/PD models for fasting plasma glucose (FPG)",
    "and glycated hemoglobin (HbA1c). The drug effect on FPG elimination",
    "is driven by steady-state AUC (AUCss = DOSE_EMPA_MGD * 1e6 / MW / CL)",
    "via an Emax function (Gmax, AUC50); FPG in turn drives HbA1c",
    "production with a boundary-condition baseline (HbA1climit). Pooled",
    "popPK/PD analysis of 4065 T2DM patients (PK n = 2761 active) from",
    "two phase I, four phase II, and four phase III studies (Baron 2016",
    "Diabetes Therapy)."
  )
  reference <- paste(
    "Baron KT, Macha S, Broedl UC, Nock V, Retlich S, Riggs M.",
    "Population Pharmacokinetics and Exposure-Response (Efficacy and",
    "Safety/Tolerability) of Empagliflozin in Patients with Type 2 Diabetes.",
    "Diabetes Ther. 2016;7(3):455-471. doi:10.1007/s13300-016-0174-y"
  )
  vignette <- "Baron_2016_empagliflozin"
  units <- list(
    time          = "hour",
    dosing        = "mg empagliflozin (oral, once daily)",
    concentration = "Cc in nmol/L (= nM; converted from mg/L via MW 450.91 g/mol); FPG in mmol/L; HbA1c in % (NGSP)"
  )
  paper_specific_compartments <- c("hba1c")

  covariateData <- list(
    AGE = list(
      description        = "Age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effects centred at 50 years on CL, V2/F, V3/F, ka (Table S1), on the model-predicted baseline FPG and Gmax (Table S3), and on kHbA1cout (Table S4).",
      source_name        = "AGE"
    ),
    BMI = list(
      description        = "Body mass index at baseline",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effects centred at 25 kg/m^2 on CL, V3/F, BFPG, Gmax, kHbA1cout. The single V2/F effect is centred at 20 kg/m^2 (Baron 2016 Table S1 row V2/F BMI).",
      source_name        = "BMI"
    ),
    CRCL = list(
      description        = "Modification-of-diet-in-renal-disease (MDRD) estimated glomerular filtration rate, BSA-normalised",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effects centred at 100 mL/min/1.73 m^2 on CL (Table S1), BFPG, Gmax (Table S3), and kHbA1cout (Table S4). Renamed from the source column eGFR to the canonical CRCL per covariate-columns.md.",
      source_name        = "eGFR"
    ),
    TPRO = list(
      description        = "Total serum protein",
      units              = "g/L (SI); paper reports g/dL with reference value 70 g/dL.",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effects centred at 70 g/dL on CL, V2/F, V3/F (Table S1). The model converts the canonical g/L column to g/dL inline (TPRO / 10).",
      source_name        = "TPRO"
    ),
    ALT = list(
      description        = "Alanine aminotransferase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on CL centred at 20 U/L (Table S1).",
      source_name        = "ALT"
    ),
    AST = list(
      description        = "Aspartate aminotransferase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on CL centred at 20 U/L (Table S1).",
      source_name        = "AST"
    ),
    ALP = list(
      description        = "Alkaline phosphatase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on CL centred at 70 U/L (Table S1). Renamed from the source column AP to the canonical ALP per covariate-columns.md.",
      source_name        = "AP"
    ),
    LDH = list(
      description        = "Lactate dehydrogenase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on CL centred at 160 U/L (Table S1).",
      source_name        = "LDH"
    ),
    SEXF = list(
      description        = "Sex (female indicator)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Multiplicative effects on CL, V2/F, V3/F, ka (Table S1), BFPG, Gmax (Table S3), and kHbA1cout (Table S4).",
      source_name        = "SEX (1 = female; 0 = male)"
    ),
    SMOKE_CURRENT = list(
      description        = "Current-smoker indicator (paired with SMOKE_NEVER)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = former smoker (when paired with SMOKE_NEVER = 0)",
      notes              = paste(
        "Baron 2016 Table S1 reports smoking as a 3-level categorical with NEVER smoker as",
        "the reference (CL multipliers theta_4 = 1.02 for exsmoker, theta_5 = 1.06 for",
        "current smoker). The canonical SMOKE_NEVER + SMOKE_CURRENT pairing uses FORMER",
        "smoker as the reference (both indicators = 0); the model recodes the published",
        "coefficients so the reference category becomes former smoker: see in-file comments",
        "on cl_smoke_never and cl_smoke_current."
      ),
      source_name        = "SMK_C (1 = current; 0 otherwise)"
    ),
    SMOKE_NEVER = list(
      description        = "Never-smoker indicator (paired with SMOKE_CURRENT)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = former smoker (when paired with SMOKE_CURRENT = 0)",
      notes              = "See SMOKE_CURRENT notes above for the recoding from Baron 2016's never-as-reference encoding to the canonical former-as-reference pairing.",
      source_name        = "SMK_N (1 = never; 0 otherwise)"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian)",
      notes              = "Multiplicative effects on CL, V2/F, V3/F, ka (Table S1), BFPG, Gmax (Table S3), and kHbA1cout (Table S4). Reference: non-Asian (pooled white + black + other).",
      source_name        = "ASIAN"
    ),
    RACE_BLACK = list(
      description        = "Black race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Black)",
      notes              = "Multiplicative effects on BFPG and Gmax only (Table S3). No PK or kHbA1cout effects.",
      source_name        = "BLACK"
    ),
    FPG = list(
      description        = "Observed baseline fasting plasma glucose at study entry",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-subject baseline FPG used as a CALIBRATION COVARIATE in the model-predicted BFPG equation (Baron 2016 Table S3): BFPG_i = IBFPG * exp(theta * 8 / FPG_obs) * (covariates) * exp(eta). Distinct from the model state glucose(t) which is the simulated time-varying FPG.",
      source_name        = "BFPG"
    ),
    T_DIAG_DIAB = list(
      description        = "Time since type 2 diabetes diagnosis at study entry",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effects centred at 2 years on BFPG, Gmax (Table S3) and kHbA1cout (Table S4). Baron 2016 categorises duration as <1 year / 1-5 years / >5 years (Table S2) but the structural model retains it as a continuous power-form covariate; supply the per-subject duration in years.",
      source_name        = "DUR"
    ),
    CONMED_METFORMIN = list(
      description        = "Concomitant metformin co-administration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant metformin)",
      notes              = "Multiplicative effects on BFPG, Gmax (Table S3) and kHbA1cout (Table S4).",
      source_name        = "MET"
    ),
    CONMED_SULFONYLUREA = list(
      description        = "Concomitant sulfonylurea co-administration indicator (class indicator)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant sulfonylurea)",
      notes              = "Multiplicative effects on BFPG and Gmax (Table S3); no kHbA1cout effect. Class indicator pooling glyburide / glibenclamide / glipizide / glimepiride / gliclazide / tolbutamide as a single second-generation oral antidiabetic stratum (Baron 2016 Methods / Table S2: 'background sulfonylurea' study-design cohort).",
      source_name        = "SU"
    ),
    CONMED_PIOGLITAZONE = list(
      description        = "Concomitant pioglitazone co-administration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant pioglitazone)",
      notes              = "Multiplicative effects on BFPG and Gmax (Table S3); no kHbA1cout effect.",
      source_name        = "PIO"
    ),
    DOSE_EMPA_MGD = list(
      description        = "Patient's own once-daily empagliflozin dose at the current dosing record",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Drives the steady-state AUC computation that feeds the FPG Emax: AUCss [nM*h] = DOSE_EMPA_MGD * 1e6 / MW_empa / cl, with MW_empa = 450.91 g/mol. Set to 0 for placebo arms (yields AUCss = 0 and STIM = 0). Studied doses: 1, 5, 10, 25, 50 mg QD (Table S2). Supply per dose record (constant within an inter-dose interval).",
      source_name        = "DOSE"
    )
  )

  population <- list(
    species             = "human",
    n_subjects_pk       = 2761L,
    n_subjects_pkpd     = 4065L,
    n_studies           = 10L,
    n_pk_observations   = 12503L,
    n_fpg_observations  = 25361L,
    n_hba1c_observations = 22012L,
    studies             = "2 phase I + 4 phase II + 4 phase III multi-region clinical trials (Baron 2016 Table 1).",
    age_range           = "36.0 - 76.0 years (PK/PD pooled dataset 2.5th - 97.5th percentile, Baron 2016 Table 2)",
    age_median          = "58.0 years",
    bmi_median          = "29.1 kg/m^2 (range 21.0 - 42.4)",
    sex_female_pct      = 45.5,
    race_ethnicity      = c(NonBlackAsian = 55.3, Black = 2.73, Asian = 42.0),
    smoking             = c(Never = 62.7, Former = 24.1, Current = 13.2),
    disease_state       = "Type 2 diabetes mellitus (T2DM); duration of T2DM > 5 years in 58.5% of patients (Baron 2016 Table S2).",
    dose_range          = "1 - 50 mg empagliflozin PO once daily; 40.9% on 10 mg QD and 46.0% on 25 mg QD (Baron 2016 Results).",
    baseline_fpg_median = "8.38 mmol/L (range 4.83 - 13.6)",
    baseline_hba1c_median = "7.90 % (range 6.70 - 9.80)",
    crcl_median         = "81.8 mL/min/1.73 m^2 (range 33.4 - 128)",
    tpro_median         = "72.0 g/dL (range 64.0 - 82.0)",
    co_medication       = "Metformin in 37.6%; metformin + sulfonylurea in 29.6%; no antidiabetic backbone in 32.8% (Baron 2016 Table S2).",
    notes               = "Demographics summarised in Baron 2016 Tables 1, 2, and S2. PK dataset excludes placebo subjects; PK/PD dataset includes placebo arms (n = 1469 / 36.1%)."
  )

  ini({
    # ---- Structural PK parameters (Baron 2016 Table S1; reference subject: 50-year-old, ----
    # ---- non-smoking male, non-Asian, BMI 25, eGFR 100, TPRO 70 g/dL, ALT 20, AST 20, ----
    # ---- AP 70, LDH 160) ----
    lcl <- log(10.6);   label("Apparent oral clearance CL/F (L/h) for the reference subject")  # Table S1 row TV_CL = 10.6 L/h
    lvc <- log(3.14);   label("Apparent central volume V2/F (L) for the reference subject; renamed V2 -> vc per canonical compartment-name conventions")     # Table S1 row TV_V2 = 3.14 L
    lq  <- log(6.34);   label("Apparent inter-compartmental clearance Q/F (L/h)")               # Table S1 row Q/F = 6.34 L/h
    lvp <- log(70.6);   label("Apparent peripheral volume V3/F (L) for the reference subject; renamed V3 -> vp per canonical compartment-name conventions")  # Table S1 row TV_V3 = 70.6 L
    lka <- log(0.196);  label("First-order absorption rate constant ka (1/h) for the reference subject")  # Table S1 row TV_ka = 0.196 1/h
    lalag <- fixed(log(0.500)); label("Absorption lag time ALAG1 (h; fixed per Baron 2016 Methods)")  # Table S1 row ALAG1 = 0.500 h; FIX per Methods "The absorption lag time was fixed to 0.5 h"

    # ---- PK covariate effects on CL/F (Table S1) ----
    # Power-form continuous effects exp(theta * log(cov / ref))
    e_age_cl  <- -0.241; label("Power exponent of (AGE/50) on CL/F")     # Table S1 theta_1 = -0.241
    e_bmi_cl  <-  0.320; label("Power exponent of (BMI/25) on CL/F")     # Table S1 theta_2 = 0.320
    e_egfr_cl <-  0.333; label("Power exponent of (CRCL/100) on CL/F")   # Table S1 theta_6 = 0.333
    e_tpro_cl <- -0.219; label("Power exponent of (TPRO/70 g/dL) on CL/F")  # Table S1 theta_8 = -0.219
    e_alt_cl  <-  0.0200;  label("Power exponent of (ALT/20) on CL/F")   # Table S1 theta_9 = 0.0200
    e_ast_cl  <- -0.0421;  label("Power exponent of (AST/20) on CL/F")   # Table S1 theta_10 = -0.0421
    e_ap_cl   <- -0.0488;  label("Power exponent of (ALP/70) on CL/F")   # Table S1 theta_11 = -0.0488
    e_ldh_cl  <-  0.0110;  label("Power exponent of (LDH/160) on CL/F")  # Table S1 theta_12 = 0.0110

    # Categorical effects on CL/F (multiplicative).
    # Smoking is recoded from Baron 2016's never-as-reference encoding (Table S1
    # theta_4^exsmoker = 1.02, theta_5^current smoker = 1.06) to the canonical
    # SMOKE_NEVER + SMOKE_CURRENT pair (former smoker as the implicit reference).
    cl_female      <- 0.886; label("Multiplicative effect of female sex on CL/F")  # Table S1 theta_3 = 0.886
    cl_asian       <- 0.880; label("Multiplicative effect of Asian race on CL/F")  # Table S1 theta_7 = 0.880
    cl_smoke_never   <- 0.9804; label("Multiplicative effect of never-smoker (vs former) on CL/F = 1 / theta_4_paper = 1 / 1.02")  # Recoded from Table S1 theta_4 = 1.02 (exsmoker vs never)
    cl_smoke_current <- 1.0392; label("Multiplicative effect of current-smoker (vs former) on CL/F = theta_5_paper / theta_4_paper = 1.06 / 1.02")  # Recoded from Table S1 theta_5 = 1.06 (current smoker vs never) and theta_4 = 1.02 (exsmoker vs never)

    # ---- PK covariate effects on V2/F (paper); renamed to V_central = vc per canonical conventions (Table S1) ----
    e_age_vc  <- 0.795; label("Power exponent of (AGE/50) on V2/F (now vc)")        # Table S1 theta_13 = 0.795
    e_tpro_vc <- 2.07;  label("Power exponent of (TPRO/70 g/dL) on V2/F (now vc)")  # Table S1 theta_15 = 2.07
    e_bmi_vc  <- 1.04;  label("Power exponent of (BMI/20) on V2/F (now vc); reference 20 kg/m^2 per Table S1")  # Table S1 theta_16 = 1.04, ref 20 (not 25)
    vc_female <- 1.25;  label("Multiplicative effect of female sex on V2/F (now vc)")  # Table S1 theta_14 = 1.25
    vc_asian  <- 1.27;  label("Multiplicative effect of Asian race on V2/F (now vc)")  # Table S1 theta_17 = 1.27

    # ---- PK covariate effects on V3/F (paper); renamed to V_peripheral = vp per canonical conventions (Table S1) ----
    e_age_vp  <-  0.135; label("Power exponent of (AGE/50) on V3/F (now vp)")        # Table S1 theta_18 = 0.135
    e_tpro_vp <- -0.196; label("Power exponent of (TPRO/70 g/dL) on V3/F (now vp)")  # Table S1 theta_20 = -0.196
    e_bmi_vp  <-  0.672; label("Power exponent of (BMI/25) on V3/F (now vp)")        # Table S1 theta_21 = 0.672
    vp_female <-  0.831; label("Multiplicative effect of female sex on V3/F (now vp)")  # Table S1 theta_19 = 0.831
    vp_asian  <-  0.959; label("Multiplicative effect of Asian race on V3/F (now vp)")  # Table S1 theta_22 = 0.959

    # ---- PK covariate effects on ka (Table S1) ----
    e_age_ka  <- 0.108; label("Power exponent of (AGE/50) on ka")  # Table S1 theta_23 = 0.108
    ka_female <- 1.17;  label("Multiplicative effect of female sex on ka")  # Table S1 theta_24 = 1.17
    ka_asian  <- 1.23;  label("Multiplicative effect of Asian race on ka")  # Table S1 theta_25 = 1.23

    # ---- PK inter-individual variability (Table S1; reported as omega^2 with CV% in adjacent column) ----
    # CL and V3 (peripheral, now `vp`) share a covariance; ka is independent.
    etalcl + etalvp ~ c(0.142,
                        0.0447, 0.0744)   # omega^2_CL = 0.142 (39.1% CV); cov(CL,V3) = 0.0447; omega^2_V3 = 0.0744 (27.8% CV)
    etalka ~ 0.0262                       # omega^2_ka = 0.0262 (16.3% CV)

    # ---- PK residual error (Table S1; proportional model) ----
    # Two phases were reported: Studies 1-2-5 (early phase, lower residual variability,
    # 16.9% CV) and Studies 3-4-6-7-8-9-10 (late phase, higher residual variability, 37.0% CV).
    # This implementation encodes the LATE-PHASE residual (37.0% CV) because it dominates
    # the dataset (Studies 6-10 are the registration / phase III trials and contributed the
    # majority of patients). The early-phase residual is documented in the vignette's
    # Assumptions and deviations section.
    # The Table S1 footnote describes a separate additive-on-log-scale residual added to
    # retain ~1.74% of outlier observations (CWRES 1 to 3 or > 3); this outlier residual
    # is not propagated to the simulation model.
    propSd <- sqrt(0.128); label("Proportional residual error for Cc (late-phase Studies 3-4-6-10, 37.0% CV; sqrt(omega^2 = 0.128))")  # Table S1 sigma^2_prop_Study3-4-6-10 = 0.128

    # ---- FPG model structural parameters (Table S3; reference: BFPG_obs 8 mM, AGE 50, ----
    # ---- BMI 25, eGFR 100, male, non-black non-Asian, duration 2 years, no concomitant therapy) ----
    libfpg <- log(14.2); label("Intrinsic typical baseline FPG (IBFPG) for Studies 6-10 (mmol/L)")  # Table S3 IBFPG (Studies 6-10) = 14.2
    e_obsbfpg_bfpg <- -0.497; label("Coefficient on (8 / FPG_observed) inside exp() for predicted BFPG")  # Table S3 theta_1 = -0.497

    # BFPG covariate effects
    e_age_bfpg   <- -0.100;   label("Power exponent of (AGE/50) on predicted BFPG")     # Table S3 theta_2 = -0.100
    e_bmi_bfpg   <-  0.0248;  label("Power exponent of (BMI/25) on predicted BFPG")     # Table S3 theta_3 = 0.0248
    e_egfr_bfpg  <- -0.0364;  label("Power exponent of (CRCL/100) on predicted BFPG")   # Table S3 theta_4 = -0.0364
    bfpg_female  <-  0.996;   label("Multiplicative effect of female sex on predicted BFPG")  # Table S3 theta_5 = 0.996
    bfpg_black   <-  0.977;   label("Multiplicative effect of Black race on predicted BFPG") # Table S3 theta_6 = 0.977
    bfpg_asian   <-  0.948;   label("Multiplicative effect of Asian race on predicted BFPG") # Table S3 theta_7 = 0.948
    e_dur_bfpg   <-  0.0512;  label("Power exponent of (T_DIAG_DIAB/2) on predicted BFPG")   # Table S3 theta_8 = 0.0512
    bfpg_met     <-  0.995;   label("Multiplicative effect of concomitant metformin on predicted BFPG")  # Table S3 theta_9 = 0.995
    bfpg_su      <-  1.01;    label("Multiplicative effect of concomitant sulfonylurea on predicted BFPG")  # Table S3 theta_10 = 1.01
    bfpg_pio     <-  0.999;   label("Multiplicative effect of concomitant pioglitazone on predicted BFPG")  # Table S3 theta_11 = 0.999

    # FPG turnover - kFPGout was held FIXED in the Baron 2016 final model (Table S3 RSE = N/A;
    # 95% CI = 0.0407, 0.0407, indicating fixed rather than estimated).
    lkfpgout <- fixed(log(0.0407)); label("First-order FPG elimination rate kFPGout (1/h; FIXED per Table S3)")  # Table S3 kFPGout = 0.0407 1/h; RSE N/A and bootstrap CI (0.0407, 0.0407) -> fixed

    # Gmax (maximum fractional reduction in FPG with drug)
    lgmax <- log(0.217); label("Typical maximal fractional FPG reduction Gmax (Emax)")  # Table S3 TV_Gmax = 0.217 (i.e. 21.7%)

    # Gmax covariate effects
    e_bfpg_gmax <-  1.76;    label("Power exponent of (BFPG_predicted/8) on Gmax")    # Table S3 theta_12 = 1.76
    e_age_gmax  <- -0.401;   label("Power exponent of (AGE/50) on Gmax")              # Table S3 theta_13 = -0.401
    e_bmi_gmax  <- -0.288;   label("Power exponent of (BMI/25) on Gmax")              # Table S3 theta_14 = -0.288
    e_egfr_gmax <-  0.512;   label("Power exponent of (CRCL/100) on Gmax")            # Table S3 theta_15 = 0.512
    gmax_female <-  0.890;   label("Multiplicative effect of female sex on Gmax")     # Table S3 theta_16 = 0.890
    gmax_black  <-  0.979;   label("Multiplicative effect of Black race on Gmax")     # Table S3 theta_17 = 0.979
    gmax_asian  <-  1.07;    label("Multiplicative effect of Asian race on Gmax")     # Table S3 theta_18 = 1.07
    e_dur_gmax  <-  0.0117;  label("Power exponent of (T_DIAG_DIAB/2) on Gmax")       # Table S3 theta_19 = 0.0117
    gmax_met    <-  0.931;   label("Multiplicative effect of concomitant metformin on Gmax")     # Table S3 theta_20 = 0.931
    gmax_su     <-  1.27;    label("Multiplicative effect of concomitant sulfonylurea on Gmax")  # Table S3 theta_21 = 1.27
    gmax_pio    <-  1.02;    label("Multiplicative effect of concomitant pioglitazone on Gmax")  # Table S3 theta_22 = 1.02

    # AUC50 (steady-state AUC at half-maximal Gmax)
    lauc50 <- log(703); label("AUC50: steady-state AUC at half of Gmax (nM*h)")  # Table S3 AUC50 = 703 nM*h

    # FPG IIV (Table S3)
    etalibfpg ~ 0.0165   # omega^2_BFPG = 0.0165 (12.9% CV)
    etalgmax  ~ 0.237    # omega^2_Gmax = 0.237 (51.7% CV)

    # FPG residual error (Table S3; proportional)
    propSd_glucose <- sqrt(0.0236); label("Proportional residual error for FPG (15.5% CV; sqrt(omega^2 = 0.0236))")  # Table S3 sigma^2_prop = 0.0236

    # ---- HbA1c model structural parameters (Table S4) ----
    # Baron 2016 Table S4 estimates only the ratio kHbA1cin / kHbA1cout (= 0.466) and
    # kHbA1cout itself; kHbA1cin is derived inside model() as kHbA1cin = ratio * kHbA1cout.
    hba1c_ratio <- 0.466; label("Ratio kHbA1cin / kHbA1cout (mmol/L glucose to % HbA1c production-elimination ratio)")  # Table S4 Ratio = 0.466
    hba1climit  <- 3.99;  label("HbA1c boundary condition HbA1climit, Studies 6-10 (%)")  # Table S4 HbA1climit (Studies 6-10) = 3.99 %
    lkhba1cout  <- log(1.59e-3); label("Typical first-order HbA1c elimination rate kHbA1cout (1/h)")  # Table S4 TV_kHbA1cout = 1.59e-3 1/h (~ half-life 2.6 weeks)

    # kHbA1cout covariate effects
    e_age_khba1cout    <- -0.241; label("Power exponent of (AGE/50) on kHbA1cout")     # Table S4 theta_1 = -0.241
    e_bmi_khba1cout    <-  0.328; label("Power exponent of (BMI/25) on kHbA1cout")     # Table S4 theta_2 = 0.328
    e_egfr_khba1cout   <- -0.119; label("Power exponent of (CRCL/100) on kHbA1cout")   # Table S4 theta_3 = -0.119
    khba1cout_female   <-  0.858; label("Multiplicative effect of female sex on kHbA1cout")  # Table S4 theta_4 = 0.858
    khba1cout_asian    <-  1.05;  label("Multiplicative effect of Asian race on kHbA1cout")  # Table S4 theta_5 = 1.05
    e_dur_khba1cout    <- -0.577; label("Power exponent of (T_DIAG_DIAB/2) on kHbA1cout")    # Table S4 theta_6 = -0.577
    khba1cout_met      <-  1.48;  label("Multiplicative effect of concomitant metformin on kHbA1cout")  # Table S4 theta_7 = 1.48

    # HbA1c IIV (Table S4)
    etalkhba1cout ~ 0.0242   # omega^2_kHbA1cout = 0.0242 (15.7% CV)

    # HbA1c residual error (Table S4; proportional)
    propSd_hba1c <- sqrt(0.00321); label("Proportional residual error for HbA1c (5.67% CV; sqrt(omega^2 = 0.00321))")  # Table S4 sigma^2_prop_HbA1c = 0.00321
  })

  model({
    # ---- Constants ----
    # Empagliflozin molecular weight; small molecule (BI 10773; SGLT2 inhibitor;
    # CAS 864070-44-0, molecular formula C23H27ClO7) -> 450.91 g/mol.
    mw_empa     <- 450.91
    nmol_per_mg <- 1e6 / mw_empa

    # ---- Population reference values for centred covariate effects (Baron 2016 Tables S1, S3, S4 footnotes) ----
    ref_age          <- 50      # years (Table S1 / S3 / S4 reference subject)
    ref_bmi_cl       <- 25      # kg/m^2 reference for CL, V3, BFPG, Gmax, kHbA1cout
    ref_bmi_v2       <- 20      # kg/m^2 reference for V2/F only (Table S1 row V2 ~ BMI)
    ref_egfr         <- 100     # mL/min/1.73 m^2
    ref_tpro_gdL     <- 70      # g/dL (Table S1 reference; covariate column is g/L canonical)
    ref_alt          <- 20      # U/L
    ref_ast          <- 20      # U/L
    ref_ap           <- 70      # U/L
    ref_ldh          <- 160     # U/L
    ref_bfpg_obs     <- 8       # mmol/L (Table S3 reference observed BFPG and BFPG/8 anchor for Gmax)
    ref_t_diag_diab  <- 2       # years (Table S3 / S4 reference duration anchor for the (DUR/2)^theta forms)

    # Convert the canonical g/L TPRO column to g/dL for use against the paper's g/dL reference
    tpro_gdL <- TPRO / 10

    # ---- Individual PK parameters (Baron 2016 Table S1; reference values defined above) ----
    cl <- exp(lcl + etalcl) *
          (AGE / ref_age)^e_age_cl *
          (BMI / ref_bmi_cl)^e_bmi_cl *
          (CRCL / ref_egfr)^e_egfr_cl *
          (tpro_gdL / ref_tpro_gdL)^e_tpro_cl *
          (ALT / ref_alt)^e_alt_cl *
          (AST / ref_ast)^e_ast_cl *
          (ALP / ref_ap)^e_ap_cl *
          (LDH / ref_ldh)^e_ldh_cl *
          cl_female^SEXF *
          cl_asian^RACE_ASIAN *
          cl_smoke_never^SMOKE_NEVER *
          cl_smoke_current^SMOKE_CURRENT

    vc <- exp(lvc) *
          (AGE / ref_age)^e_age_vc *
          (tpro_gdL / ref_tpro_gdL)^e_tpro_vc *
          (BMI / ref_bmi_v2)^e_bmi_vc *
          vc_female^SEXF *
          vc_asian^RACE_ASIAN

    q <- exp(lq)

    vp <- exp(lvp + etalvp) *
          (AGE / ref_age)^e_age_vp *
          (tpro_gdL / ref_tpro_gdL)^e_tpro_vp *
          (BMI / ref_bmi_cl)^e_bmi_vp *
          vp_female^SEXF *
          vp_asian^RACE_ASIAN

    ka <- exp(lka + etalka) *
          (AGE / ref_age)^e_age_ka *
          ka_female^SEXF *
          ka_asian^RACE_ASIAN

    # ---- Micro-rate constants ----
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- Steady-state AUC driver for the FPG PD model ----
    # AUCss [nmol*h/L = nM*h] = (DOSE_EMPA_MGD * 1e6 / MW_empa) [nmol/day] / cl [L/h]
    #                        (for QD dosing, dose-per-day matches the 24-h AUCss).
    AUCss <- DOSE_EMPA_MGD * nmol_per_mg / cl

    # ---- PK ODEs ----
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Absorption lag time (FIX 0.5 h per Methods)
    alag(depot) <- exp(lalag)

    # ---- Individual FPG / PD parameters ----
    # Predicted individual baseline FPG (BFPG) per Table S3:
    #   BFPG_i = IBFPG * exp(theta_1 * 8 / FPG_obs) * power-form continuous effects
    #           * multiplicative categorical effects * exp(eta).
    # The leading exp() term anchors the typical IBFPG (Studies 6-10, 14.2 mmol/L) to
    # an observed-FPG-of-8-mM reference; for an observed FPG of 8 mmol/L the term reduces
    # to exp(-0.497) = 0.608, so the typical-subject BFPG is ~ 14.2 * 0.608 = 8.6 mmol/L
    # (consistent with the pooled cohort observed median of 8.4 mmol/L; Baron 2016 Table 2).
    bfpg_indiv <- exp(libfpg + etalibfpg) *
                  exp(e_obsbfpg_bfpg * ref_bfpg_obs / FPG) *
                  (AGE / ref_age)^e_age_bfpg *
                  (BMI / ref_bmi_cl)^e_bmi_bfpg *
                  (CRCL / ref_egfr)^e_egfr_bfpg *
                  bfpg_female^SEXF *
                  bfpg_black^RACE_BLACK *
                  bfpg_asian^RACE_ASIAN *
                  (T_DIAG_DIAB / ref_t_diag_diab)^e_dur_bfpg *
                  bfpg_met^CONMED_METFORMIN *
                  bfpg_su^CONMED_SULFONYLUREA *
                  bfpg_pio^CONMED_PIOGLITAZONE

    kFPGout <- exp(lkfpgout)

    gmax_indiv <- exp(lgmax + etalgmax) *
                  (bfpg_indiv / ref_bfpg_obs)^e_bfpg_gmax *
                  (AGE / ref_age)^e_age_gmax *
                  (BMI / ref_bmi_cl)^e_bmi_gmax *
                  (CRCL / ref_egfr)^e_egfr_gmax *
                  gmax_female^SEXF *
                  gmax_black^RACE_BLACK *
                  gmax_asian^RACE_ASIAN *
                  (T_DIAG_DIAB / ref_t_diag_diab)^e_dur_gmax *
                  gmax_met^CONMED_METFORMIN *
                  gmax_su^CONMED_SULFONYLUREA *
                  gmax_pio^CONMED_PIOGLITAZONE

    auc50 <- exp(lauc50)

    # Steady-state Emax stimulation of FPG elimination (Baron 2016 Eq. for STIM)
    STIM <- gmax_indiv * AUCss / (auc50 + AUCss)

    # ---- FPG ODE (Baron 2016 Eq. 1) ----
    # d(FPG)/dt = kFPGin - kFPGout * FPG * (1 + STIM)
    # Steady-state baseline: kFPGin = bfpg_indiv * kFPGout
    kFPGin <- bfpg_indiv * kFPGout
    d/dt(glucose) <- kFPGin - kFPGout * glucose * (1 + STIM)
    glucose(0)    <- bfpg_indiv

    # ---- Individual HbA1c parameters and ODE (Baron 2016 Eq. 2) ----
    # d(HbA1c)/dt = kHbA1cin * FPG - kHbA1cout * HbA1c * (1 - HbA1climit/HbA1c)
    # Steady-state baseline: HbA1c_ss = HbA1climit + (kHbA1cin / kHbA1cout) * FPG
    kHbA1cout <- exp(lkhba1cout + etalkhba1cout) *
                 (AGE / ref_age)^e_age_khba1cout *
                 (BMI / ref_bmi_cl)^e_bmi_khba1cout *
                 (CRCL / ref_egfr)^e_egfr_khba1cout *
                 khba1cout_female^SEXF *
                 khba1cout_asian^RACE_ASIAN *
                 (T_DIAG_DIAB / ref_t_diag_diab)^e_dur_khba1cout *
                 khba1cout_met^CONMED_METFORMIN

    kHbA1cin <- hba1c_ratio * kHbA1cout
    hba1c_init <- hba1climit + hba1c_ratio * bfpg_indiv

    d/dt(hba1c) <- kHbA1cin * glucose - kHbA1cout * hba1c * (1 - hba1climit / hba1c)
    hba1c(0)    <- hba1c_init

    # ---- Observations and residual error ----
    # Empagliflozin plasma concentration in nmol/L (= nM); central state holds mg.
    # Cc [nM] = (central [mg] / vc [L]) * (1e6 / MW [g/mol])
    Cc <- central / vc * nmol_per_mg
    Cc      ~ prop(propSd)
    glucose ~ prop(propSd_glucose)
    hba1c   ~ prop(propSd_hba1c)
  })
}
