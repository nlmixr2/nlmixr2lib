Sokolov_2020_sglt_qsp <- function() {
  description <- paste(
    "QSP. Renal glucose filtration, SGLT2 / SGLT1-mediated reabsorption and",
    "urinary excretion in adults (healthy and type 2 diabetes mellitus)",
    "driven by the physiologically based PK of three SGLT2 inhibitors:",
    "dapagliflozin, empagliflozin and canagliflozin. Each gliflozin has an",
    "oral transit-chain absorption into a shared plasma volume (dapagliflozin",
    "with one peripheral compartment), non-renal plasma clearance, and",
    "glomerular filtration of the unbound fraction into the S1/S2 (proximal",
    "convoluted, SGLT2) and then S3 (proximal straight, SGLT1) tubular lumen,",
    "the bladder and the urine. Glucose is filtered at GFR x mean daily",
    "plasma glucose and reabsorbed by Michaelis-Menten SGLT2 (S1/S2) and SGLT1",
    "(S3) kinetics under competitive inhibition by the luminal free drug,",
    "with disease-specific Vmax and a T2DM factor on every Ki. Mean daily",
    "plasma glucose, eGFR and T2DM status enter as covariates.",
    "Typical-value model (fit to study-level mean data); no IIV and no",
    "residual error. 37 ODE states."
  )
  reference <- paste(
    "Sokolov V, Yakovleva T, Chu L, Tang W, Greasley PJ, Johansson S,",
    "Peskov K, Helmlinger G, Boulton DW, Penland RC.",
    "Differentiating the sodium-glucose cotransporter 1 inhibition capacity",
    "of canagliflozin vs. dapagliflozin and empagliflozin using quantitative",
    "systems pharmacology modeling.",
    "CPT Pharmacometrics Syst Pharmacol. 2020;9(4):222-229.",
    "doi:10.1002/psp4.12498.",
    "Model structure and parameter estimation first reported in Yakovleva T",
    "et al. Diabetes Obes Metab. 2019;21(12):2684-2693",
    "(doi:10.1111/dom.13858); all values here are taken from the Sokolov 2020",
    "deposited model code (Supplementary SGLT_model_code.txt).",
    sep = " "
  )
  vignette <- "Sokolov_2020_sglt_qsp"
  paper_specific_compartments <- c(
    "pct_dapagliflozin",
    "pst_dapagliflozin",
    "bladder_dapagliflozin",
    "pct_empagliflozin",
    "pst_empagliflozin",
    "bladder_empagliflozin",
    "pct_canagliflozin",
    "pst_canagliflozin",
    "bladder_canagliflozin",
    "glu_pct",
    "glu_pst",
    "glu_bladder",
    "glu_urine"
  )

  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL"
  )

  covariateData <- list(
    GLU = list(
      description = "Mean daily plasma glucose (MPG), held constant over the simulation",
      units = "mmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The model's plasma glucose Glucosepls equals the arm-level mean daily",
        "plasma glucose (IQRtools regressor MPGinitial); it is not a dynamic",
        "state and is not lowered by the drug. Filtered glucose load = GFR *",
        "GLU (mmol/h). Sokolov 2020 Table S2: drug-specific medians 7.8",
        "(dapagliflozin), 9.24 (empagliflozin), 10.43 (canagliflozin) mM;",
        "common T2DM median 9.33 mM. The code default is 5.5 mM. Supply one",
        "value per subject (a time-varying value is accepted and is carried",
        "forward between rows)."
      ),
      source_name = "MPGinitial"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate (eGFR), used as the absolute GFR",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Converted inside model() to the model's GFR in L/h as",
        "CRCL * 60 / 1000 (the deposited script plots GFR * 1000 / 60 as",
        "eGFR in mL/min/1.73 m^2), i.e. the BSA-normalised eGFR is used",
        "directly as the absolute filtration rate. Trials reporting only",
        "creatinine clearance were converted by eGFR = 0.9 * CrCl (Sokolov",
        "2020 Methods). Table S2 medians: 98.1 (dapagliflozin) and 100",
        "(empagliflozin, canagliflozin, common) mL/min/1.73 m^2. The code",
        "default GFR = 6.6 L/h corresponds to 110 mL/min/1.73 m^2. GFR drives",
        "both glucose filtration and the filtration of unbound drug."
      ),
      source_name = "GFR"
    ),
    DIS_DIAB = list(
      description = "Type 2 diabetes mellitus indicator (1 = T2DM, 0 = healthy)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy)",
      notes = paste(
        "Maps the deposited code's mutually exclusive switches as kt2d =",
        "DIS_DIAB and khs = 1 - DIS_DIAB. In T2DM the SGLT2 Vmax is 111.4",
        "(healthy 87.07) mmol/h, the SGLT1 Vmax is 140 - 111.4 = 28.6",
        "(healthy 105.6 - 87.07 = 18.53) mmol/h, and every SGLT1 and SGLT2 Ki",
        "is multiplied by f_ki_t2dm = 0.3047. All Sokolov 2020 simulations",
        "are for T2DM (DIS_DIAB = 1); the healthy branch reproduces the",
        "Yakovleva 2019 healthy-volunteer fit."
      ),
      source_name = "kt2d / khs"
    )
  )

  compartmentData <- list(
    depot_dapagliflozin = list(
      analyte = "dapagliflozin",
      units = "mmol",
      specimen = "administration site",
      verified = TRUE
    ),
    transit1_dapagliflozin = list(
      analyte = "dapagliflozin",
      units = "mmol",
      specimen = "administration site",
      verified = TRUE
    ),
    transit2_dapagliflozin = list(
      analyte = "dapagliflozin",
      units = "mmol",
      specimen = "administration site",
      verified = TRUE
    ),
    transit3_dapagliflozin = list(
      analyte = "dapagliflozin",
      units = "mmol",
      specimen = "administration site",
      verified = TRUE
    ),
    transit4_dapagliflozin = list(
      analyte = "dapagliflozin",
      units = "mmol",
      specimen = "administration site",
      verified = TRUE
    ),
    transit5_dapagliflozin = list(
      analyte = "dapagliflozin",
      units = "mmol",
      specimen = "administration site",
      verified = TRUE
    ),
    central_dapagliflozin = list(analyte = "dapagliflozin", units = "mmol", specimen = "plasma", verified = TRUE),
    peripheral1_dapagliflozin = list(analyte = "dapagliflozin", units = "mmol", specimen = "tissue", verified = TRUE),
    pct_dapagliflozin = list(analyte = "dapagliflozin", units = "mmol", specimen = "urine", verified = TRUE),
    pst_dapagliflozin = list(analyte = "dapagliflozin", units = "mmol", specimen = "urine", verified = TRUE),
    bladder_dapagliflozin = list(analyte = "dapagliflozin", units = "mmol", specimen = "urine", verified = TRUE),
    urine_dapagliflozin = list(analyte = "dapagliflozin", units = "mmol", specimen = "urine", verified = TRUE),
    depot_empagliflozin = list(
      analyte = "empagliflozin",
      units = "mmol",
      specimen = "administration site",
      verified = TRUE
    ),
    transit1_empagliflozin = list(
      analyte = "empagliflozin",
      units = "mmol",
      specimen = "administration site",
      verified = TRUE
    ),
    transit2_empagliflozin = list(
      analyte = "empagliflozin",
      units = "mmol",
      specimen = "administration site",
      verified = TRUE
    ),
    transit3_empagliflozin = list(
      analyte = "empagliflozin",
      units = "mmol",
      specimen = "administration site",
      verified = TRUE
    ),
    transit4_empagliflozin = list(
      analyte = "empagliflozin",
      units = "mmol",
      specimen = "administration site",
      verified = TRUE
    ),
    transit5_empagliflozin = list(
      analyte = "empagliflozin",
      units = "mmol",
      specimen = "administration site",
      verified = TRUE
    ),
    central_empagliflozin = list(analyte = "empagliflozin", units = "mmol", specimen = "plasma", verified = TRUE),
    pct_empagliflozin = list(analyte = "empagliflozin", units = "mmol", specimen = "urine", verified = TRUE),
    pst_empagliflozin = list(analyte = "empagliflozin", units = "mmol", specimen = "urine", verified = TRUE),
    bladder_empagliflozin = list(analyte = "empagliflozin", units = "mmol", specimen = "urine", verified = TRUE),
    urine_empagliflozin = list(analyte = "empagliflozin", units = "mmol", specimen = "urine", verified = TRUE),
    depot_canagliflozin = list(
      analyte = "canagliflozin",
      units = "mmol",
      specimen = "administration site",
      verified = TRUE
    ),
    transit1_canagliflozin = list(
      analyte = "canagliflozin",
      units = "mmol",
      specimen = "administration site",
      verified = TRUE
    ),
    transit2_canagliflozin = list(
      analyte = "canagliflozin",
      units = "mmol",
      specimen = "administration site",
      verified = TRUE
    ),
    transit3_canagliflozin = list(
      analyte = "canagliflozin",
      units = "mmol",
      specimen = "administration site",
      verified = TRUE
    ),
    transit4_canagliflozin = list(
      analyte = "canagliflozin",
      units = "mmol",
      specimen = "administration site",
      verified = TRUE
    ),
    central_canagliflozin = list(analyte = "canagliflozin", units = "mmol", specimen = "plasma", verified = TRUE),
    pct_canagliflozin = list(analyte = "canagliflozin", units = "mmol", specimen = "urine", verified = TRUE),
    pst_canagliflozin = list(analyte = "canagliflozin", units = "mmol", specimen = "urine", verified = TRUE),
    bladder_canagliflozin = list(analyte = "canagliflozin", units = "mmol", specimen = "urine", verified = TRUE),
    urine_canagliflozin = list(analyte = "canagliflozin", units = "mmol", specimen = "urine", verified = TRUE),
    glu_pct = list(analyte = "glucose", units = "mmol", specimen = "urine", verified = TRUE),
    glu_pst = list(analyte = "glucose", units = "mmol", specimen = "urine", verified = TRUE),
    glu_bladder = list(analyte = "glucose", units = "mmol", specimen = "urine", verified = TRUE),
    glu_urine = list(analyte = "glucose", units = "mmol", specimen = "urine", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = NA_integer_,
    n_arms_t2dm = c(dapagliflozin = 7L, canagliflozin = 11L, empagliflozin = 14L),
    age_range = NA_character_,
    weight_range = NA_character_,
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state = paste(
      "Healthy volunteers and adults with type 2 diabetes mellitus from",
      "published phase I / II dose-ranging trials of dapagliflozin,",
      "empagliflozin and canagliflozin (Sokolov 2020 Table S1)."
    ),
    dose_range = paste(
      "Single and repeated oral doses: dapagliflozin 2.5-100 mg,",
      "canagliflozin 25-400 mg, empagliflozin 1-100 mg (T2DM",
      "pharmacodynamic arms)."
    ),
    regions = NA_character_,
    renal_function = paste(
      "Median eGFR 98.1 (dapagliflozin trials) and 100 (empagliflozin and",
      "canagliflozin trials) mL/min/1.73 m^2 (Sokolov 2020 Table S2)."
    ),
    notes = paste(
      "The model was fit to study-level MEAN data (plasma PK profiles,",
      "24-h urinary drug excretion and 24-h urinary glucose excretion per",
      "arm), not individual data, so subject counts and demographics are",
      "not reported. Arm-level mean daily plasma glucose and eGFR were",
      "supplied as regressors. Of the 44 model parameters, 27 were taken",
      "from the literature and 17 were estimated (Sokolov 2020 Methods;",
      "the deposited Script_Fitting_Simulation.R re-estimates the six",
      "SGLT parameters)."
    )
  )

  ini({
    # All values are from the Sokolov 2020 deposited model code
    # (Supplementary SGLT_model_code.txt, 'MODEL PARAMETERS' section), cited
    # below as 'Code'. The ODEs and rate laws are Sokolov 2020 Table S3.
    # Literature-defined constants (27 of the 44 parameters) are fixed();
    # the 17 estimated parameters are free.

    # ---- Physiology (literature constants) ----
    lvc <- fixed(log(2.75)); label("Plasma volume Vpl shared as the central volume of all three drugs (L)") # Code Vpl = 2.75
    v_pct <- fixed(0.045); label("S1/S2 proximal convoluted tubule lumen volume Vlumen1 (L)") # Code Vlumen1 = 0.045
    v_pst <- fixed(0.01944); label("S3 proximal straight tubule lumen volume Vlumen2 (L)") # Code Vlumen2 = 0.01944
    v_bladder <- fixed(0.2); label("Bladder volume Vbladder (L)") # Code Vbladder = 0.2
    q_pct <- fixed(2.7); label("Tubular fluid flow from S1/S2 to S3, Qlumen (L/h)") # Code Qlumen = 2.7
    q_pst <- fixed(0.72); label("Tubular fluid flow from S3 to the bladder, Qbladder (L/h)") # Code Qbladder = 0.72
    q_urine <- fixed(0.055); label("Urine flow from the bladder, Qurine (L/h)") # Code Qurine = 0.055

    # ---- Dapagliflozin PK ----
    lfdepot_dapagliflozin <- fixed(log(0.78)); label("Dapagliflozin oral bioavailability Fdapa (fraction)") # Code Fdapa = 0.78
    fu_dapagliflozin <- fixed(0.086); label("Dapagliflozin unbound fraction in plasma fupdapa (fraction)") # Code fupdapa = 0.086
    mw_dapagliflozin <- fixed(408.87); label("Dapagliflozin molecular weight MWdapa (g/mol)") # Code MWdapa = 408.87
    lktr_dapagliflozin <- log(10.38); label("Dapagliflozin transit rate constant k_tr_d (1/h)") # Code k_tr_d = 10.38
    lka_dapagliflozin <- log(0.5071); label("Dapagliflozin absorption rate constant kabsdapa (1/h)") # Code kabsdapa = 0.5071
    lcl_dapagliflozin <- log(13.07); label("Dapagliflozin non-renal plasma clearance CLdapapls (L/h)") # Code CLdapapls = 13.07
    lq_dapagliflozin <- log(13.65); label("Dapagliflozin intercompartmental clearance Qdapa (L/h)") # Code Qdapa = 13.65
    lvp_dapagliflozin <- log(98.57); label("Dapagliflozin peripheral volume Vprfdapa (L)") # Code Vprfdapa = 98.57

    # ---- Empagliflozin PK ----
    lfdepot_empagliflozin <- fixed(log(0.78)); label("Empagliflozin oral bioavailability Fempa (fraction)") # Code Fempa = 0.78
    fu_empagliflozin <- fixed(0.22); label("Empagliflozin unbound fraction in plasma fupempa (fraction)") # Code fupempa = 0.22
    mw_empagliflozin <- fixed(450.91); label("Empagliflozin molecular weight MWempa (g/mol)") # Code MWempa = 450.91
    lktr_empagliflozin <- log(19.8); label("Empagliflozin transit rate constant k_tr_e (1/h)") # Code k_tr_e = 19.8
    lka_empagliflozin <- log(0.09984); label("Empagliflozin absorption rate constant kabsempa (1/h)") # Code kabsempa = 0.09984
    lcl_empagliflozin <- log(5.78); label("Empagliflozin non-renal plasma clearance CLempapls (L/h)") # Code CLempapls = 5.78

    # ---- Canagliflozin PK ----
    lfdepot_canagliflozin <- fixed(log(0.65)); label("Canagliflozin oral bioavailability Fcana (fraction)") # Code Fcana = 0.65
    fu_canagliflozin <- fixed(0.01); label("Canagliflozin unbound fraction in plasma fupcana (fraction)") # Code fupcana = 0.01
    mw_canagliflozin <- fixed(444.5); label("Canagliflozin molecular weight MWcana (g/mol)") # Code MWcana = 444.5
    lktr_canagliflozin <- log(11.93); label("Canagliflozin transit rate constant k_tr_c (1/h)") # Code k_tr_c = 11.93
    lka_canagliflozin <- log(0.1358); label("Canagliflozin absorption rate constant kabscana (1/h)") # Code kabscana = 0.1358
    lcl_canagliflozin <- log(8.242); label("Canagliflozin non-renal plasma clearance CLcanapls (L/h)") # Code CLcanapls = 8.242

    # ---- Glucose handling (literature constants) ----
    mw_glucose <- fixed(180.156); label("Glucose molecular weight MWglucose (g/mol)") # Code MWglucose = 180.156
    km_sglt1 <- fixed(0.5); label("SGLT1 glucose Michaelis-Menten constant KmreabsSGLT1 (mmol/L)") # Code KmreabsSGLT1 = 0.5
    km_sglt2 <- fixed(4); label("SGLT2 glucose Michaelis-Menten constant KmreabsSGLT2 (mmol/L)") # Code KmreabsSGLT2 = 4
    vmax_total_healthy <- fixed(105.6); label("Total SGLT1 + SGLT2 maximum reabsorption rate in healthy subjects Vmax_hs (mmol/h)") # Code Vmax_hs = 105.6
    vmax_total_t2dm <- fixed(140); label("Total SGLT1 + SGLT2 maximum reabsorption rate in T2DM Vmax_t2d (mmol/h)") # Code Vmax_t2d = 140
    reabs_base <- fixed(39.174); label("Reference total glucose reabsorption rate for the inhibition readout reabsbase (mmol/h)") # Code reabsbase = 39.174
    ratio_sglt1_dapagliflozin <- fixed(1157); label("In vitro SGLT1:SGLT2 Ki ratio for dapagliflozin (unitless)") # Code KidapaSGLT1ex = k_sglt1*KidapaSGLT2ex*1157/1000
    ratio_sglt1_empagliflozin <- fixed(1249); label("In vitro SGLT1:SGLT2 Ki ratio for empagliflozin (unitless)") # Code KiempaSGLT1ex = k_sglt1*KiempaSGLT2ex*1249/1000
    ratio_sglt1_canagliflozin <- fixed(158); label("In vitro SGLT1:SGLT2 Ki ratio for canagliflozin (unitless)") # Code KicanaSGLT1ex = k_sglt1*KicanaSGLT2ex*158/1000

    # ---- Glucose handling (estimated) ----
    lvmax_sglt2_healthy <- log(87.07); label("SGLT2 maximum reabsorption rate in healthy subjects VmaxreabsSGLT2hs (mmol/h)") # Code VmaxreabsSGLT2hs = 87.07
    lvmax_sglt2_t2dm <- log(111.4); label("SGLT2 maximum reabsorption rate in T2DM VmaxreabsSGLT2t2d (mmol/h)") # Code VmaxreabsSGLT2t2d = 111.4
    lki_sglt2_dapagliflozin <- log(103.1); label("Dapagliflozin SGLT2 Ki in healthy subjects KidapaSGLT2ex (pmol/L)") # Code KidapaSGLT2ex = 103.1 pM
    lki_sglt2_empagliflozin <- log(638.6); label("Empagliflozin SGLT2 Ki in healthy subjects KiempaSGLT2ex (pmol/L)") # Code KiempaSGLT2ex = 638.6 pM
    lki_sglt2_canagliflozin <- log(364.3); label("Canagliflozin SGLT2 Ki in healthy subjects KicanaSGLT2ex (pmol/L)") # Code KicanaSGLT2ex = 364.3 pM
    f_ki_t2dm <- 0.3047; label("Multiplier on every SGLT1 and SGLT2 Ki in T2DM, coef (unitless)") # Code coef = 0.3047
  })

  model({
    # ---- Covariate-derived inputs ----
    # GFR in L/h from eGFR in mL/min/1.73 m^2 (deposited script, Figure 1 code).
    gfr <- CRCL * 60 / 1000
    glu_plasma <- GLU
    khs <- 1 - DIS_DIAB
    kt2d <- DIS_DIAB

    # ---- Individual parameters ----
    vc <- exp(lvc)
    fdepot_dapagliflozin <- exp(lfdepot_dapagliflozin)
    fdepot_empagliflozin <- exp(lfdepot_empagliflozin)
    fdepot_canagliflozin <- exp(lfdepot_canagliflozin)
    ktr_dapagliflozin <- exp(lktr_dapagliflozin)
    ka_dapagliflozin <- exp(lka_dapagliflozin)
    cl_dapagliflozin <- exp(lcl_dapagliflozin)
    q_dapagliflozin <- exp(lq_dapagliflozin)
    vp_dapagliflozin <- exp(lvp_dapagliflozin)
    ktr_empagliflozin <- exp(lktr_empagliflozin)
    ka_empagliflozin <- exp(lka_empagliflozin)
    cl_empagliflozin <- exp(lcl_empagliflozin)
    ktr_canagliflozin <- exp(lktr_canagliflozin)
    ka_canagliflozin <- exp(lka_canagliflozin)
    cl_canagliflozin <- exp(lcl_canagliflozin)

    # SGLT capacities: SGLT1 Vmax is the total minus the SGLT2 Vmax.
    vmax_sglt2 <- exp(lvmax_sglt2_healthy) * khs + exp(lvmax_sglt2_t2dm) * kt2d
    vmax_sglt1 <- (vmax_total_healthy - exp(lvmax_sglt2_healthy)) * khs +
      (vmax_total_t2dm - exp(lvmax_sglt2_t2dm)) * kt2d

    # Ki (mmol/L). SGLT2 Ki are estimated in pM (divide by 1e9); SGLT1 Ki are
    # the SGLT2 Ki times the in vitro ratio (nM = pM * ratio / 1000; divide
    # by 1e6). In T2DM every Ki is multiplied by f_ki_t2dm.
    ki_scale <- khs + kt2d * f_ki_t2dm
    ki_sglt2_dapagliflozin <- ki_scale * exp(lki_sglt2_dapagliflozin) / 1e9
    ki_sglt2_empagliflozin <- ki_scale * exp(lki_sglt2_empagliflozin) / 1e9
    ki_sglt2_canagliflozin <- ki_scale * exp(lki_sglt2_canagliflozin) / 1e9
    ki_sglt1_dapagliflozin <- ki_scale * exp(lki_sglt2_dapagliflozin) * ratio_sglt1_dapagliflozin / 1000 / 1e6
    ki_sglt1_empagliflozin <- ki_scale * exp(lki_sglt2_empagliflozin) * ratio_sglt1_empagliflozin / 1000 / 1e6
    ki_sglt1_canagliflozin <- ki_scale * exp(lki_sglt2_canagliflozin) * ratio_sglt1_canagliflozin / 1000 / 1e6

    # ---- Concentrations (mmol/L) ----
    cp_dapagliflozin <- central_dapagliflozin / vc
    cpt_dapagliflozin <- peripheral1_dapagliflozin / vp_dapagliflozin
    cpct_dapagliflozin <- pct_dapagliflozin / v_pct
    cpst_dapagliflozin <- pst_dapagliflozin / v_pst
    cbl_dapagliflozin <- bladder_dapagliflozin / v_bladder
    cp_empagliflozin <- central_empagliflozin / vc
    cpct_empagliflozin <- pct_empagliflozin / v_pct
    cpst_empagliflozin <- pst_empagliflozin / v_pst
    cbl_empagliflozin <- bladder_empagliflozin / v_bladder
    cp_canagliflozin <- central_canagliflozin / vc
    cpct_canagliflozin <- pct_canagliflozin / v_pct
    cpst_canagliflozin <- pst_canagliflozin / v_pst
    cbl_canagliflozin <- bladder_canagliflozin / v_bladder
    cglu_pct <- glu_pct / v_pct
    cglu_pst <- glu_pst / v_pst
    cglu_bladder <- glu_bladder / v_bladder

    # ---- Glucose reabsorption (Table S3; competitive inhibition, mmol/h) ----
    reabs_sglt2 <- cglu_pct * vmax_sglt2 / (km_sglt2 * (1 +
      cpct_dapagliflozin / ki_sglt2_dapagliflozin +
      cpct_empagliflozin / ki_sglt2_empagliflozin +
      cpct_canagliflozin / ki_sglt2_canagliflozin) + cglu_pct)
    reabs_sglt1 <- cglu_pst * vmax_sglt1 / (km_sglt1 * (1 +
      cpst_dapagliflozin / ki_sglt1_dapagliflozin +
      cpst_empagliflozin / ki_sglt1_empagliflozin +
      cpst_canagliflozin / ki_sglt1_canagliflozin) + cglu_pst)

    # ---- Dapagliflozin: depot + 5 transit states, 2-compartment plasma ----
    d/dt(depot_dapagliflozin) <- -ktr_dapagliflozin * depot_dapagliflozin
    d/dt(transit1_dapagliflozin) <- ktr_dapagliflozin * depot_dapagliflozin - ktr_dapagliflozin * transit1_dapagliflozin
    d/dt(transit2_dapagliflozin) <- ktr_dapagliflozin * transit1_dapagliflozin - ktr_dapagliflozin * transit2_dapagliflozin
    d/dt(transit3_dapagliflozin) <- ktr_dapagliflozin * transit2_dapagliflozin - ktr_dapagliflozin * transit3_dapagliflozin
    d/dt(transit4_dapagliflozin) <- ktr_dapagliflozin * transit3_dapagliflozin - ktr_dapagliflozin * transit4_dapagliflozin
    d/dt(transit5_dapagliflozin) <- ktr_dapagliflozin * transit4_dapagliflozin - ka_dapagliflozin * transit5_dapagliflozin
    d/dt(central_dapagliflozin) <- ka_dapagliflozin * transit5_dapagliflozin -
      cl_dapagliflozin * cp_dapagliflozin -
      q_dapagliflozin * (cp_dapagliflozin - cpt_dapagliflozin) -
      gfr * fu_dapagliflozin * cp_dapagliflozin
    d/dt(peripheral1_dapagliflozin) <- q_dapagliflozin * (cp_dapagliflozin - cpt_dapagliflozin)
    d/dt(pct_dapagliflozin) <- gfr * fu_dapagliflozin * cp_dapagliflozin - q_pct * cpct_dapagliflozin
    d/dt(pst_dapagliflozin) <- q_pct * cpct_dapagliflozin - q_pst * cpst_dapagliflozin
    d/dt(bladder_dapagliflozin) <- q_pst * cpst_dapagliflozin - q_urine * cbl_dapagliflozin
    d/dt(urine_dapagliflozin) <- q_urine * cbl_dapagliflozin

    # ---- Empagliflozin: depot + 5 transit states, 1-compartment plasma ----
    d/dt(depot_empagliflozin) <- -ktr_empagliflozin * depot_empagliflozin
    d/dt(transit1_empagliflozin) <- ktr_empagliflozin * depot_empagliflozin - ktr_empagliflozin * transit1_empagliflozin
    d/dt(transit2_empagliflozin) <- ktr_empagliflozin * transit1_empagliflozin - ktr_empagliflozin * transit2_empagliflozin
    d/dt(transit3_empagliflozin) <- ktr_empagliflozin * transit2_empagliflozin - ktr_empagliflozin * transit3_empagliflozin
    d/dt(transit4_empagliflozin) <- ktr_empagliflozin * transit3_empagliflozin - ktr_empagliflozin * transit4_empagliflozin
    d/dt(transit5_empagliflozin) <- ktr_empagliflozin * transit4_empagliflozin - ka_empagliflozin * transit5_empagliflozin
    d/dt(central_empagliflozin) <- ka_empagliflozin * transit5_empagliflozin -
      cl_empagliflozin * cp_empagliflozin -
      gfr * fu_empagliflozin * cp_empagliflozin
    d/dt(pct_empagliflozin) <- gfr * fu_empagliflozin * cp_empagliflozin - q_pct * cpct_empagliflozin
    d/dt(pst_empagliflozin) <- q_pct * cpct_empagliflozin - q_pst * cpst_empagliflozin
    d/dt(bladder_empagliflozin) <- q_pst * cpst_empagliflozin - q_urine * cbl_empagliflozin
    d/dt(urine_empagliflozin) <- q_urine * cbl_empagliflozin

    # ---- Canagliflozin: depot + 4 transit states, 1-compartment plasma ----
    d/dt(depot_canagliflozin) <- -ktr_canagliflozin * depot_canagliflozin
    d/dt(transit1_canagliflozin) <- ktr_canagliflozin * depot_canagliflozin - ktr_canagliflozin * transit1_canagliflozin
    d/dt(transit2_canagliflozin) <- ktr_canagliflozin * transit1_canagliflozin - ktr_canagliflozin * transit2_canagliflozin
    d/dt(transit3_canagliflozin) <- ktr_canagliflozin * transit2_canagliflozin - ktr_canagliflozin * transit3_canagliflozin
    d/dt(transit4_canagliflozin) <- ktr_canagliflozin * transit3_canagliflozin - ka_canagliflozin * transit4_canagliflozin
    d/dt(central_canagliflozin) <- ka_canagliflozin * transit4_canagliflozin -
      cl_canagliflozin * cp_canagliflozin -
      gfr * fu_canagliflozin * cp_canagliflozin
    d/dt(pct_canagliflozin) <- gfr * fu_canagliflozin * cp_canagliflozin - q_pct * cpct_canagliflozin
    d/dt(pst_canagliflozin) <- q_pct * cpct_canagliflozin - q_pst * cpst_canagliflozin
    d/dt(bladder_canagliflozin) <- q_pst * cpst_canagliflozin - q_urine * cbl_canagliflozin
    d/dt(urine_canagliflozin) <- q_urine * cbl_canagliflozin

    # ---- Glucose: filtration, reabsorption, bladder, urine (mmol) ----
    d/dt(glu_pct) <- gfr * glu_plasma - reabs_sglt2 - q_pct * cglu_pct
    d/dt(glu_pst) <- q_pct * cglu_pct - reabs_sglt1 - q_pst * cglu_pst
    d/dt(glu_bladder) <- q_pst * cglu_pst - q_urine * cglu_bladder
    d/dt(glu_urine) <- q_urine * cglu_bladder

    # ---- Drug-free glucose steady state as the initial condition ----
    # The deposited code starts the glucose states at zero and runs a 72-h
    # burn-in before dosing. Here each tubular segment starts at the
    # positive root of its drug-free steady-state mass balance,
    # inflow = Vmax*C/(Km + C) + Qout*C, so a simulation can be dosed at
    # time 0. The bladder follows from q_pst * C_pst = q_urine * C_bladder.
    gin_pct <- gfr * glu_plasma
    b_pct <- q_pct * km_sglt2 + vmax_sglt2 - gin_pct
    css_pct <- (-b_pct + sqrt(b_pct^2 + 4 * q_pct * gin_pct * km_sglt2)) / (2 * q_pct)
    gin_pst <- q_pct * css_pct
    b_pst <- q_pst * km_sglt1 + vmax_sglt1 - gin_pst
    css_pst <- (-b_pst + sqrt(b_pst^2 + 4 * q_pst * gin_pst * km_sglt1)) / (2 * q_pst)
    glu_pct(0) <- css_pct * v_pct
    glu_pst(0) <- css_pst * v_pst
    glu_bladder(0) <- q_pst * css_pst / q_urine * v_bladder

    # ---- Bioavailability: oral dose in mg -> absorbed amount in mmol ----
    f(depot_dapagliflozin) <- fdepot_dapagliflozin / mw_dapagliflozin
    f(depot_empagliflozin) <- fdepot_empagliflozin / mw_empagliflozin
    f(depot_canagliflozin) <- fdepot_canagliflozin / mw_canagliflozin

    # ---- Outputs ----
    # Plasma concentrations in ng/mL (mmol/L * g/mol * 1000).
    Cc_dapagliflozin <- cp_dapagliflozin * mw_dapagliflozin * 1000
    Cc_empagliflozin <- cp_empagliflozin * mw_empagliflozin * 1000
    Cc_canagliflozin <- cp_canagliflozin * mw_canagliflozin * 1000
    # Cumulative drug excreted in urine (mg).
    ae_dapagliflozin <- urine_dapagliflozin * mw_dapagliflozin
    ae_empagliflozin <- urine_empagliflozin * mw_empagliflozin
    ae_canagliflozin <- urine_canagliflozin * mw_canagliflozin
    # Cumulative urinary glucose excretion since time 0 (g).
    uge <- glu_urine * mw_glucose / 1000
    # Reabsorption fluxes (mmol/h) and total-inhibition readout (percent).
    reabs_total <- reabs_sglt2 + reabs_sglt1
    inh_total_pct <- 100 * (1 - reabs_total / reabs_base)
  })
}
