Ooi_2023_bevurogant <- function() {
  description <- paste(
    "Two-compartment population PK model for oral bevurogant (BI 730357, a",
    "RORgammat antagonist) in healthy adults and in patients with",
    "moderate-to-severe plaque psoriasis (Ooi 2023). Dual absorption: most of",
    "the bioavailable dose enters the depot as a zero-order input after a short",
    "lag and is then absorbed first-order, and the remaining fraction enters",
    "central as a zero-order input after a long (about 11.5 h) lag. Absolute",
    "oral bioavailability on the logit scale is study-specific and falls",
    "linearly with dose, with shifts for food, the intended commercial",
    "formulation and twice-daily dosing. Clearance is lower in patients and",
    "falls with baseline CRP and AST; allometric weight scaling. A separate",
    "pair of IV states carries the 100 ug intravenous microtracer of study",
    "1407-0033, measured as its own analyte (ROUTE_ORAL = 0 selects it on an",
    "observation record; Cc_iv is always output). Combined residual error",
    "per study / diet stratum with between-subject variability on the",
    "proportional part, the additive part scaled by the assay LLOQ, and a",
    "larger residual for trough samples."
  )
  reference <- paste(
    "Ooi QX, Wenning L, Plock N, Hoch M, Kerbusch T, Jacobs A, Zheng W, Ji Q,",
    "Gong J, Shah N, Joshi A, Jeon S, Sachs M, Kemmler M, Seeburger J,",
    "Hilger E, Pfeifer E. Bounded integer model-based analysis of psoriasis",
    "area and severity index in patients with moderate-to-severe plaque",
    "psoriasis receiving BI 730357. CPT Pharmacometrics Syst Pharmacol.",
    "2023;12(6):758-769. doi:10.1002/psp4.12948.",
    "Model equations from the Appendix S1 NONMEM control stream (run42).",
    sep = " "
  )
  vignette <- "Ooi_2023_bevurogant"
  units <- list(time = "h", dosing = "mg", concentration = "nmol/L")

  paper_specific_compartments <- c("central_iv", "peripheral1_iv")
  paper_specific_etas <- c("etaRUV")

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling referenced to 70 kg with fixed exponents 0.75 on CL and Q and 1 on Vc and Vp (Table 2, footnote b).",
      source_name = "WT"
    ),
    CRP = list(
      description = "Baseline C-reactive protein",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL centred on 0.305 mg/L, the median in healthy subjects (Table 2, footnote c; control stream 'EXP(THETA(25)*LOG(BCRP/0.305))'). The control stream sets the factor to 1 when CRP is missing (-99); impute a missing value with 0.305 mg/L to reproduce that.",
      source_name = "BCRP"
    ),
    AST = list(
      description = "Baseline aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL centred on 13 U/L, the median in healthy subjects (control stream 'EXP(THETA(28)*LOG(BAST/13))'). The control stream sets the factor to 1 when AST is missing (-99); impute a missing value with 13 U/L to reproduce that.",
      source_name = "BAST"
    ),
    DIS_PSORIASIS = list(
      description = "Moderate-to-severe plaque psoriasis patient indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy subject, studies 1407-0002, 1407-0032 and 1407-0033)",
      notes = "The control stream keys the patient effects on the study number: SWTPSO = 1 when PROJ == 14070030, the only patient study. Selects the patient CL multiplier 0.773, the patient dose slope on logit F, the patient residual-error stratum and, with all STUDY_1407_* indicators 0, the study 1407-0030 bioavailability.",
      source_name = "PROJ == 14070030"
    ),
    STUDY_1407_0002 = list(
      description = "Phase Ib multiple-dose study 1407-0002 in healthy men",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (another study; with DIS_PSORIASIS = 1, the patient study 1407-0030)",
      notes = "Selects the study-specific oral bioavailability 0.652 and, when fed, the fed healthy-subject dose slope on logit F. The ka study effect of this study is fixed to 0 in the control stream (THETA(32) 0 FIX), so it has no ka term.",
      source_name = "PROJ == 14070002"
    ),
    STUDY_1407_0032 = list(
      description = "Phase I single-dose study 1407-0032 in healthy subjects",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (another study)",
      notes = "Selects the study-specific oral bioavailability 0.429 and the multiplicative ka shift 1 + 0.575.",
      source_name = "PROJ == 14070032"
    ),
    STUDY_1407_0033 = list(
      description = "Phase I absolute-bioavailability study 1407-0033 in healthy men (50 mg oral with a 100 ug IV microtracer)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (another study)",
      notes = "Selects the study-specific oral bioavailability 0.487 and the additive logit shift -0.850 on the slow-path fraction (applied here as +0.850 on the logit of the main-path fraction). The paper did not use the slow-path shift for its simulations. The ka study effect of this study is fixed to 0 in the control stream (THETA(34) 0 FIX).",
      source_name = "PROJ == 14070033"
    ),
    FED = list(
      description = "Fed (1) vs fasted (0) administration",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted)",
      notes = "Control stream SWTFED = 1 when DIET == 1. Fed dosing shifts logit F by +0.586, multiplies ka by 1.25 and D1 by 3.52 x (dose / 50 mg)^-0.591, selects the fed dose slope on logit F in study 1407-0002, and selects the fed residual-error stratum (healthy subjects and patients).",
      source_name = "DIET"
    ),
    FORM_BEVUROGANT_ICF = list(
      description = "Intended commercial formulation (1) vs the trial formulation (0)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (trial formulation, FORM == 1)",
      notes = "Additive logit-F shift +0.287 and multiplicative ka shift 1 + 0.391 (control stream LF1FORM / KAFORM, FORM == 2). Coarse-milled and unmilled intended-commercial batches (FORM > 2) were excluded from the analysis and are not covered.",
      source_name = "FORM"
    ),
    REGI_BID = list(
      description = "Twice-daily dosing regimen indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (once-daily or single dose)",
      notes = "Empirical additive logit-F shift -0.448 (control stream LF1BID, BID == 1), introduced to describe the 200 mg twice-daily troughs.",
      source_name = "BID"
    ),
    DOSE_BEVUROGANT_MG = list(
      description = "Bevurogant dose of the most recent administration",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear effect on logit F, slope x (dose - 50 mg), with a slope per population (fasted healthy, fed healthy in study 1407-0002, patients), and a power effect (dose / 50 mg)^-0.591 on D1 under fed conditions only (control stream DEFFD1**SWTFED). The IV microtracer is not a dose of this column.",
      source_name = "DOSE"
    ),
    ROUTE_ORAL = list(
      description = "Analyte of the observation: 1 = total bevurogant (oral drug plus any IV microtracer), 0 = the IV microtracer alone",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (IV microtracer sample)",
      notes = "Carries the control stream ROUTE column on OBSERVATION records: 'IPRED = A(2)/V2' (IV states only) when ROUTE == 0, '(A(4) + A(2))/V2' otherwise; ROUTE == 0 also selects the IV residual stratum. Set 1 on every record except the IV-microtracer samples of study 1407-0033. Dose records do not read it (the dose compartment carries the route).",
      source_name = "ROUTE"
    ),
    SAMPLE_PREDOSE = list(
      description = "Trough (pre-dose) sample indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (post-dose sample)",
      notes = "Multiplies the whole residual SD (proportional and additive) by 1.27 (control stream 'IF(TROUGH==1) COMB = (PROP + ADD)*THETA(36)').",
      source_name = "TROUGH"
    ),
    LLOQ = list(
      description = "Lower limit of quantification of the bevurogant assay used for the observation",
      units = "nmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Scales every additive residual SD: additive SD = estimate x LLOQ (control stream 'ADD = EPS(2)*(LOQ/1) ; Scale additive error to LLOQ in assay'). The Appendix S1 dataset gives LOQ = 5 nmol/L for study 1407-0030; the LLOQs of the healthy-volunteer assays (oral and IV microtracer) are not reported.",
      source_name = "LOQ"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "bevurogant", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "bevurogant (oral dose)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "bevurogant (oral dose)", units = "mg", specimen = "tissue", verified = TRUE),
    central_iv = list(analyte = "bevurogant IV microtracer", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_iv = list(analyte = "bevurogant IV microtracer", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 353L,
    n_studies = 4L,
    n_observations = 7686L,
    age_range = "18-75 years",
    sex_female_pct = 26.6,
    disease_state = "Healthy subjects (n = 109) and patients with moderate-to-severe plaque psoriasis (n = 244)",
    dose_range = "25-200 mg single oral doses; 25-400 mg once daily or 200 mg twice daily for up to 28 days (healthy) or 12 weeks (patients); one 100 ug IV microtracer",
    regions = "Not reported",
    notes = "Phase I studies 1407-0002 (n = 61), 1407-0032 (n = 42) and 1407-0033 (n = 6) in healthy subjects and phase II study 1407-0030 in patients (Table 1). 73.4% men in the pooled PK population (Results). Weight, CRP and AST distributions are not reported."
  )

  ini({
    # Disposition (absolute: the IV microtracer identifies F). Reference 70 kg.
    lcl <- log(6.72); label("Clearance, healthy subject (L/h)") # Table 2 'CL (L/h)' = 6.72
    e_psoriasis_cl <- 0.773; label("Multiplier on CL for patients (unitless)") # Table 2 'Multiplier for patients on CL' = 0.773
    lvc <- log(19.5); label("Central volume of distribution (L)") # Table 2 'V c (L)' = 19.5
    lq <- log(23.1); label("Intercompartmental clearance (L/h)") # Table 2 'Q (L/h)' = 23.1
    lvp <- log(169); label("Peripheral volume of distribution (L)") # Table 2 'V p (L)' = 169
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent of WT / 70 kg on CL and Q (unitless)") # Table 2 'Allometric exponent for CL and Q' = 0.750
    e_wt_vc_vp <- fixed(1); label("Allometric exponent of WT / 70 kg on Vc and Vp (unitless)") # Table 2 'Allometric exponent for V c and V p' = 1.00
    e_crp_cl <- -0.120; label("Power exponent of CRP / 0.305 mg/L on CL (unitless)") # Table 2 'Exponent of CRP on CL' = -0.120
    e_ast_cl <- -0.232; label("Power exponent of AST / 13 U/L on CL (unitless)") # Table 2 'Exponent of AST on CL' = -0.232

    # Absolute oral bioavailability: one typical value per study (logit scale)
    logitfdepot_study0030 <- logit(0.553); label("Logit of oral bioavailability, study 1407-0030 (patients) (logit)") # Table 2 'F for study 1407-0030' = 0.553
    logitfdepot_study0002 <- logit(0.652); label("Logit of oral bioavailability, study 1407-0002 (logit)") # Table 2 'F for study 1407-0002' = 0.652
    logitfdepot_study0032 <- logit(0.429); label("Logit of oral bioavailability, study 1407-0032 (logit)") # Table 2 'F for study 1407-0032' = 0.429
    logitfdepot_study0033 <- logit(0.487); label("Logit of oral bioavailability, study 1407-0033 (logit)") # Table 2 'F for study 1407-0033' = 0.487
    e_fed_fdepot <- 0.586; label("Additive shift of fed state on logit F (logit)") # Table 2 'Additive shift of fed state on F' = 0.586
    e_icf_fdepot <- 0.287; label("Additive shift of the intended commercial formulation on logit F (logit)") # Table 2 'Additive shift of iCF on F' = 0.287
    e_bid_fdepot <- -0.448; label("Additive shift of twice-daily dosing on logit F (logit)") # Table 2 'Additive shift of twice daily dosing on F' = -0.448
    e_dose_fdepot_hvfasted <- -0.00618; label("Slope of (dose - 50 mg) on logit F, fasted healthy subjects (1/mg)") # Table 2 'Slope for dose effect on F for fasted healthy subjects' = -0.00618
    e_dose_fdepot_hvfed <- -0.00514; label("Slope of (dose - 50 mg) on logit F, fed healthy subjects (1/mg)") # Table 2 'Slope for dose effect on F for fed healthy subjects' = -0.00514
    e_dose_fdepot_pso <- -0.00517; label("Slope of (dose - 50 mg) on logit F, patients (1/mg)") # Table 2 'Slope for dose effect on F for fasted patients' = -0.00517

    # Main absorption path: zero-order into the depot after a lag, then first-order
    ltlag <- fixed(log(0.167)); label("Lag time of the main absorption path (h)") # Table 2 't lagmain (h)' = 0.167
    ld1 <- log(0.538); label("Zero-order input duration of the main absorption path, fasted (h)") # Table 2 'D 0main (h)' = 0.538
    e_fed_d1 <- 3.52; label("Multiplier for fed state on D1 (unitless)") # Table 2 'Multiplier for fed state on D 0main' = 3.52
    e_dose_d1 <- -0.591; label("Power exponent of dose / 50 mg on D1, fed state only (unitless)") # Table 2 'Power for dose effect on D 0main' = -0.591
    lka <- log(0.203); label("First-order absorption rate constant, fasted (1/h)") # Table 2 'k a (1/h)' = 0.203
    e_fed_ka <- 1.25; label("Multiplier for fed state on ka (unitless)") # Table 2 'Multiplier for fed state on ka' = 1.25
    e_icf_ka <- 0.391; label("Proportional shift of the intended commercial formulation on ka (fraction)") # Table 2 'Proportional shift of iCF on k a' = 0.391
    e_study0032_ka <- 0.575; label("Proportional shift of study 1407-0032 on ka (fraction)") # Table 2 'Proportional shift of study 1407-0032 on k a' = 0.575

    # Slow absorption path: zero-order into central after a long lag.
    # ffo = 1 - Fslow is the fraction of the bioavailable dose on the main path.
    logitffo <- -logit(0.218); label("Logit of the fraction of the bioavailable dose on the main path, 1 - Fslow (logit)") # Table 2 'F slow' = 0.218; logit(1 - x) = -logit(x)
    e_study0033_ffo <- 0.850; label("Additive shift of study 1407-0033 on logit(1 - Fslow) (logit)") # Table 2 'Additive shift of study 1407-0033 on F slow' = -0.850, sign reversed for the complement
    ltlag2 <- log(11.5); label("Lag time of the slow absorption path (h)") # Table 2 't lagslow (h)' = 11.5
    ld2 <- log(17.9); label("Zero-order input duration of the slow absorption path (h)") # Table 2 'D 0slow (h)' = 17.9

    # IIV (Table 2 SDs; control stream $OMEGA variances)
    etalcl ~ 0.173499 # Table 2 IIV CL 0.417 (SD); $OMEGA OM1_CL 0.173499
    etalvc ~ 0.545649 # Table 2 IIV V c 0.739 (SD); $OMEGA OM2_V2 0.545649
    etalka ~ 0.081342 # Table 2 IIV k a 0.285 (SD); $OMEGA OM5_KA 0.081342
    etalogitfdepot ~ 0.176443 # Table 2 IIV F 0.420 (SD, logit scale); $OMEGA OM6_LF1 0.176443
    etald1 ~ 0.576859 # Table 2 IIV D 0main 0.760 (SD); $OMEGA OM7_D1 0.576859
    etalogitffo ~ 0.599612 # Table 2 IIV F slow 0.774 (SD, logit scale); $OMEGA OM12_LABS4 0.599612; same variance on the complement
    etaRUV ~ 0.166522 # Table 2 'Proportional residual unexplained variability (CV)' 0.408 (SD); $OMEGA OM11_RUV 0.166522, on the proportional SD only

    # Residual error, one combined error per stratum. Additive SDs are in units
    # of the assay LLOQ (multiplied by LLOQ in model()).
    propSd_hvfasted <- 0.134; label("Proportional residual SD, oral, fasted, healthy subjects (fraction)") # Table 2 = 0.134; $SIGMA SG1_PROP 0.0178609
    addSd_hvfasted <- 2.03; label("Additive residual SD, oral, fasted, healthy subjects (multiples of LLOQ)") # Table 2 = 2.03; $SIGMA SG2_ADD 4.1273
    propSd_pso <- 0.203; label("Proportional residual SD, oral, fasted, patients (fraction)") # Table 2 = 0.203; $SIGMA SG5_PROP_PSORIASIS 0.0411327
    addSd_pso <- 5.34; label("Additive residual SD, oral, fasted, patients (multiples of LLOQ)") # Table 2 = 5.34; $SIGMA SG6_ADD_PSORIASIS 28.5187
    propSd_fed <- 0.173; label("Proportional residual SD, oral, fed, healthy subjects and patients (fraction)") # Table 2 = 0.173; $SIGMA SG7_PROP_FED 0.0300302
    addSd_fed <- 2.84; label("Additive residual SD, oral, fed, healthy subjects and patients (multiples of LLOQ)") # Table 2 = 2.84; $SIGMA SG8_ADD_FED 8.0833
    propSd_iv <- 0.377; label("Proportional residual SD, IV microtracer samples, healthy subjects (fraction)") # Table 2 = 0.377; $SIGMA SG3_PROP_IV 0.14224
    addSd_iv <- fixed(0.0316); label("Additive residual SD, IV microtracer samples, healthy subjects (multiples of LLOQ)") # Table 2 = 0.0316; $SIGMA SG4_ADD_IV 0.001 FIX
    e_predose_ruv <- 1.27; label("Multiplier on the residual SD for trough samples (unitless)") # Table 2 'Multiplier for Ctrough on residual unexplained variability' = 1.27
  })

  model({
    # mg -> nmol: Appendix S1 dataset AMT 93875.56 nmol for DOSE 50 mg
    # (molecular weight 532.62 g/mol)
    nmol_per_mg <- 93875.56 / 50

    # Absolute oral bioavailability (control stream TVLF1 / LF1 / FT)
    stdy_other <- STUDY_1407_0002 + STUDY_1407_0032 + STUDY_1407_0033
    logitf_study <- logitfdepot_study0030 * (1 - stdy_other) +
      logitfdepot_study0002 * STUDY_1407_0002 +
      logitfdepot_study0032 * STUDY_1407_0032 +
      logitfdepot_study0033 * STUDY_1407_0033
    # Dose slope: patients; else fed study 1407-0002; else fasted healthy
    slope_hv <- e_dose_fdepot_hvfasted * (1 - FED * STUDY_1407_0002) +
      e_dose_fdepot_hvfed * FED * STUDY_1407_0002
    slope_dose <- e_dose_fdepot_pso * DIS_PSORIASIS + slope_hv * (1 - DIS_PSORIASIS)
    logit_ft <- logitf_study +
      e_icf_fdepot * FORM_BEVUROGANT_ICF +
      e_bid_fdepot * REGI_BID +
      e_fed_fdepot * FED +
      slope_dose * (DOSE_BEVUROGANT_MG - 50) +
      etalogitfdepot
    ft <- expit(logit_ft)

    # Split between the main (depot) and slow (central) paths
    ffo <- expit(logitffo + e_study0033_ffo * STUDY_1407_0033 + etalogitffo)

    # Absorption parameters; the dose effect on D1 applies to fed doses only
    ka <- exp(lka + etalka) * e_fed_ka^FED *
      (1 + e_icf_ka * FORM_BEVUROGANT_ICF) *
      (1 + e_study0032_ka * STUDY_1407_0032)
    d1 <- exp(ld1 + etald1) * (e_fed_d1 * (DOSE_BEVUROGANT_MG / 50)^e_dose_d1)^FED
    tlag <- exp(ltlag)
    tlag2 <- exp(ltlag2)
    d2 <- exp(ld2)

    # Disposition
    cl <- exp(lcl + etalcl) * e_psoriasis_cl^DIS_PSORIASIS *
      (CRP / 0.305)^e_crp_cl * (AST / 13)^e_ast_cl * (WT / 70)^e_wt_cl_q
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    q <- exp(lq) * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp) * (WT / 70)^e_wt_vc_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Oral drug (control stream compartments PODEPO, POCENTR, POPER1)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    # IV microtracer, same disposition (IVCENTR, IVPER1)
    d/dt(central_iv) <- -kel * central_iv - k12 * central_iv + k21 * peripheral1_iv
    d/dt(peripheral1_iv) <- k12 * central_iv - k21 * peripheral1_iv

    # Each oral administration is two dose records with the same amt and
    # rate = -2: one into depot, one into central.
    f(depot) <- ft * ffo
    dur(depot) <- d1
    alag(depot) <- tlag
    f(central) <- ft * (1 - ffo)
    dur(central) <- d2
    alag(central) <- tlag2

    # Residual SD per stratum (control stream $ERROR, applied in this order):
    # oral fasted healthy by default, IV microtracer sample (ROUTE == 0),
    # then patient study (PROJ == 14070030), then fed (DIET == 1).
    iv_obs <- 1 - ROUTE_ORAL
    nonfed <- 1 - FED
    propSd_hv <- propSd_hvfasted * (1 - iv_obs) + propSd_iv * iv_obs
    addSd_hv <- addSd_hvfasted * (1 - iv_obs) + addSd_iv * iv_obs
    propSd_obs <- (propSd_hv * (1 - DIS_PSORIASIS) + propSd_pso * DIS_PSORIASIS) * nonfed +
      propSd_fed * FED
    addSd_obs <- (addSd_hv * (1 - DIS_PSORIASIS) + addSd_pso * DIS_PSORIASIS) * nonfed +
      addSd_fed * FED
    ruv_trough <- e_predose_ruv^SAMPLE_PREDOSE
    propSd_i <- propSd_obs * exp(etaRUV) * ruv_trough
    addSd_i <- addSd_obs * LLOQ * ruv_trough

    # Observation: the IV microtracer alone on IV-analyte samples
    # (ROUTE_ORAL = 0), otherwise the total of oral drug and microtracer.
    Cc_iv <- central_iv / vc * nmol_per_mg
    Cc <- (central * ROUTE_ORAL + central_iv) / vc * nmol_per_mg
    Cc ~ add(addSd_i) + prop(propSd_i) + combined2()
  })
}
