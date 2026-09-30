Bulitta_2024_contezolid <- function() {
  description <- paste(
    "Integrated population PK model for the intravenous double prodrug",
    "contezolid acefosamil (CZA), its intermediate MRX-1352, active",
    "contezolid and the inactive metabolite MRX-1320, jointly fitted to",
    "110 healthy volunteers (IV CZA 150-2400 mg; oral contezolid",
    "400-1200 mg fed/fasting) and 74 adult phase 2 patients with acute",
    "bacterial skin and skin structure infection (oral contezolid 800 mg",
    "q12h with food). CZA (amount only, not assayed) converts to",
    "MRX-1352 by Michaelis-Menten kinetics with a first-order loss;",
    "MRX-1352, contezolid and MRX-1320 each have two-compartment",
    "disposition. The MRX-1352-to-contezolid conversion clearance is",
    "auto-induced by MRX-1352 concentration (sigmoid Emax) and rises",
    "sigmoidally with time since the first dose; the MRX-1352 loss",
    "clearance shares the same auto-induction EC50 and Hill coefficient.",
    "Contezolid is converted entirely by first-order clearance to",
    "MRX-1320, which is eliminated by parallel linear and",
    "Michaelis-Menten routes. Oral contezolid passes three transit",
    "compartments (each with mean time Tlag/3) and a gut compartment",
    "with first-order absorption; lag time, absorption half-life and",
    "relative bioavailability depend on the fed state and the",
    "400/800/1200 mg dose level. Allometric weight scaling (70 kg;",
    "exponents 0.75 clearances, 1 volumes). Separate typical contezolid",
    "clearance and between-subject variability for healthy volunteers",
    "and patients. Fitted in S-ADAPT (importance sampling)."
  )
  reference <- paste(
    "Bulitta JB, Fang E, Stryjewski ME, Wang W, Atiee GJ, Stark JG,",
    "Hafkin B. Population pharmacokinetic rationale for intravenous",
    "contezolid acefosamil followed by oral contezolid dosage regimens.",
    "Antimicrob Agents Chemother. 2024;68(4):e01400-23.",
    "doi:10.1128/aac.01400-23. Preliminary results presented as",
    "Bulitta JB, Hafkin B, Fang E. 1118. Population Pharmacokinetics of",
    "Contezolid Acefosamil and Contezolid - Rationale for a Safe and",
    "Effective Loading Dose Regimen. Open Forum Infect Dis.",
    "2021;8(Suppl 1):S651. doi:10.1093/ofid/ofab466.1311."
  )
  vignette <- "Bulitta_2024_contezolid"

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling of every clearance (including the two",
        "Michaelis-Menten maximum rates) with exponent 0.75 and every",
        "volume with exponent 1.0, reference weight 70 kg (Methods,",
        "'Covariate effects'). Both exponents are fixed."
      ),
      source_name = "WT"
    ),
    DIS_CSSSI = list(
      description = paste(
        "1 = phase 2 patient with acute bacterial skin and skin structure",
        "infection (ABSSSI; study MRX-I-03), 0 = healthy volunteer",
        "(studies MRX4-002 and MRX-I-02)"
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy volunteer)",
      notes = paste(
        "Selects the contezolid total-clearance typical value AND its",
        "between-subject variability: Table 2 reports CL_Con,HV = 10.2",
        "L/h (BSV 0.234) and CL_Con,PA = 11.3 L/h (BSV 0.538) as two",
        "separately estimated parameters. ABSSSI is the same entity the",
        "register records as cSSSI."
      ),
      source_name = "healthy volunteer vs phase 2 patient"
    ),
    FED = list(
      description = "1 = oral contezolid dose taken with food, 0 = fasting",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (fed; the reference state of Table 2)",
      notes = paste(
        "Time-varying per dose record. Selects the fed or fasting lag",
        "time, absorption half-life and relative bioavailability. In",
        "MRX-I-02 the fed state was a high-fat, low-tyramine meal; all",
        "multiple-dose and phase 2 dosing was with food. Only affects",
        "oral contezolid doses."
      ),
      source_name = "fed / fasting"
    ),
    DOSE_CONTEZOLID_MG = list(
      description = "Administered oral contezolid dose (mg)",
      units = "mg",
      type = "continuous",
      reference_category = "800 mg (the reference dose of Table 2)",
      notes = paste(
        "Time-varying per dose record. Used only to pick the relative",
        "bioavailability row of Table 2, which the paper estimated for",
        "the three studied dose levels: values below 600 mg use the",
        "400 mg row, 600-1000 mg use the 800 mg row and above 1000 mg",
        "use the 1200 mg row. The paper gives no interpolation rule;",
        "the band edges are the midpoints between studied doses."
      ),
      source_name = "dose level"
    ),
    STUDY_MRX4002 = list(
      description = "1 = subject from the IV contezolid acefosamil study MRX4-002, 0 = oral contezolid studies MRX-I-02 and MRX-I-03",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (oral contezolid studies)",
      notes = paste(
        "Selects the residual-error terms only: the paper estimated",
        "separate additive and proportional residual errors for",
        "contezolid and MRX-1320 after IV CZA dosing (MRX4-002) and after",
        "oral contezolid dosing (Table 2 footnote h). MRX-1352 was only",
        "measured in MRX4-002 and has one residual error."
      ),
      source_name = "study"
    ),
    OCC = list(
      description = "Occasion index for the between-occasion variability of oral bioavailability",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "The paper includes between-occasion variability (BOV) in oral",
        "bioavailability (Methods, 'Oral contezolid'; Table 2 footnote d)",
        "but does not define the occasions. Encoded with two occasions,",
        "matching the two crossover periods (fed / fasting, days 1 and 7)",
        "of the single-dose part of MRX-I-02: OCC = 1 or 2 draws the",
        "corresponding occasion eta; any other value (e.g. 0) switches",
        "BOV off. For a single-occasion simulation set OCC = 1."
      ),
      source_name = "occasion"
    )
  )

  compartmentData <- list(
    depot_iv = list(analyte = "contezolid acefosamil", units = "mg", specimen = "not applicable", verified = TRUE),
    central_mrx1352 = list(analyte = "MRX-1352", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_mrx1352 = list(analyte = "MRX-1352", units = "mg", specimen = "plasma", verified = TRUE),
    transit1 = list(analyte = "contezolid", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "contezolid", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "contezolid", units = "mg", specimen = "administration site", verified = TRUE),
    depot = list(analyte = "contezolid", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "contezolid", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "contezolid", units = "mg", specimen = "plasma", verified = TRUE),
    central_mrx1320 = list(analyte = "MRX-1320", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_mrx1320 = list(analyte = "MRX-1320", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 184L,
    n_studies = 3L,
    age_range = paste(
      "adults; mean +/- SD 35.9 +/- 10.2 years (MRX4-002), 24.6 +/- 3.3",
      "years (MRX-I-02), 38.4 +/- 10.7 years (MRX-I-03)"
    ),
    weight_range = paste(
      "mean +/- SD 77.9 +/- 10.0 kg (MRX4-002), 71.2 +/- 12.7 kg",
      "(MRX-I-02), 82.8 +/- 21.7 kg (MRX-I-03)"
    ),
    sex_female_pct = 37.0,
    race_ethnicity = c(Caucasian = 65.8, Black = 19.6, Hispanic = 10.9, Asian = 1.6, Other = 2.2),
    disease_state = paste(
      "110 healthy volunteers (MRX4-002, MRX-I-02) and 74 adult phase 2",
      "patients with acute bacterial skin and skin structure infection",
      "(MRX-I-03)"
    ),
    dose_range = paste(
      "MRX4-002: IV CZA 150-2400 mg single doses (60- or 90-min",
      "infusion) and 600 or 900 mg q12h or 2100 mg q24h for 10 days.",
      "MRX-I-02: oral contezolid 400, 800 or 1200 mg single doses with",
      "and without food, and 800 mg q12h with food for 14 or 28 days.",
      "MRX-I-03: oral contezolid 800 mg q12h with food for 10 days."
    ),
    regions = "USA and Australia",
    notes = paste(
      "Demographics from Results paragraph 1 and Table 1. Sex counts:",
      "MRX4-002 43M/23F, MRX-I-02 22M/22F, MRX-I-03 51M/23F (68/184",
      "female). Race: 121 Caucasian, 36 Black or African American, 20",
      "Hispanic, 3 Asian, 4 other. MRX-1352 was measured only in",
      "MRX4-002; MRX-1320 in MRX4-002 and the MRX-I-02 multiple-dose",
      "cohorts; contezolid in all three studies. BLQ samples handled",
      "by the Beal M3 method."
    )
  )

  ini({
    # All values are Bulitta 2024 Table 2 ('Population mean' column) at
    # the 70 kg reference weight. Between-subject variabilities are the
    # 'apparent coefficient of variation of a normal distribution on
    # natural logarithmic scale' (Table 2 footnote c), i.e. the SD of the
    # log-scale random effect, so omega^2 = BSV^2 (same S-ADAPT convention
    # as Bulitta_2019_pefloxacin).
    #
    # Half-lives and lag times printed in minutes are converted to the
    # model time unit (h) as the rate constants k = log(2) / (T1/2 / 60)
    # and ktr = 3 / (Tlag / 60).

    # ---- Contezolid disposition -----------------------------------------
    lcl <- log(10.2)
    label("Contezolid total clearance, healthy volunteers (L/h)") # Table 2: CL_Con,HV = 10.2 L/h (RSE 5.2%)
    lcl_csssi <- log(11.3)
    label("Contezolid total clearance, ABSSSI patients (L/h)") # Table 2: CL_Con,PA = 11.3 L/h (RSE 8.9%)
    lq <- log(69.9)
    label("Contezolid distribution clearance CLd (L/h)") # Table 2: CLd_Con = 69.9 L/h (RSE 18.1%)
    lvc <- fixed(log(3.0))
    label("Contezolid central volume V1 (L)") # Table 2: V1_Con = 3.0 L (fixed); Results: 'fixed to 3.0 L ... since initial estimates were smaller'
    lvp <- log(14.1)
    label("Contezolid peripheral volume V2 (L)") # Table 2: V2_Con = 14.1 L (RSE 7.0%)

    # ---- Oral contezolid absorption ------------------------------------
    lfdepot <- log(0.640)
    label("Bioavailability of 800 mg oral contezolid with food, relative to IV CZA converted to contezolid (fraction)") # Table 2: F = 0.640 (RSE 7.2%); footnote e
    lfdepot_400fed <- log(0.924)
    label("Relative bioavailability of 400 mg fed vs 800 mg fed (ratio)") # Table 2: F_rel,400mg,fed = 0.924 (RSE 9.5%)
    lfdepot_400fasted <- log(0.606)
    label("Relative bioavailability of 400 mg fasting vs 800 mg fed (ratio)") # Table 2: F_rel,400mg,fasting = 0.606 (RSE 15.1%)
    lfdepot_800fasted <- log(0.467)
    label("Relative bioavailability of 800 mg fasting vs 800 mg fed (ratio)") # Table 2: F_rel,800mg,fasting = 0.467 (RSE 8.5%)
    lfdepot_1200fed <- log(0.820)
    label("Relative bioavailability of 1200 mg fed vs 800 mg fed (ratio)") # Table 2: F_rel,1200mg,fed = 0.820 (RSE 6.0%)
    lfdepot_1200fasted <- log(0.485)
    label("Relative bioavailability of 1200 mg fasting vs 800 mg fed (ratio)") # Table 2: F_rel,1200mg,fasting = 0.485 (RSE 8.7%)
    lmtt <- log(47.5 / 60)
    label("Oral absorption lag time Tlag with food = mean transit time of the three transit compartments (h)") # Table 2: T_lag,fed = 47.5 min (RSE 11.1%); Fig. 2 each transit step Tlag/3
    lmtt_fasted <- log(22.5 / 60)
    label("Oral absorption lag time Tlag fasting (h)") # Table 2: T_lag,fasting = 22.5 min (RSE 18.8%)
    lka <- log(log(2) / (59.3 / 60))
    label("First-order absorption rate constant from the gut with food (1/h)") # Table 2: T_1/2,abs,fed = 59.3 min (RSE 4.4%) -> ka = log(2)/(59.3/60) = 0.7013 /h
    lka_fasted <- log(log(2) / (69.4 / 60))
    label("First-order absorption rate constant from the gut, fasting (1/h)") # Table 2: T_1/2,abs,fasting = 69.4 min (RSE 7.8%) -> ka = log(2)/(69.4/60) = 0.5993 /h

    # ---- Contezolid acefosamil (CZA) -----------------------------------
    lvmax_cza <- log(802)
    label("Maximum rate of conversion of CZA to MRX-1352 (mg CZA/h)") # Table 2: Vmax_CZA = 802 mg/h (RSE 5.3%)
    lkm_cza <- log(0.960)
    label("Amount of CZA giving half of Vmax_CZA, AM50 (mg)") # Table 2: AM_50,CZA = 0.960 mg (RSE 6.6%)
    lkel_cza <- log(log(2) / (19.7 / 60))
    label("First-order loss rate constant of CZA (1/h)") # Table 2: T_1/2,loss = 19.7 min (RSE 28.0%) -> k = log(2)/(19.7/60) = 2.111 /h

    # ---- MRX-1352 ------------------------------------------------------
    lcl_mrx1352 <- log(1.46)
    label("MRX-1352-to-contezolid conversion clearance without induction at time 0, CL_1352,0 (L/h)") # Table 2: CL_1352,0 = 1.46 L/h (RSE 5.3%)
    lcl_ss_mrx1352 <- log(6.60)
    label("MRX-1352-to-contezolid conversion clearance without induction at steady state, CL_1352,SS (L/h)") # Table 2: CL_1352,SS = 6.60 L/h (RSE 6.3%)
    lcl_t50_mrx1352 <- log(116)
    label("Time past first dose of half-maximal increase in MRX-1352 conversion clearance, TC50 (h)") # Table 2: TC_50 = 116 h (RSE 4.5%)
    lcl_time_hill_mrx1352 <- log(2.72)
    label("Hill coefficient of the time-dependent MRX-1352 clearance increase, Ht (unitless)") # Table 2: Ht = 2.72 (RSE 6.8%)
    lemax_mrx1352 <- log(145)
    label("Maximum extent of auto-induction of the MRX-1352 conversion clearance, Emax (unitless)") # Table 2: Emax = 145 (RSE 19.2%)
    lec50_mrx1352 <- log(64.9)
    label("MRX-1352 concentration giving half of Emax (and of Emax_Loss), EC50 (mg/L)") # Table 2: EC_50 = 64.9 mg/L (RSE 6.0%)
    lhill_mrx1352 <- log(7.01)
    label("Hill coefficient of the MRX-1352 auto-induction, H (unitless)") # Table 2: H = 7.01 (RSE 6.6%)
    lcl_loss_mrx1352 <- log(0.00987)
    label("Loss clearance of MRX-1352 without auto-induction, CL_Loss,1352 (L/h)") # Table 2: CL_Loss,1352 = 0.00987 L/h (RSE 20.8%)
    lemax_loss_mrx1352 <- log(29400)
    label("Maximum extent of auto-induction of the MRX-1352 loss clearance, Emax_Loss (unitless)") # Table 2: Emax_Loss = 29,400 (RSE 5.7%)
    lq_mrx1352 <- log(10.2)
    label("MRX-1352 distribution clearance CLd_1352 (L/h)") # Table 2: CLd_1352 = 10.2 L/h (RSE 8.6%)
    lvc_mrx1352 <- log(7.61)
    label("MRX-1352 central volume V1_1352 (L)") # Table 2: V1_1352 = 7.61 L (RSE 5.3%)
    lvp_mrx1352 <- log(13.0)
    label("MRX-1352 peripheral volume V2_1352 (L)") # Table 2: V2_1352 = 13.0 L (RSE 5.9%)

    # ---- MRX-1320 ------------------------------------------------------
    lcl_mrx1320 <- log(0.853)
    label("MRX-1320 linear elimination clearance CL_1320 (L/h)") # Table 2: CL_1320 = 0.853 L/h (RSE 18.4%)
    lvmax_mrx1320 <- log(74.2)
    label("MRX-1320 maximum rate of saturable elimination Vmax_1320 (mg/h)") # Table 2: Vmax_1320 = 74.2 mg/h (RSE 7.4%)
    lkm_mrx1320 <- log(0.657)
    label("MRX-1320 Michaelis-Menten constant Km_1320 (mg/L)") # Table 2: Km_1320 = 0.657 mg/L (RSE 9.0%)
    lq_mrx1320 <- log(16.4)
    label("MRX-1320 distribution clearance CLd_1320 (L/h)") # Table 2: CLd_1320 = 16.4 L/h (RSE 7.1%)
    lvc_mrx1320 <- log(29.6)
    label("MRX-1320 central volume V1_1320 (L)") # Table 2: V1_1320 = 29.6 L (RSE 8.6%)
    lvp_mrx1320 <- log(67.9)
    label("MRX-1320 peripheral volume V2_1320 (L)") # Table 2: V2_1320 = 67.9 L (RSE 9.4%)

    # ---- Allometric exponents ------------------------------------------
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent on all clearances and Michaelis-Menten maximum rates (unitless)") # Methods 'Covariate effects': 'fixed to 0.75 for all clearances'
    e_wt_vc <- fixed(1)
    label("Allometric exponent on all volumes of distribution (unitless)") # Methods 'Covariate effects': 'and to 1.0 for volumes'

    # ---- Between-subject variability (omega^2 = BSV^2) -----------------
    etalcl ~ 0.054756 # Table 2: CL_Con,HV BSV 0.234 -> 0.234^2
    etalcl_csssi ~ 0.289444 # Table 2: CL_Con,PA BSV 0.538 -> 0.538^2
    etalq ~ 1.8769 # Table 2: CLd_Con BSV 1.37 -> 1.37^2
    etalvc ~ 0.079524 # Table 2: V1_Con BSV 0.282 -> 0.282^2
    etalvp ~ 0.077284 # Table 2: V2_Con BSV 0.278 -> 0.278^2
    etalfdepot ~ fixed(0.0225) # Table 2: F BSV 0.150 held constant -> 0.150^2
    etalfdepot_400fed ~ fixed(0.01) # Table 2: F_rel,400mg,fed BSV 0.100 held constant -> 0.1^2
    etalfdepot_400fasted ~ fixed(0.01) # Table 2: F_rel,400mg,fasting BSV 0.100 held constant -> 0.1^2
    etalfdepot_800fasted ~ fixed(0.01) # Table 2: F_rel,800mg,fasting BSV 0.100 held constant -> 0.1^2
    etalfdepot_1200fed ~ fixed(0.01) # Table 2: F_rel,1200mg,fed BSV 0.100 held constant -> 0.1^2
    etalfdepot_1200fasted ~ fixed(0.01) # Table 2: F_rel,1200mg,fasting BSV 0.100 held constant -> 0.1^2
    etalmtt ~ 0.398161 # Table 2: T_lag,fed BSV 0.631 -> 0.631^2
    etalmtt_fasted ~ 0.781456 # Table 2: T_lag,fasting BSV 0.884 -> 0.884^2
    etalka ~ fixed(0.0225) # Table 2: T_1/2,abs,fed BSV 0.150 held constant -> 0.150^2
    etalka_fasted ~ 0.145161 # Table 2: T_1/2,abs,fasting BSV 0.381 -> 0.381^2
    etalvmax_cza ~ 0.131044 # Table 2: Vmax_CZA BSV 0.362 -> 0.362^2
    etalkm_cza ~ 0.0225 # Table 2: AM_50,CZA BSV 0.150 (RSE 112%) -> 0.150^2
    etalkel_cza ~ 0.294849 # Table 2: T_1/2,loss BSV 0.543 -> 0.543^2
    etalcl_mrx1352 ~ 0.04 # Table 2: CL_1352,0 BSV 0.200 -> 0.200^2
    etalcl_ss_mrx1352 ~ 0.011881 # Table 2: CL_1352,SS BSV 0.109 -> 0.109^2
    etalcl_t50_mrx1352 ~ fixed(0.04) # Table 2: TC_50 BSV 0.200 held constant -> 0.200^2
    etalcl_time_hill_mrx1352 ~ fixed(0.01) # Table 2: Ht BSV 0.100 held constant -> 0.100^2
    etalemax_mrx1352 ~ 0.329476 # Table 2: Emax BSV 0.574 -> 0.574^2
    etalec50_mrx1352 ~ 0.016129 # Table 2: EC_50 BSV 0.127 -> 0.127^2
    etalhill_mrx1352 ~ fixed(0.01) # Table 2: H BSV 0.100 held constant -> 0.100^2
    etalcl_loss_mrx1352 ~ 0.744769 # Table 2: CL_Loss,1352 BSV 0.863 -> 0.863^2
    etalemax_loss_mrx1352 ~ 0.00456976 # Table 2: Emax_Loss BSV 0.0676 -> 0.0676^2
    etalq_mrx1352 ~ 0.077284 # Table 2: CLd_1352 BSV 0.278 -> 0.278^2
    etalvc_mrx1352 ~ 0.021316 # Table 2: V1_1352 BSV 0.146 -> 0.146^2
    etalvp_mrx1352 ~ 0.045369 # Table 2: V2_1352 BSV 0.213 -> 0.213^2
    etalcl_mrx1320 ~ 0.026896 # Table 2: CL_1320 BSV 0.164 -> 0.164^2
    # Table 2 footnote f: correlation 0.949 between Vmax_1320 and Km_1320;
    # covariance = 0.949 * 0.306 * 0.491 = 0.14258.
    etalvmax_mrx1320 + etalkm_mrx1320 ~ c(0.093636, 0.14258, 0.241081) # Table 2: Vmax_1320 BSV 0.306, Km_1320 BSV 0.491, footnote f r = 0.949
    etalq_mrx1320 ~ 0.071824 # Table 2: CLd_1320 BSV 0.268 -> 0.268^2
    etalvc_mrx1320 ~ 0.320356 # Table 2: V1_1320 BSV 0.566 -> 0.566^2
    etalvp_mrx1320 ~ 0.341056 # Table 2: V2_1320 BSV 0.584 -> 0.584^2

    # ---- Between-occasion variability of F (two occasions, see OCC) ----
    etaiov_lfdepot_oc1 ~ 0.018225 # Table 2: F BOV 0.135 (RSE 42.5%; footnote d) -> 0.135^2
    etaiov_lfdepot_oc2 ~ fixed(0.018225) # same variance as occasion 1 (one BOV estimate in Table 2)

    # ---- Residual error (Table 2 footnote h) ---------------------------
    addSd <- 0.00233
    label("Contezolid additive residual SD after oral contezolid (mg/L)") # Table 2 footnote h: 0.00233 mg/L (RSE 11.5%)
    propSd <- 0.536
    label("Contezolid proportional residual SD after oral contezolid (fraction)") # Table 2 footnote h: CV 0.536 (RSE 3.1%)
    addSd_iv <- 0.0392
    label("Contezolid additive residual SD after IV CZA (mg/L)") # Table 2 footnote h: 0.0392 mg/L (RSE 26.8%)
    propSd_iv <- 0.192
    label("Contezolid proportional residual SD after IV CZA (fraction)") # Table 2 footnote h: CV 0.192 (RSE 3.5%)
    addSd_mrx1352 <- 0.108
    label("MRX-1352 additive residual SD (mg/L)") # Table 2 footnote h: 0.108 mg/L (RSE 21.2%)
    propSd_mrx1352 <- 0.202
    label("MRX-1352 proportional residual SD (fraction)") # Table 2 footnote h: CV 0.202 (RSE 2.7%)
    addSd_mrx1320 <- 0.00133
    label("MRX-1320 additive residual SD after oral contezolid (mg/L)") # Table 2 footnote h: 0.00133 mg/L (RSE 24.0%)
    propSd_mrx1320 <- 0.495
    label("MRX-1320 proportional residual SD after oral contezolid (fraction)") # Table 2 footnote h: CV 0.495 (RSE 4.4%)
    addSd_mrx1320_iv <- 0.00636
    label("MRX-1320 additive residual SD after IV CZA (mg/L)") # Table 2 footnote h: 0.00636 mg/L (RSE 17.7%)
    propSd_mrx1320_iv <- 0.117
    label("MRX-1320 proportional residual SD after IV CZA (fraction)") # Table 2 footnote h: CV 0.117 (RSE 4.2%)
  })

  model({
    # Molecular masses (Methods, 'Population pharmacokinetics'): CZA
    # 552.33, MRX-1352 487.30, contezolid 408.33, MRX-1320 444.36 Da. All
    # states are in mg of their own species, so every conversion flux is
    # multiplied by the product-to-substrate mass ratio.
    mw_cza <- 552.33
    mw_mrx1352 <- 487.30
    mw_con <- 408.33
    mw_mrx1320 <- 444.36

    # Allometric size factors (70 kg reference).
    wt_cl <- (WT / 70)^e_wt_cl
    wt_v <- (WT / 70)^e_wt_vc

    # Oral-dose condition indicators.
    fasted <- 1 - FED
    d400 <- 0
    d1200 <- 0
    if (DOSE_CONTEZOLID_MG < 600) d400 <- 1
    if (DOSE_CONTEZOLID_MG > 1000) d1200 <- 1
    d800 <- 1 - d400 - d1200

    # Occasion indicators for the BOV on F.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)

    # ---- Contezolid ----------------------------------------------------
    # Healthy volunteers and patients each carry their own typical value
    # and BSV (Table 2); DIS_CSSSI selects one of the two.
    cl <- exp((lcl + etalcl) * (1 - DIS_CSSSI) + (lcl_csssi + etalcl_csssi) * DIS_CSSSI) * wt_cl
    q <- exp(lq + etalq) * wt_cl
    vc <- exp(lvc + etalvc) * wt_v
    vp <- exp(lvp + etalvp) * wt_v

    # Oral bioavailability: F(800 mg fed) times the relative bioavailability
    # of the dose-level x fed-state cell (1 for 800 mg fed), with BSV on F,
    # BSV on each F_rel and BOV on F.
    lfrel <- d400 * FED * (lfdepot_400fed + etalfdepot_400fed) +
      d400 * fasted * (lfdepot_400fasted + etalfdepot_400fasted) +
      d800 * fasted * (lfdepot_800fasted + etalfdepot_800fasted) +
      d1200 * FED * (lfdepot_1200fed + etalfdepot_1200fed) +
      d1200 * fasted * (lfdepot_1200fasted + etalfdepot_1200fasted)
    fdepot <- exp(lfdepot + etalfdepot + lfrel + oc1 * etaiov_lfdepot_oc1 + oc2 * etaiov_lfdepot_oc2)

    mtt <- exp((lmtt + etalmtt) * FED + (lmtt_fasted + etalmtt_fasted) * fasted)
    ka <- exp((lka + etalka) * FED + (lka_fasted + etalka_fasted) * fasted)
    ktr <- 3 / mtt

    # ---- CZA -----------------------------------------------------------
    vmax_cza <- exp(lvmax_cza + etalvmax_cza) * wt_cl
    km_cza <- exp(lkm_cza + etalkm_cza)
    kel_cza <- exp(lkel_cza + etalkel_cza)

    # ---- MRX-1352 ------------------------------------------------------
    cl0_mrx1352 <- exp(lcl_mrx1352 + etalcl_mrx1352) * wt_cl
    cl_ss_mrx1352 <- exp(lcl_ss_mrx1352 + etalcl_ss_mrx1352) * wt_cl
    cl_t50_mrx1352 <- exp(lcl_t50_mrx1352 + etalcl_t50_mrx1352)
    cl_time_hill_mrx1352 <- exp(lcl_time_hill_mrx1352 + etalcl_time_hill_mrx1352)
    emax_mrx1352 <- exp(lemax_mrx1352 + etalemax_mrx1352)
    ec50_mrx1352 <- exp(lec50_mrx1352 + etalec50_mrx1352)
    hill_mrx1352 <- exp(lhill_mrx1352 + etalhill_mrx1352)
    cl_loss_mrx1352 <- exp(lcl_loss_mrx1352 + etalcl_loss_mrx1352) * wt_cl
    emax_loss_mrx1352 <- exp(lemax_loss_mrx1352 + etalemax_loss_mrx1352)
    q_mrx1352 <- exp(lq_mrx1352 + etalq_mrx1352) * wt_cl
    vc_mrx1352 <- exp(lvc_mrx1352 + etalvc_mrx1352) * wt_v
    vp_mrx1352 <- exp(lvp_mrx1352 + etalvp_mrx1352) * wt_v

    # ---- MRX-1320 ------------------------------------------------------
    cl_mrx1320 <- exp(lcl_mrx1320 + etalcl_mrx1320) * wt_cl
    vmax_mrx1320 <- exp(lvmax_mrx1320 + etalvmax_mrx1320) * wt_cl
    km_mrx1320 <- exp(lkm_mrx1320 + etalkm_mrx1320)
    q_mrx1320 <- exp(lq_mrx1320 + etalq_mrx1320) * wt_cl
    vc_mrx1320 <- exp(lvc_mrx1320 + etalvc_mrx1320) * wt_v
    vp_mrx1320 <- exp(lvp_mrx1320 + etalvp_mrx1320) * wt_v

    # ---- Concentrations and time-varying clearances --------------------
    Cc <- central / vc
    Cc_mrx1352 <- central_mrx1352 / vc_mrx1352
    Cc_mrx1320 <- central_mrx1320 / vc_mrx1320

    # Eq 3-4: time-dependent base conversion clearance of MRX-1352. T is
    # time past the first dose, so the data set must start at the first
    # dose (t = 0).
    fct_cl_mrx1352 <- t^cl_time_hill_mrx1352 / (t^cl_time_hill_mrx1352 + cl_t50_mrx1352^cl_time_hill_mrx1352)
    clt_mrx1352 <- cl0_mrx1352 + fct_cl_mrx1352 * (cl_ss_mrx1352 - cl0_mrx1352)

    # Eq 1, 5, 6, 7: immediate auto-induction by the MRX-1352 concentration,
    # sharing EC50 and H between the conversion and loss clearances. The
    # concentration is floored at zero for the power only: once MRX-1352 has
    # washed out the solver can return tiny negative amounts, and a negative
    # base raised to the non-integer H is NaN, which would poison every state.
    c_ind_mrx1352 <- Cc_mrx1352
    if (c_ind_mrx1352 < 0) c_ind_mrx1352 <- 0
    hill_term <- c_ind_mrx1352^hill_mrx1352 / (c_ind_mrx1352^hill_mrx1352 + ec50_mrx1352^hill_mrx1352)
    cl_conv_mrx1352 <- (1 + emax_mrx1352 * hill_term) * clt_mrx1352
    cl_lossind_mrx1352 <- (1 + emax_loss_mrx1352 * hill_term) * cl_loss_mrx1352

    # ---- Fluxes (mg of the substrate species per h) --------------------
    rate_cza <- vmax_cza * depot_iv / (km_cza + depot_iv)
    rate_mrx1352 <- cl_conv_mrx1352 * Cc_mrx1352
    rate_con <- cl * Cc

    # ---- ODEs ----------------------------------------------------------
    d/dt(depot_iv) <- -rate_cza - kel_cza * depot_iv
    d/dt(central_mrx1352) <- rate_cza * mw_mrx1352 / mw_cza - rate_mrx1352 -
      cl_lossind_mrx1352 * Cc_mrx1352 -
      q_mrx1352 * Cc_mrx1352 + q_mrx1352 * peripheral1_mrx1352 / vp_mrx1352
    d/dt(peripheral1_mrx1352) <- q_mrx1352 * Cc_mrx1352 - q_mrx1352 * peripheral1_mrx1352 / vp_mrx1352

    d/dt(transit1) <- -ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(depot) <- ktr * transit3 - ka * depot
    d/dt(central) <- ka * depot + rate_mrx1352 * mw_con / mw_mrx1352 - rate_con -
      q * Cc + q * peripheral1 / vp
    d/dt(peripheral1) <- q * Cc - q * peripheral1 / vp

    d/dt(central_mrx1320) <- rate_con * mw_mrx1320 / mw_con - cl_mrx1320 * Cc_mrx1320 -
      vmax_mrx1320 * Cc_mrx1320 / (km_mrx1320 + Cc_mrx1320) -
      q_mrx1320 * Cc_mrx1320 + q_mrx1320 * peripheral1_mrx1320 / vp_mrx1320
    d/dt(peripheral1_mrx1320) <- q_mrx1320 * Cc_mrx1320 - q_mrx1320 * peripheral1_mrx1320 / vp_mrx1320

    # Oral contezolid enters the first transit compartment (Fig. 2, Stom1).
    f(transit1) <- fdepot

    # ---- Residual error ------------------------------------------------
    # Separate contezolid and MRX-1320 residual errors after IV CZA
    # (study MRX4-002) and after oral contezolid (Table 2 footnote h).
    addSd_cc <- addSd * (1 - STUDY_MRX4002) + addSd_iv * STUDY_MRX4002
    propSd_cc <- propSd * (1 - STUDY_MRX4002) + propSd_iv * STUDY_MRX4002
    addSd_cc_mrx1320 <- addSd_mrx1320 * (1 - STUDY_MRX4002) + addSd_mrx1320_iv * STUDY_MRX4002
    propSd_cc_mrx1320 <- propSd_mrx1320 * (1 - STUDY_MRX4002) + propSd_mrx1320_iv * STUDY_MRX4002

    Cc ~ add(addSd_cc) + prop(propSd_cc)
    Cc_mrx1352 ~ add(addSd_mrx1352) + prop(propSd_mrx1352)
    Cc_mrx1320 ~ add(addSd_cc_mrx1320) + prop(propSd_cc_mrx1320)
  })
}
