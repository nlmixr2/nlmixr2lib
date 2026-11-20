Colin_2025_remimazolam_umss <- function() {
  description <- paste(
    "Sequential population PK/PD proportional-odds model for the University",
    "of Michigan Sedation Scale (UMSS, 0 awake to 4 unarousable) produced by",
    "intravenous remimazolam, with or without fentanyl, in 31 children and",
    "adolescents aged 6 to 18 years undergoing procedural sedation (Colin",
    "2025, NCT04851717). The remimazolam / CNS 7054 pharmacokinetics of",
    "modellib('Colin_2025_remimazolam') are carried forward unchanged and",
    "fentanyl plasma concentrations come from the weight-proportional",
    "two-compartment paediatric model of Ginsberg 1996. Each drug drives its",
    "own effect compartment; the effect-site concentrations, normalised by",
    "category-specific EC50s, are summed under a simplified Minto response",
    "surface (additive interaction) and give the probabilities of UMSS >= 1,",
    ">= 2 and >= 3 through one sigmoid with a shared Hill coefficient. UMSS 4",
    "was not observed and is not modelled. Age, sex, fentanyl",
    "coadministration and procedure type were tested on the remimazolam EC50",
    "and none was retained.",
    sep = " "
  )
  reference <- paste(
    "Colin PJ, Bichajian LH, Curt VR, Koomen JV, Stohr T, Struys MMRF,",
    "Mason KP. Pharmacokinetics and Pharmacodynamics of Remimazolam for",
    "Procedural Sedation in Children and Adolescents. Anesthesiology.",
    "2025;143(2):368-382. doi:10.1097/ALN.0000000000005560. PMC12227205.",
    "Fentanyl pharmacokinetics: Ginsberg B, Howell S, Glass PS, et al.",
    "Pharmacokinetic model-driven infusion of fentanyl in children.",
    "Anesthesiology. 1996;85:1268-1275. doi:10.1097/00000542-199612000-00007.",
    "Interaction model: Minto CF, et al. Response surface model for",
    "anesthetic drug interactions. Anesthesiology. 2000;92:1603-1616.",
    "doi:10.1097/00000542-200006000-00017.",
    "Pharmacokinetic layer: see modellib('Colin_2025_remimazolam').",
    sep = " "
  )
  vignette <- "Colin_2025_remimazolam"
  units <- list(
    time = "min",
    dosing = "mg (remimazolam and fentanyl)",
    concentration = "ng/mL (remimazolam, CNS 7054 and fentanyl; the UMSS outputs are category probabilities in 0-1)"
  )
  dosing <- c("central", "central_fentanyl")

  # What each ODE state holds, verified against the Supplemental Digital
  # Content 1 'final PKPD model for UMSS' control stream ($MODEL ...
  # COMP(FENT_C) COMP(FENT_P) COMP(EFFECT) COMP(EFFECT_FENT); $DES
  # 'DADT(12)= ke0*((A(1)/V1)-A(12))', 'DADT(13)= ke01*((A(10)/V10) -
  # A(13))'). The two effect states are concentrations, not amounts.
  compartmentData <- list(
    central = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE),
    transit1 = list(
      analyte = "CNS 7054",
      units = "mg (remimazolam equivalents)",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit2 = list(
      analyte = "CNS 7054",
      units = "mg (remimazolam equivalents)",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit3 = list(
      analyte = "CNS 7054",
      units = "mg (remimazolam equivalents)",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit4 = list(
      analyte = "CNS 7054",
      units = "mg (remimazolam equivalents)",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit5 = list(
      analyte = "CNS 7054",
      units = "mg (remimazolam equivalents)",
      specimen = "not applicable",
      verified = TRUE
    ),
    central_cns7054 = list(
      analyte = "CNS 7054",
      units = "mg (remimazolam equivalents)",
      specimen = "plasma",
      verified = TRUE
    ),
    central_fentanyl = list(analyte = "fentanyl", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_fentanyl = list(analyte = "fentanyl", units = "mg", specimen = "plasma", verified = TRUE),
    effect = list(analyte = "remimazolam", units = "ng/mL", specimen = "not applicable", verified = TRUE),
    effect_fentanyl = list(analyte = "fentanyl", units = "ng/mL", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling of the remimazolam / CNS 7054 clearances (exponent 0.75) and volumes (exponent 1) to 70 kg, and linear (per-kg) scaling of every fentanyl clearance and volume of the Ginsberg 1996 model (Supplemental Digital Content 1: 'CL_FENT = 0.019*WEIGHT', 'V10 = 0.32*WEIGHT'). Cohort mean 50.6 kg (SD 21.2, range 20.5-106) (Table 3).",
      source_name = "WEIGHT"
    )
  )

  # Tested on the remimazolam EC50 by likelihood-ratio test and not retained
  # (Results, Pharmacodynamic Model for UMSS): age (dOFV -0.54, P = 0.46),
  # sex (dOFV -1.17, P = 0.28), fentanyl coadministration (dOFV -2.37,
  # P = 0.12) and study procedure (dOFV -0.89, P = 0.83, 3 df). Fentanyl
  # itself enters the model mechanistically through its effect-site
  # concentration, so no fentanyl flag column is needed: the supplement's
  # 'IF(OPIOID.EQ.1)' guards only zero a fentanyl concentration that is
  # already zero without fentanyl doses.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age.",
      units = "years",
      type = "continuous",
      notes = "Tested on the remimazolam EC50; not retained (dOFV -0.54, P = 0.46)."
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male).",
      units = "(binary)",
      type = "binary",
      notes = "Tested on the remimazolam EC50; not retained (dOFV -1.17, P = 0.28). The paper does not state its coding."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 31L,
    n_studies = 1L,
    age_range = "6-18 years",
    age_mean = "12.3 years (SD 3.6)",
    weight_range = "20.5-106 kg",
    weight_mean = "50.6 kg (SD 21.2)",
    sex_female_pct = 100 * 12 / 31,
    disease_state = "Children and adolescents (ASA physical status I-III) undergoing elective diagnostic or therapeutic procedures under procedural sedation targeted to UMSS 2-3.",
    dose_range = "Intravenous remimazolam initial bolus 125-150 ug/kg (max 5-7 mg) or initial infusion 5.0-7.5 ug/kg/min, with top-up boluses (Table 2). Fentanyl in 21 of 31 patients: 50 ug then 25 ug supplements (12-18 years) or 1 ug/kg then 0.5 ug/kg (6-11 years), 2-5 min before the first remimazolam bolus.",
    regions = "United States",
    n_observations = "1,920 UMSS scores (1,394 with and 526 without fentanyl coadministration), of which 295 recorded after rescue sedation were excluded. UMSS was scored every minute from the first dose to the start of the procedure, every 3 min during the procedure, then every minute until arousable and every 5 min until discharge.",
    biomarkers = "University of Michigan Sedation Scale: 0 awake and alert, 1 minimally sedated, 2 moderately sedated, 3 deeply sedated, 4 unarousable (Supplemental Digital Content Table 1). No UMSS 4 was observed before rescue medication.",
    notes = "Fitted in NONMEM 7.4 with the Laplacian method using the individual pharmacokinetic parameter approach (individual remimazolam / CNS 7054 parameters fixed from the PK step). Fentanyl was not measured; its concentrations were predicted from the Ginsberg 1996 model and the individual dosing history. A model with competitive antagonism by CNS 7054 gave the same objective function value with an IC50 more than 10-fold above the highest observed CNS 7054 concentration and was not retained. Estimates from Table 5."
  )

  ini({
    # ---------------------------------------------------------------------
    # Remimazolam / CNS 7054 pharmacokinetic layer. The UMSS control stream
    # (Supplemental Digital Content 1) carries the PK step's final estimates
    # as '$THETA ... FIX' and reads each patient's post hoc etas from the
    # data (IETA1..IETA9), so none of these was re-estimated. They are
    # reproduced, fixed, from modellib('Colin_2025_remimazolam') so the
    # model can be simulated end to end; see that file for Table 4 traces.
    # ---------------------------------------------------------------------
    lvc <- fixed(0.223)
    label("Remimazolam central volume V1 for 70 kg (log L)") # SDC1 UMSS stream THETA(1) = 0.223 FIX; Table 4 V1 = 1.25 L/70 kg
    lvp <- fixed(2.38)
    label("Remimazolam fast peripheral volume V2 for 70 kg (log L)") # SDC1 UMSS stream THETA(2) = 2.38 FIX; Table 4 V2 = 10.8 L/70 kg
    lvp2 <- fixed(2.51)
    label("Remimazolam slow peripheral volume V3 for 70 kg (log L)") # SDC1 UMSS stream THETA(3) = 2.51 FIX; Table 4 V3 = 12.3 L/70 kg
    lcl <- fixed(-0.359)
    label("Remimazolam elimination clearance CL for 70 kg (log L/min)") # SDC1 UMSS stream THETA(4) = -0.359 FIX; Table 4 CL = 0.70 L/min/70 kg
    lq <- fixed(0.519)
    label("Remimazolam fast distributional clearance Q2 for 70 kg (log L/min)") # SDC1 UMSS stream THETA(5) = 0.519 FIX; Table 4 Q2 = 1.68 L/min/70 kg
    lq2 <- fixed(-1.47)
    label("Remimazolam slow distributional clearance Q3 for 70 kg (log L/min)") # SDC1 UMSS stream THETA(6) = -1.47 FIX; Table 4 Q3 = 0.23 L/min/70 kg
    lvc_cns7054 <- fixed(2.5)
    label("CNS 7054 volume of distribution Vm for 70 kg (log L)") # SDC1 UMSS stream THETA(7) = 2.5 FIX; Table 4 Vm = 12.2 L/70 kg
    lcl_cns7054 <- fixed(-2.69)
    label("CNS 7054 elimination clearance CLm for 70 kg (log L/min)") # SDC1 UMSS stream THETA(8) = -2.69 FIX; Table 4 CLm = 4.07 L/h/70 kg
    lmtt <- fixed(1.38)
    label("Mean transit time of CNS 7054 formation (log min)") # SDC1 UMSS stream THETA(9) = 1.38 FIX; Table 4 MTT = 3.97 min

    # ---------------------------------------------------------------------
    # Fentanyl, Ginsberg 1996 two-compartment paediatric model, coded in the
    # UMSS control stream with every clearance and volume proportional to
    # body weight and no between-subject variability (Methods: variability
    # 'was not included in the simulations for fentanyl'). Literature
    # constants, held fixed.
    # ---------------------------------------------------------------------
    lcl_fentanyl <- fixed(log(0.019))
    label("Fentanyl elimination clearance per kg body weight (log L/min/kg)") # SDC1 UMSS stream 'CL_FENT = 0.019*WEIGHT' (Ginsberg 1996)
    lvc_fentanyl <- fixed(log(0.32))
    label("Fentanyl central volume per kg body weight (log L/kg)") # SDC1 UMSS stream 'V10 = 0.32*WEIGHT'
    lvp_fentanyl <- fixed(log(1.49))
    label("Fentanyl peripheral volume per kg body weight (log L/kg)") # SDC1 UMSS stream 'V11 = 1.49*WEIGHT'
    lq_fentanyl <- fixed(log(0.036))
    label("Fentanyl intercompartmental clearance per kg body weight (log L/min/kg)") # SDC1 UMSS stream 'Q11 = 0.036*WEIGHT'

    # ---------------------------------------------------------------------
    # UMSS pharmacodynamics, Table 5. Every parameter is estimated on the
    # log scale (SDC1: 'ke0 = EXP(THETA(14))', 'TVEC50_2 = TVEC50_1 +
    # EXP(THETA(17))'). The EC50 of each higher UMSS category is the EC50
    # of the category below plus a positive increment (main-text equations
    # 1 and 2), so the categories cannot cross.
    # ---------------------------------------------------------------------
    lke0 <- log(0.230)
    label("Remimazolam effect-site equilibration rate constant ke0 (1/min)") # Table 5 ke0RMZ = 0.230 1/min (0.20-0.27)
    lhill <- log(3.56)
    label("Hill coefficient shared by the three UMSS thresholds (unitless)") # Table 5 gamma = 3.56 (3.18-3.96)
    lec50 <- log(206)
    label("Remimazolam effect-site EC50 for UMSS >= 1 (ng/mL)") # Table 5 EC50RMZ,1 = 206 ng/mL (151-281)
    ldec50_umss2 <- log(224)
    label("Increment of the remimazolam EC50 from UMSS >= 1 to UMSS >= 2 (ng/mL)") # Table 5 dEC50RMZ,2 = 224 ng/mL (162-309)
    ldec50_umss3 <- log(347)
    label("Increment of the remimazolam EC50 from UMSS >= 2 to UMSS >= 3 (ng/mL)") # Table 5 dEC50RMZ,3 = 347 ng/mL (249-488)
    lke0_fentanyl <- log(0.125)
    label("Fentanyl effect-site equilibration rate constant (1/min)") # Table 5 ke0FENT = 0.125 1/min (0.08-0.20)
    lec50_fentanyl <- log(0.560)
    label("Fentanyl effect-site EC50 for UMSS >= 1 (ng/mL)") # Table 5 EC50FENT,1 = 0.560 ng/mL (0.49-0.73)
    ldec50_umss2_fentanyl <- log(1.60)
    label("Increment of the fentanyl EC50 from UMSS >= 1 to UMSS >= 2 (ng/mL)") # Table 5 dEC50FENT,2 = 1.60 ng/mL (1.09-2.60)
    ldec50_umss3_fentanyl <- log(4.14)
    label("Increment of the fentanyl EC50 from UMSS >= 2 to UMSS >= 3 (ng/mL)") # Table 5 dEC50FENT,3 = 4.14 ng/mL (2.59-41.7)

    # Minto 2000 interaction coefficient on U50, fixed to 0 (additive
    # interaction): SDC1 '0 FIX ; INT; rel.' (THETA(23)); Discussion:
    # 'assuming beta2,U50 = beta3,U50 = beta4,U50 = 0'.
    beta2_u50 <- fixed(0)
    label("Minto interaction coefficient on the normalised potency U50 (unitless)") # SDC1 UMSS stream THETA(23) = 0 FIX

    # ---------------------------------------------------------------------
    # Between-subject variability. PK layer from Table 4 (CV% = sqrt of the
    # NONMEM variance, Table 4 footnote), fixed here as in the UMSS stream.
    # PD: Table 5 IIV on the remimazolam EC50 = 68.2%, read on the same
    # sqrt-of-variance convention as Table 4. The one eta scales all three
    # remimazolam EC50s (SDC1 'EC50_k = TVEC50_k * EXP(ETA(1))'). The Hill
    # coefficient and the fentanyl EC50s carry an IIV 'fixed to 10%'
    # ($OMEGA 1E-2 FIX 'BSV_Gamma', 'BSV_C50_FENT').
    # ---------------------------------------------------------------------
    etalvc ~ fixed(0.01) # SDC1 PK stream $OMEGA 1E-2 FIX 'IIV_V1'
    etalvp ~ fixed(0.853776) # Table 4 IIV V2 = 92.4%; 0.924^2
    etalvp2 ~ fixed(0.01) # SDC1 PK stream $OMEGA 1E-2 FIX 'IIV_V3'
    etalcl ~ fixed(0.398161) # Table 4 IIV CL = 63.1%; 0.631^2
    etalq ~ fixed(0.01) # SDC1 PK stream $OMEGA 1E-2 FIX 'IIV_Q2'
    etalq2 ~ fixed(0.01) # SDC1 PK stream $OMEGA 1E-2 FIX 'IIV_Q3'
    etalvc_cns7054 ~ fixed(0.245025) # Table 4 IIV Vm = 49.5%; 0.495^2
    etalcl_cns7054 ~ fixed(0.417316) # Table 4 IIV CLm = 64.6%; 0.646^2
    etalmtt ~ fixed(0.01) # SDC1 PK stream $OMEGA 1E-2 FIX 'BSV_MTT'
    etalec50 ~ 0.465124 # Table 5 IIV EC50RMZ,1 = 68.2% (55.9-120.1), shrinkage 32%; 0.682^2
    etalhill ~ fixed(0.01) # SDC1 UMSS stream $OMEGA 1E-2 FIX 'BSV_Gamma'
    etalec50_fentanyl ~ fixed(0.01) # SDC1 UMSS stream $OMEGA 1E-2 FIX 'BSV_C50_FENT'

    # ---------------------------------------------------------------------
    # Concentration residual errors from the PK step, which the UMSS stream
    # also carries as FIX thetas (THETA(10)-(13)). The UMSS likelihood is
    # the exact categorical probability (SDC1 'Y = DV0*P0 + ... + DV3*P3')
    # and has no residual error; the concentration endpoints are kept so
    # the PK layer can be simulated with observation noise.
    # ---------------------------------------------------------------------
    propSd <- fixed(0.411)
    label("Proportional residual error on remimazolam (fraction)") # Table 4 PropRMZ = 41.1%; SDC1 UMSS stream THETA(10) = -0.888 FIX
    addSd <- fixed(1)
    label("Additive residual error on remimazolam (ng/mL)") # SDC1 THETA(11) = 0 FIX ';; SD = 1 ng.mL-1'
    propSd_cns7054 <- fixed(0.175)
    label("Proportional residual error on CNS 7054 (fraction)") # Table 4 PropCNS7054 = 17.5%; SDC1 UMSS stream THETA(12) = -1.74 FIX
    addSd_cns7054 <- fixed(22.9)
    label("Additive residual error on CNS 7054 (ng/mL)") # Table 4 AddCNS7054 = 22.9 ng/mL; SDC1 UMSS stream THETA(13) = 3.13 FIX
  })

  model({
    # Remimazolam / CNS 7054 (as modellib('Colin_2025_remimazolam')).
    fsize <- WT / 70

    vc <- exp(lvc + etalvc) * fsize
    vp <- exp(lvp + etalvp) * fsize
    vp2 <- exp(lvp2 + etalvp2) * fsize
    cl <- exp(lcl + etalcl) * fsize^0.75
    q <- exp(lq + etalq) * fsize^0.75
    q2 <- exp(lq2 + etalq2) * fsize^0.75
    vc_cns7054 <- exp(lvc_cns7054 + etalvc_cns7054) * fsize
    cl_cns7054 <- exp(lcl_cns7054 + etalcl_cns7054) * fsize^0.75
    mtt <- exp(lmtt + etalmtt)
    ktr <- (5 + 1) / mtt

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2
    kel_cns7054 <- cl_cns7054 / vc_cns7054

    # Fentanyl, Ginsberg 1996 (per-kg constants times body weight).
    cl_fentanyl <- exp(lcl_fentanyl) * WT
    vc_fentanyl <- exp(lvc_fentanyl) * WT
    vp_fentanyl <- exp(lvp_fentanyl) * WT
    q_fentanyl <- exp(lq_fentanyl) * WT
    kel_fentanyl <- cl_fentanyl / vc_fentanyl
    k12_fentanyl <- q_fentanyl / vc_fentanyl
    k21_fentanyl <- q_fentanyl / vp_fentanyl

    ke0 <- exp(lke0)
    ke0_fentanyl <- exp(lke0_fentanyl)

    d/dt(central) <- -k12 * central + k21 * peripheral1 - k13 * central +
      k31 * peripheral2 - kel * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2
    d/dt(transit1) <- kel * central - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(transit4) <- ktr * transit3 - ktr * transit4
    d/dt(transit5) <- ktr * transit4 - ktr * transit5
    d/dt(central_cns7054) <- ktr * transit5 - kel_cns7054 * central_cns7054
    d/dt(central_fentanyl) <- -k12_fentanyl * central_fentanyl +
      k21_fentanyl * peripheral1_fentanyl - kel_fentanyl * central_fentanyl
    d/dt(peripheral1_fentanyl) <- k12_fentanyl * central_fentanyl -
      k21_fentanyl * peripheral1_fentanyl

    # Plasma concentrations in ng/mL (mg/L * 1000). CNS 7054 carries the
    # 439.3 / 425.3 factor exactly as the supplement codes it.
    Cc <- central / vc * 1000
    Cc_cns7054 <- central_cns7054 / vc_cns7054 * 1000 * 439.3 / 425.3
    Cc_fentanyl <- central_fentanyl / vc_fentanyl * 1000

    # Effect compartments driven by the plasma concentrations; the states
    # are effect-site concentrations in ng/mL.
    d/dt(effect) <- ke0 * (Cc - effect)
    d/dt(effect_fentanyl) <- ke0_fentanyl * (Cc_fentanyl - effect_fentanyl)

    # Category-specific EC50s (main-text equations 1-2; SDC1 TVEC50_k).
    ec50_umss1 <- exp(lec50 + etalec50)
    ec50_umss2 <- (exp(lec50) + exp(ldec50_umss2)) * exp(etalec50)
    ec50_umss3 <- (exp(lec50) + exp(ldec50_umss2) + exp(ldec50_umss3)) *
      exp(etalec50)
    ec50_umss1_fentanyl <- exp(lec50_fentanyl + etalec50_fentanyl)
    ec50_umss2_fentanyl <- (exp(lec50_fentanyl) + exp(ldec50_umss2_fentanyl)) *
      exp(etalec50_fentanyl)
    ec50_umss3_fentanyl <- (exp(lec50_fentanyl) + exp(ldec50_umss2_fentanyl) +
      exp(ldec50_umss3_fentanyl)) * exp(etalec50_fentanyl)
    hill <- exp(lhill + etalhill)

    # Normalised effect-site concentrations (main-text equations 3-4).
    u_umss1 <- effect / ec50_umss1
    u_umss2 <- effect / ec50_umss2
    u_umss3 <- effect / ec50_umss3
    u_umss1_fentanyl <- effect_fentanyl / ec50_umss1_fentanyl
    u_umss2_fentanyl <- effect_fentanyl / ec50_umss2_fentanyl
    u_umss3_fentanyl <- effect_fentanyl / ec50_umss3_fentanyl

    # Minto U50 for the fentanyl fraction theta = U_FENT / (U_RMZ + U_FENT):
    # U50 = 1 - beta2 * theta + beta2 * theta^2 (SDC1 U50_k), equal to 1
    # for the fitted additive interaction. theta is set to 0 when neither
    # drug is present.
    theta_umss1 <- 0
    theta_umss2 <- 0
    theta_umss3 <- 0
    if (u_umss1 + u_umss1_fentanyl > 0) {
      theta_umss1 <- u_umss1_fentanyl / (u_umss1 + u_umss1_fentanyl)
      theta_umss2 <- u_umss2_fentanyl / (u_umss2 + u_umss2_fentanyl)
      theta_umss3 <- u_umss3_fentanyl / (u_umss3 + u_umss3_fentanyl)
    }
    u50_umss1 <- 1 - beta2_u50 * theta_umss1 + beta2_u50 * theta_umss1^2
    u50_umss2 <- 1 - beta2_u50 * theta_umss2 + beta2_u50 * theta_umss2^2
    u50_umss3 <- 1 - beta2_u50 * theta_umss3 + beta2_u50 * theta_umss3^2

    # Cumulative probabilities P(UMSS >= k), main-text equation 5 (SDC1
    # DRUG1..DRUG3 = PGE1..PGE3).
    x_umss1 <- ((u_umss1 + u_umss1_fentanyl) / u50_umss1)^hill
    x_umss2 <- ((u_umss2 + u_umss2_fentanyl) / u50_umss2)^hill
    x_umss3 <- ((u_umss3 + u_umss3_fentanyl) / u50_umss3)^hill
    prob_umss_ge1 <- x_umss1 / (1 + x_umss1)
    prob_umss_ge2 <- x_umss2 / (1 + x_umss2)
    prob_umss_ge3 <- x_umss3 / (1 + x_umss3)

    # Category probabilities (SDC1 P0..P3); P(UMSS = 0) = 1 at baseline.
    prob_umss0 <- 1 - prob_umss_ge1
    prob_umss1 <- prob_umss_ge1 - prob_umss_ge2
    prob_umss2 <- prob_umss_ge2 - prob_umss_ge3
    prob_umss3 <- prob_umss_ge3

    # Expected score, the supplement's IPRED = 3*P3 + 2*P2 + 1*P1.
    umss_expected <- 3 * prob_umss3 + 2 * prob_umss2 + 1 * prob_umss1

    Cc ~ prop(propSd) + add(addSd)
    Cc_cns7054 ~ prop(propSd_cns7054) + add(addSd_cns7054)
    # Ordinal endpoint: simulated category k is UMSS score k - 1 (category
    # 4, score 3, carries 1 - the sum of the three listed probabilities).
    umss ~ c(prob_umss0, prob_umss1, prob_umss2)
  })
}
