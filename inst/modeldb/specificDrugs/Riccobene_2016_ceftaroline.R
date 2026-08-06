Riccobene_2016_ceftaroline <- function() {
  description <- "Joint two-compartment ceftaroline fosamil (prodrug) plus three-compartment ceftaroline (active metabolite) population PK model with an algebraic epithelial-lining-fluid (ELF) partition coefficient, developed on plasma and bronchoalveolar-lavage ELF data from 50 healthy adults given 600 mg ceftaroline fosamil as a 1-h IV infusion q12h or q8h to steady state (Riccobene 2016). Ceftaroline fosamil is assumed to be converted completely to ceftaroline, so the whole prodrug elimination clearance CLcf enters the ceftaroline central compartment. ELF concentrations are not a distribution compartment: the parallel decline of plasma and ELF made an ELF compartment unidentifiable, so ELF is modelled as the ceftaroline plasma concentration scaled by a partition coefficient (0.193) carrying inter-individual variability. Body weight enters allometrically on all clearances (0.75) and volumes (1); BSA-normalised creatinine clearance below 80 mL/min/1.73 m2, age, hemodialysis status and healthy-versus-patient status act on ceftaroline clearance, and healthy-versus-patient status acts on ceftaroline central volume. An intramuscular depot (first-order ka, bioavailability fixed to 1) is carried from the upstream adult model; the ELF study itself used IV dosing only."
  reference   <- paste(
    "Riccobene TA, Pushkin R, Jandourek A, Knebel W, Khariton T.",
    "Penetration of ceftaroline into the epithelial lining fluid of healthy adult subjects.",
    "Antimicrob Agents Chemother. 2016;60(10):5849-5857. doi:10.1128/AAC.02755-15.",
    "Parameter estimates are from Supplemental Table 1 and the model equations from",
    "Supplemental Equation 1 of that article's supplemental material.",
    "Structural model and the fixed covariate effects were carried from the upstream adult",
    "ceftaroline fosamil / ceftaroline population PK analysis of 10 phase 1, 1 phase 2 and",
    "4 phase 3 studies cited as reference 14 of Riccobene 2016 (Riccobene T, Khariton T,",
    "Knebel W, O'Neal T, Ghahramani P. 2013. Abstr 23rd Eur Congr Clin Microbiol Infect Dis,",
    "abstr P902); that abstract is not in nlmixr2lib.",
    sep = " "
  )
  vignette    <- "Riccobene_2016_ceftaroline"

  units       <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling of every clearance (exponent 0.75) and every volume (exponent 1) of both ceftaroline fosamil and ceftaroline, referenced to 70 kg. Supplemental Table 1 prints the scaling terms as (WT/70)^0.75 and (WT/70)^1 under each structural parameter, and Supplemental Equation 1 writes them as log(WT/70)*0.75 and log(WT/70) additive terms inside exp(). Weight ranged 58-102 kg in the ELF study (Results, Population pharmacokinetics in the lung).",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Age at study entry.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on ceftaroline clearance, (AGE/36)^-0.278, applied to every subject (Supplemental Equation 1 COV5 = log(AGE/36)*theta13, with no conditional gate). The exponent was fixed from the upstream adult model because the ELF study enrolled nobody over 45 years (Methods, Population pharmacokinetics in the lung). Ages ranged 20-45 years.",
      source_name        = "AGE"
    ),
    CRCL = list(
      description        = "Creatinine clearance normalised by body surface area.",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on ceftaroline clearance, (CRCL/80)^0.508, active ONLY below the 80 mL/min/1.73 m^2 reference and only in subjects not on hemodialysis; above 80 the term is 1. Encoded with min(CRCL, 80) so the term saturates at the reference. The exponent was fixed from the upstream adult model because no subject in the ELF study had a CRCL below 80 (Methods, Population pharmacokinetics in the lung). Source column nCLCR in the paper text and nCRCL in the supplement.",
      source_name        = "nCRCL"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-participant indicator (1 = healthy adult, 0 = patient with an infection).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (patient with an infection; the pooled phase 2 / phase 3 infected-patient cohort of the upstream adult model)",
      notes              = "Multiplicative power-form effects 3.32^DIS_HEALTHY on ceftaroline clearance and 3.67^DIS_HEALTHY on ceftaroline central volume. The source column is PAT (1 = patient, 0 = healthy) and is re-expressed here as DIS_HEALTHY = 1 - PAT, so lcl_ceftaroline and lvc_ceftaroline carry the PATIENT typical values printed in Supplemental Table 1 and the multipliers restore the healthy-adult values. See the vignette Errata: the supplement's footnote orientation is contradicted by the paper's own Table 2 noncompartmental results, which only reproduce when both multipliers are active for the ELF study's all-healthy cohort.",
      source_name        = "PAT"
    ),
    RRT_HEMODIAL_STATUS = list(
      description        = "Intermittent-hemodialysis treatment-status indicator (1 = subject is on an intermittent hemodialysis programme, 0 = not).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not on hemodialysis)",
      notes              = "Multiplicative power-form effect 0.372^RRT_HEMODIAL_STATUS on ceftaroline clearance between dialysis sessions, and it switches off the CRCL term. Supplemental Table 1 labels the source column ESRD (end-stage renal disease, 1 = yes) and Supplemental Equation 1 scopes the term to 'dialysis patients during non-dialysis periods', so within this model the ESRD flag and the on-a-hemodialysis-programme flag are operationally the same column; the canonical RRT name is used because the term's scope is dialysis, not renal disease generally. No subject in the ELF study had end-stage renal disease, so the coefficient is fixed from the upstream adult model.",
      source_name        = "ESRD"
    ),
    RRT_HEMODIAL_ACTIVE = list(
      description        = "Hemodialysis-active indicator (1 while a dialysis session is running, 0 in the interdialytic interval and in non-dialysed subjects).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no dialysis session running)",
      notes              = "Time-varying per-session gate. Supplemental Equation 1 REPLACES ceftaroline clearance entirely during a session (CLc_i = exp(theta14) = 9.97 L/h, with no covariate or eta terms) rather than adding a dialyser arm to the body clearance, so the model composes the two arms as cl_interdialytic*(1 - RRT_HEMODIAL_ACTIVE) + cl_hemodialysis*RRT_HEMODIAL_ACTIVE. No subject in the ELF study was dialysed, so the value is fixed from the upstream adult model.",
      source_name        = "(dialysis period flag, unnamed in the supplement)"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "ceftaroline fosamil", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "ceftaroline fosamil", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "ceftaroline fosamil", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    central_ceftaroline = list(
      analyte = "ceftaroline", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1_ceftaroline = list(
      analyte = "ceftaroline", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    peripheral2_ceftaroline = list(
      analyte = "ceftaroline", units = "mg",
      specimen = "plasma", verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 50L,
    n_studies      = 1L,
    age_range      = "20-45 years (enrolled 19-45; mean 34.6 years in the q12h arm and 33.1 years in the q8h arm)",
    weight_range   = "58-102 kg",
    sex_female_pct = 100 * 8 / 50,
    race_ethnicity = c(White = 81.1, Black = 13.2, Asian = 3.8, Other = 1.9),
    disease_state  = "Healthy adult nonsmokers, body mass index 18-30 kg/m^2, no clinically significant disease; none had end-stage renal disease, none were on dialysis, none had a BSA-normalised creatinine clearance below 80 mL/min/1.73 m^2 and none were over 45 years of age.",
    dose_range     = "600 mg ceftaroline fosamil as a 1-h IV infusion, q12h or q8h for 3 days plus a single dose on day 4",
    regions        = "Single centre in the United States (Pulmonary Associates, Phoenix, AZ)",
    n_observations = "856 measurable plasma concentrations (210 ceftaroline fosamil, 646 ceftaroline) and 49 measurable ELF concentrations (6 ceftaroline fosamil, 43 ceftaroline)",
    notes          = "Phase 1, open-label, multiple-dose study; 53 subjects enrolled (Table 1 demographics: 26 q12h, 27 q8h), 50 completed and contributed to the population PK analysis (25 per arm). Race and sex percentages above are computed from the 53 enrolled subjects of Table 1 because Table 1 is the only demographic breakdown reported; the 50-subject analysis population is reported only as 42 males and 8 females (Results, Population pharmacokinetics in the lung), which is the sex split used here. Each subject underwent bronchoalveolar lavage at exactly one of five post-dose times, so ELF is a sparse, composite profile. The structural model and all covariate effects come from an upstream pooled adult analysis (10 phase 1, 1 phase 2, 4 phase 3 studies); effects that the ELF data carried no information on were fixed to their upstream values."
  )

  ini({
    # ---------------------------------------------------------------------
    # Ceftaroline fosamil (parent prodrug) -- Supplemental Table 1, Fixed
    # Effects. The table prints back-transformed values (exp(theta)) in
    # natural units; Supplemental Equation 1 writes each parameter as
    # exp(theta + log(WT/70)*exponent + eta), so the printed value is the
    # typical value at WT = 70 kg.
    # ---------------------------------------------------------------------
    lcl     <- log(238);         label("Ceftaroline fosamil clearance CLcf at 70 kg (L/h)")                    # Suppl Table 1: CLcf (theta1) = 238 L/h (95% CI 206, 274)
    lvc     <- fixed(log(6.96)); label("Ceftaroline fosamil central volume Vccf at 70 kg (L)")                  # Suppl Table 1: Vccf (theta2) = 6.96 L [FIXED]
    lq      <- log(97.4);        label("Ceftaroline fosamil intercompartmental clearance Qcf at 70 kg (L/h)")   # Suppl Table 1: Qcf (theta3) = 97.4 L/h (95% CI 57.3, 165)
    lvp     <- log(6.23);        label("Ceftaroline fosamil peripheral volume Vpcf at 70 kg (L)")               # Suppl Table 1: Vpcf (theta4) = 6.23 L (95% CI 4.26, 9.1)

    # Intramuscular route: carried from the upstream adult model. The ELF
    # study dosed intravenously only, so neither parameter is informed by
    # these data and both are fixed. Supplemental Equation 1 shows an
    # eta on ka1cf, but Supplemental Table 1 reports no variance for it,
    # so no IIV is carried here (see vignette Errata).
    lka     <- fixed(log(0.549)); label("Ceftaroline fosamil intramuscular absorption rate kacf (1/h)")         # Suppl Table 1: kacf (theta5) = 0.549 1/h [FIXED]
    lfdepot <- fixed(log(1));     label("Ceftaroline fosamil intramuscular bioavailability FIMcf (fraction)")   # Suppl Table 1: FIMcf (theta6) = 1 [FIXED]

    # ---------------------------------------------------------------------
    # Ceftaroline (active metabolite) -- Supplemental Table 1.
    #
    # IMPORTANT: lcl_ceftaroline and lvc_ceftaroline carry the PATIENT
    # typical values. Supplemental Table 1 prints CLc = 3.76 L/h and
    # Vcc = 3.18 L alongside the multiplicative terms theta16^PAT = 3.32
    # (on CLc) and theta15^PAT = 3.67 (on Vcc). Applying no multiplier to
    # this all-healthy cohort overpredicts steady-state AUC0-tau by 3.5x
    # (157 vs the 45.0-53.0 mg*h/L of the paper's own Table 2) and Cmax by
    # 2.4x; applying BOTH multipliers reproduces Table 2 in both dose arms
    # (AUC0-tau 47.3, Cmax 22.2-22.6 mg/L). Applying only one of the two
    # fails (CL only -> Cmax 29.5; Vc only -> AUC 157), so the pair is
    # identified jointly. The printed values are therefore the patient
    # reference and the multipliers restore the healthy state, which is
    # the canonical DIS_HEALTHY orientation. See vignette Errata.
    # ---------------------------------------------------------------------
    lcl_ceftaroline  <- log(3.76); label("Ceftaroline clearance CLc at 70 kg, age 36 y, in patients (L/h)")                      # Suppl Table 1: CLc (theta7) = 3.76 L/h (95% CI 3.6, 3.93)
    lvc_ceftaroline  <- log(3.18); label("Ceftaroline central volume Vcc at 70 kg in patients (L)")                              # Suppl Table 1: Vcc (theta8) = 3.18 L (95% CI 2.72, 3.73)
    lq_ceftaroline   <- log(3.71); label("Ceftaroline intercompartmental clearance to peripheral1, Q1c at 70 kg (L/h)")          # Suppl Table 1: Q1c (theta9) = 3.71 L/h (95% CI 2.66, 5.19)
    lvp_ceftaroline  <- log(7.14); label("Ceftaroline first peripheral volume Vp1c at 70 kg (L)")                                # Suppl Table 1: Vp1c (theta10) = 7.14 L (95% CI 6.13, 8.32)
    lq2_ceftaroline  <- log(18.6); label("Ceftaroline intercompartmental clearance to peripheral2, Q2c at 70 kg (L/h)")          # Suppl Table 1: Q2c (theta18) = 18.6 L/h (95% CI 13.7, 25.3)
    lvp2_ceftaroline <- log(6.76); label("Ceftaroline second peripheral volume Vp2c at 70 kg (L)")                               # Suppl Table 1: Vp2c (theta19) = 6.76 L (95% CI 5.71, 8)

    # Dialysis-active ceftaroline clearance. Supplemental Equation 1
    # REPLACES CLc with this value while a session runs, rather than
    # adding a dialyser arm on top of body clearance.
    lcl_hemodialysis_ceftaroline <- fixed(log(9.97)); label("Ceftaroline clearance during a hemodialysis session, CLdial (L/h)") # Suppl Table 1: CLdial (theta14) = 9.97 L/h [FIXED]

    # ---------------------------------------------------------------------
    # Epithelial lining fluid partition coefficient.
    # Supplemental Equation 1: PC = theta17 * exp(etaPC) and the ELF
    # observation is Cij = Chat_ij * PC * (1 + eps_p3ij), i.e. ELF is an
    # algebraic scaling of the ceftaroline plasma concentration, not a
    # distribution compartment (Results: the parallel decline of plasma
    # and ELF made an ELF compartment unidentifiable).
    # ---------------------------------------------------------------------
    lrelf <- log(0.193); label("ELF / total-plasma ceftaroline concentration ratio PCELF (unitless)")                            # Suppl Table 1: PCELF (theta17) = 0.193 (95% CI 0.161, 0.226); main text Results quotes 0.193 (0.171, 0.215)

    # ---------------------------------------------------------------------
    # Covariate effects. All are fixed: Supplemental Table 1 marks each
    # [FIXED], and the Methods state that effects the ELF data carried no
    # information on were held at their upstream-model values.
    # ---------------------------------------------------------------------
    # Shared allometric exponents: 0.75 on every clearance-type parameter
    # (CLcf, Qcf, CLc, Q1c, Q2c) and 1 on every volume (Vccf, Vpcf, Vcc,
    # Vp1c, Vp2c), as printed under each row of Supplemental Table 1.
    e_wt_cl_q  <- fixed(0.75); label("Allometric exponent shared by all clearances (unitless)")                                  # Suppl Table 1: (WT/70)^0.75 under CLcf, Qcf, CLc, Q1c, Q2c
    e_wt_vc_vp <- fixed(1);    label("Allometric exponent shared by all volumes (unitless)")                                     # Suppl Table 1: (WT/70)^1 under Vccf, Vpcf, Vcc, Vp1c, Vp2c

    e_crcl_cl_ceftaroline        <- fixed(0.508);  label("Exponent of (CRCL/80) on ceftaroline CL, below 80 mL/min/1.73 m^2")    # Suppl Table 1: (nCRCL/80)^theta11, theta11 = 0.508 [FIXED]
    e_age_cl_ceftaroline         <- fixed(-0.278); label("Exponent of (AGE/36) on ceftaroline CL")                               # Suppl Table 1: (AGE/36)^theta13, theta13 = -0.278 [FIXED]
    e_dis_healthy_cl_ceftaroline <- fixed(3.32);   label("Healthy-participant multiplier on ceftaroline CL (power-form base)")    # Suppl Table 1: theta16^PAT = 3.32 [FIXED], re-expressed onto DIS_HEALTHY = 1 - PAT
    e_dis_healthy_vc_ceftaroline <- fixed(3.67);   label("Healthy-participant multiplier on ceftaroline Vc (power-form base)")    # Suppl Table 1: theta15^PAT = 3.67 [FIXED], re-expressed onto DIS_HEALTHY = 1 - PAT
    e_rrt_hemodial_status_cl_ceftaroline <- fixed(0.372); label("Hemodialysis multiplier on interdialytic ceftaroline CL (power-form base)") # Suppl Table 1: theta12^ESRD = 0.372 [FIXED]

    # ---------------------------------------------------------------------
    # Inter-individual variability -- Supplemental Table 1,
    # "Inter-individual variability". The listed VAR / COV entries form a
    # complete 4x4 lower triangle in the order CLcf, Vccf, CLc, Vcc, plus
    # two diagonal-only entries for Vp1c and PCELF. Values are log-scale
    # variances: the printed %CV equals 100*sqrt(VAR) for every row (e.g.
    # VAR(CLcf) = 0.0604 -> 24.6%), so they are carried unchanged.
    # ---------------------------------------------------------------------
    # Row-by-row source trace for the block below (comments cannot sit
    # inside the c(...) itself -- rxode2's label parser rejects them):
    #   row 1 (CLcf): VAR(CLcf)      = 0.0604   (%CV 24.6; 95% CI -0.0377, 0.158)
    #   row 2 (Vccf): COV(CLcf,Vccf) = -0.0661  (r -0.894)
    #                 VAR(Vccf)      = 0.0906   (%CV 30.1; 95% CI -0.502, 0.684)
    #   row 3 (CLc):  COV(CLcf,CLc)  = 0.000937 (r 0.0318)
    #                 COV(Vccf,CLc)  = 0.0151   (r 0.419)
    #                 VAR(CLc)       = 0.0143   (%CV 12; 95% CI 0.0056, 0.0231)
    #   row 4 (Vcc):  COV(CLcf,Vcc)  = 0.00584  (r 0.128)
    #                 COV(Vccf,Vcc)  = 0.0112   (r 0.199)
    #                 COV(CLc,Vcc)   = 0.0157   (r 0.702)
    #                 VAR(Vcc)       = 0.0346   (%CV 18.6; 95% CI 0.00572, 0.0636)
    etalcl + etalvc + etalcl_ceftaroline + etalvc_ceftaroline ~ c(
      0.0604,
      -0.0661, 0.0906,
      0.000937, 0.0151, 0.0143,
      0.00584, 0.0112, 0.0157, 0.0346
    )
    etalvp_ceftaroline ~ 0.02                  # VAR(Vp1c)  = 0.02  (%CV 14.1; 95% CI 0.00283, 0.0371)
    etalrelf           ~ 0.128                 # VAR(PCELF) = 0.128 (%CV 35.8; 95% CI 0.0458, 0.211)

    # ---------------------------------------------------------------------
    # Residual variability -- Supplemental Table 1, "Residual
    # variability". The table reports VARIANCES; the parenthetical %CV /
    # SD columns are their square roots (e.g. propCF 0.0658 -> %CV 25.7,
    # addCF 0.00237 -> SD 0.0487), so the SDs below are sqrt(variance).
    # The reported COV(propCF, propC) = 0.00544 (r = 0.401), a correlation
    # between the two plasma proportional residuals, has no nlmixr2
    # representation and is not carried (see vignette Errata).
    # ---------------------------------------------------------------------
    propSd             <- 0.2565;        label("Ceftaroline fosamil proportional residual SD (fraction)")   # Suppl Table 1: propCF = 0.0658 variance (%CV 25.7) -> sqrt(0.0658) = 0.2565
    addSd              <- 0.0487;        label("Ceftaroline fosamil additive residual SD (mg/L)")           # Suppl Table 1: addCF  = 0.00237 variance (SD 0.0487)
    propSd_ceftaroline <- 0.0528;        label("Ceftaroline proportional residual SD (fraction)")           # Suppl Table 1: propC  = 0.00279 variance (%CV 5.28) -> sqrt(0.00279) = 0.0528
    addSd_ceftaroline  <- 0.102;         label("Ceftaroline additive residual SD (mg/L)")                   # Suppl Table 1: addC   = 0.0105 variance (SD 0.102)
    propSd_Celf        <- fixed(0.0632); label("Ceftaroline ELF proportional residual SD (fraction)")       # Suppl Table 1: propELF = 0.004 variance [FIXED] -> sqrt(0.004) = 0.0632
  })

  model({
    # ----- Reference covariate values (Supplemental Table 1 / Equation 1)
    ref_wt   <- 70   # kg
    ref_age  <- 36   # years
    ref_crcl <- 80   # mL/min/1.73 m^2

    allom_cl <- (WT / ref_wt)^e_wt_cl_q
    allom_v  <- (WT / ref_wt)^e_wt_vc_vp

    # ----- Derived covariate terms on ceftaroline clearance
    # Supplemental Equation 1: COV3 = log(nCRCL/80)*theta11 for
    # non-dialysis subjects with nCRCL < 80, and COV3 = 0 otherwise. The
    # min() saturates the term at the 80 reference so it is exactly 1
    # above it, and zeroing the exponent for dialysed subjects reproduces
    # the "non-dialysis patients" gate without branching.
    crcl_cl <- (min(CRCL, ref_crcl) / ref_crcl)^(e_crcl_cl_ceftaroline * (1 - RRT_HEMODIAL_STATUS))
    # COV5 = log(AGE/36)*theta13, applied to every subject.
    age_cl  <- (AGE / ref_age)^e_age_cl_ceftaroline

    # ----- Individual parameters: ceftaroline fosamil (parent prodrug)
    cl <- exp(lcl + etalcl) * allom_cl
    vc <- exp(lvc + etalvc) * allom_v
    q  <- exp(lq)  * allom_cl
    vp <- exp(lvp) * allom_v
    ka <- exp(lka)

    # ----- Individual parameters: ceftaroline (active metabolite)
    # Clearance switches wholesale to the dialysis-session value while
    # RRT_HEMODIAL_ACTIVE is 1 (Supplemental Equation 1: "for dialysis
    # patients during dialysis, CLc_i = exp(theta14)"), and otherwise
    # carries the covariate cascade. e_dis_healthy_* multipliers move the
    # patient-reference typicals to the healthy state.
    cl_ceftaroline <-
      exp(lcl_ceftaroline + etalcl_ceftaroline) * allom_cl * crcl_cl * age_cl *
        e_rrt_hemodial_status_cl_ceftaroline^RRT_HEMODIAL_STATUS *
        e_dis_healthy_cl_ceftaroline^DIS_HEALTHY *
        (1 - RRT_HEMODIAL_ACTIVE) +
      exp(lcl_hemodialysis_ceftaroline) * RRT_HEMODIAL_ACTIVE
    vc_ceftaroline  <- exp(lvc_ceftaroline + etalvc_ceftaroline) * allom_v *
      e_dis_healthy_vc_ceftaroline^DIS_HEALTHY
    q_ceftaroline   <- exp(lq_ceftaroline)                        * allom_cl
    vp_ceftaroline  <- exp(lvp_ceftaroline + etalvp_ceftaroline)  * allom_v
    q2_ceftaroline  <- exp(lq2_ceftaroline)                       * allom_cl
    vp2_ceftaroline <- exp(lvp2_ceftaroline)                      * allom_v

    relf <- exp(lrelf + etalrelf)

    # ----- Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    kel_ceftaroline <- cl_ceftaroline / vc_ceftaroline
    k12_ceftaroline <- q_ceftaroline  / vc_ceftaroline
    k21_ceftaroline <- q_ceftaroline  / vp_ceftaroline
    k13_ceftaroline <- q2_ceftaroline / vc_ceftaroline
    k31_ceftaroline <- q2_ceftaroline / vp2_ceftaroline

    # ----- ODE system
    # Ceftaroline fosamil is assumed to be converted completely to
    # ceftaroline, so its whole elimination clearance CLcf becomes the
    # formation flux into the ceftaroline central compartment (the paper
    # reports rapid, complete conversion and Supplemental Equation 1
    # carries no formation fraction). Mass is transferred 1:1 in mg; the
    # supplement states no molar conversion between the two analytes.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    d/dt(central_ceftaroline) <-
      kel * central -
      kel_ceftaroline * central_ceftaroline -
      k12_ceftaroline * central_ceftaroline + k21_ceftaroline * peripheral1_ceftaroline -
      k13_ceftaroline * central_ceftaroline + k31_ceftaroline * peripheral2_ceftaroline
    d/dt(peripheral1_ceftaroline) <- k12_ceftaroline * central_ceftaroline - k21_ceftaroline * peripheral1_ceftaroline
    d/dt(peripheral2_ceftaroline) <- k13_ceftaroline * central_ceftaroline - k31_ceftaroline * peripheral2_ceftaroline

    # ----- Intramuscular bioavailability
    f(depot) <- exp(lfdepot)

    # ----- Observations and residual error
    Cc             <- central             / vc
    Cc_ceftaroline <- central_ceftaroline / vc_ceftaroline
    Celf           <- Cc_ceftaroline * relf

    Cc             ~ add(addSd)             + prop(propSd)
    Cc_ceftaroline ~ add(addSd_ceftaroline) + prop(propSd_ceftaroline)
    Celf           ~ prop(propSd_Celf)
  })
}
