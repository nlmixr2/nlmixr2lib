Franzese_2022_avibactam <- function() {
  description <- paste(
    "Two-compartment IV population PK model for the avibactam component of",
    "ceftazidime-avibactam in children aged 3 months to < 18 years and adults",
    "(Franzese 2022). The adult model of Li 2019 was refitted to 14,223",
    "concentrations from 2,403 subjects after adding 153 children from one",
    "single-dose phase I study and two multiple-dose phase II studies (cIAI",
    "and cUTI). Body weight acts on CL and Q with an exponent of 0.67 and on",
    "Vc and Vp with an exponent of 1, all held constant. For AGE > 2 years,",
    "BSA-normalized creatinine clearance (capped at 150 mL/min/1.73 m^2)",
    "enters CL as a power function below 80 and a linear function above it;",
    "for AGE <= 2 years it is replaced by the Rhodin renal-maturation function",
    "of postmenstrual age (TM50 47.7 weeks, Hill 3.4). End-stage renal",
    "disease, hemodialysis, elevated APACHE II and the adult phase II cIAI",
    "study act on CL; infection type and a ventilator on the PK sampling day",
    "act on Vc. Full 4x4 OMEGA block on CL, Vc, Vp and Q; residual error",
    "switched between a phase I (combined additive + proportional), an adult",
    "phase II and a phase III (proportional) stratum. Companion to",
    "Franzese_2022_ceftazidime; the two analytes were fitted as separate",
    "models."
  )
  reference <- paste(
    "Franzese RC, McFadyen L, Watson KJ, Riccobene T, Carrothers TJ,",
    "Vourvahis M, Chan PLS, Raber S, Bradley JS, Lovern M. Population",
    "Pharmacokinetic Modeling and Probability of Pharmacodynamic Target",
    "Attainment for Ceftazidime-Avibactam in Pediatric Patients Aged 3 Months",
    "and Older. Clin Pharmacol Ther. 2022;111(3):635-645.",
    "doi:10.1002/cpt.2460. Parameter estimates are supplementary Table S3;",
    "the covariate functional forms are the final avibactam NONMEM control",
    "stream printed in the Supplementary Methods (CPT-111-635-s001), whose",
    "$THETA / $OMEGA blocks hold INITIAL values. Erratum: Clin Pharmacol",
    "Ther. 2024;115(2):373, doi:10.1002/cpt.3143 (corrects the transposed",
    "analyte column headings of Table 2 only; no parameter value is",
    "affected).",
    sep = " "
  )
  vignette <- "Franzese_2022_ceftazidime_avibactam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    central = list(analyte = "avibactam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "avibactam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference 70 kg. (WT/70)^0.67 on CL and Q (control stream TH30, shared), (WT/70)^1 on Vc and Vp (TH14, shared). Observed pediatric range 4.1-80.0 kg (Table 1).",
      source_name = "WT"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Used only as a GATE between the CRCL renal function (AGE > 2 years) and the PAGE maturation function (AGE <= 2 years). The age power terms on CL and Vc in the control stream (TH33, TH27) are held at 0.",
      source_name = "AGE"
    ),
    PAGE = list(
      description = "Postmenstrual age",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "WEEKS, matching the Rhodin 2009 constants TM50 = 47.7 weeks and Hill 3.4. Taken as postnatal age + 40 weeks where unknown (Methods). Only read when AGE <= 2 years.",
      source_name = "PMA"
    ),
    CRCL = list(
      description = "BSA-normalized creatinine clearance, NCrCL (Franzese 2022 Methods, 'Model development'): updated bedside Schwartz equation in children, Cockcroft-Gault creatinine clearance x 1.73 / BSA in adults",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Capped at 150 inside model(). For AGE > 2 years and no ESRD: (CRCL/80)^0.986 below 80 and 1 + 0.00344 * (CRCL - 80) at or above 80, so the factor is exactly 1 at the 80 mL/min/1.73 m^2 hinge, which is the reference renal function of the ini() clearance. Factor 1 for ESRD subjects and for AGE <= 2 years.",
      source_name = "CLCRN"
    ),
    RENALIMP_ESRD = list(
      description = "End-stage renal disease indicator",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = not ESRD",
      notes = "Bare multiplier 0.0674 on CL (Table S3 theta5) applied in place of the CRCL factor, for ESRD subjects not on dialysis.",
      source_name = "ESRD"
    ),
    RRT_HEMODIAL_STATUS = list(
      description = "Hemodialysis indicator (control stream DIAL = 1)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = not on dialysis",
      notes = "When 1, the typical clearance theta1 * (study and ESRD factors) is REPLACED by the dialysis clearance 21.1 L/h (Table S3 theta6); the weight, maturation, renal and APACHE factors still multiply it, as in the control stream. Dialysis subjects should also carry RENALIMP_ESRD = 1, which sets the CRCL factor to 1.",
      source_name = "DIAL"
    ),
    APACHE_II_SEV = list(
      description = "Elevated-APACHE-II stratum indicator (control stream APACHE = 2)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = APACHE level 1 or missing (-99); children have no APACHE score and take 0",
      notes = "Proportional shift (1 - 0.192) on CL (Table S3 theta15). Li 2019, the parent analysis, defines the elevated level as APACHE II > 10.",
      source_name = "APACHE"
    ),
    DIS_CIAI = list(
      description = "Complicated intra-abdominal infection indicator (control stream POP = 3)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = not cIAI; the all-zero reference across the DIS_* set is a subject with no cIAI, cUTI or nosocomial-pneumonia classification",
      notes = "Proportional shift (1 + 0.214) on Vc (Table S3 theta12) for phase III adult cIAI and pediatric phase II cIAI subjects. For the adult phase II cIAI study (STUDY_CIAI_PH2 = 1) theta12 is REPLACED by theta9, so model() applies theta12 to DIS_CIAI - STUDY_CIAI_PH2.",
      source_name = "POP = 3"
    ),
    DIS_CUTI = list(
      description = "Complicated urinary tract infection indicator (control stream POP = 2)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = not cUTI",
      notes = "Proportional shift (1 + 0.412) on Vc (Table S3 theta11) for adult and pediatric cUTI subjects in phase II or III studies.",
      source_name = "POP = 2"
    ),
    DIS_HABP = list(
      description = "Hospital-acquired bacterial pneumonia indicator (component of control stream POP = 4)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = not HABP",
      notes = "One nosocomial-pneumonia level covers HAP and VAP; the shared proportional shift (1 + 0.214) on Vc (Table S3 theta12) is applied to DIS_HABP + DIS_VABP.",
      source_name = "POP = 4"
    ),
    DIS_VABP = list(
      description = "Ventilator-associated bacterial pneumonia indicator (component of control stream POP = 4)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = not VABP",
      notes = "See DIS_HABP.",
      source_name = "POP = 4"
    ),
    MECH_VENT = list(
      description = "Ventilator present on the day of PK sampling (control stream POP5 = 1)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = no ventilator on the PK sampling day",
      notes = "Proportional shift (1 + 0.267) on Vc (Table S3 theta28; RSE 55.6%).",
      source_name = "POP5 = 1"
    ),
    STUDY_CIAI_PH2 = list(
      description = "Adult phase 2 cIAI study cohort indicator (control stream POP = 3, PHASE = 2, STDY not 15)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = any other study, including the pediatric phase II cIAI study",
      notes = "Proportional shifts (1 + 0.431) on CL (Table S3 theta10) and (1 + 2.17) on Vc (theta9). The Vc shift REPLACES the cIAI shift theta12 (Li 2019 semantics). The pediatric phase II cIAI study takes theta12 and no CL shift (control stream 'IF(STDY.EQ.15) V1POP=(1 + THETA(12))'; the study-15 CL shift TH37 is held at 0).",
      source_name = "POP = 3 and PHASE = 2"
    ),
    STUDY_CAZAVI_PHASE2 = list(
      description = "Adult phase 2 study stratum of the pooled ceftazidime-avibactam analysis (residual-error switch)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 with STUDY_CAZAVI_PHASE3 also 0 selects the phase I stratum",
      notes = "Avibactam PH2 flag 'IF (STDY.EQ.2001.OR.STDY.EQ.2002) PH2=1', i.e. the two ADULT phase 2 studies only. Pediatric phase I and phase II records fall in the phase I stratum (combined additive + proportional error) in this model; see STUDY_CAZAVI_PED_PHASE2 in the ceftazidime companion.",
      source_name = "STDY = 2001 or 2002"
    ),
    STUDY_CAZAVI_PHASE3 = list(
      description = "Adult phase 3 study stratum of the pooled ceftazidime-avibactam analysis (residual-error switch)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 with STUDY_CAZAVI_PHASE2 also 0 selects the phase I stratum",
      notes = "Avibactam PH3 flag 'IF (PHASE.EQ.3) PH3=1'.",
      source_name = "PHASE = 3"
    )
  )

  # Covariates present in the final control stream with their coefficient
  # held at zero (or one), or screened and not retained.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = male",
      notes = "Control stream CLSEX (TH26) and V1SEX (TH29) are held at 0 in the final model.",
      source_name = "SEX"
    ),
    RACE_ASIAN_OTH = list(
      description = "Non-Chinese, non-Japanese Asian race indicator",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = non-Asian",
      notes = "Control stream CLRCE for RCE = 3 (TH22) is held at 0 in the final model, as are the other race thetas TH21, TH23-TH25. Li 2019 had estimated -0.0865 on CL for this group; the refit removed it.",
      source_name = "RCE = 3"
    ),
    RENAL_ARC = list(
      description = "Augmented renal clearance indicator",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = no augmented renal clearance",
      notes = "Control stream CLCLCR2 = TH13 multiplies the supra-80 CRCL slope for ARC = 1 and is held at 1 in the final model, so ARC has no effect on any CRCL >= 80 subject (Li 2019 had estimated 0.992).",
      source_name = "ARC"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 2403,
    n_observations = 14223,
    n_pediatric_subjects = 153,
    n_pediatric_observations = 488,
    age_range = "Children 0.25-17.7 years (median 7.57; Table 1) pooled with the adult phase I-III population of Li 2019",
    weight_range = "Children 4.1-80.0 kg (median 25.0; Table 1)",
    sex_female_pct = 55.6,
    race_ethnicity = "Children: White 79.1%, Black 3.9%, Chinese 11.8%, other Asian 1.3%, American Indian / Alaska Native 0.7%, other 3.3% (Table 1)",
    disease_state = "Children with suspected or confirmed infection (phase I, n = 32), cIAI (phase II, n = 58) or cUTI (phase II, n = 63), pooled with adults with cIAI, cUTI or nosocomial pneumonia (including VAP), subjects with renal impairment / ESRD / hemodialysis, and healthy volunteers",
    renal_function = "Pediatric baseline NCrCL median 104 (range 43-489) mL/min/1.73 m^2 (Table 1)",
    dose_range = "2-hour IV infusion. Children: 12.5 mg/kg avibactam q8h (maximum 500 mg) for >= 6 months, 10 mg/kg q8h for 3 to < 6 months, halved for CrCL 30 to < 50 mL/min (Table S1). Adults: 500 mg q8h with renal adjustment. Given in a fixed 4:1 ratio with ceftazidime.",
    regions = "Global",
    studies = "Pediatric NCT01893346 (phase I, single dose), NCT02475733 (phase II cIAI), NCT02497781 (phase II cUTI), plus the adult data set of Li 2019",
    unbound_fraction = "0.92 for avibactam (Methods, 'Simulations and PK/PD targets'). The model predicts TOTAL plasma concentration.",
    notes = "Final model estimated with SAEM followed by importance sampling assisted by mode a posteriori (IMPMAP) in NONMEM 7.3. 17 pediatric avibactam concentrations with |CWRES| > 4 were excluded."
  )

  ini({
    # =====================================================================
    # Structural parameters, supplementary Table S3. Reference subject:
    # 70 kg, aged > 2 years, CRCL = 80 mL/min/1.73 m^2 (the hinge, where
    # the renal factor is 1), not ESRD, not on dialysis, APACHE not
    # elevated, no infection-type classification, not ventilated.
    # (The control stream's $THETA holds the Li 2019 initial values, e.g.
    # CL 10.2, and its $OMEGA block is a 0.1 / 0.001 initial guess.)
    # =====================================================================
    lcl <- log(10.7); label("Avibactam clearance at CRCL = 80 mL/min/1.73 m^2 and WT = 70 kg (L/h)") # Table S3: theta1 CL = 10.7 L/h (RSE 3.74%)
    lvc <- log(11.5); label("Avibactam central volume at WT = 70 kg (L)") # Table S3: theta2 Vc = 11.5 L (RSE 4.31%)
    lvp <- log(7.56); label("Avibactam peripheral volume at WT = 70 kg (L)") # Table S3: theta3 Vp = 7.56 L (RSE 14.1%)
    lq <- log(6.94); label("Avibactam intercompartmental clearance at WT = 70 kg (L/h)") # Table S3: theta4 Q = 6.94 L/h (RSE 18.5%)

    lcl_dial <- log(21.1); label("Avibactam clearance in hemodialysis patients at WT = 70 kg (L/h)") # Table S3: theta6 = 21.1 L/h (RSE 9.5%)
    e_renalimp_esrd_cl <- 0.0674; label("Multiplicative factor on avibactam CL for ESRD off dialysis (unitless)") # Table S3: theta5 = 0.0674 (RSE 23.7%)

    # Renal function (AGE > 2 years, not ESRD), CRCL capped at 150.
    e_crcl_cl_lt80 <- 0.986; label("Avibactam CRCL power exponent on CL below 80 mL/min/1.73 m^2 (unitless)") # Table S3: theta7 = 0.986 (RSE 6.34%)
    e_crcl_cl_ge80 <- 0.00344; label("Avibactam CRCL linear slope on CL at or above 80 mL/min/1.73 m^2 (per mL/min/1.73 m^2)") # Table S3: theta8 = 0.00344 (RSE 11.6%)

    # Renal maturation (AGE <= 2 years), Rhodin 2009 (paper reference 26).
    tmat50 <- fixed(47.7); label("Postmenstrual age at 50 percent renal maturation (weeks)") # Results 'Avibactam'; control stream TH31 = 47.7
    hill_mat <- fixed(3.4); label("Hill coefficient of renal maturation (unitless)") # Results 'Avibactam'; control stream TH32 = 3.4

    # Body size (Results 'Avibactam': 0.67 for CL and Q, 1 for volumes).
    e_wt_cl <- fixed(0.67); label("Allometric exponent of body weight on avibactam CL (unitless)") # control stream TH30 = 0.67 (CLWT)
    e_wt_q <- fixed(0.67); label("Allometric exponent of body weight on avibactam Q (unitless)") # control stream TH30 = 0.67 (QWT, shared with CL)
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on avibactam Vc (unitless)") # control stream TH14 = 1 (V1WT)
    e_wt_vp <- fixed(1); label("Allometric exponent of body weight on avibactam Vp (unitless)") # control stream TH14 = 1 (V2WT, shared with Vc)

    # Proportional (1 + theta) covariate shifts.
    e_ciai_ph2_cl <- 0.431; label("Proportional shift in avibactam CL for the adult phase II cIAI study (fraction)") # Table S3: theta10 = 0.431 (RSE 33.4%)
    e_apache_ii_sev_cl <- -0.192; label("Proportional shift in avibactam CL for elevated APACHE II (fraction)") # Table S3: theta15 = -0.192 (RSE 15.4%)
    e_ciai_ph2_vc <- 2.17; label("Proportional shift in avibactam Vc for the adult phase II cIAI study (fraction)") # Table S3: theta9 = 2.17 (RSE 24.8%)
    e_cuti_vc <- 0.412; label("Proportional shift in avibactam Vc for cUTI (fraction)") # Table S3: theta11 = 0.412 (RSE 19.6%)
    e_ciai_habp_vabp_vc <- 0.214; label("Proportional shift in avibactam Vc for phase III cIAI, HAP/VAP or pediatric cIAI (fraction)") # Table S3: theta12 = 0.214 (RSE 26.8%)
    e_mech_vent_vc <- 0.267; label("Proportional shift in avibactam Vc for a ventilator on the PK sampling day (fraction)") # Table S3: theta28 = 0.267 (RSE 55.6%)

    # =====================================================================
    # IIV: full 4x4 block in the control-stream ETA order CL, V1 (Vc),
    # V2 (Vp), Q. Table S3 prints variances and covariances; its BSV column
    # is sqrt(omega^2) (0.3453 -> 58.8, 1.139 -> 107, 1.156 -> 108,
    # 5.487 -> 234) and its printed correlations reproduce from these
    # values: r(Vc,CL) = 0.21, r(Vp,CL) = 0.85, r(Vp,Vc) = -0.29,
    # r(Q,CL) = 0.86, r(Q,Vc) = -0.28, r(Q,Vp) = 0.99.
    # =====================================================================
    etalcl + etalvc + etalvp + etalq ~ c(
      0.3453,
      0.1305, 1.139,
      0.5397, -0.3397, 1.156,
      1.178, -0.7016, 2.495, 5.487
    ) # Table S3 rows etaCL2 through etaQ2

    # =====================================================================
    # Residual error (thetas are SDs; $SIGMA 1). Phase I stratum:
    # sqrt(add^2 + (prop * IPRED)^2); adult phase II and phase III strata:
    # proportional only. Additive SD 43.8 ng/mL = 0.0438 mg/L.
    # =====================================================================
    propSdPhase1 <- 0.174; label("Proportional residual SD, phase I stratum including all pediatric records (fraction)") # Table S3: theta17 = 0.174
    addSdPhase1 <- 0.0438; label("Additive residual SD, phase I stratum including all pediatric records (mg/L)") # Table S3: theta18 = 43.8 ng/mL
    propSdPhase2 <- 0.498; label("Proportional residual SD, adult phase II stratum (fraction)") # Table S3: theta19 = 0.498
    propSdPhase3 <- 0.364; label("Proportional residual SD, phase III stratum (fraction)") # Table S3: theta20 = 0.364
  })

  model({
    # 1. Renal factor on CL. AGE > 2 and not ESRD: power below 80, linear
    #    at or above 80, CRCL capped at 150. ESRD: 1 (theta5 multiplies
    #    the typical CL instead). AGE <= 2: renal maturation.
    crcl_cap <- CRCL * (CRCL <= 150) + 150 * (CRCL > 150)
    older <- (AGE > 2)
    renal_lo <- (crcl_cap < 80) * (1 - RENALIMP_ESRD) * older
    renal_hi <- (crcl_cap >= 80) * (1 - RENALIMP_ESRD) * older
    renal_cl <- (crcl_cap / 80)^e_crcl_cl_lt80 * renal_lo +
      (1 + e_crcl_cl_ge80 * (crcl_cap - 80)) * renal_hi +
      (1 - renal_lo - renal_hi)
    fmat <- PAGE^hill_mat / (tmat50^hill_mat + PAGE^hill_mat)
    mat_cl <- fmat * (1 - older) + older

    # 2. Typical CL before body size: theta1 x study and ESRD factors off
    #    dialysis, or the dialysis clearance theta6 (control stream CLT1).
    clt1 <- exp(lcl) * (1 + e_ciai_ph2_cl * STUDY_CIAI_PH2) *
      (1 + (e_renalimp_esrd_cl - 1) * RENALIMP_ESRD) * (1 - RRT_HEMODIAL_STATUS) +
      exp(lcl_dial) * RRT_HEMODIAL_STATUS

    cl <- clt1 * renal_cl * mat_cl * (WT / 70)^e_wt_cl *
      (1 + e_apache_ii_sev_cl * APACHE_II_SEV) * exp(etalcl)

    # 3. Population effects on Vc. The adult phase II cIAI shift replaces
    #    the cIAI shift (control stream V1POP).
    ciai_other <- DIS_CIAI - STUDY_CIAI_PH2
    pop_vc <- 1 +
      e_cuti_vc * DIS_CUTI +
      e_ciai_ph2_vc * STUDY_CIAI_PH2 +
      e_ciai_habp_vabp_vc * (ciai_other + DIS_HABP + DIS_VABP)

    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc * pop_vc *
      (1 + e_mech_vent_vc * MECH_VENT)
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vp
    q <- exp(lq + etalq) * (WT / 70)^e_wt_q

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 4. Total plasma concentration (mg/L). Multiply by 0.92 for the free
    #    concentration used in the 50% fT > 1 mg/L target.
    Cc <- central / vc

    phase1 <- 1 - STUDY_CAZAVI_PHASE2 - STUDY_CAZAVI_PHASE3
    propSd <- propSdPhase1 * phase1 +
      propSdPhase2 * STUDY_CAZAVI_PHASE2 +
      propSdPhase3 * STUDY_CAZAVI_PHASE3
    addSd <- addSdPhase1 * phase1

    Cc ~ add(addSd) + prop(propSd)
  })
}
