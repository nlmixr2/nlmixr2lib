Franzese_2022_ceftazidime <- function() {
  description <- paste(
    "Two-compartment IV population PK model for the ceftazidime component of",
    "ceftazidime-avibactam in children aged 3 months to < 18 years and adults",
    "(Franzese 2022). The adult model of Li 2019 was refitted to 9,628",
    "concentrations from 2,130 subjects after adding 153 children from one",
    "single-dose phase I study and two multiple-dose phase II studies (cIAI",
    "and cUTI). Body weight acts on CL through a normalized Emax function",
    "(half-maximal at 53.5 kg, equal to 1 at 70 kg) and on Q, Vc and Vp",
    "through power terms with exponents of 0.67, 1 and 1 held constant. For",
    "AGE > 2 years, BSA-normalized creatinine clearance (bedside Schwartz in",
    "children, Cockcroft-Gault x 1.73 / BSA in adults, capped at 150",
    "mL/min/1.73 m^2) enters CL as a two-segment linear spline hinged at 100;",
    "for AGE <= 2 years it is replaced by the Rhodin renal-maturation function",
    "of postmenstrual age (TM50 47.7 weeks, Hill 3.4). Infection type (cIAI,",
    "nosocomial pneumonia), Asian race and a ventilator on the PK sampling",
    "day act on CL and/or Vc. Diagonal IIV on CL, Vc, Vp and Q; combined",
    "additive + proportional residual error switched between a phase I and a",
    "pooled phase II/III stratum. Companion to Franzese_2022_avibactam; the",
    "two analytes were fitted as separate models."
  )
  reference <- paste(
    "Franzese RC, McFadyen L, Watson KJ, Riccobene T, Carrothers TJ,",
    "Vourvahis M, Chan PLS, Raber S, Bradley JS, Lovern M. Population",
    "Pharmacokinetic Modeling and Probability of Pharmacodynamic Target",
    "Attainment for Ceftazidime-Avibactam in Pediatric Patients Aged 3 Months",
    "and Older. Clin Pharmacol Ther. 2022;111(3):635-645.",
    "doi:10.1002/cpt.2460. Parameter estimates are supplementary Table S2;",
    "the covariate functional forms are the final ceftazidime NONMEM control",
    "stream printed in the Supplementary Methods (CPT-111-635-s001). Erratum:",
    "Clin Pharmacol Ther. 2024;115(2):373, doi:10.1002/cpt.3143 (corrects",
    "the transposed analyte column headings of Table 2 only; no parameter",
    "value is affected).",
    sep = " "
  )
  vignette <- "Franzese_2022_ceftazidime_avibactam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    central = list(analyte = "ceftazidime", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ceftazidime", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference 70 kg. Acts on CL through the normalized Emax function WT / (53.5 + WT) * (53.5 + 70) / 70 (Supplementary Methods, 'Ceftazidime Emax function for covariate effect of body weight'), whose ceiling is 1 + 53.5/70 = 1.764 for very heavy subjects (Table S2 footnote). Acts on Q as (WT/70)^0.67 and on Vc and Vp as (WT/70)^1. Observed pediatric range 4.1-80.0 kg (Table 1).",
      source_name = "WT"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Used only as a GATE: for AGE > 2 years the renal factor on CL is the CRCL spline, and for AGE <= 2 years it is the PAGE maturation function instead (control stream 'IF (AGE.LE.2) CLPMA = ...' and 'IF (NCLCR.LT.INFLECT.AND.AGE.GT.2) ...'). Age itself carries no coefficient.",
      source_name = "AGE"
    ),
    PAGE = list(
      description = "Postmenstrual age",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "WEEKS, not the register default of months, because the Rhodin 2009 maturation constants (TM50 = 47.7 weeks, Hill 3.4) are on that scale (same convention as Germovsek_2018_meropenem.R and Riccobene_2017_ceftaroline.R). Franzese 2022 Methods: where postmenstrual age was unknown it was taken as postnatal age + 40 weeks. Only read when AGE <= 2 years; supply any positive value for older subjects.",
      source_name = "PMA"
    ),
    CRCL = list(
      description = "BSA-normalized creatinine clearance, NCrCL (Franzese 2022 Methods, 'Model development'): updated bedside Schwartz equation in children, Cockcroft-Gault creatinine clearance x 1.73 / BSA in adults",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Capped at 150 mL/min/1.73 m^2 inside model() (control stream 'IF (NCLCR.GT.150) NCLCR=150'). Enters CL as slope1 * CRCL below 100 and slope1 * 100 + slope2 * (CRCL - 100) at or above 100, with both slopes carried from an earlier model iteration (Table S2 footnote: 'obtained from a previous model iteration (MS-06)'). The lower segment passes through the origin, so the factor equals 1 at CRCL = 1/0.0103036 = 97.1, which is the reference renal function of the ini() clearance. Only read when AGE > 2 years. Observed pediatric median 104 (range 43-489) (Table 1). This is a change from the parent Li 2019 model, which used raw Cockcroft-Gault mL/min.",
      source_name = "CLCRN"
    ),
    DIS_CIAI = list(
      description = "Complicated intra-abdominal infection indicator (control stream POP = 3)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = not cIAI; the all-zero reference across the DIS_* set is a subject with no cIAI, cUTI or nosocomial-pneumonia classification (the phase I population)",
      notes = "Bare multiplier 1.33 on CL (Table S2 theta16) and 1.83 on Vc, shared with nosocomial pneumonia (theta21). Applies to adult phase II and III and to pediatric phase II cIAI subjects alike. Mutually exclusive with DIS_CUTI, DIS_HABP and DIS_VABP.",
      source_name = "POP = 3"
    ),
    DIS_CUTI = list(
      description = "Complicated urinary tract infection indicator (control stream POP = 2)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = not cUTI",
      notes = "Bare multiplier 1.49 on Vc (Table S2 theta20); no effect on CL. Applies to adult and pediatric cUTI subjects alike.",
      source_name = "POP = 2"
    ),
    DIS_HABP = list(
      description = "Hospital-acquired bacterial pneumonia indicator (component of control stream POP = 4, nosocomial pneumonia)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = not HABP",
      notes = "Franzese 2022 carries a single nosocomial-pneumonia level (POP = 4) covering HAP and VAP, so DIS_HABP and DIS_VABP share one coefficient and enter model() as their sum: bare multiplier 1.10 on CL (Table S2 theta17) and 1.83 on Vc (theta21, shared with cIAI). Kept as two columns for consistency with the register and with Li_2019_ceftazidime.R.",
      source_name = "POP = 4"
    ),
    DIS_VABP = list(
      description = "Ventilator-associated bacterial pneumonia indicator (component of control stream POP = 4, nosocomial pneumonia)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = not VABP",
      notes = "See DIS_HABP; one pooled coefficient is applied to DIS_HABP + DIS_VABP. No pediatric HAP/VAP subject was studied; the paper simulates pediatric HAP/VAP by carrying the adult covariate effects over.",
      source_name = "POP = 4"
    ),
    MECH_VENT = list(
      description = "Ventilator present on the day of PK sampling (control stream POP5 = 1)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = no ventilator on the PK sampling day",
      notes = "Proportional shift (1 + 0.202) on Vc (Table S2 theta22, 'Population effect on Vc for presence of ventilator'). Same NPv definition as Li 2019: time-fixed at the sampling day and not the same contrast as VABP vs HABP.",
      source_name = "POP5 = 1"
    ),
    RACE_ASIAN_OTH = list(
      description = "Non-Chinese, non-Japanese Asian race indicator (control stream RCE = 3, ASN)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = non-Asian (White, Black and other; the dominant grouping)",
      notes = "Proportional shift (1 - 0.136) on CL (Table S2 theta18) plus the pooled Asian shift (1 - 0.135) on Vc (theta23). Mutually exclusive with RACE_CHINESE and RACE_JAPANESE.",
      source_name = "RCE = 3"
    ),
    RACE_CHINESE = list(
      description = "Chinese (including Taiwanese) race indicator (control stream RCE = 14, CHN)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = not Chinese",
      notes = "Proportional shift (1 - 0.0844) on CL (Table S2 theta19) plus the pooled Asian shift (1 - 0.135) on Vc (theta23). 18 of the 153 children were Chinese (Table 1).",
      source_name = "RCE = 14"
    ),
    RACE_JAPANESE = list(
      description = "Japanese race indicator (control stream RCE = 13, JPN)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = not Japanese",
      notes = "Pooled Asian shift (1 - 0.135) on Vc only (Table S2 theta23); the control stream's CLRCE block assigns no clearance factor to RCE = 13.",
      source_name = "RCE = 13"
    ),
    STUDY_CAZAVI_PHASE2 = list(
      description = "Adult phase 2 study stratum of the pooled ceftazidime-avibactam analysis (residual-error switch)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 with STUDY_CAZAVI_PHASE3 and STUDY_CAZAVI_PED_PHASE2 also 0 selects the phase I stratum",
      notes = "The ceftazidime residual error pools every PHASE = 2 or 3 record into one stratum (control stream $ERROR 'IF(PHASE.EQ.2.OR.PHASE.EQ.3)'), so this model uses only the sum of the three study indicators.",
      source_name = "PHASE = 2"
    ),
    STUDY_CAZAVI_PHASE3 = list(
      description = "Adult phase 3 study stratum of the pooled ceftazidime-avibactam analysis (residual-error switch)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 with the two phase 2 indicators also 0 selects the phase I stratum",
      notes = "See STUDY_CAZAVI_PHASE2.",
      source_name = "PHASE = 3"
    ),
    STUDY_CAZAVI_PED_PHASE2 = list(
      description = "Pediatric phase 2 study indicator (cIAI NCT02475733 or cUTI NCT02497781) of the Franzese 2022 pooled ceftazidime-avibactam analysis (residual-error switch)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = any other study, including the single-dose pediatric phase I study NCT01893346",
      notes = "Selects the pooled phase II/III residual stratum for ceftazidime, exactly as an adult phase 2 record does. Kept apart from STUDY_CAZAVI_PHASE2 because the avibactam companion model places these same pediatric records in its phase I residual stratum (its PH2 flag tests the adult study numbers STDY 2001/2002, not PHASE), so one covariate set can drive both models.",
      source_name = "PHASE = 2 with STDY = 15 or 16"
    )
  )

  # Covariates whose coefficients are held at zero in the final control
  # stream, or that were screened and not retained.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = male",
      notes = "Re-evaluated in forward inclusion / backward elimination; not retained in the final ceftazidime model.",
      source_name = "SEX"
    ),
    APACHE_II_SEV = list(
      description = "Elevated-APACHE-II stratum indicator",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = APACHE II not elevated",
      notes = "Not in the final ceftazidime model (as in Li 2019); the avibactam companion retains it on CL.",
      source_name = "APACHE"
    ),
    RENALIMP_ESRD = list(
      description = "End-stage renal disease indicator",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = not ESRD",
      notes = "The ceftazidime model has no ESRD or dialysis term; the CRCL spline alone spans the renal range.",
      source_name = "ESRD"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 2130,
    n_observations = 9628,
    n_pediatric_subjects = 153,
    n_pediatric_observations = 509,
    age_range = "Children 0.25-17.7 years (median 7.57; Table 1) pooled with the adult phase I-III population of Li 2019",
    weight_range = "Children 4.1-80.0 kg (median 25.0; Table 1)",
    sex_female_pct = 55.6,
    race_ethnicity = "Children: White 79.1%, Black 3.9%, Chinese 11.8%, other Asian 1.3%, American Indian / Alaska Native 0.7%, other 3.3% (Table 1)",
    disease_state = "Children with suspected or confirmed infection (phase I, n = 32), cIAI (phase II, n = 58) or cUTI (phase II, n = 63), pooled with adults with cIAI, cUTI or nosocomial pneumonia (including VAP) and healthy volunteers",
    renal_function = "Pediatric baseline NCrCL median 104 (range 43-489) mL/min/1.73 m^2 (Table 1); one pediatric cUTI subject had NCrCL < 50",
    dose_range = "2-hour IV infusion. Children: 50-12.5 mg/kg ceftazidime-avibactam q8h (maximum 2,000-500 mg) for >= 6 months, 40-10 mg/kg q8h for 3 to < 6 months, halved for CrCL 30 to < 50 mL/min (Table S1). Adults: 2,000-500 mg q8h with renal adjustment.",
    regions = "Global",
    studies = "Pediatric NCT01893346 (phase I, single dose), NCT02475733 (phase II cIAI), NCT02497781 (phase II cUTI), plus the adult data set of Li 2019 (11 phase I, 2 phase II and 5 phase III studies)",
    unbound_fraction = "0.85 for ceftazidime (Methods, 'Simulations and PK/PD targets'). The model predicts TOTAL plasma concentration.",
    notes = "FOCE-INTER in NONMEM 7.3. Outliers with |CWRES| > 4 (30 adult and 3 pediatric ceftazidime concentrations) were excluded from the final model."
  )

  ini({
    # =====================================================================
    # Structural parameters, supplementary Table S2. They refer to a 70 kg
    # subject with CRCL = 97.1 mL/min/1.73 m^2 (where the renal spline
    # equals 1), aged > 2 years, with no infection-type classification
    # (all DIS_* = 0), non-Asian and not ventilated.
    #
    # The control stream's $THETA block carries 9.1259 for CL; every other
    # $THETA matches Table S2 to rounding. The table value 7.75 is the final
    # estimate (Table S2 footnote: it is run121, whose $PROBLEM line says the
    # weight function was re-parameterized relative to run120, the source of
    # the initial values). See the vignette for the check against Table 3.
    # =====================================================================
    lcl <- log(7.75); label("Ceftazidime clearance at CRCL = 97.1 mL/min/1.73 m^2 and WT = 70 kg (L/h)") # Table S2: theta1 CL = 7.75 L/h (RSE 1.56%)
    lvc <- log(11.2); label("Ceftazidime central volume at WT = 70 kg (L)") # Table S2: theta2 Vc = 11.2 L (RSE 3.54%)
    lq <- log(5.33); label("Ceftazidime intercompartmental clearance at WT = 70 kg (L/h)") # Table S2: theta3 Q = 5.33 L/h (RSE 6.52%)
    lvp <- log(6.52); label("Ceftazidime peripheral volume at WT = 70 kg (L)") # Table S2: theta4 Vp = 6.52 L (RSE 3.12%)

    # =====================================================================
    # Renal function (AGE > 2 years): two-segment linear spline in CRCL,
    # hinged at 100 and capped at 150. Both slopes are carried unestimated
    # from an earlier model iteration (Table S2 rows 'Slope 1' / 'Slope 2'
    # and footnote; control stream constants SLOPE1 / SLOPE2).
    # =====================================================================
    e_crcl_cl_lt100 <- fixed(0.0103036); label("Ceftazidime CRCL linear slope on CL below 100 mL/min/1.73 m^2 (per mL/min/1.73 m^2)") # Table S2 'Slope 1' = 0.01030360; control stream SLOPE1
    e_crcl_cl_ge100 <- fixed(0.00125182); label("Ceftazidime CRCL linear slope on CL at or above 100 mL/min/1.73 m^2 (per mL/min/1.73 m^2)") # Table S2 'Slope 2' = 0.00125182; control stream SLOPE2

    # =====================================================================
    # Renal maturation (AGE <= 2 years), Rhodin 2009 (paper reference 26):
    # PAGE^hill / (tmat50^hill + PAGE^hill), PAGE in weeks.
    # =====================================================================
    tmat50 <- fixed(47.7); label("Postmenstrual age at 50 percent renal maturation (weeks)") # Methods 'Model development'; control stream TH7 = 4.77E+01, Rhodin 2009
    hill_mat <- fixed(3.4); label("Hill coefficient of renal maturation (unitless)") # Methods 'Model development'; control stream TH6 = 3.4E+00, Rhodin 2009

    # =====================================================================
    # Body size. CL: normalized Emax in weight with zero intercept
    # (control stream TH13 = 0) and Hill coefficient 1 (Table S2 theta14),
    #   f(WT) = WT^h / (WT50^h + WT^h) * (WT50^h + 70^h) / 70^h,
    # which equals 1 at 70 kg (Supplementary Methods). Q, Vc and Vp: power
    # terms with exponents held constant.
    # =====================================================================
    e_wt50_cl <- 53.5; label("Body weight at the half-maximal weight effect on ceftazidime CL (kg)") # Table S2: theta15 = 53.5 kg (RSE 8.81%)
    e_wt_cl_hill <- fixed(1); label("Hill coefficient of the Emax weight effect on ceftazidime CL (unitless)") # Table S2: theta14 = 1.0, not estimated in run121
    e_wt_q <- fixed(0.67); label("Allometric exponent of body weight on ceftazidime Q (unitless)") # Results 'Ceftazidime'; control stream TH12 = 0.67
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on ceftazidime Vc (unitless)") # Results 'Ceftazidime'; control stream TH5 = 1.00
    e_wt_vp <- fixed(1); label("Allometric exponent of body weight on ceftazidime Vp (unitless)") # Results 'Ceftazidime'; control stream TH5 = 1.00 (shared with Vc)

    # =====================================================================
    # Infection-type effects. The control stream writes these as BARE
    # multipliers ('IF(POP.EQ.3) CLPOP = THETA(16)'), so 1 means no effect.
    # =====================================================================
    e_ciai_cl <- 1.33; label("Multiplicative factor on ceftazidime CL for cIAI (unitless)") # Table S2: theta16 = 1.33 (RSE 2.37%)
    e_habp_vabp_cl <- 1.1; label("Multiplicative factor on ceftazidime CL for HAP/VAP (unitless)") # Table S2: theta17 = 1.1 (RSE 2.96%)
    e_cuti_vc <- 1.49; label("Multiplicative factor on ceftazidime Vc for cUTI (unitless)") # Table S2: theta20 = 1.49 (RSE 4.57%)
    e_ciai_habp_vabp_vc <- 1.83; label("Multiplicative factor on ceftazidime Vc for cIAI or HAP/VAP (unitless)") # Table S2: theta21 = 1.83 (RSE 3.97%)

    # Proportional (1 + theta) shifts.
    e_mech_vent_vc <- 0.202; label("Proportional shift in ceftazidime Vc for a ventilator on the PK sampling day (fraction)") # Table S2: theta22 = 0.202 (RSE 33.5%)
    e_race_asian_oth_cl <- -0.136; label("Proportional shift in ceftazidime CL for non-Chinese, non-Japanese Asian race (fraction)") # Table S2: theta18 = -0.136 (RSE 20.3%)
    e_race_chinese_cl <- -0.0844; label("Proportional shift in ceftazidime CL for Chinese race (fraction)") # Table S2: theta19 = -0.0844 (RSE 29.1%)
    e_race_asian_vc <- -0.135; label("Proportional shift in ceftazidime Vc for any Asian race (fraction)") # Table S2: theta23 = -0.135 (RSE 23.2%)

    # =====================================================================
    # IIV: diagonal OMEGA, log-scale variances (control stream comment
    # 'They are VARIANCES!'). Table S2 BSV column = sqrt(exp(omega^2) - 1):
    # 0.154 -> 40.8%, 0.108 -> 33.8%, 0.203 -> 47.5%, 0.0236 -> 15.4%.
    # =====================================================================
    etalcl ~ 0.154 # Table S2 row etaCL2 = 0.154 (BSV 40.8 percent CV)
    etalvc ~ 0.108 # Table S2 row etaVc2 = 0.108 (BSV 33.8 percent CV)
    etalq ~ 0.203 # Table S2 row etaQ2 = 0.203 (BSV 47.5 percent CV)
    etalvp ~ 0.0236 # Table S2 row etaVp2 = 0.0236 (BSV 15.4 percent CV)

    # =====================================================================
    # Residual error: W = sqrt(add^2 + (prop * IPRED)^2) with $SIGMA 1, so
    # the thetas are standard deviations (control stream comment 'THETAs
    # for residuals are Standard Deviations'). Additive terms are in ng/mL
    # in the source and are divided by 1000 here for mg/L.
    # =====================================================================
    propSdPhase1 <- 0.172; label("Proportional residual SD, phase I stratum (fraction)") # Table S2: theta8 = 0.172
    addSdPhase1 <- 0.125; label("Additive residual SD, phase I stratum (mg/L)") # Table S2: theta9 = 125 ng/mL
    propSdPhase23 <- 0.374; label("Proportional residual SD, phase II/III stratum (fraction)") # Table S2: theta10 = 0.374
    addSdPhase23 <- 2.56; label("Additive residual SD, phase II/III stratum (mg/L)") # Table S2: theta11 = 2560 ng/mL
  })

  model({
    # 1. Renal factor on CL, gated by age. AGE > 2: CRCL spline capped at
    #    150. AGE <= 2: renal maturation in postmenstrual age.
    crcl_cap <- CRCL * (CRCL <= 150) + 150 * (CRCL > 150)
    renal_crcl <- e_crcl_cl_lt100 * crcl_cap * (crcl_cap < 100) +
      (e_crcl_cl_lt100 * 100 + e_crcl_cl_ge100 * (crcl_cap - 100)) * (crcl_cap >= 100)
    fmat <- PAGE^hill_mat / (tmat50^hill_mat + PAGE^hill_mat)
    older <- (AGE > 2)
    renal_cl <- renal_crcl * older + fmat * (1 - older)

    # 2. Normalized Emax weight effect on CL (1 at 70 kg).
    wt_cl <- WT^e_wt_cl_hill / (e_wt50_cl^e_wt_cl_hill + WT^e_wt_cl_hill) *
      (e_wt50_cl^e_wt_cl_hill + 70^e_wt_cl_hill) / 70^e_wt_cl_hill

    # 3. Infection-type factors (mutually exclusive indicators, bare
    #    multipliers).
    infect_cl <- (1 - DIS_CIAI - DIS_HABP - DIS_VABP) +
      e_ciai_cl * DIS_CIAI +
      e_habp_vabp_cl * (DIS_HABP + DIS_VABP)
    infect_vc <- (1 - DIS_CUTI - DIS_CIAI - DIS_HABP - DIS_VABP) +
      e_cuti_vc * DIS_CUTI +
      e_ciai_habp_vabp_vc * (DIS_CIAI + DIS_HABP + DIS_VABP)

    # 4. Individual parameters (control stream $PK).
    cl <- exp(lcl + etalcl) * wt_cl * renal_cl * infect_cl *
      (1 + e_race_asian_oth_cl * RACE_ASIAN_OTH + e_race_chinese_cl * RACE_CHINESE)
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc * infect_vc *
      (1 + e_mech_vent_vc * MECH_VENT) *
      (1 + e_race_asian_vc * (RACE_ASIAN_OTH + RACE_CHINESE + RACE_JAPANESE))
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vp
    q <- exp(lq + etalq) * (WT / 70)^e_wt_q

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 5. Total plasma concentration (mg/L). Multiply by 0.85 for the free
    #    concentration used in the 50% fT > MIC target.
    Cc <- central / vc

    phase23 <- STUDY_CAZAVI_PHASE2 + STUDY_CAZAVI_PHASE3 + STUDY_CAZAVI_PED_PHASE2
    propSd <- propSdPhase1 * (1 - phase23) + propSdPhase23 * phase23
    addSd <- addSdPhase1 * (1 - phase23) + addSdPhase23 * phase23

    Cc ~ add(addSd) + prop(propSd)
  })
}
