Jamois_2025_tocilizumab_crs <- function() {
  description <- "Joint population PK / soluble IL-6 receptor (sIL-6R) model for intravenous tocilizumab given to manage cytokine release syndrome (CRS) induced by the anti-CD20 T-cell-engaging bispecific antibodies mosunetuzumab or glofitamab in adults with relapsed/refractory B-cell non-Hodgkin lymphoma (Jamois 2025). Two-compartment tocilizumab disposition with first-order SC absorption and parallel linear and Michaelis-Menten elimination; total sIL-6R described by the quasi-steady-state approximation of target-mediated drug disposition (an indirect-response model with inhibition of elimination). The rheumatoid-arthritis model is held fixed and CRS factors are estimated on CL (first CRS event grade >= 2 only), VC/VP, Q, kint and kdeg."
  reference <- "Jamois C, Turner DC, Gibiansky L, Li F, Menuet JF, Pereira LR, Lundberg L, Guerini E, Carlile DJ, Relf J, Wei MC, Kwan A, Frey N, Grange S. Tocilizumab Dosing for Management of T Cell-Engaging Bispecific Antibody-Related CRS in Patients With R/R B-Cell NHL. Clin Pharmacol Ther. 2025;118(4):917. doi:10.1002/cpt.3751"
  vignette <- "Jamois_2025_tocilizumab"
  paper_specific_etas <- c("etaruv", "etaruv_total_target")
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    depot = list(analyte = "tocilizumab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tocilizumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "tocilizumab", units = "mg", specimen = "tissue", verified = TRUE),
    total_target = list(
      analyte = "total (unbound + tocilizumab-bound) soluble IL-6 receptor",
      units = "ng/mL",
      specimen = "serum",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power scaling on CL and Q (shared exponent 0.502), on VC and VP (shared exponent 0.750) and on kdeg (-0.0781), each as (WT/70)^theta (Table S6: LOGWT = LOG(WT/70)).",
      source_name = "WT"
    ),
    HDLC = list(
      description = "High-density lipoprotein cholesterol",
      units = "mmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect (HDLC/1.5)^-0.252 on linear CL (Table S6: LOGHDL = LOG(HDL/1.5)); carried fixed from the rheumatoid-arthritis model. Table S4 reports HDL in mmol/L (median 1.5).",
      source_name = "HDL"
    ),
    CRCL = list(
      description = "BSA-normalized creatinine clearance",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect (CRCL/100)^0.357 on the Michaelis-Menten Vmax, with CRCL truncated to [15, 150] before use (Table S6: CRCLN1 bounds, LOGCRCL = LOG(CRCLN1/100)). Source column CRCLN ('calculated creatinine clearance', Table S4, mL/min/1.73m2); the estimating equation is not stated.",
      source_name = "CRCLN"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Three roles, all carried fixed from the rheumatoid-arthritis / juvenile idiopathic arthritis model: (1) CL maturation factor 1 - 0.262*exp(-log(2)*AGE/2.44); (2) power effect (AGE/50)^-0.354 on ka; (3) hyperbolic AGE/(12.0 + AGE) scaling of log SC bioavailability. In the adult CRS population (age 30-87 years) the CL maturation factor is 1 to within 0.01%.",
      source_name = "AGE"
    ),
    CRS_GRADE_GE2 = list(
      description = "First CRS event of ASTCT grade >= 2 indicator (1 = grade >= 2, 0 = grade 1)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (first CRS event grade 1)",
      notes = "Table S6: FIRSTCRS1 = FIRSTCRS (grade of the first CRS event), set to 2 when missing; CRS2 = 1 if FIRSTCRS1 >= 2. Multiplies linear CL by exp(theta31) = 1.60 (Table S5 'CL CRS grade>1'). Supply 1 when the grade is unknown, matching the source imputation.",
      source_name = "FIRSTCRS"
    ),
    TSIL6R_BASE = list(
      description = "Observed baseline (pre-tocilizumab) serum total soluble IL-6 receptor concentration",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Table S6: BASE = observed SIL6R when available, otherwise the model-predicted typical baseline exp(theta24 + eta9) (sBASE = 29.4 ng/mL, Table S5). Missing values are coded negative in the source data set (SIL6RMISS = 1 if SIL6R < 0); supply 0 or a negative value to use the model baseline. BASE is the initial condition of total_target and scales kint and kdeg as (BASE/40)^theta.",
      source_name = "SIL6R"
    ),
    STUDY_25220_28117 = list(
      description = "Indicator for the tocilizumab SC studies NA25220 and 28117 of the rheumatoid-arthritis / JIA model",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (all other studies, including every CRS study)",
      notes = "Table S6: ST2 = 1 if STUD = 25220 or 28117; multiplies ka by exp(theta16) = 0.807 (Table S1 'ka, NA25220'). Inherited fixed from the rheumatoid-arthritis model; always 0 in the CRS population (IV dosing).",
      source_name = "STUD"
    ),
    SC_INJ_GROUP3 = list(
      description = "Indicator for subcutaneous-injection category SINJ = 3 of the rheumatoid-arthritis / JIA data set (coding not defined by the source)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (any other SINJ category, and all IV dosing)",
      notes = "Table S6: FF1 = 1; IF(SINJ.EQ.3) FF1 = THETA(17); LTPF1 = THETA(8)*AGE/(THETA(23)+AGE)*FF1. Multiplies the log SC bioavailability by theta17 = 0.692 (Table S1 'F SC,SJIT3'), i.e. raises F toward 1. Neither the paper nor the supplement defines SINJ or its level 3. Inherited fixed from the rheumatoid-arthritis model; it acts only on SC doses, so it is irrelevant for the all-IV CRS population (supply 0).",
      source_name = "SINJ"
    )
  )

  covariatesDataExcluded <- list(
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Rheumatoid-arthritis model effect LOG(ALB/35)*theta12 on VC and VP (Table S1: -0.526) was fixed to 0 in the CRS model (Table S5 'V ALB = 0 Fixed'; Table S6 $THETA 12 = 0 FIX)."
    ),
    TPRO = list(
      description = "Baseline total serum protein",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Rheumatoid-arthritis model effect LOG(PROT/75)*theta13 on VC and VP (Table S1: 0.562) was fixed to 0 in the CRS model (Table S5 'V PROT = 0 Fixed'; Table S6 $THETA 13 = 0 FIX)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 61L,
    n_observations = 506L,
    n_studies = 2L,
    age_range = "30-87 years (mean 63.0, median 64)",
    weight_range = "53-128 kg (mean 75.7, median 73.8)",
    sex_female_pct = 27.9,
    disease_state = "Relapsed/refractory B-cell non-Hodgkin lymphoma with cytokine release syndrome (ASTCT grade >= 1) after mosunetuzumab (NCT02500407, n = 16) or glofitamab (NCT03075696, n = 45), receiving tocilizumab for CRS management.",
    dose_range = "Tocilizumab IV 8 mg/kg (12 mg/kg if < 30 kg; no more than 800 mg per infusion); repeat doses permitted, up to four doses every 8 h in the external-validation studies.",
    regions = "Global (Roche-sponsored phase 1/2 studies)",
    baseline_sil6r = "Total sIL-6R median 33.4 ng/mL (range 20.3-103) (Table S4)",
    notes = "Model-development data set (Table 2 and Table S4): 61 patients, 276 tocilizumab and 230 sIL-6R samples. External validation used 6 further mosunetuzumab-treated patients from NCT03677141, NCT03671018, NCT03677154 and NCT03075696 (67 patients, 543 samples in total). Assay LLOQ 0.1 ug/mL for tocilizumab."
  )

  ini({
    # Rheumatoid-arthritis / JIA model parameters held fixed in the CRS fit
    # (Table S6 $THETA FIX; values match Table S1 to its printed precision).
    lcl <- fixed(-1.556766391); label("Linear clearance CL at the reference patient, RA model (L/day)") # Table S6 THETA(1) FIX; Table S1 CL = exp(theta1) = 0.211 L/day
    lvc <- fixed(1.213239670); label("Central volume VC, RA model (L)") # Table S6 THETA(2) FIX; Table S1 VC = exp(theta2) = 3.36 L
    lvp <- fixed(1.190255642); label("Peripheral volume VP, RA model (L)") # Table S6 THETA(3) FIX; Table S1 VP = exp(theta3) = 3.29 L
    lq <- fixed(-0.937571287); label("Intercompartmental clearance Q, RA model (L/day)") # Table S6 THETA(4) FIX; Table S1 Q = exp(theta4) = 0.392 L/day
    lvmax <- fixed(0.836626291); label("Michaelis-Menten maximum elimination rate VM (mg/L/day)") # Table S6 THETA(5) FIX; Table S1 VM = exp(theta5) = 2.31 mg/L/day
    lkm <- fixed(-0.882927120); label("Michaelis-Menten constant KM (ug/mL)") # Table S6 THETA(6) FIX; Table S1 KM = exp(theta6) = 0.414 mcg/mL
    lka <- fixed(-1.564654589); label("First-order SC absorption rate ka (1/day)") # Table S6 THETA(7) FIX; Table S1 ka = exp(theta7) = 0.209 1/day
    lfdepot <- fixed(-0.342109352); label("Log SC bioavailability at the adult asymptote (unitless)") # Table S6 THETA(8) FIX; Table S1 FSC = exp(theta8) = 0.710
    e_wt_cl <- fixed(0.502061606); label("Power exponent of WT/70 on CL and Q (unitless)") # Table S6 THETA(9) FIX; Table S1 'CL WT, Q WT' = 0.502
    e_wt_vc <- fixed(0.749708593); label("Power exponent of WT/70 on VC and VP (unitless)") # Table S6 THETA(10) FIX; Table S1 'VC,WT, VP,WT' = 0.750
    e_hdlc_cl <- fixed(-0.252391011); label("Power exponent of HDLC/1.5 on CL (unitless)") # Table S6 THETA(11) FIX; Table S1 'CL HDL' = -0.252
    e_crcl_vmax <- fixed(0.356767267); label("Power exponent of CRCL/100 on VM (unitless)") # Table S6 THETA(14) FIX; Table S1 'VM,CRCLN' = 0.357
    e_age_ka <- fixed(-0.354390651); label("Power exponent of AGE/50 on ka (unitless)") # Table S6 THETA(15) FIX; Table S1 'ka,AGE' = -0.354
    e_study_25220_28117_ka <- fixed(-0.214487419); label("Log factor on ka for studies 25220/28117 (unitless)") # Table S6 THETA(16) FIX; Table S1 'ka, NA25220' = exp(theta16) = 0.807
    e_sc_inj_group3_lfdepot <- fixed(0.691991210); label("Multiplicative factor on log SC bioavailability for SINJ = 3 (unitless)") # Table S6 THETA(17) FIX; Table S1 'F SC,SJIT3' = 0.692
    sdL <- fixed(0.998582423); label("Log-scale residual SD at low tocilizumab concentrations SDL (unitless)") # Table S6 THETA(18) FIX; Table S1 'sigma SDL' = 0.999
    sdH <- fixed(0.116924070); label("Log-scale residual SD at high tocilizumab concentrations SDH (unitless)") # Table S6 THETA(19) FIX; Table S1 'sigma SDH' = 0.117
    sd50 <- fixed(2.189973116); label("Tocilizumab concentration at which the residual SD is midway between SDL and SDH (ug/mL)") # Table S6 THETA(20) FIX; Table S1 'sigma SD50' = 2.19
    e_age_cl <- fixed(0.261574358); label("Maximum fractional reduction of CL at birth (unitless)") # Table S6 THETA(21) FIX; Table S1 'CL maturation' = 0.262
    thalf_age_cl <- fixed(2.440289736); label("Half-life of CL maturation with age (years)") # Table S6 THETA(22) FIX; Table S1 't1/2,matur' = 2.44 years
    age50_fdepot <- fixed(12.03537464); label("Age at half the adult log-bioavailability reduction (years)") # Table S6 THETA(23) FIX; Table S1 'FSC,Age50' = 12.0 years
    lkint <- fixed(-2.615842819); label("Elimination rate of tocilizumab-bound sIL-6R kint at BASE = 40 ng/mL, RA model (1/day)") # Table S6 THETA(25) FIX; Table S1 kint = exp(theta25) = 0.0731
    lkdeg <- fixed(0.290626168); label("Degradation rate of unbound sIL-6R kdeg at BASE = 40 ng/mL and WT = 70 kg, RA model (1/day)") # Table S6 THETA(26) FIX; Table S1 kdeg = exp(theta26) = 1.34
    lkss <- fixed(-1.422316790); label("Quasi-steady-state constant KSS (ug/mL)") # Table S6 THETA(27) FIX; Table S1 KSS = exp(theta27) = 0.241
    e_base_kint <- fixed(0.192887381); label("Power exponent of BASE/40 on kint (unitless)") # Table S6 THETA(28) FIX; Table S1 'kint,sBASE' = 0.193
    e_base_kdeg <- fixed(-0.527820349); label("Power exponent of BASE/40 on kdeg (unitless)") # Table S6 THETA(29) FIX; Table S1 'kdeg,sBASE' = -0.528
    e_wt_kdeg <- fixed(-0.078069784); label("Power exponent of WT/70 on kdeg (unitless)") # Table S6 THETA(30) FIX; Table S1 'kdeg,WT' = -0.0781

    # Parameters estimated in the CRS population (Table S5).
    lrbase <- log(29.4); label("Typical baseline total sIL-6R when no observed baseline is available (ng/mL)") # Table S5 sBASE = exp(theta24) = 29.4
    e_crs_grade_ge2_cl <- log(1.60); label("Log factor on CL for a first CRS event of grade >= 2 (unitless)") # Table S5 'CL CRS grade>1' = exp(theta31) = 1.60
    e_crs_vc <- log(1.15); label("Log factor on VC and VP in patients with CRS (unitless)") # Table S5 'VC,CRS = VP,CRS' = exp(theta32) = 1.15
    e_crs_q <- log(3.59); label("Log factor on Q in patients with CRS (unitless)") # Table S5 'Q CRS' = exp(theta33) = 3.59
    e_crs_kint <- log(2.09); label("Log factor on kint in patients with CRS (unitless)") # Table S5 'kint,CRS' = exp(theta34) = 2.09
    e_crs_kdeg <- log(1.76); label("Log factor on kdeg in patients with CRS (unitless)") # Table S5 'Kdeg,CRS' = exp(theta35) = 1.76

    # Inter-individual variability. Etas on Q, VM, KM, kint and KSS are
    # 0 FIX in Table S6 and are omitted.
    etalcl ~ 0.291 # Table S5 omega2 CL = 0.291
    etalvc ~ 0.0681 # Table S5 omega2 'VC,VP' = 0.0681; one eta shared by VC and VP (Table S6 V3 uses ETA(2))
    etalka ~ fixed(0.205326751) # Table S6 OMEGA 7 FIX; Table S1 omega2 ka = 0.205
    etalfdepot ~ fixed(0.008794452) # Table S6 OMEGA 8 FIX; Table S1 omega2 FSC = 0.00879
    etalrbase ~ fixed(0.381151468) # Table S6 OMEGA 9 FIX; Table S1 omega2 sBASE = 0.381
    etalkdeg ~ fixed(0.012983219) # Table S6 OMEGA 11 FIX; Table S1 omega2 kdeg = 0.013
    etaruv ~ fixed(0.127645999) # Table S6 OMEGA 13 FIX; Table S1 omega2 EPS,TCZ = 0.128; IIV on the tocilizumab residual-error magnitude
    etaruv_total_target ~ fixed(0.353774965) # Table S6 OMEGA 14 FIX; Table S1 omega2 EPS,sIL6R = 0.354; IIV on the sIL-6R residual-error magnitude

    # sIL-6R residual error: Y = log(RTOT) + exp(ETA(14)) * EPS(2), SIGMA(2,2) = 0.015.
    expSd_total_target <- fixed(sqrt(0.014993567)); label("Log-scale residual SD of total sIL-6R (unitless)") # Table S6 $SIGMA 2 = 0.014993567 FIX (variance); Table S1 sigma2 sIL6R = 0.015
  })

  model({
    # ---- Covariate transforms (Table S6 $PK) ------------------------------
    logwt <- log(WT / 70)
    loghdl <- log(HDLC / 1.5)
    crcl_trunc <- min(max(CRCL, 15), 150)
    logcrcl <- log(crcl_trunc / 100)
    logage <- log(AGE / 50)
    logclage <- log(1 - e_age_cl * exp(-log(2) * AGE / thalf_age_cl))

    # ---- Individual PK parameters (Table S6 MU_1 .. MU_8) -----------------
    cl <- exp(lcl + e_wt_cl * logwt + e_hdlc_cl * loghdl + logclage + e_crs_grade_ge2_cl * CRS_GRADE_GE2 + etalcl)
    vc <- exp(lvc + e_wt_vc * logwt + etalvc) * exp(e_crs_vc)
    vp <- exp(lvp + e_wt_vc * logwt + etalvc) * exp(e_crs_vc)
    q <- exp(lq + e_wt_cl * logwt + e_crs_q)
    vmax <- exp(lvmax + e_crcl_vmax * logcrcl)
    km <- exp(lkm)
    ka <- exp(lka + e_age_ka * logage + e_study_25220_28117_ka * STUDY_25220_28117 + etalka)
    ff_sc_inj <- 1 + (e_sc_inj_group3_lfdepot - 1) * SC_INJ_GROUP3
    fdepot <- exp(lfdepot * AGE / (age50_fdepot + AGE) * ff_sc_inj + etalfdepot)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- sIL-6R parameters (Table S6 MU_9 .. MU_12) -----------------------
    # Observed baseline when supplied (> 0), otherwise the model baseline.
    rbase_model <- exp(lrbase + etalrbase)
    if (TSIL6R_BASE > 0) {
      rbase <- TSIL6R_BASE
    } else {
      rbase <- rbase_model
    }
    kint <- exp(lkint + e_crs_kint) * (rbase / 40)^e_base_kint
    kdeg <- exp(lkdeg + e_wt_kdeg * logwt + etalkdeg) * exp(e_crs_kdeg) * (rbase / 40)^e_base_kdeg
    kss <- exp(lkss)

    Cc <- central / vc
    # Tocilizumab-bound sIL-6R from the QSS relation (Table S6 $DES RC).
    bound_target <- total_target * Cc / (kss + Cc)

    # ---- ODEs (Table S6 $DES) ---------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - vmax * central / (km + Cc) - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(total_target) <- kdeg * (rbase - total_target) - (kint - kdeg) * bound_target
    f(depot) <- fdepot
    total_target(0) <- rbase

    # ---- Derived outputs ---------------------------------------------------
    free_target <- total_target - bound_target
    ro <- 100 * Cc / (kss + Cc) # sIL-6R occupancy (%)

    # ---- Residual error (Table S6 $ERROR) ---------------------------------
    # Tocilizumab: Y = log(C) + W * EPS(1), W = (SDL-(SDL-SDH)*C/(SD50+C))*exp(ETA(13)), SIGMA(1,1) = 1 FIX.
    W <- (sdL - (sdL - sdH) * Cc / (sd50 + Cc)) * exp(etaruv)
    Cc ~ lnorm(W)
    # sIL-6R: Y = log(RTOT) + exp(ETA(14)) * EPS(2).
    W_total_target <- expSd_total_target * exp(etaruv_total_target)
    total_target ~ lnorm(W_total_target)
  })
}
