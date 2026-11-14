Jamois_2025_tocilizumab_ra <- function() {
  description <- "Joint population PK / soluble IL-6 receptor (sIL-6R) model for intravenous and subcutaneous tocilizumab in adults with rheumatoid arthritis (RA) and children with systemic or polyarticular-course juvenile idiopathic arthritis (Jamois 2025 Table S1; the internally updated reference model that the CRS model adapts). Two-compartment tocilizumab disposition with first-order SC absorption and age-dependent SC bioavailability, parallel linear and Michaelis-Menten elimination, CL maturation with age; total sIL-6R described by the quasi-steady-state approximation of target-mediated drug disposition (an indirect-response model with inhibition of elimination)."
  reference <- "Jamois C, Turner DC, Gibiansky L, Li F, Menuet JF, Pereira LR, Lundberg L, Guerini E, Carlile DJ, Relf J, Wei MC, Kwan A, Frey N, Grange S. Tocilizumab Dosing for Management of T Cell-Engaging Bispecific Antibody-Related CRS in Patients With R/R B-Cell NHL. Clin Pharmacol Ther. 2025;118(4):917. doi:10.1002/cpt.3751"
  vignette <- "Jamois_2025_tocilizumab"
  paper_specific_etas <- c("etaruv", "etaruv_total_target")
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")
  dosing <- c("depot", "central")

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
      notes = "Power effect (HDLC/1.5)^-0.252 on linear CL (Table S1 'CL HDL'; Table S6: LOGHDL = LOG(HDL/1.5)).",
      source_name = "HDL"
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect (ALB/35)^-0.526 on VC and VP (Table S1 'V ALB'; Table S6: LOGALB = LOG(ALB/35)).",
      source_name = "ALB"
    ),
    TPRO = list(
      description = "Baseline total serum protein",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect (TPRO/75)^0.562 on VC and VP (Table S1 'V PROT'; Table S6: LOGPROT = LOG(PROT/75)).",
      source_name = "PROT"
    ),
    CRCL = list(
      description = "BSA-normalized creatinine clearance",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect (CRCL/100)^0.357 on the Michaelis-Menten Vmax, with CRCL truncated to [15, 150] before use (Table S6: CRCLN1 bounds, LOGCRCL = LOG(CRCLN1/100)). Source column CRCLN ('calculated creatinine clearance'); the estimating equation is not stated.",
      source_name = "CRCLN"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Three roles: (1) CL maturation factor 1 - 0.262*exp(-log(2)*AGE/2.44); (2) power effect (AGE/50)^-0.354 on ka; (3) hyperbolic AGE/(12.0 + AGE) scaling of log SC bioavailability (Table S1; Table S6 LOGCLAGE, MU_7, LTPF1).",
      source_name = "AGE"
    ),
    TSIL6R_BASE = list(
      description = "Observed baseline (pre-tocilizumab) serum total soluble IL-6 receptor concentration",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Table S6: BASE = observed SIL6R when available, otherwise the model-predicted typical baseline exp(theta24 + eta9) (sBASE = 47.6 ng/mL in the RA / JIA model, Table S1). Missing values are coded negative in the source data set; supply 0 or a negative value to use the model baseline. BASE is the initial condition of total_target and scales kint and kdeg as (BASE/40)^theta.",
      source_name = "SIL6R"
    ),
    STUDY_25220_28117 = list(
      description = "Indicator for the tocilizumab SC studies NA25220 and 28117",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (all other studies)",
      notes = "Table S6: ST2 = 1 if STUD = 25220 or 28117; multiplies ka by exp(theta16) = 0.807 (Table S1 'ka, NA25220').",
      source_name = "STUD"
    ),
    SC_INJ_GROUP3 = list(
      description = "Indicator for subcutaneous-injection category SINJ = 3 (coding not defined by the source)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (any other SINJ category, and all IV dosing)",
      notes = "Table S6: FF1 = 1; IF(SINJ.EQ.3) FF1 = THETA(17); LTPF1 = THETA(8)*AGE/(THETA(23)+AGE)*FF1. Multiplies the log SC bioavailability by theta17 = 0.692 (Table S1 'F SC,SJIT3'), i.e. raises F toward 1. Neither the paper nor the supplement defines SINJ or its level 3. Acts only on SC doses.",
      source_name = "SINJ"
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = NA_integer_,
    age_range = "Adults and children (pediatric sJIA and pcJIA); range not reported",
    weight_range = "not reported",
    disease_state = "Rheumatoid arthritis (adults), systemic juvenile idiopathic arthritis and polyarticular-course juvenile idiopathic arthritis (children).",
    dose_range = "Intravenous and subcutaneous tocilizumab (regimens not reported); SC studies include NA25220 and 28117.",
    regions = "Global (Roche-sponsored studies)",
    notes = "Internal (otherwise unpublished) update of the tocilizumab PK-sIL-6R model of the co-authors, fit to data from multiple RA, sJIA and pcJIA studies; Jamois 2025 publishes only the parameter estimates (Table S1). The number of studies, subjects and samples and the demographic distributions are not reported."
  )

  ini({
    lcl <- -1.556766391; label("Linear clearance CL at the reference patient (L/day)") # Table S6 THETA(1); Table S1 CL = exp(theta1) = 0.211 L/day
    lvc <- 1.213239670; label("Central volume VC (L)") # Table S6 THETA(2); Table S1 VC = exp(theta2) = 3.36 L
    lvp <- 1.190255642; label("Peripheral volume VP (L)") # Table S6 THETA(3); Table S1 VP = exp(theta3) = 3.29 L
    lq <- -0.937571287; label("Intercompartmental clearance Q (L/day)") # Table S6 THETA(4); Table S1 Q = exp(theta4) = 0.392 L/day
    lvmax <- 0.836626291; label("Michaelis-Menten maximum elimination rate VM (mg/L/day)") # Table S6 THETA(5); Table S1 VM = exp(theta5) = 2.31 mg/L/day
    lkm <- -0.882927120; label("Michaelis-Menten constant KM (ug/mL)") # Table S6 THETA(6); Table S1 KM = exp(theta6) = 0.414 mcg/mL
    lka <- -1.564654589; label("First-order SC absorption rate ka (1/day)") # Table S6 THETA(7); Table S1 ka = exp(theta7) = 0.209 1/day
    lfdepot <- -0.342109352; label("Log SC bioavailability at the adult asymptote (unitless)") # Table S6 THETA(8); Table S1 FSC = exp(theta8) = 0.710
    e_wt_cl <- 0.502061606; label("Power exponent of WT/70 on CL and Q (unitless)") # Table S6 THETA(9); Table S1 'CL WT, Q WT' = 0.502
    e_wt_vc <- 0.749708593; label("Power exponent of WT/70 on VC and VP (unitless)") # Table S6 THETA(10); Table S1 'VC,WT, VP,WT' = 0.750
    e_hdlc_cl <- -0.252391011; label("Power exponent of HDLC/1.5 on CL (unitless)") # Table S6 THETA(11); Table S1 'CL HDL' = -0.252
    e_alb_vc <- -0.526; label("Power exponent of ALB/35 on VC and VP (unitless)") # Table S1 'V ALB' theta12 = -0.526 (Table S6 fixes it to 0 for the CRS fit)
    e_tpro_vc <- 0.562; label("Power exponent of TPRO/75 on VC and VP (unitless)") # Table S1 'V PROT' theta13 = 0.562 (Table S6 fixes it to 0 for the CRS fit)
    e_crcl_vmax <- 0.356767267; label("Power exponent of CRCL/100 on VM (unitless)") # Table S6 THETA(14); Table S1 'VM,CRCLN' = 0.357
    e_age_ka <- -0.354390651; label("Power exponent of AGE/50 on ka (unitless)") # Table S6 THETA(15); Table S1 'ka,AGE' = -0.354
    e_study_25220_28117_ka <- -0.214487419; label("Log factor on ka for studies 25220/28117 (unitless)") # Table S6 THETA(16); Table S1 'ka, NA25220' = exp(theta16) = 0.807
    e_sc_inj_group3_lfdepot <- 0.691991210; label("Multiplicative factor on log SC bioavailability for SINJ = 3 (unitless)") # Table S6 THETA(17); Table S1 'F SC,SJIT3' = 0.692
    sdL <- 0.998582423; label("Log-scale residual SD at low tocilizumab concentrations SDL (unitless)") # Table S6 THETA(18); Table S1 'sigma SDL' = 0.999
    sdH <- 0.116924070; label("Log-scale residual SD at high tocilizumab concentrations SDH (unitless)") # Table S6 THETA(19); Table S1 'sigma SDH' = 0.117
    sd50 <- 2.189973116; label("Tocilizumab concentration at which the residual SD is midway between SDL and SDH (ug/mL)") # Table S6 THETA(20); Table S1 'sigma SD50' = 2.19
    e_age_cl <- 0.261574358; label("Maximum fractional reduction of CL at birth (unitless)") # Table S6 THETA(21); Table S1 'CL maturation' = 0.262
    thalf_age_cl <- 2.440289736; label("Half-life of CL maturation with age (years)") # Table S6 THETA(22); Table S1 't1/2,matur' = 2.44 years
    age50_fdepot <- 12.03537464; label("Age at half the adult log-bioavailability reduction (years)") # Table S6 THETA(23); Table S1 'FSC,Age50' = 12.0 years
    lrbase <- 3.862452507; label("Typical baseline total sIL-6R when no observed baseline is available (ng/mL)") # Table S6 $THETA 24 trailing comment 3.862452507; Table S1 sBASE = exp(theta24) = 47.6
    lkint <- -2.615842819; label("Elimination rate of tocilizumab-bound sIL-6R kint at BASE = 40 ng/mL (1/day)") # Table S6 THETA(25); Table S1 kint = exp(theta25) = 0.0731
    lkdeg <- 0.290626168; label("Degradation rate of unbound sIL-6R kdeg at BASE = 40 ng/mL and WT = 70 kg (1/day)") # Table S6 THETA(26); Table S1 kdeg = exp(theta26) = 1.34
    lkss <- -1.422316790; label("Quasi-steady-state constant KSS (ug/mL)") # Table S6 THETA(27); Table S1 KSS = exp(theta27) = 0.241
    e_base_kint <- 0.192887381; label("Power exponent of BASE/40 on kint (unitless)") # Table S6 THETA(28); Table S1 'kint,sBASE' = 0.193
    e_base_kdeg <- -0.527820349; label("Power exponent of BASE/40 on kdeg (unitless)") # Table S6 THETA(29); Table S1 'kdeg,sBASE' = -0.528
    e_wt_kdeg <- -0.078069784; label("Power exponent of WT/70 on kdeg (unitless)") # Table S6 THETA(30); Table S1 'kdeg,WT' = -0.0781

    # Inter-individual variability (Table S1; Table S6 carries the more precise
    # values of the omegas it holds fixed in the CRS fit).
    etalcl ~ 0.0673 # Table S1 omega2 CL = 0.0673
    etalvc + etalvp ~ c(0.0398, 0.0345, 0.159) # Table S1 omega2 VC = 0.0398, Omega(2,3) = 0.0345 (R = 0.434), omega2 VP = 0.159
    etalq ~ fixed(0.01) # Table S1 omega2 Q = 0.01
    etalvmax ~ fixed(0.01) # Table S1 omega2 VM = 0.01
    etalkm ~ fixed(0.01) # Table S1 omega2 KM = 0.01
    etalka ~ 0.205326751 # Table S6 OMEGA 7; Table S1 omega2 ka = 0.205
    etalfdepot ~ 0.008794452 # Table S6 OMEGA 8; Table S1 omega2 FSC = 0.00879
    etalrbase ~ 0.381151468 # Table S6 OMEGA 9; Table S1 omega2 sBASE = 0.381
    etalkint ~ 0.0186 # Table S1 omega2 kint = 0.0186
    etalkdeg ~ 0.012983219 # Table S6 OMEGA 11; Table S1 omega2 kdeg = 0.013
    etalkss ~ 0.0994 # Table S1 omega2 KSS = 0.0994
    etaruv ~ 0.127645999 # Table S6 OMEGA 13; Table S1 omega2 EPS,TCZ = 0.128; IIV on the tocilizumab residual-error magnitude
    etaruv_total_target ~ 0.353774965 # Table S6 OMEGA 14; Table S1 omega2 EPS,sIL6R = 0.354; IIV on the sIL-6R residual-error magnitude

    # sIL-6R residual error: Y = log(RTOT) + exp(ETA(14)) * EPS(2), SIGMA(2,2) = 0.015.
    expSd_total_target <- sqrt(0.014993567); label("Log-scale residual SD of total sIL-6R (unitless)") # Table S6 $SIGMA 2 = 0.014993567 (variance); Table S1 sigma2 sIL6R = 0.015
  })

  model({
    # ---- Covariate transforms (Table S6 $PK) ------------------------------
    logwt <- log(WT / 70)
    loghdl <- log(HDLC / 1.5)
    logalb <- log(ALB / 35)
    logprot <- log(TPRO / 75)
    crcl_trunc <- min(max(CRCL, 15), 150)
    logcrcl <- log(crcl_trunc / 100)
    logage <- log(AGE / 50)
    logclage <- log(1 - e_age_cl * exp(-log(2) * AGE / thalf_age_cl))

    # ---- Individual PK parameters (Table S6 MU_1 .. MU_8, CRS terms removed)
    cl <- exp(lcl + e_wt_cl * logwt + e_hdlc_cl * loghdl + logclage + etalcl)
    vc <- exp(lvc + e_wt_vc * logwt + e_alb_vc * logalb + e_tpro_vc * logprot + etalvc)
    vp <- exp(lvp + e_wt_vc * logwt + e_alb_vc * logalb + e_tpro_vc * logprot + etalvp)
    q <- exp(lq + e_wt_cl * logwt + etalq)
    vmax <- exp(lvmax + e_crcl_vmax * logcrcl + etalvmax)
    km <- exp(lkm + etalkm)
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
    kint <- exp(lkint + etalkint) * (rbase / 40)^e_base_kint
    kdeg <- exp(lkdeg + e_wt_kdeg * logwt + etalkdeg) * (rbase / 40)^e_base_kdeg
    kss <- exp(lkss + etalkss)

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
