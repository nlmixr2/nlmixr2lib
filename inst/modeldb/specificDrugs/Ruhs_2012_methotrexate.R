Ruhs_2012_methotrexate <- function() {
  description <- "PK/PD model of methotrexate (MTX) and homocysteine (HCY) after high-dose MTX treatment in children with acute lymphoblastic leukemia (Ruhs 2012). Two-compartment IV PK for MTX with linear BSA scaling on CL, V1, Q, V2 (theta values reported per m^2 BSA in the paper) and a power effect of the age- and gender-adjusted serum creatinine ratio (CREAT_REF / CREAT) on CL (Eq. 1); coupled to a single-compartment indirect response model for HCY where MTX inhibits the HCY elimination rate kout via an inverse Emax function (Emax fixed to 1) and the typical HCY baseline depends linearly on age. IOV reported in the paper on CL (17.15% CV) and HCYBL (23.83% CV) across the window and four consolidation HDMTX administrations is not encoded in the model file (only between-subject IIV is carried)."
  reference <- "Ruhs H, Becker A, Drescher A, Panetta JC, Pui CH, Relling MV, Jaehde U. Population PK/PD Model of Homocysteine Concentrations after High-Dose Methotrexate Treatment in Patients with Acute Lymphoblastic Leukemia. PLoS ONE. 2012;7(9):e46015. doi:10.1371/journal.pone.0046015"
  vignette <- "Ruhs_2012_methotrexate"
  units <- list(
    time          = "hour",
    dosing        = "mg",
    concentration = "umol/L",
    notes         = "Dose entered as MTX amount in mg into the central compartment. The observation Cc is reported in umol/L (uM) via the conversion (central/vc) * 1000 / MW_MTX (MW_MTX = 454.44 g/mol) so that EC50 and HCY remain in their published uM units."
  )

  covariateData <- list(
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear multiplicative scaling on CL, V1, Q, V2 (theta values in Table 2 are reported per m^2 BSA; the implicit BSA exponent is 1). Index dataset median 0.83 m^2 (range 0.40-2.97) per Ruhs 2012 Table 1. The paper does not state which BSA formula was used (DuBois, Mosteller, Haycock); record as unspecified.",
      source_name        = "BSA"
    ),
    CREAT = list(
      description        = "Measured serum creatinine",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Patient's measured serum creatinine (Ruhs 2012 Table 1: median 0.4 mg/dL, range 0.1-1.2). Enters the model only via the ratio (CREAT_REF / CREAT)^theta_CR on CL (Eq. 1).",
      source_name        = "CCR"
    ),
    CREAT_REF = list(
      description        = "Age- and gender-adjusted reference serum creatinine",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Externally-computed expected normal SCR for the individual (denoted CCR,adj in Ruhs 2012). The paper cites its reference [23] for the maturation-dependent adjustment to account for age and gender but does not give the explicit formula in the main text. Users must compute CREAT_REF from age and sex before passing to the model. When no covariate value can be derived, set CREAT_REF = CREAT so the renal-function factor evaluates to 1 (no adjustment).",
      source_name        = "CCR,adj"
    ),
    AGE = list(
      description        = "Age at study entry",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at study entry. Enters the model only as a linear additive effect on the typical HCY baseline: HCYBL_tv = theta_BL + theta_BL,AGE * AGE, with theta_BL,AGE = 0.116 uM/year (Ruhs 2012 Table 2). Index dataset median 5.42 years (range 1.03-18.85).",
      source_name        = "AGE"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 494L,
    n_index         = 331L,
    n_evaluation    = 163L,
    n_centres       = 2L,
    n_mtx_obs       = 6722L,
    n_hcy_obs       = 2567L,
    age_range       = "1.03-18.85 years",
    age_median      = "5.42 years (index dataset)",
    weight_range    = "7.8-160.1 kg",
    weight_median   = "22.0 kg (index dataset)",
    bsa_range       = "0.40-2.97 m^2",
    bsa_median      = "0.83 m^2",
    creat_range     = "0.1-1.2 mg/dL",
    sex_female_pct  = 42.9,
    race_ethnicity  = NA,
    disease_state   = "Newly diagnosed acute lymphoblastic leukemia (ALL); low-risk (LR) and standard/high-risk (SHR) subgroups stratified by leukocyte count, DNA index, ETV6-RUNX1 / BCR-ABL1 status and minimal residual disease",
    dose_range      = "Window therapy: 1000 mg/m^2 over 4 h (1 g/m^2) or over 24 h (200 mg/m^2 IV bolus followed by 800 mg/m^2 over 24 h). Consolidation HDMTX every other week for 4 doses, individualized to target MTX Cpss = 33 uM (LR; median 2653 mg/m^2) or 65 uM (SHR; median 4654 mg/m^2); 10% of the individualized dose given over 1 h as loading and the remaining 90% over 23 h. Folinate rescue starting 42 h after start of infusion.",
    regions         = "USA (St. Jude Children's Research Hospital, Memphis TN; Cook Children's Medical Center, Fort Worth TX)",
    notes           = "Demographics from Ruhs 2012 Table 1. The PK model was developed on the index dataset (331 patients) and externally evaluated on a separate evaluation dataset (163 patients). The paper estimated interoccasion variability (IOV) of 17.15% CV on CL and 23.83% CV on HCYBL across the window and four consolidation HDMTX administrations; that IOV is not encoded in this model file (only between-subject IIV is carried)."
  )

  ini({
    # ---------------------------------------------------------------
    # MTX PK structural parameters (Ruhs 2012 Table 2; per m^2 BSA).
    # ---------------------------------------------------------------
    lcl <- log(6.68);  label("Typical MTX clearance per m^2 BSA (L/h/m^2)")                 # Ruhs 2012 Table 2 theta_CL = 6.68
    lvc <- log(18.6);  label("Typical MTX central volume V1 per m^2 BSA (L/m^2)")           # Ruhs 2012 Table 2 theta_V1 = 18.6
    lq  <- log(0.161); label("Typical MTX inter-compartmental clearance per m^2 BSA (L/h/m^2)")  # Ruhs 2012 Table 2 theta_Q = 0.161
    lvp <- log(3.09);  label("Typical MTX peripheral volume V2 per m^2 BSA (L/m^2)")        # Ruhs 2012 Table 2 theta_V2 = 3.09

    # Power exponent on the renal-function factor (CREAT_REF / CREAT) on CL.
    e_creat_cl <- 0.314; label("Power exponent on (CREAT_REF / CREAT) for MTX CL (unitless)")  # Ruhs 2012 Table 2 theta_CR = 0.314; Eq. 1

    # ---------------------------------------------------------------
    # HCY PD structural parameters (Ruhs 2012 Table 2).
    # HCYBL_tv = theta_BL + theta_BL,AGE * AGE so theta_BL is the
    # zero-age intercept; log-normal IIV is applied to this typical
    # value.
    # ---------------------------------------------------------------
    lhcybl      <- log(4.88);  label("Typical HCY baseline at age 0 y (uM)")                # Ruhs 2012 Table 2 theta_BL = 4.88 (intercept of the linear age effect)
    lkout       <- log(0.027); label("Typical HCY elimination rate constant kout (1/h)")    # Ruhs 2012 Table 2 theta_kout = 0.027
    lec50       <- log(0.648); label("Typical MTX EC50 for inhibition of HCY elimination (uM)")  # Ruhs 2012 Table 2 theta_EC50 = 0.648
    emax        <- fixed(1);   label("MTX maximum fractional inhibition of HCY elimination (unitless; fixed)")  # Ruhs 2012 Table 2 theta_Emax = 1 (fixed)
    e_age_hcybl <- 0.116;      label("Linear effect of age on typical HCY baseline (uM/year)")  # Ruhs 2012 Table 2 theta_BL,AGE = 0.116

    # ---------------------------------------------------------------
    # IIV (Ruhs 2012 Table 2 % CV; omega^2 = log(1 + CV^2)).
    # No correlations reported; diagonal Omega.
    # ---------------------------------------------------------------
    etalcl    ~ 0.02089  # log(1 + 0.1453^2); Table 2 CL    IIV 14.53% CV
    etalvc    ~ 0.03199  # log(1 + 0.1803^2); Table 2 V1    IIV 18.03% CV
    etalq     ~ 0.12584  # log(1 + 0.3661^2); Table 2 Q     IIV 36.61% CV
    etalhcybl ~ 0.03528  # log(1 + 0.1895^2); Table 2 HCYBL IIV 18.95% CV
    etalkout  ~ 0.10472  # log(1 + 0.3323^2); Table 2 kout  IIV 33.23% CV
    etalec50  ~ 0.25285  # log(1 + 0.5363^2); Table 2 EC50  IIV 53.63% CV

    # ---------------------------------------------------------------
    # Residual error (Ruhs 2012 Table 2).
    #   MTX: exponential 0.235 (log-scale SD) -> nlmixr2 prop(propSd)
    #        in the linear small-sigma limit (CV ~ 23.5% either way).
    #   HCY: combined additive 0.911 uM + exponential 0.165 ->
    #        nlmixr2 prop(propSd_HCY) + add(addSd_HCY) in the small-
    #        sigma limit (the simulation-pattern equivalent of the
    #        NONMEM additive + exponential error model used here is
    #        natively supported by nlmixr2 only as additive +
    #        proportional; the difference is < 2% in implied %CV at
    #        the reported sigma values).
    # ---------------------------------------------------------------
    propSd     <- 0.235; label("MTX proportional residual SD (fraction; equivalent to NONMEM exponential error = 0.235 in the small-sigma limit)")  # Ruhs 2012 Table 2 PK exponential error = 0.235
    addSd_HCY  <- 0.911; label("HCY additive residual SD (uM)")                              # Ruhs 2012 Table 2 PD additive error = 0.911
    propSd_HCY <- 0.165; label("HCY proportional residual SD (fraction; equivalent to NONMEM exponential error = 0.165 in the small-sigma limit)")  # Ruhs 2012 Table 2 PD exponential error = 0.165
  })

  model({
    # ---------------------------------------------------------------
    # Constants. MTX molecular weight is used to convert the central
    # MTX amount (mg, the dosing unit) into a uM concentration so that
    # the inhibition function (driven by EC50 in uM) and the published
    # HCY baseline / EC50 use matching units.
    # ---------------------------------------------------------------
    MW_MTX <- 454.44  # methotrexate molecular weight (g/mol = mg/mmol)

    # ---------------------------------------------------------------
    # 1. Renal-function factor on CL (Ruhs 2012 Eq. 1):
    #      CL_i = theta_CL * BSA * (CCR,adj / CCR)^theta_CR
    # CREAT_REF = CCR,adj (age/gender-adjusted reference SCR).
    # ---------------------------------------------------------------
    f_scr <- (CREAT_REF / CREAT)^e_creat_cl

    # ---------------------------------------------------------------
    # 2. Individual MTX PK parameters. BSA scaling is linear (the
    # paper reports thetas per m^2; no estimated exponent on BSA).
    # No IIV on V2 (none reported in Table 2).
    # ---------------------------------------------------------------
    cl <- exp(lcl + etalcl) * BSA * f_scr
    vc <- exp(lvc + etalvc) * BSA
    q  <- exp(lq  + etalq)  * BSA
    vp <- exp(lvp)          * BSA

    # ---------------------------------------------------------------
    # 3. Two-compartment IV micro-constants.
    # ---------------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ---------------------------------------------------------------
    # 4. MTX disposition (IV infusion expected on the central
    # compartment; no depot). Doses in mg -> central amount in mg.
    # ---------------------------------------------------------------
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # MTX plasma concentration in uM (mg/L * 1000 / MW_g/mol).
    Cc <- (central / vc) * 1000 / MW_MTX

    # ---------------------------------------------------------------
    # 5. HCY indirect response (Ruhs 2012 Eq. 2 steady state + Eq. 3
    # dynamics). Inverse-Emax inhibition of kout by MTX:
    #   d/dt(CHCY) = kin - kout * (1 - Emax * Cc/(EC50 + Cc)) * CHCY
    # with kin = kout * HCYBL_i so that CHCY(0) = HCYBL_i is a steady
    # state. The wording in the paper ("Emax is the remaining fraction
    # to which kout will be reduced at the maximal effect of MTX") is
    # consistent with the published numerics: at Emax = 1 and Cc = 33 uM
    # (LR Cpss) the factor (1 - Cc/(EC50 + Cc)) = 1 - 33/33.648 = 0.0193,
    # matching the reported 1.93% residual kout; at 65 uM (SHR Cpss)
    # the factor is 0.0099, matching the reported 0.99%.
    # ---------------------------------------------------------------
    hcybl_tv <- exp(lhcybl) + e_age_hcybl * AGE
    hcybl    <- hcybl_tv * exp(etalhcybl)
    kout     <- exp(lkout + etalkout)
    ec50     <- exp(lec50 + etalec50)
    kin      <- kout * hcybl

    d/dt(effect) <- kin - kout * (1 - emax * Cc / (ec50 + Cc)) * effect

    effect(0) <- hcybl

    HCY <- effect

    Cc  ~ prop(propSd)
    HCY ~ prop(propSd_HCY) + add(addSd_HCY)
  })
}
