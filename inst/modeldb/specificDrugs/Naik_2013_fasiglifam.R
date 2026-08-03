Naik_2013_fasiglifam <- function() {
  description <- paste(
    "Joint population PK / FPG / HbA1c model of fasiglifam (TAK-875, a novel",
    "GPR40 / free-fatty-acid-receptor-1 agonist) in adults with type 2",
    "diabetes mellitus (T2DM) inadequately controlled on metformin (Naik",
    "2013 CPT PSP). Structural PK is a two-compartment model with",
    "first-order absorption and linear elimination (Q/F, V2/F, ka FIXED to",
    "values estimated in an earlier multiple-rising-dose study because the",
    "phase-2 sparse-sampling design could not identify them). Sex is an",
    "additive covariate on CL/F (males ~41 pct higher). Drug effect on",
    "fasting plasma glucose (FPG) is a semi-mechanistic indirect-response",
    "model with Emax stimulation of Kout; Emax depends exponentially on",
    "baseline FPG (BFPG) and on baseline aspartate aminotransferase (AST",
    "= SGOT). HbA1c is driven by FPG via a first-order production /",
    "elimination pair, with a placebo factor LIFE(t) = 1 - MPL * (1 -",
    "exp(-ln(2) * t / HL_pl)) that scales HbA1c production down over the",
    "12-week trial (MPL half-life FIXED to 720 h from graphical analysis).",
    "Baseline HbA1c depends linearly on disease duration (T_DIAG_DIAB)",
    "and exponentially on BFPG. MPL has an additive sex effect (males",
    "have larger placebo response). Data: 1211 PK samples from 286 T2DM",
    "patients and 2710 FPG + 1381 HbA1c samples from 346 T2DM patients",
    "on 6.25-200 mg oral once-daily fasiglifam for 12 weeks.",
    "Fasiglifam was subsequently withdrawn from phase-3 development",
    "in December 2013 for hepatotoxicity signals; this popPK/PD model",
    "characterises the phase-2 dose-selection analysis that preceded",
    "that finding."
  )
  reference <- paste(
    "Naik H, Lu J, Cao C, Pfister M, Vakilynejad M, Leifke E.",
    "Pharmacometric Approaches to Guide Dose Selection of the Novel",
    "GPR40 Agonist TAK-875 in Subjects With Type 2 Diabetes Mellitus.",
    "CPT: Pharmacometrics & Systems Pharmacology (2013) 2, e22;",
    "doi:10.1038/psp.2012.23. Supplement PSP-2012-0026-T contains",
    "the three NONMEM control streams reproduced here (s07 PK,",
    "s06 FPG, s05 HbA1c). Q/F, V2/F, and ka were FIXED per the",
    "phase-2 sparse-sampling rationale to values estimated in the",
    "phase-1 multiple-rising-dose study of Leifke E et al.,",
    "Clin. Pharmacol. Ther. 92, 29-39 (2012),",
    "doi:10.1038/clpt.2012.43.",
    sep = " "
  )
  vignette <- "Naik_2013_fasiglifam"
  paper_specific_compartments <- c("glucose", "Hba1c")

  units <- list(
    time          = "h",
    dosing        = "mg fasiglifam (TAK-875, oral, once daily)",
    concentration = paste(
      "Cc in mg/L (equivalent to ug/mL; EC50 = 3.16 ug/mL = 3.16 mg/L);",
      "FPG in mg/dL; HbA1c in %"
    )
  )

  covariateData <- list(
    SEXF = list(
      description        = "Sex (female indicator)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female is the reference / typical subject)",
      notes              = paste(
        "Additive covariate effects. Naik 2013 reports SEX = 1 for male",
        "and SEX = 0 for female (source column SEX); the canonical SEXF",
        "uses SEXF = 1 for female. The model applies the paper's",
        "coefficients to (1 - SEXF) so that SEXF = 1 (female) reproduces",
        "the paper's typical CL/F = 0.75 L/h and MPL = 0.0590 whereas",
        "SEXF = 0 (male) yields CL/F = 0.75 + 0.31 = 1.06 L/h and",
        "MPL = 0.0590 + 0.0363 = 0.0953."
      ),
      source_name        = "SEX (1 = male, 0 = female in Naik 2013 supplement)"
    ),
    FPG = list(
      description        = "Observed baseline fasting plasma glucose at study entry",
      units              = "mg/dL (paper convention; 1 mmol/L = 18.02 mg/dL)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline FPG (per-subject, time-fixed) is used as a covariate",
        "in three places: exponential effect on Emax of the FPG",
        "indirect-response model (Table 2 row 'BFPG on Emax' = 0.00746",
        "per (mg/dL); centred at 163.5); exponential effect on baseline",
        "HbA1c BLA1 (Table 2 row 'BFPG on BLA1' = 0.00181 %/(mg/dL);",
        "centred at 163.5); anchor for the steady-state KIG =",
        "BLA1 * KA1C / FPG in the HbA1c ODE. Distinct from the modelled",
        "glucose(t) state. Reference value 163.5 mg/dL is the typical",
        "subject in Naik 2013 Figure 2 tornado plot."
      ),
      source_name        = "BFPG (paper narrative) / FPG (source dataset column, baseline value)"
    ),
    AST = list(
      description        = "Serum aspartate aminotransferase (baseline)",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Exponential effect on Emax of the FPG indirect-response model",
        "(Table 2 row 'AST on Emax' = 0.00731 per (U/L); centred at",
        "21 U/L). Naik 2013 supplement control stream s06 encodes this",
        "as EXP(theta5 * (FPG - 163.5) + theta6 * (SGOT - 21)) whereas",
        "the Discussion prose describes it as 'linear'. The control",
        "stream is the definitive source and the exponential form is",
        "used here (see vignette Errata). Naik 2013 Discussion also",
        "notes that removing seven patients with AST >= 70 U/L removed",
        "the statistical significance of this effect (P >= 0.005);",
        "the coefficient is retained here as reported for the full",
        "cohort. Reference value 21 U/L is the typical subject in",
        "Naik 2013 Figure 2 tornado plot."
      ),
      source_name        = "SGOT (serum glutamic-oxaloacetic transaminase; legacy name for AST)"
    ),
    T_DIAG_DIAB = list(
      description        = "Time since T2DM diagnosis at study entry",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Linear (additive) effect on baseline HbA1c BLA1 (Table 2 row",
        "'DD on BLA1' = 0.0133 %/year; centred at 4.61 years). Naik",
        "2013 Discussion notes 'Disease duration had no significant",
        "effect as a covariate on PK or efficacy parameters if",
        "controlled for BLA1'; the reported effect is retained in the",
        "BLA1 equation as encoded in the supplement control stream",
        "s05. Reference value 4.61 years is the typical subject in",
        "Naik 2013 Figure 2 tornado plot. Cohort range 0.39 - 14.75",
        "years (Naik 2013 Discussion)."
      ),
      source_name        = "DD (disease duration)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Age at study entry (screened, not retained in final model)",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened against PK and efficacy parameters via forward-",
        "addition + backward-elimination stepwise covariate selection",
        "(Wahlby 2001 protocol). Not retained in the final PK, PK-FPG,",
        "or PK-HbA1c models. Cohort mean 51.5 +/- 10.6 years (range",
        "21 - 79) per Naik 2013 Table 1."
      ),
      source_name        = "AGE"
    ),
    WT = list(
      description        = "Body weight at baseline (screened, not retained in final model)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened; not retained. Cohort mean 86.0 kg (range 49.5 -",
        "172.7) per Naik 2013 Table 1."
      ),
      source_name        = "WT"
    ),
    BMI = list(
      description        = "Body mass index at baseline (screened, not retained in final model)",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened; not retained. Cohort mean 31.7 +/- 5.1 kg/m^2 per",
        "Naik 2013 Table 1."
      ),
      source_name        = "BMI"
    ),
    CRCL = list(
      description        = "Creatinine clearance (screened, not retained in final model)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened; not retained. Cohort mean 129.7 +/- 44.0 mL/min per",
        "Naik 2013 Table 1. Renal excretion of fasiglifam is minimal",
        "(CLr <= 0.003 L/h; Naik 2013 Introduction)."
      ),
      source_name        = "CRCL"
    ),
    ALB = list(
      description        = "Serum albumin (screened, not retained in final model)",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened; not retained. Cohort mean 41.8 +/- 3.3 g/L per",
        "Naik 2013 Table 1."
      ),
      source_name        = "ALB (source column ALBU)"
    ),
    CONMED_METFORMIN = list(
      description        = "Concomitant metformin (screened, not retained in final model)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no metformin)",
      notes              = paste(
        "Screened; not retained. 76.3 pct of the PK/PD cohort was on",
        "background metformin at baseline (Naik 2013 Table 1). Naik",
        "2013 Discussion attributes the lack-of-effect to prior",
        "metformin steady state at study entry."
      ),
      source_name        = "METF (paper reports metformin dose 1696 +/- 442 mg for 273 subjects; the binary Yes/No indicator was screened)"
    )
  )

  population <- list(
    species                 = "human",
    n_subjects_pk           = 286L,
    n_subjects_pkpd         = 346L,
    n_pk_observations       = 1211L,
    n_fpg_observations      = 2710L,
    n_hba1c_observations    = 1381L,
    n_studies               = 1L,
    study_design            = paste(
      "Phase 2, randomised, double-blind, placebo- and active-comparator",
      "(glimepiride)-controlled, parallel-group, multicentre study",
      "(Burant 2012 Lancet). Seven treatment groups: placebo, TAK-875",
      "6.25, 12.5, 25, 50, 100, 200 mg once-daily for 12 weeks. All",
      "patients were on stable metformin monotherapy for >= 8 weeks",
      "before screening, or on diet + exercise alone for >= 8 weeks."
    ),
    age_range               = "21 - 79 years (mean 51.5 +/- 10.6)",
    weight_range            = "49.5 - 172.7 kg (mean 86.0 kg)",
    sex_female_pct          = 53.5,
    race_ethnicity          = c(Caucasian = 82.4, AfricanAmerican = 10.1, Asian = 3.2, Other = 4.3),
    ethnicity_hispanic_pct  = 66.8,
    disease_state           = paste(
      "Type 2 diabetes mellitus (T2DM) inadequately controlled on stable",
      "metformin monotherapy or on diet + exercise alone. Baseline FPG",
      "mean 170.3 +/- 51.1 mg/dL; baseline HbA1c mean 8.4 +/- 0.93 pct.",
      "Disease duration mean 5.7 +/- 4.9 years (range 0.39 - 14.75)."
    ),
    dose_range              = paste(
      "Fasiglifam 6.25, 12.5, 25, 50, 100, and 200 mg orally once daily",
      "for 12 weeks. Placebo arm was included in the PK-efficacy",
      "dataset but not in the PK-only dataset."
    ),
    regions                 = "Not tabulated in Naik 2013; the parent Burant 2012 phase-2 trial was US / Latin America / Europe multicentre.",
    notes                   = paste(
      "Data from 346 T2DM patients (PK-efficacy dataset) and 286 T2DM",
      "patients (PK dataset). Naik 2013 Table 1 tabulates baseline",
      "demographics and clinical laboratory values. Population-PK-PD",
      "modelling was performed sequentially (PK first, then PK-FPG",
      "with fixed PK parameters, then PK-HbA1c with fixed PK-FPG",
      "parameters); the joint dynamics are reproduced here in a",
      "single nlmixr2 file. The subsequent phase-3 development of",
      "fasiglifam was terminated in December 2013 for signals of",
      "hepatotoxicity (an outcome that post-dates this popPK/PD",
      "analysis)."
    )
  )

  ini({
    # ================================================================
    # PK model (Naik 2013 Table 2, 'Pharmacokinetic analysis' block;
    # supplement control stream s07 PSP-2012-0026-T-s07.doc). Two-
    # compartment oral, first-order absorption + linear elimination.
    # Reference typical subject: female (SEXF = 1).
    # ================================================================
    lcl <- log(0.75)             ; label("Apparent oral clearance CL/F (L/h) for typical female subject")   # Table 2 CL/F = 0.75 L/h (7.46 pct RSE); s07 TH1
    cl_male_delta <- 0.31        ; label("Additive shift on CL/F for male (L/h); paper reports males have ~41 pct higher clearance than females")  # Table 2 'Gender on CL/F' = 0.31 L/h (27.2 pct RSE); s07 TH6 additive form
    lvc <- log(5.86)             ; label("Apparent central volume V1/F (L)")                                  # Table 2 V1/F = 5.86 L (11.2 pct RSE); s07 TH2
    lq  <- fixed(log(0.833))     ; label("Apparent intercompartmental clearance Q/F (L/h; per Naik 2013 Discussion because phase-2 sparse sampling could not identify Q, V2, ka; value from Leifke 2012 multiple-rising-dose study)")  # Table 2 Q/F = 0.833 L/h FIXED; s07 TH3 FIX
    lvp <- fixed(log(23.7))      ; label("Apparent peripheral volume V2/F (L; per Naik 2013 Discussion; value from Leifke 2012)")   # Table 2 V2/F = 23.7 L FIXED; s07 TH4 FIX
    lka <- fixed(log(0.075))     ; label("First-order absorption rate constant ka (1/h; per Naik 2013 Discussion; value from Leifke 2012)")  # Table 2 Ka = 0.075 1/h FIXED; s07 TH5 FIX

    etalcl ~ 0.488               # Table 2 omega^2 on CL/F = 0.488 (7.45 pct RSE); s07 OMEGA1 EXP
    propSd <- sqrt(0.152)        ; label("Proportional residual SD for fasiglifam plasma concentration (fraction)")  # Table 2 sigma^2 (exponential error label; proportional error per Methods) = 0.152 (3.67 pct RSE); s07 SIGMA1

    # ================================================================
    # FPG model (Naik 2013 Table 2, 'Pharmacokinetic-FPG analysis'
    # block; supplement control stream s06 PSP-2012-0026-T-s06.doc).
    # Indirect response with Emax stimulation of Kout.
    # Reference typical subject: BFPG = 163.5 mg/dL, AST = 21 U/L.
    # ================================================================
    lbl_fpg <- log(164)          ; label("Model-predicted baseline FPG BL (mg/dL) for the typical subject")   # Table 2 BL = 164 mg/dL (1.5 pct RSE); s06 TH1
    lkout_fpg <- log(0.00542)    ; label("First-order rate constant for FPG elimination Kout (1/h)")           # Table 2 Kout = 0.00542 1/h (12.5 pct RSE); s06 TH2
    lemax_fpg <- log(0.366)      ; label("Maximum fractional stimulation of Kout by fasiglifam Emax for typical subject (unitless fraction)")  # Table 2 Emax = 0.366 (16.2 pct RSE); s06 TH3
    lec50_fpg <- log(3.16)       ; label("Fasiglifam plasma concentration producing half-maximal Emax EC50 (mg/L = ug/mL)")  # Table 2 EC50 = 3.16 ug/mL (32.3 pct RSE); s06 TH4
    e_bfpg_emax <- 0.00746       ; label("Exponential coefficient of (FPG - 163.5) on Emax (1/(mg/dL))")        # Table 2 'BFPG on Emax' = 0.00746 (13.1 pct RSE); s06 TH5
    e_ast_emax <- 0.00731        ; label("Exponential coefficient of (AST - 21) on Emax (1/(U/L))")            # Table 2 'AST on Emax' = 0.00731 (44.2 pct RSE); s06 TH6 -- exponential per supplement control stream s06, despite paper text 'increased linearly with AST' (see vignette Errata)

    etalbl_fpg ~ 0.0610          # Table 2 omega^2 on BL = 0.0610 (9.3 pct RSE); s06 OMEGA1 EXP
    etalemax_fpg ~ 0.0854        # Table 2 omega^2 on Emax = 0.0854 (53.9 pct RSE); s06 OMEGA2 EXP
    propSd_glucose <- sqrt(0.0187) ; label("Proportional residual SD for FPG (fraction; additive on log scale per Methods = proportional in linear space)")  # Table 2 sigma^2 (additive on log-transformed scale) = 0.0187 (1.39 pct RSE); s06 SIGMA1

    # ================================================================
    # HbA1c model (Naik 2013 Table 2, 'Pharmacokinetic-efficacy
    # analysis' block; supplement control stream s05
    # PSP-2012-0026-T-s05.doc). First-order HbA1c production driven
    # by FPG(t) times the exponential-decay placebo factor LIFE(t) =
    # 1 - MPL * (1 - exp(-ln(2) * t / HL_pl)); first-order removal
    # by KA1C. KIG = BLA1 * KA1C / FPG is the steady-state anchor.
    # Reference typical subject: female, DD = 4.61 y, BFPG = 163.5 mg/dL.
    # ================================================================
    bla1 <- 8.25                 ; label("Baseline HbA1c BLA1 (pct) for the typical subject")               # Table 2 BLA1 = 8.25 pct (0.5 pct RSE); s05 TH1
    e_dd_bla1 <- 0.0133          ; label("Linear (additive) coefficient of (T_DIAG_DIAB - 4.61) on BLA1 (pct / year)")   # Table 2 'DD on BLA1' = 0.0133 pct/year (55.6 pct RSE); s05 TH6
    e_bfpg_bla1 <- 0.00181       ; label("Exponential coefficient of (FPG - 163.5) on BLA1 (1/(mg/dL))")     # Table 2 'BFPG on BLA1' = 0.00181 pct/(mg/dL) (5.1 pct RSE); s05 TH7 -- entered in supplement s05 as EXP(TH7 * (BLI - 163.5)) so 0.00181 is an exponential coefficient (see vignette Errata)
    lka1c <- log(0.00052)        ; label("First-order rate constant for HbA1c elimination KA1C (1/h)")       # Table 2 KA1C = 0.00052 1/h (11.3 pct RSE); s05 TH2
    hl_placebo <- fixed(720)     ; label("Placebo-effect half-life HL (h; per Naik 2013 Methods: '720 h observed based upon graphical analysis of data')")   # Table 2 HL = 720 h FIXED; s05 TH3 FIX
    mpl <- 0.0590                ; label("Maximum placebo response MPL for typical female subject (unitless)")   # Table 2 MPL = 0.0590 (21.7 pct RSE); s05 TH4
    mpl_male_delta <- 0.0363     ; label("Additive shift on MPL for male (unitless)")                        # Table 2 'Gender on MPL' = 0.0363 (32.1 pct RSE); s05 TH5

    # omega^2 BLA1 = 0.0057 (9.4 pct RSE) s05 OMEGA(diag) EXP; block-2:
    # omega^2 KA1C = 0.95 (17.6 pct RSE), cov(MPL, KA1C) = 0.506
    # (35.0 pct RSE), omega^2 MPL = 0.305 (44.9 pct RSE) -- Table 2
    # 'Pharmacokinetic-efficacy analysis' rows + s05 OMEGA + OMEGA
    # BLOCK(2).
    etabla1 ~ 0.0057
    etalka1c + etampl ~ c(0.95,
                          0.506, 0.305)
    propSd_Hba1c <- sqrt(0.00164); label("Proportional residual SD for HbA1c (fraction; additive on log scale per Methods = proportional in linear space)")   # Table 2 sigma^2 (additive on log-transformed scale) = 0.00164 (3.3 pct RSE); s05 SIGMA1
  })

  model({
    # ================================================================
    # Reference values for centred covariate effects (Naik 2013
    # Figure 2 tornado-plot reference-line: 52-year-old female with
    # AST = 21 U/L, BFPG = 163.5 mg/dL, DD = 4.61 years).
    # ================================================================
    ref_bfpg <- 163.5       # mg/dL
    ref_ast  <- 21.0        # U/L
    ref_dd   <- 4.61        # years

    # ================================================================
    # Individual PK parameters. NONMEM control stream s07 encodes
    #   TVCL = THETA(1) + THETA(6) * SEX_MALE
    #   CL   = TVCL * EXP(ETA(1))
    # where SEX_MALE = 1 for male, 0 for female. The canonical SEXF
    # is 1 for female, so SEX_MALE = 1 - SEXF.
    # ================================================================
    cl_typ <- exp(lcl) + cl_male_delta * (1.0 - SEXF)
    cl <- cl_typ * exp(etalcl)
    vc <- exp(lvc)
    q  <- exp(lq)
    vp <- exp(lvp)
    ka <- exp(lka)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ================================================================
    # PK ODEs (2-compartment oral).
    # ================================================================
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Observed concentration Cc in mg/L (equivalent to ug/mL); EC50
    # is stated as 3.16 ug/mL = 3.16 mg/L so no unit-conversion factor
    # is needed between Cc and EC50.
    Cc <- central / vc

    # ================================================================
    # Individual FPG-model parameters (s06). Emax is exponentially
    # modulated by (FPG - 163.5) and (AST - 21). Baseline FPG (the
    # covariate `FPG`, per-subject observed BFPG) enters both here
    # and in the HbA1c KIG steady-state anchor below.
    # ================================================================
    bl_fpg_typ <- exp(lbl_fpg)
    bl_fpg_i   <- bl_fpg_typ * exp(etalbl_fpg)

    kout_fpg <- exp(lkout_fpg)
    kin_fpg  <- kout_fpg * bl_fpg_i

    emax_typ <- exp(lemax_fpg) *
                exp(e_bfpg_emax * (FPG - ref_bfpg) +
                    e_ast_emax  * (AST - ref_ast))
    emax_i <- emax_typ * exp(etalemax_fpg)

    ec50 <- exp(lec50_fpg)

    # Stimulation of Kout by fasiglifam (Naik 2013 Methods Eq. block):
    #   dFPG/dt = KIN - KOUT * (1 + Emax * Cc / (EC50 + Cc)) * FPG
    stim_fpg <- 1.0 + emax_i * Cc / (ec50 + Cc)

    d/dt(glucose) <- kin_fpg - kout_fpg * stim_fpg * glucose
    glucose(0)    <- bl_fpg_i

    # ================================================================
    # Individual HbA1c-model parameters (s05).
    # BLA1 = (8.25 + 0.0133 * (DD - 4.61)) * exp(0.00181 * (FPG - 163.5))
    # MPL  = 0.0590 + 0.0363 * (1 - SEXF)      # additive male shift
    # KIG  = BLA1 * KA1C / FPG                 # steady-state anchor
    # LIFE(t) = 1 - MPL * (1 - exp(-ln(2) * t / HL_pl))
    # dHbA1c/dt = KIG * LIFE(t) * glucose - KA1C * HbA1c
    # ================================================================
    bla1_typ_dd  <- bla1 + e_dd_bla1 * (T_DIAG_DIAB - ref_dd)
    bla1_typ_fpg <- bla1_typ_dd * exp(e_bfpg_bla1 * (FPG - ref_bfpg))
    bla1_i       <- bla1_typ_fpg * exp(etabla1)

    ka1c_i <- exp(lka1c + etalka1c)

    mpl_typ <- mpl + mpl_male_delta * (1.0 - SEXF)
    mpl_i   <- mpl_typ * exp(etampl)

    kig_i <- bla1_i * ka1c_i / FPG

    # Placebo factor LIFE(t): at t = 0, LIFE = 1 (no placebo effect);
    # as t -> Inf, LIFE -> 1 - MPL. HL_pl is the half-life of the
    # placebo-approach.
    life_t <- 1.0 - mpl_i * (1.0 - exp(-log(2.0) * t / hl_placebo))

    d/dt(Hba1c) <- kig_i * life_t * glucose - ka1c_i * Hba1c
    Hba1c(0)    <- bla1_i

    # ================================================================
    # Residual error (proportional in linear space; Methods:
    # 'additive on log-transformed scale' == proportional on
    # untransformed scale).
    # ================================================================
    Cc      ~ prop(propSd)
    glucose ~ prop(propSd_glucose)
    Hba1c   ~ prop(propSd_Hba1c)
  })
}
