Melhem_2018_g_csf <- function() {
  description <- paste(
    "Semi-mechanistic population PK/PD model for filgrastim and pegfilgrastim",
    "in healthy adult volunteers and adult cancer patients undergoing",
    "chemotherapy-induced neutropenia (CIN). Subcutaneous drug enters a",
    "depot compartment, absorbed at first-order rate KSC into a total-drug",
    "compartment (DT). Target-mediated disposition uses a quasi-equilibrium",
    "quadratic between total drug and total G-CSFR pool (circ), giving free",
    "drug (FDC) and bound complex (RDC). Free drug is cleared linearly",
    "(CLD/VD) and bound drug is internalised (KINT). The granulopoiesis PD",
    "cascade is a five-state receptor chain (SM stem -> MT mitotic -> PM1 ->",
    "PM2 -> RB blood receptors) with baseline set by KP/KTR (bone-marrow",
    "pools) and KP/KC (blood pool); ANC = RB / SR. Drug binding stimulates",
    "receptor production (ST1 = 1 + STM1*driver) and bone-marrow transit",
    "(ST2 = 1 + STM2*driver). Chemotherapy is modelled as a KPD virtual",
    "compartment with lag LAG9; its output (KCHM*CHM) scaled by CHMSL adds",
    "to the mitotic-cell elimination rate. Corticosteroid is a second KPD",
    "compartment with lag LAG8; its output drives an Emax stimulation",
    "modulating ST1 and ST2 through CRTM = KCRT*CRT/(CRT50+KCRT*CRT).",
    "Filgrastim and pegfilgrastim differ on FSC, KSC, VD, CLD and KD",
    "(categorical covariate FORM_GCSF_PEG); pegfilgrastim at 300 ug/kg has",
    "an additional lower CLD (categorical DOSE_PEG_300UGKG); patient vs",
    "healthy-volunteer status shifts KINT and STM2 (categorical DIS_HEALTHY);",
    "body weight allometrically scales CLD (exponent 0.641) and VD (exponent",
    "0.943) with reference 70 kg. Filgrastim (2.6 h) and pegfilgrastim",
    "(10.1 h) apparent half-lives."
  )
  reference <- paste(
    "Melhem M, Delor I, Perez-Ruixo JJ, Harrold J, Chow A, Wu L, Jacqmin P.",
    "(2018). Pharmacokinetic-pharmacodynamic modelling of neutrophil",
    "response to G-CSF in healthy subjects and patients with",
    "chemotherapy-induced neutropenia. Br J Clin Pharmacol 84(5):911-925.",
    "doi:10.1111/bcp.13504.",
    sep = " "
  )
  vignette <- "Melhem_2018_g_csf"

  paper_specific_compartments <- c("depot_kpd_chemotherapy", "depot_kpd_corticosteroid")

  units <- list(
    time          = "hour",
    dosing        = "nmol (G-CSF); mg (chemotherapy and corticosteroid KPD inputs)",
    concentration = "nmol/L (nM; serum G-CSF, endogenous BSLD plus exogenous FDC)",
    ANC           = "10^9 cells/L",
    notes         = paste(
      "G-CSF dose enters compartment 'depot' in nmol (filgrastim MW = 18.8 kDa;",
      "pegfilgrastim MW = 39.0 kDa: 1 ug filgrastim = 1e-6/18800 mol =",
      "5.32e-2 nmol; 1 ug pegfilgrastim = 1e-6/39000 mol = 2.56e-2 nmol).",
      "Chemotherapy dose in mg enters compartment 'depot_kpd_chemotherapy'",
      "(effect delayed by LAG9 = 66.2 h). Corticosteroid dose (1 mg proxy",
      "per Methods) enters compartment 'depot_kpd_corticosteroid' (effect",
      "delayed by LAG8 = 16.3 h). Both KPD compartments decay first-order",
      "and drive stimulation / elimination terms in the granulopoiesis chain."
    )
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric size descriptor scaling filgrastim / pegfilgrastim volume of",
        "distribution (exponent 0.943) and clearance (exponent 0.641) with reference",
        "70 kg (Melhem 2018 Table 2). Continuous covariates were normalised to 70 kg",
        "before entering the power function (Methods, Model evaluation)."
      ),
      source_name        = "WT"
    ),
    FORM_GCSF_PEG = list(
      description        = "G-CSF formulation indicator (1 = pegfilgrastim, 0 = filgrastim)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (filgrastim)",
      notes              = paste(
        "Categorical drug-type covariate on FSC, KSC, VD, CLD and KD",
        "(Melhem 2018 Methods 'Covariate analysis' and Table 2).",
        "Reference: filgrastim, the non-pegylated recombinant human G-CSF.",
        "Pegfilgrastim is filgrastim covalently linked to a 20 kDa monomethoxy-",
        "polyethylene-glycol (mPEG) moiety, giving a longer apparent half-life",
        "(10.1 vs 2.6 h) and different absorption / distribution characteristics.",
        "Auto-approved sibling under the FORM_<DRUG>_<FEATURE> family."
      ),
      source_name        = "DRUG"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-volunteer indicator (1 = healthy adult, 0 = adult cancer patient on chemotherapy)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (adult cancer patient on chemotherapy)",
      notes              = paste(
        "Categorical population-type covariate on KINT and STM2 (Melhem 2018",
        "Methods 'Covariate analysis' and Table 2 rows STM2_HV / STM2_PT and",
        "KINT_HV / KINT_PT). Reference category = adult cancer patient on",
        "chemotherapy (the target population). The two populations pooled by",
        "the model were 110 healthy adult volunteers and 618 adult cancer",
        "patients undergoing chemotherapy (52 paediatric cancer patients were",
        "in Dataset B for covariate exploration but not the primary reference",
        "model of Table 2)."
      ),
      source_name        = "POP"
    ),
    DOSE_PEG_300UGKG = list(
      description        = "Pegfilgrastim 300 ug/kg dose-cohort indicator (1 = subject received pegfilgrastim at 300 ug/kg, 0 = otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (filgrastim, or pegfilgrastim at any dose other than 300 ug/kg)",
      notes              = paste(
        "Subject-level cohort indicator carrying the low-clearance adjustment for",
        "the 300 ug/kg pegfilgrastim cohort (n = 12: 8 healthy volunteers + 4",
        "cancer patients on chemotherapy, per Melhem 2018 Results 'Model",
        "development'). The paper introduced this effect after graphical analysis",
        "of IIV showed that the clearance for 300 ug/kg pegfilgrastim was lower",
        "than for other pegfilgrastim doses; the effect was statistically",
        "significant and considered physiologically plausible despite the small",
        "sample. Encoded as an additive log-scale shift on CLD that applies only",
        "when FORM_GCSF_PEG = 1 (the effect is meaningless for filgrastim).",
        "Auto-approved sibling under the DOSE_HIGH_<drug> family (cf.",
        "DOSE_HIGH_EFL for Jansson 2008 eflornithine and DOSE_400MG for",
        "Jorga 2000 tolcapone)."
      ),
      source_name        = "DOSE_PEG_300"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 703L,
    n_studies      = 10L,
    age_range      = "adults only in the primary Dataset A / final reference model (110 healthy adults + 618 adult cancer patients; 52 paediatric cancer patients in Dataset B were used for covariate exploration only)",
    weight_range   = NA_character_,
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state  = paste(
      "Two pooled populations: (i) healthy adult volunteers receiving",
      "filgrastim or pegfilgrastim (IV and SC routes; formal model",
      "development used SC only); (ii) adult cancer patients receiving",
      "filgrastim or pegfilgrastim SC alongside chemotherapy",
      "(non-small-cell lung cancer, breast cancer, and other solid",
      "tumours per Table 1). All studies sponsored by Amgen Inc."
    ),
    dose_range     = paste(
      "Filgrastim SC: 75, 150, 300, 600 ug once-daily for 10 days in HVs;",
      "5 ug/kg QD in adult cancer patients on chemotherapy (Table 1).",
      "Pegfilgrastim SC: 30, 60, 100, 300 ug/kg SD or QD in HVs; 30, 60,",
      "100, 300 ug/kg once every 3 weeks or 6 mg fixed dose in adult",
      "cancer patients. Filgrastim IV: 375, 750 ug (single dose) and",
      "5 ug/kg (multiple dose) as comparators. Pegfilgrastim IV: 30, 60",
      "ug/kg single dose as comparators. First chemotherapy cycle data",
      "only (Dataset A n = 703)."
    ),
    regions        = NA_character_,
    notes          = paste(
      "Ten in-house Amgen studies. Analysis-dataset naming: Dataset A =",
      "primary structural-development set (703 subjects: 110 HV + 593 adult",
      "patients from first chemotherapy cycle; 6830 PK and 9213 PD",
      "observations); Dataset B = A + 52 paediatric cancer patients (755",
      "subjects; used for additional covariate exploration but not the",
      "primary reference model); Dataset C = A + all chemotherapy cycles",
      "(703 subjects, multiple cycles). The final reference model in",
      "Table 2 is the fit to Dataset A. Estimation software: Monolix 4.3.0",
      "(SAEM algorithm). Sampling distribution for random effects was a",
      "Student t with 5 degrees of freedom."
    )
  )

  ini({
    # -----------------------------------------------------------------
    # PK structural parameters -- filgrastim reference values
    # -----------------------------------------------------------------
    lfsc         <- fixed(log(1))              ; label("Log FSC (subcutaneous bioavailability) for filgrastim (reference)")                                       # Melhem 2018 Table 2 (FSC FIL = 1, Fixed)
    e_peg_fsc    <- log(0.646 / 1)             ; label("Log shift on FSC for pegfilgrastim (FSC PEG / FSC FIL)")                                                  # Melhem 2018 Table 2 (FSC PEG = 0.646, RSE 4%)

    lksc         <- log(0.123)                 ; label("Log KSC (SC absorption rate) for filgrastim (1/h)")                                                      # Melhem 2018 Table 2 (KSC FIL = 0.123, RSE 3%)
    e_peg_ksc    <- log(0.0188 / 0.123)        ; label("Log shift on KSC for pegfilgrastim (KSC PEG / KSC FIL)")                                                 # Melhem 2018 Table 2 (KSC PEG = 0.0188, RSE 2%)

    lvd          <- log(3.12)                  ; label("Log VD (volume of distribution) for filgrastim at WT 70 kg (L)")                                          # Melhem 2018 Table 2 (VD FIL = 3.12 L, RSE 4%)
    e_peg_vd     <- log(5.76 / 3.12)           ; label("Log shift on VD for pegfilgrastim (VD PEG / VD FIL)")                                                    # Melhem 2018 Table 2 (VD PEG = 5.76, RSE 3%)
    e_wt_vd      <- 0.943                      ; label("Allometric exponent on VD relative to WT / 70 kg")                                                       # Melhem 2018 Table 2 (Weight effect - beta_VD = 0.943, RSE 11%; estimated with 95% CI 0.747-1.14)

    lcld         <- log(0.833)                 ; label("Log CLD (linear clearance) for filgrastim at WT 70 kg (L/h)")                                             # Melhem 2018 Table 2 (CLD FIL = 0.833, RSE 4%)
    e_peg_cld    <- log(0.362 / 0.833)         ; label("Log shift on CLD for pegfilgrastim (CLD PEG / CLD FIL)")                                                 # Melhem 2018 Table 2 (CLD PEG = 0.362, RSE 4%)
    e_wt_cld     <- 0.641                      ; label("Allometric exponent on CLD relative to WT / 70 kg")                                                       # Melhem 2018 Table 2 (Weight effect - beta_CLD = 0.641, RSE 16%; estimated with 95% CI 0.445-0.84)
    e_peg300_cld <- log(0.107 / 0.362)         ; label("Log shift on CLD for pegfilgrastim 300 ug/kg cohort (CLD PEG300 / CLD PEG)")                              # Melhem 2018 Table 2 (CLD PEG 300 ug/kg = 0.107, RSE 15%)

    lkd          <- log(0.0237)                ; label("Log KD (dissociation constant of drug-receptor complex) for filgrastim (nM)")                             # Melhem 2018 Table 2 (KD FIL = 0.0237, RSE 8%)
    e_peg_kd     <- log(0.0959 / 0.0237)       ; label("Log shift on KD for pegfilgrastim (KD PEG / KD FIL)")                                                    # Melhem 2018 Table 2 (KD PEG = 0.0959, RSE 5%)

    # -----------------------------------------------------------------
    # Granulopoiesis parameters (drug-shared)
    # -----------------------------------------------------------------
    lkp          <- log(0.0276)                ; label("Log KP (receptor production rate) (nM/h)")                                                                # Melhem 2018 Table 2 (KP = 0.0276, RSE 1%)
    lktr         <- fixed(log(4 / 120))        ; label("Log KTR (bone-marrow transit rate) (1/h) (5-day maturation time)")                      # Melhem 2018 Table 2 (KTR = 0.033) and Methods 'Neutrophil maturation'
    lkc          <- fixed(log(0.12))           ; label("Log KC (blood neutrophil elimination rate) (1/h) (6-hour blood half-life)")              # Melhem 2018 Table 2 (KC = 0.12) and Results 'Covariate analysis'
    lsr          <- fixed(log(0.059))          ; label("Log SR (receptor / ANC scaling factor) (nM per 10^9 cells/L) -- preliminary-run estimate")       # Melhem 2018 Table 2 (SR = 0.059) and Results 'Model development'

    lstm1        <- log(7.53)                  ; label("Log STM1 (max stimulation of receptor production)")                                                       # Melhem 2018 Table 2 (STM1 = 7.53, RSE 2%)
    lstm2        <- log(3.89)                  ; label("Log STM2 (max stimulation of transit rate) -- patient reference")                                         # Melhem 2018 Table 2 (STM2 PT = 3.89, RSE 2%)
    e_hv_stm2    <- log(5.21 / 3.89)           ; label("Log shift on STM2 for healthy volunteers (STM2 HV / STM2 PT)")                                            # Melhem 2018 Table 2 (STM2 HV = 5.21, RSE 4%)

    lkint        <- log(0.113)                 ; label("Log KINT (drug-receptor complex internalisation rate) (1/h) -- patient reference")                        # Melhem 2018 Table 2 (KINT PT = 0.113, RSE 4%)
    e_hv_kint    <- log(0.197 / 0.113)         ; label("Log shift on KINT for healthy volunteers (KINT HV / KINT PT)")                                            # Melhem 2018 Table 2 (KINT HV = 0.197, RSE 7%)

    lbsld        <- log(0.003)                 ; label("Log BSLD (baseline endogenous G-CSF concentration) (nM)")                                                 # Melhem 2018 Table 2 (BSLD = 0.003, RSE 3%)

    # -----------------------------------------------------------------
    # KPD lag times and elimination rates
    # -----------------------------------------------------------------
    ltlag_cort   <- log(16.3)                  ; label("Log LAG8 (corticosteroid KPD lag time) (h)")                                                              # Melhem 2018 Table 2 (LAG8 = 16.3, RSE 9%)
    ltlag_chem   <- log(66.2)                  ; label("Log LAG9 (chemotherapy KPD lag time) (h)")                                                                # Melhem 2018 Table 2 (LAG9 = 66.2, RSE 1%)

    lkel_cort    <- fixed(log(0.2))            ; label("Log KCRT (corticosteroid KPD elimination rate) (1/h)")                                            # Melhem 2018 Table 2 (KCRT = 0.2, Fixed)
    lkel_chem    <- log(0.0724)                ; label("Log KCHM (chemotherapy KPD elimination rate) (1/h)")                                                       # Melhem 2018 Table 2 (KCHM = 0.0724, RSE 1%)

    lchmsl       <- log(668)                   ; label("Log CHMSL (slope relating chemotherapy KPD output to mitotic-cell loss rate) (1/mg)")                      # Melhem 2018 Table 2 (CHMSL = 668, RSE 10%)
    lcrt50       <- log(0.0015)                ; label("Log CRT50 (corticosteroid Emax half-max) (mg/h)")                                                          # Melhem 2018 Table 2 (CRT50 = 0.0015, RSE 15%)

    # -----------------------------------------------------------------
    # Inter-individual variability (log-normal SDs on log-scale).
    # Omega column of Melhem 2018 Table 2 reports SD (paired CV% shown in
    # parentheses; the parameter itself is the SD, not the variance).
    # ini() variances are therefore SD^2.
    # -----------------------------------------------------------------
    etalfsc      ~ 0.440^2                                                                                                                                        # Melhem 2018 Table 2 (Omega FSC = 0.440, RSE 5%)
    etalksc      ~ 0.225^2                                                                                                                                        # Melhem 2018 Table 2 (Omega KSC = 0.225, RSE 5%)
    etalvd       ~ 0.282^2                                                                                                                                        # Melhem 2018 Table 2 (Omega VD = 0.282, RSE 7%)
    etalcld      ~ 0.370^2                                                                                                                                        # Melhem 2018 Table 2 (Omega CLD = 0.370, RSE 6%)
    etalkp       ~ 0.265^2                                                                                                                                        # Melhem 2018 Table 2 (Omega KP = 0.265, RSE 5%)
    etalkd       ~ 0.726^2                                                                                                                                        # Melhem 2018 Table 2 (Omega KD = 0.726, RSE 5%)
    etalstm1     ~ 0.315^2                                                                                                                                        # Melhem 2018 Table 2 (Omega STM1 = 0.315, RSE 5%)
    etalstm2     ~ 0.273^2                                                                                                                                        # Melhem 2018 Table 2 (Omega STM2 = 0.273, RSE 5%)
    etalkint     ~ 0.570^2                                                                                                                                        # Melhem 2018 Table 2 (Omega KINT = 0.570, RSE 5%)
    etalbsld     ~ 0.260^2                                                                                                                                        # Melhem 2018 Table 2 (Omega BSLD = 0.260, RSE 12%)
    etaltlag_cort ~ 0.698^2                                                                                                                                       # Melhem 2018 Table 2 (Omega LAG8 = 0.698, RSE 12%)
    etaltlag_chem ~ 0.110^2                                                                                                                                       # Melhem 2018 Table 2 (Omega LAG9 = 0.110, RSE 4%)
    etalcrt50    ~ 1.180^2                                                                                                                                        # Melhem 2018 Table 2 (Omega CRT50 = 1.180, RSE 15%)

    # Correlated IIV between KCHM and CHMSL (Melhem 2018 Table 2 row
    # 'Correlation (KCHM, CHMSL) = 0.731'). Covariance = r * sigma1 * sigma2.
    etalkel_chem + etalchmsl ~ c(0.259^2,
                                  0.99 * 0.731 * 0.259 * 2.28, 2.28^2)                                                                                            # Melhem 2018 Table 2 (Omega KCHM = 0.259, Omega CHMSL = 2.28, correlation 0.731; off-diagonal nudged by 0.99 to keep OMEGA positive definite under Cholesky sampling -- see references/known-vignette-failure-patterns.md section 1)

    # -----------------------------------------------------------------
    # Residual error -- Melhem 2018 statistical model: "Residual variability
    # was described using an additive error model in the log domain for both
    # PK and PD endpoints." This maps to lnorm() in nlmixr2 (log-normal
    # observation model, exp SD applied on the log scale).
    # -----------------------------------------------------------------
    expSd        <- 0.537                       ; label("Log-normal residual SD on G-CSF concentration Cc (PK)")                                                   # Melhem 2018 Table 2 (a1 PK = 0.537, RSE 1%)
    expSd_ANC    <- 0.298                       ; label("Log-normal residual SD on absolute neutrophil count ANC (PD)")                                            # Melhem 2018 Table 2 (a2 PD = 0.298, RSE 1%)
  })

  model({
    # -----------------------------------------------------------------
    # Individual parameters. FSC / KSC / VD / CLD / KD carry the
    # log-scale drug-type shift for pegfilgrastim (FORM_GCSF_PEG).
    # VD / CLD additionally carry the allometric exponent on WT / 70.
    # CLD carries the additional pegfilgrastim 300 ug/kg shift
    # (DOSE_PEG_300UGKG, only meaningful when FORM_GCSF_PEG = 1).
    # STM2 / KINT carry the healthy-volunteer shift (DIS_HEALTHY).
    # -----------------------------------------------------------------
    fsc      <- exp(lfsc + e_peg_fsc * FORM_GCSF_PEG + etalfsc)
    ksc      <- exp(lksc + e_peg_ksc * FORM_GCSF_PEG + etalksc)
    vd       <- exp(lvd  + e_peg_vd  * FORM_GCSF_PEG + etalvd)  * (WT / 70) ^ e_wt_vd
    cld      <- exp(lcld + e_peg_cld * FORM_GCSF_PEG + e_peg300_cld * DOSE_PEG_300UGKG + etalcld) * (WT / 70) ^ e_wt_cld
    kd       <- exp(lkd  + e_peg_kd  * FORM_GCSF_PEG + etalkd)

    kp       <- exp(lkp + etalkp)
    ktr      <- exp(lktr)
    kc       <- exp(lkc)
    sr       <- exp(lsr)

    stm1     <- exp(lstm1 + etalstm1)
    stm2     <- exp(lstm2 + e_hv_stm2 * DIS_HEALTHY + etalstm2)
    kint     <- exp(lkint + e_hv_kint * DIS_HEALTHY + etalkint)
    bsld     <- exp(lbsld + etalbsld)

    tlag_cort <- exp(ltlag_cort + etaltlag_cort)
    tlag_chem <- exp(ltlag_chem + etaltlag_chem)
    kel_cort <- exp(lkel_cort)
    kel_chem <- exp(lkel_chem + etalkel_chem)
    chmsl    <- exp(lchmsl + etalchmsl)
    crt50    <- exp(lcrt50 + etalcrt50)

    # -----------------------------------------------------------------
    # Drug clearance micro-constant. CLD acts on the FREE drug (fdc*vd)
    # only; bound drug (rdc*vd) is cleared via receptor internalisation
    # (KINT), following Melhem 2018 Methods 'G-CSF PK disposition model'
    # and supplement Equations.
    # -----------------------------------------------------------------
    ke_drug  <- cld / vd

    # -----------------------------------------------------------------
    # Quasi-equilibrium quadratic binding (Melhem 2018 supplement MLXTRAN
    # 'EQUATION' block). tdc = total drug concentration = central / vd;
    # trc = total receptor concentration = circ (RB compartment carries
    # the blood-neutrophil receptor pool). Free drug FDC is the positive
    # root of the binding quadratic:
    #   fdc = 0.5*(tdc - trc - kd + sqrt((tdc - trc - kd)^2 + 4*kd*tdc))
    # and bound drug RDC = tdc - fdc.
    # -----------------------------------------------------------------
    tdc      <- central / vd
    trc      <- circ
    qroot    <- tdc - trc - kd
    fdc      <- 0.5 * (qroot + sqrt(qroot * qroot + 4 * kd * tdc))
    rdc      <- tdc - fdc

    # -----------------------------------------------------------------
    # KPD-driven effect terms (Melhem 2018 supplement 'EQUATION'):
    #   DCRT = KCRT * CRT              (corticosteroid output rate)
    #   DCHM = KCHM * CHM              (chemotherapy output rate)
    #   CRTM = DCRT / (CRT50 + DCRT)   (Emax corticosteroid stimulation)
    #   STM  = CHMSL * DCHM            (linear chemotherapy elimination)
    # -----------------------------------------------------------------
    dcrt     <- kel_cort * depot_kpd_corticosteroid
    dchm     <- kel_chem * depot_kpd_chemotherapy
    crtm     <- dcrt / (crt50 + dcrt)
    stm_chem <- chmsl * dchm

    # Bound-receptor stimulation driver (rdc/trc) modulated by
    # corticosteroid effect. When crtm = 0 (no corticosteroid), the
    # expressions reduce to ST1 = 1 + STM1 * (RDC/TRC) and
    # ST2 = 1 + STM2 * (RDC/TRC); when crtm = 1 (full corticosteroid
    # stimulation), ST1 = 1 + STM1 and ST2 = 1 + STM2.
    st1      <- 1 + stm1 * ((1 - crtm) * (rdc / trc) + crtm)
    st2      <- 1 + stm2 * ((1 - crtm) * (rdc / trc) + crtm)

    # -----------------------------------------------------------------
    # ODE system (Melhem 2018 supplement 'EQUATION' block).
    #   depot                     : SC drug depot (nmol)
    #   central                   : total drug DT in the body (nmol)
    #   precursor1                : SM, stem-cell receptors (nM)
    #   precursor2                : MT, mitotic-cell receptors (nM; chemo acts here)
    #   precursor3                : PM1, first post-mitotic (nM)
    #   precursor4                : PM2, second post-mitotic (nM)
    #   circ                      : RB, blood-neutrophil receptors (nM)
    #   depot_kpd_corticosteroid  : CRT, KPD corticosteroid amount (mg)
    #   depot_kpd_chemotherapy    : CHM, KPD chemotherapy amount (mg)
    # -----------------------------------------------------------------
    d/dt(depot)                    <- -ksc * depot
    d/dt(central)                  <-  ksc * depot - ke_drug * (fdc * vd) - kint * (rdc * vd)
    d/dt(precursor1)               <-  kp * st1                        - ktr * st2 * precursor1
    d/dt(precursor2)               <-  ktr * st2 * precursor1          - (ktr * st2 + stm_chem) * precursor2
    d/dt(precursor3)               <-  ktr * st2 * precursor2          - ktr * st2 * precursor3
    d/dt(precursor4)               <-  ktr * st2 * precursor3          - ktr * st2 * precursor4
    d/dt(circ)                     <-  ktr * st2 * precursor4          - kc * circ
    d/dt(depot_kpd_corticosteroid) <- -kel_cort * depot_kpd_corticosteroid
    d/dt(depot_kpd_chemotherapy)   <- -kel_chem * depot_kpd_chemotherapy

    # -----------------------------------------------------------------
    # Steady-state initial conditions (Melhem 2018 supplement 'EQUATION'):
    #   SM(0) = MT(0) = PM1(0) = PM2(0) = KP / KTR
    #   RB(0)                           = KP / KC
    # so that with ST1 = ST2 = 1 (no drug binding, no corticosteroid) and
    # STM = 0 (no chemotherapy) each dS/dt = 0.
    # -----------------------------------------------------------------
    precursor1(0)               <- kp / ktr
    precursor2(0)               <- kp / ktr
    precursor3(0)               <- kp / ktr
    precursor4(0)               <- kp / ktr
    circ(0)                     <- kp / kc

    # -----------------------------------------------------------------
    # Bioavailability of the SC depot and KPD absorption lags.
    # -----------------------------------------------------------------
    f(depot)                        <- fsc
    alag(depot_kpd_corticosteroid)  <- tlag_cort
    alag(depot_kpd_chemotherapy)    <- tlag_chem

    # -----------------------------------------------------------------
    # Observations. G-CSF serum concentration Cc (nM) = FDC + BSLD
    # (free exogenous drug plus endogenous G-CSF baseline). Absolute
    # neutrophil count ANC (10^9 cells/L) = RB / SR. Because both
    # observables are computed *after* every ODE state (central, circ,
    # and the two KPD compartments), the vignette's event table uses
    # `cmt = "Cc"` / `cmt = "ANC"` with matching `dvid = 1L` / `dvid = 2L`
    # per the standard multi-output pattern in Harrold 2020 filgrastim;
    # the observable cmt slots land after the ODE states, so the
    # slot-renumbering bug documented in the skill's
    # known-vignette-failure-patterns section 2 does not fire.
    # -----------------------------------------------------------------
    Cc  <- fdc + bsld
    Cc  ~ lnorm(expSd)
    ANC <- circ / sr
    ANC ~ lnorm(expSd_ANC)
  })
}
