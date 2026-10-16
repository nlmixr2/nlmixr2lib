AsinPrieto_2021_hepatitisB_qsp <- function() {
  description <- paste(
    "QSP. Quantitative systems pharmacology model of acute hepatitis B",
    "virus (HBV) infection in adults across liver, plasma and lymph-node",
    "compartments: viral dynamics (hepatocyte infection, HBV and HBsAg",
    "synthesis, debris hepatocytes and ALT release), the innate response",
    "(dendritic cells, IFN-alpha from plasmacytoid DC, NK activation,",
    "IFN-gamma, TRAIL+ NK lysis, HBsAg inhibition of IFN-alpha), the",
    "cellular adaptive response (CTL priming in lymph node, CTL",
    "proliferation, IFN-gamma production and cytolysis in liver), the",
    "humoral response (B-cell activation, a five-stage plasmablast",
    "chain, plasma cells, anti-HBs/anti-HBc with antigen complexation)",
    "and the immunotolerant response (Treg-driven CTL exhaustion).",
    "36 ODE states. Infection starts from an HBV inoculum in the liver",
    "(HBV0 = 0.67 virion/mL, about 1000 copies in 1500 mL); with",
    "HBV0 = 0 the system holds at its uninfected baseline. All values are",
    "literature-derived, homeostasis-derived or calibrated constants;",
    "there is no IIV and no residual-error model. The memory-CTL branch",
    "(an output-only sink) is omitted and the short- and long-lived",
    "plasma-cell pools are merged because their split ratio and the",
    "long-lived death rate are not reported.",
    sep = " "
  )
  reference <- paste(
    "Asin-Prieto E, Parra-Guillen ZP, Gomez Mantilla JD, Vandenbossche J,",
    "Stuyckens K, Woot de Trixhe X, Perez-Ruixo JJ, Troconiz IF.",
    "A quantitative systems pharmacology model for acute viral hepatitis B.",
    "Comput Struct Biotechnol J. 2021;19:4997-5007.",
    "doi:10.1016/j.csbj.2021.08.052.",
    "Equations from Supplementary Appendix 1 (equations 1-47); parameter",
    "values from Table S1 and initial conditions from Table S2 of the",
    "supplementary tables.",
    sep = " "
  )
  vignette <- "AsinPrieto_2021_hepatitisB_qsp"

  units <- list(
    time = "day",
    dosing = "virion/mL",
    concentration = "IU/mL"
  )

  # Every state is a mechanistic species of the HBV QSP model; none maps
  # onto a canonical PK compartment. Suffixes name the anatomical
  # compartment: _LV = liver, _PL = plasma, _LN = lymph node.
  paper_specific_compartments <- c(
    # ---- liver (Appendix 1 eq 2-18) ----
    "iHep_LV",
    "dHep_LV",
    "HBV_LV",
    "HBsAg_LV",
    "DC_LV",
    "DCa_LV",
    "IFNa_LV",
    "NK_LV",
    "NKa_LV",
    "IFNg_LV",
    "CTLa_LV",
    "CTLexh_LV",
    "Th0_LV",
    "Treg_LV",
    "ALT_LV",
    # ---- plasma (eq 19-34) ----
    "ALT_PL",
    "HBV_PL",
    "HBsAg_PL",
    "DC_PL",
    "NK_PL",
    "CTLa_PL",
    "PC_PL",
    "antiHBs_PL",
    "antiHBc_PL",
    "cpxHBV_PL",
    "cpxHBsAg_PL",
    # ---- lymph node (eq 35-47) ----
    "DCa_LN",
    "CTL_LN",
    "CTLa_LN",
    "Bcell_LN",
    "PB1_LN",
    "PB2_LN",
    "PB3_LN",
    "PB4_LN",
    "PB5_LN",
    "PC_LN"
  )

  # What each ODE state holds, in what units, in what matrix. Transcribed
  # from Table S2 (Model Species: compartment, definition, units).
  compartmentData <- list(
    iHep_LV = list(analyte = "infected hepatocytes", units = "cell/mL", specimen = "tissue", verified = TRUE),
    dHep_LV = list(analyte = "debris (dead) hepatocytes", units = "cell/mL", specimen = "tissue", verified = TRUE),
    HBV_LV = list(analyte = "hepatitis B virions", units = "virion/mL", specimen = "tissue", verified = TRUE),
    HBsAg_LV = list(
      analyte = "hepatitis B surface antigen",
      units = "molecule/mL",
      specimen = "tissue",
      verified = TRUE
    ),
    DC_LV = list(analyte = "naive dendritic cells", units = "cell/mL", specimen = "tissue", verified = TRUE),
    DCa_LV = list(analyte = "activated dendritic cells", units = "cell/mL", specimen = "tissue", verified = TRUE),
    IFNa_LV = list(analyte = "interferon-alpha", units = "pg/mL", specimen = "tissue", verified = TRUE),
    NK_LV = list(analyte = "naive natural killer cells", units = "cell/mL", specimen = "tissue", verified = TRUE),
    NKa_LV = list(analyte = "activated natural killer cells", units = "cell/mL", specimen = "tissue", verified = TRUE),
    IFNg_LV = list(analyte = "interferon-gamma", units = "pg/mL", specimen = "tissue", verified = TRUE),
    CTLa_LV = list(
      analyte = "HBV-specific activated cytotoxic T lymphocytes",
      units = "cell/mL",
      specimen = "tissue",
      verified = TRUE
    ),
    CTLexh_LV = list(
      analyte = "exhausted cytotoxic T lymphocytes",
      units = "cell/mL",
      specimen = "tissue",
      verified = TRUE
    ),
    Th0_LV = list(analyte = "regulatory T-cell precursors", units = "cell/mL", specimen = "tissue", verified = TRUE),
    Treg_LV = list(analyte = "regulatory T cells", units = "cell/mL", specimen = "tissue", verified = TRUE),
    ALT_LV = list(analyte = "alanine aminotransferase", units = "IU/mL", specimen = "tissue", verified = TRUE),
    ALT_PL = list(analyte = "alanine aminotransferase", units = "IU/mL", specimen = "plasma", verified = TRUE),
    HBV_PL = list(analyte = "hepatitis B virions", units = "virion/mL", specimen = "plasma", verified = TRUE),
    HBsAg_PL = list(
      analyte = "hepatitis B surface antigen",
      units = "molecule/mL",
      specimen = "plasma",
      verified = TRUE
    ),
    DC_PL = list(analyte = "naive dendritic cells", units = "cell/mL", specimen = "plasma", verified = TRUE),
    NK_PL = list(analyte = "naive natural killer cells", units = "cell/mL", specimen = "plasma", verified = TRUE),
    CTLa_PL = list(
      analyte = "HBV-specific activated cytotoxic T lymphocytes",
      units = "cell/mL",
      specimen = "plasma",
      verified = TRUE
    ),
    PC_PL = list(
      analyte = "plasma cells (short-lived plus long-lived)",
      units = "cell/mL",
      specimen = "plasma",
      verified = TRUE
    ),
    antiHBs_PL = list(analyte = "anti-HBs antibody", units = "molecule/mL", specimen = "plasma", verified = TRUE),
    antiHBc_PL = list(analyte = "anti-HBc antibody", units = "molecule/mL", specimen = "plasma", verified = TRUE),
    cpxHBV_PL = list(analyte = "HBV:anti-HBs complex", units = "molecule/mL", specimen = "plasma", verified = TRUE),
    cpxHBsAg_PL = list(analyte = "HBsAg:anti-HBs complex", units = "molecule/mL", specimen = "plasma", verified = TRUE),
    DCa_LN = list(analyte = "activated dendritic cells", units = "cell/mL", specimen = "lymph", verified = TRUE),
    CTL_LN = list(
      analyte = "cytotoxic T-lymphocyte precursors",
      units = "cell/mL",
      specimen = "lymph",
      verified = TRUE
    ),
    CTLa_LN = list(
      analyte = "HBV-specific activated cytotoxic T lymphocytes",
      units = "cell/mL",
      specimen = "lymph",
      verified = TRUE
    ),
    Bcell_LN = list(analyte = "naive B cells", units = "cell/mL", specimen = "lymph", verified = TRUE),
    PB1_LN = list(
      analyte = "plasmablasts (maturation stage 1)",
      units = "cell/mL",
      specimen = "lymph",
      verified = TRUE
    ),
    PB2_LN = list(
      analyte = "plasmablasts (maturation stage 2)",
      units = "cell/mL",
      specimen = "lymph",
      verified = TRUE
    ),
    PB3_LN = list(
      analyte = "plasmablasts (maturation stage 3)",
      units = "cell/mL",
      specimen = "lymph",
      verified = TRUE
    ),
    PB4_LN = list(
      analyte = "plasmablasts (maturation stage 4)",
      units = "cell/mL",
      specimen = "lymph",
      verified = TRUE
    ),
    PB5_LN = list(
      analyte = "plasmablasts (maturation stage 5)",
      units = "cell/mL",
      specimen = "lymph",
      verified = TRUE
    ),
    PC_LN = list(
      analyte = "plasma cells (short-lived plus long-lived)",
      units = "cell/mL",
      specimen = "lymph",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 49,
    n_studies = 4,
    age_range = "adults; not reported by the source cohorts",
    disease_state = paste(
      "Acute, self-resolving hepatitis B virus infection in",
      "immunocompetent adults (horizontal transmission)."
    ),
    dose_range = paste(
      "No drug. Infection is initiated by an HBV inoculum of 0.67",
      "virion/mL in the liver (Table S2), about 1000 HBV DNA copies in",
      "1500 mL, the minimum needle-stick inoculum cited in the Discussion."
    ),
    notes = paste(
      "The model was not fitted to patient-level data. Parameters come",
      "from the literature, from homeostasis derivations or from manual",
      "calibration (Table S1); only the two ALT parameters were",
      "calibrated against clinical ALT data. Typical-value predictions",
      "were compared with digitised acute-HBV data from 4 published",
      "cohorts (Table 1): Webster 2000 (n = 5), Dunn 2009 (n = 21",
      "sampled), Fisicaro 2009 (n = 2) and Chulanov 2003 (n = 21);",
      "n_subjects sums those cohort sizes. The chronicity analysis",
      "(500 virtual subjects, 30% log-normal variability on 11",
      "parameters) is an exploratory simulation, not a fitted",
      "variability model, and is not encoded as IIV."
    )
  )

  ini({
    # =====================================================================
    # All values are literature-derived, homeostasis-derived or calibrated
    # constants from Table S1 (Definition, estimates, ranges and sources of
    # the different model parameters), not maximum-likelihood estimates.
    # They are exposed as typical values without fixed() so downstream
    # users can recalibrate (same convention as
    # BoivinChampeaux_2026_hepatitisB_qsp).
    # Time unit is days throughout.
    # =====================================================================

    # ---- Volumes (Table S1) ---------------------------------------------
    V_LV <- 1500 ; label("Liver volume (mL)")                                                   # Table S1 V_LV
    V_PL <- 3000 ; label("Plasma volume (mL)")                                                  # Table S1 V_PL
    V_LN <- 730  ; label("Lymph nodes volume (mL)")                                             # Table S1 V_LN

    # ---- Inoculum (Table S2 HBV_LV initial condition) --------------------
    HBV0 <- 0.67 ; label("Initial liver HBV inoculum (virion/mL); 0 = uninfected baseline")    # Table S2 HBV_LV = 0.67 virions/mL

    # ---- Synthesis, proliferation, recruitment & activation (Table S1) ---
    kinf          <- 3e-10     ; label("Hepatocyte infection rate constant (mL/virion/day)")               # Table S1 k_inf
    ksyn_HBV      <- 315       ; label("HBV synthesis rate constant by infected hepatocytes (virion/cell/day)") # Table S1 k_syn_HBV
    ksyn_HBsAg    <- 45165     ; label("HBsAg synthesis rate constant by infected hepatocytes (molecule/cell/day)") # Table S1 k_syn_HBsAg
    kact_DC       <- 0.2       ; label("Naive DC activation rate constant (1/day)")                         # Table S1 k_act_DC
    ksyn_IFNa     <- 459.2     ; label("IFN-alpha independent production (pg/mL/day)")                     # Table S1 k_syn_IFNa
    ksyn_IFNa_DC  <- 0.1       ; label("IFN-alpha produced by the pDC fraction of activated DC (pg/cell/day)") # Table S1 k_syn_IFNa_DC*
    kact_NK       <- 1.82e-5   ; label("Naive NK activation rate constant (1/day)")                         # Table S1 k_act_NK
    ksyn_IFNg_NK  <- 0.00334   ; label("IFN-gamma production by activated NK (pg/cell/day)")                # Table S1 k_syn_IFNg_NK
    ksyn_IFNg_CTL <- 0.0034    ; label("IFN-gamma production by liver IFN-gamma+ CTL (pg/cell/day)")        # Table S1 k_syn_IFNg_CTL
    kprol_CTL     <- 0.57      ; label("Proliferation of activated CTL in liver and lymph node (1/day)")    # Table S1 k_prol_CTL
    kact_Treg     <- 0.0256    ; label("Activation of Treg in the liver (mL/cell/day)")                     # Table S1 k_act_Treg
    kprol_Treg    <- 0.7       ; label("Proliferation of Treg in the liver (1/day)")                        # Table S1 k_prol_Treg
    ksyn_ALT      <- 8.1e-8    ; label("ALT synthesis by debris hepatocytes (IU/cell/day)")                 # Table S1 k_syn_ALT
    ksyn_DC       <- 1.076e4   ; label("Naive DC daily production (cell/mL/day)")                           # Table S1 k_syn_DC
    ksyn_NK       <- 1.185e5   ; label("Naive NK daily production (cell/mL/day)")                           # Table S1 k_syn_NK
    kprod_antiHBs <- 1.2e9     ; label("Anti-HBs production per plasma cell (molecule/cell/day)")           # Table S1 k_prod_antiHBs
    kprod_antiHBc <- 1.2e9     ; label("Anti-HBc production per plasma cell (molecule/cell/day)")           # Table S1 k_prod_antiHBc
    ksyn_CTL      <- 1000      ; label("CTL production in the lymph node (cell/mL/day)")                    # Table S1 k_syn_CTL
    krec_CTL      <- 0.1       ; label("CTL recruitment by activated DC in the lymph node (cell/cell/day)") # Table S1 k_rec_CTL
    kact_CTL      <- 1e-4      ; label("CTL activation by activated DC in the lymph node (mL/cell/day)")    # Table S1 k_act_CTL
    ksyn_Bcell    <- 200       ; label("Naive B-cell daily production (cell/mL/day)")                       # Table S1 k_syn_Bcell
    kact_Bcell    <- 1         ; label("Activation of B cells (1/day)")                                     # Table S1 k_act_Bcell

    # ---- Degradation / disappearance (Table S1) --------------------------
    kdeath_Hep    <- 0.0039    ; label("Natural death of healthy and infected hepatocytes (1/day)")         # Table S1 k_death_Hep
    kdeath_dHep   <- 13        ; label("Disappearance of debris hepatocytes (1/day)")                       # Table S1 k_death_dHep
    kdeg_HBV      <- 0.7       ; label("Natural degradation of HBV (1/day)")                                # Table S1 k_deg_HBV
    kdeath_DC     <- 0.231     ; label("Natural death of activated and naive DC (1/day)")                   # Table S1 k_death_DC
    kdeg_IFNa     <- 5.74      ; label("Degradation of IFN-alpha (1/day)")                                  # Table S1 k_deg_IFNa
    kdeath_NK     <- 0.069     ; label("Natural death of activated and naive NK (1/day)")                   # Table S1 k_death_NK
    kdeg_IFNg     <- 33.4      ; label("Degradation of IFN-gamma (1/day)")                                  # Table S1 k_deg_IFNg
    kexh_CTLa     <- 0.403     ; label("Exhaustion rate constant of activated CTL (1/day)")                 # Table S1 k_ex_CTL*
    kdeath_CTLa   <- 0.33      ; label("Natural death of activated CTL (1/day)")                            # Table S1 k_death_CTL*
    kdeath_CTLexh <- 0.033     ; label("Natural death of exhausted CTL (1/day)")                            # Table S1 k_death_CTLexh
    kdeath_Treg   <- 0.1       ; label("Natural death of Treg (1/day)")                                     # Table S1 k_death_Treg
    kdeg_ALT      <- 0.09      ; label("Degradation of ALT (1/day)")                                        # Table S1 k_deg_ALT
    kdeg_HBsAg    <- 0.0835    ; label("Natural degradation of HBsAg (1/day)")                              # Table S1 k_deg_HBsAg
    kdeg_antiHBs  <- 0.033     ; label("Natural degradation of anti-HBs (1/day)")                           # Table S1 k_deg_antiHBs
    kdeg_cpxHBV   <- 2.7       ; label("Natural degradation of the HBV:anti-HBs complex (1/day)")           # Table S1 k_deg_complexHBV
    kdeg_cpxHBsAg <- 2.7       ; label("Natural degradation of the HBsAg:anti-HBs complex (1/day)")         # Table S1 k_deg_complexantiHBs
    kdeg_antiHBc  <- 0.033     ; label("Natural degradation of anti-HBc (1/day)")                           # Table S1 k_deg_antiHBc
    kdeath_CTL    <- 0.102     ; label("Natural death of CTL precursors in the lymph node (1/day)")         # Table S1 k_death_CTL
    kdeath_Bcell  <- 0.2       ; label("Natural death of naive B cells (1/day)")                            # Table S1 k_death_Bcell
    MTT_PB        <- 50        ; label("Mean transit time of plasmablasts (day)")                           # Table S1 MTT_PB
    N_PB          <- 5         ; label("Number of plasmablast transit compartments (unitless)")             # Table S1 N_PB
    kdeath_PB     <- 0.1       ; label("Natural death of plasmablasts (1/day)")                             # Table S1 k_death_PB
    kdeath_PC     <- 0.1       ; label("Natural death of plasma cells in plasma and lymph node (1/day)")    # Table S1 k_death_PC

    # ---- Maximum concentrations and ratios (Table S1) --------------------
    Hep_tot       <- 1.33e8    ; label("Total liver hepatocytes (cell/mL)")                                 # Table S1 Hep_tot
    pDC_ratio     <- 0.36      ; label("Fraction of pDC among activated DC (unitless)")                     # Table S1 pDC_RATIO
    CTLa_ratio    <- 0.65      ; label("Fraction of activated CTL that are IFN-gamma+ (unitless)")          # Table S1 CTL*_RATIO
    NK_ratio      <- 0.5       ; label("Fraction of activated NK that are TRAIL+ (unitless)")               # Table S1 NK_RATIO
    CTLa_max_LV   <- 2e6       ; label("Maximal concentration of activated CTL in the liver (cell/mL)")     # Table S1 CTL*_max_LV
    Treg_max_LV   <- 2e6       ; label("Maximal concentration of Treg in the liver (cell/mL)")              # Table S1 Treg_max_LV
    CTLa_max_LN   <- 1e5       ; label("Maximal concentration of activated CTL in the lymph node (cell/mL)") # Table S1 CTL*_max_LN

    # ---- Distribution and binding (Table S1) -----------------------------
    kLV_PL_HBV    <- 1         ; label("HBV distribution from liver to plasma (1/day)")                     # Table S1 k_LV_PL_HBV
    kPL_LV_HBV    <- 1         ; label("HBV distribution from plasma to liver (1/day)")                     # Table S1 k_PL_LV_HBV
    kLV_PL_HBsAg  <- 1         ; label("HBsAg distribution from liver to plasma (1/day)")                   # Table S1 k_LV_PL_HBsAg
    kPL_LV_HBsAg  <- 1         ; label("HBsAg distribution from plasma to liver (1/day)")                   # Table S1 k_PL_LV_HBsAg
    kPL_LV_DC     <- 0.0599    ; label("DC distribution from plasma to liver (1/day)")                      # Table S1 k_PL_LV_DC
    kLV_LN_DCa    <- 0.9       ; label("Activated DC distribution from liver to lymph node (1/day)")        # Table S1 k_LV_LN_DC*
    kPL_LV_NK     <- 0.575     ; label("NK distribution from plasma to liver (1/day)")                      # Table S1 k_PL_LV_NK
    kPL_LV_CTLa   <- 0.9       ; label("Activated CTL distribution from plasma to liver (1/day)")           # Table S1 k_PL_LV_CTL*
    kLV_PL_ALT    <- 0.09      ; label("ALT distribution from liver to plasma (1/day)")                     # Table S1 k_LV_PL_ALT
    kLN_PL_CTLa   <- 0.9       ; label("Activated CTL distribution from lymph node to plasma (1/day)")      # Table S1 k_LN_PL_CTL*
    kon_antiHBs   <- 1e-12     ; label("Association rate of HBV and HBsAg with anti-HBs (mL/molecule/day)") # Table S1 k_on_antiHBs
    koff_antiHBs  <- 10        ; label("Dissociation rate of anti-HBs complexes (1/day)")                   # Table S1 k_off_antiHBs
    kLN_PL_PC     <- 0.9       ; label("Plasma-cell distribution from lymph node to plasma (1/day)")        # Table S1 k_LN_PL_PC

    # ---- Activity constants (Table S1) -----------------------------------
    klytic_CTL      <- 3.7e-5    ; label("Cytopathic effect of CTL on infected hepatocytes (mL/cell/day)")  # Table S1 k_lytic_CTL
    klytic_NK       <- 7.536e-7  ; label("Cytopathic effect of NK on infected hepatocytes (mL/cell/day)")   # Table S1 k_lytic_NK
    IFNa50_HBV      <- 100       ; label("IFN-alpha inhibiting 50% of HBV synthesis (pg/mL)")               # Table S1 IFNa_50_HBV
    IFNg50_HBV      <- 531.7     ; label("IFN-gamma inhibiting 50% of HBV synthesis (pg/mL)")               # Table S1 IFNg_50_HBV
    HBV50_DC        <- 10        ; label("HBV triggering 50% of maximal DC activation (virion/mL)")         # Table S1 HBV_50_DC
    HBsAg50_IFNa    <- 2.2e11    ; label("Liver HBsAg inhibiting 50% of DC-induced IFN-alpha synthesis (molecule/mL)") # Table S1 HBsAg_LV_50_IFNa
    SLP_IFNa_NK     <- 5.25      ; label("Fractional change in NK activation per unit IFN-alpha (mL/pg)")   # Table S1 SLP_IFNa_NK
    HBV50_CTL_prol  <- 10        ; label("HBV triggering 50% of maximal CTL proliferation (virion/mL)")     # Table S1 HBV_50_CTL_prol
    gamma_CTL       <- 4         ; label("Shape parameter of HBV effect on CTL proliferation (unitless)")   # Table S1 gamma_CTL
    Treg50_CTL_exh  <- 1e6       ; label("Treg triggering 50% of maximal CTL exhaustion (cell/mL)")         # Table S1 Treg_50_CTL_exh
    CTL50_Treg_prol <- 1e6       ; label("Activated CTL triggering 50% of maximal Treg proliferation (cell/mL)") # Table S1 CTL_50_Treg_prol
    DC50_CTL_prol   <- 1         ; label("Activated DC triggering 50% of maximal CTL proliferation in lymph node (cell/mL)") # Table S1 DC_50_CTL_prol
    DC50_Bcell_prol <- 1         ; label("Activated DC triggering 50% of maximal B-cell activation in lymph node (cell/mL)") # Table S1 DC_50_Bcell_prol
  })

  model({
    # ---- Algebraic species (Appendix 1 eq 1, 11, 13) ---------------------
    Hep_LV <- Hep_tot - iHep_LV                                                   # eq 1
    NKtrail_LV <- NK_ratio * NKa_LV                                               # eq 11
    CTLifng_LV <- CTLa_ratio * CTLa_LV                                            # eq 13

    # Non-cytolytic inhibition of HBV and HBsAg synthesis by liver IFN-alpha
    # and IFN-gamma (eq 4, 5).
    inh_IFN <- (1 - IFNa_LV / (IFNa_LV + IFNa50_HBV)) *
      (1 - IFNg_LV / (IFNg_LV + IFNg50_HBV))

    # Main text section 2.2: viral replication is negligible while fewer
    # than 1 hepatocyte is infected. Read literally as a whole-liver cell
    # count (iHep_LV * V_LV < 1 cell) gating HBV synthesis; not latched.
    repl_HBV <- 1
    if (iHep_LV * V_LV < 1) {
      repl_HBV <- 0
    }

    # Hill enablers (exponent 4 as printed in eq 6, 7, 14, 17, 37, 40, 41)
    act_DC_HBV <- HBV_LV^4 / (HBV_LV^4 + HBV50_DC^4)                              # eq 6, 7
    prol_CTL_HBV <- HBV_LV^gamma_CTL / (HBV_LV^gamma_CTL + HBV50_CTL_prol^gamma_CTL) # eq 14
    exh_Treg <- Treg_LV^4 / (Treg_LV^4 + Treg50_CTL_exh^4)                        # eq 14, 15
    prol_Treg_CTL <- CTLa_LV^4 / (CTLa_LV^4 + CTL50_Treg_prol^4)                  # eq 17
    prol_CTL_DC <- DCa_LN^4 / (DCa_LN^4 + DC50_CTL_prol^4)                        # eq 37
    act_Bcell_DC <- DCa_LN^4 / (DCa_LN^4 + DC50_Bcell_prol^4)                     # eq 40, 41
    ktr_PB <- N_PB / MTT_PB                                                       # eq 41-47

    # ---- Liver (eq 2-18) ----------------------------------------------------
    d/dt(iHep_LV) <- kinf * HBV_LV * Hep_LV - kdeath_Hep * iHep_LV -
      klytic_CTL * CTLa_LV * iHep_LV - klytic_NK * NKtrail_LV * iHep_LV          # eq 2
    d/dt(dHep_LV) <- kdeath_Hep * (Hep_LV + iHep_LV) +
      klytic_CTL * CTLa_LV * iHep_LV + klytic_NK * NKtrail_LV * iHep_LV -
      kdeath_dHep * dHep_LV                                                       # eq 3
    d/dt(HBV_LV) <- ksyn_HBV * repl_HBV * iHep_LV * inh_IFN -
      kinf * HBV_LV * Hep_LV - kLV_PL_HBV * HBV_LV +
      kPL_LV_HBV * HBV_PL * V_PL / V_LV - kdeg_HBV * HBV_LV                       # eq 4
    d/dt(HBsAg_LV) <- ksyn_HBsAg * iHep_LV * inh_IFN - kLV_PL_HBsAg * HBsAg_LV +
      kPL_LV_HBsAg * HBsAg_PL * V_PL / V_LV - kdeg_HBsAg * HBsAg_LV              # eq 5
    d/dt(DC_LV) <- kPL_LV_DC * DC_PL * V_PL / V_LV - kact_DC * act_DC_HBV * DC_LV -
      kdeath_DC * DC_LV                                                           # eq 6
    d/dt(DCa_LV) <- kact_DC * act_DC_HBV * DC_LV - kdeath_DC * DCa_LV -
      kLV_LN_DCa * DCa_LV                                                         # eq 7
    # eq 8: HBsAg inhibits the whole IFN-alpha synthesis (basal + pDC).
    d/dt(IFNa_LV) <- (ksyn_IFNa + ksyn_IFNa_DC * DCa_LV * pDC_ratio) *
      (1 - HBsAg_LV / (HBsAg_LV + HBsAg50_IFNa)) - kdeg_IFNa * IFNa_LV            # eq 8
    d/dt(NK_LV) <- kPL_LV_NK * NK_PL * V_PL / V_LV -
      kact_NK * NK_LV * (1 + SLP_IFNa_NK * IFNa_LV) - kdeath_NK * NK_LV           # eq 9
    d/dt(NKa_LV) <- kact_NK * NK_LV * (1 + SLP_IFNa_NK * IFNa_LV) -
      kdeath_NK * NKa_LV                                                          # eq 10
    # eq 12: basal IFN-gamma production is carried by ksyn_IFNg_NK * NKa_LV.
    d/dt(IFNg_LV) <- ksyn_IFNg_NK * NKa_LV + ksyn_IFNg_CTL * CTLifng_LV -
      kdeg_IFNg * IFNg_LV                                                         # eq 12
    # eq 14: liver CTL* death uses kdeath_CTLa (Table S1 k_death_CTL*).
    d/dt(CTLa_LV) <- kPL_LV_CTLa * CTLa_PL * V_PL / V_LV +
      kprol_CTL * CTLa_LV * (1 - (CTLa_LV + CTLexh_LV) / CTLa_max_LV) * prol_CTL_HBV -
      kexh_CTLa * CTLa_LV * exh_Treg - kdeath_CTLa * CTLa_LV                      # eq 14
    d/dt(CTLexh_LV) <- kexh_CTLa * CTLa_LV * exh_Treg - kdeath_CTLexh * CTLexh_LV # eq 15
    d/dt(Th0_LV) <- -kact_Treg * Th0_LV * DCa_LV                                  # eq 16
    d/dt(Treg_LV) <- kact_Treg * Th0_LV * DCa_LV +
      kprol_Treg * Treg_LV * (1 - Treg_LV / Treg_max_LV) * prol_Treg_CTL -
      kdeath_Treg * Treg_LV                                                       # eq 17
    d/dt(ALT_LV) <- ksyn_ALT * dHep_LV - kLV_PL_ALT * ALT_LV - kdeg_ALT * ALT_LV  # eq 18

    # ---- Plasma (eq 19-34) -------------------------------------------------
    d/dt(ALT_PL) <- kLV_PL_ALT * ALT_LV * V_LV / V_PL - kdeg_ALT * ALT_PL          # eq 19
    d/dt(HBV_PL) <- kLV_PL_HBV * HBV_LV * V_LV / V_PL - kPL_LV_HBV * HBV_PL -
      kdeg_HBV * HBV_PL - kon_antiHBs * HBV_PL * antiHBs_PL +
      koff_antiHBs * cpxHBV_PL                                                    # eq 20
    d/dt(HBsAg_PL) <- kLV_PL_HBsAg * HBsAg_LV * V_LV / V_PL -
      kPL_LV_HBsAg * HBsAg_PL - kdeg_HBsAg * HBsAg_PL -
      kon_antiHBs * HBsAg_PL * antiHBs_PL + koff_antiHBs * cpxHBsAg_PL            # eq 22
    d/dt(DC_PL) <- ksyn_DC - kPL_LV_DC * DC_PL - kdeath_DC * DC_PL                # eq 23
    # eq 25: the death term acts on NK_PL (footnote 10 homeostasis).
    d/dt(NK_PL) <- ksyn_NK - kPL_LV_NK * NK_PL - kdeath_NK * NK_PL                # eq 25
    d/dt(CTLa_PL) <- kLN_PL_CTLa * CTLa_LN * V_LN / V_PL - kPL_LV_CTLa * CTLa_PL -
      kdeath_CTLa * CTLa_PL                                                       # eq 27
    # eq 29 + eq 30 summed: PC = SPC + LPC with k_death_LPC = k_death_PC.
    d/dt(PC_PL) <- kLN_PL_PC * PC_LN * V_LN / V_PL - kdeath_PC * PC_PL            # eq 29, 30
    d/dt(antiHBs_PL) <- kprod_antiHBs * PC_PL - kdeg_antiHBs * antiHBs_PL -
      kon_antiHBs * (HBV_PL + HBsAg_PL) * antiHBs_PL +
      koff_antiHBs * (cpxHBV_PL + cpxHBsAg_PL)                                    # eq 31
    d/dt(antiHBc_PL) <- kprod_antiHBc * PC_PL - kdeg_antiHBc * antiHBc_PL         # eq 32
    d/dt(cpxHBV_PL) <- kon_antiHBs * HBV_PL * antiHBs_PL - koff_antiHBs * cpxHBV_PL -
      kdeg_cpxHBV * cpxHBV_PL                                                     # eq 33
    d/dt(cpxHBsAg_PL) <- kon_antiHBs * HBsAg_PL * antiHBs_PL -
      koff_antiHBs * cpxHBsAg_PL - kdeg_cpxHBsAg * cpxHBsAg_PL                    # eq 34

    # ---- Lymph node (eq 35-47) ---------------------------------------------
    d/dt(DCa_LN) <- kLV_LN_DCa * DCa_LV * V_LV / V_LN - kdeath_DC * DCa_LN         # eq 35
    d/dt(CTL_LN) <- ksyn_CTL + krec_CTL * DCa_LN - kdeath_CTL * CTL_LN -
      kact_CTL * CTL_LN * DCa_LN                                                  # eq 36
    # eq 37: lymph-node efflux acts on CTLa_LN (mass balance with eq 27).
    d/dt(CTLa_LN) <- kact_CTL * CTL_LN * DCa_LN +
      kprol_CTL * CTLa_LN * (1 - CTLa_LN / CTLa_max_LN) * prol_CTL_DC -
      kLN_PL_CTLa * CTLa_LN - kdeath_CTLa * CTLa_LN                               # eq 37
    d/dt(Bcell_LN) <- ksyn_Bcell - kact_Bcell * Bcell_LN * act_Bcell_DC -
      kdeath_Bcell * Bcell_LN                                                     # eq 40
    d/dt(PB1_LN) <- kact_Bcell * Bcell_LN * act_Bcell_DC - ktr_PB * PB1_LN -
      kdeath_PB * PB1_LN                                                          # eq 41
    d/dt(PB2_LN) <- ktr_PB * (2 * PB1_LN - PB2_LN) - kdeath_PB * PB2_LN          # eq 42
    d/dt(PB3_LN) <- ktr_PB * (2 * PB2_LN - PB3_LN) - kdeath_PB * PB3_LN          # eq 43
    d/dt(PB4_LN) <- ktr_PB * (2 * PB3_LN - PB4_LN) - kdeath_PB * PB4_LN          # eq 44
    d/dt(PB5_LN) <- ktr_PB * (2 * PB4_LN - PB5_LN) - kdeath_PB * PB5_LN          # eq 45
    # eq 46 + eq 47 summed: the PC_RATIO split cancels.
    d/dt(PC_LN) <- 2 * ktr_PB * PB5_LN - kLN_PL_PC * PC_LN - kdeath_PC * PC_LN    # eq 46, 47

    # ---- Reported outputs ------------------------------------------------
    HBVDNA <- HBV_PL / 5                                                          # eq 21; 0.2 IU/copy (Table S2 footnote 33)
    HBsAg_ngmL <- HBsAg_PL / 3.011e11 * 1000                                      # 3.011e11 molecule/ug (Table S1 footnote 2)
    ALT_UL <- ALT_PL * 1000                                                       # IU/mL to U/L
    IFNa_PL <- IFNa_LV * V_LV / V_PL                                              # eq 24
    IFNg_PL <- IFNg_LV * V_LV / V_PL                                              # eq 26

    # ---- Initial conditions (Table S2) -------------------------------------
    iHep_LV(0) <- 0
    dHep_LV(0) <- 40000
    HBV_LV(0) <- HBV0
    HBsAg_LV(0) <- 0
    DC_LV(0) <- 19200
    DCa_LV(0) <- 0
    IFNa_LV(0) <- 80
    NK_LV(0) <- 2.7e6
    NKa_LV(0) <- 3e5
    IFNg_LV(0) <- 30
    CTLa_LV(0) <- 0
    CTLexh_LV(0) <- 0
    Th0_LV(0) <- 2e6
    Treg_LV(0) <- 0
    ALT_LV(0) <- 0.018
    ALT_PL(0) <- 0.009
    HBV_PL(0) <- 0
    HBsAg_PL(0) <- 0
    DC_PL(0) <- 37000
    NK_PL(0) <- 1.84e5
    CTLa_PL(0) <- 0
    PC_PL(0) <- 0
    antiHBs_PL(0) <- 0
    antiHBc_PL(0) <- 0
    cpxHBV_PL(0) <- 0
    cpxHBsAg_PL(0) <- 0
    DCa_LN(0) <- 0
    # Table S2 footnote 35 defines CTL_LN(0) = k_syn_CTL / k_death_CTL
    # (9804 cell/mL); the table prints it rounded to 1e4.
    CTL_LN(0) <- ksyn_CTL / kdeath_CTL
    CTLa_LN(0) <- 0
    Bcell_LN(0) <- 1000
    PB1_LN(0) <- 0
    PB2_LN(0) <- 0
    PB3_LN(0) <- 0
    PB4_LN(0) <- 0
    PB5_LN(0) <- 0
    PC_LN(0) <- 0
  })
}
