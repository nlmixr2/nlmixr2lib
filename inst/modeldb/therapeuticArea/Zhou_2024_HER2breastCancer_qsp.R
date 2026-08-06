Zhou_2024_HER2breastCancer_qsp <- function() {
  description <- paste(
    "QSP. Quantitative systems pharmacology model of HER2-positive",
    "metastatic breast cancer coupling the ErbB receptor network to",
    "PI3K/AKT, Ras/MAPK and BTK signal transduction and to tumour",
    "growth, with five clinically approved therapeutics built in:",
    "the antibody-drug conjugates T-DM1 and T-DXd, the tyrosine",
    "kinase inhibitors lapatinib and pyrotinib, and capecitabine",
    "(with its 5'DFCR / 5'DFUR / 5-FU metabolite cascade).",
    "54 ODE states and 79 reactions: 18 PK states (ADC, naked",
    "antibody, released payload, drug-antibody ratio, and the two",
    "TKI and capecitabine disposition chains) plus 36 tumour-cell",
    "states (EGFR / HER2 / HER3 / HER4 monomers, EGF- and",
    "NRG1-driven homo- and heterodimers and their phosphorylated",
    "forms, Raf / ERK / PI3K / AKT / BTK phospho-cycles, an",
    "AKT-driven negative feedback on HER3 synthesis, receptor-",
    "mediated ADC binding, internalisation and payload release, and",
    "a logistic tumour-growth state). Default parameterisation is",
    "the in vivo SKBR3 mouse xenograft with T-DM1; in vitro and",
    "T-DXd, BT-474, NCI-N87 and ZR-75-1 alternates are given in the",
    "per-parameter comments and in population$notes.",
    "All parameters are literature constants or the authors'",
    "pattern-search calibration estimates; the paper reports no",
    "standard errors, no inter-individual variability and no",
    "residual-error model, so the model is deterministic.",
    sep = " "
  )
  reference <- paste(
    "Zhou YT, Chu JH, Zhao SH, Li GL, Fu ZY, Zhang SJ, Gao XH, Ma W,",
    "Shen K, Gao Y, Li W, Yin YM, Zhao C.",
    "Quantitative systems pharmacology modeling of HER2-positive",
    "metastatic breast cancer for translational efficacy evaluation",
    "and combination assessment across therapeutic modalities.",
    "Acta Pharmacol Sin. 2024;45(6):1287-1304.",
    "doi:10.1038/s41401-024-01232-9.",
    "Model compartments, species, initial conditions, parameters,",
    "reactions and the full ODE system are taken from Supplementary",
    "Table S1 (file 41401_2024_1232_MOESM2_ESM.xlsx, sheets",
    "Compartments / Species / Parameters / Reactions / Equations).",
    sep = " "
  )
  vignette <- "Zhou_2024_HER2breastCancer_qsp"

  # Every state below is a paper-mechanistic species from Supplementary
  # Table S1 sheet Species. rxode2 identifiers replace the supplement's
  # compartment-qualified names (e.g. Cell.[NRG1:E2:E3_p] -> NRG1_E2_E3p,
  # Vp_ADC.ADC -> adc_peripheral).
  paper_specific_compartments <- c(
    # --- pharmacokinetics -------------------------------------------
    "adc_central", "adc_peripheral", "mab_central", "mab_peripheral",
    "pl_central", "pl_peripheral", "dar",
    "cap_depot", "cap_central", "dfcr", "dfur", "fu",
    "lap_depot", "lap_central", "lap_peripheral",
    "pyr_depot", "pyr_central", "pyr_peripheral",
    # --- ErbB receptor network --------------------------------------
    "EGF", "NRG1", "E1", "E2", "E3", "E4",
    "EGF_E1", "EGF_E1_E2", "EGF_E1_E2p",
    "E2_E2", "E2_E2p", "E2_E3", "E2_E3p",
    "NRG1_E3", "NRG1_E2_E3", "NRG1_E2_E3p",
    "NRG1_E4", "NRG1_E2_E4", "NRG1_E2_E4p",
    # --- downstream signalling --------------------------------------
    "Raf", "Raf_p", "ERK", "ERK_p", "PI3K", "PI3K_p", "AKT", "AKT_p",
    "feedback", "BTK", "BTK_p",
    # --- tumour ADC handling and growth -----------------------------
    "ADC_ex", "ADC_E2", "ADC_E2_int", "PL_cell", "PL_ex", "Cells"
  )

  units <- list(
    time          = "day",
    dosing        = "nmol/L",
    concentration = "nM"
  )

  # Every state is a concentration (nmol/L) except dar (molecules of
  # payload per antibody) and Cells (tumour volume in mm3 in vivo, or a
  # relative cell number in vitro).  Source: Supplementary Table S1
  # sheets Species and Compartments.
  compartmentData <- list(
    adc_central    = list(analyte = "conjugated ADC (T-DM1 or T-DXd)", units = "nmol/L", specimen = "plasma", verified = TRUE),
    adc_peripheral = list(analyte = "conjugated ADC (T-DM1 or T-DXd)", units = "nmol/L", specimen = "tissue", verified = TRUE),
    mab_central    = list(analyte = "naked trastuzumab antibody", units = "nmol/L", specimen = "plasma", verified = TRUE),
    mab_peripheral = list(analyte = "naked trastuzumab antibody", units = "nmol/L", specimen = "tissue", verified = TRUE),
    pl_central     = list(analyte = "released payload (DM1 or DXd)", units = "nmol/L", specimen = "plasma", verified = TRUE),
    pl_peripheral  = list(analyte = "released payload (DM1 or DXd)", units = "nmol/L", specimen = "tissue", verified = TRUE),
    dar            = list(analyte = "drug-antibody ratio", units = "molecules per antibody", specimen = "not applicable", verified = TRUE),
    cap_depot      = list(analyte = "capecitabine", units = "nmol/L", specimen = "administration site", verified = TRUE),
    cap_central    = list(analyte = "capecitabine", units = "nmol/L", specimen = "plasma", verified = TRUE),
    dfcr           = list(analyte = "5'-deoxy-5-fluorocytidine (5'DFCR)", units = "nmol/L", specimen = "plasma", verified = TRUE),
    dfur           = list(analyte = "5'-deoxy-5-fluorouridine (5'DFUR)", units = "nmol/L", specimen = "plasma", verified = TRUE),
    fu             = list(analyte = "5-fluorouracil", units = "nmol/L", specimen = "plasma", verified = TRUE),
    lap_depot      = list(analyte = "lapatinib", units = "nmol/L", specimen = "administration site", verified = TRUE),
    lap_central    = list(analyte = "lapatinib", units = "nmol/L", specimen = "plasma", verified = TRUE),
    lap_peripheral = list(analyte = "lapatinib", units = "nmol/L", specimen = "tissue", verified = TRUE),
    pyr_depot      = list(analyte = "pyrotinib", units = "nmol/L", specimen = "administration site", verified = TRUE),
    pyr_central    = list(analyte = "pyrotinib", units = "nmol/L", specimen = "plasma", verified = TRUE),
    pyr_peripheral = list(analyte = "pyrotinib", units = "nmol/L", specimen = "tissue", verified = TRUE),
    EGF            = list(analyte = "epidermal growth factor", units = "nmol/L", specimen = "tumor", verified = TRUE),
    NRG1           = list(analyte = "neuregulin-1 (heregulin)", units = "nmol/L", specimen = "tumor", verified = TRUE),
    E1             = list(analyte = "free EGFR (ErbB1)", units = "nmol/L", specimen = "tumor", verified = TRUE),
    E2             = list(analyte = "free HER2 (ErbB2)", units = "nmol/L", specimen = "tumor", verified = TRUE),
    E3             = list(analyte = "free HER3 (ErbB3)", units = "nmol/L", specimen = "tumor", verified = TRUE),
    E4             = list(analyte = "free HER4 (ErbB4)", units = "nmol/L", specimen = "tumor", verified = TRUE),
    EGF_E1         = list(analyte = "EGF-bound EGFR", units = "nmol/L", specimen = "tumor", verified = TRUE),
    EGF_E1_E2      = list(analyte = "EGF:EGFR/HER2 heterodimer", units = "nmol/L", specimen = "tumor", verified = TRUE),
    EGF_E1_E2p     = list(analyte = "phospho-EGF:EGFR/HER2 heterodimer", units = "nmol/L", specimen = "tumor", verified = TRUE),
    E2_E2          = list(analyte = "HER2/HER2 homodimer", units = "nmol/L", specimen = "tumor", verified = TRUE),
    E2_E2p         = list(analyte = "phospho-HER2/HER2 homodimer", units = "nmol/L", specimen = "tumor", verified = TRUE),
    E2_E3          = list(analyte = "HER2/HER3 heterodimer", units = "nmol/L", specimen = "tumor", verified = TRUE),
    E2_E3p         = list(analyte = "phospho-HER2/HER3 heterodimer", units = "nmol/L", specimen = "tumor", verified = TRUE),
    NRG1_E3        = list(analyte = "NRG1-bound HER3", units = "nmol/L", specimen = "tumor", verified = TRUE),
    NRG1_E2_E3     = list(analyte = "NRG1:HER3/HER2 heterodimer", units = "nmol/L", specimen = "tumor", verified = TRUE),
    NRG1_E2_E3p    = list(analyte = "phospho-NRG1:HER3/HER2 heterodimer", units = "nmol/L", specimen = "tumor", verified = TRUE),
    NRG1_E4        = list(analyte = "NRG1-bound HER4", units = "nmol/L", specimen = "tumor", verified = TRUE),
    NRG1_E2_E4     = list(analyte = "NRG1:HER4/HER2 heterodimer", units = "nmol/L", specimen = "tumor", verified = TRUE),
    NRG1_E2_E4p    = list(analyte = "phospho-NRG1:HER4/HER2 heterodimer", units = "nmol/L", specimen = "tumor", verified = TRUE),
    Raf            = list(analyte = "inactive Raf", units = "nmol/L", specimen = "tumor", verified = TRUE),
    Raf_p          = list(analyte = "phospho-Raf", units = "nmol/L", specimen = "tumor", verified = TRUE),
    ERK            = list(analyte = "inactive ERK", units = "nmol/L", specimen = "tumor", verified = TRUE),
    ERK_p          = list(analyte = "phospho-ERK", units = "nmol/L", specimen = "tumor", verified = TRUE),
    PI3K           = list(analyte = "inactive PI3K", units = "nmol/L", specimen = "tumor", verified = TRUE),
    PI3K_p         = list(analyte = "phospho-PI3K", units = "nmol/L", specimen = "tumor", verified = TRUE),
    AKT            = list(analyte = "inactive AKT", units = "nmol/L", specimen = "tumor", verified = TRUE),
    AKT_p          = list(analyte = "phospho-AKT", units = "nmol/L", specimen = "tumor", verified = TRUE),
    feedback       = list(analyte = "AKT-driven negative-feedback protein acting on HER3 synthesis", units = "nmol/L", specimen = "tumor", verified = TRUE),
    BTK            = list(analyte = "inactive BTK", units = "nmol/L", specimen = "tumor", verified = TRUE),
    BTK_p          = list(analyte = "phospho-BTK", units = "nmol/L", specimen = "tumor", verified = TRUE),
    ADC_ex         = list(analyte = "conjugated ADC in tumour extracellular space", units = "nmol/L", specimen = "tumor", verified = TRUE),
    ADC_E2         = list(analyte = "HER2-bound ADC", units = "nmol/L", specimen = "tumor", verified = TRUE),
    ADC_E2_int     = list(analyte = "internalised ADC:HER2 complex", units = "nmol/L", specimen = "tumor", verified = TRUE),
    PL_cell        = list(analyte = "intracellular released payload", units = "nmol/L", specimen = "tumor", verified = TRUE),
    PL_ex          = list(analyte = "extracellular released payload in the tumour", units = "nmol/L", specimen = "tumor", verified = TRUE),
    Cells          = list(analyte = "tumour burden", units = "mm3 in vivo or relative cell number in vitro", specimen = "tumor", verified = TRUE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list()

  population <- list(
    species        = paste(
      "in vitro (SKBR3 primary; BT-474, NCI-N87 and ZR-75-1",
      "alternates) and mouse (BALB/c nude SKBR3 and BT-474 CDX,",
      "KPL4 CDX and a HER2-overexpressing PDX xenograft)"
    ),
    n_subjects     = NA_integer_,
    n_studies      = NA_integer_,
    disease_state  = paste(
      "HER2-positive breast cancer. Cell-level calibration used SKBR3",
      "(HER2 IHC 3+, 1.5e6 HER2 receptors per cell); generality was",
      "shown in BT-474 (IHC 3+), NCI-N87 (IHC 2+, 0.5e6 receptors) and",
      "ZR-75-1 (IHC 1+, 0.1e6 receptors). In vivo translation used",
      "SKBR3, KPL4 and BT-474 mouse xenografts plus a",
      "HER2-overexpressing PDX for T-DXd."
    ),
    dose_range     = paste(
      "In vivo mouse regimens reproduced by the paper: lapatinib",
      "100 mg/kg PO daily, pyrotinib 10-30 mg/kg PO daily,",
      "capecitabine 200-400 mg/kg PO daily, T-DM1 0.2-20 mg/kg IV",
      "single dose, T-DXd 10 mg/kg IV single dose. The authors'",
      "in-house SKBR3 experiment treated from a starting tumour",
      "volume of about 80 mm3 for 14 days."
    ),
    regions        = "Preclinical (Nanjing Medical University plus digitised literature datasets)",
    notes          = paste(
      "Scenario-specific parameter sets from Supplementary Table S1.",
      "T-DXd instead of T-DM1: Vc_PL = 0.0416 L, Vp_ADC = 0.000896 L,",
      "Vp_PL = 0.02 L, CL_ADC = 0.000138 L/day, CL_PL = 1.61088 L/day,",
      "dar(0) = 8, kdec_ADC_1 = kdec_ADC_2 = 0.023 /day,",
      "kcl_ADC_2 = 0.1605 /day/molecule, kout_PL = 0,",
      "kper_PL = 0.0001 /min, p3 = 1e-14 /day, and in vivo",
      "w14 = 35, n14 = 1, km14 = 0.066.",
      "In vitro (SKBR3): umax = 0.1 /h, dmax = 0.01 /h,",
      "Cellsmax = 12 relative cells, Cells(0) = 1, km13 = 1,",
      "n13 = 0.1, km10 = 40000, n10 = 1, w15 = 2.3, and",
      "w14 / n14 / km14 = 3.2 / 2 / 0.07 (T-DM1) or",
      "2.66 / 1 / 0.066 (T-DXd).",
      "BT-474 CDX in vivo: umax = 0.014 /h, dmax = 0.0001 /h,",
      "km13 = 1.3, n13 = 0.3, w14 = 500, n14 = 0.8, km14 = 0.35.",
      "Other cell lines in vitro (E2(0) / E3(0) / E2_E2(0) /",
      "E2_E2p(0) / E2_E3(0) / E2_E3p(0) / ks2 / ks3 / feedback(0)):",
      "BT-474 as SKBR3; NCI-N87 = 264.8 / 52.85 / 107.9 / 168.1 /",
      "10.77 / 2.503 / 3.057 / 20.9 / 2.739;",
      "ZR-75-1 = 92.8 / 60.72 / 13.25 / 20.64 / 4.334 / 1.007 /",
      "0.46 / 19 / 0.08967, with the Raf / ERK / PI3K / AKT",
      "phospho-fractions in Supplementary Table S1 sheet Species.",
      "Per-cell-line in vitro growth and drug-potency parameters",
      "(km13 / n13 / km10 / w15, then T-DM1 w14 / n14 / km14):",
      "SKBR3 = 1 / 0.1 / 40000 / 2.3, then 3.2 / 2 / 0.07",
      "(T-DXd in SKBR3 = 2.66 / 1 / 0.066);",
      "BT-474 = 1.3 / 0.3 / 40000 / 2.3, then 1.2 / 2.7 / 0.04;",
      "NCI-N87 = 0.0741 / 0.1658 / 1000 / 5.7615, then",
      "3.7855 / 1 / 0.6279;",
      "ZR-75-1 = 0.0675 / 0.1025 / 127090 / 9.8864, then",
      "7.5476 / 1 / 5.0835. n10 = 1 in every calibrated scenario."
    )
  )

  ini({
    # =================================================================
    # NOTE ON UNITS.  Supplementary Table S1 reports cell-level rate
    # constants per MINUTE, PK rate constants per DAY, and the tumour
    # proliferation / death rates per HOUR.  The model time unit is DAY,
    # so per-minute values are multiplied by 1440 and per-hour values by
    # 24.  The arithmetic is left visible so each supplement value can
    # be read directly off the expression.
    # =================================================================

    # -----------------------------------------------------------------
    # Compartment volumes (Supplementary Table S1, sheet Compartments)
    # -----------------------------------------------------------------
    Vc      <- fixed(0.00086)  ; label("Central volume for ADC, mAb, lapatinib and pyrotinib (L)")   # Compartments: Vc, PMID 27029797 / 23315145
    Vc_PL   <- fixed(0.066)    ; label("Central volume for released payload (L)")                    # Compartments: Vc_PL, T-DM1 (PMID 27029797); T-DXd 0.0416
    Vp_ADC  <- fixed(0.001896) ; label("Peripheral volume for ADC and mAb (L)")                      # Compartments: Vp_ADC, T-DM1 (PMID 27029797); T-DXd 0.000896
    Vp_PL   <- fixed(0.0402)   ; label("Peripheral volume for released payload (L)")                 # Compartments: Vp_PL, T-DM1 (PMID 27029797); T-DXd 0.02
    Vp_lap  <- fixed(0.0139)   ; label("Peripheral volume for lapatinib (L)")                        # Compartments: Vp_lap, fitted
    Vp_pyr  <- fixed(0.06)     ; label("Peripheral volume for pyrotinib (L)")                        # Compartments: Vp_pyr, fitted
    V_cap   <- fixed(0.04)     ; label("Apparent distribution volume for capecitabine (L)")          # Compartments: V_cap, fitted
    V_met   <- fixed(0.03)     ; label("Apparent distribution volume for capecitabine metabolites (L)")  # Compartments: V_met, fitted
    Cell    <- fixed(1e-12)    ; label("Volume of a single tumour cell (L)")                         # Compartments: Cell, PMID 28374319

    # -----------------------------------------------------------------
    # ADC and payload disposition (Table S1 sheet Parameters; 1/day)
    # -----------------------------------------------------------------
    CL_ADC     <- fixed(0.000186) ; label("Central clearance of ADC (L/day)")                        # Parameters: CL_ADC, T-DM1 (PMID 27029797); T-DXd 0.000138
    CLD_ADC    <- fixed(0.00236)  ; label("Distributional clearance of ADC (L/day)")                 # Parameters: CLD_ADC, PMID 27029797
    CL_PL      <- fixed(0.2258)   ; label("Central clearance of released payload (L/day)")           # Parameters: CL_PL, T-DM1 (PMID 27029797); T-DXd 1.61088
    CLD_PL     <- fixed(3.1)      ; label("Distributional clearance of released payload (L/day)")    # Parameters: CLD_PL, PMID 27029797
    kdec_ADC_1 <- fixed(0.241)    ; label("Systemic ADC deconjugation rate constant (1/day)")        # Parameters: kdec_ADC_1, T-DM1 (PMID 27029797); T-DXd 0.023
    kdec_ADC_2 <- fixed(0.241)    ; label("Systemic ADC deconjugation rate constant releasing payload (1/day/molecule)")  # Parameters: kdec_ADC_2, T-DM1; T-DXd 0.023
    kcl_ADC_2  <- fixed(0.2163)   ; label("ADC elimination rate constant releasing payload (1/day/molecule)")             # Parameters: kcl_ADC_2, T-DM1; T-DXd 0.1605

    # -----------------------------------------------------------------
    # Capecitabine and metabolites (Table S1 sheet Parameters; 1/day)
    # The supplement names these ka_cap, CL10, CL12, k10, k12, k23, k34
    # and k40; the k-constants carry a _cap suffix here so they cannot
    # be confused with the canonical two-compartment micro-constants.
    # -----------------------------------------------------------------
    ka_cap  <- fixed(50)   ; label("Absorption rate constant of capecitabine (1/day)")               # Parameters: ka_cap, PMID 16284918
    CL10    <- fixed(1)    ; label("Apparent elimination clearance of capecitabine (L/day)")         # Parameters: CL10, fitted
    CL12    <- fixed(15)   ; label("Apparent metabolic clearance of capecitabine to 5'DFCR (L/day)") # Parameters: CL12, fitted
    k23_cap <- fixed(70)   ; label("Conversion rate constant of 5'DFCR to 5'DFUR (1/day)")           # Parameters: k23, fitted
    k34_cap <- fixed(200)  ; label("Conversion rate constant of 5'DFUR to 5-FU (1/day)")             # Parameters: k34, fitted
    k40_cap <- fixed(3000) ; label("Elimination rate constant of 5-FU (1/day)")                      # Parameters: k40, fitted

    # -----------------------------------------------------------------
    # Lapatinib and pyrotinib disposition (Table S1 Parameters; 1/day)
    # -----------------------------------------------------------------
    ka_lap  <- fixed(5.688)  ; label("Absorption rate constant of lapatinib (1/day)")                # Parameters: ka_lap, PMID 23315145
    CLD_lap <- fixed(50)     ; label("Distributional clearance of lapatinib (L/day)")                # Parameters: CLD_lap, fitted
    CL_lap  <- fixed(1.6)    ; label("Central clearance of lapatinib (L/day)")                       # Parameters: CL_lap, fitted
    ka_pyr  <- fixed(7.9375) ; label("Absorption rate constant of pyrotinib (1/day)")                # Parameters: ka_pyr, fitted
    CLD_pyr <- fixed(50)     ; label("Distributional clearance of pyrotinib (L/day)")                # Parameters: CLD_pyr, fitted
    CL_pyr  <- fixed(2.25)   ; label("Central clearance of pyrotinib (L/day)")                       # Parameters: CL_pyr, fitted

    # -----------------------------------------------------------------
    # ErbB receptor synthesis and degradation (per minute -> per day)
    # ks1 and ks4 are initial assignments ks1 = kd1*E1(0), ks4 = kd4*E4(0)
    # -----------------------------------------------------------------
    ks2 <- fixed(10.29 * 1440) ; label("Synthesis rate of HER2 (nmol/L/day)")                        # Parameters: ks2 = 10.29 nM/min (SKBR3); NCI-N87 3.057, ZR-75-1 0.46
    ks3 <- fixed(45 * 1440)    ; label("Maximum synthesis rate of HER3 (nmol/L/day)")                # Parameters: ks3 = 45 nM/min (SKBR3); NCI-N87 20.9, ZR-75-1 19
    kd1 <- fixed(0.0018 * 1440); label("Degradation rate constant of EGFR (1/day)")                  # Parameters: kd1 = 0.0018 /min, PMID 8617810
    kd2 <- fixed(0.0013 * 1440); label("Degradation rate constant of HER2 (1/day)")                  # Parameters: kd2 = 0.0013 /min, PMID 34725188
    kd3 <- fixed(0.31 * 1440)  ; label("Degradation rate constant of HER3 (1/day)")                  # Parameters: kd3 = 0.31 /min, fitted
    kd4 <- fixed(0.0021 * 1440); label("Degradation rate constant of HER4 (1/day)")                  # Parameters: kd4 = 0.0021 /min, PMID 8617810

    # -----------------------------------------------------------------
    # Ligand binding, dimerisation and phosphorylation (per min -> day)
    # -----------------------------------------------------------------
    k1on  <- fixed(0.0618 * 1440); label("EGF binding to EGFR association rate (L/nmol/day)")        # Parameters: k1on = 0.0618 /nM/min, PMID 30796217
    k1off <- fixed(0.34 * 1440)  ; label("EGF:EGFR dissociation rate constant (1/day)")              # Parameters: k1off = 0.34 /min, PMID 30796217
    k2on  <- fixed(0.6 * 1440)   ; label("EGF:EGFR to HER2 dimerisation rate (L/nmol/day)")          # Parameters: k2on = 0.6 /nM/min, PMID 16533841 and fitting
    k2off <- fixed(1.2 * 1440)   ; label("EGF:EGFR/HER2 dissociation rate constant (1/day)")         # Parameters: k2off = 1.2 /min
    k3on  <- fixed(0.07 * 1440)  ; label("EGF:EGFR/HER2 phosphorylation rate constant (1/day)")      # Parameters: k3on = 0.07 /min
    k3off <- fixed(0.02 * 1440)  ; label("EGF:EGFR/HER2 dephosphorylation rate constant (1/day)")    # Parameters: k3off = 0.02 /min
    k4on  <- fixed(0.2 * 1440)   ; label("HER2 homodimerisation rate constant (L/nmol/day)")         # Parameters: k4on = 0.2 /nM/min
    k4off <- fixed(130 * 1440)   ; label("HER2/HER2 dissociation rate constant (1/day)")             # Parameters: k4off = 130 /min
    k5on  <- fixed(0.48 * 1440)  ; label("HER2/HER2 phosphorylation rate constant (1/day)")          # Parameters: k5on = 0.48 /min
    k5off <- fixed(0.3 * 1440)   ; label("HER2/HER2 dephosphorylation rate constant (1/day)")        # Parameters: k5off = 0.3 /min
    k6on  <- fixed(0.1 * 1440)   ; label("HER2 to HER3 dimerisation rate constant (L/nmol/day)")     # Parameters: k6on = 0.1 /nM/min
    k6off <- fixed(130 * 1440)   ; label("HER2/HER3 dissociation rate constant (1/day)")             # Parameters: k6off = 130 /min
    k7on  <- fixed(0.56 * 1440)  ; label("HER2/HER3 phosphorylation rate constant (1/day)")          # Parameters: k7on = 0.56 /min
    k7off <- fixed(2.4 * 1440)   ; label("HER2/HER3 dephosphorylation rate constant (1/day)")        # Parameters: k7off = 2.4 /min
    k8on  <- fixed(0.354 * 1440) ; label("NRG1 binding to HER3 association rate (L/nmol/day)")       # Parameters: k8on = 0.354 /nM/min, PMID 10745024 and fitting
    k8off <- fixed(0.54 * 1440)  ; label("NRG1:HER3 dissociation rate constant (1/day)")             # Parameters: k8off = 0.54 /min
    k9on  <- fixed(0.6 * 1440)   ; label("NRG1:HER3 to HER2 dimerisation rate (L/nmol/day)")         # Parameters: k9on = 0.6 /nM/min
    k9off <- fixed(0.2 * 1440)   ; label("NRG1:HER3/HER2 dissociation rate constant (1/day)")        # Parameters: k9off = 0.2 /min
    k10on <- fixed(2.3 * 1440)   ; label("NRG1:HER3/HER2 phosphorylation rate constant (1/day)")     # Parameters: k10on = 2.3 /min
    k10off<- fixed(0.03 * 1440)  ; label("NRG1:HER3/HER2 dephosphorylation rate constant (1/day)")   # Parameters: k10off = 0.03 /min
    k11on <- fixed(0.15 * 1440)  ; label("NRG1 binding to HER4 association rate (L/nmol/day)")       # Parameters: k11on = 0.15 /nM/min, PMID 10745024 and fitting
    k11off<- fixed(0.024 * 1440) ; label("NRG1:HER4 dissociation rate constant (1/day)")             # Parameters: k11off = 0.024 /min, fitted
    k12on <- fixed(0.6 * 1440)   ; label("NRG1:HER4 to HER2 dimerisation rate (L/nmol/day)")         # Parameters: k12on = 0.6 /nM/min, fitted
    k12off<- fixed(1.2 * 1440)   ; label("NRG1:HER4/HER2 dissociation rate constant (1/day)")        # Parameters: k12off = 1.2 /min, fitted
    k13on <- fixed(0.96 * 1440)  ; label("NRG1:HER4/HER2 phosphorylation rate constant (1/day)")     # Parameters: k13on = 0.96 /min, fitted
    k13off<- fixed(0.3 * 1440)   ; label("NRG1:HER4/HER2 dephosphorylation rate constant (1/day)")   # Parameters: k13off = 0.3 /min, fitted

    # -----------------------------------------------------------------
    # Dimer degradation (per minute -> per day)
    # -----------------------------------------------------------------
    kd12    <- fixed(1 * 1440)      ; label("TKI-induced degradation rate of EGF:EGFR/HER2 (1/day)")      # Parameters: kd12 = 1 /min, fitted
    kd22    <- fixed(1 * 1440)      ; label("TKI-induced degradation rate of HER2/HER2 (1/day)")          # Parameters: kd22 = 1 /min, fitted
    kd23_1  <- fixed(1 * 1440)      ; label("TKI-induced degradation rate of NRG1:HER3/HER2 (1/day)")     # Parameters: kd23_1 = 1 /min, fitted
    kd23_2  <- fixed(0.125 * 1440)  ; label("TKI-induced degradation rate of HER2/HER3 (1/day)")          # Parameters: kd23_2 = 0.125 /min, fitted
    kd24    <- fixed(1 * 1440)      ; label("TKI-induced degradation rate of NRG1:HER4/HER2 (1/day)")     # Parameters: kd24 = 1 /min, fitted
    kd12p   <- fixed(0.01 * 1440)   ; label("Degradation rate of phospho-EGF:EGFR/HER2 (1/day)")          # Parameters: kd12p = 0.01 /min, fitted
    kd22p   <- fixed(0.008 * 1440)  ; label("Degradation rate of phospho-HER2/HER2 (1/day)")              # Parameters: kd22p = 0.008 /min, fitted
    kd23p_1 <- fixed(0.0003 * 1440) ; label("Degradation rate of phospho-NRG1:HER3/HER2 (1/day)")         # Parameters: kd23p_1 = 0.0003 /min, fitted
    kd23p_2 <- fixed(0.009 * 1440)  ; label("Degradation rate of phospho-HER2/HER3 (1/day)")              # Parameters: kd23p_2 = 0.009 /min, fitted
    kd24p   <- fixed(0.008 * 1440)  ; label("Degradation rate of phospho-NRG1:HER4/HER2 (1/day)")         # Parameters: kd24p = 0.008 /min, fitted

    # -----------------------------------------------------------------
    # Downstream signalling phospho-cycles (per minute -> per day)
    # -----------------------------------------------------------------
    kon_PI3K  <- fixed(1 * 1440)      ; label("PI3K phosphorylation rate constant (L/nmol/day)")      # Parameters: kon_PI3K = 1 /nM/min, fitted
    koff_PI3K <- fixed(28 * 1440)     ; label("PI3K dephosphorylation rate constant (1/day)")         # Parameters: koff_PI3K = 28 /min, fitted
    kon_AKT   <- fixed(10 * 1440)     ; label("AKT phosphorylation rate constant (1/day)")            # Parameters: kon_AKT = 10 /min, fitted
    koff_AKT  <- fixed(0.4 * 1440)    ; label("AKT dephosphorylation rate constant (1/day)")          # Parameters: koff_AKT = 0.4 /min, fitted
    kon_Raf   <- fixed(1 * 1440)      ; label("Raf phosphorylation rate constant (L/nmol/day)")       # Parameters: kon_Raf = 1 /nM/min, fitted
    koff_Raf  <- fixed(100 * 1440)    ; label("Raf dephosphorylation rate constant (1/day)")          # Parameters: koff_Raf = 100 /min, fitted
    kon_ERK   <- fixed(50 * 1440)     ; label("ERK phosphorylation rate constant (1/day)")            # Parameters: kon_ERK = 50 /min, fitted
    koff_ERK  <- fixed(8 * 1440)      ; label("ERK dephosphorylation rate constant (1/day)")          # Parameters: koff_ERK = 8 /min, fitted
    kon_BTK   <- fixed(1 * 1440)      ; label("BTK phosphorylation rate constant (1/day)")            # Parameters: kon_BTK = 1 /min, fitted
    koff_BTK  <- fixed(0.0005 * 1440) ; label("BTK dephosphorylation rate constant (1/day)")          # Parameters: koff_BTK = 0.0005 /min, fitted
    kfb1      <- fixed(1 * 1440)      ; label("Activation rate of the HER3 negative-feedback protein (1/day)")    # Parameters: kfb1 = 1 /min, fitted
    kfb2      <- fixed(0.001 * 1440)  ; label("Inactivation rate of the HER3 negative-feedback protein (1/day)")  # Parameters: kfb2 = 0.001 /min, fitted

    # -----------------------------------------------------------------
    # Tumour ADC binding, internalisation and payload handling
    # -----------------------------------------------------------------
    kon_ADC     <- fixed(0.00614 * 1440) ; label("ADC association rate to HER2 (L/nmol/day)")            # Parameters: kon_ADC = 0.00614 /nM/min, PMID 26912181
    koff_ADC    <- fixed(0.00023 * 1440) ; label("ADC:HER2 dissociation rate constant (1/day)")          # Parameters: koff_ADC = 0.00023 /min, PMID 26912181
    kint_ADC    <- fixed(0.0015 * 1440)  ; label("Internalisation rate constant of the ADC:HER2 complex (1/day)")  # Parameters: kint_ADC = 0.0015 /min, PMID 26912181
    kdeg_ADC_1  <- fixed(0.00063 * 1440) ; label("Proteasomal degradation rate constant of internalised ADC (1/day)")  # Parameters: kdeg_ADC_1 = 0.00063 /min, PMID 26912181
    kdeg_ADC_2  <- fixed(0.00063 * 1440) ; label("Proteasomal degradation rate constant releasing payload (1/day/molecule)")  # Parameters: kdeg_ADC_2 = 0.00063 /min, PMID 26912181
    kout_PL     <- fixed(0.00025 * 1440) ; label("Efflux rate constant of intracellular payload (1/day)")     # Parameters: kout_PL = 0.00025 /min for T-DM1 (PMID 26912181); 0 assumed for T-DXd
    kper_PL     <- fixed(0 * 1440)       ; label("Bidirectional payload permeation rate constant (1/day)")    # Parameters: kper_PL = 0 assumed for T-DM1; 0.0001 /min fitted for T-DXd
    kdec_ADC_3  <- fixed(0.000372 * 1440); label("Extracellular ADC deconjugation rate constant (1/day)")     # Parameters: kdec_ADC_3 = 0.000372 /min, PMID 27029797
    kdec_ADC_4  <- fixed(0.000372 * 1440); label("Extracellular ADC deconjugation rate releasing payload (1/day/molecule)")  # Parameters: kdec_ADC_4 = 0.000372 /min, PMID 27029797

    # -----------------------------------------------------------------
    # Hill functions: HER3 negative feedback and TKI actions
    # -----------------------------------------------------------------
    kE3  <- fixed(1)      ; label("Hill coefficient for the HER3 negative feedback")                 # Parameters: kE3, fitted
    kmE3 <- fixed(10)     ; label("Hill constant for the HER3 negative feedback")                    # Parameters: kmE3, fitted
    km1  <- fixed(1000)   ; label("Hill constant, lapatinib inhibition of EGF:EGFR/HER2 phosphorylation")   # Parameters: km1, fitted
    km2  <- fixed(1000)   ; label("Hill constant, lapatinib-induced degradation of EGF:EGFR/HER2")   # Parameters: km2, fitted
    km3  <- fixed(5000)   ; label("Hill constant, lapatinib inhibition of HER2/HER2 phosphorylation")# Parameters: km3, fitted
    km4  <- fixed(500)    ; label("Hill constant, lapatinib-induced degradation of HER2/HER2")       # Parameters: km4, fitted
    km5  <- fixed(2500)   ; label("Hill constant, lapatinib inhibition of HER2/HER3 phosphorylation")# Parameters: km5, fitted
    km6  <- fixed(500)    ; label("Hill constant, lapatinib-induced degradation of HER2/HER3")       # Parameters: km6, fitted
    km7  <- fixed(150000) ; label("Hill constant, lapatinib inhibition of NRG1:HER3/HER2 phosphorylation")  # Parameters: km7, fitted
    km8  <- fixed(100000) ; label("Hill constant, lapatinib-induced degradation of NRG1:HER3/HER2")  # Parameters: km8, fitted
    km9  <- fixed(0.5)    ; label("Hill constant, phospho-Raf activation of ERK")                    # Parameters: km9, fitted
    km10 <- fixed(40000)  ; label("Hill constant, 5-FU stimulation of tumour cell death")            # Parameters: km10 = 40000 (SKBR3 and in vivo); NCI-N87 1000, ZR-75-1 127090
    km11 <- fixed(10)     ; label("Hill constant, phospho-PI3K activation of AKT")                   # Parameters: km11, fitted
    km12 <- fixed(1)      ; label("Hill constant, PTEN inhibition of AKT activation")                # Parameters: km12, fitted
    km13 <- fixed(1)      ; label("Hill constant, phospho-ERK / AKT / BTK stimulation of tumour growth")  # Parameters: km13 = 1 (SKBR3 in vivo); BT-474 1.3, NCI-N87 0.0741, ZR-75-1 0.0675
    km14 <- fixed(0.07)   ; label("Hill constant, payload stimulation of tumour cell death")         # Parameters: km14 = 0.07 in vivo/in vitro (T-DM1, SKBR3); T-DXd 0.066; in vivo BT-474 CDX 0.35; in vitro T-DM1 BT-474 0.04, NCI-N87 0.6279, ZR-75-1 5.0835
    km15 <- fixed(11)     ; label("Hill constant, pyrotinib inhibition of EGF:EGFR/HER2 phosphorylation")   # Parameters: km15, fitted
    km16 <- fixed(217)    ; label("Hill constant, pyrotinib inhibition of HER2/HER2 phosphorylation")# Parameters: km16, fitted
    km17 <- fixed(1)      ; label("Hill constant, pyrotinib inhibition of HER2/HER3 phosphorylation")# Parameters: km17, fitted
    km18 <- fixed(2)      ; label("Hill constant, pyrotinib inhibition of NRG1:HER3/HER2 phosphorylation")  # Parameters: km18, fitted
    km19 <- fixed(19)     ; label("Hill constant, pyrotinib inhibition of NRG1:HER4/HER2 phosphorylation")  # Parameters: km19, fitted
    km20 <- fixed(2272)   ; label("Hill constant, pyrotinib-induced degradation of NRG1:HER4/HER2")  # Parameters: km20, fitted
    km21 <- fixed(80)     ; label("Hill constant, phospho-NRG1 dimer activation of BTK")             # Parameters: km21, fitted
    km22 <- fixed(20000)  ; label("Hill constant, pyrotinib-induced degradation of EGF:EGFR/HER2")   # Parameters: km22, fitted
    km23 <- fixed(430)    ; label("Hill constant, pyrotinib-induced degradation of HER2/HER2")       # Parameters: km23, fitted
    km24 <- fixed(24)     ; label("Hill constant, pyrotinib-induced degradation of HER2/HER3")       # Parameters: km24, fitted
    km25 <- fixed(61986)  ; label("Hill constant, pyrotinib-induced degradation of NRG1:HER3/HER2")  # Parameters: km25, fitted

    n1  <- fixed(1)   ; label("Hill coefficient, lapatinib inhibition of EGF:EGFR/HER2 phosphorylation")    # Parameters: n1, fitted
    n2  <- fixed(1)   ; label("Hill coefficient, lapatinib-induced degradation of EGF:EGFR/HER2")     # Parameters: n2, fitted
    n3  <- fixed(2)   ; label("Hill coefficient, lapatinib inhibition of HER2/HER2 phosphorylation")  # Parameters: n3, fitted
    n4  <- fixed(1)   ; label("Hill coefficient, lapatinib-induced degradation of HER2/HER2")         # Parameters: n4, fitted
    n5  <- fixed(2)   ; label("Hill coefficient, lapatinib inhibition of HER2/HER3 phosphorylation")  # Parameters: n5, fitted
    n6  <- fixed(1)   ; label("Hill coefficient, lapatinib-induced degradation of HER2/HER3")         # Parameters: n6, fitted
    n7  <- fixed(2)   ; label("Hill coefficient, lapatinib inhibition of NRG1:HER3/HER2 phosphorylation")   # Parameters: n7, fitted
    n8  <- fixed(2)   ; label("Hill coefficient, lapatinib-induced degradation of NRG1:HER3/HER2")    # Parameters: n8, fitted
    n9  <- fixed(2)   ; label("Hill coefficient, phospho-Raf activation of ERK")                      # Parameters: n9, fitted
    n10 <- fixed(1)   ; label("Hill coefficient, 5-FU stimulation of tumour cell death")              # Parameters: n10 = 1 in every calibrated scenario
    n11 <- fixed(3)   ; label("Hill coefficient, phospho-PI3K activation of AKT")                     # Parameters: n11, fitted
    n12 <- fixed(3)   ; label("Hill coefficient, PTEN inhibition of AKT activation")                  # Parameters: n12, fitted
    n13 <- fixed(0.1) ; label("Hill coefficient, phospho-ERK / AKT / BTK stimulation of tumour growth")     # Parameters: n13 = 0.1 (SKBR3 in vivo); BT-474 0.3, NCI-N87 0.1658, ZR-75-1 0.1025
    n14 <- fixed(2)   ; label("Hill coefficient, payload stimulation of tumour cell death")           # Parameters: n14 = 2 (T-DM1, in vivo and in vitro SKBR3); T-DXd 1; in vivo BT-474 CDX 0.8; in vitro T-DM1 BT-474 2.7, NCI-N87 1, ZR-75-1 1
    n15 <- fixed(1)   ; label("Hill coefficient, pyrotinib inhibition of EGF:EGFR/HER2 phosphorylation")    # Parameters: n15, fitted
    n16 <- fixed(3)   ; label("Hill coefficient, pyrotinib inhibition of HER2/HER2 phosphorylation")  # Parameters: n16, fitted
    n17 <- fixed(3)   ; label("Hill coefficient, pyrotinib inhibition of HER2/HER3 phosphorylation")  # Parameters: n17, fitted
    n18 <- fixed(1)   ; label("Hill coefficient, pyrotinib inhibition of NRG1:HER3/HER2 phosphorylation")   # Parameters: n18, fitted
    n19 <- fixed(1)   ; label("Hill coefficient, pyrotinib inhibition of NRG1:HER4/HER2 phosphorylation")   # Parameters: n19, fitted
    n20 <- fixed(1)   ; label("Hill coefficient, pyrotinib-induced degradation of NRG1:HER4/HER2")    # Parameters: n20, fitted
    n21 <- fixed(3)   ; label("Hill coefficient, phospho-NRG1 dimer activation of BTK")               # Parameters: n21, fitted
    n22 <- fixed(1)   ; label("Hill coefficient, pyrotinib-induced degradation of EGF:EGFR/HER2")     # Parameters: n22, fitted
    n23 <- fixed(1)   ; label("Hill coefficient, pyrotinib-induced degradation of HER2/HER2")         # Parameters: n23, fitted
    n24 <- fixed(1)   ; label("Hill coefficient, pyrotinib-induced degradation of HER2/HER3")         # Parameters: n24, fitted
    n25 <- fixed(1)   ; label("Hill coefficient, pyrotinib-induced degradation of NRG1:HER3/HER2")    # Parameters: n25, fitted

    # -----------------------------------------------------------------
    # Pathway weights
    # -----------------------------------------------------------------
    w1  <- fixed(30)    ; label("Weight, phospho-EGF:EGFR/HER2 activation of PI3K")                   # Parameters: w1, fitted
    w2  <- fixed(0.01)  ; label("Weight, phospho-HER2/HER2 activation of PI3K")                       # Parameters: w2, fitted
    w3  <- fixed(1)     ; label("Weight, phospho-HER2/HER3 activation of PI3K")                       # Parameters: w3, fitted
    w4  <- fixed(8)     ; label("Weight, phospho-NRG1:HER3/HER2 activation of PI3K")                  # Parameters: w4, fitted
    w5  <- fixed(100)   ; label("Weight, phospho-NRG1:HER4/HER2 activation of PI3K")                  # Parameters: w5, fitted
    w6  <- fixed(30)    ; label("Weight, phospho-EGF:EGFR/HER2 activation of Raf")                    # Parameters: w6, fitted
    w7  <- fixed(0.01)  ; label("Weight, phospho-HER2/HER2 activation of Raf")                        # Parameters: w7, fitted
    w8  <- fixed(0.3)   ; label("Weight, phospho-HER2/HER3 activation of Raf")                        # Parameters: w8, fitted
    w9  <- fixed(2)     ; label("Weight, phospho-NRG1:HER3/HER2 activation of Raf")                   # Parameters: w9, fitted
    w10 <- fixed(40)    ; label("Weight, phospho-NRG1:HER4/HER2 activation of Raf")                   # Parameters: w10, fitted
    w11 <- fixed(0.5)   ; label("Weight, phospho-ERK stimulation of tumour growth")                   # Parameters: w11, fitted
    w12 <- fixed(8)     ; label("Weight, phospho-AKT stimulation of tumour growth")                   # Parameters: w12, fitted
    w13 <- fixed(1)     ; label("Weight, phospho-BTK stimulation of tumour growth")                   # Parameters: w13, fitted
    w14 <- fixed(50)    ; label("Maximum payload stimulation of tumour cell death")                   # Parameters: w14 = 50 (T-DM1 in vivo); T-DXd 35; in vivo BT-474 CDX 500; in vitro T-DM1 SKBR3 3.2, T-DXd SKBR3 2.66, BT-474 1.2, NCI-N87 3.7855, ZR-75-1 7.5476
    w15 <- fixed(120)   ; label("Maximum 5-FU stimulation of tumour cell death")                      # Parameters: w15 = 120 in vivo; in vitro 2.3 (SKBR3 and BT-474), NCI-N87 5.7615, ZR-75-1 9.8864
    w16 <- fixed(9)     ; label("Weight, phospho-NRG1:HER3/HER2 activation of BTK")                   # Parameters: w16, fitted
    w17 <- fixed(0.001) ; label("Weight, phospho-NRG1:HER4/HER2 activation of BTK")                   # Parameters: w17, fitted

    # -----------------------------------------------------------------
    # Tumour partition coefficients, switches and growth parameters
    # -----------------------------------------------------------------
    p1 <- fixed(0.04)    ; label("Lapatinib tumour partition coefficient")                            # Parameters: p1, fitted
    p2 <- fixed(0.5)     ; label("Pyrotinib tumour partition coefficient")                            # Parameters: p2, fitted
    p3 <- fixed(3.5e-13) ; label("Rate of ADC entering the tumour microenvironment (1/day)")          # Parameters: p3 = 3.5e-13 (T-DM1); T-DXd 1e-14
    p4 <- fixed(35)      ; label("5-FU tumour partition coefficient")                                 # Parameters: p4, fitted
    s0 <- fixed(1)       ; label("Dimensionalising divisor for Hill-function arguments (nmol/L)")     # Parameters: s0
    Rx <- fixed(1)       ; label("Tumour-present switch, 1 when a tumour exists and 0 otherwise")     # Parameters: Rx
    PTEN <- fixed(1)     ; label("PTEN level (nmol/L), constant")                                     # Species: Cell.PTEN initial amount 1 nM, no ODE

    umax     <- fixed(0.015 * 24)  ; label("Maximum tumour proliferation rate (1/day)")               # Parameters: umax = 0.015 /h in vivo; BT-474 CDX 0.014 /h; in vitro 0.1 /h
    dmax     <- fixed(0.0005 * 24) ; label("Basal tumour death rate (1/day)")                         # Parameters: dmax = 0.0005 /h in vivo; BT-474 CDX 0.0001 /h; in vitro 0.01 /h
    Cellsmax <- fixed(2000)        ; label("Carrying capacity of the tumour (mm3)")                   # Parameters: Vmax = 2000 mm3 in vivo; 12 relative cells in vitro

    # -----------------------------------------------------------------
    # Initial conditions that a user is expected to vary
    # -----------------------------------------------------------------
    dar0    <- fixed(3.5) ; label("Drug-antibody ratio at dosing (molecules per antibody)")           # Species: DAR.DAR = 3.5 (T-DM1); 8 (T-DXd)
    cells0  <- fixed(80)  ; label("Tumour volume at the start of treatment (mm3)")                    # Species: Cells is unfixed in vivo; 80 mm3 in the authors in-house SKBR3 experiment
    E1_0    <- fixed(249)   ; label("Steady-state free EGFR (nmol/L)")                                # Species: Cell.E1 = 249 nM (150000 receptors per SKBR3 cell)
    E2_0    <- fixed(500.4) ; label("Steady-state free HER2 (nmol/L)")                                # Species: Cell.E2 = 500.4 nM (SKBR3 and BT-474); NCI-N87 264.8; ZR-75-1 92.8
    E3_0    <- fixed(44.83) ; label("Steady-state free HER3 (nmol/L)")                                # Species: Cell.E3 = 44.83 nM (SKBR3 and BT-474); NCI-N87 52.85; ZR-75-1 60.72
    E4_0    <- fixed(3)     ; label("Steady-state free HER4 (nmol/L)")                                # Species: Cell.E4 = 3 nM (2000 receptors per SKBR3 cell)
    E2_E2_0 <- fixed(385.1) ; label("Steady-state HER2/HER2 homodimer (nmol/L)")                      # Species: Cell.[E2:E2] = 385.1 nM (SKBR3); NCI-N87 107.9; ZR-75-1 13.25
    E2_E2p_0<- fixed(600.2) ; label("Steady-state phospho-HER2/HER2 homodimer (nmol/L)")              # Species: Cell.[E2:E2_p] = 600.2 nM (SKBR3); NCI-N87 168.1; ZR-75-1 20.64
    E2_E3_0 <- fixed(17.25) ; label("Steady-state HER2/HER3 heterodimer (nmol/L)")                    # Species: Cell.[E2:E3] = 17.25 nM (SKBR3); NCI-N87 10.77; ZR-75-1 4.334
    E2_E3p_0<- fixed(4.011) ; label("Steady-state phospho-HER2/HER3 heterodimer (nmol/L)")            # Species: Cell.[E2:E3_p] = 4.011 nM (SKBR3); NCI-N87 2.503; ZR-75-1 1.007
    Raf_0   <- fixed(0.9328); label("Steady-state inactive Raf (nmol/L)")                             # Species: Cell.Raf = 0.9328 (SKBR3); NCI-N87 0.9763; ZR-75-1 0.9949
    Rafp_0  <- fixed(0.0672); label("Steady-state phospho-Raf (nmol/L)")                              # Species: Cell.Raf_p = 0.0672 (SKBR3); NCI-N87 0.0237; ZR-75-1 0.0051
    ERK_0   <- fixed(0.947) ; label("Steady-state inactive ERK (nmol/L)")                             # Species: Cell.ERK = 0.947 (SKBR3); NCI-N87 0.993; ZR-75-1 0.9997
    ERKp_0  <- fixed(0.053) ; label("Steady-state phospho-ERK (nmol/L)")                              # Species: Cell.ERK_p = 0.053 (SKBR3); NCI-N87 0.007; ZR-75-1 0.0003
    PI3K_0  <- fixed(0.7366); label("Steady-state inactive PI3K (nmol/L)")                            # Species: Cell.PI3K = 0.7366 (SKBR3); NCI-N87 0.87; ZR-75-1 0.9584
    PI3Kp_0 <- fixed(0.2634); label("Steady-state phospho-PI3K (nmol/L)")                             # Species: Cell.PI3K_p = 0.2634 (SKBR3); NCI-N87 0.13; ZR-75-1 0.0416
    AKT_0   <- fixed(0.9777); label("Steady-state inactive AKT (nmol/L)")                             # Species: Cell.AKT = 0.9777 (SKBR3); NCI-N87 0.9973; ZR-75-1 0.9999
    AKTp_0  <- fixed(0.0223); label("Steady-state phospho-AKT (nmol/L)")                              # Species: Cell.AKT_p = 0.0223 (SKBR3); NCI-N87 0.0027; ZR-75-1 0.0001
    fb_0    <- fixed(22.3)  ; label("Steady-state HER3 negative-feedback protein (nmol/L)")           # Species: Cell.feedback = 22.3 (SKBR3); NCI-N87 2.739; ZR-75-1 0.08967
    BTK_0   <- fixed(1)     ; label("Steady-state inactive BTK (nmol/L)")                             # Species: Cell.BTK = 1 nM
  })

  model({
    # =================================================================
    # Derived rate constants (Table S1 sheet Parameters, the entries
    # flagged as -defined by an initial assignment-)
    # =================================================================
    k12_ADC   <- CLD_ADC / Vc
    k21_ADC   <- CLD_ADC / Vp_ADC
    k12_PL    <- CLD_PL / Vc_PL
    k21_PL    <- CLD_PL / Vp_PL
    kcl_PL    <- CL_PL / Vc_PL
    kcl_ADC_1 <- CL_ADC / Vc
    k10_cap   <- CL10 / V_cap
    k12_cap   <- CL12 / V_cap
    k12_lap   <- CLD_lap / Vc
    k21_lap   <- CLD_lap / Vp_lap
    kcl_lap   <- CL_lap / Vc
    k12_pyr   <- CLD_pyr / Vc
    k21_pyr   <- CLD_pyr / Vp_pyr
    kcl_pyr   <- CL_pyr / Vc
    ks1       <- kd1 * E1_0
    ks4       <- kd4 * E4_0

    # =================================================================
    # Repeated assignments.  Supplementary Table S1 sheet Species notes
    # give TTmAb explicitly.  The three tumour drug concentrations are
    # the partition-coefficient assignments described in Methods
    # (-drug exposures in the tumor are proportional to that in the
    # blood, as described by tumor partition coefficients for TKIs and
    # capecitabine-) using p1, p2 and p4; the assignment itself is not
    # printed as an equation, see the vignette Errata.
    # =================================================================
    TTmAb     <- adc_central + mab_central
    lap_tumor <- p1 * lap_central
    pyr_tumor <- p2 * pyr_central
    fu_tumor  <- p4 * fu

    # =================================================================
    # Reaction fluxes v1-v79 (Table S1 sheet Reactions, column Fluxes)
    # Every flux carries its compartment-volume factor exactly as
    # printed, and each ODE divides by its own compartment volume.
    # =================================================================

    # ---- ADC, naked antibody and payload disposition ----------------
    v1  <- (k12_ADC * adc_central) * Vc - (k21_ADC * adc_peripheral) * Vp_ADC
    v2  <- (k12_ADC * mab_central) * Vc - (k21_ADC * mab_peripheral) * Vp_ADC
    v3  <- (k12_PL * pl_central) * Vc_PL - (k21_PL * pl_peripheral) * Vp_PL
    v4  <- (kcl_PL * pl_central) * Vc_PL
    v5  <- (kcl_ADC_1 * adc_central) * Vc
    v6  <- (kcl_ADC_2 * dar * adc_central) * Vc
    v7  <- (kdec_ADC_1 * adc_central) * Vc
    v8  <- (kdec_ADC_2 * dar * adc_central) * Vc
    v9  <- (kcl_ADC_1 * mab_central) * Vc
    v10 <- kdec_ADC_1 * dar

    # ---- capecitabine and its metabolite cascade --------------------
    v11 <- (ka_cap * cap_depot) * V_cap
    v12 <- (k12_cap * cap_central) * V_cap
    v13 <- (k23_cap * dfcr) * V_met
    v14 <- (k34_cap * dfur) * V_met
    v15 <- (k40_cap * fu) * V_met
    v16 <- (k10_cap * cap_central) * V_cap

    # ---- lapatinib and pyrotinib ------------------------------------
    v17 <- (ka_lap * lap_depot) * Vc
    v18 <- (k12_lap * lap_central) * Vc - (k21_lap * lap_peripheral) * Vp_lap
    v19 <- (kcl_lap * lap_central) * Vc
    v20 <- (ka_pyr * pyr_depot) * Vc
    v21 <- (k12_pyr * pyr_central) * Vc - (k21_pyr * pyr_peripheral) * Vp_pyr
    v22 <- (kcl_pyr * pyr_central) * Vc

    # ---- ErbB receptor turnover -------------------------------------
    v23 <- ks1 * Cell
    v24 <- (kd1 * E1) * Cell
    v25 <- ks2 * Cell
    v26 <- (kd2 * E2) * Cell
    v27 <- (ks3 * (1 - (feedback / s0)^kE3 / (kmE3 + (feedback / s0)^kE3))) * Cell
    v28 <- (kd3 * E3) * Cell
    v29 <- ks4 * Cell
    v30 <- (kd4 * E4) * Cell

    # ---- EGF / EGFR / HER2 arm --------------------------------------
    v31 <- (k1on * EGF * E1) * Cell - (k1off * EGF_E1) * Cell
    v32 <- (k2on * EGF_E1 * E2) * Cell - (k2off * EGF_E1_E2) * Cell
    v33 <- (k3on *
              (1 - (lap_tumor / s0)^n1 / (km1 + (lap_tumor / s0)^n1)) *
              (1 - (pyr_tumor / s0)^n15 / (km15 + (pyr_tumor / s0)^n15)) *
              EGF_E1_E2 - k3off * EGF_E1_E2p) * Cell
    v34 <- (kd12 * ((lap_tumor / s0)^n2 / (km2 + (lap_tumor / s0)^n2) +
                      (pyr_tumor / s0)^n22 / (km22 + (pyr_tumor / s0)^n22)) *
              EGF_E1_E2) * Cell
    v35 <- (kd12p * EGF_E1_E2p) * Cell

    # ---- HER2 homodimer arm -----------------------------------------
    v36 <- (k4on * E2 * E2) * Cell - (k4off * E2_E2) * Cell
    v37 <- (k5on *
              (1 - (lap_tumor / s0)^n3 / (km3 + (lap_tumor / s0)^n3)) *
              (1 - (pyr_tumor / s0)^n16 / (km16 + (pyr_tumor / s0)^n16)) *
              E2_E2 - k5off * E2_E2p) * Cell
    v38 <- (kd22 * ((lap_tumor / s0)^n4 / (km4 + (lap_tumor / s0)^n4) +
                      (pyr_tumor / s0)^n23 / (km23 + (pyr_tumor / s0)^n23)) *
              E2_E2) * Cell
    v39 <- (kd22p * E2_E2p) * Cell

    # ---- ligand-independent HER2/HER3 arm ---------------------------
    v40 <- (k6on * E2 * E3) * Cell - (k6off * E2_E3) * Cell
    v41 <- (k7on *
              (1 - (lap_tumor / s0)^n5 / (km5 + (lap_tumor / s0)^n5)) *
              (1 - (pyr_tumor / s0)^n17 / (km17 + (pyr_tumor / s0)^n17)) *
              E2_E3 - k7off * E2_E3p) * Cell
    v42 <- (kd23_2 * ((lap_tumor / s0)^n6 / (km6 + (lap_tumor / s0)^n6) +
                        (pyr_tumor / s0)^n24 / (km24 + (pyr_tumor / s0)^n24)) *
              E2_E3) * Cell
    v43 <- (kd23p_2 * E2_E3p) * Cell

    # ---- NRG1 / HER3 / HER2 arm -------------------------------------
    v44 <- (k8on * NRG1 * E3) * Cell - (k8off * NRG1_E3) * Cell
    v45 <- (k9on * NRG1_E3 * E2) * Cell - (k9off * NRG1_E2_E3) * Cell
    v46 <- (k10on *
              (1 - (lap_tumor / s0)^n7 / (km7 + (lap_tumor / s0)^n7)) *
              (1 - (pyr_tumor / s0)^n18 / (km18 + (pyr_tumor / s0)^n18)) *
              NRG1_E2_E3 - k10off * NRG1_E2_E3p) * Cell
    v47 <- (kd23_1 * ((lap_tumor / s0)^n8 / (km8 + (lap_tumor / s0)^n8) +
                        (pyr_tumor / s0)^n25 / (km25 + (pyr_tumor / s0)^n25)) *
              NRG1_E2_E3) * Cell
    v48 <- (kd23p_1 * NRG1_E2_E3p) * Cell

    # ---- NRG1 / HER4 / HER2 arm (lapatinib has no HER4 activity) ----
    v49 <- (k11on * NRG1 * E4) * Cell - (k11off * NRG1_E4) * Cell
    v50 <- (k12on * NRG1_E4 * E2) * Cell - (k12off * NRG1_E2_E4) * Cell
    v51 <- (k13on *
              (1 - (pyr_tumor / s0)^n19 / (km19 + (pyr_tumor / s0)^n19)) *
              NRG1_E2_E4 - k13off * NRG1_E2_E4p) * Cell
    v52 <- (kd24 * ((pyr_tumor / s0)^n20 / (km20 + (pyr_tumor / s0)^n20)) *
              NRG1_E2_E4) * Cell
    v53 <- (kd24p * NRG1_E2_E4p) * Cell

    # ---- Ras / MAPK, PI3K / AKT, feedback and BTK --------------------
    v54 <- (kon_Raf * Raf * (w6 * EGF_E1_E2p + w7 * E2_E2p + w8 * E2_E3p +
                               w9 * NRG1_E2_E3p + w10 * NRG1_E2_E4p)) * Cell
    v55 <- (koff_Raf * Raf_p) * Cell
    v56 <- (kon_ERK * ERK * ((Raf_p / s0)^n9 / (km9 + (Raf_p / s0)^n9))) * Cell
    v57 <- (koff_ERK * ERK_p) * Cell
    v58 <- (kon_PI3K * PI3K * (w1 * EGF_E1_E2p + w2 * E2_E2p + w3 * E2_E3p +
                                 w4 * NRG1_E2_E3p + w5 * NRG1_E2_E4p)) * Cell
    v59 <- (koff_PI3K * PI3K_p) * Cell
    v60 <- (kon_AKT * AKT * ((PI3K_p / s0)^n11 / (km11 + (PI3K_p / s0)^n11)) *
              (1 - ((PTEN / s0)^n12 / (km12 + (PTEN / s0)^n12)))) * Cell
    v61 <- (koff_AKT * AKT_p) * Cell
    v62 <- (kfb1 * AKT_p) * Cell
    v63 <- (kfb2 * feedback) * Cell
    v64 <- (kon_BTK * BTK *
              (w16 * NRG1_E2_E3p / s0 + w17 * NRG1_E2_E4p / s0)^n21 /
              (km21 + (w16 * NRG1_E2_E3p / s0 + w17 * NRG1_E2_E4p / s0)^n21)) * Cell
    v65 <- (koff_BTK * BTK_p) * Cell

    # ---- ADC trafficking in the tumour ------------------------------
    v66 <- (Rx * p3 * adc_central) * Vc
    v67 <- (kon_ADC * ADC_ex * E2) * Cell - (koff_ADC * ADC_E2) * Cell
    v68 <- (kint_ADC * ADC_E2) * Cell
    v69 <- (kdeg_ADC_1 * ADC_E2_int) * Cell
    v70 <- (kdeg_ADC_2 * dar * ADC_E2_int) * Cell
    v71 <- (kout_PL * PL_cell) * Cell
    v72 <- (kper_PL * PL_cell) * Cell - (kper_PL * PL_ex) * Cell
    v73 <- (kdec_ADC_3 * ADC_ex) * Cell
    v74 <- (kdec_ADC_4 * dar * ADC_ex) * Cell
    v75 <- (kdec_ADC_3 * ADC_E2) * Cell
    v76 <- (kdec_ADC_4 * dar * ADC_E2) * Cell
    v77 <- (kcl_PL * PL_ex) * Cell

    # ---- tumour growth and death ------------------------------------
    growthDrive <- w11 * ERK_p / s0 + w12 * AKT_p / s0 + w13 * BTK_p / s0
    v78 <- umax * growthDrive^n13 / (km13 + growthDrive^n13) *
      (1 - Cells / Cellsmax) * Cells
    v79 <- dmax * (1 +
                     w14 * (PL_cell / s0)^n14 / (km14 + (PL_cell / s0)^n14) +
                     w15 * (fu_tumor / s0)^n10 / (km10 + (fu_tumor / s0)^n10)) * Cells

    # =================================================================
    # ODE system (Table S1 sheet Equations)
    # =================================================================
    d/dt(adc_central)    <- 1 / Vc * (-v1 - v5 - v7 - v66)
    d/dt(mab_central)    <- 1 / Vc * (-v2 - v9 + v7)
    d/dt(adc_peripheral) <- 1 / Vp_ADC * (v1)
    d/dt(mab_peripheral) <- 1 / Vp_ADC * (v2)
    d/dt(pl_central)     <- 1 / Vc_PL * (-v3 - v4 + v6 + v8)
    d/dt(pl_peripheral)  <- 1 / Vp_PL * (v3)
    d/dt(dar)            <- -v10
    d/dt(cap_depot)      <- 1 / V_cap * (-v11)
    d/dt(cap_central)    <- 1 / V_cap * (v11 - v12 - v16)
    d/dt(dfcr)           <- 1 / V_met * (v12 - v13)
    d/dt(dfur)           <- 1 / V_met * (v13 - v14)
    d/dt(fu)             <- 1 / V_met * (v14 - v15)
    d/dt(lap_depot)      <- 1 / Vc * (-v17)
    d/dt(lap_central)    <- 1 / Vc * (v17 - v19 - v18)
    d/dt(lap_peripheral) <- 1 / Vp_lap * (v18)
    d/dt(pyr_depot)      <- 1 / Vc * (-v20)
    d/dt(pyr_central)    <- 1 / Vc * (v20 - v21 - v22)
    d/dt(pyr_peripheral) <- 1 / Vp_pyr * (v21)

    d/dt(EGF)  <- 1 / Cell * (-v31)
    d/dt(NRG1) <- 1 / Cell * (-v44 - v49)
    d/dt(E1)   <- 1 / Cell * (-v31 + v23 - v24)
    d/dt(E2)   <- 1 / Cell * (-v32 - 2 * v36 - v40 - v45 - v50 - v67 + v25 - v26)
    d/dt(E3)   <- 1 / Cell * (-v40 - v44 + v27 - v28)
    d/dt(E4)   <- 1 / Cell * (-v49 + v29 - v30)

    d/dt(EGF_E1)     <- 1 / Cell * (v31 - v32)
    d/dt(EGF_E1_E2)  <- 1 / Cell * (v32 - v33 - v34)
    d/dt(EGF_E1_E2p) <- 1 / Cell * (v33 - v35)
    d/dt(E2_E2)      <- 1 / Cell * (v36 - v37 - v38)
    d/dt(E2_E2p)     <- 1 / Cell * (v37 - v39)
    d/dt(E2_E3)      <- 1 / Cell * (v40 - v41 - v42)
    d/dt(E2_E3p)     <- 1 / Cell * (v41 - v43)
    d/dt(NRG1_E3)     <- 1 / Cell * (v44 - v45)
    d/dt(NRG1_E2_E3)  <- 1 / Cell * (v45 - v46 - v47)
    d/dt(NRG1_E2_E3p) <- 1 / Cell * (v46 - v48)
    d/dt(NRG1_E4)     <- 1 / Cell * (v49 - v50)
    d/dt(NRG1_E2_E4)  <- 1 / Cell * (v50 - v51 - v52)
    d/dt(NRG1_E2_E4p) <- 1 / Cell * (v51 - v53)

    d/dt(Raf)      <- 1 / Cell * (-v54 + v55)
    d/dt(Raf_p)    <- 1 / Cell * (v54 - v55)
    d/dt(ERK)      <- 1 / Cell * (-v56 + v57)
    d/dt(ERK_p)    <- 1 / Cell * (v56 - v57)
    d/dt(PI3K)     <- 1 / Cell * (-v58 + v59)
    d/dt(PI3K_p)   <- 1 / Cell * (v58 - v59)
    d/dt(AKT)      <- 1 / Cell * (-v60 + v61)
    d/dt(AKT_p)    <- 1 / Cell * (v60 - v61)
    d/dt(feedback) <- 1 / Cell * (v62 - v63)
    d/dt(BTK)      <- 1 / Cell * (-v64 + v65)
    d/dt(BTK_p)    <- 1 / Cell * (v64 - v65)

    d/dt(ADC_ex)     <- 1 / Cell * (-v67 - v73 + v66)
    d/dt(ADC_E2)     <- 1 / Cell * (v67 - v68 - v75)
    d/dt(ADC_E2_int) <- 1 / Cell * (v68 - v69)
    d/dt(PL_cell)    <- 1 / Cell * (-v71 + v70 - v72)
    d/dt(PL_ex)      <- 1 / Cell * (v71 + v76 + v74 - v77 + v72)
    d/dt(Cells)      <- v78 - v79

    # =================================================================
    # Initial conditions (Table S1 sheet Species, SKBR3 / in vivo)
    # The listed values are the drug-free steady state the authors
    # reached by pre-simulating the ligand-free network.
    # =================================================================
    dar(0)      <- dar0
    E1(0)       <- E1_0
    E2(0)       <- E2_0
    E3(0)       <- E3_0
    E4(0)       <- E4_0
    E2_E2(0)    <- E2_E2_0
    E2_E2p(0)   <- E2_E2p_0
    E2_E3(0)    <- E2_E3_0
    E2_E3p(0)   <- E2_E3p_0
    Raf(0)      <- Raf_0
    Raf_p(0)    <- Rafp_0
    ERK(0)      <- ERK_0
    ERK_p(0)    <- ERKp_0
    PI3K(0)     <- PI3K_0
    PI3K_p(0)   <- PI3Kp_0
    AKT(0)      <- AKT_0
    AKT_p(0)    <- AKTp_0
    feedback(0) <- fb_0
    BTK(0)      <- BTK_0
    Cells(0)    <- cells0

    # =================================================================
    # Readouts.  The paper reports no residual-error model, no IIV and
    # no parameter standard errors, so the model is deterministic.
    # =================================================================
    Cc        <- adc_central   # conjugated ADC in plasma (nmol/L)
    Cc_tab    <- TTmAb         # total antibody, ADC plus naked mAb (nmol/L)
    tumorSize <- Cells         # tumour volume in vivo (mm3) or relative cell number in vitro
  })
}
