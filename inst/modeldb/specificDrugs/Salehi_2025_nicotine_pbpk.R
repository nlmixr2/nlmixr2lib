Salehi_2025_nicotine_pbpk <- function() {
  description <- "PBPK (whole-body, MCSim/deSolve). Nicotine, cotinine and nicotine glucuronide disposition in adults who use tobacco-free nicotine pouches, moist smokeless tobacco (MST) or intravenous nicotine (Salehi 2025). A buccal-cavity tissue permeation front-end solves Fick's second law across the buccal epithelium as a 20-slab method-of-lines diffusion chain (buccal_slab1..buccal_slab20) with a released-nicotine flux boundary at the saliva interface and a perfect-sink boundary at the effective tissue depth; the sink-face flux is the blood uptake rate feeding the perfused buccal submucosa, and the balance of released nicotine that is swallowed enters the gut. The whole-body disposition model is inherited from Rostami 2022 with hepatic and urinary clearances reduced by 50 percent per Salehi 2025, and carries eleven flow-limited nicotine tissues plus a nicotine-driven heart-rate feedback on cardiac output, a five-tissue cotinine sub-model and a one-compartment nicotine-glucuronide sub-model. Deterministic typical-value model: the sources report no between-subject variability and no residual-error model. The inhalation (cigarette) route of the parent Rostami model is NOT implemented because its CFD-derived deposition fractions are unreported in every on-disk source; see the vignette Errata."
  reference <- paste(
    "Salehi A, Sarkar MA, Smith JH, Rostami AA. (2025).",
    "Physiologically Based Pharmacokinetic Modeling to Predict Nicotine Pharmacokinetics of Nicotine Pouches Under Naturalistic Use Conditions.",
    "J Clin Pharmacol 65(10):1297-1309. doi:10.1002/jcph.70038. PMCID PMC12484417.",
    "Correction added 2025-06-27 (final sentence of Methods; prose only, no parameter impact).",
    "Whole-body disposition structure and physiological / chemical parameters inherited from",
    "Rostami AA, Campbell JL, Pithawalla YB, et al. (2022) Sci Rep 12:2436, doi:10.1038/s41598-022-06209-4",
    "(with Author Correction, Sci Rep 12:2966, doi:10.1038/s41598-022-07016-7 -- bibliography fix only).",
    sep = " "
  )
  vignette <- "Salehi_2025_nicotine_pbpk"

  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ng/mL",
    amount        = "mg",
    weight        = "kg"
  )

  # Nicotine is released from the product into the mouth (`depot`) and is
  # dosed there; the intravenous validation arm doses arterial blood
  # directly. Both sit outside the depot/central default, so the dosing
  # targets are declared explicitly.
  dosing <- c("depot", "a_arterial")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Scales every tissue volume, blood flow and clearance.",
        "Rostami 2022 Table 1 uses a 73.0 kg reference adult;",
        "Salehi 2025 Figure 2 reports study-matched mean body weights of",
        "~92 kg (TP1 / on! cohort, which also supplied the 2 g MST data) and",
        "~97 kg (TP2 / on! PLUS cohort, which also supplied the 4 g MST data)."
      ),
      source_name        = "BW"
    )
  )

  # Every ODE state. `analyte` / `specimen` checked against the Salehi 2025
  # supplementary MCSim listing and Rostami 2022 Tables 1-2.
  compartmentData <- list(
    depot                  = list(analyte = "nicotine", units = "mg", specimen = "administration site", verified = TRUE),
    buccal_slab1           = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    buccal_slab2           = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    buccal_slab3           = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    buccal_slab4           = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    buccal_slab5           = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    buccal_slab6           = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    buccal_slab7           = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    buccal_slab8           = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    buccal_slab9           = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    buccal_slab10          = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    buccal_slab11          = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    buccal_slab12          = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    buccal_slab13          = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    buccal_slab14          = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    buccal_slab15          = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    buccal_slab16          = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    buccal_slab17          = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    buccal_slab18          = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    buccal_slab19          = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    buccal_slab20          = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    a_buccal               = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    a_conducting_airway    = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    a_transitional_airway  = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    a_pulmonary            = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    a_venous               = list(analyte = "nicotine", units = "mg", specimen = "whole blood", verified = TRUE),
    a_arterial             = list(analyte = "nicotine", units = "mg", specimen = "whole blood", verified = TRUE),
    a_heart                = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    a_brain                = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    a_liver                = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    a_skin                 = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    a_muscle               = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    a_fat                  = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    a_rapidly_perfused     = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    a_slowly_perfused      = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    a_gut                  = list(analyte = "nicotine", units = "mg", specimen = "administration site", verified = TRUE),
    gut_lumen              = list(analyte = "nicotine", units = "mg", specimen = "faeces", verified = TRUE),
    a_hepatic              = list(analyte = "nicotine", units = "mg", specimen = "tissue", verified = TRUE),
    a_urine                = list(analyte = "nicotine", units = "mg", specimen = "urine", verified = TRUE),
    effect                 = list(analyte = "nicotine", units = "mg/L", specimen = "not applicable", verified = TRUE),
    a_venous_cot           = list(analyte = "cotinine", units = "mg", specimen = "whole blood", verified = TRUE),
    a_arterial_cot         = list(analyte = "cotinine", units = "mg", specimen = "whole blood", verified = TRUE),
    a_liver_cot            = list(analyte = "cotinine", units = "mg", specimen = "tissue", verified = TRUE),
    a_muscle_cot           = list(analyte = "cotinine", units = "mg", specimen = "tissue", verified = TRUE),
    a_fat_cot              = list(analyte = "cotinine", units = "mg", specimen = "tissue", verified = TRUE),
    a_rapidly_perfused_cot = list(analyte = "cotinine", units = "mg", specimen = "tissue", verified = TRUE),
    a_slowly_perfused_cot  = list(analyte = "cotinine", units = "mg", specimen = "tissue", verified = TRUE),
    a_hepatic_cot          = list(analyte = "cotinine", units = "mg", specimen = "tissue", verified = TRUE),
    a_urine_cot            = list(analyte = "cotinine", units = "mg", specimen = "urine", verified = TRUE),
    central_gluc           = list(analyte = "nicotine glucuronide", units = "mg", specimen = "plasma", verified = TRUE),
    a_urine_gluc           = list(analyte = "nicotine glucuronide", units = "mg", specimen = "urine", verified = TRUE)
  )

  population <- list(
    species       = "human",
    age_range     = "adults who use tobacco products (21+ years)",
    weight_median = "92 kg (on! / TP1 cohort) and 97 kg (on! PLUS / TP2 cohort); Rostami 2022 reference adult 73.0 kg",
    disease_state = "healthy adults who smoke cigarettes or use moist smokeless tobacco",
    dose_range    = paste(
      "on! nicotine pouches 2 / 4 / 8 mg (30 min controlled use);",
      "on! PLUS nicotine pouches 6 / 9 / 12 mg (45 min controlled use);",
      "own-brand MST 2 g (30 min) and 4 g (45 min) at ~10 mg nicotine/g;",
      "intravenous nicotine 2 ug/kg/min over 30 min (Figure S3 validation arm)"
    ),
    regions       = "United States",
    notes         = paste(
      "The PBPK model is not fitted to individual-level data. Clinical PK profiles used",
      "for the base-use-case regressions are cohort means adapted from Liu 2022, Rensch 2021",
      "and associated conference presentations (Salehi 2025 'Clinical Data');",
      "naturalistic use patterns are from the Becker 2024 actual-use study",
      "(median ~6 pouches/day, ~12-13 min in mouth) and the PATH study (~15 cigarettes/day).",
      "Mean nicotine extraction was ~61% (SD ~21%) for on! 8 mg over 30 min and",
      "67-72% (SD 28-32%) for on! PLUS over 45 min."
    )
  )

  ini({
    # =====================================================================
    # BUCCAL-CAVITY TISSUE PERMEATION MODEL (Salehi 2025 Equations 1, 2, 5)
    # The (D, l) pair is treated by the authors as a property of buccal
    # tissue itself -- independent of product and formulation -- and is
    # fixed once optimized. Both were obtained by regression against 2 g
    # MST clinical PK plus the Rostami 2022 tissue-uptake measurement, so
    # they are estimated quantities, not literature constants.
    # =====================================================================
    DBUC   <- 1.2e-5   ; label("Nicotine diffusivity in buccal tissue (cm^2/s)")                    # Salehi 2025 Results 'Tissue Permeation Model Parameters': Dopt = 1.2e-5 cm^2/s
    LBUC   <- 0.175    ; label("Effective buccal tissue thickness, perfect-sink depth (cm)")        # Salehi 2025 Results: l_opt = 1.75 mm

    # =====================================================================
    # PRODUCT-SPECIFIC RELEASE AND UPTAKE (Salehi 2025 Equations 2, 6;
    # Table 2 'Adjusted parameter'). Defaults below are the on! (TP1) 8 mg
    # 30-minute base use case. Change these five for another product /
    # nicotine level / use scenario -- see the model vignette.
    # =====================================================================
    KREL   <- 1.8832   ; label("First-order nicotine release rate constant from the product (1/h)") # Salehi 2025 Methods 'Clinical Data': ~61% of on! 8 mg extracted over 30 min, exponential interpolation A0*(1-exp(-k*t)) => k = -log(1-0.61)/0.5 h = 1.8832
    TUSE   <- 0.5      ; label("Product use duration, time held in the mouth (h)")                  # Salehi 2025 Table 2 / Figure 2: on! base use case is 30 min
    FTISM  <- 0.46     ; label("Mean buccal tissue uptake fraction of released nicotine (unitless)")# Salehi 2025 Table 2, on! 8 mg row: f_tissue-bar = 0.46
    FTISE  <- 0.46     ; label("Buccal tissue uptake fraction at end of use (unitless)")            # Salehi 2025 Table 2: f_e reported only for on! PLUS; constant-uptake products set f_e = f_tissue-bar so Equation 6 collapses to the constant model
    FSWAL  <- fixed(1) ; label("Fraction of non-absorbed released nicotine that is swallowed (unitless)") # Salehi 2025 Methods: pouches are not expectorated so the balance transfers to the GI tract; set to 0 for MST, where 'transfer to the GI tract during MST use is negligible due to expectoration'

    # =====================================================================
    # GASTROINTESTINAL UPTAKE
    # =====================================================================
    KA1    <- 0.8      ; label("Gastrointestinal uptake rate constant (1/h)")                       # Salehi 2025 Results 'GI Tract Uptake Rate Constant (KA1)': optimized to 0.8 /h (Teeguarden 2013 reported 1.34 /h for a combined BC+GI pseudo-compartment)
    FA     <- fixed(0.67); label("Oral bioavailability of swallowed nicotine (unitless)")           # Rostami 2022 Table 2: FA = 0.67 (Teeguarden et al.)

    # =====================================================================
    # METABOLISM AND EXCRETION
    # Salehi 2025 Methods 'Nicotine PBPK Model': "the hepatic and urinary
    # clearance and metabolic parameters adjusted (reduced by 50%) to
    # improve the predictions of the terminal phase". Each value below is
    # therefore half of the Rostami 2022 Table 2 entry.
    # =====================================================================
    CLMC   <- fixed(1.35)  ; label("Hepatic nicotine metabolic clearance constant (L/h/kg)")                 # Rostami 2022 Table 2 CLMC = 2.70 L/h/kg^0.75, halved per Salehi 2025. NOTE: applied linearly in body weight, not BW^0.75 -- see model() and vignette Errata
    CLKC   <- fixed(0.21)  ; label("Urinary nicotine clearance constant (L/h/kg^0.75)")                      # Rostami 2022 Table 2 CLKC = 0.42, halved per Salehi 2025
    CLLMC  <- fixed(0.07)  ; label("Hepatic cotinine metabolic clearance constant (L/h/kg^0.75)")            # Rostami 2022 Table 2 CLLMC = 0.14, halved per Salehi 2025
    CLKMC  <- fixed(0.0125); label("Urinary cotinine clearance constant (L/h/kg^0.75)")                      # Rostami 2022 Table 2 CLKMC = 0.025, halved per Salehi 2025
    FNC    <- fixed(0.80)  ; label("Fraction of metabolized nicotine converted to cotinine (unitless)")      # Rostami 2022 Table 2 FNC = 0.80 (the supplementary MCSim listings of both papers instead carry FNC = 0.72; the published table value is adopted -- see vignette Errata)
    FNG    <- fixed(0.10)  ; label("Fraction of metabolized nicotine converted to nicotine glucuronide (unitless)") # Salehi 2025 supplementary MCSim listing: FNG = 0.10
    VGDC   <- fixed(1.0)   ; label("Nicotine glucuronide distribution volume constant (L/kg)")               # Salehi 2025 supplementary MCSim listing: VGDC = 1.0
    GCLC   <- fixed(0.1)   ; label("Urinary nicotine glucuronide clearance constant (L/h/kg^0.75)")          # Salehi 2025 supplementary MCSim listing: GCLC = 0.1

    # =====================================================================
    # PHYSIOLOGY -- Rostami 2022 Table 1
    # =====================================================================
    QCC    <- fixed(16)      ; label("Cardiac output constant (L/h/kg^0.75)")                # Rostami 2022 Table 1 QCC = 16 (ICRP)
    HRO    <- fixed(61.1)    ; label("Basal heart rate (beats/min)")                         # Rostami 2022 Table 1 HRO = 61.1

    VHC    <- fixed(0.0044)  ; label("Heart volume (fraction of body weight)")               # Rostami 2022 Table 1 VHC
    VBC    <- fixed(0.02)    ; label("Brain volume (fraction of body weight)")               # Rostami 2022 Table 1 VBC
    VFC    <- fixed(0.258)   ; label("Fat volume (fraction of body weight)")                 # Rostami 2022 Table 1 VFC
    VLC    <- fixed(0.024)   ; label("Liver volume (fraction of body weight)")               # Rostami 2022 Table 1 VLC
    VSKC   <- fixed(0.042)   ; label("Skin volume (fraction of body weight)")                # Rostami 2022 Table 1 VSKC
    VMC    <- fixed(0.34)    ; label("Muscle volume (fraction of body weight)")              # Rostami 2022 Table 1 VMC
    VABC   <- fixed(0.02)    ; label("Arterial blood volume (fraction of body weight)")      # Rostami 2022 Table 1 VABC
    VVBC   <- fixed(0.05)    ; label("Venous blood volume (fraction of body weight)")        # Rostami 2022 Table 1 VVBC
    VRC    <- fixed(0.03)    ; label("Rapidly perfused volume (fraction of body weight)")    # Rostami 2022 Table 1 VRC
    VSLOWC <- fixed(0.08)    ; label("Slowly perfused volume (fraction of body weight)")     # Rostami 2022 Table 1 VSLOWC

    QFC    <- fixed(0.068)   ; label("Fat blood flow (fraction of cardiac output)")          # Rostami 2022 Table 1 QFC
    QBC    <- fixed(0.12)    ; label("Brain blood flow (fraction of cardiac output)")        # Rostami 2022 Table 1 QBC
    QHC    <- fixed(0.04)    ; label("Heart blood flow (fraction of cardiac output)")        # Rostami 2022 Table 1 QHC
    QSKC   <- fixed(0.05)    ; label("Skin blood flow (fraction of cardiac output)")         # Rostami 2022 Table 1 QSKC
    QMC    <- fixed(0.14)    ; label("Muscle blood flow (fraction of cardiac output)")       # Rostami 2022 Table 1 QMC
    QLC    <- fixed(0.26)    ; label("Liver blood flow (fraction of cardiac output)")        # Rostami 2022 Table 1 QLC
    QRC    <- fixed(0.19)    ; label("Rapidly perfused blood flow (fraction of cardiac output)")  # Rostami 2022 Table 1 QRC
    QSC    <- fixed(0.08)    ; label("Slowly perfused blood flow (fraction of cardiac output)")   # Rostami 2022 Table 1 QSC
    QBUC   <- fixed(0.0215)  ; label("Buccal cavity blood flow (fraction of cardiac output)")     # Rostami 2022 Table 1 QBUC (Corley et al.)
    QCAC   <- fixed(0.025)   ; label("Conducting airway blood flow (fraction of cardiac output)") # Rostami 2022 Table 1 QCAC (Corley et al.)
    QTAC   <- fixed(0.007)   ; label("Transitional airway blood flow (fraction of cardiac output)")# Rostami 2022 Table 1 QTAC (Corley et al.)

    SABU   <- fixed(103.10)  ; label("Buccal cavity surface area (cm^2)")                    # Rostami 2022 Table 1 SABU
    SACA   <- fixed(199.50)  ; label("Conducting airway surface area (cm^2)")                # Rostami 2022 Table 1 SACA
    SATA   <- fixed(163.60)  ; label("Transitional airway surface area (cm^2)")              # Rostami 2022 Table 1 SATA
    SAPUL  <- fixed(540000)  ; label("Pulmonary surface area (cm^2)")                        # Rostami 2022 Table 1 SAPUL
    WTBU   <- fixed(0.0065)  ; label("Buccal cavity epithelium width (cm)")                  # Rostami 2022 Table 1 WTBU
    WTCA   <- fixed(0.0065)  ; label("Conducting airway epithelium width (cm)")              # Rostami 2022 Table 1 WTCA
    WTTA   <- fixed(0.0065)  ; label("Transitional airway epithelium width (cm)")            # Rostami 2022 Table 1 WTTA
    WTPUL  <- fixed(0.000036); label("Pulmonary epithelium width (cm)")                      # Rostami 2022 Table 1 WTPUL

    # =====================================================================
    # PARTITION COEFFICIENTS -- Rostami 2022 Table 2
    # =====================================================================
    PLU    <- fixed(0.90)  ; label("Nicotine lung:blood partition coefficient (unitless)")            # Rostami 2022 Table 2 PLU
    PF     <- fixed(0.80)  ; label("Nicotine fat:blood partition coefficient (unitless)")             # Rostami 2022 Table 2 PF
    PBR    <- fixed(3.00)  ; label("Nicotine brain:blood partition coefficient (unitless)")           # Rostami 2022 Table 2 PBR
    PL     <- fixed(7.50)  ; label("Nicotine liver:blood partition coefficient (unitless)")           # Rostami 2022 Table 2 PL
    PH     <- fixed(1.60)  ; label("Nicotine heart:blood partition coefficient (unitless)")           # Rostami 2022 Table 2 PH
    PSK    <- fixed(1.50)  ; label("Nicotine skin:blood partition coefficient (unitless)")            # Rostami 2022 Table 2 PSK
    PM     <- fixed(1.50)  ; label("Nicotine muscle:blood partition coefficient (unitless)")          # Rostami 2022 Table 2 PM
    PR     <- fixed(7.50)  ; label("Nicotine rapidly-perfused:blood partition coefficient (unitless)")# Rostami 2022 Table 2 PR (set to liver)
    PS     <- fixed(1.50)  ; label("Nicotine slowly-perfused:blood partition coefficient (unitless)") # Rostami 2022 Table 2 PS (set to muscle)
    PML    <- fixed(2.00)  ; label("Cotinine liver:blood partition coefficient (unitless)")           # Rostami 2022 Table 2 PML
    PMM    <- fixed(1.50)  ; label("Cotinine muscle:blood partition coefficient (unitless)")          # Rostami 2022 Table 2 PMM
    PMR    <- fixed(1.50)  ; label("Cotinine rapidly-perfused:blood partition coefficient (unitless)")# Rostami 2022 Table 2 PMR
    PMS    <- fixed(1.00)  ; label("Cotinine slowly-perfused:blood partition coefficient (unitless)") # Rostami 2022 Table 2 PMS
    PMF    <- fixed(0.50)  ; label("Cotinine fat:blood partition coefficient (unitless)")             # Rostami 2022 Table 2 PMF

    # =====================================================================
    # SATURABLE TISSUE BINDING -- Salehi 2025 supplementary MCSim listing.
    # Both affinities are zero in the published code, which degenerates the
    # Langmuir isotherm to a hard-threshold sink of capacity BM*V. The
    # capacities are ~1.6e-3 mg (heart) and ~5e-5 mg (lung), i.e. well
    # under 0.1% of a 2-12 mg product dose. Implemented as coded.
    # =====================================================================
    BMLURC <- fixed(0.00235); label("Lung nicotine binding capacity per unit tissue volume (mg/L)")   # Salehi 2025 supplementary MCSim listing: BMLURC = 0.00235
    BMHRC  <- fixed(0.00427); label("Heart nicotine binding capacity per unit tissue volume (mg/L)")  # Salehi 2025 supplementary MCSim listing: BMHRC = 0.00427
    KBLU   <- fixed(0)      ; label("Lung nicotine binding affinity (mg/L)")                          # Salehi 2025 supplementary MCSim listing: KBLU = 0
    KBH    <- fixed(0)      ; label("Heart nicotine binding affinity (mg/L)")                         # Salehi 2025 supplementary MCSim listing: KBH = 0

    # =====================================================================
    # HEART-RATE / CARDIAC-OUTPUT FEEDBACK -- Rostami 2022 Table 2
    # =====================================================================
    SPD    <- fixed(933.66); label("Nicotine concentration-heart rate slope (beats/min per mg/L)")   # Rostami 2022 Table 2 S = 933.66 (Teeguarden et al.)
    KANT   <- fixed(1.6617); label("First-order rate of loss of nicotine tolerance (1/h)")           # Rostami 2022 Table 2 KANT = 1.6617
    CANT50 <- fixed(0.0152); label("Nicotine tolerance half-maximal concentration (mg/L)")           # Rostami 2022 Table 2 CANT50 = 0.0152

    # =====================================================================
    # MOLECULAR WEIGHTS (physical constants; declared as 0 placeholders in
    # the published MCSim listing, whose .in input files were not published)
    # =====================================================================
    MW     <- fixed(162.23); label("Nicotine molecular weight (g/mol)")               # physical constant; PubChem CID 89594. Declared MW = 0 in the published MCSim listing
    MWCOT  <- fixed(178)   ; label("Cotinine molecular weight (g/mol)")               # Salehi 2025 supplementary MCSim listing: hard-coded literal 178 in the cotinine formation term
    MWGLU  <- fixed(338)   ; label("Nicotine glucuronide molecular weight (g/mol)")   # Salehi 2025 supplementary MCSim listing: hard-coded literal 338 in the glucuronide formation term (the MWG = 338.4 declaration is not the value used)
  })

  model({
    # =====================================================================
    # 1. TISSUE VOLUMES (L)
    #    Epithelial volumes are surface area x width, converted cm^3 -> L.
    #    Rostami 2022 renormalizes every body-weight fraction by VTOTALC so
    #    the explicit epithelial volumes do not double-count body mass.
    #    The nasal and perfused-lung compartments of the parent model are
    #    omitted: their blood flows QNC and QLNGC are absent from Rostami
    #    2022 Table 1 (whose flow fractions already sum to ~1.0) and are 0
    #    in the published code, so both compartments are inert.
    # =====================================================================
    vca     <- SACA * WTCA / 1000
    vta     <- SATA * WTTA / 1000
    vpul    <- SAPUL * WTPUL / 1000
    vbu     <- SABU * WTBU / 1000
    vlu     <- vca + vta + vpul
    vtotalc <- (vlu + vbu) / WT + VHC + VBC + VFC + VSKC + VMC + VABC + VVBC + VSLOWC + VRC + VLC
    vf      <- VFC * WT / vtotalc
    vb      <- VBC * WT / vtotalc
    vh      <- VHC * WT / vtotalc
    vr      <- VRC * WT / vtotalc
    vl      <- VLC * WT / vtotalc
    vsk     <- VSKC * WT / vtotalc
    vm      <- VMC * WT / vtotalc
    vs      <- VSLOWC * WT / vtotalc
    vab     <- VABC * WT / vtotalc
    vvb     <- VVBC * WT / vtotalc
    vgd     <- VGDC * WT
    bmh     <- BMHRC * vh
    bmlu    <- BMLURC * vlu

    # =====================================================================
    # 2. SCALED CLEARANCES (L/h)
    #    CLM is assigned twice in the published MCSim listing -- first as
    #    CLMC*BW^0.75, then four lines later as CLMC*BW. MCSim's generated
    #    C keeps the LAST assignment, so hepatic metabolic clearance is
    #    linear in body weight while every other clearance is allometric.
    #    Reproduced as coded; see vignette Errata.
    # =====================================================================
    clk  <- CLKC * WT^0.75
    clm  <- CLMC * WT
    clkm <- CLKMC * WT^0.75
    cllm <- CLLMC * WT^0.75
    gcl  <- GCLC * WT^0.75

    # =====================================================================
    # 3. CARDIAC OUTPUT WITH NICOTINE HEART-RATE FEEDBACK
    #    Rostami 2022: heart rate rises with venous nicotine and is damped
    #    by an acquired-tolerance state; cardiac output is heart rate x
    #    stroke volume, and every tissue flow is a fixed fraction of it.
    # =====================================================================
    qci     <- QCC * WT^0.75
    sv      <- qci / (HRO * 60)
    cvvb    <- a_venous / vvb
    eout    <- HRO + SPD * cvvb / (1 + effect / CANT50)
    qc      <- eout * 60 * sv
    qtotalc <- QFC + QBC + QRC + QLC + QHC + QSKC + QMC + QSC + QBUC + QTAC + QCAC
    qf   <- QFC * qc / qtotalc
    qb   <- QBC * qc / qtotalc
    qr   <- QRC * qc / qtotalc
    ql   <- QLC * qc / qtotalc
    qh   <- QHC * qc / qtotalc
    qsk  <- QSKC * qc / qtotalc
    qm   <- QMC * qc / qtotalc
    qs   <- QSC * qc / qtotalc
    qbu  <- QBUC * qc / qtotalc
    qca  <- QCAC * qc / qtotalc
    qta  <- QTAC * qc / qtotalc

    # Tolerance state: rises toward the current venous concentration and is
    # rectified so tolerance is never lost faster than it is gained
    # (published code: DCANT = (KANT*(CVVB-CANT) < 0) ? 0 : KANT*(CVVB-CANT)).
    d/dt(effect) <- max(KANT * (cvvb - effect), 0)

    # =====================================================================
    # 4. PRODUCT RELEASE AND BUCCAL UPTAKE PARTITION
    #    Salehi 2025 approximates release by differentiating the exponential
    #    interpolation A0*(1 - exp(-k*t)) of pre- and post-use product
    #    nicotine content, which is exactly first-order loss from `depot`.
    #    The product is removed from the mouth at the end of the use period
    #    by zeroing `depot` with a replacement event (evid = 5) in the event
    #    table, so any nicotine still unreleased leaves with the product.
    #
    #    Equation 6 (linear tissue uptake) is written exactly as printed.
    #    It requires the time elapsed since the product entered the mouth,
    #    so the event table MUST carry a `depot` record; for a pure
    #    intravenous simulation include a placeholder `amt = 0` depot dose
    #    at time 0, otherwise tad(depot) is NA and every state becomes NA.
    #    For constant-uptake products FTISE == FTISM and the ramp term is
    #    algebraically zero.
    # =====================================================================
    rrel <- KREL * depot
    ftis <- ((FTISE - FTISM) / 2) * (tad(depot) / TUSE) + (FTISE + FTISM) / 2
    jin  <- ftis * rrel
    gitr <- FSWAL * (1 - ftis) * rrel
    d/dt(depot) <- -rrel

    # =====================================================================
    # 5. BUCCAL TISSUE DIFFUSION (Salehi 2025 Equations 1, 2, 5)
    #    Fick's second law dC/dt = D d2C/dx2 across the epithelium, solved
    #    as a mass-conserving finite-volume method of lines over 20 equal
    #    slabs. Slab 1 receives the tissue-uptake flux at the saliva
    #    interface (Equation 2); the perfect-sink boundary C(t, l) = 0 at
    #    the far face gives the outflow `jout` over a half-cell distance,
    #    and that flux is the blood uptake rate R_blood of Equation 5.
    #    AUC is converged to four decimal places from 10 slabs upward;
    #    20 slabs puts Cmax within 0.03% of an 80-slab solution.
    #    Slab volumes are cm^3 = mL and slab concentrations mg/mL, so the
    #    conductance `dcond` (cm^2/h x cm^2 / cm = mL/h) yields mg/h fluxes.
    # =====================================================================
    hstep <- LBUC / 20
    vslab <- SABU * hstep
    dcond <- (DBUC * 3600) * SABU / hstep

    cbu1  <- buccal_slab1 / vslab
    cbu2  <- buccal_slab2 / vslab
    cbu3  <- buccal_slab3 / vslab
    cbu4  <- buccal_slab4 / vslab
    cbu5  <- buccal_slab5 / vslab
    cbu6  <- buccal_slab6 / vslab
    cbu7  <- buccal_slab7 / vslab
    cbu8  <- buccal_slab8 / vslab
    cbu9  <- buccal_slab9 / vslab
    cbu10 <- buccal_slab10 / vslab
    cbu11 <- buccal_slab11 / vslab
    cbu12 <- buccal_slab12 / vslab
    cbu13 <- buccal_slab13 / vslab
    cbu14 <- buccal_slab14 / vslab
    cbu15 <- buccal_slab15 / vslab
    cbu16 <- buccal_slab16 / vslab
    cbu17 <- buccal_slab17 / vslab
    cbu18 <- buccal_slab18 / vslab
    cbu19 <- buccal_slab19 / vslab
    cbu20 <- buccal_slab20 / vslab

    jf1  <- dcond * (cbu1 - cbu2)
    jf2  <- dcond * (cbu2 - cbu3)
    jf3  <- dcond * (cbu3 - cbu4)
    jf4  <- dcond * (cbu4 - cbu5)
    jf5  <- dcond * (cbu5 - cbu6)
    jf6  <- dcond * (cbu6 - cbu7)
    jf7  <- dcond * (cbu7 - cbu8)
    jf8  <- dcond * (cbu8 - cbu9)
    jf9  <- dcond * (cbu9 - cbu10)
    jf10 <- dcond * (cbu10 - cbu11)
    jf11 <- dcond * (cbu11 - cbu12)
    jf12 <- dcond * (cbu12 - cbu13)
    jf13 <- dcond * (cbu13 - cbu14)
    jf14 <- dcond * (cbu14 - cbu15)
    jf15 <- dcond * (cbu15 - cbu16)
    jf16 <- dcond * (cbu16 - cbu17)
    jf17 <- dcond * (cbu17 - cbu18)
    jf18 <- dcond * (cbu18 - cbu19)
    jf19 <- dcond * (cbu19 - cbu20)
    jout <- 2 * dcond * cbu20

    d/dt(buccal_slab1)  <- jin - jf1
    d/dt(buccal_slab2)  <- jf1 - jf2
    d/dt(buccal_slab3)  <- jf2 - jf3
    d/dt(buccal_slab4)  <- jf3 - jf4
    d/dt(buccal_slab5)  <- jf4 - jf5
    d/dt(buccal_slab6)  <- jf5 - jf6
    d/dt(buccal_slab7)  <- jf6 - jf7
    d/dt(buccal_slab8)  <- jf7 - jf8
    d/dt(buccal_slab9)  <- jf8 - jf9
    d/dt(buccal_slab10) <- jf9 - jf10
    d/dt(buccal_slab11) <- jf10 - jf11
    d/dt(buccal_slab12) <- jf11 - jf12
    d/dt(buccal_slab13) <- jf12 - jf13
    d/dt(buccal_slab14) <- jf13 - jf14
    d/dt(buccal_slab15) <- jf14 - jf15
    d/dt(buccal_slab16) <- jf15 - jf16
    d/dt(buccal_slab17) <- jf16 - jf17
    d/dt(buccal_slab18) <- jf17 - jf18
    d/dt(buccal_slab19) <- jf18 - jf19
    d/dt(buccal_slab20) <- jf19 - jout

    rblood <- jout

    # =====================================================================
    # 6. NICOTINE TISSUE CONCENTRATIONS (mg/L)
    #    The lung and heart venous-outflow concentrations invert a Langmuir
    #    binding isotherm; the remaining tissues are flow-limited. Fat uses
    #    a MULTIPLICATION by its partition coefficient where every other
    #    tissue divides -- reproduced as coded; see vignette Errata.
    # =====================================================================
    ca    <- a_arterial / vab
    ch    <- a_heart / vh
    cbr   <- a_brain / vb
    cr    <- a_rapidly_perfused / vr
    cliv  <- a_liver / vl
    csk   <- a_skin / vsk
    cmu   <- a_muscle / vm
    cslow <- a_slowly_perfused / vs
    cfa   <- a_fat / vf

    cvh   <- 0.5 * (sqrt(((bmh + vh * PH * KBH - a_heart) / (vh * PH))^2 +
                           4 * a_heart * KBH / (vh * PH)) -
                      (bmh + vh * PH * KBH - a_heart) / (vh * PH))
    cvpul <- 0.5 * (sqrt(((bmlu + vpul * PLU * KBLU - a_pulmonary) / (vpul * PLU))^2 +
                           4 * a_pulmonary * KBLU / (vpul * PLU)) -
                      (bmlu + vpul * PLU * KBLU - a_pulmonary) / (vpul * PLU))
    cvbr  <- cbr / PBR
    cvr   <- cr / PR
    cvl   <- cliv / PL
    cvsk  <- csk / PSK
    cvm   <- cmu / PM
    cvfa  <- cfa * PF
    cvs   <- cslow / PS
    cvbus <- (a_buccal / vbu) / PLU
    cvcas <- (a_conducting_airway / vca) / PLU
    cvtas <- (a_transitional_airway / vta) / PLU

    # =====================================================================
    # 7. NICOTINE MASS BALANCE
    # =====================================================================
    rametl <- clm * cvl

    d/dt(a_venous) <- qm * cvm + qh * cvh + qf * cvfa + qb * cvbr + qs * cvs +
      qsk * cvsk + qr * cvr + ql * cvl + qta * cvtas + qbu * cvbus +
      qca * cvcas - qc * cvvb
    d/dt(a_arterial)            <- qc * (cvpul - ca) - clk * ca
    d/dt(a_urine)               <- clk * ca
    d/dt(a_pulmonary)           <- qc * (cvvb - cvpul)
    d/dt(a_buccal)              <- qbu * (ca - cvbus) + rblood
    d/dt(a_conducting_airway)   <- qca * (ca - cvcas)
    d/dt(a_transitional_airway) <- qta * (ca - cvtas)
    d/dt(a_gut)                 <- FA * gitr - KA1 * a_gut
    d/dt(gut_lumen)             <- (1 - FA) * gitr
    d/dt(a_liver)               <- ql * (ca - cvl) - rametl + KA1 * a_gut
    d/dt(a_hepatic)             <- rametl
    d/dt(a_heart)               <- qh * (ca - cvh)
    d/dt(a_brain)               <- qb * (ca - cvbr)
    d/dt(a_rapidly_perfused)    <- qr * (ca - cvr)
    d/dt(a_skin)                <- qsk * (ca - cvsk)
    d/dt(a_muscle)              <- qm * (ca - cvm)
    d/dt(a_fat)                 <- qf * (ca - cvfa)
    d/dt(a_slowly_perfused)     <- qs * (ca - cvs)

    # =====================================================================
    # 8. COTININE SUB-MODEL
    #    Rostami 2022 lumps cotinine into five tissue groups: rapidly
    #    perfused + brain, liver, muscle + heart, fat, and slowly perfused
    #    + skin + buccal + airways. Cotinine is formed in the liver from
    #    the metabolized-nicotine flux, mass-corrected by the cotinine /
    #    nicotine molecular-weight ratio.
    # =====================================================================
    cmr  <- a_rapidly_perfused_cot / (vr + vb)
    cml  <- a_liver_cot / vl
    cmm  <- a_muscle_cot / (vm + vh)
    cmf  <- a_fat_cot / vf
    cms  <- a_slowly_perfused_cot / (vs + vsk + vbu + vca + vta)
    cvmr <- cmr / PMR
    cvml <- cml / PML
    cvmm <- cmm / PMM
    cvmf <- cmf / PMF
    cvms <- cms / PMS
    cvbm <- a_venous_cot / vvb
    cam  <- a_arterial_cot / vab

    d/dt(a_venous_cot) <- (qr + qb) * cvmr + ql * cvml + (qm + qh) * cvmm +
      qf * cvmf + (qs + qsk + qbu + qca + qta) * cvms - qc * cvbm
    d/dt(a_arterial_cot)         <- qc * (cvbm - cam) - clkm * cam
    d/dt(a_urine_cot)            <- clkm * cam
    d/dt(a_rapidly_perfused_cot) <- (qr + qb) * (cam - cvmr)
    d/dt(a_liver_cot)            <- ql * (cam - cvml) + FNC * rametl * (MWCOT / MW) - cllm * cvml
    d/dt(a_hepatic_cot)          <- cllm * cvml
    d/dt(a_muscle_cot)           <- (qm + qh) * (cam - cvmm)
    d/dt(a_fat_cot)              <- qf * (cam - cvmf)
    d/dt(a_slowly_perfused_cot)  <- (qs + qsk + qbu + qca + qta) * (cam - cvms)

    # =====================================================================
    # 9. NICOTINE GLUCURONIDE SUB-MODEL (one compartment, renal clearance)
    # =====================================================================
    cgb <- central_gluc / vgd
    d/dt(central_gluc) <- FNG * rametl * (MWGLU / MW) - gcl * cgb
    d/dt(a_urine_gluc) <- gcl * cgb

    # =====================================================================
    # 10. OUTPUTS -- mg/L converted to ng/mL. Salehi 2025 reports venous
    #     plasma nicotine; the arterial output reproduces the arterial /
    #     venous separation of Figure S3. No residual-error model is
    #     reported in any on-disk source, so none is declared.
    # =====================================================================
    Cc      <- cvvb * 1000
    Cart    <- ca * 1000
    Cc_cot  <- cvbm * 1000
    Cc_gluc <- cgb * 1000
    HR      <- eout
  })
}
