Sun_2026_tilmicosin_pbpk <- function() {
  description <- paste(
    "Veterinary (pig). PBPK (whole-body, seventeen-state, Berkeley Madonna",
    "10.6.1) for the macrolide tilmicosin after oral administration in swine,",
    "developed to optimise the dosage regimen against Pasteurella multocida",
    "and to project the withdrawal interval for edible tissues (Sun et al.",
    "2026, J Agric Food Chem 74:4754-4766). Oral absorption is a lumped",
    "four-segment luminal transit chain: stomach -> duodenum -> remaining",
    "small intestine -> large intestine -> feces, with absorption into the",
    "portal circulation from the duodenum (Ka1) and, far more slowly, from",
    "the remaining small intestine (Ka2). Absorbed drug enters the liver, so",
    "the model is first-pass. Liver, kidney, muscle and a lumped",
    "rest-of-body compartment are perfusion-limited and well stirred; the",
    "lung is permeability-limited and resolves into pulmonary blood,",
    "pulmonary interstitial fluid (PIF) and pulmonary tissue, because PIF is",
    "the site of action for respiratory pathogens and was measured directly",
    "by microdialysis in conscious pigs in this study. Elimination is",
    "hepatobiliary (the major route) plus renal (minor); biliary output is",
    "returned to the duodenum, so the model carries enterohepatic",
    "recirculation. The eight parameters that Sun 2026 Table S4 carried",
    "through the 1000-animal Monte Carlo pop-PBPK analysis carry",
    "between-animal variability here; every other parameter is a fixed",
    "physiological or fitted chemical-specific constant. The model was",
    "hand-calibrated in Berkeley Madonna against microdialysis and digitised",
    "literature data, so no standard errors and no residual-error model are",
    "reported -- the propSd term is a placeholder. Three main-text versus",
    "supplement-code contradictions were resolved numerically against the",
    "paper's own printed claims rather than reproduced as coded; see the",
    "vignette Errata."
  )
  reference <- paste(
    "Sun L, Zhang C, Mi K, Wang H, Pan Y, Tao Y, Huang L.",
    "Dose optimization of tilmicosin against Pasteurella multocida in swine",
    "by physiologically based pharmacokinetic-pharmacodynamic model.",
    "J Agric Food Chem. 2026;74(8):4754-4766.",
    "doi:10.1021/acs.jafc.5c11368.",
    "Model equations transcribed from the Supporting Information (Berkeley",
    "Madonna code) and main-text equations 3-5; parameter values from",
    "Table 2 and Supporting Information Tables S3 and S4.",
    sep = " "
  )
  vignette <- "Sun_2026_tilmicosin"

  # Time in hours, doses in mg (the paper prescribes mg/kg; the vignette
  # multiplies by WT), concentrations in ug/mL. States hold amounts in mg
  # and volumes are in L, so amount/volume is mg/L, which is numerically
  # identical to the ug/mL the paper reports (and to the mg/kg used for
  # the maximum residue limits in liver, kidney and muscle).
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against the Supporting Information code
  # comments, each of which names its state's matrix and gives mg.
  compartmentData <- list(
    stomach            = list(analyte = "tilmicosin", units = "mg", specimen = "administration site", verified = TRUE),
    duodenum           = list(analyte = "tilmicosin", units = "mg", specimen = "administration site", verified = TRUE),
    a_small_intestine  = list(analyte = "tilmicosin", units = "mg", specimen = "administration site", verified = TRUE),
    a_large_intestine  = list(analyte = "tilmicosin", units = "mg", specimen = "administration site", verified = TRUE),
    a_oral_absorbed    = list(analyte = "tilmicosin", units = "mg", specimen = "not applicable", verified = TRUE),
    liver              = list(analyte = "tilmicosin", units = "mg", specimen = "tissue", verified = TRUE),
    a_bile             = list(analyte = "tilmicosin", units = "mg", specimen = "bile", verified = TRUE),
    venous             = list(analyte = "tilmicosin", units = "mg", specimen = "plasma", verified = TRUE),
    arterial           = list(analyte = "tilmicosin", units = "mg", specimen = "plasma", verified = TRUE),
    kidney             = list(analyte = "tilmicosin", units = "mg", specimen = "tissue", verified = TRUE),
    urine              = list(analyte = "tilmicosin", units = "mg", specimen = "urine", verified = TRUE),
    muscle             = list(analyte = "tilmicosin", units = "mg", specimen = "tissue", verified = TRUE),
    other              = list(analyte = "tilmicosin", units = "mg", specimen = "tissue", verified = TRUE),
    vp_lung            = list(analyte = "tilmicosin", units = "mg", specimen = "tissue", verified = TRUE),
    is_lung            = list(analyte = "tilmicosin", units = "mg", specimen = "tissue", verified = TRUE),
    int_lung           = list(analyte = "tilmicosin", units = "mg", specimen = "tissue", verified = TRUE),
    a_feces            = list(analyte = "tilmicosin", units = "mg", specimen = "faeces", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Sun 2026 Table S3 and the Supporting Information code both fix",
        "BW = 40 kg. Every organ volume and blood flow is a fixed fraction",
        "of BW and both elimination clearances are per-kg constants, so WT",
        "is carried here as a covariate exactly as BW enters the published",
        "code. Note that the model's 40 kg reference does NOT match the",
        "20 +/- 1 kg pigs of the authors' own microdialysis experiment; the",
        "calibration and validation datasets span 10-40 mg/kg doses in pigs",
        "the paper does not weigh individually. Because volumes, flows,",
        "clearances and the mg/kg dose all scale linearly in WT, predicted",
        "CONCENTRATIONS are invariant to WT in this model and only the",
        "amounts change; the lung transfer terms are the sole exception,",
        "since K2 and K3 divide by weight-scaled sub-compartment volumes",
        "while K1 and K4 multiply an amount directly. See the vignette",
        "Errata."
      ),
      source_name        = "BW"
    )
  )

  population <- list(
    species        = "pig (crossbred Landrace x Large White x Duroc)",
    n_subjects     = 3L,
    n_studies      = 5L,
    age_range      = "nursery / grower swine",
    weight_range   = "20 +/- 1 kg (microdialysis experiment); model reference BW = 40 kg",
    disease_state  = "healthy",
    dose_range     = "40 mg/kg oral tilmicosin once daily for 3 consecutive days (the optimised regimen); 10, 15, 20, 40, 50 and 60 mg/kg regimens were simulated",
    regions        = "China",
    organism       = "Pasteurella multocida isolate ZJWZ-A, isolated in Wuhan, Hubei Province, China in 2023 and deposited in the National Reference Laboratory for Veterinary Drug Residues at Huazhong Agricultural University; tilmicosin MIC 8 ug/mL by CLSI broth microdilution",
    notes          = paste(
      "The three pigs contributing pulmonary interstitial fluid (Sun 2026",
      "Table S7, columns Lung-1 to Lung-3) are the animals of the authors'",
      "own microdialysis experiment, the first reported collection of PIF",
      "from conscious pigs. The structural model was additionally calibrated",
      "and validated against four previously published swine datasets (Sun",
      "2026 Table 1): 20 animals dosed 20 mg/kg twice daily for 15 days",
      "(liver, kidney, muscle) and 5 animals given a single 40 mg/kg dose",
      "(plasma) for calibration; 40 animals dosed 10 mg/kg three times daily",
      "for 15 days (liver, kidney), 20 animals dosed 15 mg/kg once daily for",
      "3 days (muscle) and 5 animals given a single 20 mg/kg dose (plasma)",
      "for validation. All data were from healthy animals. The pop-PBPK",
      "layer simulates 1000 virtual animals by Monte Carlo sampling of the",
      "eight influential parameters in Table S4; no individual-level fitting",
      "was performed and no between-animal variance was estimated from data",
      "-- the coefficients of variation are literature defaults (20% for",
      "partition coefficients and transport rate constants, 30% otherwise)",
      "as described in Sun 2026 'Monte Carlo analysis'."
    )
  )

  ini({
    # ---------------------------------------------------------------
    # Lumped oral absorption chain (Sun 2026 Table 2, "Absorption"; all
    # "estimated"). Sun 2026 'PBPK Model Parametrization': the oral
    # absorption parameters were estimated by visually fitting the model
    # to the 20 mg/kg data set of Shen et al. (ref 26) in Berkeley
    # Madonna 10.6.1.
    # ---------------------------------------------------------------
    lkst <- log(0.8)
    label("Gastric emptying rate constant Kst, stomach to duodenum (1/h)")                  # Sun 2026 Table 2 (Kst = 0.8 /h, estimated); SI code `Kst = 0.8`

    lkt_duodenum <- log(2.2)
    label("Duodenal transit rate constant Kd, duodenum to remaining small intestine (1/h)") # Sun 2026 Table 2 (Kd = 2.2 /h, estimated); SI code `Kd = 2.2`

    lkt_ileocolic <- log(2.0)
    label("Ileocolic transit rate constant Kint, small to large intestine (1/h)")           # Sun 2026 Table 2 (Kint = 2.0 /h, estimated); SI code `Kint = 2`

    lka_duodenum <- log(0.8)
    label("Duodenal absorption rate constant Ka1 (1/h)")                                    # Sun 2026 Table 2 (Ka1 = 0.8 /h, estimated); SI code `Ka1 = 0.8`

    lka_small_intestine <- log(0.007)
    label("Remaining small intestinal absorption rate constant Ka2 (1/h)")                  # Sun 2026 Table 2 (Ka2 = 0.007 /h, estimated); SI code `Ka2 = 0.007`

    lkf <- log(0.05)
    label("Fecal excretion rate constant Kf, large intestine to feces (1/h)")               # Sun 2026 Table 2 (Kf = 0.05 /h, estimated, described in the Figure 1 legend as the "fecal excretion rate"); SI code `Kf = 0.05`

    # ---------------------------------------------------------------
    # Tissue-to-plasma partition coefficients (Sun 2026 Table 2, all
    # "estimated"). Sun 2026 'PBPK Model Parametrization': "In the
    # absence of pig-specific data for TIL, tissue/plasma partition
    # coefficients were initially assigned values reported for
    # tulathromycin in goats and were subsequently optimized by visually
    # aligning model simulations with measured TIL plasma
    # concentrations."
    # ---------------------------------------------------------------
    lkp_liver <- log(20)
    label("Liver-to-plasma partition coefficient PL (unitless)")                            # Sun 2026 Table 2 (PL = 20, estimated); SI code `PL = 20`

    lkp_kidney <- log(23)
    label("Kidney-to-plasma partition coefficient PK (unitless)")                           # Sun 2026 Table 2 (PK = 23, estimated); SI code `PK = 23`

    lkp_muscle <- log(1.8)
    label("Muscle-to-plasma partition coefficient PM (unitless)")                           # Sun 2026 Table 2 (PM = 1.8, estimated); SI code `PM = 1.8`

    lkp_lung <- log(2.2)
    label("Lung-to-plasma partition coefficient PLu (unitless)")                            # Sun 2026 Table 2 (PLu = 2.2, estimated); SI code `PLU = 2.2`

    lkp_other <- log(6.24)
    label("Rest-of-body-to-plasma partition coefficient PR (unitless)")                     # Sun 2026 Table 2 (PR = 6.24, estimated); SI code `PR = 6.24`

    # ---------------------------------------------------------------
    # Permeability-limited lung sub-model (Sun 2026 Table 2, "Transport
    # Parameters in Lung Subcompartment"; all "estimated").
    #
    # RESOLVED CONTRADICTION 1 (vignette Errata E1). Sun 2026 Table 2
    # and the SI code disagree about which numeric value belongs to
    # which direction; the two readings are transposed WITHIN each
    # bidirectional pair. The SI code's own `K1..K4` declarations are
    # dead code -- its ODEs reference KBE / KEB / KET / KTE, which are
    # never declared anywhere -- so the only surviving name-to-slot map
    # is main-text equations 3-5, which use K2 for pulmonary blood ->
    # PIF, K1 for PIF -> pulmonary blood, K4 for PIF -> tissue and K3
    # for tissue -> PIF. Feeding the SI's declared VALUES through that
    # main-text SLOT map is the reading shipped below, and it is the
    # only one of the sixteen possible readings that reproduces the
    # paper's own Figure 3B claim ("all predicted values except one at
    # the 48 h time point fell within the 2-fold error range"). See the
    # vignette for the full sweep.
    # ---------------------------------------------------------------
    lk_blood_isf <- log(2.6355773)
    label("Pulmonary blood to interstitial fluid transfer K2 (1/h)")                        # Sun 2026 eq 3 slot K2 (blood -> PIF) carrying the SI code's declared value `K2 = 2.6355773`; Table 2 prints 0.03 for this direction. Errata E1

    lk_isf_blood <- log(0.03)
    label("Pulmonary interstitial fluid to blood transfer K1 (1/h)")                        # Sun 2026 eq 3 slot K1 (PIF -> blood) carrying the SI code's declared value `K1 = 0.03`; Table 2 prints 2.63 for this direction. Errata E1

    lk_isf_tissue <- log(0.3)
    label("Pulmonary interstitial fluid to lung tissue transfer K4 (1/h)")                  # Sun 2026 eq 4 slot K4 (PIF -> tissue) carrying the SI code's declared value `K4 = 0.3`; Table 2 prints 9.65 for this direction. Errata E1

    lk_tissue_isf <- log(9.645148)
    label("Lung tissue to pulmonary interstitial fluid transfer K3 (1/h)")                  # Sun 2026 eq 4 slot K3 (tissue -> PIF) carrying the SI code's declared value `K3 = 9.645148`; Table 2 prints 0.30 for this direction. Errata E1

    lfb_lung <- log(0.200)
    label("Fraction of tilmicosin bound to pulmonary tissue protein Pt (unitless)")         # Sun 2026 Table 2 (Pt = 0.200, "pulmonary tissue protein binding"); SI code `Pt = 0.2`. The SI lung equations call this quantity `PBtissue`, which is never assigned in the code; see vignette Errata E5

    # ---------------------------------------------------------------
    # Plasma protein binding (Sun 2026 Table 2, "pulmonary plasma
    # protein binding"). Table 2 reports Pb = 0.188 as the BOUND
    # fraction; the SI code uses it as `(1 - PB)`, so the unbound
    # fraction is 1 - 0.188. Encoded on the canonical unbound scale.
    # Applied only in the lung sub-model, exactly as coded: unlike the
    # sibling Mi 2023 cefquinome model, Sun 2026 perfuses the systemic
    # tissues with TOTAL rather than unbound arterial drug.
    # ---------------------------------------------------------------
    fu <- fixed(0.812)
    label("Fraction of tilmicosin unbound in pulmonary plasma (unitless) = 1 - Pb")         # Sun 2026 Table 2 (Pb = 0.188 bound); SI code `Pb = 0.188`, `Free = 1-PB`, used as `(1 - PB)` in the RLUB / RLUE equations

    # ---------------------------------------------------------------
    # Elimination (Sun 2026 Table 2, "Elimination"; both "estimated").
    #
    # RESOLVED CONTRADICTION 3 (vignette Errata E3). The SI code
    # declares `KbileC = 0.002` and `KurineC = 0.126`, i.e. swapped
    # relative to Table 2, and swaps their unit comments alongside. Three
    # independent documents agree against the code: Table 2 (hepatic
    # clearance 0.126, renal clearance 0.002), Table S4 (Monte Carlo
    # means 0.126 for KbileC and 0.002 for KurineC) and the main-text
    # prose "Following excretion in bile, the majority are excreted in
    # feces, while a minor portion is eliminated via renal excretion in
    # urine." The Table 2 assignment is shipped.
    # ---------------------------------------------------------------
    lcl_nonren <- log(0.126)
    label("Hepatobiliary (non-renal) clearance KbileC (L/h/kg), applied to liver venous concentration") # Sun 2026 Table 2 (KbileC = 0.126 L/h/kg) + Table S4 (mean 0.126); the SI code declares 0.002. Errata E3

    lcl_renal <- log(0.002)
    label("Renal clearance KurineC (L/h/kg), applied to kidney venous concentration")       # Sun 2026 Table 2 (KurineC = 0.002 L/h/kg) + Table S4 (mean 0.002); the SI code declares 0.126. Errata E3

    # ---------------------------------------------------------------
    # Cardiac output. Held fixed (literature swine physiology, Sun 2026
    # Table S3) but living in ini() for parity with the other scaling
    # constants; every other physiological constant is written directly
    # in model(). Sun 2026 did NOT carry any physiological parameter
    # through the Monte Carlo analysis -- Table S4 lists only
    # chemical-specific parameters -- so unlike Mi 2023 this one has no
    # between-animal variability.
    # ---------------------------------------------------------------
    qcc <- fixed(4.944)
    label("Cardiac output QCC (L/h/kg)")                                                    # Sun 2026 Table S3 (QCC = 4.944 L/h/kg); SI code `QCC = 4.944`

    # ---------------------------------------------------------------
    # Between-animal variability, from the Monte Carlo distributions of
    # Sun 2026 Table S4. These are NOT estimated variances -- they are
    # the literature default coefficients of variation the authors
    # assumed (Sun 2026 'Monte Carlo analysis': "The default coefficient
    # of variation (CV) for partition coefficients (PCs) and transport
    # constant rates is assumed to be 20%. The default CVs of other
    # parameters were set as 30%."), so all are fixed().
    #
    # All eight Table S4 parameters are lognormal, so each takes a
    # log-scale eta with variance log(1 + CV^2). Table S4's printed
    # Upper / Lower bounds are the 2.5 / 97.5 percentiles; see the
    # vignette source-trace table for the arithmetic.
    #
    # Table S4 parameterises each lognormal so that its ARITHMETIC MEAN
    # equals the Table 2 point estimate, which puts the median about
    # 1-2% below it. Here the Table 2 value is used as the typical value
    # (i.e. the median), so that a typical-value simulation reproduces
    # the deterministic Table 2 model and Figures 3-4 exactly. See
    # vignette Errata E6.
    # ---------------------------------------------------------------
    etalkp_liver ~ fixed(0.03922071)
    # Sun 2026 Table S4 (PL, lognormal, CV 0.2): log(1 + 0.20^2) = 0.03922071
    etalkp_kidney ~ fixed(0.03922071)
    # Sun 2026 Table S4 (PK, lognormal, CV 0.2)
    etalkp_muscle ~ fixed(0.03922071)
    # Sun 2026 Table S4 (PM, lognormal, CV 0.2)
    etalkp_lung ~ fixed(0.03922071)
    # Sun 2026 Table S4 (PLU, lognormal, CV 0.2)

    etalkt_duodenum ~ fixed(0.0861777)
    # Sun 2026 Table S4 (Kd, lognormal, CV 0.3): log(1 + 0.30^2) = 0.0861777
    etalka_duodenum ~ fixed(0.0861777)
    # Sun 2026 Table S4 (Ka1, lognormal, CV 0.3)
    etalcl_renal ~ fixed(0.0861777)
    # Sun 2026 Table S4 (KurineC, lognormal, CV 0.3)
    etalcl_nonren ~ fixed(0.0861777)
    # Sun 2026 Table S4 (KbileC, lognormal, CV 0.3)

    # ---------------------------------------------------------------
    # Sun 2026 calibrated the model by hand in Berkeley Madonna 10.6.1
    # by minimising the sum of squared residuals and reports no
    # residual-error model, no standard errors and no objective
    # function. Model performance is reported only as MAPE (32.7-54.4%
    # in PIF, muscle, liver and kidney; 90.8% in plasma; Table S8
    # narrative). nlmixr2 model definitions require a residual-error
    # term, so propSd below is a fixed placeholder for syntactic
    # completeness only and must NOT be read as an estimate. Same
    # convention as Mi_2023_cefquinome_pbpk.
    # ---------------------------------------------------------------
    propSd <- fixed(0.10)
    label("Proportional residual error placeholder (fraction)")                             # not reported in Sun 2026; placeholder only
  })

  model({
    # =================================================================
    # Swine physiology (Sun 2026 Table S3, refs 31-32; SI code
    # {Physiological Parameters}). Blood flows are fractions of cardiac
    # output; volumes are fractions of body weight.
    # =================================================================
    q_co <- qcc * WT                                  # L/h; SI code `QC = QCC * BW`

    # Fractional flows. The rest-of-body fraction is computed as the
    # complement exactly as the code does. It evaluates to 0.3565,
    # matching the 0.3565 printed in Table S3.
    qlc_frac <- 0.2725                                # Sun 2026 Table S3 (QLC = 0.2725)
    qkc_frac <- 0.12                                  # Sun 2026 Table S3 (QKC = 0.1200)
    qmc_frac <- 0.251                                 # Sun 2026 Table S3 (QMC = 0.2510)
    qrc_frac <- 1 - qlc_frac - qkc_frac - qmc_frac    # SI code `QRC = 1 - QLC - QKC - QMC`

    q_liver <- qlc_frac * q_co
    q_kidney <- qkc_frac * q_co
    q_muscle <- qmc_frac * q_co
    q_other <- qrc_frac * q_co
    # SI code `QLUC = 1`: the whole cardiac output perfuses the lung.

    # Organ volumes (L). SI code
    # `VRC = 1 - VLC - VKC - VMC - VartC - VLUC - VvenC`, which
    # evaluates to 0.5013, matching the 0.5013 printed in Table S3.
    v_liver <- 0.0247 * WT                            # Sun 2026 Table S3 (VLC = 0.0247)
    v_kidney <- 0.004 * WT                            # Sun 2026 Table S3 (VKC = 0.0040)
    v_muscle <- 0.4 * WT                              # Sun 2026 Table S3 (VMC = 0.4000)
    v_arterial <- 0.044 * WT                          # Sun 2026 Table S3 (VartC = 0.0440)
    v_venous <- 0.016 * WT                            # Sun 2026 Table S3 (VvenC = 0.0160)
    v_lung <- 0.01 * WT                               # Sun 2026 Table S3 (VLUC = 0.0100)
    v_other <- (1 - 0.0247 - 0.004 - 0.4 - 0.044 - 0.01 - 0.016) * WT   # SI code `VRC`

    # Lung sub-compartment volumes (L). SI code
    # `VLUB = VLU * FVBLU`, `VLUE = VLU * FVILU`, `VLUT = VLU - VLUE - VLUB`.
    #
    # RESOLVED CONTRADICTION 2 (vignette Errata E2). Table S3 gives
    # FVBLU = 0.0200 and FVILU = 0.4850; the SI code declares 0.262 and
    # 0.188. The code's two values are exactly the cefquinome numbers of
    # the same laboratory's Mi 2023 model (shipped here as
    # Mi_2023_cefquinome_pbpk, `v_vp_lung <- v_lung * 0.262` and
    # `v_is_lung <- v_lung * 0.188`), i.e. a copy-paste leftover from
    # the shared Berkeley Madonna template. Table S3's values are
    # shipped, and they are required for the Figure 3B match.
    v_vp_lung <- v_lung * 0.02                        # Sun 2026 Table S3 (FVBLU = 0.0200); the SI code declares 0.262. Errata E2
    v_is_lung <- v_lung * 0.485                       # Sun 2026 Table S3 (FVILU = 0.4850); the SI code declares 0.188. Errata E2
    v_int_lung <- v_lung - v_is_lung - v_vp_lung      # SI code `VLUT = VLU - VLUE - VLUB`; note this evaluates to 0.495 of lung, not the 0.0100 printed for FVTLU in Table S3. Errata E2

    # =================================================================
    # Individual chemical-specific parameters
    # =================================================================
    kst <- exp(lkst)
    kt_duodenum <- exp(lkt_duodenum + etalkt_duodenum)
    kt_ileocolic <- exp(lkt_ileocolic)
    ka_duodenum <- exp(lka_duodenum + etalka_duodenum)
    ka_small_intestine <- exp(lka_small_intestine)
    kf <- exp(lkf)

    kp_liver <- exp(lkp_liver + etalkp_liver)
    kp_kidney <- exp(lkp_kidney + etalkp_kidney)
    kp_muscle <- exp(lkp_muscle + etalkp_muscle)
    kp_lung <- exp(lkp_lung + etalkp_lung)
    kp_other <- exp(lkp_other)

    k_blood_isf <- exp(lk_blood_isf)
    k_isf_blood <- exp(lk_isf_blood)
    k_isf_tissue <- exp(lk_isf_tissue)
    k_tissue_isf <- exp(lk_tissue_isf)
    fb_lung <- exp(lfb_lung)

    # Per-kg clearances scaled to the individual (L/h). SI code
    # `Kbile = KbileC * BW` and `Kurine = KurineC * BW`.
    cl_nonren <- exp(lcl_nonren + etalcl_nonren) * WT
    cl_renal <- exp(lcl_renal + etalcl_renal) * WT

    # =================================================================
    # Concentrations (ug/mL == mg/L). States hold amounts in mg.
    # =================================================================
    c_venous <- venous / v_venous                     # SI code `CV = AV/Vven`
    c_arterial <- arterial / v_arterial               # SI code `CA = AA/Vart`

    # Emergent venous concentrations leaving each perfusion-limited
    # organ, SI code `CVL = AL/(VL*PL)`. Note that `CVK`, `CVM` and
    # `CVR` are USED by the SI code but never assigned in it; they take
    # the standard flow-limited form A/(V*P) shown for `CVL`. See
    # vignette Errata E5.
    cv_liver <- liver / (v_liver * kp_liver)
    cv_kidney <- kidney / (v_kidney * kp_kidney)
    cv_muscle <- muscle / (v_muscle * kp_muscle)
    cv_other <- other / (v_other * kp_other)

    # Lung. The three sub-compartments equilibrate with the systemic
    # circulation through a single lumped lung:plasma partition
    # coefficient applied to the TOTAL lung amount, SI code
    # `CVLU = (ALUE+ALUT+ALUB)/(VLU*PLU)` and `CLU = (ALUE+ALUT+ALUB)/VLU`.
    a_lung_total <- vp_lung + is_lung + int_lung
    cv_lung <- a_lung_total / (v_lung * kp_lung)
    c_lung <- a_lung_total / v_lung
    c_is_lung <- is_lung / v_is_lung                  # SI code `CLUE = ALUE/VLUE`

    # =================================================================
    # Fluxes (mg/h)
    # =================================================================
    # Absorption out of the two absorbing luminal segments, SI code
    # `RAO = Ka1 * ADI + Ka2 * AOSI`.
    r_absorb <- ka_duodenum * duodenum + ka_small_intestine * a_small_intestine
    # Hepatobiliary output, SI code `Rmet = Kbile * CVL`. Despite the
    # code's `met` naming, Sun 2026 describes this route as biliary
    # excretion, and the flux is returned to the duodenum below rather
    # than removed from the system.
    r_bile <- cl_nonren * cv_liver
    r_renal <- cl_renal * cv_kidney                   # SI code `Rurine = Kurine * CVK`
    r_feces <- kf * a_large_intestine                 # SI code `Rf = Kf * ALI`

    # Lung permeability fluxes, SI code RLUB / RLUE / RLUT and main-text
    # equations 3-5. NOTE the dimensional asymmetry, printed identically
    # in both documents and reproduced here as published: the inbound
    # terms multiply a CONCENTRATION (amount / volume) while the
    # outbound terms multiply an AMOUNT, so the four constants cannot
    # all carry the "/h" units Table 2 assigns them. See vignette
    # Errata E4. Only unbound drug crosses out of pulmonary blood and
    # out of pulmonary tissue, so those two terms carry fu and
    # (1 - Pt) respectively, while interstitial fluid is treated as
    # entirely unbound (it is what microdialysis samples).
    r_blood_isf <- k_blood_isf * (vp_lung * fu / v_vp_lung)
    # SI code `KBE*(ALUB*(1-PB)/VLUB)`; main text eq 3 `K2 x (ALUB x (1-Pb)/VLUB)`
    r_isf_blood <- k_isf_blood * is_lung
    # SI code `KEB*ALUE`; main text eq 3 `K1 x ALUE`
    r_isf_tissue <- k_isf_tissue * is_lung
    # SI code `KET*ALUE`; main text eq 4 `K4 x ALUE`
    r_tissue_isf <- k_tissue_isf * (int_lung * (1 - fb_lung) / v_int_lung)
    # SI code `KTE*(ALUT*(1-PBtissue)/VLUT)`; main text eq 4 `K3 x (ALUT x (1-Pt)/VLUT)`

    # =================================================================
    # ODE system (SI code)
    # =================================================================
    # Lumped luminal transit chain. The oral dose enters `stomach`.
    d/dt(stomach) <- -kst * stomach
    # SI code `RAST = RDOSEoral - Kst * AST`

    # The duodenum receives gastric emptying AND the biliary output, so
    # the model carries enterohepatic recirculation.
    d/dt(duodenum) <-
      kst * stomach -
      kt_duodenum * duodenum -
      ka_duodenum * duodenum +
      r_bile
    # SI code `RDI = Kst * AST - Kd * ADI - Ka1 * ADI + Rmet`

    # SI code `ROSI = Kdm * ADI - Ka2 * AOSI - Kint * AOSI`. `Kdm` is
    # used here but never declared anywhere in the deposited code; it is
    # the duodenal transit constant `Kd`, which is the only rate
    # constant that can drain the duodenum into this segment and is what
    # main-text Figure 1 and Table 2 describe. See vignette Errata E5.
    d/dt(a_small_intestine) <-
      kt_duodenum * duodenum -
      ka_small_intestine * a_small_intestine -
      kt_ileocolic * a_small_intestine

    d/dt(a_large_intestine) <- kt_ileocolic * a_small_intestine - r_feces
    # SI code `RLI = Kint * AOSI - Rf`

    # Cumulative absorbed amount, the input side of the published mass
    # balance (SI code `Input = AAO`). Bookkeeping only; no feedback.
    d/dt(a_oral_absorbed) <- r_absorb
    # SI code `RAO = Ka1 * ADI + Ka2 * AOSI`

    # Liver. Absorbed drug enters here rather than directly into venous
    # blood, so oral tilmicosin undergoes first pass.
    d/dt(liver) <- q_liver * (c_arterial - cv_liver) - r_bile + r_absorb
    # SI code `RL = QL * (CA - CVL) - Rmet + RAO`

    # Cumulative biliary output. Carried so the published mass balance
    # can be checked; the same flux is recycled into the duodenum above,
    # so this state is a record of biliary output, not a closed sink.
    d/dt(a_bile) <- r_bile
    # SI code `d/dt(Amet) = Rmet`

    # Venous blood receives every perfusion-limited organ's venous
    # outflow. Unlike the sibling Mi 2023 model there is no absorption
    # term here, because absorbed drug enters the liver.
    d/dt(venous) <-
      q_liver * cv_liver +
      q_kidney * cv_kidney +
      q_muscle * cv_muscle +
      q_other * cv_other -
      q_co * c_venous
    # SI code `RV = (QL*CVL + QK*CVK + QM*CVM + QR*CVR) - QC*CV`

    # Arterial blood is fed by the lung and drained by the tissues.
    # Sun 2026 perfuses the tissues with TOTAL arterial drug (`CA`),
    # not the unbound fraction; reproduced as coded.
    d/dt(arterial) <- q_co * (cv_lung - c_arterial)
    # SI code `RA = QC * (CVLU - CA)`

    # Perfusion-limited, well-stirred organs.
    d/dt(kidney) <- q_kidney * (c_arterial - cv_kidney) - r_renal
    # SI code `RK = QK * (CA - CVK) - Rurine`
    d/dt(muscle) <- q_muscle * (c_arterial - cv_muscle)
    # SI code `RM = QM * (CA - CVM)`
    d/dt(other) <- q_other * (c_arterial - cv_other)
    # SI code `RR = QR * (CA - CVR)`

    # Permeability-limited lung: pulmonary blood -> interstitial fluid
    # -> tissue.
    d/dt(vp_lung) <- q_co * (c_venous - cv_lung) - r_blood_isf + r_isf_blood
    # SI code `RLUB`; main text eq 3
    d/dt(is_lung) <- r_blood_isf - r_isf_blood - r_isf_tissue + r_tissue_isf
    # SI code `RLUE`; main text eq 4
    d/dt(int_lung) <- r_isf_tissue - r_tissue_isf
    # SI code `RLUT`; main text eq 5

    # Terminal excretion sinks.
    d/dt(urine) <- r_renal                            # SI code `d/dt(Aurine) = Rurine`
    d/dt(a_feces) <- r_feces                          # SI code `d/dt(Af) = Rf`

    # =================================================================
    # Observations
    # =================================================================
    # Sun 2026 compares model output against plasma, pulmonary
    # interstitial fluid, liver, kidney and muscle. The SI code makes no
    # blood-to-plasma correction, so the venous blood concentration is
    # the paper's plasma prediction.
    Cc <- c_venous
    # Free drug in pulmonary interstitial fluid, the microdialysis
    # measurement and the driver of the PIF AUC/MIC target.
    Cpif <- c_is_lung
    Cliver <- liver / v_liver
    Ckidney <- kidney / v_kidney
    Cmuscle <- muscle / v_muscle
    Clung <- c_lung

    Cc ~ prop(propSd)
  })
}
