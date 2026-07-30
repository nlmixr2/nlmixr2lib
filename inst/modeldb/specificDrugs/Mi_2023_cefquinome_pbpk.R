Mi_2023_cefquinome_pbpk <- function() {
  description <- paste(
    "Veterinary (pig). PBPK (whole-body, six-compartment, Berkeley",
    "Madonna 10.1.3) for the cephalosporin cefquinome after intramuscular",
    "injection in swine, developed to optimise the dosage regimen against",
    "respiratory-tract pathogens and to predict the withdrawal interval",
    "for edible tissues (Mi et al. 2023, PLoS Comput Biol 19:e1011331).",
    "Plasma (venous and arterial blood pools), liver, kidney, muscle and a",
    "lumped rest-of-body compartment are perfusion-limited and well",
    "stirred; the lung is permeability-limited and resolves into vascular",
    "blood, interstitial fluid and tissue sub-compartments, because lung",
    "interstitial fluid is the site of action for respiratory pathogens",
    "and was measured directly by microdialysis in this study. Only",
    "unbound drug crosses between the lung sub-compartments and only",
    "unbound arterial drug perfuses the tissues. Intramuscular absorption",
    "is a two-step process: 90% of the injected dose is immediately",
    "available at the injection site and is absorbed into venous blood at",
    "first-order rate Kim, while 10% is released from a slow depot at",
    "first-order rate Kdiss into the same injection-site pool. Elimination",
    "is renal (from kidney) plus hepatobiliary (from liver). The 11",
    "parameters that Mi 2023 Table 4 carried through the 1000-animal Monte",
    "Carlo pop-PBPK analysis carry between-animal variability here; every",
    "other parameter is a fixed physiological or fitted chemical-specific",
    "constant. The model was hand-calibrated in Berkeley Madonna against",
    "digitised literature data, so no standard errors and no residual-error",
    "model are reported -- the propSd term is a placeholder. Several",
    "Table 3 / supplement-code discrepancies are reproduced as coded; see",
    "the vignette Errata."
  )
  reference <- paste(
    "Mi K, Sun L, Hou Y, Cai X, Zhou K, Ma W, Xu X, Pan Y, Liu Z, Huang L.",
    "A physiologically based pharmacokinetic model to optimize the dosage",
    "regimen and withdrawal time of cefquinome in pigs. PLoS Comput Biol.",
    "2023;19(8):e1011331. doi:10.1371/journal.pcbi.1011331.",
    "Model equations transcribed from S1 Text (Berkeley Madonna code);",
    "parameter values from Table 3 and Table 4.",
    sep = " "
  )
  vignette <- "Mi_2023_cefquinome_pbpk"

  # The `bile` state is the cumulative amount eliminated by the
  # hepatobiliary route (`Amet` in the S1 Text code, whose own comment
  # reads "the amount of drug metabolized in liver"). It is a terminal
  # bookkeeping sink with no feedback into the dynamics; it exists so the
  # published mass-balance identity (S1 Text, {Mass balance equations})
  # can be checked. There is no canonical compartment for a biliary
  # elimination sink in inst/references/compartment-names.md, so it is
  # declared here rather than silently extending the register.
  paper_specific_compartments <- c("bile")

  # Time in hours, doses in mg (the paper prescribes mg/kg; the vignette
  # multiplies by WT), concentrations in ug/mL. States hold amounts in mg
  # and volumes are in L, so amount/volume is mg/L, which is numerically
  # identical to the ug/mL the paper reports (and to the PPM used for the
  # maximum residue limits in liver and kidney).
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Mi 2023 Table 3 fixes BW = 25 kg (a nursery pig, the class most",
        "susceptible to respiratory-tract pathogens). Every organ volume",
        "and blood flow in the S1 Text code is a fixed fraction of BW, and",
        "the renal and hepatobiliary clearances are per-kg constants, so",
        "WT is carried here as a covariate exactly as BW enters the",
        "published code. The lung permeability terms KBI / KIB / KIT / KTI",
        "are NOT scaled by body weight but are divided by weight-scaled",
        "sub-compartment volumes, so the lung kinetics are not",
        "weight-invariant. The authors state explicitly that the dosage",
        "regimen and withdrawal intervals derived from this model do not",
        "transfer to market-age swine (about 90-100 kg) and must be",
        "re-assessed. Calibration and validation data span 15-45 kg."
      ),
      source_name        = "BW"
    )
  )

  population <- list(
    species        = "pig (crossbred Landrace x Large White x Duroc)",
    n_subjects     = 4L,
    n_studies      = 5L,
    age_range      = "nursery / grower swine",
    weight_range   = "20 +/- 2 kg (microdialysis experiment); 15-45 kg across the calibration and validation datasets; model reference BW = 25 kg",
    disease_state  = "healthy",
    dose_range     = "2 mg/kg intramuscular cefquinome sulfate (label dose); extra-label regimens of 3, 4 and 5 mg/kg once or twice daily were simulated",
    regions        = "China",
    notes          = paste(
      "The four crossbred pigs are the animals of the authors' own",
      "microdialysis experiment (Mi 2023 Materials and methods,",
      "'Animals'), which supplied the lung interstitial-fluid data. The",
      "structural model was additionally calibrated and validated against",
      "four previously published swine datasets digitised with",
      "WebPlotDigitizer (Mi 2023 Table 2): Li/Wu ref [19] plasma (n = 5,",
      "25 kg, single 2 mg/kg IM) and Zhang/Li ref [18] liver and kidney",
      "(n = 40, 45 kg, five doses over 24 h) for calibration; Mi/Li ref",
      "[17] plasma (n = 6, 15 kg) and Xu/Yang ref [21] liver and kidney",
      "(n = 40, 30 kg) for validation. All data were from healthy",
      "animals. The pop-PBPK layer simulates 1000 virtual animals by",
      "Monte Carlo sampling of the 11 influential parameters in Table 4;",
      "no individual-level fitting was performed and no between-animal",
      "variance was estimated from data -- the coefficients of variation",
      "are literature defaults (20% for partition coefficients, 30%",
      "otherwise) as described in Mi 2023 'Pop-PBPK model'."
    )
  )

  ini({
    # ---------------------------------------------------------------
    # Intramuscular absorption (Mi 2023 Table 3, "Absorption rate
    # constant"; all three "Model fitting"). Mi 2023 Discussion:
    # "Kdiss, Frac, and Kim were optimized as 0.05, 0.1 and 7 per h".
    # ---------------------------------------------------------------
    lka <- log(7)
    label("Absorption rate constant Kim, injection site to venous blood (1/h)")            # Mi 2023 Table 3 (Kim = 7 /h, Model fitting); S1 Text `Kim = 7`

    lkdiss <- log(0.05)
    label("Dissolution rate constant Kdiss, slow depot to injection site (1/h)")           # Mi 2023 Table 3 (Kdiss = 0.05 /h, Model fitting); S1 Text `Kdiss = 0.05`

    frac <- 0.1
    label("Fraction of the intramuscular dose entering the slow-release depot (unitless)") # Mi 2023 Table 3 (Frac = 0.1, Model fitting); S1 Text `Rppgim = Rinputim*Frac`

    # ---------------------------------------------------------------
    # Tissue-to-plasma partition coefficients (Mi 2023 Table 3).
    # PL, PK and PLU are credited to ref [46] and are held fixed; PM and
    # PR are credited to "Model fitting". Note that the Mi 2023
    # "Model parameterization and calibration" prose says all partition
    # coefficients were adjusted against swine data; Table 3's per-row
    # Reference column is the more specific record and is followed here.
    # See vignette Errata E1.
    # ---------------------------------------------------------------
    lkp_liver <- fixed(log(6))
    label("Liver-to-plasma partition coefficient PL (unitless)")                           # Mi 2023 Table 3 (PL = 6, ref [46]); S1 Text `PL = 6`

    lkp_kidney <- fixed(log(15.2))
    label("Kidney-to-plasma partition coefficient PK (unitless)")                          # Mi 2023 Table 3 (PK = 15.2, ref [46]); S1 Text `PK = 15.2`

    lkp_lung <- fixed(log(1.5))
    label("Lung-to-plasma partition coefficient PLU (unitless)")                           # Mi 2023 Table 3 (PLU = 1.5, ref [46]); S1 Text `PLU = 1.5`

    lkp_muscle <- log(0.1)
    label("Muscle-to-plasma partition coefficient PM (unitless)")                          # Mi 2023 Table 3 (PM = 0.1, Model fitting); S1 Text `PM = 0.1`

    lkp_other <- log(0.1)
    label("Rest-of-body-to-plasma partition coefficient PR (unitless)")                    # Mi 2023 Table 3 (PR = 0.1, Model fitting); S1 Text `PR = 0.1`

    # ---------------------------------------------------------------
    # Permeability-limited lung sub-model (Mi 2023 Table 3, all "Model
    # fitting"). Mi 2023 Discussion: "KBI, KIB, KIT, and KTI were
    # optimized as 0.110, 0.052, 3.56, and 2.60 per h". Table 3 labels
    # these "/h", but in the S1 Text equations each multiplies a
    # CONCENTRATION to give an amount rate, so they act dimensionally as
    # clearances (L/h). Reproduced exactly as coded; see vignette
    # Errata E2.
    # ---------------------------------------------------------------
    lk_blood_isf <- log(0.110)
    label("Lung vascular blood to interstitial fluid transfer KBI (L/h, printed as 1/h)")  # Mi 2023 Table 3 (KBI = 0.110); S1 Text `KBI = 0.110`

    lk_isf_blood <- log(0.052)
    label("Lung interstitial fluid to vascular blood transfer KIB (L/h, printed as 1/h)")  # Mi 2023 Table 3 (KIB = 0.052); S1 Text `KIB = 0.052`

    lk_isf_tissue <- log(3.56)
    label("Lung interstitial fluid to lung tissue transfer KIT (L/h, printed as 1/h)")     # Mi 2023 Table 3 (KIT = 3.56); S1 Text `KIT = 3.56`

    lk_tissue_isf <- log(2.60)
    label("Lung tissue to interstitial fluid transfer KTI (L/h, printed as 1/h)")          # Mi 2023 Table 3 (KTI = 2.60); S1 Text `KTI = 2.6`

    lfb_lung <- log(0.30)
    label("Fraction of cefquinome bound to lung tissue protein PT (unitless)")             # Mi 2023 Table 3 (PT = 0.3, Model fitting); Discussion "PT was defined as 0.3 for the final model". The S1 Text lung equations call this quantity `PBtissue`, which is never assigned in the code; see vignette Errata E3

    # ---------------------------------------------------------------
    # Plasma protein binding (Mi 2023 Table 3, ref [31]). Table 3
    # reports PB = 0.188 as the BOUND fraction; the S1 Text code uses it
    # as `CAfree = CA*(1-PB)`, so the unbound fraction is 1 - 0.188.
    # Encoded on the canonical unbound scale.
    # ---------------------------------------------------------------
    fu <- fixed(0.812)
    label("Fraction of cefquinome unbound in plasma (unitless) = 1 - PB")                  # Mi 2023 Table 3 (PB = 0.188 bound, ref [31]); S1 Text `PB = 0.188`, `CAfree = CA*(1-PB)`

    # ---------------------------------------------------------------
    # Elimination (Mi 2023 Table 3, both "Model fitting"). Mi 2023
    # Discussion: "The elimination parameters KurineC and KblieC were
    # optimized as 0.3 and 0.01 L/h/kg". Both are per-kg clearances and
    # are multiplied by WT in model().
    # ---------------------------------------------------------------
    lcl_renal <- log(0.3)
    label("Renal clearance KurineC (L/h/kg), applied to kidney venous concentration")      # Mi 2023 Table 3 (KurineC = 0.3 L/h/kg); S1 Text `KurineC = 0.3`

    lcl_nonren <- log(0.01)
    label("Hepatobiliary (non-renal) clearance KbileC (L/h/kg), applied to liver venous concentration") # Mi 2023 Table 3 (KbileC = 0.01 L/h/kg); S1 Text `KbileC = 0.01`

    # ---------------------------------------------------------------
    # Physiological parameters that Mi 2023 Table 4 carried through the
    # Monte Carlo analysis. They are literature swine physiology (ref
    # [29]) and so are fixed, but they need to live in ini() to carry
    # their between-animal variability. Every other physiological
    # constant is written directly in model().
    # ---------------------------------------------------------------
    qcc <- fixed(4.944)
    label("Cardiac output QCC (L/h/kg)")                                                   # Mi 2023 Table 3 (QCC = 4.944 L/h/kg, ref [29]); S1 Text `QCC = 4.944`

    qkc <- fixed(0.1398)
    label("Fraction of cardiac output to the kidney QKC (unitless)")                       # Mi 2023 Table 3 (QKC = 0.1398, ref [29]); S1 Text `QKC = 0.1398`

    # ---------------------------------------------------------------
    # Between-animal variability, from the Monte Carlo distributions of
    # Mi 2023 Table 4. These are NOT estimated variances -- they are the
    # literature default coefficients of variation the authors assumed
    # (Mi 2023 "Pop-PBPK model": "The default coefficient of variation
    # (CV) for partition coefficients (PCs) and transport constant rates
    # are assumed to be 20%. The default CVs of other parameters were set
    # as 30%."), so all are fixed().
    #
    # Normal parameters (QCC, QKC) take an additive eta with
    # variance (mean * CV)^2. Lognormal parameters take a log-scale eta
    # with variance log(1 + CV^2). Both reproduce the 2.5 / 97.5
    # percentile bounds printed in Table 4; see the vignette source-trace
    # table for the arithmetic.
    #
    # Table 4 parameterises each lognormal so that its ARITHMETIC MEAN
    # equals the Table 3 point estimate, which puts the median about
    # 1-2% below it. Here the Table 3 value is used as the typical value
    # (i.e. the median), so that a typical-value simulation reproduces
    # the deterministic Table 3 model and Figures 3-4 exactly. See
    # vignette Errata E4.
    # ---------------------------------------------------------------
    etaqcc ~ fixed(2.199882)
    # Mi 2023 Table 4 (QCC, Normal, CV 0.3): SD = 4.944 * 0.3 = 1.4832; 1.4832^2 = 2.199882
    etaqkc ~ fixed(0.001758964)
    # Mi 2023 Table 4 (QKC, Normal, CV 0.3): SD = 0.1398 * 0.3 = 0.04194; 0.04194^2 = 0.001758964

    etalkp_liver ~ fixed(0.03922071)
    # Mi 2023 Table 4 (PL, Lognormal, CV 0.2): log(1 + 0.20^2) = 0.03922071
    etalkp_kidney ~ fixed(0.03922071)
    # Mi 2023 Table 4 (PK, Lognormal, CV 0.2)
    etalkp_muscle ~ fixed(0.03922071)
    # Mi 2023 Table 4 (PM, Lognormal, CV 0.2)
    etalkp_other ~ fixed(0.03922071)
    # Mi 2023 Table 4 (PR, Lognormal, CV 0.2)
    etalkp_lung ~ fixed(0.03922071)
    # Mi 2023 Table 4 (PLU, Lognormal, CV 0.2)

    etalk_isf_tissue ~ fixed(0.0861777)
    # Mi 2023 Table 4 (KIT, Lognormal, CV 0.3): log(1 + 0.30^2) = 0.0861777
    etalk_tissue_isf ~ fixed(0.0861777)
    # Mi 2023 Table 4 (KTI, Lognormal, CV 0.3)
    etalcl_renal ~ fixed(0.0861777)
    # Mi 2023 Table 4 (KurineC, Lognormal, CV 0.3)
    etalfb_lung ~ fixed(0.0861777)
    # Mi 2023 Table 4 (PT, Lognormal, CV 0.3)

    # ---------------------------------------------------------------
    # Mi 2023 calibrated the model by hand in Berkeley Madonna 10.1.3
    # against digitised literature data and reports no residual-error
    # model, no standard errors and no objective function. nlmixr2 model
    # definitions require a residual-error term, so propSd below is a
    # fixed placeholder for syntactic completeness only and must NOT be
    # read as an estimate. Same convention as
    # Kang_2023_artesunate_hamster_pbpk and An_2012_mitoxantrone_*_pbpk.
    # ---------------------------------------------------------------
    propSd <- fixed(0.10)
    label("Proportional residual error placeholder (fraction)")                            # not reported in Mi 2023; placeholder only
  })

  model({
    # =================================================================
    # Swine physiology (Mi 2023 Table 3, ref [29]; S1 Text
    # {Physiological Parameters}). Blood flows are fractions of cardiac
    # output; volumes are fractions of body weight.
    # =================================================================
    q_co <- (qcc + etaqcc) * WT   # L/h; S1 Text `QC = QCC * BW`

    # Fractional flows. The rest-of-body fraction is computed as the
    # complement exactly as the code does, so it tracks the sampled QKC.
    # S1 Text `QRC = 1 - QLC - QKC - QMC`. Note this evaluates to 0.3025,
    # not the 0.3055 printed in Mi 2023 Table 3; see vignette Errata E5.
    qlc_frac <- 0.3053                                # Mi 2023 Table 3 (QLC = 0.3053)
    qmc_frac <- 0.2524                                # Mi 2023 Table 3 (QMC = 0.2524)
    qkc_frac <- qkc + etaqkc                          # Mi 2023 Table 3 / Table 4 (QKC = 0.1398)
    qrc_frac <- 1 - qlc_frac - qkc_frac - qmc_frac    # S1 Text `QRC = 1 - QLC - QKC - QMC`

    q_liver <- qlc_frac * q_co
    q_kidney <- qkc_frac * q_co
    q_muscle <- qmc_frac * q_co
    q_other <- qrc_frac * q_co
    # S1 Text `QLUC = 1`: the whole cardiac output perfuses the lung.

    # Organ volumes (L). S1 Text
    # `VRC = 1 - VLC - VKC - VMC - VbloodC - VLUC` with VbloodC = 0.06
    # (= VartC 0.016 + VvenC 0.044). This evaluates to 0.4966, not the
    # 0.232 printed in Mi 2023 Table 3; see vignette Errata E5.
    v_liver <- 0.0294 * WT                            # Mi 2023 Table 3 (VLC = 0.0294)
    v_kidney <- 0.004 * WT                            # Mi 2023 Table 3 (VKC = 0.004)
    v_muscle <- 0.4 * WT                              # Mi 2023 Table 3 (VMC = 0.4)
    v_venous <- 0.044 * WT                            # Mi 2023 Table 3 (VvenC = 0.044)
    v_arterial <- 0.016 * WT                          # Mi 2023 Table 3 (VartC = 0.016)
    v_lung <- 0.01 * WT                               # Mi 2023 Table 3 (VLUC = 0.01)
    v_other <- (1 - 0.0294 - 0.004 - 0.4 - 0.06 - 0.01) * WT   # S1 Text `VRC`

    # Lung sub-compartment volumes (L). S1 Text
    # `VLUB = VLU * FVBLU`, `VLUI = VLU * FVILU`, `VLUT = VLU - VLUI - VLUB`.
    v_vp_lung <- v_lung * 0.262                       # Mi 2023 Table 3 (VLUB = 0.262 of lung, ref [30])
    v_is_lung <- v_lung * 0.188                       # Mi 2023 Table 3 (VLUI = 0.188 of lung, ref [8])
    v_int_lung <- v_lung - v_is_lung - v_vp_lung      # Mi 2023 Table 3 (VLUT = 0.55 of lung, Calculated)

    # =================================================================
    # Individual chemical-specific parameters
    # =================================================================
    ka <- exp(lka)
    kdiss <- exp(lkdiss)

    kp_liver <- exp(lkp_liver + etalkp_liver)
    kp_kidney <- exp(lkp_kidney + etalkp_kidney)
    kp_muscle <- exp(lkp_muscle + etalkp_muscle)
    kp_other <- exp(lkp_other + etalkp_other)
    kp_lung <- exp(lkp_lung + etalkp_lung)

    k_blood_isf <- exp(lk_blood_isf)
    k_isf_blood <- exp(lk_isf_blood)
    k_isf_tissue <- exp(lk_isf_tissue + etalk_isf_tissue)
    k_tissue_isf <- exp(lk_tissue_isf + etalk_tissue_isf)
    fb_lung <- exp(lfb_lung + etalfb_lung)

    # Per-kg clearances scaled to the individual (L/h).
    cl_renal <- exp(lcl_renal + etalcl_renal) * WT
    cl_nonren <- exp(lcl_nonren) * WT

    # =================================================================
    # Concentrations (ug/mL == mg/L). States hold amounts in mg.
    # =================================================================
    c_venous <- venous / v_venous                     # S1 Text `CV = AV/Vven`
    c_arterial <- arterial / v_arterial               # S1 Text `CA = AA/Vart`
    c_arterial_free <- c_arterial * fu                # S1 Text `CAfree = CA*(1-PB)`

    # Emergent venous concentrations leaving each perfusion-limited
    # organ, S1 Text `CVL = AL/(VL*PL)` and siblings.
    cv_liver <- liver / (v_liver * kp_liver)
    cv_kidney <- kidney / (v_kidney * kp_kidney)
    cv_muscle <- muscle / (v_muscle * kp_muscle)
    cv_other <- other / (v_other * kp_other)

    # Lung. The three sub-compartments equilibrate with the systemic
    # circulation through a single lumped lung:plasma partition
    # coefficient applied to the TOTAL lung amount, S1 Text
    # `CVLU = (ALUI+ALUT+ALUB)/(VLU*PLU)` and `CLU = (ALUI+ALUT+ALUB)/VLU`.
    a_lung_total <- vp_lung + is_lung + int_lung
    cv_lung <- a_lung_total / (v_lung * kp_lung)
    c_lung <- a_lung_total / v_lung
    c_is_lung <- is_lung / v_is_lung                  # S1 Text `CLUI = ALUI/VLUI`

    # =================================================================
    # Fluxes (mg/h)
    # =================================================================
    r_absorb <- ka * depot                            # S1 Text `Rim = Kim*Amtsiteim`
    r_renal <- cl_renal * cv_kidney                   # S1 Text `Rurine = Kurine*CVK`
    r_bile <- cl_nonren * cv_liver                    # S1 Text `Rmet = Kbile*CVL`

    # Lung permeability fluxes. Each transfer constant multiplies a
    # sub-compartment CONCENTRATION (not an amount), exactly as printed
    # in the S1 Text code; only unbound drug crosses, so the vascular and
    # tissue terms carry fu and (1 - PT) respectively, while interstitial
    # fluid is treated as entirely unbound (it is what microdialysis
    # samples). See vignette Errata E2.
    r_blood_isf <- k_blood_isf * (vp_lung * fu / v_vp_lung)
    # S1 Text `KBI*(ALUB*(1-PB)/VLUB)`
    r_isf_blood <- k_isf_blood * is_lung / v_is_lung
    # S1 Text `KIB*ALUI/VLUI`
    r_isf_tissue <- k_isf_tissue * is_lung / v_is_lung
    # S1 Text `KIT*ALUI/VLUI`
    r_tissue_isf <- k_tissue_isf * (int_lung * (1 - fb_lung) / v_int_lung)
    # S1 Text `KTI*(ALUT*(1-PBtissue)/VLUT)`

    # =================================================================
    # ODE system (S1 Text, {CEQ distribution in each compartment})
    # =================================================================
    # Intramuscular absorption. `depot` is the injection-site pool
    # (S1 Text `Amtsiteim`) and `depot2` the slow-release depot
    # (S1 Text `DOSEppgim`). The dose is split between them by the
    # bioavailability terms below rather than by an in-ODE input rate.
    d/dt(depot2) <- -kdiss * depot2
    d/dt(depot) <- kdiss * depot2 - r_absorb

    # Venous blood receives every organ's venous outflow plus the
    # absorbed intramuscular drug, S1 Text
    # `RV = (QL*CVL+QK*CVK+QM*CVM+QR*CVR+Rim)-QC*CV`.
    d/dt(venous) <-
      q_liver * cv_liver +
      q_kidney * cv_kidney +
      q_muscle * cv_muscle +
      q_other * cv_other +
      r_absorb -
      q_co * c_venous

    # Arterial blood is fed by the lung and drained by the tissues, which
    # take only unbound drug, S1 Text `RA = QC*(CVLu-CAfree)`.
    d/dt(arterial) <- q_co * (cv_lung - c_arterial_free)

    # Perfusion-limited, well-stirred organs.
    d/dt(liver) <- q_liver * (c_arterial_free - cv_liver) - r_bile
    d/dt(kidney) <- q_kidney * (c_arterial_free - cv_kidney) - r_renal
    d/dt(muscle) <- q_muscle * (c_arterial_free - cv_muscle)
    d/dt(other) <- q_other * (c_arterial_free - cv_other)

    # Permeability-limited lung: vascular blood -> interstitial fluid ->
    # tissue, S1 Text `RLUB`, `RLUI`, `RLUT`.
    d/dt(vp_lung) <- q_co * (c_venous - cv_lung) - r_blood_isf + r_isf_blood
    d/dt(is_lung) <- r_blood_isf - r_isf_blood - r_isf_tissue + r_tissue_isf
    d/dt(int_lung) <- r_isf_tissue - r_tissue_isf

    # Terminal excretion sinks, carried so the published mass balance can
    # be checked, S1 Text `d/dt(Aurine)` and `d/dt(Amet)`.
    d/dt(urine) <- r_renal
    d/dt(bile) <- r_bile

    # Split the intramuscular dose: 1 - Frac is immediately available at
    # the injection site, Frac goes to the slow-release depot. S1 Text
    # `Rpenim = Rinputim*(1-Frac)` and `Rppgim = Rinputim*Frac`. Dose the
    # same amount to both compartments and let these fractions apportion
    # it.
    f(depot) <- 1 - frac
    f(depot2) <- frac

    # =================================================================
    # Observations
    # =================================================================
    # Mi 2023 compares model output against plasma, lung interstitial
    # fluid, kidney and liver. The S1 Text code makes no
    # blood-to-plasma correction, so the venous blood concentration is
    # the paper's plasma prediction.
    Cc <- c_venous
    Cliver <- liver / v_liver
    Ckidney <- kidney / v_kidney
    Cmuscle <- muscle / v_muscle
    Clung <- c_lung
    # Free drug in lung interstitial fluid, i.e. the pulmonary epithelial
    # lining fluid (PELF) concentration that microdialysis measures and
    # that drives the %fT>MIC target.
    Cisf_lung <- c_is_lung
    # Unbound arterial concentration, the other %fT>MIC driver
    # (S1 Text Eq 5 evaluates T>MIC on CAfree).
    Cfree <- c_arterial_free

    Cc ~ prop(propSd)
  })
}
