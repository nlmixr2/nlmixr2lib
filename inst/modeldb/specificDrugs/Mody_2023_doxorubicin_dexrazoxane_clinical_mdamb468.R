Mody_2023_doxorubicin_dexrazoxane_clinical_mdamb468 <- function() {
  description <- "In vitro to clinical translation (MDA-MB-468 human breast cancer cell line; triple-negative breast cancer) of the doxorubicin (DOX) + dexrazoxane (DEX) anticancer pharmacodynamic model (Mody 2023 Methods, 'In vitro-in vivo translation of PD responses'; Figure 6). The in-vitro cellular PD model of Mody_2023_doxorubicin_dexrazoxane_mdamb468 is driven here by clinically relevant TUMOUR-SITE concentrations instead of cell-culture-medium concentrations. Clinical plasma PK is reproduced from the companion paper (Mody 2023 Sci Rep 13:3100, doi:10.1038/s41598-023-29964-4): DOX is a 3-compartment mammillary model with linear elimination and DEX a 2-compartment mammillary model parameterised by rate constants, both given as 15-min IV infusions for a typical 1.8 m^2 body-surface-area subject. Tumour-site concentrations are obtained by scaling the simulated plasma concentrations by a total tumour:plasma AUC ratio, kp_tumor_dox = 57.1 (digitised from He 2018 human tumour distribution data) and kp_tumor_dexrazoxane, which is ARBITRARY because no DEX tumour-distribution data exist (the paper simulates 0.1, 1 and 10). Plasma concentrations in mg/L are converted to uM at the PK-PD interface using literature molar masses. The paper's simulated regimen is DOX 50 mg/m^2 Q3W with DEX 500 mg/m^2 (a 10:1 DEX:DOX ratio, previously identified as maximally cardioprotective) over three cycles. For MDA-MB-468 the paper predicts a COMPARABLE AUEC of percent cell viability across the DOX and DOX + DEX groups, consistent with the additive in-vitro interaction. An arbitrary 10% CV inter-individual variability is applied to the PD parameters and to both tumour:plasma ratios, exactly as the paper's 500-subject population simulation does; these are simulation conventions, not fitted variances, and are wrapped in fixed(). Sibling models: Mody_2023_doxorubicin_dexrazoxane_clinical_jimt1, Mody_2023_doxorubicin_dexrazoxane_mdamb468 (the in-vitro fit these PD parameters come from) and Mody_2023_doxorubicin_dexrazoxane_static_mdamb468."
  reference <- paste(
    "Mody H, Vaidya TR, Lezeau J, Taha K, Ait-Oudhia S (2023).",
    "In vitro to clinical translation of combinatorial effects of doxorubicin",
    "and dexrazoxane in breast cancer: a mechanism-based pharmacokinetic/",
    "pharmacodynamic modeling approach.",
    "Frontiers in Pharmacology 14:1239141.",
    "doi:10.3389/fphar.2023.1239141.",
    "Clinical plasma PK parameters are from the companion paper",
    "Mody H, Vaidya TR, Ait-Oudhia S (2023), Scientific Reports 13:3100,",
    "doi:10.1038/s41598-023-29964-4, Table 2.",
    sep = " "
  )
  vignette <- "Mody_2023_doxorubicin_dexrazoxane_breast_cancer"

  units <- list(
    time          = "h",
    dosing        = "mg (IV 15-min infusion; mg/m^2 doses pre-multiplied by the 1.8 m^2 typical body surface area in the event table)",
    concentration = "mg/L (plasma); uM (tumour-site PD driver, derived); % (cell viability, the PD readout)"
  )

  compartmentData <- list(
    central_dox              = list(analyte = "doxorubicin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_dox          = list(analyte = "doxorubicin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2_dox          = list(analyte = "doxorubicin", units = "mg", specimen = "plasma", verified = TRUE),
    central_dexrazoxane      = list(analyte = "dexrazoxane", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_dexrazoxane  = list(analyte = "dexrazoxane", units = "mg", specimen = "plasma", verified = TRUE),
    transit1_dox             = list(analyte = "doxorubicin cell-death signal", units = "1/h", specimen = "not applicable", verified = TRUE),
    transit2_dox             = list(analyte = "doxorubicin cell-death signal", units = "1/h", specimen = "not applicable", verified = TRUE),
    transit3_dox             = list(analyte = "doxorubicin cell-death signal", units = "1/h", specimen = "not applicable", verified = TRUE),
    transit1_dexrazoxane     = list(analyte = "dexrazoxane cell-death signal", units = "1/h", specimen = "not applicable", verified = TRUE),
    transit2_dexrazoxane     = list(analyte = "dexrazoxane cell-death signal", units = "1/h", specimen = "not applicable", verified = TRUE),
    transit3_dexrazoxane     = list(analyte = "dexrazoxane cell-death signal", units = "1/h", specimen = "not applicable", verified = TRUE),
    viability                = list(analyte = "cell viability", units = "%", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species          = "human (simulated) driving an in vitro MDA-MB-468 human breast cancer cell PD model",
    n_subjects       = 500L,
    n_studies        = 1L,
    disease_state    = "Simulated typical adult cancer patient with a body surface area of 1.8 m^2; the PD layer is the MDA-MB-468 (triple-negative breast cancer) cell model established in vitro.",
    dose_range       = "DOX 50 mg/m^2 (90 mg at 1.8 m^2) Q3W as a 15-min IV infusion, alone or with DEX 500 mg/m^2 (900 mg at 1.8 m^2), over three cycles (9 weeks).",
    dosing_frequency = "Q3W, three cycles, both drugs given simultaneously as 15-min IV infusions.",
    notes            = "Population simulation of 500 subjects with an arbitrary 10% inter-individual variability on the PD parameters and on the tumour:plasma ratios (paper Methods). Efficacy was summarised as the area under the percent-cell-viability effect curve (AUEC) over the three cycles for DOX alone and for DOX + DEX at kp_tumor_dexrazoxane of 0.1, 1 and 10. The clinical plasma PK parameters are typical values for a 1.8 m^2 subject and carry no covariates."
  )

  ini({
    # ================================================================
    # DOXORUBICIN CLINICAL PLASMA PK -- 3-compartment, IV
    # ================================================================
    # Companion paper Mody 2023 Sci Rep Table 2 (top). The underlying
    # PK model is Kontny et al. 2013; the companion reproduces the
    # typical values for a 1.8 m^2 body-surface-area subject, which is
    # the form used here. Fixed: this paper did not re-estimate them.
    lcl_dox  <- fixed(log(53.3)) ; label("Doxorubicin clearance CL (L/h per 1.8 m^2)")                          # Sci Rep Table 2: 53.3 (4% RSE)
    lvc_dox  <- fixed(log(17.7)) ; label("Doxorubicin central volume Vc (L per 1.8 m^2)")                       # Sci Rep Table 2: 17.7 (8% RSE)
    lq_dox   <- fixed(log(58.7)) ; label("Doxorubicin inter-compartmental clearance Q2 (L/h per 1.8 m^2)")      # Sci Rep Table 2: 58.7 (8% RSE)
    lvp_dox  <- fixed(log(1830)) ; label("Doxorubicin first peripheral volume V2 (L per 1.8 m^2)")              # Sci Rep Table 2: 1830 (7% RSE)
    lq2_dox  <- fixed(log(21.8)) ; label("Doxorubicin inter-compartmental clearance Q3 (L/h per 1.8 m^2)")      # Sci Rep Table 2: 21.8 (13% RSE)
    lvp2_dox <- fixed(log(71.6)) ; label("Doxorubicin second peripheral volume V3 (L per 1.8 m^2)")             # Sci Rep Table 2: 71.6 (15% RSE)

    # ================================================================
    # DEXRAZOXANE CLINICAL PLASMA PK -- 2-compartment, IV
    # ================================================================
    # Companion paper Sci Rep Table 2 (bottom), fitted to the phase 1
    # data of Earhart et al. 1982. Reported directly as rate constants
    # (all three estimated as exactly 1 /h with distinct %RSEs of 4, 3
    # and 8), so the rate-constant parameterisation is kept rather than
    # being converted to CL/Q macro-constants.
    lkel_dexrazoxane <- fixed(log(1))    ; label("Dexrazoxane elimination rate constant kel (1/h)")                       # Sci Rep Table 2: 1 (4% RSE)
    lk12_dexrazoxane <- fixed(log(1))    ; label("Dexrazoxane central-to-peripheral rate constant k12 (1/h)")             # Sci Rep Table 2: 1 (3% RSE)
    lk21_dexrazoxane <- fixed(log(1))    ; label("Dexrazoxane peripheral-to-central rate constant k21 (1/h)")             # Sci Rep Table 2: 1 (8% RSE)
    lvc_dexrazoxane  <- fixed(log(14.6)) ; label("Dexrazoxane central volume Vc (L)")                                     # Sci Rep Table 2: 14.6 (5% RSE)

    # ================================================================
    # TUMOUR:PLASMA PENETRATION
    # ================================================================
    lkp_tumor_dox <- fixed(log(57.1)) ; label("Total tumour:plasma AUC ratio for DOX (unitless)")   # Results: AUCtumor/AUCplasma = 57.1, digitised from He 2018
    # ARBITRARY. The paper states no human tumour-distribution data
    # exist for DEX and simulates 0.1, 1 and 10 as scenarios. The
    # middle scenario is shipped as the default; change it to run the
    # other two. The paper found the choice barely affects the AUEC.
    lkp_tumor_dexrazoxane <- fixed(log(1)) ; label("Total tumour:plasma AUC ratio for DEX (unitless; ARBITRARY, paper simulates 0.1, 1 and 10)")   # Methods: arbitrary values 0.1, 1 and 10

    # ================================================================
    # PLASMA mg/L -> uM CONVERSION (literature molar masses, FIXED)
    # ================================================================
    # NOT from the paper. The PD parameters were fitted in uM while the
    # clinical PK is in mg/L, so the interface needs a molar mass:
    #   c_uM = C (mg/L) * 1000 / MW (g/mol)
    # DOX free base MW 543.52 g/mol (PubChem CID 31703) -> 1.8397
    # DEX free base MW 268.27 g/mol (PubChem CID 71384) -> 3.7275
    lconv_dox         <- fixed(log(1.8397)) ; label("DOX mg/L to uM conversion factor (unitless; literature molar mass, not from the paper)")   # PubChem CID 31703, MW 543.52 g/mol
    lconv_dexrazoxane <- fixed(log(3.7275)) ; label("DEX mg/L to uM conversion factor (unitless; literature molar mass, not from the paper)")   # PubChem CID 71384, MW 268.27 g/mol

    # ================================================================
    # CELLULAR-LEVEL PD PARAMETERS -- Mody 2023 Table 2, MDA-MB-468 column
    # ================================================================
    lrbase    <- fixed(log(100))    ; label("Baseline percent cell viability R0 (%)")                          # Table 2: 100 (Fixed)
    lkgrow    <- log(0.0167)       ; label("First-order cell growth rate constant kg (1/h)")                  # Table 2, MDA-MB-468
    lkmax_dox <- log(1.16)    ; label("Maximal DOX killing rate constant Smax,DOX (1/h)")                # Table 2, MDA-MB-468
    lkc50_dox <- log(0.761)    ; label("DOX concentration for half-maximal killing rate SC50,DOX (uM)")   # Table 2, MDA-MB-468
    lktr_dox  <- log(0.0275)     ; label("DOX cell-death transit rate constant 1/tau_DOX (1/h)")            # Table 2, MDA-MB-468
    lkmax_dexrazoxane <- log(0.29) ; label("Maximal DEX killing rate constant Smax,DEX (1/h)")           # Table 2, MDA-MB-468
    lkc50_dexrazoxane <- log(17.5) ; label("DEX concentration for half-maximal killing rate SC50,DEX (uM)")  # Table 2, MDA-MB-468
    lktr_dexrazoxane  <- log(0.0182)  ; label("DEX cell-death transit rate constant 1/tau_DEX (1/h)")       # Table 2, MDA-MB-468
    lpsi <- log(1) ; label("DOX-DEX interaction parameter psi (unitless; >1 antagonistic, <1 synergistic, =1 additive)")  # Table 2, MDA-MB-468: ~1 (additive)

    # ================================================================
    # ARBITRARY SIMULATION IIV (paper Methods, not fitted)
    # ================================================================
    # "Population simulations were conducted for 500 subjects with an
    # arbitrary inter-individual variability (IIV) of 10% introduced on
    # the PD parameters and the Fac1 and Fac2 parameters." These are
    # simulation conventions, so every one is fixed(); 0.01 is the
    # log-scale variance giving an approximately 10% CV.
    etalkgrow            ~ fixed(0.01)
    etalkmax_dox         ~ fixed(0.01)
    etalkc50_dox         ~ fixed(0.01)
    etalktr_dox          ~ fixed(0.01)
    etalkmax_dexrazoxane ~ fixed(0.01)
    etalkc50_dexrazoxane ~ fixed(0.01)
    etalktr_dexrazoxane  ~ fixed(0.01)
    etalpsi              ~ fixed(0.01)
    etalkp_tumor_dox         ~ fixed(0.01)
    etalkp_tumor_dexrazoxane ~ fixed(0.01)

    # ================================================================
    # RESIDUAL ERROR
    # ================================================================
    # Neither paper reports a residual-error magnitude for plasma DOX,
    # plasma DEX or cell viability. Encoded as fixed(0), not invented.
    addSd_dox         <- fixed(0) ; label("Additive residual SD on plasma DOX (mg/L; ZERO - not reported in source)")
    addSd_dexrazoxane <- fixed(0) ; label("Additive residual SD on plasma DEX (mg/L; ZERO - not reported in source)")
    addSd_viability      <- fixed(0) ; label("Additive residual SD on percent cell viability (%; ZERO - not reported in source)")
  })

  model({
    cl_dox  <- exp(lcl_dox)
    vc_dox  <- exp(lvc_dox)
    q_dox   <- exp(lq_dox)
    vp_dox  <- exp(lvp_dox)
    q2_dox  <- exp(lq2_dox)
    vp2_dox <- exp(lvp2_dox)
    kel_dexrazoxane <- exp(lkel_dexrazoxane)
    k12_dexrazoxane <- exp(lk12_dexrazoxane)
    k21_dexrazoxane <- exp(lk21_dexrazoxane)
    vc_dexrazoxane  <- exp(lvc_dexrazoxane)

    kp_tumor_dox         <- exp(lkp_tumor_dox         + etalkp_tumor_dox)
    kp_tumor_dexrazoxane <- exp(lkp_tumor_dexrazoxane + etalkp_tumor_dexrazoxane)
    conv_dox             <- exp(lconv_dox)
    conv_dexrazoxane     <- exp(lconv_dexrazoxane)

    rbase            <- exp(lrbase)
    kgrow            <- exp(lkgrow            + etalkgrow)
    kmax_dox         <- exp(lkmax_dox         + etalkmax_dox)
    kc50_dox         <- exp(lkc50_dox         + etalkc50_dox)
    ktr_dox          <- exp(lktr_dox          + etalktr_dox)
    kmax_dexrazoxane <- exp(lkmax_dexrazoxane + etalkmax_dexrazoxane)
    kc50_dexrazoxane <- exp(lkc50_dexrazoxane + etalkc50_dexrazoxane)
    ktr_dexrazoxane  <- exp(lktr_dexrazoxane  + etalktr_dexrazoxane)
    psi              <- exp(lpsi              + etalpsi)

    # ---- DOX micro-constants --------------------------------------
    kel_dox <- cl_dox / vc_dox
    k12_dox <- q_dox  / vc_dox
    k21_dox <- q_dox  / vp_dox
    k13_dox <- q2_dox / vc_dox
    k31_dox <- q2_dox / vp2_dox

    # ---- Clinical plasma PK ---------------------------------------
    d/dt(central_dox)     <- -kel_dox * central_dox -
                              k12_dox * central_dox + k21_dox * peripheral1_dox -
                              k13_dox * central_dox + k31_dox * peripheral2_dox
    d/dt(peripheral1_dox) <-  k12_dox * central_dox - k21_dox * peripheral1_dox
    d/dt(peripheral2_dox) <-  k13_dox * central_dox - k31_dox * peripheral2_dox

    d/dt(central_dexrazoxane)     <- -kel_dexrazoxane * central_dexrazoxane -
                                      k12_dexrazoxane * central_dexrazoxane +
                                      k21_dexrazoxane * peripheral1_dexrazoxane
    d/dt(peripheral1_dexrazoxane) <-  k12_dexrazoxane * central_dexrazoxane -
                                      k21_dexrazoxane * peripheral1_dexrazoxane

    Cc_dox         <- central_dox         / vc_dox
    Cc_dexrazoxane <- central_dexrazoxane / vc_dexrazoxane

    # ---- Plasma -> tumour site -> uM PD driver --------------------
    ctumor_dox         <- Cc_dox         * kp_tumor_dox         * conv_dox
    ctumor_dexrazoxane <- Cc_dexrazoxane * kp_tumor_dexrazoxane * conv_dexrazoxane

    # ---- Cellular PD (Mody 2023 Eqs 4a-4i) ------------------------
    kkill_dox         <- kmax_dox         * ctumor_dox         / (kc50_dox * psi   + ctumor_dox)
    kkill_dexrazoxane <- kmax_dexrazoxane * ctumor_dexrazoxane / (kc50_dexrazoxane + ctumor_dexrazoxane)

    d/dt(transit1_dox) <- ktr_dox * (kkill_dox    - transit1_dox)
    d/dt(transit2_dox) <- ktr_dox * (transit1_dox - transit2_dox)
    d/dt(transit3_dox) <- ktr_dox * (transit2_dox - transit3_dox)

    d/dt(transit1_dexrazoxane) <- ktr_dexrazoxane * (kkill_dexrazoxane    - transit1_dexrazoxane)
    d/dt(transit2_dexrazoxane) <- ktr_dexrazoxane * (transit1_dexrazoxane - transit2_dexrazoxane)
    d/dt(transit3_dexrazoxane) <- ktr_dexrazoxane * (transit2_dexrazoxane - transit3_dexrazoxane)

    d/dt(viability) <- kgrow * viability -
                       (transit3_dox + transit3_dexrazoxane) * viability
    viability(0)    <- rbase

    Cc_dox         ~ add(addSd_dox)
    Cc_dexrazoxane ~ add(addSd_dexrazoxane)
    viability      ~ add(addSd_viability)
  })
}
