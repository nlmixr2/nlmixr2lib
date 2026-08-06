Mody_2023_doxorubicin_dexrazoxane_mdamb468 <- function() {
  description <- "In vitro (MDA-MB-468 human breast cancer cell line; triple-negative breast cancer). Cellular-level pharmacodynamic (PD) model for the single-agent and combined anticancer effects of doxorubicin (DOX) and dexrazoxane (DEX) on percent cell viability over 96 h (Mody 2023 Eqs 3, 3a-3e and 4a-4i; Table 2, MDA-MB-468 column). Both drugs degrade in the cell culture medium by first-order kinetics, so their driving concentrations are dosable one-compartment states conc_dox and conc_dexrazoxane that decay with kdeg (Mody 2023 Sci Rep Eqs 1-2, rate constants restated in this paper's Results). Each drug stimulates cell death through a capacity-limited Hill function whose signal is delayed by a three-compartment signal-transduction transit chain (transit1_<drug> to transit3_<drug>, first-order rate ktr = 1/tau). The two terminal transit signals act ADDITIVELY on the loss term of an exponentially growing viability state (Eq 4i). The drug-drug interaction parameter psi multiplies the DOX half-maximal killing concentration kc50_dox ONLY (Eq 4a); psi > 1 is antagonistic, psi < 1 synergistic and psi = 1 additive. Estimated psi is 1 (additive) in MDA-MB-468. Setting psi = 1 recovers the paper's single-agent DOX equations (Eq 3a-3e) exactly, which is how the single-agent arms were fitted. No IIV or residual error was estimated (Table 2 reports point estimates with %RSE only); an arbitrary 10% CV IIV is introduced at simulation time per the paper's Methods. Sibling models: Mody_2023_doxorubicin_dexrazoxane_jimt1 (the other cell line), Mody_2023_doxorubicin_dexrazoxane_static_mdamb468 (the 72-h static concentration-response fit) and Mody_2023_doxorubicin_dexrazoxane_clinical_mdamb468 (the same PD model driven by clinical tumour-site concentrations)."
  reference <- paste(
    "Mody H, Vaidya TR, Lezeau J, Taha K, Ait-Oudhia S (2023).",
    "In vitro to clinical translation of combinatorial effects of doxorubicin",
    "and dexrazoxane in breast cancer: a mechanism-based pharmacokinetic/",
    "pharmacodynamic modeling approach.",
    "Frontiers in Pharmacology 14:1239141.",
    "doi:10.3389/fphar.2023.1239141.",
    sep = " "
  )
  vignette <- "Mody_2023_doxorubicin_dexrazoxane_breast_cancer"

  units <- list(
    time          = "h",
    dosing        = "uM (drug concentration applied to the culture medium at t = 0)",
    concentration = "uM (medium drug concentration); % (cell viability, the PD readout)"
  )

  compartmentData <- list(
    conc_dox                = list(analyte = "doxorubicin", units = "uM", specimen = "administration site", verified = TRUE),
    conc_dexrazoxane        = list(analyte = "dexrazoxane", units = "uM", specimen = "administration site", verified = TRUE),
    transit1_dox            = list(analyte = "doxorubicin cell-death signal", units = "1/h", specimen = "not applicable", verified = TRUE),
    transit2_dox            = list(analyte = "doxorubicin cell-death signal", units = "1/h", specimen = "not applicable", verified = TRUE),
    transit3_dox            = list(analyte = "doxorubicin cell-death signal", units = "1/h", specimen = "not applicable", verified = TRUE),
    transit1_dexrazoxane    = list(analyte = "dexrazoxane cell-death signal", units = "1/h", specimen = "not applicable", verified = TRUE),
    transit2_dexrazoxane    = list(analyte = "dexrazoxane cell-death signal", units = "1/h", specimen = "not applicable", verified = TRUE),
    transit3_dexrazoxane    = list(analyte = "dexrazoxane cell-death signal", units = "1/h", specimen = "not applicable", verified = TRUE),
    viability               = list(analyte = "cell viability", units = "%", specimen = "not applicable", verified = TRUE)
  )

  # conc_dox / conc_dexrazoxane are the cell-culture-medium drug
  # concentration states. `conc_<drug>` follows the in-vitro precedent of
  # Kroemer_2024_ceftazidime_avibactam_fosfomycin_hfim.R; the bare `conc`
  # prefix is not a registered canonical compartment, so both are declared
  # paper-specific here rather than registered.
  paper_specific_compartments <- c("conc_dox", "conc_dexrazoxane")

  covariateData <- list()

  population <- list(
    species        = "in vitro (MDA-MB-468 human breast cancer cell line)",
    n_subjects     = NA_integer_,
    n_studies      = 1L,
    disease_state  = "triple-negative breast cancer. MDA-MB-468 cells seeded at 5 x 10^3 cells per well (100 uL) of a 96-well plate, incubated overnight, then exposed to DOX, DEX or their combination and assayed for viability with CCK-8 (absorbance 450 nm).",
    dose_range     = "DOX 0.005-1 uM, DEX 6.25-200 uM, and 36 DOX x DEX combinations of those levels.",
    cell_line      = "MDA-MB-468 (ATCC, Manassas, VA); triple-negative breast cancer.",
    culture        = "DMEM with 10% fetal bovine serum and 1% penicillin/streptomycin; 37 C, humidified 5% CO2; passaged at confluency with 0.25% trypsin / 2.21 nM EDTA.",
    notes          = "No between-well or between-experiment random effects were estimated; Table 2 reports point estimates with %RSE only. Experiments were run in at least triplicate against matched vehicle controls. The paper introduces an arbitrary 10% CV IIV on the PD parameters only at the clinical-simulation stage (see Mody_2023_doxorubicin_dexrazoxane_clinical_mdamb468)."
  )

  ini({
    # ================================================================
    # MEDIUM DEGRADATION (Mody 2023 Sci Rep Eqs 1-2; rate constants
    # restated verbatim in this paper's Results, "Time course effects
    # of DOX and DEX as single agents and in combination")
    # ================================================================
    # Carried from the companion paper and held FIXED here: this paper
    # did not re-estimate them ("previously estimated at ... were used").
    lkdeg_dox         <- fixed(log(0.022))  ; label("First-order degradation rate constant of DOX in cell culture medium (1/h)")   # Results: 0.022 (+/-0.0004) 1/h
    lkdeg_dexrazoxane <- fixed(log(0.054))  ; label("First-order degradation rate constant of DEX in cell culture medium (1/h)")   # Results: 0.054 (+/-0.0016) 1/h

    # ================================================================
    # CELLULAR-LEVEL PD PARAMETERS -- Mody 2023 Table 2, MDA-MB-468 column
    # ================================================================
    lrbase    <- fixed(log(100))    ; label("Baseline percent cell viability R0 (%)")                                        # Table 2: 100 (Fixed)
    lkgrow    <- log(0.0167)       ; label("First-order cell growth rate constant kg (1/h)")                                # Table 2: 0.0167 (1.45% RSE)

    lkmax_dox <- log(1.16)    ; label("Maximal DOX killing rate constant Smax,DOX (1/h)")                              # Table 2: 1.16 (3.55% RSE)
    lkc50_dox <- log(0.761)    ; label("DOX concentration for half-maximal killing rate SC50,DOX (uM)")                 # Table 2: 0.761 (7.48% RSE)
    lktr_dox  <- log(0.0275)     ; label("DOX cell-death transit rate constant 1/tau_DOX (1/h)")                          # Table 2: 0.0275 (3.17% RSE); tau_DOX = 36.4 h

    lkmax_dexrazoxane <- log(0.29) ; label("Maximal DEX killing rate constant Smax,DEX (1/h)")                         # Table 2: 0.29 (48.5% RSE)
    lkc50_dexrazoxane <- log(17.5) ; label("DEX concentration for half-maximal killing rate SC50,DEX (uM)")            # Table 2: 17.5 (23.7% RSE)
    lktr_dexrazoxane  <- log(0.0182)  ; label("DEX cell-death transit rate constant 1/tau_DEX (1/h)")                     # Table 2: 0.0182 (19.4% RSE); tau_DEX = 54.9 h

    # Table 2 prints psi only to one significant figure, as "~1"
    # (%RSE 0.542); no more precise value exists in the paper or its
    # supplement. Encoded as the printed value. psi = 1 recovers the
    # single-agent DOX model of Eqs 3a-3e.
    lpsi <- log(1) ; label("DOX-DEX interaction parameter psi (unitless; >1 antagonistic, <1 synergistic, =1 additive)")  # Table 2: ~1 (0.542% RSE)

    # ================================================================
    # RESIDUAL ERROR
    # ================================================================
    # Mody 2023 Table 2 reports point estimates with %RSE only; no
    # residual-error magnitude is reported anywhere in the paper or its
    # supplement. Encoded as fixed(0) rather than an invented value.
    addSd_viability <- fixed(0) ; label("Additive residual SD on percent cell viability (%; ZERO - not reported in source)")
  })

  model({
    kdeg_dox         <- exp(lkdeg_dox)
    kdeg_dexrazoxane <- exp(lkdeg_dexrazoxane)
    rbase            <- exp(lrbase)
    kgrow            <- exp(lkgrow)
    kmax_dox         <- exp(lkmax_dox)
    kc50_dox         <- exp(lkc50_dox)
    ktr_dox          <- exp(lktr_dox)
    kmax_dexrazoxane <- exp(lkmax_dexrazoxane)
    kc50_dexrazoxane <- exp(lkc50_dexrazoxane)
    ktr_dexrazoxane  <- exp(lktr_dexrazoxane)
    psi              <- exp(lpsi)

    # ---- Medium drug degradation (Sci Rep Eqs 1-2) ----------------
    # The applied concentration enters as a bolus "dose" (in uM) into
    # each state at t = 0, so the state IS the medium concentration.
    d/dt(conc_dox)         <- -kdeg_dox         * conc_dox
    d/dt(conc_dexrazoxane) <- -kdeg_dexrazoxane * conc_dexrazoxane

    # ---- Capacity-limited stimulation of cell death ---------------
    # Eq 4a: psi scales the DOX half-maximal killing concentration ONLY.
    # Eq 4e: the DEX arm carries no interaction term.
    kkill_dox         <- kmax_dox         * conc_dox         / (kc50_dox * psi     + conc_dox)
    kkill_dexrazoxane <- kmax_dexrazoxane * conc_dexrazoxane / (kc50_dexrazoxane   + conc_dexrazoxane)

    # ---- Three-compartment signal-transduction transit chains -----
    # Eqs 4b-4d (DOX) and 4f-4h (DEX); all chains start empty.
    d/dt(transit1_dox) <- ktr_dox * (kkill_dox     - transit1_dox)
    d/dt(transit2_dox) <- ktr_dox * (transit1_dox  - transit2_dox)
    d/dt(transit3_dox) <- ktr_dox * (transit2_dox  - transit3_dox)

    d/dt(transit1_dexrazoxane) <- ktr_dexrazoxane * (kkill_dexrazoxane    - transit1_dexrazoxane)
    d/dt(transit2_dexrazoxane) <- ktr_dexrazoxane * (transit1_dexrazoxane - transit2_dexrazoxane)
    d/dt(transit3_dexrazoxane) <- ktr_dexrazoxane * (transit2_dexrazoxane - transit3_dexrazoxane)

    # ---- Exponential growth with additive drug-induced death ------
    # Eq 4i. With conc_dexrazoxane = 0 this reduces to Eq 3e for DOX
    # alone (and vice versa).
    d/dt(viability) <- kgrow * viability -
                       (transit3_dox + transit3_dexrazoxane) * viability
    viability(0)    <- rbase

    viability ~ add(addSd_viability)
  })
}
