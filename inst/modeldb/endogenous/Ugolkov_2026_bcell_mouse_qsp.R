Ugolkov_2026_bcell_mouse_qsp <- function() {
  description <- paste(
    "QSP. Preclinical (mouse). Ugolkov 2026 quantitative systems pharmacology",
    "model of the T cell-dependent B cell immune response to antigen exposure.",
    "Twenty ODEs and 31 parameters spanning two coupled layers. A homeostasis",
    "sub-model (7 ODEs) carries immature B cells generated in the bone marrow",
    "through transitional type 1 (T1) B cells in bone, blood and spleen into",
    "naive B cells recirculating between spleen, blood and lymph nodes. An",
    "activation layer (13 ODEs) drives antibody-secreting cell (ASC)",
    "generation from an empirical antigen forcing function passed through a",
    "six-compartment transit chain (mean transit time 18 days), generates ASC",
    "in spleen and lymph nodes in proportion to the local antigen-specific",
    "naive B cell precursor pool (frequency 1 in 10^6), and distributes ASC",
    "via blood into a saturable bone-marrow survival niche and into",
    "non-lymphoid peripheral tissues. ASC in spleen, lymph nodes, bone marrow",
    "and blood drive antigen-specific IgG production in plasma. ASC",
    "elimination is deliberately absent: the source reports ASC levels stable",
    "to 300 days post-immunization. Deterministic: the source reports no IIV",
    "and tabulates no residual-error magnitudes. Of the 31 parameters, 19 are",
    "fixed from experimental data, 7 are derived analytically from the",
    "steady-state cell counts, and 5 (k9, k12, k13, kmat_3, Vmax) were",
    "estimated (all RSE <= 24%). Two Table 1 values are NOT reproducible and",
    "the deposited RxODE code is used instead; see the vignette Errata.",
    sep = " "
  )
  reference <- paste(
    "Ugolkov Y, Volkova A, Helmlinger G, Peskov K, Sokolov V (2026).",
    "Quantitative systems pharmacology model of B cell immune response in",
    "mouse. Front Immunol 17:1745710. doi:10.3389/fimmu.2026.1745710.",
    "PMCID: PMC13171755.",
    "Parameter values from Table 1; model equations from Supplementary",
    "Equations 1-28; initial conditions, the antigen forcing function",
    "constants and the mean transit time from the Supplementary RxODE model",
    "code (Supplementary Material DataSheet1.pdf).",
    sep = " "
  )
  vignette <- "Ugolkov_2026_bcell_mouse_qsp"

  # Every state is a paper-mechanistic B cell / ASC / antigen / antibody
  # species; none maps onto a canonical PK compartment role. Names are the
  # lower-case snake_case forms of the Supplementary Equations' state names
  # (ImmBone -> imm_bone, ASCSpleen -> asc_spleen, and so on).
  paper_specific_compartments <- c(
    "imm_bone", "t1_bone", "t1_blood", "t1_spleen",
    "naive_blood", "naive_spleen", "naive_ln",
    "antigen0", "antigen1", "antigen2", "antigen3", "antigen4", "antigen5",
    "antigen",
    "asc_spleen", "asc_ln", "asc_bone", "asc_blood", "asc_peripheral",
    "igg_blood"
  )

  units <- list(
    time          = "day",
    dosing        = "none (no drug is administered; antigen exposure enters as the built-in empirical forcing function of Supplementary Equation 15, not as a dosing event)",
    concentration = "cell (ASC states, absolute cell counts) and 10^6 cell (B cell states); IgG in pg"
  )

  covariateData <- list()

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Supplementary Equations 8-28, the
  # Table 1 "Description" and "Units" columns, and the Methods narrative.
  #
  # NOTE the two amount scales, which the source code carries implicitly via
  # its "#cell / 10^6" comments: the seven homeostasis B cell states are in
  # MILLIONS of cells, while the five ASC states are ABSOLUTE cell counts.
  # The conversion happens exactly once, in the ASC generation terms of
  # Supplementary Equations 22-23 -- see the model() block.
  compartmentData <- list(
    imm_bone       = list(analyte = "immature B cells",             units = "10^6 cell", specimen = "tissue",         verified = TRUE),
    t1_bone        = list(analyte = "transitional type 1 B cells",  units = "10^6 cell", specimen = "tissue",         verified = TRUE),
    t1_blood       = list(analyte = "transitional type 1 B cells",  units = "10^6 cell", specimen = "whole blood",    verified = TRUE),
    t1_spleen      = list(analyte = "transitional type 1 B cells",  units = "10^6 cell", specimen = "tissue",         verified = TRUE),
    naive_blood    = list(analyte = "naive B cells",                units = "10^6 cell", specimen = "whole blood",    verified = TRUE),
    naive_spleen   = list(analyte = "naive B cells",                units = "10^6 cell", specimen = "tissue",         verified = TRUE),
    naive_ln       = list(analyte = "naive B cells",                units = "10^6 cell", specimen = "lymph",          verified = TRUE),
    antigen0       = list(analyte = "antigen",                      units = NA_character_, specimen = "not applicable", verified = TRUE),
    antigen1       = list(analyte = "antigen",                      units = NA_character_, specimen = "not applicable", verified = TRUE),
    antigen2       = list(analyte = "antigen",                      units = NA_character_, specimen = "not applicable", verified = TRUE),
    antigen3       = list(analyte = "antigen",                      units = NA_character_, specimen = "not applicable", verified = TRUE),
    antigen4       = list(analyte = "antigen",                      units = NA_character_, specimen = "not applicable", verified = TRUE),
    antigen5       = list(analyte = "antigen",                      units = NA_character_, specimen = "not applicable", verified = TRUE),
    antigen        = list(analyte = "antigen",                      units = NA_character_, specimen = "not applicable", verified = TRUE),
    asc_spleen     = list(analyte = "antibody-secreting cells",     units = "cell",      specimen = "tissue",         verified = TRUE),
    asc_ln         = list(analyte = "antibody-secreting cells",     units = "cell",      specimen = "lymph",          verified = TRUE),
    asc_bone       = list(analyte = "antibody-secreting cells",     units = "cell",      specimen = "tissue",         verified = TRUE),
    asc_blood      = list(analyte = "antibody-secreting cells",     units = "cell",      specimen = "whole blood",    verified = TRUE),
    asc_peripheral = list(analyte = "antibody-secreting cells",     units = "cell",      specimen = "tissue",         verified = TRUE),
    igg_blood      = list(analyte = "antigen-specific IgG",         units = "pg",        specimen = "plasma",         verified = TRUE)
  )

  population <- list(
    species       = "mouse",
    n_subjects    = NA_integer_,
    n_studies     = 21L,
    disease_state = paste(
      "Healthy, non-genetically-modified laboratory mice. The homeostasis",
      "sub-model describes unimmunized steady-state B cell subpopulations;",
      "the activation layer describes a primary T cell-dependent response to",
      "a single antigen exposure.",
      sep = " "
    ),
    dose_range    = paste(
      "No drug. Immunogens across the calibration and validation studies were",
      "virus-like particles, LCMV, Plasmodium chabaudi, dengue virus 2, sheep",
      "red blood cells, keyhole limpet hemocyanin, NP-ova/CFA, chicken",
      "ovalbumin, NP18-chicken gamma-globulin, tetanus toxin heavy chain",
      "fragment and a 7-valent pneumococcal conjugate vaccine, given",
      "intraperitoneally, intravenously or intradermally",
      "(Supplementary Table S1).",
      sep = " "
    ),
    notes         = paste(
      "Aggregated across 21 published mouse studies identified by a systematic",
      "literature search (Supplementary Table S1); no individual-level data",
      "were used, so there are no per-animal demographics and n_subjects is",
      "not reported. Strains span C57BL/6 (and C57BL/6J, CBA/J, CBA/N) and",
      "BALB/c (and BALB/cJRj); age 5-24 weeks where reported; both sexes.",
      "Calibration targeted ASC dynamics in spleen, lymph nodes, bone marrow",
      "and blood (references 33-38 of the main text); external validation used",
      "independent plasma IgG kinetics from 5 studies, normalized to each",
      "study's maximum (RMSE 0.253).",
      "Two tissue-scaling conventions are baked into the calibration data and",
      "are reproduced here: bone marrow ASC counts are two-femur counts",
      "multiplied by 7.9 (femur = 12.7% of total murine marrow, Katz 1968),",
      "and lymph node ASC counts are single-node counts multiplied by 22, the",
      "number of murine lymph nodes.",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # Ugolkov 2026 Table 1 gives all 31 parameters. Time is in DAYS.
    #
    # TWO AMOUNT SCALES (see compartmentData above): the homeostasis B cell
    # states and their steady-state anchors are in 10^6 cells; the ASC states
    # and IgG are absolute cells / pg. Table 1 prints the anchors as
    # "1.97 * 10^6 cell" etc.; the deposited RxODE code writes them as 1.97
    # with the comment "#cell / 10^6". The value below is the code's scale.
    #
    # ESTIMATED vs FIXED: the Methods state that 19 parameters were fixed
    # from experimental data, 7 were derived analytically from steady-state
    # cell counts (those 7 are computed in model(), not listed here), and 5
    # were estimated. Only the 5 estimated parameters are log-transformed,
    # matching the paper's own statement that "all calibrated parameters were
    # log-transformed for fitting".
    # =====================================================================

    # --- Steady-state cell-count anchors (Table 1, "Algebraic steady-state
    #     parameters"). These are observed steady-state counts, and they are
    #     also the initial conditions of the homeostasis states.
    bl_imm_bone     <- fixed(1.97)  ; label("Immature B cells in bone marrow at steady state (10^6 cell)")       # Table 1 ImmBone_ss = 1.97 * 10^6 cell, calculated from ref 54; code line ImmBone_ss = 1.97
    bl_t1_bone      <- fixed(0.76)  ; label("Transitional T1 B cells in bone marrow at steady state (10^6 cell)")# Table 1 T1Bone_ss = 0.76 * 10^6 cell, calculated from ref 54
    bl_t1_blood     <- fixed(0.004) ; label("Transitional T1 B cells in blood at steady state (10^6 cell)")      # Table 1 T1Blood_ss = 0.004 * 10^6 cell, calculated from ref 55
    bl_t1_spleen    <- fixed(4.40)  ; label("Transitional T1 B cells in spleen at steady state (10^6 cell)")     # Table 1 T1Spleen_ss = 4.40 * 10^6 cell, calculated from ref 56
    bl_naive_spleen <- fixed(24.30) ; label("Naive B cells in spleen at steady state (10^6 cell)")               # Table 1 NaiveSpleen_ss = 24.30 * 10^6 cell, calculated from refs 4, 56
    bl_naive_blood  <- fixed(4.64)  ; label("Naive B cells in blood at steady state (10^6 cell)")                # Table 1 NaiveBlood_ss = 4.64 * 10^6 cell, calculated from ref 57
    bl_naive_ln     <- fixed(1.90)  ; label("Naive B cells in lymph nodes at steady state (10^6 cell)")          # Table 1 NaiveLN_ss = 1.90 * 10^6 cell, calculated from ref 58

    # --- Non-activated B cell migration rate constants that were NOT derived
    #     from steady state (Table 1, "Rate constant of migration"). The other
    #     four (k1, k3, k5, k7) are computed in model() from the anchors above.
    k2 <- fixed(1.01)  ; label("T1 B cell migration, spleen to blood (1/day)")        # Table 1 k2 = 1.01, calculated from ref 59
    k4 <- fixed(10.08) ; label("Naive B cell migration, spleen to blood (1/day)")     # Table 1 k4 = 10.08, calculated from ref 60
    k6 <- fixed(4.9)   ; label("Naive B cell migration, blood to lymph nodes (1/day)")# Table 1 k6 = 4.9, calculated from ref 60

    # --- Other homeostasis rate constants (Table 1, "Other rate constants").
    ksyn_imm <- fixed(20)   ; label("Immature B cell generation rate in bone marrow (10^6 cell/day)")     # Table 1 ksyn_imm = 20 * 10^6 cell/day, calculated from ref 61; code ksyn_imm = 20 "#(cell / 10^6)/day"
    kdeg_ln  <- fixed(0.01) ; label("Naive B cell elimination rate constant in lymph nodes (1/day)")
    # ^^ ERRATUM. Table 1 prints kdeg_ln = 0.05 (calculated from ref 62); the
    #    deposited RxODE model code sets kdeg_ln = 0.01. The code value is
    #    used because it is the value that reproduces the paper's own
    #    published outputs -- see the vignette Errata section, which shows
    #    that 0.01 (with the k7 derivation below) reproduces every Figure 2
    #    Cmax / Tmax within 8%, whereas the Table 1 pair collapses the naive
    #    B cell pool and misses the ASC Cmax values by roughly 90%.

    # --- ASC migration rate constants that were NOT estimated
    #     (Table 1, "Rate constant of ASC migration").
    k8  <- fixed(2.52)  ; label("ASC migration, spleen to blood (1/day)")            # Table 1 k8 = 2.52, calculated from ref 63
    k10 <- fixed(143.60); label("ASC migration, blood to lymph nodes (1/day)")       # Table 1 k10 = 143.60, assumed on ref 64
    k11 <- fixed(2.72)  ; label("ASC migration, lymph nodes to blood (1/day)")       # Table 1 k11 = 2.72, calculated from ref 63
    k14 <- fixed(0.02)  ; label("ASC migration, peripheral tissues to blood (1/day)")# Table 1 k14 = 0.02, assumed on refs 50, 51; Results: set equal to k12 because ASC persist in these compartments

    # --- IgG production and elimination (Table 1, "Other rate constants").
    ksin_igg <- fixed(22.00); label("IgG synthesis rate constant per ASC (pg/cell/day)")# Table 1 ksin_igg = 22.00, calculated from ref 65
    kdeg_igg <- fixed(0.12) ; label("IgG elimination rate constant (1/day)")            # Table 1 kdeg_igg = 0.12, calculated from ref 66

    # --- Bone marrow survival-niche saturation (Table 1, "Other parameters").
    khalf <- fixed(0.10) ; label("Half-saturation constant of ASC influx into bone marrow (cell)")  # Table 1 Khalf = 0.10; Methods: fixed at 0.1 to reproduce the rapid marrow accumulation of refs 33-38

    # --- Antigen transit chain. Not in Table 1; from the Results narrative
    #     and the deposited code.
    mtt <- fixed(18) ; label("Mean transit time of the antigen transit chain (day)")  # Results: "A mean transit time of 18 days was selected to match the observed average time to reach the maximum ASC counts (Tmax) (17.4 +/- 9.45 days)"; code MTT = 18

    # --- The 5 ESTIMATED parameters (Table 1; RSE column). All were fitted on
    #     the log scale, so they are the only log-transformed entries here.
    lk9      <- log(8.53)  ; label("ASC migration, blood to spleen (1/day)")                       # Table 1 k9 = 8.53, RSE 13.14%, estimated from refs 33-38
    lk12     <- log(0.02)  ; label("ASC migration, bone marrow to blood (1/day)")                  # Table 1 k12 = 0.02, RSE 23.23%, estimated from refs 33-38
    lk13     <- log(16.26) ; label("ASC migration, blood to peripheral tissues (1/day)")           # Table 1 k13 = 16.26, RSE 23.55%, estimated from refs 33-38
    lkmat_3  <- log(294.19); label("Effective ASC generation rate constant from naive B cells (1/day)") # Table 1 kmat_3 = 294.19, RSE 11.39%, estimated from refs 33-38
    lvmax    <- log(490.60); label("Maximum ASC influx into the bone marrow survival niche (cell/day)") # Table 1 Vmax = 490.60, RSE 19.51%, estimated from refs 33-38; code Vmax = 490.6. The inter-study heterogeneity section quotes 492.8 as "the typical calibrated value" and a 26.66-836 range across the six datasets refitted individually.
  })

  model({
    # ===================================================================
    # 1. Estimated parameters back-transformed from the log scale.
    # ===================================================================
    k9     <- exp(lk9)
    k12    <- exp(lk12)
    k13    <- exp(lk13)
    kmat_3 <- exp(lkmat_3)
    vmax   <- exp(lvmax)

    # ===================================================================
    # 2. The seven parameters derived analytically from steady state.
    #    Results: "Seven parameters (k1, k3, k5, k7, kmat_1, kmat_2,
    #    kdeg_spl) were derived analytically using steady-state
    #    calculations. For each relevant equation, the time derivative was
    #    set to zero, and the resulting algebraic expression was solved for
    #    the unknown parameter."
    #    Reproduced here exactly as Supplementary Equations 1-7 and the
    #    deposited RxODE code write them.
    # ===================================================================
    kmat_1 <- ksyn_imm / bl_imm_bone                                             # Supp Eq 1  -> 10.15 (Table 1 kmat_1 = 10.15)
    k1     <- (kmat_1 * bl_imm_bone) / bl_t1_bone                                # Supp Eq 2  -> 26.32 (Table 1 k1 = 26.32)
    k3     <- (k1 * bl_t1_bone + k2 * bl_t1_spleen) / bl_t1_blood                # Supp Eq 3  -> 6111  (Table 1 k3 = 6111)
    kmat_2 <- (k3 * bl_t1_blood - k2 * bl_t1_spleen) / bl_t1_spleen              # Supp Eq 4  -> 4.55  (Table 1 kmat_2 = 4.55)

    # Supp Eq 5. Two notes on this one line, both detailed in the vignette
    # Errata:
    #  (a) The degradation term uses the SPLEEN anchor, not the lymph-node
    #      anchor that a strict steady-state solve of Supp Eq 14 would give.
    #      That is what BOTH Supp Eq 5 and the deposited code write, so it is
    #      reproduced rather than silently corrected. With kdeg_ln = 0.01 it
    #      leaves the lymph-node balance off by ~1%, and the naive lymph-node
    #      pool drifts to 1.02x baseline over 350 days.
    #  (b) It does NOT reproduce the Table 1 value of k7 = 62.62; it gives
    #      11.84. Table 1's 62.62 is recoverable only by substituting the
    #      spleen anchor into the INFLUX term as well
    #      ((4.9 * 24.30 - 0.05 * 1.90) / 1.90 = 62.62), which is not what any
    #      printed equation says, and which breaks the model's published
    #      behaviour.
    k7 <- (k6 * bl_naive_blood - kdeg_ln * bl_naive_spleen) / bl_naive_ln        # Supp Eq 5 / code line k7 = (k6 * NaiveBlood_ss - kdeg_ln * NaiveSpleen)/(NaiveLN_ss)

    k5 <- (k7 * bl_naive_ln - k6 * bl_naive_blood + k4 * bl_naive_spleen) /
      bl_naive_blood                                                             # Supp Eq 6
    kdeg_spl <- (kmat_2 * bl_t1_spleen - k4 * bl_naive_spleen + k5 * bl_naive_blood) /
      bl_naive_spleen                                                            # Supp Eq 7

    # Antigen transit rate constant. The chain has six transfers
    # (antigen0 -> antigen1 -> ... -> antigen5 -> antigen), and the paper's
    # code sets klag = 7 / MTT. The canonical bare name for this quantity is
    # ktr; the source calls it klag.
    ktr <- 7 / mtt                                                               # code line klag = 7 / MTT

    # ===================================================================
    # 3. Homeostasis sub-model. Supplementary Equations 8-14.
    #    States in 10^6 cells.
    # ===================================================================
    d/dt(imm_bone)     <- ksyn_imm - kmat_1 * imm_bone                           # Supp Eq 8
    d/dt(t1_bone)      <- kmat_1 * imm_bone - k1 * t1_bone                       # Supp Eq 9
    d/dt(t1_blood)     <- k1 * t1_bone + k2 * t1_spleen - k3 * t1_blood          # Supp Eq 10
    d/dt(t1_spleen)    <- k3 * t1_blood - k2 * t1_spleen - kmat_2 * t1_spleen    # Supp Eq 11
    d/dt(naive_blood)  <- k4 * naive_spleen - k5 * naive_blood +
      k7 * naive_ln - k6 * naive_blood                                           # Supp Eq 12
    d/dt(naive_spleen) <- kmat_2 * t1_spleen - kdeg_spl * naive_spleen -
      k4 * naive_spleen + k5 * naive_blood                                       # Supp Eq 13
    d/dt(naive_ln)     <- k6 * naive_blood - k7 * naive_ln - kdeg_ln * naive_ln  # Supp Eq 14

    # ===================================================================
    # 4. Antigen forcing function and transit chain.
    #    Supplementary Equations 15-21.
    #
    #    Eq 15 is the derivative of a pre-defined empirical antigen-kinetics
    #    curve fitted to the viral-load time course of Castillo-Mendez 2007
    #    (Supplementary Figure 1). The three constants 393.12, 66.98 and
    #    189.1 are printed as bare literals in Eq 15 and in the deposited
    #    code; they are unnamed shape constants of that curve, so they are
    #    reproduced here as literals rather than promoted to named
    #    parameters.
    #
    #    tclamp guards the 1/t, 1/t^2 and log(t) singularity at t = 0 so the
    #    system can be integrated from time zero. The guard is numerically
    #    inert: exp(189.1 - 393.12/t - 66.98*log(t)) underflows to exactly 0
    #    for t <= 1e-3, which is also the true limit of Eq 15 as t -> 0+.
    # ===================================================================
    tclamp <- max(t, 1e-3)
    d/dt(antigen0) <- (393.12 / tclamp^2 - 66.98 / tclamp) *
      exp(189.1 - 393.12 / tclamp - 66.98 * log(tclamp))                         # Supp Eq 15
    d/dt(antigen1) <- ktr * antigen0 - ktr * antigen1                            # Supp Eq 16
    d/dt(antigen2) <- ktr * antigen1 - ktr * antigen2                            # Supp Eq 17
    d/dt(antigen3) <- ktr * antigen2 - ktr * antigen3                            # Supp Eq 18
    d/dt(antigen4) <- ktr * antigen3 - ktr * antigen4                            # Supp Eq 19
    d/dt(antigen5) <- ktr * antigen4 - ktr * antigen5                            # Supp Eq 20
    d/dt(antigen)  <- ktr * antigen5 - ktr * antigen                             # Supp Eq 21

    # ===================================================================
    # 5. ASC generation, distribution and the bone marrow survival niche.
    #    Supplementary Equations 22-26.
    #
    #    Eqs 22-23 are where the two amount scales meet. The paper writes the
    #    activated precursor pool as NaiveSpleen / 10^6 and NaiveLN * 22 /
    #    10^6, with the naive counts in absolute cells and the 10^6 being the
    #    physiological antigen-specific precursor frequency of 1 in 10^6
    #    (Results, citing Abbott 2018). Because naive_spleen and naive_ln are
    #    already carried in units of 10^6 cells, the division by 10^6 is
    #    already applied and the terms below correctly read as absolute
    #    numbers of activated precursors -- which is exactly how the
    #    deposited code writes them. The 22 is the number of murine lymph
    #    nodes, converting a single-node count to the whole-body pool.
    #
    #    Naive B cells are deliberately NOT depleted by activation: at a
    #    precursor frequency of 1 in 10^6 the loss is negligible (Methods).
    #    ASC elimination is likewise absent: ASC levels are stable to 300
    #    days post-immunization (Methods, citing ref 36).
    # ===================================================================
    d/dt(asc_spleen) <- kmat_3 * antigen * naive_spleen -
      k8 * asc_spleen + k9 * asc_blood                                           # Supp Eq 22
    d/dt(asc_ln)     <- kmat_3 * antigen * (naive_ln * 22) -
      k11 * asc_ln + k10 * asc_blood                                             # Supp Eq 23
    d/dt(asc_bone)   <- (vmax * asc_blood) / (khalf + asc_blood) -
      k12 * asc_bone                                                             # Supp Eq 24
    d/dt(asc_blood)  <- k11 * asc_ln + k8 * asc_spleen + k12 * asc_bone +
      k14 * asc_peripheral -
      (k9 + k10 + k13 + vmax / (khalf + asc_blood)) * asc_blood                  # Supp Eq 25
    d/dt(asc_peripheral) <- k13 * asc_blood - k14 * asc_peripheral               # Supp Eq 26

    # ===================================================================
    # 6. Antigen-specific IgG. Supplementary Equations 27-28.
    #    ASC in peripheral non-lymphoid tissues are excluded from the pool
    #    driving systemic IgG: they are tissue-adapted and secrete locally
    #    (Methods, citing refs 39, 40).
    # ===================================================================
    ascTotal <- asc_spleen + asc_bone + asc_blood + asc_ln                       # Supp Eq 27
    d/dt(igg_blood) <- ksin_igg * ascTotal - kdeg_igg * igg_blood                # Supp Eq 28

    # ===================================================================
    # 7. Initial conditions (deposited RxODE model code). The homeostasis
    #    states start at their observed steady-state counts; every antigen,
    #    ASC and IgG state starts at zero (rxode2's default, written out
    #    here only for the homeostasis states as the source does).
    # ===================================================================
    imm_bone(0)     <- bl_imm_bone
    t1_bone(0)      <- bl_t1_bone
    t1_blood(0)     <- bl_t1_blood
    t1_spleen(0)    <- bl_t1_spleen
    naive_blood(0)  <- bl_naive_blood
    naive_spleen(0) <- bl_naive_spleen
    naive_ln(0)     <- bl_naive_ln

    # No residual-error model. Ugolkov 2026 minimized a -2LL whose per-variable
    # sigma was computed from the residuals themselves, and tabulates no
    # residual-error magnitude for any output; no IIV is reported either (the
    # reported 95% CIs are parameter-uncertainty bands from the
    # variance-covariance matrix, explicitly NOT prediction intervals).
  })
}
