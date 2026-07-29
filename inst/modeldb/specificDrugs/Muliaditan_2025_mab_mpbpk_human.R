Muliaditan_2025_mab_mpbpk_human <- function() {
  description <- paste(
    "Human-scaled. Translational minimal physiologically based",
    "pharmacokinetic (mPBPK) model for transferrin-receptor (TfR)",
    "mediated brain delivery of monoclonal antibodies, projected",
    "forward from the cynomolgus monkey fit",
    "Muliaditan_2025_mab_mpbpk_nhp by replacing the Bloomingdale",
    "2017 NHP physiology with the human physiology (Muliaditan 2025",
    "Supplementary Table S1, human column), allometrically scaling",
    "the bsAb-TfR internalization rate kint by (70/6.2)^(-0.25) =",
    "0.546 (paper Methods: standard rate-constant exponent -0.25),",
    "and recalibrating the luminal BCSFB unbound TfR baseline",
    "uTFR0_BCSFB to be 3-fold higher than the NHP estimate",
    "(0.256 -> 0.768 nM) per the paper Results. The other TfR-",
    "related parameters (TfRpt, uTFR0_BBB, FACQ_BECF, TfRtotn,",
    "ktrans, kdeg_uTfR_BBB, kdeg_uTfR_BCSFB, FACBR) are assumed",
    "identical to the NHP estimates per paper Methods. Per-compound",
    "TfR binding parameters (kon_T, koff_T) MUST be set per simulated",
    "antibody from biophysical measurements (paper Table S2); the",
    "default ini() encodes kon_T = 0 (non-TfR control IgG). For",
    "trontinemab in human, the paper reports KD,TfR = 131 nM",
    "(versus 249 nM in NHP) with kon_T = 1.0548 nM^-1 h^-1 and",
    "koff_T = 138.24 h^-1 (Table S2, Grimm 2023 column for human).",
    "Clinical validation in the paper was against single ascending",
    "doses 0.1-7.2 mg/kg IV trontinemab in healthy human subjects",
    "(NCT04023994; Grimm 2023)."
  )
  reference <- paste(
    "Muliaditan M, van Steeg TJ, Avery LB, Sun W, Hammond TR, Hijdra D,",
    "Choi SL, Pillai N, Leksa NC, Mavroudis PD.",
    "Translational minimal physiologically based pharmacokinetic model",
    "for transferrin receptor-mediated brain delivery of antibodies.",
    "mAbs. 2025;17(1):2515414. doi:10.1080/19420862.2025.2515414.",
    "Human physiology from Bloomingdale 2017 (Muliaditan 2025 Supp",
    "Table S1, human column); TfR parameters carried from the NHP fit",
    "with kint allometrically scaled and uTFR0_BCSFB recalibrated.",
    "See modellib('Muliaditan_2025_mab_mpbpk_nhp') for the primary",
    "NHP fit.",
    sep = " "
  )
  vignette <- "Muliaditan_2025_mab_mpbpk"
  units <- list(
    time          = "hour",
    dosing        = "nmol (convert mg via dose_nmol = dose_mg / mw_kda * 1000)",
    concentration = "nmol/L"
  )

  covariateData <- list(
    MIX_FAST_ELIM = list(
      description        = paste(
        "Mixture-model class indicator for fast vs slow bsAb-TfR complex",
        "internalization rate (kint). 1 = POP1 fast eliminator",
        "(allometrically scaled human kint = 0.0179 h^-1);",
        "0 = POP2 slow eliminator (allometrically scaled human kint =",
        "0.0125 * (70/6.2)^(-0.25) = 0.00683 h^-1)."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "Population probability of MIX_FAST_ELIM = 1 (POP1, fast)",
        "is 0.437 (Muliaditan 2025 Table 2, Fraction POP1, RSE 30%).",
        "For typical-value simulation set MIX_FAST_ELIM = 0 (slow,",
        "dominant 56.3% subpopulation). The mixture-fraction was",
        "estimated in NHP and assumed to translate to humans by the",
        "paper. Allometric scaling applied to both POP1 and POP2 kint",
        "values uses the (70/6.2)^(-0.25) = 0.546 factor (paper",
        "Methods, exponent -0.25)."
      ),
      source_name        = "MIXTURE"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = NA_integer_,
    n_studies       = 1L,
    weight_range    = "Bloomingdale 2017 70-kg human reference subject for the physiology",
    disease_state   = paste(
      "Healthy adult subjects (trontinemab Phase 1 single ascending",
      "dose study NCT04023994; Grimm 2023). Plasma exposure metrics",
      "(AUC0-168h, AUC0-inf, Cmax) and CSF concentrations at study",
      "days 3 and 5 were used for clinical validation in the paper",
      "(Muliaditan 2025 Methods / Figure 5)."
    ),
    dose_range      = "0.1, 0.4, 1.2, 3.6, 7.2 mg/kg single IV infusion (trontinemab single ascending dose; Grimm 2023)",
    notes           = paste(
      "Human PK projection only - no human dataset was fit. The",
      "translated mPBPK parameters are kept identical to the NHP",
      "fit except for (a) Bloomingdale 2017 70-kg human physiology",
      "from Supp Table S1, (b) allometric scaling of kint with",
      "exponent -0.25, and (c) 3x higher uTFR0_BCSFB recalibrated",
      "to improve CSF prediction accuracy (paper Results)."
    )
  )

  ini({
    # ============================================================
    # Population structural parameters: TfRpt, uTFR0_BBB, FACQ_BECF
    # are carried from the NHP fit (assumed identical in humans per
    # Muliaditan 2025 Methods). uTFR0_BCSFB is the NHP value (0.256)
    # multiplied by 3 = 0.768 nM per paper Results recalibration.
    # ============================================================
    ltfrpt          <- log(1672);    label("Total plasma TfR concentration TfRpt (nM)")                              # Muliaditan 2025 Table 2 (NHP, assumed same in human)
    lutfr0_bbb      <- log(175);     label("Baseline luminal unbound TfR on BBB, uTFR0_BBB (nM)")                    # Muliaditan 2025 Table 2 (NHP, assumed same in human)
    lutfr0_bcsfb    <- log(0.768);   label("Baseline luminal unbound TfR on BCSFB, uTFR0_BCSFB (nM, 3x NHP value)")  # Muliaditan 2025 Results (0.256 * 3 = 0.768; recalibrated for CSF prediction)
    lkint_pop1      <- log(0.0179);  label("Allometrically scaled kint POP1, fast subpop (h^-1)")                    # Muliaditan 2025 Methods (0.0329 * (70/6.2)^(-0.25) = 0.0179)
    lkint_pop2      <- log(0.0125 * (6.2 / 70)^(0.25)); label("Allometrically scaled kint POP2, slow subpop (h^-1)") # Muliaditan 2025 Methods (0.0125 * 0.546 = 0.00683; same allometric exponent)
    lfacq_becf      <- log(0.00814); label("Correction factor on Q_BECF for ISF -> CSF distribution (unitless)")     # Muliaditan 2025 Table 2 (NHP, assumed same in human)

    # ============================================================
    # Fixed structural parameters (NHP and human identical).
    # ============================================================
    ltfrtotn         <- fixed(log(559));  label("Total neuronal TfR in brain ISF TfRtotn (nM)")        # Muliaditan 2025 Table 2 / Chang 2022 (559 nM, FIXED)
    lktrans          <- fixed(log(6));    label("Transcytosis rate ktrans across BBB and BCSFB (h^-1)") # Muliaditan 2025 Table 2 / Chang 2022 (6 h^-1, FIXED)
    lkdeg_utfr_bbb   <- fixed(log(20));   label("Degradation rate of luminal unbound TfR on BBB (h^-1)") # Muliaditan 2025 Table 2 (20 h^-1, FIXED)
    lkdeg_utfr_bcsfb <- fixed(log(20));   label("Degradation rate of luminal unbound TfR on BCSFB (h^-1)") # Muliaditan 2025 Table 2 (20 h^-1, FIXED)
    lfacbr           <- fixed(log(0.05)); label("Brain-prediction correction factor FAC_BPRED (unitless)") # Muliaditan 2025 Table 2 (0.05, FIXED)

    # ============================================================
    # Per-compound TfR binding parameters (set per simulated antibody;
    # Muliaditan 2025 Supplementary Table S2). Default = non-TfR
    # control IgG (kon_T = 0). To simulate trontinemab in HUMAN:
    #   kon_T  = 1.0548   nM^-1 h^-1, koff_T = 138.24 h^-1
    #   (Muliaditan 2025 Table S2 trontinemab human row;
    #    KD,TfR = 131 nM versus 249 nM in NHP).
    # ============================================================
    lkon_t   <- fixed(log(1e-12)); label("TfR association rate kon_T (nM^-1 h^-1) - per-compound input; 1e-12 placeholder for non-TfR control") # Muliaditan 2025 Table S2; default = non-TfR control
    lkoff_t  <- fixed(log(1));     label("TfR dissociation rate koff_T (h^-1) - per-compound input")                                              # Muliaditan 2025 Table S2; default placeholder
    mw_kda   <- fixed(150);        label("Antibody molecular weight (kDa) for mg -> nmol dose conversion") # default 150 kDa IgG (trontinemab is 194 kDa per Table S2)

    # ============================================================
    # Residual error - additive on the log-transformed concentration.
    # Carried from the NHP fit (paper does not re-estimate for human).
    # ============================================================
    expSd        <- sqrt(0.256); label("Log-scale residual SD for plasma total drug (unitless)")     # Muliaditan 2025 Table 2 (NHP estimate)
    expSd_Ccsf   <- sqrt(0.893); label("Log-scale residual SD for CSF unbound drug (unitless)")      # Muliaditan 2025 Table 2 (NHP estimate)
    expSd_Cbrain <- sqrt(0.426); label("Log-scale residual SD for whole-brain homogenate (unitless)")# Muliaditan 2025 Table 2 (NHP estimate)
  })

  model({
    # Resolve parameters from log space.
    tfrpt        <- exp(ltfrpt)
    utfr0_bbb    <- exp(lutfr0_bbb)
    utfr0_bcsfb  <- exp(lutfr0_bcsfb)
    kint_pop1    <- exp(lkint_pop1)
    kint_pop2    <- exp(lkint_pop2)
    facq_becf    <- exp(lfacq_becf)
    tfrtotn      <- exp(ltfrtotn)
    ktrans       <- exp(lktrans)
    kdeg_utfr_bbb   <- exp(lkdeg_utfr_bbb)
    kdeg_utfr_bcsfb <- exp(lkdeg_utfr_bcsfb)
    facbr        <- exp(lfacbr)
    kon_t        <- exp(lkon_t)
    koff_t       <- exp(lkoff_t)

    kint <- kint_pop2 + (kint_pop1 - kint_pop2) * MIX_FAST_ELIM
    krec_utfr       <- kint
    ksyn_utfr_bbb   <- kdeg_utfr_bbb   * utfr0_bbb
    ksyn_utfr_bcsfb <- kdeg_utfr_bcsfb * utfr0_bcsfb

    # ============================================================
    # Bloomingdale 2017 human physiology (Muliaditan 2025 Supp Table
    # S1, human column). QT uses the volume-weighted-average value
    # from the supplement footnote (a) rather than the original
    # Bloomingdale QT.
    # ============================================================
    vp           <- 3.13            # plasma                          # Muliaditan 2025 Supp Table S1 human
    v_tissue_v   <- 1.68            # tissue vascular VTv             # Muliaditan 2025 Supp Table S1 human
    v_tissue_e   <- 0.335           # tissue endosomal VTe            # Muliaditan 2025 Supp Table S1 human
    v_tissue_i   <- 11.1            # tissue interstitial VTi         # Muliaditan 2025 Supp Table S1 human
    v_brain_v    <- 0.0319          # brain vascular VBv              # Muliaditan 2025 Supp Table S1 human
    v_bbb_e      <- 0.00659         # BBB endosomal V_BE_BBB          # Muliaditan 2025 Supp Table S1 human
    v_bcsfb_e    <- 6.59e-4         # BCSFB endosomal V_BE_BCSFB      # Muliaditan 2025 Supp Table S1 human
    v_brain_i    <- 0.261           # brain interstitial VBi          # Muliaditan 2025 Supp Table S1 human
    v_csf        <- 0.143           # CSF VCSF                        # Muliaditan 2025 Supp Table S1 human
    v_lymph      <- 0.274           # lymph VL                        # Muliaditan 2025 Supp Table S1 human

    qt           <- 22.1            # tissue plasma flow              # Muliaditan 2025 Supp Table S1 human (volume-weighted average)
    qb           <- 21.5            # brain plasma flow               # Muliaditan 2025 Supp Table S1 human
    lt           <- 0.321           # tissue lymphatic flow           # Muliaditan 2025 Supp Table S1 human
    lb           <- 0.0345          # brain lymphatic flow            # Muliaditan 2025 Supp Table S1 human
    qb_ecf       <- 0.0105          # brain ECF flow                  # Muliaditan 2025 Supp Table S1 human
    qb_csf       <- 0.024           # brain CSF flow                  # Muliaditan 2025 Supp Table S1 human

    k_clupt      <- 0.55            # tissue endosomal uptake rate    # Muliaditan 2025 Supp Table S1 human
    k_clupb      <- 0.0195          # brain endosomal uptake rate     # Muliaditan 2025 Supp Table S1 human
    clup_t       <- 0.184           # tissue uptake clearance CLUPT   # Muliaditan 2025 Supp Table S1 human
    clup_bbb     <- 1.29e-4         # BBB uptake clearance CLUPBBB    # Muliaditan 2025 Supp Table S1 human
    clup_bcsfb   <- 1.29e-5         # BCSFB uptake clearance CLUPBCSFB# Muliaditan 2025 Supp Table S1 human

    sig_bbb      <- 1               # BBB reflection coefficient      # Muliaditan 2025 Supp Table S1 human
    sig_bcsfb    <- 0.9974          # BCSFB reflection coefficient    # Muliaditan 2025 Supp Table S1 human
    sig_tv       <- 0.9233          # tissue vascular reflection      # Muliaditan 2025 Supp Table S1 human
    sig_tl       <- 0.2             # tissue lymph reflection         # Muliaditan 2025 Supp Table S1 human
    sig_isf      <- 0.2             # brain ISF reflection            # Muliaditan 2025 Supp Table S1 human
    sig_csf      <- 0.2             # CSF reflection                  # Muliaditan 2025 Supp Table S1 human

    fr_tissue    <- 0.715           # tissue FcRn fractional recycling # Muliaditan 2025 Supp Table S1 human
    fr_brain     <- 0.715           # brain FcRn fractional recycling  # Muliaditan 2025 Supp Table S1 human

    fcrn_baseline_M <- 4.98e-5      # mol/L                           # Muliaditan 2025 Supp Table S1 human
    fcrn_baseline   <- fcrn_baseline_M * 1e9

    kon_fcrn_M  <- 5.59e8           # M^-1 h^-1                       # Muliaditan 2025 Supp Table S1 human (differs from NHP 7.92e8)
    kon_fcrn    <- kon_fcrn_M * 1e-9
    koff_fcrn   <- 23.9             # h^-1                            # Muliaditan 2025 Supp Table S1 human (differs from NHP 46.8)
    kdeg_endo   <- 26.6             # h^-1                            # Muliaditan 2025 Supp Table S1 human

    # Dose unit conversion (mg -> nmol).
    f(central) <- 1000 / mw_kda

    # Concentrations from amounts.
    c_plasma          <- central          / vp
    c_tissue_vasc     <- tissue_vasc      / v_tissue_v
    c_tissue_endo_u   <- tissue_endo_u    / v_tissue_e
    c_tissue_endo_b   <- tissue_endo_b    / v_tissue_e
    c_tissue_isf      <- tissue_isf       / v_tissue_i
    c_brain_vascular      <- brain_vascular       / v_brain_v
    c_bbb_endo_u      <- bbb_endo_u       / v_bbb_e
    c_bbb_endo_b      <- bbb_endo_b       / v_bbb_e
    c_brain_isf       <- brain_isf        / v_brain_i
    c_bcsfb_endo_u    <- bcsfb_endo_u     / v_bcsfb_e
    c_bcsfb_endo_b    <- bcsfb_endo_b     / v_bcsfb_e
    c_csf             <- csf              / v_csf
    c_lymph           <- lymph            / v_lymph
    c_tissue_fcrn     <- tissue_fcrn      / v_tissue_e
    c_bbb_fcrn        <- bbb_fcrn         / v_bbb_e
    c_bcsfb_fcrn      <- bcsfb_fcrn       / v_bcsfb_e
    c_complex_plasma  <- complex_plasma   / vp
    c_delta_utfr_bbb  <- delta_utfr_bbb   / v_brain_v
    c_complex_bbb_lum <- complex_bbb_lum  / v_brain_v
    c_complex_bbb_abl <- complex_bbb_abl  / v_brain_i
    c_utfr_bbb_abl    <- utfr_bbb_abl     / v_brain_i
    c_delta_utfr_bcsfb<- delta_utfr_bcsfb / v_brain_v
    c_complex_bcsfb_lum <- complex_bcsfb_lum / v_brain_v
    c_complex_bcsfb_abl <- complex_bcsfb_abl / v_csf
    c_utfr_bcsfb_abl    <- utfr_bcsfb_abl    / v_csf
    c_complex_neuron    <- complex_neuron    / v_brain_i

    utfr_bbb_lum    <- utfr0_bbb   + c_delta_utfr_bbb
    utfr_bcsfb_lum  <- utfr0_bcsfb + c_delta_utfr_bcsfb
    utfr_plasma     <- tfrpt   - c_complex_plasma
    utfr_neuron     <- tfrtotn - c_complex_neuron

    tfrb_plasma_rate    <- -kon_t * c_plasma     * utfr_plasma  + koff_t * c_complex_plasma
    tfrb_bbb_lum_rate   <- -kon_t * c_brain_vascular * utfr_bbb_lum    + koff_t * c_complex_bbb_lum
    tfrb_bcsfb_lum_rate <- -kon_t * c_brain_vascular * utfr_bcsfb_lum  + koff_t * c_complex_bcsfb_lum
    tfrb_bbb_abl_rate   <- -kon_t * c_brain_isf  * c_utfr_bbb_abl  + koff_t * c_complex_bbb_abl
    tfrb_neuron_rate    <- -kon_t * c_brain_isf  * utfr_neuron     + koff_t * c_complex_neuron
    tfrb_bcsfb_abl_rate <- -kon_t * c_csf        * c_utfr_bcsfb_abl+ koff_t * c_complex_bcsfb_abl

    # ODE system (same structural form as NHP file).
    d/dt(central) <- (qt - lt) * c_tissue_vasc + (qb - lb) * c_brain_vascular +
                     (lt + lb) * c_lymph -
                     qt * c_plasma - qb * c_plasma +
                     tfrb_plasma_rate * vp

    d/dt(tissue_vasc) <- qt * c_plasma - (qt - lt) * c_tissue_vasc -
                         (1 - sig_tv) * lt * c_tissue_vasc -
                         clup_t * c_tissue_vasc +
                         clup_t * fr_tissue * c_tissue_endo_b

    d/dt(tissue_endo_u) <- clup_t * (c_tissue_vasc + c_tissue_isf) -
                           (kon_fcrn * c_tissue_endo_u * c_tissue_fcrn -
                            koff_fcrn * c_tissue_endo_b) * v_tissue_e -
                           kdeg_endo * c_tissue_endo_u * v_tissue_e

    d/dt(tissue_endo_b) <- (kon_fcrn * c_tissue_endo_u * c_tissue_fcrn -
                            koff_fcrn * c_tissue_endo_b) * v_tissue_e -
                           clup_t * c_tissue_endo_b

    d/dt(tissue_isf) <- (1 - sig_tv) * lt * c_tissue_vasc -
                        (1 - sig_tl) * lt * c_tissue_isf +
                        clup_t * (1 - fr_tissue) * c_tissue_endo_b -
                        clup_t * c_tissue_isf

    d/dt(brain_vascular) <- qb * c_plasma - (qb - lb) * c_brain_vascular -
                        (1 - sig_bbb)   * qb_ecf * c_brain_vascular -
                        (1 - sig_bcsfb) * qb_csf * c_brain_vascular -
                        k_clupb * c_brain_vascular * v_brain_v +
                        clup_bbb   * fr_brain * c_bbb_endo_b +
                        clup_bcsfb * fr_brain * c_bcsfb_endo_b +
                        (tfrb_bbb_lum_rate + tfrb_bcsfb_lum_rate) * v_brain_v

    d/dt(bbb_endo_u) <- clup_bbb * (c_brain_vascular + c_brain_isf) -
                        (kon_fcrn * c_bbb_endo_u * c_bbb_fcrn -
                         koff_fcrn * c_bbb_endo_b) * v_bbb_e -
                        kdeg_endo * c_bbb_endo_u * v_bbb_e

    d/dt(bbb_endo_b) <- (kon_fcrn * c_bbb_endo_u * c_bbb_fcrn -
                         koff_fcrn * c_bbb_endo_b) * v_bbb_e -
                        clup_bbb * c_bbb_endo_b

    d/dt(brain_isf) <- (1 - sig_bbb)  * qb_ecf * c_brain_vascular -
                       (1 - sig_isf) * qb_ecf * c_brain_isf +
                       clup_bbb * (1 - fr_brain) * c_bbb_endo_b -
                       clup_bbb * c_brain_isf -
                       qb_ecf * facq_becf * c_brain_isf +
                       qb_ecf * c_csf +
                       (tfrb_bbb_abl_rate + tfrb_neuron_rate) * v_brain_i

    d/dt(bcsfb_endo_u) <- clup_bcsfb * (c_brain_vascular + c_csf) -
                          (kon_fcrn * c_bcsfb_endo_u * c_bcsfb_fcrn -
                           koff_fcrn * c_bcsfb_endo_b) * v_bcsfb_e -
                          kdeg_endo * c_bcsfb_endo_u * v_bcsfb_e

    d/dt(bcsfb_endo_b) <- (kon_fcrn * c_bcsfb_endo_u * c_bcsfb_fcrn -
                           koff_fcrn * c_bcsfb_endo_b) * v_bcsfb_e -
                          clup_bcsfb * c_bcsfb_endo_b

    d/dt(csf) <- (1 - sig_bcsfb) * qb_csf * c_brain_vascular -
                 clup_bcsfb * c_csf +
                 clup_bcsfb * (1 - fr_brain) * c_bcsfb_endo_b +
                 qb_ecf * facq_becf * c_brain_isf -
                 (1 - sig_csf) * qb_csf * c_csf -
                 qb_ecf * c_csf +
                 tfrb_bcsfb_abl_rate * v_csf

    d/dt(lymph) <- (1 - sig_tl)  * lt     * c_tissue_isf +
                   (1 - sig_csf) * qb_csf * c_csf +
                   (1 - sig_isf) * qb_ecf * c_brain_isf -
                   (lt + lb) * c_lymph

    d/dt(tissue_fcrn) <- -(kon_fcrn * c_tissue_endo_u * c_tissue_fcrn -
                           koff_fcrn * c_tissue_endo_b) * v_tissue_e +
                          clup_t * c_tissue_endo_b

    d/dt(bbb_fcrn) <- -(kon_fcrn * c_bbb_endo_u * c_bbb_fcrn -
                        koff_fcrn * c_bbb_endo_b) * v_bbb_e +
                       clup_bbb * c_bbb_endo_b

    d/dt(bcsfb_fcrn) <- -(kon_fcrn * c_bcsfb_endo_u * c_bcsfb_fcrn -
                          koff_fcrn * c_bcsfb_endo_b) * v_bcsfb_e +
                         clup_bcsfb * c_bcsfb_endo_b

    d/dt(complex_plasma) <- -tfrb_plasma_rate * vp - kint * complex_plasma

    d/dt(delta_utfr_bbb) <- (ksyn_utfr_bbb -
                             kdeg_utfr_bbb * utfr_bbb_lum +
                             tfrb_bbb_lum_rate +
                             krec_utfr * (v_brain_i / v_brain_v) * c_utfr_bbb_abl) * v_brain_v

    d/dt(complex_bbb_lum) <- (-tfrb_bbb_lum_rate -
                              ktrans * c_complex_bbb_lum -
                              kint   * c_complex_bbb_lum) * v_brain_v

    d/dt(complex_bbb_abl) <- ( ktrans * (v_brain_v / v_brain_i) * c_complex_bbb_lum -
                              tfrb_bbb_abl_rate -
                              kint   * c_complex_bbb_abl) * v_brain_i

    d/dt(utfr_bbb_abl) <- ( tfrb_bbb_abl_rate -
                            krec_utfr * c_utfr_bbb_abl) * v_brain_i

    d/dt(delta_utfr_bcsfb) <- (ksyn_utfr_bcsfb -
                               kdeg_utfr_bcsfb * utfr_bcsfb_lum +
                               tfrb_bcsfb_lum_rate +
                               krec_utfr * (v_csf / v_brain_v) * c_utfr_bcsfb_abl) * v_brain_v

    d/dt(complex_bcsfb_lum) <- (-tfrb_bcsfb_lum_rate -
                                ktrans * c_complex_bcsfb_lum -
                                kint   * c_complex_bcsfb_lum) * v_brain_v

    d/dt(complex_bcsfb_abl) <- ( ktrans * (v_brain_v / v_csf) * c_complex_bcsfb_lum -
                                tfrb_bcsfb_abl_rate -
                                kint   * c_complex_bcsfb_abl) * v_csf

    d/dt(utfr_bcsfb_abl) <- ( tfrb_bcsfb_abl_rate -
                              krec_utfr * c_utfr_bcsfb_abl) * v_csf

    d/dt(complex_neuron) <- (-tfrb_neuron_rate - kint * c_complex_neuron) * v_brain_i

    # Initial conditions: free FcRn pools start at the baseline level.
    tissue_fcrn(0) <- fcrn_baseline * v_tissue_e
    bbb_fcrn(0)    <- fcrn_baseline * v_bbb_e
    bcsfb_fcrn(0)  <- fcrn_baseline * v_bcsfb_e

    # Observations.
    Cc      <- c_plasma + c_complex_plasma             # plasma total drug
    Ccsf    <- c_csf                                   # CSF unbound

    homo_num <- (c_bbb_endo_u   + c_bbb_endo_b)   * v_bbb_e   +
                 c_brain_isf                      * v_brain_i +
                (c_bcsfb_endo_u + c_bcsfb_endo_b) * v_bcsfb_e
    homo_den <- v_bbb_e + v_brain_i + v_bcsfb_e
    Cbrain   <- facbr * homo_num / homo_den

    Cc     ~ lnorm(expSd)
    Ccsf   ~ lnorm(expSd_Ccsf)
    Cbrain ~ lnorm(expSd_Cbrain)
  })
}
