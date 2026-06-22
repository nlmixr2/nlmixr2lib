Muliaditan_2025_mab_mpbpk_nhp <- function() {
  description <- paste(
    "Preclinical (cynomolgus monkey). Translational minimal",
    "physiologically based pharmacokinetic (mPBPK) model for",
    "transferrin-receptor (TfR) mediated brain delivery of monoclonal",
    "antibodies in non-human primates. 26-compartment NONMEM ADVAN8",
    "structure combining the Bloomingdale 2017 mAb mPBPK framework",
    "(plasma, tissue vascular / endosomal / interstitial / FcRn,",
    "brain vascular, BBB endosomal (unbound + FcRn-bound), brain ISF,",
    "BCSFB endosomal (unbound + FcRn-bound), CSF, lymph, with FcRn",
    "recycling) and the Chang 2022 whole-body-plasma + brain-vascular",
    "+ ISF + neuronal TfR binding with empirical TMDD-style",
    "elimination of the bsAb-TfR complex (kint). Two TfR binding sites",
    "on the brain barriers (luminal BBB, luminal BCSFB) transcytose",
    "bound complex into the abluminal side (brain ISF and CSF",
    "respectively), where it may dissociate or be degraded.",
    "Parameters were fit to 395 plasma, 81 CSF, and 102 brain mean",
    "concentrations digitised from eight literature studies in",
    "cynomolgus monkey (7 non-TfR mAbs + 10 anti-TfR bsAbs with KD,TfR",
    "36-1900 nM). The kint mixture (POP1 fast 0.0329 h^-1 fraction",
    "0.437; POP2 slow 0.0125 h^-1) is selected per subject via the",
    "MIX_FAST_ELIM covariate (1 = POP1, 0 = POP2). Per-compound TfR",
    "binding parameters (kon_T, koff_T) are NOT estimated population",
    "values - they are set per simulated antibody from biophysical",
    "measurements (paper Table S2); the default ini() encodes kon_T =",
    "0 (non-TfR control IgG). Plasma observation is TOTAL drug (free +",
    "TfR-bound complex), CSF observation is unbound CBCSF, brain",
    "observation is whole-brain homogenate (simplified volume-weighted",
    "average across BBB + ISF + BCSFB endosomal spaces, scaled by the",
    "estimated FACBR correction factor 0.05). No inter-individual",
    "variability was estimated (dataset was mean digitised profiles)."
  )
  reference <- paste(
    "Muliaditan M, van Steeg TJ, Avery LB, Sun W, Hammond TR, Hijdra D,",
    "Choi SL, Pillai N, Leksa NC, Mavroudis PD.",
    "Translational minimal physiologically based pharmacokinetic model",
    "for transferrin receptor-mediated brain delivery of antibodies.",
    "mAbs. 2025;17(1):2515414. doi:10.1080/19420862.2025.2515414.",
    "Bloomingdale physiology adapted from Bloomingdale et al. 2017",
    "J Pharmacokinet Pharmacodyn (ref 25); TfR-binding structure adapted",
    "from Chang et al. 2022 (ref 18); per-compound TfR binding parameters",
    "from Muliaditan 2025 Supplementary Table S2.",
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
        "(kint = 0.0329 h^-1); 0 = POP2 slow eliminator (kint =",
        "0.0125 h^-1). Not a measured covariate - this is the",
        "$MIXTURE class assignment from the Muliaditan 2025 NONMEM run."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "Population probability of MIX_FAST_ELIM = 1 (POP1, fast",
        "internalization) is 0.437 per Muliaditan 2025 Table 2",
        "(Fraction POP1, RSE 30%, 95% CI 0.180-0.695). For typical-",
        "value simulation set MIX_FAST_ELIM = 0 (slow, dominant",
        "subpopulation 56.3%). For population simulation, draw",
        "MIX_FAST_ELIM ~ Bernoulli(0.437) per subject. The two",
        "populations were not correlated with KD,TfR or data source",
        "in the original fit (paper Discussion)."
      ),
      source_name        = "MIXTURE"
    )
  )

  population <- list(
    species         = "cynomolgus monkey",
    n_subjects      = NA_integer_,
    n_studies       = 8L,
    n_compounds     = 17L,
    n_plasma_obs    = 395L,
    n_csf_obs       = 81L,
    n_brain_obs     = 102L,
    weight_range    = "Bloomingdale 2017 6.2-kg cynomolgus reference subject for the physiology",
    disease_state   = "healthy non-human primate (cynomolgus monkey)",
    dose_range      = "1-100 mg/kg IV bolus or IV infusion (training data: 2, 10, 30, 50 mg/kg single IV; Edavettal 2022 also includes a weekly RD regimen)",
    notes           = paste(
      "Training dataset (paper Table 1) covers 7 non-TfR mAbs (Control",
      "IgG, anti-Tau-IgG, Lu-AF82422, gantenerumab, anti-BACE1, Yu 2014",
      "Control IgG, Kariolis 2020 Control IgG) and 10 anti-TfR bsAbs",
      "(TfR-K-mut, TfR-J-wt, TfR-J-mut, trontinemab, anti-TfR/BACE1 at",
      "3 affinities, ATV35.21.16:BACE1, Anti-TfR1-BACE1, Anti-TfR2-BACE1)",
      "with KD,TfR ranging from 0 to 1900 nM. Inter-individual",
      "variability was not evaluated - all observations are mean",
      "literature digitisations except Kariolis 2020 individual",
      "profiles."
    )
  )

  ini({
    # ============================================================
    # Estimated structural parameters (Muliaditan 2025 Table 2).
    # NHP NONMEM ADVAN8 run, FOCE-I; no inter-individual variability;
    # additive residual error on natural-log-transformed concentrations.
    # ============================================================

    # Total whole-body-minus-brain TfR concentration in plasma (nM).
    ltfrpt          <- log(1672);    label("Total plasma TfR concentration TfRpt (nM)")              # Muliaditan 2025 Table 2 (1672 nM; RSE 38%)
    # Baseline luminal unbound TfR on the BBB (nM).
    lutfr0_bbb      <- log(175);     label("Baseline luminal unbound TfR on BBB, uTFR0_BBB (nM)")    # Muliaditan 2025 Table 2 (175 nM; RSE 34%)
    # Baseline luminal unbound TfR on the BCSFB (nM).
    lutfr0_bcsfb    <- log(0.256);   label("Baseline luminal unbound TfR on BCSFB, uTFR0_BCSFB (nM)")# Muliaditan 2025 Table 2 (0.256 nM; RSE 47%)
    # Fast subpopulation internalization rate (h^-1).
    lkint_pop1      <- log(0.0329);  label("bsAb-TfR complex internalization rate kint POP1 (h^-1, fast subpop)") # Muliaditan 2025 Table 2 (0.0329 h^-1; RSE 14%)
    # Slow subpopulation internalization rate (h^-1).
    lkint_pop2      <- log(0.0125);  label("bsAb-TfR complex internalization rate kint POP2 (h^-1, slow subpop)") # Muliaditan 2025 Table 2 (0.0125 h^-1; RSE 8.4%)
    # Correction factor restricting unbound drug distribution from brain
    # ISF to CSF (unitless multiplier on Q_BECF for the ISF -> CSF flux).
    lfacq_becf      <- log(0.00814); label("Correction factor on Q_BECF for ISF -> CSF distribution (unitless)") # Muliaditan 2025 Table 2 (0.00814; RSE 29%)

    # ============================================================
    # Fixed structural parameters (Muliaditan 2025 Table 2 a-footer).
    # ============================================================
    # Total neuronal TfR in brain ISF, fixed from Chang 2022.
    ltfrtotn        <- fixed(log(559));  label("Total neuronal TfR in brain ISF TfRtotn (nM)")        # Muliaditan 2025 Table 2 / Chang 2022 (559 nM, FIXED)
    # Transcytosis rate of bsAb-TfR complex across BBB and BCSFB, fixed.
    lktrans         <- fixed(log(6));    label("Transcytosis rate ktrans across BBB and BCSFB (h^-1)") # Muliaditan 2025 Table 2 / Chang 2022 (6 h^-1, FIXED)
    # Luminal-BBB unbound TfR degradation rate, fixed.
    lkdeg_utfr_bbb  <- fixed(log(20));   label("Degradation rate of luminal unbound TfR on BBB (h^-1)") # Muliaditan 2025 Table 2 / Chang 2022 (20 h^-1, FIXED)
    # Luminal-BCSFB unbound TfR degradation rate, fixed to 20 (paper
    # text: in contrast to Chang 2022's 1.42 h^-1, this paper fixed
    # to 20 h^-1 to prevent total BCSFB TfR from decreasing).
    lkdeg_utfr_bcsfb <- fixed(log(20));  label("Degradation rate of luminal unbound TfR on BCSFB (h^-1, set to BBB value)") # Muliaditan 2025 Table 2 (20 h^-1, FIXED; contrast Chang 2022's 1.42)
    # Brain-prediction correction factor (FAC_BPRED = FACBR in supp eqs).
    lfacbr          <- fixed(log(0.05)); label("Brain-prediction correction factor FAC_BPRED (unitless)") # Muliaditan 2025 Table 2 (0.05, FIXED)

    # ============================================================
    # Per-compound TfR binding parameters. NOT estimated in the
    # population fit - they are biophysical inputs taken from each
    # compound's binding assay (Muliaditan 2025 Supplementary Table S2).
    # Default encoding here represents a non-TfR control IgG
    # (kon_T = 0, no TfR binding). To simulate a specific anti-TfR
    # bsAb, override these via the ini() argument to rxSolve():
    #   trontinemab (NHP): kon_T = 0.21168 nM^-1 h^-1, koff_T = 52.56 h^-1
    #   TfR-J-wt:          kon_T = 0.5832  nM^-1 h^-1, koff_T = 21.06 h^-1
    #   ATV35.21.16:BACE1: kon_T = 0.5193  nM^-1 h^-1, koff_T = 986.67 h^-1
    #   See vignette and Muliaditan 2025 Table S2.
    # ============================================================
    lkon_t   <- fixed(log(1e-12)); label("TfR association rate kon_T (nM^-1 h^-1) - per-compound input; 1e-12 placeholder for non-TfR control") # Muliaditan 2025 Table S2; default = non-TfR control
    lkoff_t  <- fixed(log(1));     label("TfR dissociation rate koff_T (h^-1) - per-compound input")                                              # Muliaditan 2025 Table S2; default placeholder

    # Antibody molecular weight for dose unit conversion (kDa).
    # Default 150 kDa (typical IgG). Override via ini() for non-150-kDa
    # constructs (e.g., trontinemab 194 kDa per Muliaditan 2025 Table S2).
    mw_kda   <- fixed(150);        label("Antibody molecular weight (kDa) for mg -> nmol dose conversion") # Muliaditan 2025 Table S2 typical IgG MW

    # ============================================================
    # Residual error - additive on the log-transformed concentration
    # (paper Methods: yi,k = PRED + eps, var(eps) = sigma^2). Maps to
    # nlmixr2 lnorm() with expSd = sqrt(sigma^2).
    # ============================================================
    expSd        <- sqrt(0.256); label("Log-scale residual SD for plasma total drug (unitless)")     # Muliaditan 2025 Table 2, sigma^2_plasma = 0.256
    expSd_Ccsf   <- sqrt(0.893); label("Log-scale residual SD for CSF unbound drug (unitless)")      # Muliaditan 2025 Table 2, sigma^2_CSF = 0.893
    expSd_Cbrain <- sqrt(0.426); label("Log-scale residual SD for whole-brain homogenate (unitless)")# Muliaditan 2025 Table 2, sigma^2_brain = 0.426
  })

  model({
    # ============================================================
    # 1. Resolve estimated/fixed parameters from log space.
    # ============================================================
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

    # Mixture selection: MIX_FAST_ELIM = 1 -> POP1 (fast); 0 -> POP2 (slow).
    # Default reference category 0 = POP2 slow (dominant 56.3% subpop).
    kint <- kint_pop2 + (kint_pop1 - kint_pop2) * MIX_FAST_ELIM
    # Recycling rate constant of unbound TfR from brain ISF/CSF back to
    # luminal BBB/BCSFB is FIXED equal to kint (paper Methods / Table 2).
    krec_utfr <- kint

    # Synthesis rates derived from steady-state balance (no drug):
    #   d[uTfR]/dt = 0 -> Ksyn = Kdeg * uTfR_baseline.
    ksyn_utfr_bbb   <- kdeg_utfr_bbb   * utfr0_bbb
    ksyn_utfr_bcsfb <- kdeg_utfr_bcsfb * utfr0_bcsfb

    # ============================================================
    # 2. Fixed NHP physiology from Bloomingdale 2017 (Muliaditan 2025
    #    Supplementary Table S1). Plasma flow QT here uses the
    #    volume-weighted-average value reported in the supplement
    #    footnote (a) rather than the original Bloomingdale QT.
    # ============================================================
    # Volumes (L)
    vp           <- 0.187           # plasma                                                  # Muliaditan 2025 Supp Table S1
    v_tissue_v   <- 0.151           # tissue vascular VTv                                     # Muliaditan 2025 Supp Table S1
    v_tissue_e   <- 0.0286          # tissue endosomal VTe                                    # Muliaditan 2025 Supp Table S1
    v_tissue_i   <- 0.976           # tissue interstitial VTi                                 # Muliaditan 2025 Supp Table S1
    v_brain_v    <- 0.00207         # brain vascular VBv                                      # Muliaditan 2025 Supp Table S1
    v_bbb_e      <- 4.27e-4         # BBB endosomal V_BE_BBB                                  # Muliaditan 2025 Supp Table S1
    v_bcsfb_e    <- 4.27e-5         # BCSFB endosomal V_BE_BCSFB                              # Muliaditan 2025 Supp Table S1
    v_brain_i    <- 0.0169          # brain interstitial VBi                                  # Muliaditan 2025 Supp Table S1
    v_csf        <- 0.00926         # CSF VCSF                                                # Muliaditan 2025 Supp Table S1
    v_lymph      <- 0.0251          # lymph VL                                                # Muliaditan 2025 Supp Table S1

    # Flows (L/h) - QT uses the volume-weighted average (Table S1 footnote a).
    qt           <- 2.9             # tissue plasma flow                                      # Muliaditan 2025 Supp Table S1 footnote a (volume-weighted avg of non-brain tissues)
    qb           <- 1.51            # brain plasma flow                                       # Muliaditan 2025 Supp Table S1
    lt           <- 0.0419          # tissue lymphatic flow                                   # Muliaditan 2025 Supp Table S1
    lb           <- 0.00369         # brain lymphatic flow                                    # Muliaditan 2025 Supp Table S1
    qb_ecf       <- 0.00123         # brain ECF flow                                          # Muliaditan 2025 Supp Table S1
    qb_csf       <- 0.00246         # brain CSF flow                                          # Muliaditan 2025 Supp Table S1

    # Uptake clearance rate constants
    k_clupt      <- 0.55            # tissue endosomal uptake rate constant                   # Muliaditan 2025 Supp Table S1
    k_clupb      <- 0.0195          # brain endosomal uptake rate constant                    # Muliaditan 2025 Supp Table S1
    clup_t       <- 0.0157          # tissue uptake clearance CLUPT                           # Muliaditan 2025 Supp Table S1
    clup_bbb     <- 8.35e-6         # BBB uptake clearance CLUPBBB                            # Muliaditan 2025 Supp Table S1
    clup_bcsfb   <- 8.35e-7         # BCSFB uptake clearance CLUPBCSFB                        # Muliaditan 2025 Supp Table S1

    # Reflection coefficients (sigma)
    sig_bbb      <- 1               # BBB reflection coefficient                              # Muliaditan 2025 Supp Table S1
    sig_bcsfb    <- 0.9974          # BCSFB reflection coefficient                            # Muliaditan 2025 Supp Table S1
    sig_tv       <- 0.9239          # tissue vascular reflection coefficient                  # Muliaditan 2025 Supp Table S1
    sig_tl       <- 0.2             # tissue lymph reflection coefficient                     # Muliaditan 2025 Supp Table S1
    sig_isf      <- 0.2             # brain ISF reflection coefficient                        # Muliaditan 2025 Supp Table S1
    sig_csf      <- 0.2             # CSF reflection coefficient                              # Muliaditan 2025 Supp Table S1

    # Fractional recycling
    fr_tissue    <- 0.715           # tissue FcRn fractional recycling FR                     # Muliaditan 2025 Supp Table S1
    fr_brain     <- 0.715           # brain FcRn fractional recycling FR_B                    # Muliaditan 2025 Supp Table S1

    # FcRn baseline concentration (initial condition in endosomes).
    fcrn_baseline_M <- 4.98e-5      # mol/L                                                   # Muliaditan 2025 Supp Table S1
    # Convert to nM to match the in-model concentration units.
    fcrn_baseline   <- fcrn_baseline_M * 1e9   # = 49800 nM

    # FcRn binding rate constants.
    # Source values in M^-1 h^-1 (kon) and h^-1 (koff); convert kon to nM^-1 h^-1.
    kon_fcrn_M  <- 7.92e8           # M^-1 h^-1                                               # Muliaditan 2025 Supp Table S1
    kon_fcrn    <- kon_fcrn_M * 1e-9  # nM^-1 h^-1
    koff_fcrn   <- 46.8             # h^-1                                                    # Muliaditan 2025 Supp Table S1
    kdeg_endo   <- 26.6             # h^-1 endosomal antibody degradation                     # Muliaditan 2025 Supp Table S1

    # ============================================================
    # 3. Dose unit conversion (mg -> nmol). f(central) = 1000/mw_kda
    #    so a dose entered as mg becomes nmol inside the compartment.
    # ============================================================
    f(central) <- 1000 / mw_kda

    # ============================================================
    # 4. Concentrations (nM) derived from compartment amounts (nmol).
    #    Compartments hold amount in nmol; c_X = amount_X / V_X.
    # ============================================================
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
    c_delta_utfr_bbb  <- delta_utfr_bbb   / v_brain_v   # relative-level state on BBB
    c_complex_bbb_lum <- complex_bbb_lum  / v_brain_v
    c_complex_bbb_abl <- complex_bbb_abl  / v_brain_i
    c_utfr_bbb_abl    <- utfr_bbb_abl     / v_brain_i
    c_delta_utfr_bcsfb<- delta_utfr_bcsfb / v_brain_v   # relative-level state on BCSFB
    c_complex_bcsfb_lum <- complex_bcsfb_lum / v_brain_v
    c_complex_bcsfb_abl <- complex_bcsfb_abl / v_csf
    c_utfr_bcsfb_abl    <- utfr_bcsfb_abl    / v_csf
    c_complex_neuron    <- complex_neuron    / v_brain_i

    # Total unbound TfR on BBB / BCSFB luminal surfaces:
    #   total = baseline + delta
    utfr_bbb_lum    <- utfr0_bbb   + c_delta_utfr_bbb
    utfr_bcsfb_lum  <- utfr0_bcsfb + c_delta_utfr_bcsfb

    # Free TfR available in plasma and on neurons (conservation):
    utfr_plasma     <- tfrpt   - c_complex_plasma
    utfr_neuron     <- tfrtotn - c_complex_neuron

    # ============================================================
    # 5. TfR binding flux helpers (nM/h, concentration rates).
    #    Multiply by compartment volume to get mass rates (nmol/h).
    # ============================================================
    # Plasma: free drug + free plasma TfR <-> complex_plasma
    tfrb_plasma_rate <- -kon_t * c_plasma     * utfr_plasma  + koff_t * c_complex_plasma
    # Luminal BBB: free brain vascular drug + free luminal-BBB TfR <-> complex_bbb_lum
    tfrb_bbb_lum_rate <- -kon_t * c_brain_vascular * utfr_bbb_lum    + koff_t * c_complex_bbb_lum
    # Luminal BCSFB (uses the same brain vascular drug pool):
    tfrb_bcsfb_lum_rate <- -kon_t * c_brain_vascular * utfr_bcsfb_lum + koff_t * c_complex_bcsfb_lum
    # Abluminal BBB: free brain ISF drug + abluminal BBB free TfR <-> complex_bbb_abl
    tfrb_bbb_abl_rate <- -kon_t * c_brain_isf * c_utfr_bbb_abl   + koff_t * c_complex_bbb_abl
    # Neuron: free brain ISF drug + free neuronal TfR <-> complex_neuron
    tfrb_neuron_rate  <- -kon_t * c_brain_isf * utfr_neuron      + koff_t * c_complex_neuron
    # Abluminal BCSFB: free CSF drug + abluminal BCSFB free TfR <-> complex_bcsfb_abl
    tfrb_bcsfb_abl_rate <- -kon_t * c_csf     * c_utfr_bcsfb_abl + koff_t * c_complex_bcsfb_abl

    # ============================================================
    # 6. ODE system (paper Supplementary "Model equations"), translated
    #    from concentration rates (NONMEM DADT in nM/h) to mass rates
    #    (rxode2 d/dt in nmol/h) by multiplying the local concentration-
    #    rate terms by the compartment volume where appropriate.
    # ============================================================

    # (1) Plasma free drug.
    d/dt(central) <- (qt - lt) * c_tissue_vasc + (qb - lb) * c_brain_vascular +
                     (lt + lb) * c_lymph -
                     qt * c_plasma - qb * c_plasma +
                     tfrb_plasma_rate * vp

    # (2) Tissue vascular.
    d/dt(tissue_vasc) <- qt * c_plasma - (qt - lt) * c_tissue_vasc -
                         (1 - sig_tv) * lt * c_tissue_vasc -
                         clup_t * c_tissue_vasc +
                         clup_t * fr_tissue * c_tissue_endo_b

    # (3) Tissue endosomal unbound.
    d/dt(tissue_endo_u) <- clup_t * (c_tissue_vasc + c_tissue_isf) -
                           (kon_fcrn * c_tissue_endo_u * c_tissue_fcrn -
                            koff_fcrn * c_tissue_endo_b) * v_tissue_e -
                           kdeg_endo * c_tissue_endo_u * v_tissue_e

    # (4) Tissue endosomal FcRn-bound.
    d/dt(tissue_endo_b) <- (kon_fcrn * c_tissue_endo_u * c_tissue_fcrn -
                            koff_fcrn * c_tissue_endo_b) * v_tissue_e -
                           clup_t * c_tissue_endo_b

    # (5) Tissue interstitial.
    d/dt(tissue_isf) <- (1 - sig_tv) * lt * c_tissue_vasc -
                        (1 - sig_tl) * lt * c_tissue_isf +
                        clup_t * (1 - fr_tissue) * c_tissue_endo_b -
                        clup_t * c_tissue_isf

    # (6) Brain vascular.
    d/dt(brain_vascular) <- qb * c_plasma - (qb - lb) * c_brain_vascular -
                        (1 - sig_bbb)   * qb_ecf * c_brain_vascular -
                        (1 - sig_bcsfb) * qb_csf * c_brain_vascular -
                        k_clupb * c_brain_vascular * v_brain_v +
                        clup_bbb   * fr_brain * c_bbb_endo_b +
                        clup_bcsfb * fr_brain * c_bcsfb_endo_b +
                        (tfrb_bbb_lum_rate + tfrb_bcsfb_lum_rate) * v_brain_v

    # (7) BBB endosomal unbound.
    d/dt(bbb_endo_u) <- clup_bbb * (c_brain_vascular + c_brain_isf) -
                        (kon_fcrn * c_bbb_endo_u * c_bbb_fcrn -
                         koff_fcrn * c_bbb_endo_b) * v_bbb_e -
                        kdeg_endo * c_bbb_endo_u * v_bbb_e

    # (8) BBB endosomal FcRn-bound.
    d/dt(bbb_endo_b) <- (kon_fcrn * c_bbb_endo_u * c_bbb_fcrn -
                         koff_fcrn * c_bbb_endo_b) * v_bbb_e -
                        clup_bbb * c_bbb_endo_b

    # (9) Brain ISF.
    d/dt(brain_isf) <- (1 - sig_bbb)  * qb_ecf * c_brain_vascular -
                       (1 - sig_isf) * qb_ecf * c_brain_isf +
                       clup_bbb * (1 - fr_brain) * c_bbb_endo_b -
                       clup_bbb * c_brain_isf -
                       qb_ecf * facq_becf * c_brain_isf +
                       qb_ecf * c_csf +
                       (tfrb_bbb_abl_rate + tfrb_neuron_rate) * v_brain_i

    # (10) BCSFB endosomal unbound.
    d/dt(bcsfb_endo_u) <- clup_bcsfb * (c_brain_vascular + c_csf) -
                          (kon_fcrn * c_bcsfb_endo_u * c_bcsfb_fcrn -
                           koff_fcrn * c_bcsfb_endo_b) * v_bcsfb_e -
                          kdeg_endo * c_bcsfb_endo_u * v_bcsfb_e

    # (11) BCSFB endosomal FcRn-bound.
    d/dt(bcsfb_endo_b) <- (kon_fcrn * c_bcsfb_endo_u * c_bcsfb_fcrn -
                           koff_fcrn * c_bcsfb_endo_b) * v_bcsfb_e -
                          clup_bcsfb * c_bcsfb_endo_b

    # (12) CSF unbound.
    d/dt(csf) <- (1 - sig_bcsfb) * qb_csf * c_brain_vascular -
                 clup_bcsfb * c_csf +
                 clup_bcsfb * (1 - fr_brain) * c_bcsfb_endo_b +
                 qb_ecf * facq_becf * c_brain_isf -
                 (1 - sig_csf) * qb_csf * c_csf -
                 qb_ecf * c_csf +
                 tfrb_bcsfb_abl_rate * v_csf

    # (13) Lymph.
    d/dt(lymph) <- (1 - sig_tl)  * lt     * c_tissue_isf +
                   (1 - sig_csf) * qb_csf * c_csf +
                   (1 - sig_isf) * qb_ecf * c_brain_isf -
                   (lt + lb) * c_lymph

    # (14) FcRn free in tissue endosome - mirror of (4) sign-flipped.
    d/dt(tissue_fcrn) <- -(kon_fcrn * c_tissue_endo_u * c_tissue_fcrn -
                           koff_fcrn * c_tissue_endo_b) * v_tissue_e +
                          clup_t * c_tissue_endo_b

    # (15) FcRn free in BBB endosome.
    d/dt(bbb_fcrn) <- -(kon_fcrn * c_bbb_endo_u * c_bbb_fcrn -
                        koff_fcrn * c_bbb_endo_b) * v_bbb_e +
                       clup_bbb * c_bbb_endo_b

    # (16) FcRn free in BCSFB endosome.
    d/dt(bcsfb_fcrn) <- -(kon_fcrn * c_bcsfb_endo_u * c_bcsfb_fcrn -
                          koff_fcrn * c_bcsfb_endo_b) * v_bcsfb_e +
                         clup_bcsfb * c_bcsfb_endo_b

    # (17) Plasma bsAb-TfR complex.
    d/dt(complex_plasma) <- -tfrb_plasma_rate * vp - kint * complex_plasma

    # (18) Delta luminal unbound TfR on BBB (relative level, units nM).
    #   State represents (uTfR_actual - uTFR0_BBB); compartment volume
    #   is treated as v_brain_v so that the dimensional rxode2 amount
    #   d/dt(.) is (nM/h * L) = nmol/h.
    d/dt(delta_utfr_bbb) <- (ksyn_utfr_bbb -
                             kdeg_utfr_bbb * utfr_bbb_lum +
                             tfrb_bbb_lum_rate +
                             krec_utfr * (v_brain_i / v_brain_v) * c_utfr_bbb_abl) * v_brain_v

    # (19) Luminal BBB bsAb-TfR complex.
    d/dt(complex_bbb_lum) <- (-tfrb_bbb_lum_rate -
                              ktrans * c_complex_bbb_lum -
                              kint   * c_complex_bbb_lum) * v_brain_v

    # (20) Abluminal BBB bsAb-TfR complex.
    d/dt(complex_bbb_abl) <- ( ktrans * (v_brain_v / v_brain_i) * c_complex_bbb_lum -
                              tfrb_bbb_abl_rate -
                              kint   * c_complex_bbb_abl) * v_brain_i

    # (21) Abluminal BBB free TfR.
    d/dt(utfr_bbb_abl) <- ( tfrb_bbb_abl_rate -
                            krec_utfr * c_utfr_bbb_abl) * v_brain_i

    # (22) Delta luminal unbound TfR on BCSFB (relative level).
    d/dt(delta_utfr_bcsfb) <- (ksyn_utfr_bcsfb -
                               kdeg_utfr_bcsfb * utfr_bcsfb_lum +
                               tfrb_bcsfb_lum_rate +
                               krec_utfr * (v_csf / v_brain_v) * c_utfr_bcsfb_abl) * v_brain_v

    # (23) Luminal BCSFB bsAb-TfR complex.
    d/dt(complex_bcsfb_lum) <- (-tfrb_bcsfb_lum_rate -
                                ktrans * c_complex_bcsfb_lum -
                                kint   * c_complex_bcsfb_lum) * v_brain_v

    # (24) Abluminal BCSFB bsAb-TfR complex.
    d/dt(complex_bcsfb_abl) <- ( ktrans * (v_brain_v / v_csf) * c_complex_bcsfb_lum -
                                tfrb_bcsfb_abl_rate -
                                kint   * c_complex_bcsfb_abl) * v_csf

    # (25) Abluminal BCSFB free TfR.
    d/dt(utfr_bcsfb_abl) <- ( tfrb_bcsfb_abl_rate -
                              krec_utfr * c_utfr_bcsfb_abl) * v_csf

    # (26) Neuronal bsAb-TfR complex.
    d/dt(complex_neuron) <- (-tfrb_neuron_rate - kint * c_complex_neuron) * v_brain_i

    # ============================================================
    # 7. Initial conditions. All compartments start at 0 except the
    #    free FcRn pools (paper Supplement: "initial condition for all
    #    compartments is 0, except for unbound FcRn (A14-A16)").
    # ============================================================
    tissue_fcrn(0) <- fcrn_baseline * v_tissue_e
    bbb_fcrn(0)    <- fcrn_baseline * v_bbb_e
    bcsfb_fcrn(0)  <- fcrn_baseline * v_bcsfb_e

    # ============================================================
    # 8. Observations.
    #   - Plasma observation = TOTAL drug concentration (free + TfR-
    #     bound complex), nM, per paper Methods.
    #   - CSF observation = unbound CBCSF, nM.
    #   - Brain observation = whole-brain homogenate concentration,
    #     nM, computed as a volume-weighted average across the BBB
    #     endosomal compartments, brain ISF, and BCSFB endosomal
    #     compartments, with the FACBR correction. The published
    #     equation (Muliaditan 2025 Supplement "Predicted drug
    #     concentration in whole brain homogenate") also includes a
    #     ventricular-CSF contribution A(12)*(VCSFLV+VCSFTFV) divided
    #     by (VBtotal + VCSFLV + VCSFTFV); those Bloomingdale
    #     ventricular-CSF volume parameters (VCSFLV, VCSFTFV, VBtotal)
    #     are not reported on disk in the source materials, so the
    #     brain homogenate here is the simplified volume-weighted
    #     average over the (BBB endosomal + brain ISF + BCSFB
    #     endosomal) compartments. The simplification is documented
    #     in vignette Errata.
    # ============================================================
    Cc      <- c_plasma + c_complex_plasma             # plasma total drug
    Ccsf    <- c_csf                                   # CSF unbound

    homo_num <- (c_bbb_endo_u   + c_bbb_endo_b)   * v_bbb_e   +
                 c_brain_isf                      * v_brain_i +
                (c_bcsfb_endo_u + c_bcsfb_endo_b) * v_bcsfb_e
    homo_den <- v_bbb_e + v_brain_i + v_bcsfb_e
    Cbrain   <- facbr * homo_num / homo_den            # whole-brain homogenate

    Cc     ~ lnorm(expSd)
    Ccsf   ~ lnorm(expSd_Ccsf)
    Cbrain ~ lnorm(expSd_Cbrain)
  })
}
