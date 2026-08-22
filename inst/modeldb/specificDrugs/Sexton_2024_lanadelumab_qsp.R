Sexton_2024_lanadelumab_qsp <- function() {
  description <- "QSP. Sexton 2024 plasma kallikrein-kinin system (KKS) model of hereditary angioedema (HAE) due to C1-inhibitor deficiency: 61 ODE states spanning a vascular (plasma) space, an endothelial-cell proximal space and a gC1q-R / bradykinin-B2 receptor surface, coupled to a one-compartment subcutaneous popPK for lanadelumab. HAE attacks are driven by a transient fold-increase in the FXII autoactivation rate delivered as a dose into the 'trigger' compartment; an attack is scored when plasma bradykinin exceeds 20 pM. Exogenous C1-INH may be given as a dose into the 'c1inh' compartment (nM increment), reproducing the published flux_C1Inh_inj term."
  reference <- paste(
    "Sexton D, Nguyen HQ, Juethner S, Luo H, Zhang Z, Jasper P, Zhu AZX.",
    "A quantitative systems pharmacology model of plasma kallikrein-kinin system",
    "dysregulation in hereditary angioedema.",
    "J Pharmacokinet Pharmacodyn. 2024;51(6):721-733.",
    "doi:10.1007/s10928-024-09919-6.",
    "Species from Supplementary Table S4; governing equations from Supplementary Table S5",
    "and the published C source (Electronic Supplementary Material MOESM2,",
    "HAE_Contact_System_Model.c); parameter values from Supplementary Table S6 and the",
    "published R parameter file (MOESM5, Param_HAE.r); initial conditions from MOESM3",
    "(Init_HAE.r); lanadelumab popPK, body-weight model and attack-severity distribution",
    "from MOESM7 (SimInput_HAE.r) and MOESM6 (virtual-patient sampler).",
    sep = " "
  )
  vignette <- "Sexton_2024_lanadelumab_qsp"
  units <- list(
    time          = "h",
    dosing        = "mg (lanadelumab SC into depot); nM increment (exogenous C1-INH into c1inh); unitless normalized attack severity (into trigger)",
    concentration = "ug/mL (lanadelumab plasma concentration Cc); nM (all KKS species); number/cell (surface receptor species)"
  )

  # Non-canonical, paper-mechanistic states. Names follow Supplementary Table S4
  # with the paper's "_in_plasma" suffix dropped for the vascular space, a "px_"
  # prefix for the endothelial proximal space, a "_gc1qr" / "_bdkrb2" suffix for
  # surface-receptor species, and a "_deg" suffix for the degradation sinks.
  paper_specific_compartments <- c(
    "trigger",
    # vascular (plasma) space -- nM
    "fxii", "fxiia", "c1inh", "hmwk", "chmwk", "prekal", "kal", "bk",
    "prekal_hmwk", "kal_hmwk", "kal_chmwk", "c1inh_fxiia", "c1inh_kal",
    "c1inh_kal_hmwk", "lana_kal", "lana_kal_hmwk",
    # endothelial proximal space -- nM
    "px_fxii", "px_fxiia", "px_c1inh", "px_prekal_hmwk", "px_kal_hmwk",
    "px_c1inh_fxiia", "px_c1inh_kal_hmwk", "px_kal_chmwk", "px_bk",
    "px_lana_kal_hmwk",
    # endothelial cell surface -- number/cell
    "gc1qr", "bdkrb2", "fxii_gc1qr", "fxiia_gc1qr", "prekal_hmwk_gc1qr",
    "kal_hmwk_gc1qr", "kal_chmwk_gc1qr", "c1inh_fxiia_gc1qr",
    "c1inh_kal_hmwk_gc1qr", "lana_kal_hmwk_gc1qr", "bk_bdkrb2",
    # degradation sinks (Table S4 lists these as model species; they carry the
    # mass-balance bookkeeping used by the vignette's mass-balance gate)
    "prekal_deg", "kal_deg", "hmwk_deg", "chmwk_deg", "c1inh_deg", "fxii_deg",
    "fxiia_deg", "bk_deg", "c1inh_kal_deg", "c1inh_kal_hmwk_deg",
    "c1inh_fxiia_deg", "gc1qr_deg", "bdkrb2_deg", "bk_bdkrb2_deg",
    "prekal_hmwk_gc1qr_deg", "kal_hmwk_gc1qr_deg", "kal_chmwk_gc1qr_deg",
    "fxii_gc1qr_deg", "fxiia_gc1qr_deg", "c1inh_fxiia_gc1qr_deg",
    "c1inh_kal_hmwk_gc1qr_deg", "lana_kal_hmwk_gc1qr_deg"
  )

  # Analyte / units / specimen per Supplementary Table S4 (verified against that table).
  compartmentData <- list(
    depot   = list(analyte = "lanadelumab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "lanadelumab", units = "mg", specimen = "plasma", verified = TRUE),
    trigger = list(analyte = "normalized HAE attack severity (fold-increase driver for FXII autoactivation)", units = "unitless", specimen = "not applicable", verified = TRUE),

    fxii           = list(analyte = "factor XII", units = "nM", specimen = "plasma", verified = TRUE),
    fxiia          = list(analyte = "activated factor XII", units = "nM", specimen = "plasma", verified = TRUE),
    c1inh          = list(analyte = "C1 esterase inhibitor", units = "nM", specimen = "plasma", verified = TRUE),
    hmwk           = list(analyte = "high molecular weight kininogen", units = "nM", specimen = "plasma", verified = TRUE),
    chmwk          = list(analyte = "cleaved high molecular weight kininogen", units = "nM", specimen = "plasma", verified = TRUE),
    prekal         = list(analyte = "prekallikrein", units = "nM", specimen = "plasma", verified = TRUE),
    kal            = list(analyte = "kallikrein", units = "nM", specimen = "plasma", verified = TRUE),
    bk             = list(analyte = "bradykinin", units = "nM", specimen = "plasma", verified = TRUE),
    prekal_hmwk    = list(analyte = "prekallikrein-HMWK complex", units = "nM", specimen = "plasma", verified = TRUE),
    kal_hmwk       = list(analyte = "kallikrein-HMWK complex", units = "nM", specimen = "plasma", verified = TRUE),
    kal_chmwk      = list(analyte = "kallikrein-cHMWK complex", units = "nM", specimen = "plasma", verified = TRUE),
    c1inh_fxiia    = list(analyte = "C1-INH-FXIIa complex", units = "nM", specimen = "plasma", verified = TRUE),
    c1inh_kal      = list(analyte = "C1-INH-kallikrein complex", units = "nM", specimen = "plasma", verified = TRUE),
    c1inh_kal_hmwk = list(analyte = "C1-INH-kallikrein-HMWK complex", units = "nM", specimen = "plasma", verified = TRUE),
    lana_kal       = list(analyte = "lanadelumab-kallikrein complex", units = "nM", specimen = "plasma", verified = TRUE),
    lana_kal_hmwk  = list(analyte = "lanadelumab-kallikrein-HMWK complex", units = "nM", specimen = "plasma", verified = TRUE),

    px_fxii            = list(analyte = "factor XII (endothelial proximal space)", units = "nM", specimen = "tissue", verified = TRUE),
    px_fxiia           = list(analyte = "activated factor XII (endothelial proximal space)", units = "nM", specimen = "tissue", verified = TRUE),
    px_c1inh           = list(analyte = "C1 esterase inhibitor (endothelial proximal space)", units = "nM", specimen = "tissue", verified = TRUE),
    px_prekal_hmwk     = list(analyte = "prekallikrein-HMWK complex (endothelial proximal space)", units = "nM", specimen = "tissue", verified = TRUE),
    px_kal_hmwk        = list(analyte = "kallikrein-HMWK complex (endothelial proximal space)", units = "nM", specimen = "tissue", verified = TRUE),
    px_c1inh_fxiia     = list(analyte = "C1-INH-FXIIa complex (endothelial proximal space)", units = "nM", specimen = "tissue", verified = TRUE),
    px_c1inh_kal_hmwk  = list(analyte = "C1-INH-kallikrein-HMWK complex (endothelial proximal space)", units = "nM", specimen = "tissue", verified = TRUE),
    px_kal_chmwk       = list(analyte = "kallikrein-cHMWK complex (endothelial proximal space)", units = "nM", specimen = "tissue", verified = TRUE),
    px_bk              = list(analyte = "bradykinin (endothelial proximal space)", units = "nM", specimen = "tissue", verified = TRUE),
    px_lana_kal_hmwk   = list(analyte = "lanadelumab-kallikrein-HMWK complex (endothelial proximal space)", units = "nM", specimen = "tissue", verified = TRUE),

    gc1qr                = list(analyte = "gC1q receptor complex (endothelial cell surface)", units = "number/cell", specimen = "tissue", verified = TRUE),
    bdkrb2               = list(analyte = "bradykinin B2 receptor (endothelial cell surface)", units = "number/cell", specimen = "tissue", verified = TRUE),
    fxii_gc1qr           = list(analyte = "FXII-gC1qR complex (endothelial cell surface)", units = "number/cell", specimen = "tissue", verified = TRUE),
    fxiia_gc1qr          = list(analyte = "FXIIa-gC1qR complex (endothelial cell surface)", units = "number/cell", specimen = "tissue", verified = TRUE),
    prekal_hmwk_gc1qr    = list(analyte = "prekallikrein-HMWK-gC1qR complex (endothelial cell surface)", units = "number/cell", specimen = "tissue", verified = TRUE),
    kal_hmwk_gc1qr       = list(analyte = "kallikrein-HMWK-gC1qR complex (endothelial cell surface)", units = "number/cell", specimen = "tissue", verified = TRUE),
    kal_chmwk_gc1qr      = list(analyte = "kallikrein-cHMWK-gC1qR complex (endothelial cell surface)", units = "number/cell", specimen = "tissue", verified = TRUE),
    c1inh_fxiia_gc1qr    = list(analyte = "C1-INH-FXIIa-gC1qR complex (endothelial cell surface)", units = "number/cell", specimen = "tissue", verified = TRUE),
    c1inh_kal_hmwk_gc1qr = list(analyte = "C1-INH-kallikrein-HMWK-gC1qR complex (endothelial cell surface)", units = "number/cell", specimen = "tissue", verified = TRUE),
    lana_kal_hmwk_gc1qr  = list(analyte = "lanadelumab-kallikrein-HMWK-gC1qR complex (endothelial cell surface)", units = "number/cell", specimen = "tissue", verified = TRUE),
    bk_bdkrb2            = list(analyte = "bradykinin-B2 receptor complex (endothelial cell surface)", units = "number/cell", specimen = "tissue", verified = TRUE),

    prekal_deg               = list(analyte = "degraded prekallikrein", units = "nM", specimen = "plasma", verified = TRUE),
    kal_deg                  = list(analyte = "degraded kallikrein", units = "nM", specimen = "plasma", verified = TRUE),
    hmwk_deg                 = list(analyte = "degraded HMWK", units = "nM", specimen = "plasma", verified = TRUE),
    chmwk_deg                = list(analyte = "degraded cHMWK", units = "nM", specimen = "plasma", verified = TRUE),
    c1inh_deg                = list(analyte = "degraded C1-INH", units = "nM", specimen = "plasma", verified = TRUE),
    fxii_deg                 = list(analyte = "degraded FXII", units = "nM", specimen = "plasma", verified = TRUE),
    fxiia_deg                = list(analyte = "degraded FXIIa", units = "nM", specimen = "plasma", verified = TRUE),
    bk_deg                   = list(analyte = "degraded bradykinin", units = "nM", specimen = "plasma", verified = TRUE),
    c1inh_kal_deg            = list(analyte = "degraded C1-INH-kallikrein complex", units = "nM", specimen = "plasma", verified = TRUE),
    c1inh_kal_hmwk_deg       = list(analyte = "degraded C1-INH-kallikrein-HMWK complex", units = "nM", specimen = "plasma", verified = TRUE),
    c1inh_fxiia_deg          = list(analyte = "degraded C1-INH-FXIIa complex", units = "nM", specimen = "plasma", verified = TRUE),
    gc1qr_deg                = list(analyte = "degraded gC1q receptor (endothelial cell surface)", units = "number/cell", specimen = "tissue", verified = TRUE),
    bdkrb2_deg               = list(analyte = "degraded bradykinin B2 receptor (endothelial cell surface)", units = "number/cell", specimen = "tissue", verified = TRUE),
    bk_bdkrb2_deg            = list(analyte = "degraded bradykinin-B2 receptor complex (endothelial cell surface)", units = "number/cell", specimen = "tissue", verified = TRUE),
    prekal_hmwk_gc1qr_deg    = list(analyte = "degraded prekallikrein-HMWK-gC1qR complex (endothelial cell surface)", units = "number/cell", specimen = "tissue", verified = TRUE),
    kal_hmwk_gc1qr_deg       = list(analyte = "degraded kallikrein-HMWK-gC1qR complex (endothelial cell surface)", units = "number/cell", specimen = "tissue", verified = TRUE),
    kal_chmwk_gc1qr_deg      = list(analyte = "degraded kallikrein-cHMWK-gC1qR complex (endothelial cell surface)", units = "number/cell", specimen = "tissue", verified = TRUE),
    fxii_gc1qr_deg           = list(analyte = "degraded FXII-gC1qR complex (endothelial cell surface)", units = "number/cell", specimen = "tissue", verified = TRUE),
    fxiia_gc1qr_deg          = list(analyte = "degraded FXIIa-gC1qR complex (endothelial cell surface)", units = "number/cell", specimen = "tissue", verified = TRUE),
    c1inh_fxiia_gc1qr_deg    = list(analyte = "degraded C1-INH-FXIIa-gC1qR complex (endothelial cell surface)", units = "number/cell", specimen = "tissue", verified = TRUE),
    c1inh_kal_hmwk_gc1qr_deg = list(analyte = "degraded C1-INH-kallikrein-HMWK-gC1qR complex (endothelial cell surface)", units = "number/cell", specimen = "tissue", verified = TRUE),
    lana_kal_hmwk_gc1qr_deg  = list(analyte = "degraded lanadelumab-kallikrein-HMWK-gC1qR complex (endothelial cell surface)", units = "number/cell", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NA_character_,
      notes              = paste(
        "Allometric size descriptor on lanadelumab CL/F and V/F, normalised to 70 kg",
        "(MOESM6 virtual-patient sampler line 47:",
        "value = mean * (body_weight / 70)^WT_exponent * exp(eta)).",
        "The virtual population samples WT log-normally with mean 81.1 kg and",
        "CV 28.1%, truncated to 36.7-178 kg (MOESM7 weight_Param_info)."
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    species       = "human",
    n_subjects    = 1000,
    n_studies     = 8,
    disease_state = "Hereditary angioedema (HAE) due to C1-inhibitor deficiency (type I / II), simulated during remission and during acute attacks; healthy-control behaviour is obtained by raising ksyn_c1inh from 11.883 to 39.608 nM/h (Param_HAE.r line 8).",
    weight_range  = "36.7-178 kg (virtual population; log-normal mean 81.1 kg, CV 28.1%)",
    dose_range    = "Lanadelumab 30, 100, 150, 300 and 400 mg SC Q2W or Q4W; fixed-dose C1-INH 1000 U IV twice weekly",
    regions       = "Not stated (multinational phase 1a/1b/3 lanadelumab and C1-INH programmes)",
    notes         = paste(
      "Virtual cohort of 1000 patients differing in lanadelumab PK, attack frequency,",
      "attack timing and attack severity. Baseline attack rate 3.5 attacks/month (HELP,",
      "NCT02586805); attacks are Poisson-timed with a normally distributed severity.",
      "Calibration and validation data come from the eight trials in Supplementary Table S1",
      "(lanadelumab NCT01923207, NCT02093923, NCT02586805, NCT02741596; Cinryze",
      "NCT00289211, NCT01005888, NCT00438815, NCT00462709).",
      "The system parameters are identical for healthy subjects, HAE patients in remission",
      "and HAE patients during an attack, except for the C1-INH synthesis rate."
    )
  )

  ini({
    # ==================================================================
    # Lanadelumab popPK -- one-compartment SC, apparent (CL/F, V/F).
    # Point estimates and IIV variances inherited from the population PK
    # analysis cited as reference [15]; reproduced in MOESM7
    # (SimInput_HAE.r popPK_Param_info, lines 40-42). Not re-estimated here.
    # ==================================================================
    lka <- fixed(log(0.0179)); label("Absorption rate constant after SC dosing (ka, 1/h)")            # MOESM7 line 42 (PK.ka mean 0.0179)
    lcl <- fixed(log(0.0249)); label("Apparent clearance (CL/F, L/h) at 70 kg")                       # MOESM7 line 41 (PK.CL_F mean 0.0249)
    lvc <- fixed(log(12.8));   label("Apparent central volume (V/F, L) at 70 kg")                     # MOESM7 line 40 (PK.V_F mean 12.8)
    e_wt_cl <- fixed(0.891);   label("Power exponent on (WT/70) for CL/F (unitless)")                 # MOESM7 line 41 (WT_exponent)
    e_wt_vc <- fixed(0.717);   label("Power exponent on (WT/70) for V/F (unitless)")                  # MOESM7 line 40 (WT_exponent)
    mw_lana <- fixed(150000);  label("Lanadelumab molecular weight (g/mol)")                          # MOESM7 line 17 (PK_MW)

    # ==================================================================
    # HAE attack driver. An attack multiplies the FXII autoactivation rate
    # by (1 + aa_fold_max * severity) for aa_duration hours, where severity
    # is the normalized fold sampled from N(0.35, 0.1) truncated to
    # (0.2, 1.0]. The reported mean (SD) 4.2-fold (1.2-fold) increase is
    # aa_fold_max * 0.35 (aa_fold_max * 0.1).
    # ==================================================================
    aa_fold_max <- fixed(12);  label("Maximum fold-increase in FXII autoactivation rate during an attack (unitless)")  # Param_HAE.r line 11 (AA_AA_Fold_Max)
    aa_duration <- fixed(12);  label("Duration of elevated FXII autoactivation during an attack (h)")                  # Param_HAE.r line 12 (AA_AA_t_duration)
    bk_threshold <- fixed(0.020); label("Plasma bradykinin threshold defining an HAE attack (nM)")                     # MOESM7 line 31 (BK_threshold); main text 20 pM

    # ==================================================================
    # Physiologic volumes and vascular <-> proximal exchange
    # ==================================================================
    v_proximal <- fixed(8e-15);    label("Proximal-space volume per endothelial cell (L)")                    # Param_HAE.r line 92 (20 nm x 400 um2)
    v_medium <- fixed(1.23e-12);   label("Plasma volume per endothelial cell (L)")                            # Param_HAE.r line 93 (3.126 L / 2.54e12 cells)
    circ_rate <- fixed(0.01);      label("Inverse vascular circulation time (1/s; 100 s circulation time)")   # Param_HAE.r line 94; Table S7 assumption 3
    bk_exchange_ratio <- fixed(30); label("Fold-faster proximal exchange for bradykinin vs other species (unitless)")  # Param_HAE.r line 95 (calibrated)
    avogadro <- fixed(6.02e23);    label("Avogadro constant (number/mol)")                                    # Param_HAE.r line 16

    # ==================================================================
    # Species half-lives (h) -- degradation rate constants are log(2)/t1/2
    # ==================================================================
    thalf_prekal <- fixed(24);          label("Prekallikrein half-life (h)")                                  # Param_HAE.r line 24 (Labcorp.com)
    thalf_kal <- fixed(0.0833333);      label("Kallikrein half-life (h; 5 min)")                              # Param_HAE.r line 25 (Cumming 1984)
    thalf_hmwk <- fixed(157);           label("HMWK half-life (h)")                                           # Param_HAE.r line 26 (Weidmann 2017)
    thalf_chmwk <- fixed(11.2142857);   label("cHMWK half-life (h; calibrated as 157/14)")                     # Param_HAE.r line 27 (calibrated)
    thalf_c1inh <- fixed(42);           label("C1-INH half-life (h)")                                         # Param_HAE.r line 28 (Weidmann 2017)
    thalf_bound_c1inh <- fixed(0.05);   label("Inhibitor-bound C1-INH complex half-life (h; 3 min)")           # Param_HAE.r line 29 (assumed as A2M, Jensen 2004)
    thalf_fxii <- fixed(60);            label("FXII half-life (h)")                                           # Param_HAE.r line 30 (Weidmann 2017)
    thalf_fxiia <- fixed(0.0833333);    label("FXIIa half-life (h; 5 min, assumed same as kallikrein)")        # Param_HAE.r line 31
    thalf_bk <- fixed(0.0125);          label("Bradykinin half-life (h; 45 s)")                               # Param_HAE.r line 32 (Cyr 2001)
    thalf_gc1qr <- fixed(2);            label("gC1q receptor half-life (h)")                                  # Param_HAE.r line 33 (typical receptor turnover)
    thalf_bdkrb2 <- fixed(2);           label("Bradykinin B2 receptor half-life (h)")                          # Param_HAE.r line 34 (typical receptor turnover)

    # ==================================================================
    # Synthesis rates
    # ==================================================================
    ksyn_c1inh <- fixed(11.883);    label("C1-INH synthesis rate for HAE patients (nM/h); healthy controls use 39.608")  # Param_HAE.r line 8
    ksyn_prekal <- fixed(41.58883); label("Prekallikrein synthesis rate (nM/h)")                               # Param_HAE.r line 54 (calibrated)
    ksyn_fxii <- fixed(10.83042);   label("FXII synthesis rate (nM/h)")                                        # Param_HAE.r line 55 (calibrated)
    ksyn_hmwk <- fixed(39.93322);   label("HMWK synthesis rate (nM/h)")                                        # Param_HAE.r line 56 (calibrated)

    # ==================================================================
    # Binding affinities (nM)
    # ==================================================================
    kd_prekal_hmwk <- fixed(12);     label("Kd for prekallikrein-HMWK binding (nM)")                           # Param_HAE.r line 60 (Bock 1985)
    kd_kal_hmwk <- fixed(15);        label("Kd for kallikrein-HMWK binding (nM)")                              # Param_HAE.r line 61 (Bock 1985)
    kd_kal_chmwk <- fixed(72);       label("Kd for kallikrein-cHMWK binding (nM)")                             # Param_HAE.r line 62 (Kenniston 2014)
    kd_c1inh_kal <- fixed(150);      label("Kd for C1-INH-kallikrein binding (nM)")                            # Param_HAE.r line 63 (Harpel 1985)
    kd_c1inh_fxiia <- fixed(1720);   label("Kd for C1-INH-FXIIa binding (nM)")                                 # Param_HAE.r line 64 (Bjorkqvist 2015)
    kd_hmwk_gc1qr <- fixed(10.35);   label("Kd for HMWK-gC1qR binding (nM)")                                   # Param_HAE.r line 65 (Fernando 2003)
    kd_chmwk_gc1qr <- fixed(10.35);  label("Kd for cHMWK-gC1qR binding (nM; assumed same as HMWK)")            # Param_HAE.r line 66
    kd_fxii_gc1qr <- fixed(144);     label("Kd for FXII-gC1qR binding (nM)")                                   # Param_HAE.r line 67 (Reddigari 1993)
    kd_fxiia_gc1qr <- fixed(144);    label("Kd for FXIIa-gC1qR binding (nM; assumed same as FXII)")            # Param_HAE.r line 68
    kd_bk_bdkrb2 <- fixed(0.5);      label("Kd for bradykinin-B2 receptor binding (nM)")                       # Param_HAE.r line 69 (Paquet 1999)
    kd_kal_lana <- fixed(0.12);      label("Kd for lanadelumab-kallikrein binding (nM)")                       # MOESM7 line 15 (Kenniston 2014)

    # ==================================================================
    # Association rate constants, 1/(nM*h)
    # (published as 1/(M*h) in Table S6; the tabulated numbers are already
    # divided by the 1e9 M -> nM scaling factor, so they are 1/(nM*h))
    # ==================================================================
    kon_hmwk_gc1qr <- fixed(0.4428);    label("kon for HMWK-gC1qR binding (1/(nM*h); 1.23e5 1/(M*s))")         # Param_HAE.r line 71 (Pixley 2011)
    kon_complex_ratio <- fixed(10);     label("Fold-higher kon for soluble prekallikrein/kallikrein-HMWK binding vs HMWK-gC1qR (unitless)")  # Param_HAE.r lines 76-78 (lower steric constraint)
    kon_c1inh_kal <- fixed(0.0612);     label("kon for C1-INH-kallikrein binding (1/(nM*h); 1.7e4 1/(M*s))")   # Param_HAE.r line 80 (Drouet 2018)
    kon_c1inh_fxiia <- fixed(0.1332);   label("kon for C1-INH-FXIIa binding (1/(nM*h); 3.7e4 1/(M*s))")        # Param_HAE.r line 81 (Drouet 2018)
    kon_bk_bdkrb2 <- fixed(36);         label("kon for bradykinin-B2 receptor binding (1/(nM*h); 1e7 1/(M*s))") # Param_HAE.r line 82 (assumed fast)
    kon_kal_lana <- fixed(12.096);      label("kon for lanadelumab-kallikrein binding (1/(nM*h); 3.36e6 1/(M*s))")  # MOESM7 line 16 (Sexton 2013)

    # ==================================================================
    # Enzymatic conversion. Note that the published model applies FXII
    # autoactivation as a first-order process in FXII_gC1qR (Table S5 E37/E39
    # and MOESM2 line 370); the tabulated Km_FXII_AutoActivation = 110 nM is
    # declared but never used. See vignette Errata.
    # ==================================================================
    kcat_fxii_autoact <- fixed(0.0475); label("kcat for FXII autoactivation (1/h)")                            # Param_HAE.r line 85 (calibrated)
    km_fxii_act <- fixed(510);          label("Km for kallikrein-mediated FXII activation (nM)")               # Param_HAE.r line 86 (Tankersley 1984)
    kcat_fxii_act <- fixed(15);         label("kcat for kallikrein-mediated FXII activation (1/h)")            # Param_HAE.r line 87 (calibrated)
    km_prekal_act <- fixed(91);         label("Km for FXIIa-mediated prekallikrein activation (nM)")           # Param_HAE.r line 88 (Tankersley 1984)
    kcat_prekal_act <- fixed(18);       label("kcat for FXIIa-mediated prekallikrein activation (1/h)")        # Param_HAE.r line 89 (calibrated)
    kcat_hmwk_cleave <- fixed(394.74);  label("kcat for kallikrein-mediated HMWK cleavage (1/h)")              # Param_HAE.r line 90 (calibrated)

    # ==================================================================
    # Initial conditions / receptor densities
    # ==================================================================
    init_fxii <- fixed(375);      label("Initial plasma FXII concentration (nM)")                              # Init_HAE.r line 10
    init_prekal <- fixed(450);    label("Initial plasma prekallikrein concentration (nM)")                      # Init_HAE.r line 8
    init_hmwk <- fixed(670);      label("Initial plasma HMWK concentration (nM)")                               # Init_HAE.r line 9
    init_c1inh <- fixed(720);     label("Initial plasma C1-INH concentration for HAE patients (nM); healthy controls use 2400")  # Init_HAE.r line 7
    init_gc1qr <- fixed(100000);  label("Initial gC1q receptor density (number/cell)")                          # Init_HAE.r line 11 (Mahdi 2001)
    init_bdkrb2 <- fixed(100000); label("Initial bradykinin B2 receptor density (number/cell)")                 # Init_HAE.r line 12 (Paquet 1999)

    pct_chmwk_basal <- fixed(35.6); label("Basal percent cHMWK used to normalise the %cHMWK readout (%)")       # MOESM7 line 18 (CS_prct_HK2Chain_basal)

    # ==================================================================
    # Between-subject variability on lanadelumab PK (variances on the log
    # scale), inherited from the population PK analysis reference [15]
    # ==================================================================
    etalka ~ fixed(0.58041)   # MOESM7 line 42 (Eta_Cov)
    etalcl ~ fixed(0.09575)   # MOESM7 line 41 (Eta_Cov)
    etalvc ~ fixed(0.08399)   # MOESM7 line 40 (Eta_Cov)

    propSd <- fixed(0.001); label("Residual proportional error on lanadelumab concentration (nominal placeholder; not published)")
  })

  model({
    # ------------------------------------------------------------------
    # Derived rate constants
    # ------------------------------------------------------------------
    kdeg_prekal <- log(2) / thalf_prekal
    kdeg_kal <- log(2) / thalf_kal
    kdeg_hmwk <- log(2) / thalf_hmwk
    kdeg_chmwk <- log(2) / thalf_chmwk
    kdeg_c1inh <- log(2) / thalf_c1inh
    kdeg_bound_c1inh <- log(2) / thalf_bound_c1inh
    kdeg_fxii <- log(2) / thalf_fxii
    kdeg_fxiia <- log(2) / thalf_fxiia
    kdeg_bk <- log(2) / thalf_bk
    kdeg_gc1qr <- log(2) / thalf_gc1qr
    kdeg_bdkrb2 <- log(2) / thalf_bdkrb2

    # Receptor synthesis balances turnover at the initial density
    ksyn_gc1qr <- kdeg_gc1qr * init_gc1qr
    ksyn_bdkrb2 <- kdeg_bdkrb2 * init_bdkrb2

    # Soluble prekallikrein / kallikrein binding to HMWK is 10x faster than
    # HMWK binding to the surface receptor complex (lower steric constraint)
    kon_prekal_hmwk <- kon_complex_ratio * kon_hmwk_gc1qr
    kon_kal_hmwk <- kon_complex_ratio * kon_hmwk_gc1qr
    kon_kal_chmwk <- kon_complex_ratio * kon_hmwk_gc1qr
    kon_chmwk_gc1qr <- kon_hmwk_gc1qr
    kon_fxii_gc1qr <- kon_hmwk_gc1qr
    kon_fxiia_gc1qr <- kon_hmwk_gc1qr

    # koff = Kd * kon for every reversible binding event
    koff_prekal_hmwk <- kd_prekal_hmwk * kon_prekal_hmwk
    koff_kal_hmwk <- kd_kal_hmwk * kon_kal_hmwk
    koff_kal_chmwk <- kd_kal_chmwk * kon_kal_chmwk
    koff_c1inh_kal <- kd_c1inh_kal * kon_c1inh_kal
    koff_c1inh_fxiia <- kd_c1inh_fxiia * kon_c1inh_fxiia
    koff_hmwk_gc1qr <- kd_hmwk_gc1qr * kon_hmwk_gc1qr
    koff_chmwk_gc1qr <- kd_chmwk_gc1qr * kon_chmwk_gc1qr
    koff_fxii_gc1qr <- kd_fxii_gc1qr * kon_fxii_gc1qr
    koff_fxiia_gc1qr <- kd_fxiia_gc1qr * kon_fxiia_gc1qr
    koff_bk_bdkrb2 <- kd_bk_bdkrb2 * kon_bk_bdkrb2
    koff_kal_lana <- kd_kal_lana * kon_kal_lana

    # Vascular <-> proximal exchange. The published distribution clearance is
    # CLD = v_proximal / circulation_time; the C source writes each exchange
    # flux as k12*C_plasma*v_medium - k21*C_prox*v_proximal, which reduces
    # exactly to CLD*(C_plasma - C_prox). Dividing by the respective volume
    # gives the k12 / k21 forms used below.
    cld <- circ_rate * 3600 * v_proximal
    cld_bk <- bk_exchange_ratio * cld
    k12 <- cld / v_medium
    k21 <- cld / v_proximal
    k12_bk <- cld_bk / v_medium
    k21_bk <- cld_bk / v_proximal

    # Converts a surface species (number/cell) to a proximal-space
    # concentration (nM)
    num_to_conc <- 1e9 / (avogadro * v_proximal)

    # ------------------------------------------------------------------
    # Lanadelumab PK (one-compartment, first-order SC absorption)
    # ------------------------------------------------------------------
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc

    Cc <- central / vc                 # ug/mL (= mg/L)
    lana_nm <- central / vc / mw_lana * 1e6   # nM

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl / vc * central

    # ------------------------------------------------------------------
    # HAE attack driver. Each attack is delivered as a dose into `trigger`
    # whose amount is the normalized severity; the fold-increase applies for
    # aa_duration hours after that dose. Before the first attack tad() is
    # missing, which the comparison below resolves to the remission state.
    # ------------------------------------------------------------------
    tad_trigger <- tad(trigger)
    if (tad_trigger < aa_duration) {
      fold_fxii <- 1 + aa_fold_max * podo(trigger)
    } else {
      fold_fxii <- 1
    }
    d/dt(trigger) <- 0

    # ------------------------------------------------------------------
    # Vascular-space reversible binding fluxes (nM/h). MOESM2 Rxn 12-17, 22-23.
    # ------------------------------------------------------------------
    bind_prekal_hmwk <- kon_prekal_hmwk * prekal * hmwk - koff_prekal_hmwk * prekal_hmwk
    bind_kal_hmwk <- kon_kal_hmwk * kal * hmwk - koff_kal_hmwk * kal_hmwk
    bind_c1inh_kal <- kon_c1inh_kal * c1inh * kal - koff_c1inh_kal * c1inh_kal
    bind_c1inh_kal_hmwk <- kon_c1inh_kal * c1inh * kal_hmwk - koff_c1inh_kal * c1inh_kal_hmwk
    bind_c1inh_fxiia <- kon_c1inh_fxiia * c1inh * fxiia - koff_c1inh_fxiia * c1inh_fxiia
    bind_kal_chmwk <- kon_kal_chmwk * kal * chmwk - koff_kal_chmwk * kal_chmwk
    bind_lana_kal <- kon_kal_lana * kal * lana_nm - koff_kal_lana * lana_kal
    bind_lana_kal_hmwk <- kon_kal_lana * kal_hmwk * lana_nm - koff_kal_lana * lana_kal_hmwk

    # ------------------------------------------------------------------
    # Endothelial-surface fluxes (number/cell/h). MOESM2 surface Rxn 1-29.
    # ------------------------------------------------------------------
    sf_fxii_bind <- kon_fxii_gc1qr * px_fxii * gc1qr - koff_fxii_gc1qr * fxii_gc1qr
    sf_fxiia_bind <- kon_fxiia_gc1qr * px_fxiia * gc1qr - koff_fxiia_gc1qr * fxiia_gc1qr
    sf_prekal_hmwk_bind <- kon_hmwk_gc1qr * px_prekal_hmwk * gc1qr - koff_hmwk_gc1qr * prekal_hmwk_gc1qr
    sf_kal_hmwk_bind <- kon_hmwk_gc1qr * px_kal_hmwk * gc1qr - koff_hmwk_gc1qr * kal_hmwk_gc1qr
    sf_kal_chmwk_bind <- kon_chmwk_gc1qr * px_kal_chmwk * gc1qr - koff_chmwk_gc1qr * kal_chmwk_gc1qr
    sf_c1inh_fxiia_bind <- kon_fxiia_gc1qr * px_c1inh_fxiia * gc1qr - koff_fxiia_gc1qr * c1inh_fxiia_gc1qr
    sf_c1inh_kal_hmwk_bind <- kon_hmwk_gc1qr * px_c1inh_kal_hmwk * gc1qr - koff_hmwk_gc1qr * c1inh_kal_hmwk_gc1qr
    sf_lana_kal_hmwk_bind <- kon_hmwk_gc1qr * px_lana_kal_hmwk * gc1qr - koff_hmwk_gc1qr * lana_kal_hmwk_gc1qr
    sf_bk_bdkrb2 <- kon_bk_bdkrb2 * px_bk * bdkrb2 - koff_bk_bdkrb2 * bk_bdkrb2

    # C1-INH inhibition on the surface (Table S7 assumption 6)
    sf_c1inh_fxiia_gc1qr <- kon_c1inh_fxiia * px_c1inh * fxiia_gc1qr - koff_c1inh_fxiia * c1inh_fxiia_gc1qr
    sf_c1inh_kal_hmwk_gc1qr <- kon_c1inh_kal * px_c1inh * kal_hmwk_gc1qr - koff_c1inh_kal * c1inh_kal_hmwk_gc1qr

    # Lanadelumab inhibition of surface-bound kallikrein (Table S7 assumption 6)
    sf_lana_bind_surface <- kon_kal_lana * kal_hmwk_gc1qr * lana_nm - koff_kal_lana * lana_kal_hmwk_gc1qr

    # Enzymatic conversions. Autoactivation is first-order in FXII_gC1qR
    # (Table S5 E37/E39); the two facilitated activations are Michaelis-Menten
    # in the substrate expressed as a proximal-space concentration.
    sf_autoact <- fold_fxii * kcat_fxii_autoact * fxii_gc1qr
    sf_fxii_act <- kcat_fxii_act * kal_hmwk_gc1qr * fxii_gc1qr * num_to_conc /
      (km_fxii_act + fxii_gc1qr * num_to_conc)
    sf_prekal_act <- kcat_prekal_act * fxiia_gc1qr * prekal_hmwk_gc1qr * num_to_conc /
      (km_prekal_act + prekal_hmwk_gc1qr * num_to_conc)
    sf_hmwk_cleave <- kcat_hmwk_cleave * kal_hmwk_gc1qr

    # ------------------------------------------------------------------
    # Vascular (plasma) space -- nM. MOESM2 diff_*_in_plasma.
    # Exogenous C1-INH is administered as a dose into `c1inh`, reproducing
    # the flux_C1Inh_inj_nmol_per_hr term of Table S5 E6.
    # ------------------------------------------------------------------
    d/dt(fxii) <- ksyn_fxii - kdeg_fxii * fxii - k12 * (fxii - px_fxii)
    d/dt(fxiia) <- -kdeg_fxiia * fxiia - bind_c1inh_fxiia - k12 * (fxiia - px_fxiia)
    d/dt(c1inh) <- ksyn_c1inh - kdeg_c1inh * c1inh - bind_c1inh_kal -
      bind_c1inh_kal_hmwk - bind_c1inh_fxiia - k12 * (c1inh - px_c1inh)
    d/dt(hmwk) <- ksyn_hmwk - kdeg_hmwk * hmwk - bind_prekal_hmwk - bind_kal_hmwk
    d/dt(prekal) <- ksyn_prekal - kdeg_prekal * prekal - bind_prekal_hmwk
    d/dt(kal) <- -kdeg_kal * kal - bind_kal_hmwk - bind_c1inh_kal - bind_kal_chmwk - bind_lana_kal
    d/dt(chmwk) <- -bind_kal_chmwk - kdeg_chmwk * chmwk
    d/dt(bk) <- -kdeg_bk * bk - k12_bk * (bk - px_bk)
    d/dt(prekal_hmwk) <- bind_prekal_hmwk - k12 * (prekal_hmwk - px_prekal_hmwk)
    d/dt(kal_hmwk) <- bind_kal_hmwk - bind_c1inh_kal_hmwk - bind_lana_kal_hmwk -
      k12 * (kal_hmwk - px_kal_hmwk)
    d/dt(kal_chmwk) <- bind_kal_chmwk - k12 * (kal_chmwk - px_kal_chmwk)
    d/dt(c1inh_fxiia) <- bind_c1inh_fxiia - kdeg_bound_c1inh * c1inh_fxiia -
      k12 * (c1inh_fxiia - px_c1inh_fxiia)
    d/dt(c1inh_kal) <- bind_c1inh_kal - kdeg_bound_c1inh * c1inh_kal
    d/dt(c1inh_kal_hmwk) <- bind_c1inh_kal_hmwk - kdeg_bound_c1inh * c1inh_kal_hmwk -
      k12 * (c1inh_kal_hmwk - px_c1inh_kal_hmwk)
    d/dt(lana_kal) <- bind_lana_kal
    d/dt(lana_kal_hmwk) <- bind_lana_kal_hmwk - k12 * (lana_kal_hmwk - px_lana_kal_hmwk)

    # ------------------------------------------------------------------
    # Endothelial proximal space -- nM. MOESM2 diff_CS_<species>.
    # ------------------------------------------------------------------
    d/dt(px_fxii) <- k21 * (fxii - px_fxii) - sf_fxii_bind * num_to_conc
    d/dt(px_fxiia) <- k21 * (fxiia - px_fxiia) - sf_fxiia_bind * num_to_conc
    d/dt(px_c1inh) <- k21 * (c1inh - px_c1inh) -
      (sf_c1inh_fxiia_gc1qr + sf_c1inh_kal_hmwk_gc1qr) * num_to_conc
    d/dt(px_prekal_hmwk) <- k21 * (prekal_hmwk - px_prekal_hmwk) - sf_prekal_hmwk_bind * num_to_conc
    d/dt(px_kal_hmwk) <- k21 * (kal_hmwk - px_kal_hmwk) - sf_kal_hmwk_bind * num_to_conc
    d/dt(px_c1inh_fxiia) <- k21 * (c1inh_fxiia - px_c1inh_fxiia) - sf_c1inh_fxiia_bind * num_to_conc
    d/dt(px_c1inh_kal_hmwk) <- k21 * (c1inh_kal_hmwk - px_c1inh_kal_hmwk) -
      sf_c1inh_kal_hmwk_bind * num_to_conc
    d/dt(px_kal_chmwk) <- k21 * (kal_chmwk - px_kal_chmwk) - sf_kal_chmwk_bind * num_to_conc
    d/dt(px_bk) <- k21_bk * (bk - px_bk) + (sf_hmwk_cleave - sf_bk_bdkrb2) * num_to_conc
    d/dt(px_lana_kal_hmwk) <- k21 * (lana_kal_hmwk - px_lana_kal_hmwk) -
      sf_lana_kal_hmwk_bind * num_to_conc

    # ------------------------------------------------------------------
    # Endothelial cell surface -- number/cell. MOESM2 diff_CS_*_gC1qR.
    # ------------------------------------------------------------------
    d/dt(gc1qr) <- ksyn_gc1qr - kdeg_gc1qr * gc1qr - sf_fxii_bind - sf_prekal_hmwk_bind -
      sf_fxiia_bind - sf_kal_hmwk_bind - sf_kal_chmwk_bind - sf_c1inh_fxiia_bind -
      sf_c1inh_kal_hmwk_bind - sf_lana_kal_hmwk_bind
    d/dt(bdkrb2) <- ksyn_bdkrb2 - kdeg_bdkrb2 * bdkrb2 - sf_bk_bdkrb2
    d/dt(fxii_gc1qr) <- sf_fxii_bind - sf_autoact - sf_fxii_act - kdeg_gc1qr * fxii_gc1qr
    d/dt(fxiia_gc1qr) <- sf_autoact + sf_fxii_act + sf_fxiia_bind - sf_c1inh_fxiia_gc1qr -
      kdeg_gc1qr * fxiia_gc1qr
    d/dt(prekal_hmwk_gc1qr) <- sf_prekal_hmwk_bind - sf_prekal_act - kdeg_gc1qr * prekal_hmwk_gc1qr
    d/dt(kal_hmwk_gc1qr) <- sf_prekal_act - sf_hmwk_cleave + sf_kal_hmwk_bind -
      sf_c1inh_kal_hmwk_gc1qr - kdeg_gc1qr * kal_hmwk_gc1qr - sf_lana_bind_surface
    d/dt(kal_chmwk_gc1qr) <- sf_hmwk_cleave + sf_kal_chmwk_bind - kdeg_gc1qr * kal_chmwk_gc1qr
    d/dt(c1inh_fxiia_gc1qr) <- sf_c1inh_fxiia_gc1qr + sf_c1inh_fxiia_bind -
      kdeg_gc1qr * c1inh_fxiia_gc1qr
    d/dt(c1inh_kal_hmwk_gc1qr) <- sf_c1inh_kal_hmwk_gc1qr + sf_c1inh_kal_hmwk_bind -
      kdeg_gc1qr * c1inh_kal_hmwk_gc1qr
    d/dt(lana_kal_hmwk_gc1qr) <- sf_lana_bind_surface + sf_lana_kal_hmwk_bind -
      kdeg_gc1qr * lana_kal_hmwk_gc1qr
    d/dt(bk_bdkrb2) <- sf_bk_bdkrb2 - kdeg_bdkrb2 * bk_bdkrb2

    # ------------------------------------------------------------------
    # Degradation sinks (Table S4). These carry no feedback; they exist so
    # that cumulative degraded mass can be audited against synthesis.
    # ------------------------------------------------------------------
    d/dt(prekal_deg) <- kdeg_prekal * prekal
    d/dt(kal_deg) <- kdeg_kal * kal
    d/dt(hmwk_deg) <- kdeg_hmwk * hmwk
    d/dt(chmwk_deg) <- kdeg_chmwk * chmwk
    d/dt(c1inh_deg) <- kdeg_c1inh * c1inh
    d/dt(fxii_deg) <- kdeg_fxii * fxii
    d/dt(fxiia_deg) <- kdeg_fxiia * fxiia
    d/dt(bk_deg) <- kdeg_bk * bk
    d/dt(c1inh_kal_deg) <- kdeg_bound_c1inh * c1inh_kal
    d/dt(c1inh_kal_hmwk_deg) <- kdeg_bound_c1inh * c1inh_kal_hmwk
    d/dt(c1inh_fxiia_deg) <- kdeg_bound_c1inh * c1inh_fxiia
    d/dt(gc1qr_deg) <- kdeg_gc1qr * gc1qr
    d/dt(bdkrb2_deg) <- kdeg_bdkrb2 * bdkrb2
    d/dt(bk_bdkrb2_deg) <- kdeg_bdkrb2 * bk_bdkrb2
    d/dt(prekal_hmwk_gc1qr_deg) <- kdeg_gc1qr * prekal_hmwk_gc1qr
    d/dt(kal_hmwk_gc1qr_deg) <- kdeg_gc1qr * kal_hmwk_gc1qr
    d/dt(kal_chmwk_gc1qr_deg) <- kdeg_gc1qr * kal_chmwk_gc1qr
    d/dt(fxii_gc1qr_deg) <- kdeg_gc1qr * fxii_gc1qr
    d/dt(fxiia_gc1qr_deg) <- kdeg_gc1qr * fxiia_gc1qr
    d/dt(c1inh_fxiia_gc1qr_deg) <- kdeg_gc1qr * c1inh_fxiia_gc1qr
    d/dt(c1inh_kal_hmwk_gc1qr_deg) <- kdeg_gc1qr * c1inh_kal_hmwk_gc1qr
    d/dt(lana_kal_hmwk_gc1qr_deg) <- kdeg_gc1qr * lana_kal_hmwk_gc1qr

    # ------------------------------------------------------------------
    # Initial conditions (Init_HAE.r). All other states start at zero.
    # ------------------------------------------------------------------
    fxii(0) <- init_fxii
    prekal(0) <- init_prekal
    hmwk(0) <- init_hmwk
    c1inh(0) <- init_c1inh
    gc1qr(0) <- init_gc1qr
    bdkrb2(0) <- init_bdkrb2

    # ------------------------------------------------------------------
    # Readouts (MOESM2 yout[0..9])
    # ------------------------------------------------------------------
    bk_pm <- bk * 1000                                              # plasma bradykinin in pM
    total_prekal <- prekal + prekal_hmwk                            # nM
    total_hmwk <- hmwk + prekal_hmwk + kal_hmwk + c1inh_kal_hmwk    # nM
    pct_free_prekal <- prekal / (prekal + prekal_hmwk) * 100        # %
    pct_chmwk <- chmwk / (chmwk + hmwk) * 100                       # %
    pct_chmwk_norm <- chmwk / (chmwk + hmwk) * 100 / pct_chmwk_basal * 100
    ro_bdkrb2 <- bk_bdkrb2 / (bk_bdkrb2 + bdkrb2) * 100             # % B2 receptor occupancy
    attack_on <- (bk > bk_threshold)                                # 1 while an attack is scored

    Cc ~ prop(propSd)
  })
}
