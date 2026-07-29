Selvaggio_2021_mRNA_vaccine_QSP <- function() {
  description <- "QSP. Preclinical (rhesus macaque, calibration). Multicompartmental quantitative systems pharmacology model of mRNA vaccine immunogenicity: kinetics of a lipid nanoparticle (LNP) carrying mRNA at the injection site, uptake and antigen expression by dendritic cells (myeloid + plasmacytoid), monocytes and neutrophils, migration to the draining lymph node, T-cell and B-cell adaptive immune response, and antibody production. Adapted and extended from the Chen 2014 immunogenicity backbone (adaptive immune response retained unvaried; Selvaggio 2021 supplement 'Model equations'). Innate-cell recruitment / uptake / migration parameters fitted to Liang 2017 rhesus macaque data (fluorescent-protein-encoding modified mRNA-LNP). Time in days; LNP in number of molecules; cell states in cell counts; antibody in pmole. Delivered as a mean-field extraction with a single mean-affinity B-cell subgroup (Chen 2014 encodes 17 affinity-subclones for the BCR distribution; the single-subgroup extraction with mean Ka is a documented deviation, see vignette Assumptions) and an explicit low-antibody closure for the receptor-competition equation (paper Eq for total antigen requires a Newton solve on Ag_f per timestep; not natively expressible in rxode2)."
  reference <- "Selvaggio G, Leonardelli L, Lofano G, Fresnay S, Parolo S, Medini D, Siena E, Marchetti L. A quantitative systems pharmacology approach to support mRNA vaccine development and optimization. CPT Pharmacometrics Syst Pharmacol. 2021;10(12):1448-1451. doi:10.1002/psp4.12721. Adaptive-immune backbone from Chen X, Hickling T, Vicini P. A mechanistic, multiscale mathematical model of immunogenicity for therapeutic proteins: Part 1-theoretical model. CPT Pharmacometrics Syst Pharmacol. 2014;3:e133. doi:10.1038/psp.2014.30. Innate-cell parameters fitted to macaque data from Liang F, et al. Efficient targeting and activation of antigen-presenting cells in vivo after modified mRNA vaccine administration in rhesus macaques. Mol Ther. 2017;25:2635-2647. doi:10.1016/j.ymthe.2017.08.006."
  vignette <- "Selvaggio_2021_mRNA_vaccine_QSP"

  paper_specific_compartments <- c(
    "lnp",
    "mdc_is", "mdc_is_lnp", "mdc_is_ag", "mdc_is_agl", "mdc_is_agm", "mdc_is_agh",
    "pdc_is", "pdc_is_lnp", "pdc_is_ag", "pdc_is_agl", "pdc_is_agm", "pdc_is_agh",
    "mc_is",  "mc_is_lnp",  "mc_is_ag",
    "np_is",  "np_is_lnp",  "np_is_ag",
    "mdc_ln_lnp", "mdc_ln_ag", "mdc_ln_agl", "mdc_ln_agm", "mdc_ln_agh",
    "pdc_ln_lnp", "pdc_ln_ag", "pdc_ln_agl", "pdc_ln_agm", "pdc_ln_agh",
    "mc_ln_lnp",  "mc_ln_ag",
    "np_ln_lnp",  "np_ln_ag",
    "nt", "ant", "mta", "mt", "ft",
    "nb", "anb", "amb", "mb", "sp", "lp", "ab"
  )

  units <- list(time = "day", dosing = "number of LNP molecules", concentration = "pmole antibody (Ab)")

  covariateData <- list()

  population <- list(
    species        = "rhesus macaque (Macaca mulatta) - calibration cohort; model framework intended for human mRNA vaccine QSP simulation",
    n_subjects     = NA_integer_,
    n_studies      = 1,
    disease_state  = "Naive; single 50 ug modified mRNA-LNP intramuscular injection (fluorescent-protein-encoding mRNA-LNP)",
    dose_range     = "50 ug single-dose mRNA-LNP (Selvaggio 2021 sensitivity analysis Table 1); Liang 2017 macaque data with fluorescent-protein-encoding modified mRNA-LNP",
    regions        = "Italy / United States (COSBI + GlaxoSmithKline Vaccines collaboration)",
    notes          = "Innate-cell kinetics fitted to Liang 2017 macaque data (see Selvaggio 2021 supplement Figure S1 - migration from injection site to lymph node, model fit vs experimental cell counts). Adaptive immune response parameters (T/B cells, antibody kinetics) carried unvaried from Chen 2014 human therapeutic-protein immunogenicity model. Selvaggio 2021 states the model was developed as a QSP framework to support mRNA vaccine development; it is a mean-field deterministic mechanistic model without population-level IIV or residual error."
  )

  ini({
    # NOTE ON PARAMETER STYLE. Every value is fixed at the paper's reported
    # estimate (Selvaggio 2021 Table S1; adaptive-immune values inherited
    # from Chen 2014 as noted). The nlmixr2 UI parser rejects compound
    # expressions containing two or more ini-level parameters, so at the top
    # of model() every fixed parameter is aliased to a local variable
    # (suffix "_v"); the ODEs then use only those local aliases. This is a
    # mechanical restatement, not a change to the model.

    # ==========================================================================
    # LNP (vaccine vector) - Table S1
    # ==========================================================================
    k_deg_lnp <- fixed(3.7697);     label("LNP tissue degradation rate constant (1/day)")                   # Table S1, fitted

    # ==========================================================================
    # Myeloid dendritic cells (mDC) - Table S1
    # ==========================================================================
    mdc0        <- fixed(4.2881e3); label("Initial number of mature mDC at IS (#cells)")                    # Table S1, ref [2]
    k_dt_mdc    <- fixed(0.02310);  label("Death rate of mature mDC (1/day)")                               # Table S1, ref [1]
    k_rc_mdc    <- fixed(3.4865e-12); label("Recruitment rate of mDC as function of LNP (1/(day*#LNP))")    # Table S1, fitted
    k_up_mdc    <- fixed(3.0466e-17); label("mDC LNP uptake rate (1/(day*#LNP))")                           # Table S1, fitted
    k_exp_mdc   <- fixed(16.0766);  label("mDC antigen expression rate (1/day)")                            # Table S1, fitted
    k_is2ln_mdc_lnp <- fixed(3.1487);  label("Migration rate IS to LN, mDC with LNP (1/day)")               # Table S1, fitted
    k_is2ln_mdc_ag  <- fixed(0.1754);  label("Migration rate IS to LN, mDC with Ag (1/day)")                # Table S1, fitted
    k_ln2bl_mdc_ag  <- fixed(29.4478); label("Migration rate LN to BL, mDC with Ag (1/day)")                # Table S1, fitted
    k_is2bl_mdc     <- fixed(0.6603);  label("Migration rate IS to BL, mDC (1/day)")                        # Table S1, fitted

    # ==========================================================================
    # Plasmacytoid dendritic cells (pDC) - Table S1
    # ==========================================================================
    # pDC baseline initial number is zero (Table S1); hard-coded in model()
    # since a zero-value ini() parameter would break the log/exp workflow
    # used for structural rates.
    k_dt_pdc    <- fixed(0.02310);  label("Death rate of mature pDC (1/day)")                               # Table S1, ref [1]
    k_rc_pdc    <- fixed(7.4591e-14); label("Recruitment rate of pDC as function of LNP (1/(day*#LNP))")    # Table S1, fitted
    k_up_pdc    <- fixed(7.9513e-18); label("pDC LNP uptake rate (1/(day*#LNP))")                           # Table S1, fitted
    k_exp_pdc   <- fixed(5.0280);   label("pDC antigen expression rate (1/day)")                            # Table S1, fitted
    k_is2ln_pdc_lnp <- fixed(1.4633);  label("Migration rate IS to LN, pDC with LNP (1/day)")               # Table S1, fitted
    k_is2ln_pdc_ag  <- fixed(0.3504);  label("Migration rate IS to LN, pDC with Ag (1/day)")                # Table S1, fitted
    k_ln2bl_pdc_ag  <- fixed(0.2567);  label("Migration rate LN to BL, pDC with Ag (1/day)")                # Table S1, fitted
    k_is2bl_pdc     <- fixed(0.02605); label("Migration rate IS to BL, pDC (1/day)")                        # Table S1, fitted (0.02605; the printed '02605' reads as 0.02605 by context)

    # ==========================================================================
    # DC MHC and antigen-transport constants - Table S1
    # ==========================================================================
    mhc0     <- fixed(95.1e3);      label("Number of MHC molecules per mature DC (#molecules)")             # Table S1, ref [1]
    k_tr     <- fixed(28.8);        label("Antigen+MHC transport rate cytoplasm to membrane (1/day)")       # Table S1, ref [1]

    # ==========================================================================
    # Monocytes (MC) - Table S1
    # ==========================================================================
    mc0         <- fixed(297.9);    label("Initial number of monocytes at IS (#cells)")                     # Table S1, ref [2]
    k_dt_mc     <- fixed(0.6931);   label("Death rate of monocytes (1/day)")                                # Table S1, ref [1]
    k_rc_mc     <- fixed(3.6350e-11); label("Recruitment rate of MC as function of LNP (1/(day*#LNP))")     # Table S1, fitted
    k_up_mc     <- fixed(2.9225e-17); label("MC LNP uptake rate (1/(day*#LNP))")                            # Table S1, fitted
    k_exp_mc    <- fixed(6.8631);   label("MC antigen expression rate (1/day)")                             # Table S1, fitted
    k_is2ln_mc_lnp <- fixed(1.4867);  label("Migration rate IS to LN, MC with LNP (1/day)")                 # Table S1, fitted
    k_is2ln_mc_ag  <- fixed(0.9636);  label("Migration rate IS to LN, MC with Ag (1/day)")                  # Table S1, fitted
    k_ln2bl_mc_ag  <- fixed(27.8383); label("Migration rate LN to BL, MC with Ag (1/day)")                  # Table S1, fitted
    k_is2bl_mc     <- fixed(2.5318e-5); label("Migration rate IS to BL, MC (1/day)")                        # Table S1, fitted

    # ==========================================================================
    # Neutrophils (NP) - Table S1
    # ==========================================================================
    # NP baseline initial number is zero (Table S1); hard-coded in model().
    k_dt_np     <- fixed(2.3765);   label("Death rate of neutrophils (1/day)")                              # Table S1, ref [1]
    k_rc_np     <- fixed(8.2702e-11); label("Recruitment rate of NP as function of LNP (1/(day*#LNP))")     # Table S1, fitted
    k_up_np     <- fixed(1.5771e-17); label("NP LNP uptake rate (1/(day*#LNP))")                            # Table S1, fitted
    k_exp_np    <- fixed(2.1845);   label("NP antigen expression rate (1/day)")                             # Table S1, fitted
    k_is2ln_np_lnp <- fixed(25.0563); label("Migration rate IS to LN, NP with LNP (1/day)")                 # Table S1, fitted
    k_is2ln_np_ag  <- fixed(0.0277);  label("Migration rate IS to LN, NP with Ag (1/day)")                  # Table S1, fitted
    k_ln2bl_np_ag  <- fixed(734.3855); label("Migration rate LN to BL, NP with Ag (1/day)")                 # Table S1, fitted
    k_is2bl_np     <- fixed(2.14e-5); label("Migration rate IS to BL, NP (1/day)")                          # Table S1, fitted

    # ==========================================================================
    # T cells - Table S1
    # ==========================================================================
    nt0        <- fixed(1.445e3);   label("Initial number of naive T cells in LN (#cells)")                 # Table S1, ref [1]
    k_dt_nt    <- fixed(0.0029);    label("Death rate of naive T cells (1/day)")                            # Table S1, ref [1]
    k_dt_at    <- fixed(0.18);      label("Death rate of active T cells (1/day)")                           # Table S1, ref [1]
    k_dt_mt    <- fixed(2.7937e-4); label("Death rate of memory T cells (1/day)")                           # Table S1, ref [1]
    k_dt_ft    <- fixed(0.18);      label("Death rate of functional T cells (1/day)")                       # Table S1, ref [1]
    k_act_nt   <- fixed(1.5);       label("Maximum activation rate of naive T cells (1/day)")               # Table S1, ref [1]
    k_act_mt   <- fixed(1.5);       label("Maximum activation rate of memory T cells (1/day)")              # Table S1, ref [1]
    k_prol_at  <- fixed(0.5793);    label("Maximum proliferation rate of active T cells (1/day)")           # Table S1, ref [1]
    f1_t       <- fixed(0.5);       label("Differentiation fraction, active T to memory T (unitless)")      # Table S1, ref [1]

    # ==========================================================================
    # B cells (single-subgroup mean-affinity extraction) - Table S1
    # ==========================================================================
    # Selvaggio 2021 encodes 17 B-cell affinity subclones (J=17) with Ka_j
    # spanning 3.91e-1 to 2.56e-4 pM^-1 (per Chen 2014); this extraction uses
    # a single subgroup with geometric-mean Ka to keep the ODE system
    # tractable in rxode2 (17-clone dynamics would multiply the B-cell block
    # by 17x and require an implicit Newton solve for the receptor-
    # competition equation Ag = Ag_f * [1 + 2 * sum(K_a_i * Ab_i / (1 +
    # K_a_i * Ag_f)) + sum(K_a_i * BCR_i / (1 + K_a_i * Ag_f))], which is
    # not natively expressible in rxode2). See vignette Assumptions.
    nb0        <- fixed(5200);      label("Initial number of naive B cells (single subgroup, #cells)")      # Table S1, ref [1]; per-subgroup value used for the single-subgroup extraction
    brn        <- fixed(75e3);      label("Number of BCR per B cell (#molecules/cell)")                     # Table S1, ref [1]
    ka_mean    <- fixed(0.010005);  label("Geometric-mean association rate constant across 17 Chen 2014 B-cell subclones (1/pM)")  # sqrt(3.91e-1 * 2.56e-4) = 0.010005; Table S1 Ka_j range 3.91e-1 to 2.56e-4
    kr_bcr     <- fixed(1);         label("Occupied BCR number at 50% naive-B activation (unitless)")       # Table S1 K_R, ref [1]
    cc_n       <- fixed(10);        label("Carrying capacity for FT to activate a naive B target (unitless)")  # Table S1 CC_N, ref [1]
    cc_m       <- fixed(100);       label("Carrying capacity for FT to activate a memory B target (unitless)") # Table S1 CC_M, ref [1]
    k_prol_anb <- fixed(0.3333);    label("Max proliferation rate, ANB from naive B (1/day)")               # Table S1, ref [1]
    k_prol_amb <- fixed(0.7273);    label("Max proliferation rate, AMB from memory B (1/day)")              # Table S1, ref [1]
    k_dt_ab_b  <- fixed(0.2518);    label("Death rate of active B cells (1/day)")                           # Table S1, ref [1]
    k_dt_mb    <- fixed(7.8278e-5); label("Death rate of memory B cells (1/day)")                           # Table S1, ref [1]
    k_dt_lp    <- fixed(0.0050);    label("Death rate of long-lived plasma cells (1/day)")                  # Table S1, ref [1]
    k_dt_sp    <- fixed(0.2310);    label("Death rate of short-lived plasma cells (1/day)")                 # Table S1, ref [1]
    g1_b       <- fixed(0.5);       label("Fraction of activated B cells differentiating to memory B (unitless)")     # Table S1, ref [1]
    g2_b       <- fixed(0.4);       label("Fraction of activated B cells differentiating to short-lived plasma (unitless)")  # Table S1, ref [1]
    k_act_nb   <- fixed(3);         label("Max activation rate of naive B cells (1/day)")                    # Table S1, ref [1]
    k_act_mb   <- fixed(3);         label("Max activation rate of memory B cells (1/day)")                   # Table S1, ref [1]

    # ==========================================================================
    # Antibody - Table S1
    # ==========================================================================
    k_prod_ab <- fixed(8.64e8);    label("Antibody secretion rate per plasma cell (#molecules/(cell*day))") # Table S1, ref [1]
    k_deg_ab  <- fixed(0.030130435); label("Antibody degradation rate (1/day)")                             # Table S1, ref [1]

    # Lymph node volume (fixed literature default; needed to convert BCR/Ab
    # per-cell counts to lymph-node concentrations for the receptor-competition
    # closure). Not tabulated in Selvaggio 2021 or Chen 2014; a canonical
    # lymph-node reference volume of 0.5 mL is used here to complete the
    # closure. See vignette Assumptions.
    v_ln    <- fixed(5e-4);        label("Lymph node volume (L; fixed literature default)")                 # Assumption

    # Nominal residual placeholder. The paper is a deterministic mechanistic
    # QSP and does not report a residual-error model; this small additive SD
    # is included only to satisfy the nlmixr2 UI observation-declaration
    # requirement (an estimable parameter is needed to bind the mu-reference
    # parser). Users can suppress residual noise via rxode2::zeroRe() in
    # typical-value simulations. See vignette Assumptions.
    addSd_Ab <- 0.01;              label("Nominal additive residual SD on Ab (pmole)")                       # placeholder
  })

  model({
    # ==========================================================================
    # 0. Local aliases (nlmixr2 UI mu-reference workaround)
    # Every ini-level fixed parameter is aliased to a "_v" local before use;
    # subsequent expressions combine only local variables.
    # ==========================================================================
    lnp_mw <- 535       # LNP molecular weight, g/mol (Table S1)
    n_avog <- 6.022e23  # Avogadro constant, #molecules/mol (Table S1)

    k_deg_lnp_v       <- k_deg_lnp
    mdc0_v            <- mdc0
    k_dt_mdc_v        <- k_dt_mdc
    k_rc_mdc_v        <- k_rc_mdc
    k_up_mdc_v        <- k_up_mdc
    k_exp_mdc_v       <- k_exp_mdc
    k_is2ln_mdc_lnp_v <- k_is2ln_mdc_lnp
    k_is2ln_mdc_ag_v  <- k_is2ln_mdc_ag
    k_ln2bl_mdc_ag_v  <- k_ln2bl_mdc_ag
    k_is2bl_mdc_v     <- k_is2bl_mdc

    k_dt_pdc_v        <- k_dt_pdc
    k_rc_pdc_v        <- k_rc_pdc
    k_up_pdc_v        <- k_up_pdc
    k_exp_pdc_v       <- k_exp_pdc
    k_is2ln_pdc_lnp_v <- k_is2ln_pdc_lnp
    k_is2ln_pdc_ag_v  <- k_is2ln_pdc_ag
    k_ln2bl_pdc_ag_v  <- k_ln2bl_pdc_ag
    k_is2bl_pdc_v     <- k_is2bl_pdc

    mhc0_v            <- mhc0
    k_tr_v            <- k_tr

    mc0_v             <- mc0
    k_dt_mc_v         <- k_dt_mc
    k_rc_mc_v         <- k_rc_mc
    k_up_mc_v         <- k_up_mc
    k_exp_mc_v        <- k_exp_mc
    k_is2ln_mc_lnp_v  <- k_is2ln_mc_lnp
    k_is2ln_mc_ag_v   <- k_is2ln_mc_ag
    k_ln2bl_mc_ag_v   <- k_ln2bl_mc_ag
    k_is2bl_mc_v      <- k_is2bl_mc

    k_dt_np_v         <- k_dt_np
    k_rc_np_v         <- k_rc_np
    k_up_np_v         <- k_up_np
    k_exp_np_v        <- k_exp_np
    k_is2ln_np_lnp_v  <- k_is2ln_np_lnp
    k_is2ln_np_ag_v   <- k_is2ln_np_ag
    k_ln2bl_np_ag_v   <- k_ln2bl_np_ag
    k_is2bl_np_v      <- k_is2bl_np

    nt0_v             <- nt0
    k_dt_nt_v         <- k_dt_nt
    k_dt_at_v         <- k_dt_at
    k_dt_mt_v         <- k_dt_mt
    k_dt_ft_v         <- k_dt_ft
    k_act_nt_v        <- k_act_nt
    k_act_mt_v        <- k_act_mt
    k_prol_at_v       <- k_prol_at
    f1_t_v            <- f1_t

    nb0_v             <- nb0
    brn_v             <- brn
    ka_mean_v         <- ka_mean
    kr_bcr_v          <- kr_bcr
    cc_n_v            <- cc_n
    cc_m_v            <- cc_m
    k_prol_anb_v      <- k_prol_anb
    k_prol_amb_v      <- k_prol_amb
    k_dt_ab_b_v       <- k_dt_ab_b
    k_dt_mb_v         <- k_dt_mb
    k_dt_lp_v         <- k_dt_lp
    k_dt_sp_v         <- k_dt_sp
    g1_b_v            <- g1_b
    g2_b_v            <- g2_b
    k_act_nb_v        <- k_act_nb
    k_act_mb_v        <- k_act_mb

    k_prod_ab_v       <- k_prod_ab
    k_deg_ab_v        <- k_deg_ab

    v_ln_v            <- v_ln

    # ==========================================================================
    # 1. Physiological baseline recruitment (steady-state closure)
    # k_br^X = X0 * (k_dt^X + k_IS2BL^X); Selvaggio 2021 Table S1 note.
    # For pDC and NP the baseline count is zero, so k_br is zero.
    # ==========================================================================
    k_br_mdc <- mdc0_v * (k_dt_mdc_v + k_is2bl_mdc_v)
    k_br_pdc <- 0.0
    k_br_mc  <- mc0_v  * (k_dt_mc_v  + k_is2bl_mc_v)
    k_br_np  <- 0.0

    # ==========================================================================
    # 2. LNP kinetics (state in #molecules)
    # Only naive (non-LNP-carrying, non-Ag-carrying) cells at the IS take up
    # LNP - consistent with the paper's individual-cell equations that source
    # uptake from the naive form only.
    # ==========================================================================
    d/dt(lnp) <- -k_up_mdc_v * lnp * mdc_is -
                  k_up_pdc_v * lnp * pdc_is -
                  k_up_mc_v  * lnp * mc_is -
                  k_up_np_v  * lnp * np_is -
                  k_deg_lnp_v * lnp

    # ==========================================================================
    # 3. Myeloid DC compartments (Selvaggio 2021 supplement Model equations)
    # ==========================================================================
    d/dt(mdc_is)     <- +k_br_mdc + k_rc_mdc_v * lnp - k_dt_mdc_v * mdc_is -
                         k_up_mdc_v * lnp * mdc_is - k_is2bl_mdc_v * mdc_is
    d/dt(mdc_is_lnp) <- +k_up_mdc_v * lnp * mdc_is - k_dt_mdc_v * mdc_is_lnp -
                         k_exp_mdc_v * mdc_is_lnp - k_is2ln_mdc_lnp_v * mdc_is_lnp
    d/dt(mdc_is_ag)  <- +k_exp_mdc_v * mdc_is_lnp - k_tr_v * mdc_is_ag -
                         k_dt_mdc_v * mdc_is_ag - k_is2ln_mdc_ag_v * mdc_is_ag
    d/dt(mdc_is_agl) <- +k_tr_v * mdc_is_ag - k_tr_v * mdc_is_agl -
                         k_dt_mdc_v * mdc_is_agl - k_is2ln_mdc_ag_v * mdc_is_agl
    d/dt(mdc_is_agm) <- +k_tr_v * mdc_is_agl - k_tr_v * mdc_is_agm -
                         k_dt_mdc_v * mdc_is_agm - k_is2ln_mdc_ag_v * mdc_is_agm
    d/dt(mdc_is_agh) <- +k_tr_v * mdc_is_agm - k_dt_mdc_v * mdc_is_agh -
                         k_is2ln_mdc_ag_v * mdc_is_agh

    d/dt(mdc_ln_lnp) <- +k_is2ln_mdc_lnp_v * mdc_is_lnp - k_dt_mdc_v * mdc_ln_lnp -
                         k_exp_mdc_v * mdc_ln_lnp
    d/dt(mdc_ln_ag)  <- +k_is2ln_mdc_ag_v * mdc_is_ag - k_dt_mdc_v * mdc_ln_ag +
                         k_exp_mdc_v * mdc_ln_lnp - k_ln2bl_mdc_ag_v * mdc_ln_ag
    d/dt(mdc_ln_agl) <- +k_is2ln_mdc_ag_v * mdc_is_agl - k_dt_mdc_v * mdc_ln_agl +
                         k_tr_v * mdc_ln_ag - k_ln2bl_mdc_ag_v * mdc_ln_agl
    d/dt(mdc_ln_agm) <- +k_is2ln_mdc_ag_v * mdc_is_agm - k_dt_mdc_v * mdc_ln_agm +
                         k_tr_v * mdc_ln_agl - k_ln2bl_mdc_ag_v * mdc_ln_agm
    d/dt(mdc_ln_agh) <- +k_is2ln_mdc_ag_v * mdc_is_agh - k_dt_mdc_v * mdc_ln_agh +
                         k_tr_v * mdc_ln_agm - k_ln2bl_mdc_ag_v * mdc_ln_agh

    # ==========================================================================
    # 4. Plasmacytoid DC compartments (same shape as mDC)
    # ==========================================================================
    d/dt(pdc_is)     <- +k_br_pdc + k_rc_pdc_v * lnp - k_dt_pdc_v * pdc_is -
                         k_up_pdc_v * lnp * pdc_is - k_is2bl_pdc_v * pdc_is
    d/dt(pdc_is_lnp) <- +k_up_pdc_v * lnp * pdc_is - k_dt_pdc_v * pdc_is_lnp -
                         k_exp_pdc_v * pdc_is_lnp - k_is2ln_pdc_lnp_v * pdc_is_lnp
    d/dt(pdc_is_ag)  <- +k_exp_pdc_v * pdc_is_lnp - k_tr_v * pdc_is_ag -
                         k_dt_pdc_v * pdc_is_ag - k_is2ln_pdc_ag_v * pdc_is_ag
    d/dt(pdc_is_agl) <- +k_tr_v * pdc_is_ag - k_tr_v * pdc_is_agl -
                         k_dt_pdc_v * pdc_is_agl - k_is2ln_pdc_ag_v * pdc_is_agl
    d/dt(pdc_is_agm) <- +k_tr_v * pdc_is_agl - k_tr_v * pdc_is_agm -
                         k_dt_pdc_v * pdc_is_agm - k_is2ln_pdc_ag_v * pdc_is_agm
    d/dt(pdc_is_agh) <- +k_tr_v * pdc_is_agm - k_dt_pdc_v * pdc_is_agh -
                         k_is2ln_pdc_ag_v * pdc_is_agh

    d/dt(pdc_ln_lnp) <- +k_is2ln_pdc_lnp_v * pdc_is_lnp - k_dt_pdc_v * pdc_ln_lnp -
                         k_exp_pdc_v * pdc_ln_lnp
    d/dt(pdc_ln_ag)  <- +k_is2ln_pdc_ag_v * pdc_is_ag - k_dt_pdc_v * pdc_ln_ag +
                         k_exp_pdc_v * pdc_ln_lnp - k_ln2bl_pdc_ag_v * pdc_ln_ag
    d/dt(pdc_ln_agl) <- +k_is2ln_pdc_ag_v * pdc_is_agl - k_dt_pdc_v * pdc_ln_agl +
                         k_tr_v * pdc_ln_ag - k_ln2bl_pdc_ag_v * pdc_ln_agl
    d/dt(pdc_ln_agm) <- +k_is2ln_pdc_ag_v * pdc_is_agm - k_dt_pdc_v * pdc_ln_agm +
                         k_tr_v * pdc_ln_agl - k_ln2bl_pdc_ag_v * pdc_ln_agm
    d/dt(pdc_ln_agh) <- +k_is2ln_pdc_ag_v * pdc_is_agh - k_dt_pdc_v * pdc_ln_agh +
                         k_tr_v * pdc_ln_agm - k_ln2bl_pdc_ag_v * pdc_ln_agh

    # ==========================================================================
    # 5. Monocytes (only naive/LNP/Ag states; no membrane-level Ag stratification)
    # ==========================================================================
    d/dt(mc_is)      <- +k_br_mc + k_rc_mc_v * lnp - k_dt_mc_v * mc_is -
                         k_up_mc_v * lnp * mc_is - k_is2bl_mc_v * mc_is
    d/dt(mc_is_lnp)  <- +k_up_mc_v * lnp * mc_is - k_dt_mc_v * mc_is_lnp -
                         k_exp_mc_v * mc_is_lnp - k_is2ln_mc_lnp_v * mc_is_lnp
    d/dt(mc_is_ag)   <- +k_exp_mc_v * mc_is_lnp - k_dt_mc_v * mc_is_ag -
                         k_is2ln_mc_ag_v * mc_is_ag
    d/dt(mc_ln_lnp)  <- +k_is2ln_mc_lnp_v * mc_is_lnp - k_dt_mc_v * mc_ln_lnp -
                         k_exp_mc_v * mc_ln_lnp
    d/dt(mc_ln_ag)   <- +k_is2ln_mc_ag_v * mc_is_ag - k_dt_mc_v * mc_ln_ag +
                         k_exp_mc_v * mc_ln_lnp - k_ln2bl_mc_ag_v * mc_ln_ag

    # ==========================================================================
    # 6. Neutrophils
    # ==========================================================================
    d/dt(np_is)      <- +k_br_np + k_rc_np_v * lnp - k_dt_np_v * np_is -
                         k_up_np_v * lnp * np_is - k_is2bl_np_v * np_is
    d/dt(np_is_lnp)  <- +k_up_np_v * lnp * np_is - k_dt_np_v * np_is_lnp -
                         k_exp_np_v * np_is_lnp - k_is2ln_np_lnp_v * np_is_lnp
    d/dt(np_is_ag)   <- +k_exp_np_v * np_is_lnp - k_dt_np_v * np_is_ag -
                         k_is2ln_np_ag_v * np_is_ag
    d/dt(np_ln_lnp)  <- +k_is2ln_np_lnp_v * np_is_lnp - k_dt_np_v * np_ln_lnp -
                         k_exp_np_v * np_ln_lnp
    d/dt(np_ln_ag)   <- +k_is2ln_np_ag_v * np_is_ag - k_dt_np_v * np_ln_ag +
                         k_exp_np_v * np_ln_lnp - k_ln2bl_np_ag_v * np_ln_ag

    # ==========================================================================
    # 7. Antigen presented in the lymph node
    # Ag = MHC0 * [0.1 * (mDC_LN^AgL + pDC_LN^AgL) + 0.5 * (mDC_LN^AgM +
    # pDC_LN^AgM) + 0.9 * (mDC_LN^AgH + pDC_LN^AgH)] (Selvaggio 2021 supp).
    # ==========================================================================
    ag <- mhc0_v * (0.1 * (mdc_ln_agl + pdc_ln_agl) +
                    0.5 * (mdc_ln_agm + pdc_ln_agm) +
                    0.9 * (mdc_ln_agh + pdc_ln_agh))

    # ==========================================================================
    # 8. Free-antigen closure (explicit low-antibody approximation)
    # Paper Eq: Ag = Ag_f * [1 + 2 * sum_i(K_a_i * Ab_i / (1 + K_a_i * Ag_f)) +
    #                       sum_i(K_a_i * BCR_i / (1 + K_a_i * Ag_f))]
    # This implicit equation for Ag_f given Ag is not natively expressible in
    # rxode2 (would require Newton iteration per timestep). We use the
    # low-antibody limit:
    #   Ag_f ~ Ag / (1 + 2 * K_a * Ab_conc + K_a * BCR_conc)
    # where concentrations are per-lymph-node-volume in pM. See vignette
    # Assumptions.
    # ==========================================================================
    ab_conc  <- ab / (v_ln_v * n_avog) * 1e12
    bcr_tot  <- brn_v * (nb + anb + amb + mb)
    bcr_conc <- bcr_tot / (v_ln_v * n_avog) * 1e12

    denom_agf <- 1 + 2 * ka_mean_v * ab_conc + ka_mean_v * bcr_conc
    ag_f <- ag / denom_agf

    # ==========================================================================
    # 9. Receptor occupancy and B-cell activation functions
    # ==========================================================================
    ro  <- ka_mean_v * ag_f / (1 + ka_mean_v * ag_f)
    r   <- ro * brn_v
    f_b <- r / (r + kr_bcr_v)
    g_b <- (1 - ro) * f_b
    h_b <- (r - kr_bcr_v) / (r + kr_bcr_v)

    # ==========================================================================
    # 10. T-cell distribution / effectiveness functions
    # Small numerical floor (1e-12) added to denominators to avoid division-by-
    # zero when all antigen-presenting states are zero pre-dose.
    # ==========================================================================
    sum_apc_ge_m <- (mdc_ln_agm + mdc_ln_agh) + (pdc_ln_agm + pdc_ln_agh)
    sum_apc_ge_l <- (mdc_ln_agl + mdc_ln_agm + mdc_ln_agh) +
                    (pdc_ln_agl + pdc_ln_agm + pdc_ln_agh)
    t_pool <- nt + ant + mta + mt

    dn <- sum_apc_ge_m / (sum_apc_ge_m + t_pool + 1e-12)
    en <- dn * ((mdc_ln_agh - mdc_ln_agm) + (pdc_ln_agh - pdc_ln_agm)) /
          (sum_apc_ge_m + 1e-12)
    dm <- sum_apc_ge_l / (sum_apc_ge_l + t_pool + 1e-12)
    em <- dm * (((mdc_ln_agm - mdc_ln_agl) + (mdc_ln_agh)) +
                ((pdc_ln_agm - pdc_ln_agl) + (pdc_ln_agh))) /
          (sum_apc_ge_l + 1e-12)

    # ==========================================================================
    # 11. T-cell ODEs (Selvaggio 2021 supplement)
    # State 'mta' corresponds to the paper's AMT (active memory T); renamed
    # from AMT because 'amt' is a reserved rxode2 dose-amount keyword.
    # ==========================================================================
    d/dt(nt)  <- +k_dt_nt_v * nt0_v - k_dt_nt_v * nt - k_act_nt_v * dn * nt
    d/dt(ant) <- +k_act_nt_v * dn * nt - k_dt_at_v * ant - k_prol_at_v * en * ant
    d/dt(mta) <- +k_act_mt_v * dm * mt + k_prol_at_v * em * mta - k_dt_at_v * mta
    d/dt(mt)  <- +k_prol_at_v * (1 - en) * f1_t_v * ant +
                  k_prol_at_v * (1 - em) * f1_t_v * mta -
                  k_dt_mt_v * mt - k_act_mt_v * dm * mt
    d/dt(ft)  <- +k_prol_at_v * (1 - en) * (1 - f1_t_v) * ant +
                  k_prol_at_v * (1 - em) * (1 - f1_t_v) * mta - k_dt_ft_v * ft

    # ==========================================================================
    # 12. B-cell carrying capacities
    # ==========================================================================
    p_n <- (cc_n_v * ft) / (cc_n_v * ft + nb + anb + amb + mb + 1e-12)
    p_m <- (cc_m_v * ft) / (cc_m_v * ft + nb + anb + amb + mb + 1e-12)

    # ==========================================================================
    # 13. B-cell + antibody ODEs (single-subgroup extraction)
    # ==========================================================================
    d/dt(nb)  <- -k_act_nb_v * f_b * p_n * nb
    d/dt(anb) <- +k_act_nb_v * g_b * p_n * nb +
                  k_prol_anb_v * h_b * p_n * anb - k_dt_ab_b_v * anb
    d/dt(amb) <- +k_act_mb_v * g_b * p_m * mb +
                  k_prol_amb_v * h_b * p_m * amb - k_dt_ab_b_v * amb
    d/dt(mb)  <- +k_prol_anb_v * (1 - h_b) * g1_b_v * p_n * anb +
                  k_prol_amb_v * (1 - h_b) * g1_b_v * p_m * amb -
                  k_act_mb_v * f_b * p_m * mb - k_dt_mb_v * mb
    d/dt(sp)  <- +k_prol_anb_v * (1 - h_b) * g2_b_v * p_n * anb +
                  k_prol_amb_v * (1 - h_b) * g2_b_v * p_m * amb - k_dt_sp_v * sp
    d/dt(lp)  <- +k_prol_anb_v * (1 - h_b) * (1 - g1_b_v - g2_b_v) * p_n * anb +
                  k_prol_amb_v * (1 - h_b) * (1 - g1_b_v - g2_b_v) * p_m * amb -
                  k_dt_lp_v * lp
    d/dt(ab)  <- +k_prod_ab_v * (sp + lp) / n_avog * 1e12 - k_deg_ab_v * ab

    # ==========================================================================
    # 14. Initial conditions
    # (pDC and NP baseline counts are zero per Table S1; those states use the
    # default zero initial condition.)
    # ==========================================================================
    lnp(0)     <- 0
    mdc_is(0)  <- mdc0_v
    mc_is(0)   <- mc0_v
    nt(0)      <- nt0_v
    nb(0)      <- nb0_v

    # ==========================================================================
    # 15. Observation variable
    # ==========================================================================
    Ab <- ab
    Ab ~ add(addSd_Ab)
  })
}
