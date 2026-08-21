Ivanova_2024_synucleinopathy_qsp <- function() {
  description <- paste(
    "QSP. Preclinical (mouse). Whole-brain quantitative systems pharmacology",
    "model of alpha-synuclein (aSyn) pathology in Parkinson's-disease-like mice,",
    "with a passive anti-aSyn immunotherapy layer (MEDI1341). 27 ODEs track three",
    "aSyn conformers (monomer, oligomer, fibril) across six brain matrices --",
    "non-striatal brain cells (BC) and striatal/substantia-nigra neurons (STR),",
    "non-striatal and striatal microglia (MG, MGstr), non-striatal and striatal",
    "brain interstitial fluid (BIF, BIFstr) -- plus monomer and oligomer in CSF.",
    "Intracellular aggregation is reversible oligomerization, irreversible",
    "oligomer maturation into fibrils, and fibril-templated (secondary-nucleation)",
    "growth from monomer. Aggregated aSyn drives microglial release of a composite",
    "'inflammatory mediator' variable (IM), which in turn accelerates aggregation",
    "and dopaminergic neurodegeneration; the surviving fraction of tyrosine-",
    "hydroxylase-positive neurons (perc_TH) feeds back on striatal aSyn synthesis.",
    "The antibody layer is a two-compartment plasma PK with saturable",
    "(Fc-receptor-like) peripheral and brain transport, equilibrium aSyn binding in",
    "BIF, and four selectable mechanisms of action toggled by the sw_mAb_* switches:",
    "block of neuronal uptake, enhanced microglial uptake, enhanced microglial",
    "degradation, and block of aggregation. Scenario switches (SW_TG, SW_VIRUS,",
    "SW_OVEREXPR_STR, SW_PFF, SW_IM_INH) reproduce the transgenic, viral-vector,",
    "preformed-fibril (PFF) and anti-inflammatory arms of Table 1. Defaults are the",
    "untreated wild-type mouse. aSyn states are CONCENTRATIONS in nM; antibody",
    "states are AMOUNTS in mg.",
    sep = " "
  )
  reference <- paste(
    "Ivanova O, Karelina T. Quantitative systems pharmacology model of",
    "alpha-synuclein pathology in Parkinson's disease-like mouse for",
    "investigation of passive immunotherapy mechanisms.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(10):1798-1809.",
    "doi:10.1002/psp4.13223.",
    "Model structure, rate laws and every parameter value are taken from the",
    "Supporting Information: the Heta-compiler model export",
    "(PSP4-13-1798-s001.xlsx; sheets Compartment / Species / Reaction / Record /",
    "Const / Process / TimeSwitcher) and the supplementary methods document",
    "(PSP4-13-1798-s002.docx, equations S1-S28).",
    sep = " "
  )
  vignette <- "Ivanova_2024_synucleinopathy_qsp"

  # Paper-mechanistic states: aSyn conformer x brain matrix, the composite
  # inflammatory-mediator variable, and the neurodegeneration log-survival
  # accumulator. None of these map onto a canonical nlmixr2lib compartment role.
  # The antibody states DO map onto canonicals and use them
  # (depot / central / peripheral1 / isf).
  paper_specific_compartments <- c(
    "aSyn_CSF", "aSynO_CSF",
    "aSyn_BIF", "aSynO_BIF", "aSynFm_BIF",
    "aSyn_BIFstr", "aSynO_BIFstr", "aSynFm_BIFstr",
    "aSyn_BC", "aSynO_BC", "aSynFm_BC",
    "aSyn_STR", "aSynO_STR", "aSynFm_STR",
    "aSyn_MG", "aSynO_MG", "aSynFm_MG",
    "aSyn_MGstr", "aSynO_MGstr", "aSynFm_MGstr",
    "IM_BIF", "IM_BIFstr", "l_prob"
  )

  units <- list(time = "h", dosing = "mg", concentration = "nM")

  # buildModelDb() auto-detects only compartments literally named `depot` /
  # `central`, which would advertise this model as antibody-dosed only. The
  # intrastriatal preformed-fibril (PFF) challenge of Table 1 is administered as
  # a dose into the striatal-ISF oligomer pool (TimeSwitcher!evt_PFFstr adds
  # Dose_PFFstr nM to aSynO_BIFstr), so it is declared explicitly.
  dosing <- c("depot", "central", "aSynO_BIFstr")

  # Issue #482. NOTE: the 20 aSyn states and the two IM states are Heta
  # CONCENTRATION states (nM / dimensionless), not amounts -- the flux
  # expressions carry the compartment volume explicitly and each d/dt is
  # divided by its volume. The four antibody states are amounts in mg.
  # analyte + specimen verified against the Species sheet `notes` column of
  # PSP4-13-1798-s001.xlsx and against supplementary equations S1-S28.
  compartmentData <- list(
    aSyn_CSF       = list(analyte = "alpha-synuclein (monomer)",   units = "nM", specimen = "CSF",        verified = TRUE),
    aSynO_CSF      = list(analyte = "alpha-synuclein (oligomer)",  units = "nM", specimen = "CSF",        verified = TRUE),
    aSyn_BIF       = list(analyte = "alpha-synuclein (monomer)",   units = "nM", specimen = "brain ISF",  verified = TRUE),
    aSynO_BIF      = list(analyte = "alpha-synuclein (oligomer)",  units = "nM", specimen = "brain ISF",  verified = TRUE),
    aSynFm_BIF     = list(analyte = "alpha-synuclein (fibril)",    units = "nM", specimen = "brain ISF",  verified = TRUE),
    aSyn_BIFstr    = list(analyte = "alpha-synuclein (monomer)",   units = "nM", specimen = "brain ISF",  verified = TRUE),
    aSynO_BIFstr   = list(analyte = "alpha-synuclein (oligomer)",  units = "nM", specimen = "brain ISF",  verified = TRUE),
    aSynFm_BIFstr  = list(analyte = "alpha-synuclein (fibril)",    units = "nM", specimen = "brain ISF",  verified = TRUE),
    aSyn_BC        = list(analyte = "alpha-synuclein (monomer)",   units = "nM", specimen = "tissue",     verified = TRUE),
    aSynO_BC       = list(analyte = "alpha-synuclein (oligomer)",  units = "nM", specimen = "tissue",     verified = TRUE),
    aSynFm_BC      = list(analyte = "alpha-synuclein (fibril)",    units = "nM", specimen = "tissue",     verified = TRUE),
    aSyn_STR       = list(analyte = "alpha-synuclein (monomer)",   units = "nM", specimen = "tissue",     verified = TRUE),
    aSynO_STR      = list(analyte = "alpha-synuclein (oligomer)",  units = "nM", specimen = "tissue",     verified = TRUE),
    aSynFm_STR     = list(analyte = "alpha-synuclein (fibril)",    units = "nM", specimen = "tissue",     verified = TRUE),
    aSyn_MG        = list(analyte = "alpha-synuclein (monomer)",   units = "nM", specimen = "tissue",     verified = TRUE),
    aSynO_MG       = list(analyte = "alpha-synuclein (oligomer)",  units = "nM", specimen = "tissue",     verified = TRUE),
    aSynFm_MG      = list(analyte = "alpha-synuclein (fibril)",    units = "nM", specimen = "tissue",     verified = TRUE),
    aSyn_MGstr     = list(analyte = "alpha-synuclein (monomer)",   units = "nM", specimen = "tissue",     verified = TRUE),
    aSynO_MGstr    = list(analyte = "alpha-synuclein (oligomer)",  units = "nM", specimen = "tissue",     verified = TRUE),
    aSynFm_MGstr   = list(analyte = "alpha-synuclein (fibril)",    units = "nM", specimen = "tissue",     verified = TRUE),
    IM_BIF         = list(analyte = "inflammatory mediators (composite)", units = "dimensionless", specimen = "brain ISF", verified = TRUE),
    IM_BIFstr      = list(analyte = "inflammatory mediators (composite)", units = "dimensionless", specimen = "brain ISF", verified = TRUE),
    l_prob         = list(analyte = "not applicable",              units = "dimensionless", specimen = "not applicable", verified = TRUE),
    depot          = list(analyte = "MEDI1341",                    units = "mg", specimen = "administration site", verified = TRUE),
    central        = list(analyte = "MEDI1341",                    units = "mg", specimen = "plasma",     verified = TRUE),
    peripheral1    = list(analyte = "MEDI1341",                    units = "mg", specimen = "plasma",     verified = TRUE),
    isf            = list(analyte = "MEDI1341",                    units = "mg", specimen = "brain ISF",  verified = TRUE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list()

  population <- list(
    species       = "mouse (C57BL/6 wild-type; Thy1/line-61 and A53T-M83 transgenic; rAAV2/7 and lentiviral A53T vector; TLR2 knockout)",
    n_subjects    = NA_integer_,
    n_studies     = 6L,
    disease_state = paste(
      "Parkinson's-disease-like alpha-synucleinopathy induced by transgenic or",
      "viral-vector aSyn overexpression, intrastriatal preformed-fibril (PFF)",
      "injection, or their combination; wild-type mice are the untreated",
      "reference state.",
      sep = " "
    ),
    dose_range = paste(
      "Anti-aSyn antibody 10-50 mg/kg IV weekly (default 10 mg/kg, i.e. 0.25 mg",
      "in a 25 g mouse; the paper's prediction runs use 20 mg/kg starting 1 week",
      "post-insult). Intrastriatal aSyn PFF challenge 5 ug (Dose_PFFstr = 93984 nM",
      "into striatal interstitial fluid); 12.5 ug in the A53T + PFF arm.",
      sep = " "
    ),
    regions = "Preclinical; model calibrated and validated against published mouse and rat datasets",
    notes = paste(
      "No individual-level cohort: this is a deterministic typical-value QSP model",
      "calibrated against aggregate published data. The six experimental",
      "PD-like mouse models are enumerated in Table 1 of the paper (Thy1/line-61",
      "and A53T-M83 transgenics; A53T + 12.5 ug PFF; rAAV2/7 and lentiviral A53T",
      "vectors; wild-type + 5 ug PFF). Body weight is fixed at 0.025 kg and every",
      "compartment volume is derived from it (Const sheet BW / BRC / BIFC / BCC /",
      "MGC / STRC / CSFC / PLC). Antibody PK was fitted to MEDI1341 rat data",
      "(Figure S2f) and the mouse-scaled parameters in the Const sheet were reused",
      "unchanged for prasinezumab and cinpanemab, whose rodent PK is unpublished.",
      sep = " "
    )
  )

  ini({
    # =================================================================
    # All values are from the Heta model export
    # PSP4-13-1798-s001.xlsx (Supporting Information of Ivanova 2024).
    # Sheet + row identifier is given after each value.
    #
    # These are mechanistic constants reported as point values with no
    # uncertainty, no IIV and no residual error, so they are kept on the
    # LINEAR scale under the paper's own symbol names (see
    # parameter-names.md "Endogenous / mechanistic parameters"). The four
    # antibody PK parameters that map onto canonical names are the
    # exception and are log-transformed (lcl / lvc / lvp / lq).
    #
    # Every value is wrapped in fixed(): the paper reports a single
    # calibrated QSP parameterisation and estimates nothing from the
    # user's data.
    # =================================================================

    # ---------------- Anatomy / compartment fractions ----------------
    BW   <- fixed(0.025);  label("Mouse body weight (kg)")                                  # Const!BW
    BRC  <- fixed(0.017);  label("Brain mass fraction of body weight (unitless)")            # Const!BRC
    BIFC <- fixed(0.19);   label("Interstitial-fluid fraction of brain (unitless)")          # Const!BIFC
    CSFC <- fixed(0.0017); label("CSF fraction of body weight (unitless)")                   # Const!CSFC
    BCC  <- fixed(0.7);    label("Brain-cell fraction of brain (unitless)")                  # Const!BCC
    MGC  <- fixed(0.03);   label("Microglia fraction of brain (unitless)")                   # Const!MGC
    STRC <- fixed(0.047);  label("Striatum + substantia nigra fraction of brain (unitless)") # Const!STRC
    PLC  <- fixed(0.04);   label("Plasma fraction of body weight (unitless)")                # Const!PLC

    # ---------------- aSyn aggregation ----------------
    k_olig_aSyn_BC   <- fixed(0.005);    label("aSyn oligomerization rate constant (L/nmol/h)")        # Const!k_olig_aSyn_BC
    k_mat_aSynFm_BC  <- fixed(0.000144); label("Oligomer maturation into fibril rate constant (1/h)")  # Const!k_mat_aSynFm_BC
    k_gr_aSynFm_BC   <- fixed(0.0025);   label("De novo fibril growth from monomer (L/nmol/h)")        # Const!k_gr_aSynFm_BC
    k_diss_aSynO     <- fixed(0.0191);   label("Oligomer dissociation rate constant (1/h)")            # Const!k_diss_aSynO
    Kd_F             <- fixed(1000);     label("Monomer binding to fibril ends, Kd (nM)")              # Const!Kd_F; AIC-selected over 40-7000 nM (Supp. methods)
    Eff              <- fixed(1000);     label("Fibril templating effect on oligomerization (unitless)")# Const!Eff
    N                <- fixed(714);      label("aSyn monomers per fibril (unitless)")                  # Const!N; 10000 kDa fibril / 14 kDa monomer (Polinski 2018)

    # ---------------- aSyn synthesis / clearance ----------------
    k_syn_aSyn_STR0    <- fixed(2.5);      label("Striatal-neuron aSyn synthesis rate (nM/h)")             # Record!k_syn_aSyn_STR start_
    BC_STR_coef        <- fixed(4);        label("Striatal / non-striatal aSyn synthesis ratio (unitless)")# Const!BC_STR_coef
    k_deg_aSyn_BC      <- fixed(2.801515); label("Monomer degradation in neurons (1/h)")                   # Const!k_deg_aSyn_BC
    k_upt_aSyn_BIF_BC  <- fixed(99.18175); label("Monomer uptake by neurons (1/h)")                        # Const!k_upt_aSyn_BIF_BC
    k_sec_aSyn_BC_BIF  <- fixed(1.5);      label("Monomer secretion by neurons (1/h)")                     # Const!k_sec_aSyn_BC_BIF
    k_upt_aSynO_BIF_BC <- fixed(2);        label("Oligomer uptake by neurons, Vmax (nmol/L/h)")            # Const!k_upt_aSynO_BIF_BC
    k_deg_aSynO_BC     <- fixed(0.25);     label("Oligomer degradation in neurons, Vmax (nmol/L/h)")       # Const!k_deg_aSynO_BC
    k_sec_aSynO_BC_BIF <- fixed(0.15);     label("Oligomer secretion by neurons (1/h)")                    # Const!k_sec_aSynO_BC_BIF
    k_upt_aSynO_BIF_MG <- fixed(0.0699);   label("Oligomer uptake by microglia (1/h)")                     # Const!k_upt_aSynO_BIF_MG
    k_deg_aSynO_MG     <- fixed(0.1);      label("Oligomer degradation in microglia (1/h)")                # Const!k_deg_aSynO_MG
    k_upt_aSynFm_BIF_MG <- fixed(0.0699);  label("Fibril uptake by microglia (1/h)")                       # Record!k_upt_aSynFm_BIF_MG
    EC50_aSynO_upt_BC  <- fixed(10);       label("Km for aggregate uptake by neurons (nM)")                # Const!EC50_aSynO_upt_BC
    EC50_aSynF_deg     <- fixed(30);       label("Km for aggregate degradation in neurons (nM)")           # Const!EC50_aSynF_deg
    # k_sec_aSynFm_STR_BIF is NOT a Const entry: Record!k_sec_aSynFm_STR_BIF is
    # assigned `k_sec_aSynFm_BC_BIF`, so it is derived in model() below rather
    # than given a numeric value here.

    # ---------------- Bulk flows ----------------
    Q_CSF_BIF <- fixed(1.1565e-05); label("Bulk flow CSF to brain ISF (L/h)")   # Const!Q_CSF_BIF
    Q_BIF_CSF <- fixed(1.1565e-05); label("Bulk flow brain ISF to CSF (L/h)")   # Const!Q_BIF_CSF
    Q_CSF_PL  <- fixed(1.25e-05);   label("Bulk flow CSF to plasma (L/h)")      # Const!Q_CSF_PL

    # ---------------- Microglial activation (inflammatory mediators) ----------------
    k_el_IM        <- fixed(0.01);      label("IM elimination scaling factor (unitless)")             # Const!k_el_IM
    Emax_aSynFm_IM <- fixed(785.2809);  label("Emax of fibrils on IM secretion (unitless)")           # Const!Emax_aSynFm_IM
    EC50_aSynFm_IM <- fixed(7420);      label("EC50 of fibrils on IM secretion (nM)")                 # Const!EC50_aSynFm_IM
    Nh_aSynFm_IM   <- fixed(3);         label("Hill coefficient, fibrils on IM secretion (unitless)") # Const!Nh_aSynFm_IM
    Emax_aSynO_IM  <- fixed(609.6);     label("Emax of oligomers on IM secretion (unitless)")         # Const!Emax_aSynO_IM
    EC50_aSynO_IM  <- fixed(5);         label("EC50 of oligomers on IM secretion (nM)")               # Const!EC50_aSynO_IM
    Nh_aSynO_IM    <- fixed(2);         label("Hill coefficient, oligomers on IM secretion (unitless)")# Const!Nh_aSynO_IM
    Emax_IM_aggr   <- fixed(8.39);      label("Emax of IM on aSyn aggregation (unitless)")            # Const!Emax_IM_aggr
    EC50_IM_aggr   <- fixed(74.6);      label("EC50 of IM on aSyn aggregation (nM)")                  # Const!EC50_IM_aggr
    Nh_IM_aggr     <- fixed(2);         label("Hill coefficient, IM on aggregation (unitless)")       # Const!Nh_IM_aggr

    # ---------------- Dopaminergic neurodegeneration ----------------
    k0_neurodeg         <- fixed(1e-05); label("Baseline dopaminergic neurodegeneration hazard (1/h)")   # Const!k0_neurodeg
    Emax_aSynO_neurodeg <- fixed(25.9);  label("Emax of oligomers on neurodegeneration (unitless)")      # Const!Emax_aSynO_neurodeg
    EC50_aSynO_neurodeg <- fixed(10);    label("EC50 of aggregates on neurodegeneration (nM)")           # Const!EC50_aSynO_neurodeg
    k_neurodeg_IM       <- fixed(67.2);  label("Microglial neurotoxicity coefficient (unitless)")        # Const!k_neurodeg_IM
    EC50_neurodeg_IM    <- fixed(203);   label("EC50 of IM on neurodegeneration (unitless)")             # Const!EC50_neurodeg_IM

    # ---------------- Reference / baseline levels used by binding and reporting ----------------
    asyn_M        <- fixed(0.0011);     label("Baseline monomeric aSyn in CSF (nM)")            # Const!asyn_M
    asyn_O        <- fixed(0.00232);    label("Baseline oligomeric aSyn in brain ISF (nM)")     # Const!asyn_O
    aSynM_BC_BL   <- fixed(0.458);      label("Baseline monomeric aSyn in neurons (nM)")        # Const!aSynM_BC_BL
    aSynO_BC_BL   <- fixed(0.04702547); label("Baseline oligomeric aSyn in neurons (nM)")       # Const!aSynO_BC_BL
    IM_BIF_BL     <- fixed(1.78);       label("Baseline whole-brain IM level (unitless)")       # Const!IM_BIF_BL
    IM_BIFstr_BL  <- fixed(0.9435);     label("Baseline striatal IM level (unitless)")          # Const!IM_BIFstr_BL
    sol_aSyn_BL   <- fixed(98);         label("Baseline soluble brain aSyn normaliser")         # Const!sol_aSyn_BL; see vignette Errata (scale not reproducible)
    insol_aSyn_BL <- fixed(1.895e-06);  label("Baseline insoluble brain aSyn normaliser (nmol)")# Const!insol_aSyn_BL

    # ---------------- Antibody PK (MEDI1341; fitted to rat data, mouse-scaled) ----------------
    lcl <- fixed(log(1.5267575e-05)); label("Antibody clearance from plasma (L/h)")            # Const!CL_Ig
    lvc <- fixed(log(0.0006875));     label("Antibody central volume (L)")                     # Const!Vc_Ig
    lvp <- fixed(log(0.0004244));     label("Antibody peripheral volume (L)")                  # Const!Vp_Ig
    lq  <- fixed(log(0.29231425));    label("Antibody plasma-periphery flow (L/h)")            # Const!Q_Ig
    # No canonical name exists for a plasma-to-brain flow or a brain clearance,
    # so these keep the paper's symbols on the linear scale.
    Q_PL_BR_Ig <- fixed(0.00075);      label("Antibody flow plasma to brain (L/h)")             # Const!Q_PL_BR_Ig
    CL_BR_Ig   <- fixed(1.544873e-05); label("Antibody clearance from brain (L/h)")             # Const!CL_BR_Ig
    Kd_Fc      <- fixed(10000);        label("Fc-receptor-like saturation constant (mg/L)")     # Const!Kd_Fc
    Mr_BP      <- fixed(150000);       label("Antibody molar mass (g/mol)")                     # Const!Mr_BP
    nM_mg      <- fixed(1e6);          label("mg-to-nmol unit conversion factor (unitless)")    # Const!nM_mg
    Dose_Ig0   <- fixed(10);           label("Reference antibody dose (mg/kg)")                 # Const!Dose_Ig0

    # ---------------- Antibody target binding and mechanism-of-action switches ----------------
    Kd_mon  <- fixed(8); label("Antibody affinity for aSyn monomer (nM)")     # Const!Kd_mon; recalibrated in vivo (Figure S3)
    Kd_MEDI <- fixed(8); label("Antibody affinity for aggregated aSyn (nM)")  # Const!Kd_MEDI; recalibrated in vivo (Figure S3)
    Emax_mAb_upt <- fixed(800); label("Emax of antibody on microglial aggregate uptake (unitless)")      # Const!Emax_mAb_upt
    Emax_mAb_deg <- fixed(200); label("Emax of antibody on microglial aggregate degradation (unitless)") # Const!Emax_mAb_deg
    # MoA switches. 1 = mechanism active, 0 = inactive. Const-sheet defaults are
    # reproduced here. Paper MoA definitions (Results, "Passive immunotherapy
    # simulation"): MoA1 = block neuronal uptake only (sw_mAb_BC_upt = 1, rest 0);
    # MoA2 = MoA1 + microglial uptake and degradation; MoA3 = block aggregation
    # only (sw_mAb_aggr = 1, rest 0); MoA4 = MoA3 + microglial uptake and
    # degradation. With no antibody dosed every switch is inert because the free
    # fractions all equal 1.
    sw_mAb_BC_upt   <- fixed(1); label("Switch: antibody blocks aSyn uptake by neurons")          # Const!sw_mAb_BC_upt
    sw_mAb_MG_upt   <- fixed(1); label("Switch: antibody enhances aSyn uptake by microglia")      # Const!sw_mAb_MG_upt
    sw_mAb_MG_deg   <- fixed(1); label("Switch: antibody enhances aSyn degradation by microglia") # Const!sw_mAb_MG_deg
    sw_mAb_aggr     <- fixed(1); label("Switch: antibody blocks intracellular aSyn aggregation")  # Const!sw_mAb_aggr
    sw_mAb_IM       <- fixed(0); label("Switch: antibody dampens microglial activation")          # Const!sw_mAb_IM
    sw_mAb_neurodeg <- fixed(0); label("Switch: antibody acts directly on neurodegeneration")     # Const!sw_mAb_neurodeg

    # ---------------- PD-like-model scenario switches ----------------
    # The Heta export encodes each experimental arm as a TimeSwitcher whose
    # `active` flag is FALSE by default (wild-type mouse). Each SW_* below turns
    # the corresponding TimeSwitcher on; the switch itself is an encoding device
    # introduced here and is not a paper parameter.
    SW_TG            <- fixed(0); label("Switch: transgenic aSyn overexpression from t = 0")    # TimeSwitcher!evt_Tg
    SW_VIRUS         <- fixed(0); label("Switch: viral aSyn overexpression at virus_time")      # TimeSwitcher!evt_virus
    SW_OVEREXPR_STR  <- fixed(0); label("Switch: striatal viral overexpression")                # TimeSwitcher!evt_overexpr_str
    SW_PFF           <- fixed(0); label("Switch: intrastriatal PFF challenge")                  # TimeSwitcher!evt_PFFstr
    SW_IM_INH        <- fixed(0); label("Switch: anti-inflammatory treatment")                  # TimeSwitcher!evt_IM_inh

    Tg_over            <- fixed(1.7);   label("Transgenic aSyn overexpression fold (unitless)")            # Const!Tg_over
    AAV                <- fixed(3);     label("Viral aSyn overexpression fold (unitless)")                 # Const!AAV
    virus_time         <- fixed(2000);  label("Time of viral vector injection (h)")                        # Const!virus_time
    rate_overexpr_str  <- fixed(10);    label("Striatal viral overexpression fold (unitless)")             # Const!rate_overexpr_str
    time_overexpr_str  <- fixed(0);     label("Time of striatal viral overexpression onset (h)")           # Const!time_overexpr_str
    time_PFFstr        <- fixed(2000);  label("Time of intrastriatal PFF injection (h)")                   # Const!time_PFFstr
    time_PFF_flow_off  <- fixed(2048);  label("End of the accelerated post-PFF clearance window (h)")      # TimeSwitcher!evt_PFF_flow start = 2048 h = time_PFFstr + 48 h (Table 1)
    coef_flow          <- fixed(4);     label("Fold increase in aggregate diffusion / clearance after PFF (unitless)") # Const!coef_flow
    time_IM_inh        <- fixed(2168);  label("Time of anti-inflammatory treatment start (h)")             # Const!time_IM_inh
    IM_inh_rate        <- fixed(0.5);   label("Residual fraction of IM synthesis under anti-inflammatory treatment (unitless)") # Const!IM_inh_rate; 0.5 = 50% inhibition, 0.1 = 90%
    Dose_PFFstr        <- fixed(93984); label("Intrastriatal PFF challenge dose (nM into aSynO_BIFstr)")   # Const!Dose_PFFstr; 5 ug / (BIFstr * 14 ug/nmol)
    PFF_ug             <- fixed(5);     label("Mass of aSyn PFF injected into the striatum (ug)")          # Const!PFF_ug; 12.5 ug in the A53T + PFF arm
    aSyn_MW            <- fixed(14);    label("aSyn monomer molar mass (ug/nmol, i.e. kDa)")               # Supp. methods, aSyn aggregation: 14 kDa monomer
  })

  model({
    # =================================================================
    # 1. Compartment volumes (L). Heta Compartment sheet, evaluated in
    #    dependency order: each line consumes the volumes above it.
    # =================================================================
    BR     <- BRC * BW                                  # Record!BR
    PL     <- PLC * BW                                  # Compartment!PL
    CSF    <- CSFC * BW                                 # Compartment!CSF
    BIFstr <- BIFC * STRC * BR                          # Compartment!BIFstr
    STR    <- BCC * (STRC * BR - BIFstr)                # Compartment!STR
    BIF    <- BIFC * (BR - STR - BIFstr)                # Compartment!BIF
    BC     <- BCC * (BR - BIF - STR - BIFstr)           # Compartment!BC
    MG     <- MGC * (BR - BIF - STRC * BR)              # Compartment!MG
    MGstr  <- MGC * STRC * BR                           # Compartment!MGstr

    # =================================================================
    # 2. Derived rate constants (Heta Record sheet). Several constants are
    #    tied to others by explicit modelling assumptions (Supporting
    #    Information, "Model assumptions").
    # =================================================================
    Q_BIFstr_BIF        <- Q_BIF_CSF / 100              # Record!Q_BIFstr_BIF; striatal/non-striatal ISF exchange 100x below BIF-CSF flow
    k_syn_aSyn_BC       <- k_syn_aSyn_STR0 / BC_STR_coef # Record!k_syn_aSyn_BC
    k_deg_aSyn_MG       <- k_deg_aSyn_BC                # Record!k_deg_aSyn_MG; monomer handling assumed non-specific
    k_upt_aSyn_BIF_MG   <- k_upt_aSyn_BIF_BC            # Record!k_upt_aSyn_BIF_MG
    k_sec_aSyn_MG_BIF   <- k_sec_aSyn_BC_BIF            # Record!k_sec_aSyn_MG_BIF
    k_upt_aSynFm_BIF_BC <- k_upt_aSynO_BIF_BC           # Record!k_upt_aSynFm_BIF_BC
    k_deg_aSynFm_BC     <- k_deg_aSynO_BC / 10          # Record!k_deg_aSynFm_BC
    k_deg_aSynFm_MG     <- k_deg_aSynO_MG / 10          # Record!k_deg_aSynFm_MG
    k_sec_aSynO_MG_BIF  <- k_deg_aSynO_MG / 3           # Record!k_sec_aSynO_MG_BIF; secretion/degradation ratio 1:3
    k_sec_aSynFm_MG_BIF <- k_deg_aSynFm_MG / 3          # Record!k_sec_aSynFm_MG_BIF
    k_sec_aSynFm_BC_BIF <- k_deg_aSynFm_BC / 3          # Record!k_sec_aSynFm_BC_BIF
    k_sec_aSynFm_STR_BIF <- k_sec_aSynFm_BC_BIF         # Record!k_sec_aSynFm_STR_BIF; assigned equal to the non-striatal constant
    k_diss_aSynFm       <- k_diss_aSynO / 100           # Record!k_diss_aSynFm; see vignette Errata (prose says /10 and "=")
    k_neurodeg_aSynFm   <- Emax_aSynO_neurodeg          # Record!k_neurodeg_aSynFm
    k0_IM               <- Q_BIF_CSF * k_el_IM / BIF    # Record!k0_IM; sets the healthy-brain IM level to 1
    aSynF_BC            <- aSynFm_BC / N                # Record!aSynF_BC; fibril-end concentration
    aSynF_STR           <- aSynFm_STR / N               # Record!aSynF_STR

    # =================================================================
    # 3. Scenario switches. The Heta TimeSwitchers apply a one-off
    #    multiplicative assignment at a fixed time; the equivalent
    #    algebraic step functions are written out here.
    # =================================================================
    tg_mult <- 1
    if (SW_TG > 0.5) {
      tg_mult <- Tg_over                                # TimeSwitcher!evt_Tg at t = 0
    }
    virus_mult <- 1
    if (SW_VIRUS > 0.5 && t >= virus_time) {
      virus_mult <- AAV                                 # TimeSwitcher!evt_virus
    }
    k_syn_aSyn_STR <- k_syn_aSyn_STR0 * tg_mult * virus_mult

    aSyn_over_STR <- 1
    if (SW_OVEREXPR_STR > 0.5 && t >= time_overexpr_str) {
      aSyn_over_STR <- rate_overexpr_str                # TimeSwitcher!evt_overexpr_str
    }

    # Accelerated aggregate diffusion / clearance in the 48 h after PFF
    # injection (Table 1: "increased elimination of oligomeric aSyn from the
    # brain 48 h after injection").
    coef_PFFstr <- 0
    if (SW_PFF > 0.5 && t >= time_PFFstr && t < time_PFF_flow_off) {
      coef_PFFstr <- coef_flow                          # Record!coef_PFFstr
    }

    sw_IM <- 1
    if (SW_IM_INH > 0.5 && t >= time_IM_inh) {
      sw_IM <- IM_inh_rate                              # TimeSwitcher!evt_IM_inh
    }

    # =================================================================
    # 4. Antibody exposure, target binding and free fractions.
    #    AB_complex_* solve the 1:1 equilibrium-binding quadratic; each
    #    Free_* is the unbound fraction of the corresponding aSyn pool.
    # =================================================================
    cl <- exp(lcl)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q  <- exp(lq)

    Ig_nM        <- nM_mg * central / vc / Mr_BP        # Record!Ig_nM; plasma antibody (nM)
    Ig_nM_ISFtot <- nM_mg * isf / BIF / Mr_BP           # Record!Ig_nM_ISFtot; total brain-ISF antibody (nM)
    mAb_tot      <- Ig_nM_ISFtot                        # Record!mAb_tot; site-of-action concentration assumed equal to BIF

    AB_complex_Asyn <- -0.5 * (-Kd_MEDI - asyn_O - mAb_tot +
      sqrt(Kd_MEDI^2 + 2 * Kd_MEDI * (mAb_tot + asyn_O) + (mAb_tot - asyn_O)^2))       # Record!AB_complex_Asyn
    Free_asyn <- 1 - AB_complex_Asyn / asyn_O                                          # Record!Free_asyn

    AB_complex_aSynM <- -0.5 * (-Kd_mon - asyn_M - mAb_tot +
      sqrt(Kd_mon^2 + 2 * Kd_mon * (mAb_tot + asyn_M) + (mAb_tot - asyn_M)^2))         # Record!AB_complex_aSynM
    Free_aSynM <- 1 - AB_complex_aSynM / asyn_M                                        # Record!Free_aSynM

    AB_complex_AsynO_BC <- -0.5 * (-Kd_MEDI - aSynO_BC_BL - mAb_tot +
      sqrt(Kd_MEDI^2 + 2 * Kd_MEDI * (mAb_tot + aSynO_BC_BL) + (mAb_tot - aSynO_BC_BL)^2)) # Record!AB_complex_AsynO_BC
    Free_aSynO_BC <- 1 - sw_mAb_aggr * AB_complex_AsynO_BC / aSynO_BC_BL               # Record!Free_aSynO_BC

    AB_complex_AsynM_BC <- -0.5 * (-Kd_mon - aSynM_BC_BL - mAb_tot +
      sqrt(Kd_mon^2 + 2 * Kd_mon * (mAb_tot + aSynM_BC_BL) + (mAb_tot - aSynM_BC_BL)^2)) # Record!AB_complex_AsynM_BC
    Free_aSynM_BC <- 1 - sw_mAb_aggr * AB_complex_AsynM_BC / aSynM_BC_BL              # Record!Free_aSynM_BC

    # =================================================================
    # 5. Inflammation and neurodegeneration feedbacks.
    # =================================================================
    aggr_coef_IM <- 1 + Emax_IM_aggr * (IM_BIF / (EC50_IM_aggr + IM_BIF))^Nh_IM_aggr   # Record!aggr_coef_IM
    aggr_coef_IMstr <- 1 + Emax_IM_aggr * (IM_BIFstr / EC50_IM_aggr)^Nh_IM_aggr /
      (1 + (IM_BIFstr / EC50_IM_aggr)^Nh_IM_aggr)                                      # Record!aggr_coef_IMstr

    TH_neurodeg <- 100 * (1 - exp(l_prob))   # Record!TH_neurodeg; % dopaminergic neurons lost
    perc_TH     <- 100 - TH_neurodeg         # Record!perc_TH; % viable TH+ neurons

    # =================================================================
    # 6. Reaction rate laws, verbatim from the Heta Reaction sheet
    #    (fluxes in nmol/h; the compartment volume is carried explicitly).
    # =================================================================
    V_tr_aSyn_CSF_BIF <- Q_CSF_BIF * (aSyn_CSF - aSyn_BIF)
    V_tr_aSynO_CSF_BIF <- Q_CSF_BIF * (aSynO_CSF - aSynO_BIF) * (1 + coef_PFFstr)
    V_tr_aSyn_BIFstr_BIF <- Q_BIFstr_BIF * (aSyn_BIFstr - aSyn_BIF)
    V_tr_aSynO_BIFstr_BIF <- Q_BIFstr_BIF * (aSynO_BIFstr - aSynO_BIF) * (1 + coef_PFFstr)
    V_tr_aSynFm_BIFstr_BIF <- Q_BIFstr_BIF * (aSynFm_BIFstr - aSynFm_BIF)
    V_tr_aSyn_BIFstr_CSF <- Q_BIF_CSF * STRC * (aSyn_BIFstr - aSyn_CSF)
    V_tr_aSynO_BIFstr_CSF <- Q_BIF_CSF * STRC * (aSynO_BIFstr - aSynO_CSF) * (1 + coef_PFFstr)
    V_el_aSynO_CSF <- Q_CSF_PL * aSynO_CSF * (0 + coef_PFFstr)

    V_syn_aSyn_BC <- BC * k_syn_aSyn_BC * aggr_coef_IM
    V_upt_aSyn_BIF_BC <- BC * k_upt_aSyn_BIF_BC * aSyn_BIF * Free_aSynM^sw_mAb_BC_upt
    V_upt_aSynO_BIF_BC <- BC * (k_upt_aSynO_BIF_BC * aSynO_BIF * Free_asyn^sw_mAb_BC_upt) /
      (aSynO_BIF * Free_asyn^sw_mAb_BC_upt + EC50_aSynO_upt_BC)
    V_upt_aSynFm_BIF_BC <- BC * (k_upt_aSynFm_BIF_BC * aSynFm_BIF * Free_asyn^sw_mAb_BC_upt /
      (aSynFm_BIF * Free_asyn^sw_mAb_BC_upt + EC50_aSynO_upt_BC))
    V_sec_aSyn_BC_BIF <- BC * k_sec_aSyn_BC_BIF * aSyn_BC
    V_sec_aSynO_BC_BIF <- BC * k_sec_aSynO_BC_BIF * aSynO_BC
    V_sec_aSynFm_BC_BIF <- BC * k_sec_aSynFm_BC_BIF * aSynFm_BC
    V_deg_aSyn_BC <- BC * k_deg_aSyn_BC * aSyn_BC
    V_deg_aSynO_BC <- BC * k_deg_aSynO_BC * aSynO_BC / (EC50_aSynF_deg + aSynO_BC)
    V_deg_aSynFm_BC <- BC * k_deg_aSynFm_BC * aSynFm_BC / (EC50_aSynF_deg + aSynFm_BC)
    V_olig_aSyn_BC <- BC * k_olig_aSyn_BC * Free_aSynM_BC^2 * aSyn_BC^2 * aggr_coef_IM *
      (1 + Free_aSynO_BC * Eff * aSynF_BC / (aSyn_BC + Kd_F))
    V_mat_aSynFm_BC <- BC * k_mat_aSynFm_BC * aSynO_BC
    V_diss_aSynO_BC <- BC * k_diss_aSynO * aSynO_BC
    V_gr_aSynFm_BC <- BC * k_gr_aSynFm_BC * aSyn_BC * Free_aSynM_BC * Free_aSynO_BC *
      aSynF_BC * aggr_coef_IM
    V_diss_aSynFm_BC <- BC * k_diss_aSynFm * aSynFm_BC

    V_syn_aSyn_STR <- aSyn_over_STR * STR * k_syn_aSyn_STR * perc_TH / 100
    V_sec_aSyn_STR_BIFstr <- STR * k_sec_aSyn_BC_BIF * aSyn_STR
    V_sec_aSynO_STR_BIFstr <- STR * k_sec_aSynO_BC_BIF * aSynO_STR
    V_sec_aSynFm_STR_BIFstr <- STR * k_sec_aSynFm_BC_BIF * aSynFm_STR
    V_sec_aSynFm_STR_BIF <- STR * k_sec_aSynFm_STR_BIF * aSynFm_STR
    V_upt_aSyn_BIFstr_STR <- STR * k_upt_aSyn_BIF_BC * aSyn_BIFstr * Free_aSynM^sw_mAb_BC_upt
    V_upt_aSynO_BIFstr_STR <- STR * (k_upt_aSynO_BIF_BC * aSynO_BIFstr * Free_asyn^sw_mAb_BC_upt) /
      (aSynO_BIFstr * Free_asyn^sw_mAb_BC_upt + EC50_aSynO_upt_BC)
    V_upt_aSynFm_BIFstr_STR <- STR * (k_upt_aSynFm_BIF_BC * aSynFm_BIFstr * Free_asyn^sw_mAb_BC_upt /
      (aSynFm_BIFstr * Free_asyn^sw_mAb_BC_upt + EC50_aSynO_upt_BC))
    V_deg_aSyn_STR <- STR * k_deg_aSyn_BC * aSyn_STR * perc_TH / 100
    V_deg_aSynO_STR <- STR * k_deg_aSynO_BC * aSynO_STR / (EC50_aSynF_deg + aSynO_STR)
    V_deg_aSynFm_STR <- STR * k_deg_aSynFm_BC * aSynFm_STR / (EC50_aSynF_deg + aSynFm_STR)
    V_olig_aSyn_STR <- STR * k_olig_aSyn_BC * Free_aSynM_BC^2 * aSyn_STR^2 * aggr_coef_IMstr *
      (1 + Free_aSynO_BC * Eff * aSynF_STR / (aSyn_STR + Kd_F))
    V_mat_aSynFm_STR <- STR * k_mat_aSynFm_BC * aSynO_STR
    V_gr_aSynFm_STR <- STR * k_gr_aSynFm_BC * aSyn_STR * Free_aSynM_BC * Free_aSynO_BC *
      aSynF_STR * aggr_coef_IMstr
    V_diss_aSynO_STR <- STR * k_diss_aSynO * aSynO_STR
    V_diss_aSynFm_STR <- STR * k_diss_aSynFm * aSynFm_STR

    V_upt_aSyn_BIF_MG <- MG * k_upt_aSyn_BIF_MG * (aSyn_BIF * Free_aSynM^sw_mAb_MG_upt +
      Emax_mAb_upt * sw_mAb_MG_upt * aSyn_BIF * AB_complex_aSynM / asyn_M)
    V_upt_aSynO_BIF_MG <- MG * k_upt_aSynO_BIF_MG * (aSynO_BIF * Free_asyn^sw_mAb_MG_upt +
      Emax_mAb_upt * sw_mAb_MG_upt * aSynO_BIF * AB_complex_Asyn / asyn_O)
    V_upt_aSynFm_BIF_MG <- MG * k_upt_aSynFm_BIF_MG * (aSynFm_BIF * Free_asyn^sw_mAb_MG_upt +
      Emax_mAb_upt * sw_mAb_MG_upt * aSynFm_BIF * AB_complex_Asyn / asyn_O)
    V_sec_aSyn_MG_BIF <- MG * k_sec_aSyn_MG_BIF * aSyn_MG
    V_sec_aSynO_MG_BIF <- MG * k_sec_aSynO_MG_BIF * aSynO_MG
    V_sec_aSynFm_MG_BIF <- MG * k_sec_aSynFm_MG_BIF * aSynFm_MG
    V_deg_aSyn_MG <- MG * k_deg_aSyn_MG * (aSyn_MG * Free_aSynM^sw_mAb_MG_deg +
      Emax_mAb_deg * sw_mAb_MG_deg * aSyn_MG * AB_complex_aSynM / asyn_M)
    V_deg_aSynO_MG <- MG * k_deg_aSynO_MG * (aSynO_MG * Free_asyn^sw_mAb_MG_deg +
      Emax_mAb_deg * sw_mAb_MG_deg * aSynO_MG * AB_complex_Asyn / asyn_O)
    V_deg_aSynFm_MG <- MG * k_deg_aSynFm_MG * (aSynFm_MG * Free_asyn^sw_mAb_MG_deg +
      Emax_mAb_deg * sw_mAb_MG_deg * aSynFm_MG * AB_complex_Asyn / asyn_O)

    V_upt_aSyn_BIFstr_MGstr <- MGstr * k_upt_aSyn_BIF_MG * (aSyn_BIFstr * Free_aSynM^sw_mAb_MG_upt +
      Emax_mAb_upt * sw_mAb_MG_upt * aSyn_BIFstr * AB_complex_aSynM / asyn_M)
    V_upt_aSynO_BIFstr_MGstr <- MGstr * k_upt_aSynO_BIF_MG * (aSynO_BIFstr * Free_asyn^sw_mAb_MG_upt +
      Emax_mAb_upt * sw_mAb_MG_upt * aSynO_BIFstr * AB_complex_Asyn / asyn_O)
    V_upt_aSynFm_BIFstr_MGstr <- MGstr * k_upt_aSynFm_BIF_MG * (aSynFm_BIFstr * Free_asyn^sw_mAb_MG_upt +
      Emax_mAb_upt * sw_mAb_MG_upt * aSynFm_BIFstr * AB_complex_Asyn / asyn_O)
    V_sec_aSyn_MGstr_BIFstr <- MGstr * k_sec_aSyn_MG_BIF * aSyn_MGstr
    V_sec_aSynO_MGstr_BIFstr <- MGstr * k_sec_aSynO_MG_BIF * aSynO_MGstr
    V_sec_aSynFm_MGstr_BIFstr <- MGstr * k_sec_aSynFm_MG_BIF * aSynFm_MGstr
    V_deg_aSyn_MGstr <- MGstr * k_deg_aSyn_MG * (aSyn_MGstr * Free_aSynM^sw_mAb_MG_deg +
      Emax_mAb_deg * sw_mAb_MG_deg * aSyn_MGstr * AB_complex_aSynM / asyn_M)
    V_deg_aSynO_MGstr <- MGstr * k_deg_aSynO_MG * (aSynO_MGstr * Free_asyn^sw_mAb_MG_deg +
      Emax_mAb_deg * sw_mAb_MG_deg * aSynO_MGstr * AB_complex_Asyn / asyn_O)
    V_deg_aSynFm_MGstr <- MGstr * k_deg_aSynFm_MG * (aSynFm_MGstr * Free_asyn^sw_mAb_MG_deg +
      Emax_mAb_deg * sw_mAb_MG_deg * aSynFm_MGstr * AB_complex_Asyn / asyn_O)

    # Process-sheet rate laws. Heta Processes contribute directly to d/dt and are
    # NOT divided by a compartment volume (confirmed two ways: this is what makes
    # the healthy-brain IM level equal 1 as stated in the Supporting Information,
    # and what keeps l_prob in a physiological range).
    V_tr_IM_BIF_CSF <- k_el_IM * Q_BIF_CSF * IM_BIF / BIF
    V_tr_IM_BIFstr_CSF <- k_el_IM * Q_BIF_CSF * STRC * IM_BIFstr / BIFstr
    V_sec_IM_MG_BIF <- sw_IM * k0_IM *
      (1 + Emax_aSynFm_IM * ((aSynFm_BC + aSynFm_BIF) * Free_asyn^sw_mAb_IM / EC50_aSynFm_IM)^Nh_aSynFm_IM +
         Emax_aSynO_IM * ((aSynO_BC + aSynO_BIF) * Free_asyn^sw_mAb_IM / EC50_aSynO_IM)^Nh_aSynO_IM) /
      (1 + ((aSynFm_BC + aSynFm_BIF) * Free_asyn^sw_mAb_IM / EC50_aSynFm_IM)^Nh_aSynFm_IM +
         ((aSynO_BC + aSynO_BIF) * Free_asyn^sw_mAb_IM / EC50_aSynO_IM)^Nh_aSynO_IM)
    V_sec_IM_MGstr_BIFstr <- sw_IM * k0_IM *
      (1 + Emax_aSynFm_IM * ((aSynFm_STR + aSynFm_BIFstr) * Free_asyn^sw_mAb_IM / EC50_aSynFm_IM)^Nh_aSynFm_IM +
         Emax_aSynO_IM * ((aSynO_STR + aSynO_BIFstr) * Free_asyn^sw_mAb_IM / EC50_aSynO_IM)^Nh_aSynO_IM) /
      (1 + ((aSynFm_STR + aSynFm_BIFstr) * Free_asyn^sw_mAb_IM / EC50_aSynFm_IM)^Nh_aSynFm_IM +
         ((aSynO_STR + aSynO_BIFstr) * Free_asyn^sw_mAb_IM / EC50_aSynO_IM)^Nh_aSynO_IM)
    V_haz_deg <- -k0_neurodeg *
      (1 + aSynFm_STR * Free_asyn^sw_mAb_neurodeg * k_neurodeg_aSynFm / EC50_aSynO_neurodeg +
         k_neurodeg_IM * (IM_BIFstr / EC50_neurodeg_IM) +
         Emax_aSynO_neurodeg * aSynO_STR * Free_asyn^sw_mAb_neurodeg / EC50_aSynO_neurodeg) /
      (1 + aSynFm_STR * Free_asyn^sw_mAb_neurodeg / EC50_aSynO_neurodeg +
         aSynO_STR * Free_asyn^sw_mAb_neurodeg / EC50_aSynO_neurodeg +
         (IM_BIFstr / EC50_neurodeg_IM))

    # Antibody disposition. V_PL_Ig / V_BR_Ig carry a saturable
    # Fc-receptor-like denominator exactly as published.
    V_abs_Ig   <- depot
    V_el_Ig    <- cl * central / vc
    V_PL_Ig    <- q * (central - vc / vp * peripheral1) /
      (central / vc + peripheral1 / vp + Kd_Fc) / vc
    V_BR_Ig    <- Q_PL_BR_Ig * central / (central / vc + Kd_Fc) / vc
    V_Ig_CL_BR <- CL_BR_Ig * isf / BIF

    # =================================================================
    # 7. ODE system. Supplementary equations S1-S28; stoichiometry taken
    #    from the `actors` column of the Heta Reaction / Process sheets.
    # =================================================================
    # CSF (S1, S3)
    d/dt(aSyn_CSF) <- (-V_tr_aSyn_CSF_BIF + V_tr_aSyn_BIFstr_CSF) / CSF
    d/dt(aSynO_CSF) <- (-V_tr_aSynO_CSF_BIF + V_tr_aSynO_BIFstr_CSF - V_el_aSynO_CSF) / CSF

    # Non-striatal brain interstitial fluid (S2, S4, S8)
    d/dt(aSyn_BIF) <- (V_tr_aSyn_CSF_BIF + V_tr_aSyn_BIFstr_BIF - V_upt_aSyn_BIF_BC +
      V_sec_aSyn_BC_BIF - V_upt_aSyn_BIF_MG + V_sec_aSyn_MG_BIF) / BIF
    d/dt(aSynO_BIF) <- (V_tr_aSynO_CSF_BIF + V_tr_aSynO_BIFstr_BIF - V_upt_aSynO_BIF_BC +
      V_sec_aSynO_BC_BIF - V_upt_aSynO_BIF_MG + V_sec_aSynO_MG_BIF) / BIF
    d/dt(aSynFm_BIF) <- (V_tr_aSynFm_BIFstr_BIF - V_upt_aSynFm_BIF_BC + V_sec_aSynFm_BC_BIF +
      V_sec_aSynFm_STR_BIF - V_upt_aSynFm_BIF_MG + V_sec_aSynFm_MG_BIF) / BIF

    # Striatal interstitial fluid (S5, S6, S7)
    d/dt(aSyn_BIFstr) <- (-V_tr_aSyn_BIFstr_BIF - V_tr_aSyn_BIFstr_CSF + V_sec_aSyn_STR_BIFstr -
      V_upt_aSyn_BIFstr_STR - V_upt_aSyn_BIFstr_MGstr + V_sec_aSyn_MGstr_BIFstr) / BIFstr
    d/dt(aSynO_BIFstr) <- (-V_tr_aSynO_BIFstr_BIF - V_tr_aSynO_BIFstr_CSF + V_sec_aSynO_STR_BIFstr -
      V_upt_aSynO_BIFstr_STR - V_upt_aSynO_BIFstr_MGstr + V_sec_aSynO_MGstr_BIFstr) / BIFstr
    d/dt(aSynFm_BIFstr) <- (-V_tr_aSynFm_BIFstr_BIF + V_sec_aSynFm_STR_BIFstr -
      V_upt_aSynFm_BIFstr_STR - V_upt_aSynFm_BIFstr_MGstr + V_sec_aSynFm_MGstr_BIFstr) / BIFstr

    # Non-striatal brain cells (S9, S10, S11)
    d/dt(aSyn_BC) <- (V_syn_aSyn_BC + V_upt_aSyn_BIF_BC - V_sec_aSyn_BC_BIF - V_deg_aSyn_BC -
      V_olig_aSyn_BC + V_diss_aSynO_BC - V_gr_aSynFm_BC + V_diss_aSynFm_BC) / BC
    d/dt(aSynO_BC) <- (V_upt_aSynO_BIF_BC - V_sec_aSynO_BC_BIF - V_deg_aSynO_BC +
      V_olig_aSyn_BC - V_mat_aSynFm_BC - V_diss_aSynO_BC) / BC
    d/dt(aSynFm_BC) <- (V_upt_aSynFm_BIF_BC - V_sec_aSynFm_BC_BIF - V_deg_aSynFm_BC +
      V_mat_aSynFm_BC + V_gr_aSynFm_BC - V_diss_aSynFm_BC) / BC

    # Striatal / substantia-nigra neurons (S12, S13, S14)
    d/dt(aSyn_STR) <- (V_syn_aSyn_STR - V_sec_aSyn_STR_BIFstr + V_upt_aSyn_BIFstr_STR -
      V_deg_aSyn_STR - V_olig_aSyn_STR - V_gr_aSynFm_STR + V_diss_aSynO_STR +
      V_diss_aSynFm_STR) / STR
    d/dt(aSynO_STR) <- (-V_sec_aSynO_STR_BIFstr + V_upt_aSynO_BIFstr_STR - V_deg_aSynO_STR +
      V_olig_aSyn_STR - V_mat_aSynFm_STR - V_diss_aSynO_STR) / STR
    d/dt(aSynFm_STR) <- (-V_sec_aSynFm_STR_BIFstr - V_sec_aSynFm_STR_BIF +
      V_upt_aSynFm_BIFstr_STR - V_deg_aSynFm_STR + V_mat_aSynFm_STR + V_gr_aSynFm_STR -
      V_diss_aSynFm_STR) / STR

    # Non-striatal microglia (S15, S16, S17)
    d/dt(aSyn_MG) <- (V_upt_aSyn_BIF_MG - V_sec_aSyn_MG_BIF - V_deg_aSyn_MG) / MG
    d/dt(aSynO_MG) <- (V_upt_aSynO_BIF_MG - V_sec_aSynO_MG_BIF - V_deg_aSynO_MG) / MG
    d/dt(aSynFm_MG) <- (V_upt_aSynFm_BIF_MG - V_sec_aSynFm_MG_BIF - V_deg_aSynFm_MG) / MG

    # Striatal microglia (S18, S19, S20)
    d/dt(aSyn_MGstr) <- (V_upt_aSyn_BIFstr_MGstr - V_sec_aSyn_MGstr_BIFstr - V_deg_aSyn_MGstr) / MGstr
    d/dt(aSynO_MGstr) <- (V_upt_aSynO_BIFstr_MGstr - V_sec_aSynO_MGstr_BIFstr - V_deg_aSynO_MGstr) / MGstr
    d/dt(aSynFm_MGstr) <- (V_upt_aSynFm_BIFstr_MGstr - V_sec_aSynFm_MGstr_BIFstr - V_deg_aSynFm_MGstr) / MGstr

    # Inflammatory mediators (S21, S22) and neurodegeneration hazard (S28)
    d/dt(IM_BIF) <- V_sec_IM_MG_BIF - V_tr_IM_BIF_CSF
    d/dt(IM_BIFstr) <- V_sec_IM_MGstr_BIFstr - V_tr_IM_BIFstr_CSF
    d/dt(l_prob) <- V_haz_deg

    # Antibody (S24, S25, S26, S27)
    d/dt(depot) <- -V_abs_Ig
    d/dt(central) <- V_abs_Ig - V_el_Ig - V_PL_Ig - V_BR_Ig
    d/dt(peripheral1) <- V_PL_Ig
    d/dt(isf) <- V_BR_Ig - V_Ig_CL_BR

    # =================================================================
    # 8. Initial conditions. Species sheet `assignments.start_`: every aSyn
    #    pool starts at a numerically negligible 1e-20 nM and the wild-type
    #    baseline is reached by integrating forward (the paper's insults
    #    start at t = 2000 h).
    # =================================================================
    aSyn_CSF(0) <- 1e-20
    aSynO_CSF(0) <- 1e-20
    aSyn_BIF(0) <- 1e-20
    aSynO_BIF(0) <- 1e-20
    aSynFm_BIF(0) <- 1e-20
    aSyn_BIFstr(0) <- 1e-20
    aSynO_BIFstr(0) <- 1e-20
    aSynFm_BIFstr(0) <- 1e-20
    aSyn_BC(0) <- 1e-20
    aSynO_BC(0) <- 1e-20
    aSynFm_BC(0) <- 1e-20
    aSyn_STR(0) <- 1e-20
    aSynO_STR(0) <- 1e-20
    aSynFm_STR(0) <- 1e-20
    aSyn_MG(0) <- 1e-20
    aSynO_MG(0) <- 1e-20
    aSynFm_MG(0) <- 1e-20
    aSyn_MGstr(0) <- 1e-20
    aSynO_MGstr(0) <- 1e-20
    aSynFm_MGstr(0) <- 1e-20

    # =================================================================
    # 9. Reported outputs (Heta Record sheet).
    # =================================================================
    Cc <- Ig_nM                                        # plasma antibody (nM)
    Ig_nM_ISF <- Ig_nM_ISFtot - AB_complex_Asyn        # Record!Ig_nM_ISF; free brain-ISF antibody (nM)
    Ig_nM_BR  <- nM_mg * isf / BR / Mr_BP              # Record!Ig_nM_BR
    Dose_Ig   <- Dose_Ig0 * BW                         # Record!Dose_Ig; mg per administration

    aSyn_BR <- (aSyn_STR * STR + aSyn_BIF * BIF + aSyn_BC * BC + aSyn_MG * MG +
      aSyn_MGstr * MGstr + aSyn_BIFstr * BIFstr) / BR                            # Record!aSyn_BR
    aSynO_BR <- (aSynO_STR * STR + aSynO_BIFstr * BIFstr + aSynO_MGstr * MGstr +
      aSynO_BIF * BIF + aSynO_BC * BC + aSynO_MG * MG) / BR                      # Record!aSynO_BR
    aSynFm_BR <- (aSynFm_STR * STR + aSynFm_BC * BC + aSynFm_MG * MG +
      aSynFm_BIF * BIF) / BR                                                     # Record!aSynFm_BR
    aSyn_tot_BR <- ((aSyn_STR + aSynO_STR + aSynFm_STR) * STR +
      (aSyn_BIF + aSynO_BIF + aSynFm_BIF) * BIF +
      (aSyn_BIFstr + aSynO_BIFstr + aSynFm_BIFstr) * BIFstr +
      (aSyn_BC + aSynO_BC + aSynFm_BC) * BC +
      (aSyn_MG + aSynO_MG + aSynFm_MG) * MG +
      (aSyn_MGstr + aSynO_MGstr + aSynFm_MGstr) * MGstr) / BR                    # Record!aSyn_tot_BR
    sol_aSyn_BR <- ((aSyn_STR + aSynO_STR) * STR + (aSyn_BIF + aSynO_BIF) * BIF +
      (aSyn_BIFstr + aSynO_BIFstr) * BIFstr + (aSyn_BC + aSynO_BC) * BC +
      (aSyn_MG + aSynO_MG) * MG + (aSyn_MGstr + aSynO_MGstr) * MGstr) / BR       # Record!sol_aSyn_BR
    aggr_aSyn_BR <- ((aSynO_STR + aSynFm_STR) * STR + (aSynO_BIF + aSynFm_BIF) * BIF +
      (aSynO_BIFstr + aSynFm_BIFstr) * BIFstr + (aSynO_BC + aSynFm_BC) * BC +
      (aSynO_MG + aSynFm_MG) * MG + (aSynO_MGstr + aSynFm_MGstr) * MGstr) / BR   # Record!aggr_aSyn_BR
    aggr_aSyn_STR <- aSynO_STR + aSynFm_STR                                      # Record!aggr_aSyn_STR

    aSyn_tot_BIF <- aSyn_BIF + aSynO_BIF + aSynFm_BIF                            # Record!aSyn_tot_BIF
    aSyn_tot_BIFstr <- aSyn_BIFstr + aSynO_BIFstr + aSynFm_BIFstr                # Record!aSyn_tot_BIFstr
    aSyn_tot_STR <- ((aSyn_STR + aSynO_STR + aSynFm_STR) * STR +
      (aSyn_MG + aSynO_MG + aSynFm_MG) * MG +
      (aSyn_BIFstr + aSynO_BIFstr + aSynFm_BIFstr) * BIFstr) / (BR * STRC)       # Record!aSyn_tot_STR
    sol_aSyn_STR <- aSyn_STR + aSynO_STR + aSyn_BIFstr + aSynO_BIFstr +
      aSyn_MGstr + aSynO_MGstr                                                   # Record!sol_aSyn_STR
    aSyn_tot_CSF <- aSyn_CSF + aSynO_CSF                                         # Record!aSyn_tot_CSF
    aSynO_CSF_perc <- aSynO_CSF * 100 / (aSyn_CSF + aSynO_CSF)                   # Record!aSynO_CSF_perc
    aSynO_BR_perc <- aSynO_BR * 100 / aSyn_tot_BR                                # Record!aSynO_BR_perc
    aSynFm_BR_perc <- aSynFm_BR * 100 / aSyn_tot_BR                              # Record!aSynFm_BR_perc

    # Axonal aSyn: no axonal compartment, so the average of the striatal
    # intraneuronal and striatal interstitial pools is used (Results, Figure 3c).
    aSyn_axons <- (aSyn_tot_STR * STR + aSyn_tot_BIFstr * BIFstr) / (STR + BIFstr) # Record!aSyn_axons

    insol_aSyn_BR <- aSynFm_STR * STR + aSynFm_BIF * BIF + aSynFm_BIFstr * BIFstr +
      aSynFm_BC * BC + aSynFm_MG * MG + aSynFm_MGstr * MGstr                     # Record!insol_aSyn_BR
    perc_insol_aSyn_BR <- 100 * insol_aSyn_BR / insol_aSyn_BL                    # Record!perc_insol_aSyn_BR
    perc_sol_aSyn_BR <- 100 * ((aSyn_STR + aSynO_STR) * STR +
      (aSyn_BIF + aSynO_BIF) * BIF + (aSyn_BIFstr + aSynO_BIFstr) * BIFstr +
      (aSyn_BC + aSynO_BC) * BC + (aSyn_MG + aSynO_MG) * MG +
      (aSyn_MGstr + aSynO_MGstr) * MGstr) / sol_aSyn_BL                          # Record!perc_sol_aSyn_BR

    IM_BIF_tot <- IM_BIF + IM_BIFstr                                             # Record!IM_BIF_tot
    perc_IM_tot <- 100 * (IM_BIF + IM_BIFstr) / IM_BIF_BL                        # Record!perc_IM_tot
    perc_IM_BIFstr <- 100 * IM_BIFstr / IM_BIFstr_BL                             # Record!perc_IM_BIFstr

    Month <- t / 30 / 24                                                         # Record!Month
    DPI <- (t - time_PFFstr) / 24                                                # Record!DPI; days post-injection

    # Dose amounts, surfaced so an event table can be built from the model
    # itself. Dose_PFFstr_nM is the value evt_PFFstr adds to aSynO_BIFstr;
    # aSyn_dose_BIFstr_ug is the paper's own conversion of the injected mass
    # into a striatal-ISF concentration. The two agree to 0.13%, which
    # independently confirms the derived compartment volumes above.
    Dose_PFFstr_nM <- Dose_PFFstr                                                # Const!Dose_PFFstr
    aSyn_dose_BIFstr_ug <- PFF_ug / (BIFstr * aSyn_MW)                           # Record!aSyn_dose_BIFstr_ug
  })
}
