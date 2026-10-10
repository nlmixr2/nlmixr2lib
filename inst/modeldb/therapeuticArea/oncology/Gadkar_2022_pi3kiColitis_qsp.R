Gadkar_2022_pi3kiColitis_qsp <- function() {
  description <- paste(
    "QSP. Mechanistic model of colonic mucosal inflammation and epithelial",
    "barrier damage under prolonged PI3K-inhibitor treatment, used to",
    "predict the onset of diarrhea and colitis. A two-compartment oral",
    "(or IV bolus) PK model drives unbound mucosal and blood",
    "concentrations; the drug inhibits each PI3K isoform (alpha, beta,",
    "gamma, delta) through its own IC50, and the isoform activities",
    "modulate epithelial-cell formation and apoptosis, APC, effector-T-cell,",
    "Treg and neutrophil recruitment, activation, differentiation and",
    "proliferation, and IL-10 / IL-12 production. 21 ODE states: 3 PK",
    "states, 14 mucosal cell and cytokine states, an acute-injury",
    "stimulus (Activating_Event), a cumulative risk state (RISK), a",
    "tissue-exposure AUC and circulating natural Tregs. Outputs",
    "are the epithelial quality, the GI risk score (diarrhea onset when it",
    "exceeds 0.6) and RISK (colitis onset when it exceeds 200). Default",
    "parameters are the authors' baseline virtual patient (SimBiology",
    "variant 'VP0 (baseline VP)') with taselisib PK and potency; the six",
    "other inhibitors (idelalisib, duvelisib, umbralisib, alpelisib,",
    "pictilisib, copanlisib) differ only in the PK, protein-binding and",
    "IC50 parameters listed in population$notes. The published incidences",
    "come from a prevalence-weighted virtual population of 704 parameter",
    "sets, which is not carried here; the model is deterministic.",
    sep = " "
  )
  reference <- paste(
    "Gadkar K, Friedrich C, Hurez V, Ruiz M-L, Dickmann L, Jolly MK,",
    "Schutt L, Jin J, Ware JA, Ramanujan S. Quantitative systems",
    "pharmacology model-based investigation of adverse gastrointestinal",
    "events associated with prolonged treatment with PI3-kinase",
    "inhibitors. CPT Pharmacometrics Syst Pharmacol. 2022;11(5):616-627.",
    "doi:10.1002/psp4.12749.",
    "Equations: Supplement S10 (PSP4-11-616-s002.docx). Parameter values,",
    "dosing and per-drug variants: the authors' deposited SimBiology",
    "project GI_risk_model.sbproj and Parameters.xlsx (PSP4-11-616-s001.zip).",
    sep = " "
  )
  vignette <- "Gadkar_2022_pi3kiColitis_qsp"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Mucosal states keep the names of the deposited SimBiology species
  # (Supplement S10 'ODEs'); the drug states use the canonical PK names.
  paper_specific_compartments <- c(
    "T_Eff_act",
    "APC_act",
    "Treg_act",
    "APC_inact",
    "IL_12",
    "IL_10",
    "T_Eff_inact",
    "Treg_inact",
    "T_Eff_CKs",
    "T_quies",
    "Inf_CKs",
    "Activating_Event",
    "Neutrophils",
    "PI3K_Inh_Tissue_AUC",
    "EChealthy",
    "ECact",
    "RISK",
    "nTreg_circ"
  )

  compartmentData <- list(
    depot = list(analyte = "PI3K inhibitor", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "PI3K inhibitor", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "PI3K inhibitor", units = "mg", specimen = "tissue", verified = TRUE),
    T_Eff_act = list(
      analyte = "activated effector CD4 T cells",
      units = "cells/mm2",
      specimen = "tissue",
      verified = TRUE
    ),
    APC_act = list(
      analyte = "activated antigen-presenting cells",
      units = "cells/mm2",
      specimen = "tissue",
      verified = TRUE
    ),
    Treg_act = list(
      analyte = "activated regulatory T cells",
      units = "cells/mm2",
      specimen = "tissue",
      verified = TRUE
    ),
    APC_inact = list(
      analyte = "inactive antigen-presenting cells",
      units = "cells/mm2",
      specimen = "tissue",
      verified = TRUE
    ),
    IL_12 = list(analyte = "IL-12 (normalized to its EC50)", units = "AU", specimen = "tissue", verified = TRUE),
    IL_10 = list(analyte = "IL-10 (normalized to its EC50)", units = "AU", specimen = "tissue", verified = TRUE),
    T_Eff_inact = list(
      analyte = "inactive effector CD4 T cells",
      units = "cells/mm2",
      specimen = "tissue",
      verified = TRUE
    ),
    Treg_inact = list(
      analyte = "inactive regulatory T cells",
      units = "cells/mm2",
      specimen = "tissue",
      verified = TRUE
    ),
    T_Eff_CKs = list(
      analyte = "effector T-cell cytokines (normalized to their EC50)",
      units = "AU",
      specimen = "tissue",
      verified = TRUE
    ),
    T_quies = list(analyte = "quiescent CD4 T cells", units = "cells/mm2", specimen = "tissue", verified = TRUE),
    Inf_CKs = list(
      analyte = "inflammatory cytokines (normalized to their EC50)",
      units = "AU",
      specimen = "tissue",
      verified = TRUE
    ),
    Activating_Event = list(
      analyte = "acute epithelial injury stimulus",
      units = "AU",
      specimen = "not applicable",
      verified = TRUE
    ),
    Neutrophils = list(analyte = "neutrophils", units = "cells/mm2", specimen = "tissue", verified = TRUE),
    PI3K_Inh_Tissue_AUC = list(
      analyte = "cumulative unbound mucosal PI3K-inhibitor exposure",
      units = "ng*h/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    EChealthy = list(
      analyte = "healthy intestinal epithelial cells",
      units = "AU",
      specimen = "tissue",
      verified = TRUE
    ),
    ECact = list(analyte = "activated intestinal epithelial cells", units = "AU", specimen = "tissue", verified = TRUE),
    RISK = list(
      analyte = "cumulative GI risk (colitis driver)",
      units = "AU",
      specimen = "not applicable",
      verified = TRUE
    ),
    nTreg_circ = list(
      analyte = "circulating natural regulatory T cells",
      units = "AU",
      specimen = "whole blood",
      verified = TRUE
    )
  )

  covariateData <- list()

  population <- list(
    species = "human (virtual population)",
    n_subjects = 704L,
    n_studies = NA_integer_,
    disease_state = paste(
      "Patients with cancer (solid tumours or B-cell malignancies)",
      "receiving PI3K inhibitors. Calibrated to the published diarrhea",
      "and colitis incidences of six inhibitors and to individual onset",
      "times for taselisib (proprietary trial data: about 50% diarrhea,",
      "8-10% colitis, onset 3-200 and 80-200 days); copanlisib was held",
      "out for validation."
    ),
    dose_range = paste(
      "Deposited regimens: taselisib 4 mg QD, idelalisib 150 mg BID,",
      "duvelisib 25 mg BID, umbralisib 800 mg QD, alpelisib 300 mg QD,",
      "pictilisib 340 mg QD (all oral), copanlisib 60 mg IV on days 1, 8",
      "and 15 of a 28-day cycle."
    ),
    regions = NA_character_,
    notes = paste(
      "The 704 virtual patients are prevalence-weighted parameter sets",
      "varying 88 mechanistic parameters (authors' VPOP.mat); this file",
      "carries only the baseline virtual patient 'VP0 (baseline VP)'.",
      "Per-drug parameter sets from the deposited SimBiology variants",
      "(CL L/h, Vc L, Q L/h, Vp L, ka 1/h, lag h, bound fraction, and",
      "unbound IC50 alpha / beta / delta / gamma in ng/mL):",
      "taselisib 5, 200, 10, 50, 0.7, 0.2, 0.76, 0.13 / 4 / 0.055 / 0.45;",
      "idelalisib 14.88, 22.65, 11.82, 72.97, 0.482, 0.247, 0.84,",
      "340.64 / 234.71 / 1.04 / 36.97;",
      "duvelisib 5, 200, 10, 50, 0.7, 0.2, 0.7625,",
      "667.81 / 35.43 / 1.04 / 11.26;",
      "umbralisib (TGR-1202) 8.7, 550, 10, 160, 0.8, 0.2, 0.7625,",
      "101200 / 506 / 10.12 / 485;",
      "alpelisib 11.5, 118, 0 (one-compartment), 1, 0.784, 0.489, 0.9,",
      "2.21 / 510 / 110 / 128;",
      "pictilisib 34.12, 441.4, 8.145, 228.1, 2.065, 0.2, 0.95,",
      "1.38 / 37.17 / 0.83 / 22.02;",
      "copanlisib (IV bolus into central) 38.4, 317, 113, 815, -, -, 0.84,",
      "0.276 / 2.05 / 0.387 / 3.54.",
      "The deposit also carries a GDC-0077 (inavolisib) variant that the",
      "paper does not use: CL 5, Vc 200, Q 10, Vp 50, ka 0.7, lag 0.2,",
      "bound fraction 0.37, IC50 0.0184 / 53.36 / 6.3 / 10.3 ng/mL."
    )
  )

  ini({
    # -----------------------------------------------------------------
    # PI3K-inhibitor PK, binding and isoform potency. Defaults are the
    # taselisib values, which equal the 'VP0 (baseline VP)' variant of
    # the deposited SimBiology project (also Parameters.xlsx tab 'All
    # parameters'). Other drugs: population$notes and Table S5.
    # Deposited IC50s are unbound concentrations in ng/mL (Table S3 nM
    # times molecular weight / 1000; taselisib alpha 0.29 nM -> 0.13 ng/mL).
    # -----------------------------------------------------------------
    lka <- fixed(log(0.7)); label("Absorption rate constant (1/hr)")  # Table S5 taselisib ka 0.70 1/hr; PI3Kinh_Ka
    lcl <- fixed(log(5)); label("Clearance (L/hr)")  # Table S5 taselisib CL 5.0 L/hr; PI3Kinh_CL
    lvc <- fixed(log(200)); label("Central volume of distribution (L)")  # Table S5 taselisib Vc 200.0 L; PI3Kinh_Vc
    lq <- fixed(log(10)); label("Intercompartmental clearance (L/hr)")  # Table S5 taselisib Q 10.0 L/hr; PI3Kinh_Q
    lvp <- fixed(log(50)); label("Peripheral volume of distribution (L)")  # Table S5 taselisib Vp 50.0 L; PI3Kinh_Vp
    lfdepot <- fixed(log(1)); label("Fraction of the depot absorbed (unitless)")  # deposit PI3Kinh_F = 1 for every drug
    ltlag <- fixed(log(0.2)); label("Oral absorption lag time (hr)")  # deposit PI3Kinh_lag (dose LagParameterName)
    ld1 <- fixed(log(0.1)); label("Duration of the zero-order input into the depot (hr)")  # deposit PI3Kinh_duration (dose DurationParameterName)
    fu <- fixed(1 - 0.76); label("Unbound fraction in plasma (unitless)")  # 1 - PI3Kinh_bound_fraction; Table S5 taselisib f 0.76
    lkp_gut <- fixed(log(0.5)); label("Mucosa to unbound-plasma concentration ratio (unitless)")  # deposit PI3Kinh_partition (Mucosa.PI3K_Inhibitor rule)
    lic50_pi3ka <- fixed(log(0.13)); label("Unbound IC50 for PI3K-alpha (ng/mL)")  # deposit PI3Kinh_PI3Ka_IC50; Table S3 0.29 nM
    lic50_pi3kb <- fixed(log(4)); label("Unbound IC50 for PI3K-beta (ng/mL)")  # deposit PI3Kinh_PI3Kb_IC50; 30x alpha (Table S3 footnote)
    lic50_pi3kg <- fixed(log(0.45)); label("Unbound IC50 for PI3K-gamma (ng/mL)")  # deposit PI3Kinh_PI3Kg_IC50; Table S3 0.97 nM
    lic50_pi3kd <- fixed(log(0.055)); label("Unbound IC50 for PI3K-delta (ng/mL)")  # deposit PI3Kinh_PI3Kd_IC50; Table S3 0.12 nM
    imax_pi3ka <- fixed(1); label("Maximal inhibition of PI3K-alpha (fraction)")  # deposit PI3Kinh_PI3Ka_Imax
    imax_pi3kb <- fixed(1); label("Maximal inhibition of PI3K-beta (fraction)")  # deposit PI3Kinh_PI3Kb_Imax
    imax_pi3kg <- fixed(1); label("Maximal inhibition of PI3K-gamma (fraction)")  # deposit PI3Kinh_PI3Kg_Imax
    imax_pi3kd <- fixed(1); label("Maximal inhibition of PI3K-delta (fraction)")  # deposit PI3Kinh_PI3Kd_Imax
    Epi_act_event_lag <- fixed(0.2); label("Lag time of an acute-injury dose into Activating_Event (hr)")  # deposit dose Epi_activation_single
    Epi_act_event_dur <- fixed(0.1); label("Duration of an acute-injury dose into Activating_Event (hr)")  # deposit dose Epi_activation_single

    # -----------------------------------------------------------------
    # Mucosal immunology and epithelium. Every value below is the
    # deposited baseline virtual patient (SimBiology variant 'VP0
    # (baseline VP)' = Parameters.xlsx tab 'All parameters', in that
    # tab's row order); rate constants are per hour. Rounded values of
    # the cell half-lives, cytokine clearances (2.77 /h = 15-min
    # half-life) and baseline densities appear in Table S1.
    # -----------------------------------------------------------------
    Tquies_per_APC_diff          <- fixed(1); label("Quiescent T cells served per activated APC for differentiation (cells/cell)")
    TEff_per_APC_activation      <- fixed(10); label("Effector T cells served per activated APC for activation (cells/cell)")
    Treg_per_APC_activation      <- fixed(10); label("Tregs served per activated APC for activation (cells/cell)")
    TEff_per_APC_prolif          <- fixed(5); label("Effector T cells served per activated APC for proliferation (cells/cell)")
    Treg_per_APC_prolif          <- fixed(5); label("Tregs served per activated APC for proliferation (cells/cell)")
    Epi_Qual_Colitis_w           <- fixed(0.2); label("Weight of epithelial damage (1 - EQ) in the GI risk score (unitless)")
    APC_recruitment_rate_k       <- fixed(0.192068635784105); label("APC recruitment rate constant from blood (cells/mm2/hr per unit source)")
    InfCK_APC_recruit_EC50       <- fixed(1); label("Inflammatory cytokines for half-maximal stimulation of APC recruitment (normalized AU)")
    InfCK_APC_recruit_Emax       <- fixed(20); label("Maximal stimulation of APC recruitment by inflammatory cytokines (unitless)")
    PI3Kg_APC_recruit_EC50       <- fixed(1); label("PI3K-gamma activity for half-maximal stimulation of APC recruitment (activity, unitless)")
    PI3Kg_APC_recruit_Emax       <- fixed(2); label("Maximal stimulation of APC recruitment by PI3K-gamma activity (unitless)")
    TEffCK_APC_recruit_EC50      <- fixed(1); label("Effector T-cell cytokines for half-maximal stimulation of APC recruitment (normalized AU)")
    TEffCK_APC_recruit_Emax      <- fixed(20); label("Maximal stimulation of APC recruitment by effector T-cell cytokines (unitless)")
    InfCK_Tquies_recruit_EC50    <- fixed(1); label("Inflammatory cytokines for half-maximal stimulation of quiescent T-cell recruitment (normalized AU)")
    InfCK_Tquies_recruit_Emax    <- fixed(3); label("Maximal stimulation of quiescent T-cell recruitment by inflammatory cytokines (unitless)")
    PI3Kg_Tquies_recruit_EC50    <- fixed(1); label("PI3K-gamma activity for half-maximal stimulation of quiescent T-cell recruitment (activity, unitless)")
    PI3Kg_Tquies_recruit_Emax    <- fixed(2); label("Maximal stimulation of quiescent T-cell recruitment by PI3K-gamma activity (unitless)")
    TEffCK_Tquies_recruit_EC50   <- fixed(1); label("Effector T-cell cytokines for half-maximal stimulation of quiescent T-cell recruitment (normalized AU)")
    TEffCK_Tquies_recruit_Emax   <- fixed(3); label("Maximal stimulation of quiescent T-cell recruitment by effector T-cell cytokines (unitless)")
    Tquies_recruitment_rate_k    <- fixed(0.461535337124289); label("Quiescent CD4 T-cell recruitment rate constant from blood (cells/mm2/hr per unit source)")
    APC_activation_rate_k        <- fixed(0.0113778992440558); label("APC activation rate constant (1/hr)")
    Treg_APC_act_Imax            <- fixed(0.9); label("Maximal inhibition of APC activation by activated Tregs (fraction)")
    Treg_APC_act_IC50            <- fixed(20); label("Activated Treg density for half-maximal inhibition of APC activation (cells/mm2)")
    APCact_IL10_prod_rate_k      <- fixed(0.01136); label("IL-10 production rate constant per activated APC (AU/cell/hr)")
    PI3Kd_APCIL10_prod_rate_EC50 <- fixed(1); label("PI3K-delta activity for half-maximal stimulation of IL-10 production by activated APCs (activity, unitless)")
    PI3Kd_APCIL10_prod_rate_Emax <- fixed(3); label("Maximal stimulation of IL-10 production by activated APCs by PI3K-delta activity (unitless)")
    Treg_IL10_prod_rate_k        <- fixed(0.0172); label("IL-10 production rate constant per activated Treg (AU/cell/hr)")
    APCact_IL12_prod_rate_k      <- fixed(0.024); label("IL-12 production rate constant per activated APC (AU/cell/hr)")
    PI3Kd_APCIL12_prod_rate_Imax <- fixed(0.75); label("Maximal inhibition of IL-12 production by activated APCs by PI3K-delta activity (fraction)")
    PI3Kd_APCIL12_prod_rate_IC50 <- fixed(1); label("PI3K-delta activity for half-maximal inhibition of IL-12 production by activated APCs (activity, unitless)")
    PI3Kg_TEff_activation_EC50   <- fixed(1); label("PI3K-gamma activity for half-maximal stimulation of effector T-cell activation (activity, unitless)")
    PI3Kg_TEff_activation_Emax   <- fixed(1); label("Maximal stimulation of effector T-cell activation by PI3K-gamma activity (unitless)")
    TEffCK_TEff_activation_EC50  <- fixed(0.5); label("Effector T-cell cytokines for half-maximal stimulation of effector T-cell activation (normalized AU)")
    TEffCK_TEff_activation_Emax  <- fixed(20); label("Maximal stimulation of effector T-cell activation by effector T-cell cytokines (unitless)")
    T_Eff_act_rate_k             <- fixed(0.000119342254855249); label("Effector T-cell activation rate constant (1/hr)")
    IL10_Treg_activation_EC50    <- fixed(0.5); label("IL-10 for half-maximal stimulation of Treg activation (normalized AU)")
    IL10_Treg_activation_Emax    <- fixed(20); label("Maximal stimulation of Treg activation by IL-10 (unitless)")
    PI3Kd_Treg_activation_EC50   <- fixed(1); label("PI3K-delta activity for half-maximal stimulation of Treg activation (activity, unitless)")
    PI3Kd_Treg_activation_Emax   <- fixed(1); label("Maximal stimulation of Treg activation by PI3K-delta activity (unitless)")
    Treg_act_rate_k              <- fixed(6.41512998007073e-05); label("Treg activation rate constant (1/hr)")
    PI3Kd_TEff_prolif_EC50       <- fixed(1); label("PI3K-delta activity for half-maximal stimulation of effector T-cell proliferation (activity, unitless)")
    PI3Kd_TEff_prolif_Emax       <- fixed(0.5); label("Maximal stimulation of effector T-cell proliferation by PI3K-delta activity (unitless)")
    PI3Kg_TEff_prolif_EC50       <- fixed(1); label("PI3K-gamma activity for half-maximal stimulation of effector T-cell proliferation (activity, unitless)")
    PI3Kg_TEff_prolif_Emax       <- fixed(1); label("Maximal stimulation of effector T-cell proliferation by PI3K-gamma activity (unitless)")
    TEffCK_TEff_prolif_EC50      <- fixed(1); label("Effector T-cell cytokines for half-maximal stimulation of effector T-cell proliferation (normalized AU)")
    TEffCK_TEff_prolif_Emax      <- fixed(20); label("Maximal stimulation of effector T-cell proliferation by effector T-cell cytokines (unitless)")
    T_Eff_prolif_rate_k          <- fixed(3.68663234689859e-05); label("Effector T-cell proliferation rate constant (1/hr)")
    IL_10_clearance_rate_k       <- fixed(2.77258872223978); label("IL-10 clearance rate constant (1/hr; 15-min half-life)")
    IL_12_clearance_rate_k       <- fixed(2.77258872223978); label("IL-12 clearance rate constant (1/hr; 15-min half-life)")
    TEffAct_clearance_rate_k     <- fixed(0.008); label("Activated effector T-cell clearance rate constant (1/hr)")
    TEffCK_prod_rate_k           <- fixed(0.0338); label("Effector T-cell cytokine production rate constant per activated effector T cell (AU/cell/hr)")
    T_eff_CK_clear_rate_k        <- fixed(2.77258872223978); label("Effector T-cell cytokine clearance rate constant (1/hr; 15-min half-life)")
    IL10_TEff_diff_IC50          <- fixed(1); label("IL-10 for half-maximal inhibition of effector T-cell differentiation (normalized AU)")
    IL10_TEff_diff_Imax          <- fixed(0.5); label("Maximal inhibition of effector T-cell differentiation by IL-10 (fraction)")
    IL12_TEff_diff_EC50          <- fixed(1); label("IL-12 for half-maximal stimulation of effector T-cell differentiation (normalized AU)")
    IL12_TEff_diff_Emax          <- fixed(20); label("Maximal stimulation of effector T-cell differentiation by IL-12 (unitless)")
    PI3Kd_TEff_diff_EC50         <- fixed(1); label("PI3K-delta activity for half-maximal stimulation of effector T-cell differentiation (activity, unitless)")
    PI3Kd_TEff_diff_Emax         <- fixed(0.5); label("Maximal stimulation of effector T-cell differentiation by PI3K-delta activity (unitless)")
    PI3Kg_TEff_diff_EC50         <- fixed(1); label("PI3K-gamma activity for half-maximal stimulation of effector T-cell differentiation (activity, unitless)")
    PI3Kg_TEff_diff_Emax         <- fixed(1); label("Maximal stimulation of effector T-cell differentiation by PI3K-gamma activity (unitless)")
    T_Eff_diff_rate_k            <- fixed(0.000108189414026834); label("Quiescent-to-effector T-cell differentiation rate constant (1/hr)")
    Treg_diff_rate_k             <- fixed(6.32621016963091e-05); label("Quiescent-to-Treg differentiation rate constant (1/hr)")
    IL10_Treg_diff_EC50          <- fixed(1); label("IL-10 for half-maximal stimulation of Treg differentiation (normalized AU)")
    IL10_Treg_diff_Emax          <- fixed(20); label("Maximal stimulation of Treg differentiation by IL-10 (unitless)")
    IL12_Treg_diff_IC50          <- fixed(1); label("IL-12 for half-maximal inhibition of Treg differentiation (normalized AU)")
    IL12_Treg_diff_Imax          <- fixed(0.5); label("Maximal inhibition of Treg differentiation by IL-12 (fraction)")
    PI3Kd_Treg_diff_EC50         <- fixed(1); label("PI3K-delta activity for half-maximal stimulation of Treg differentiation (activity, unitless)")
    PI3Kd_Treg_diff_Emax         <- fixed(1); label("Maximal stimulation of Treg differentiation by PI3K-delta activity (unitless)")
    PI3Kg_Treg_diff_EC50         <- fixed(1); label("PI3K-gamma activity for half-maximal stimulation of Treg differentiation (activity, unitless)")
    PI3Kg_Treg_diff_Emax         <- fixed(1); label("Maximal stimulation of Treg differentiation by PI3K-gamma activity (unitless)")
    PI3Kb_EC_prolif_EC50         <- fixed(1); label("PI3K-beta inhibition (1 - activity) for half-maximal reduction of epithelial-cell formation (unitless)")
    PI3Kb_EC_prolif_Emax         <- fixed(2); label("Maximal reduction of epithelial-cell formation by PI3K-beta inhibition (fraction)")
    PI3Ka_ECquies_apop_Imax      <- fixed(0.5); label("Maximal inhibition of healthy epithelial-cell apoptosis by PI3K-alpha activity (fraction)")
    PI3Ka_ECquies_apop_IC50      <- fixed(1); label("PI3K-alpha activity for half-maximal inhibition of healthy epithelial-cell apoptosis (activity, unitless)")
    Epi_quies_apoptosis_rate_k   <- fixed(0.0101337307099407); label("Healthy epithelial-cell apoptosis rate constant (1/hr)")
    APCact_clearance_rate_k      <- fixed(0.0289); label("Activated APC clearance rate constant (1/hr)")
    APCinact_clearance_rate_k    <- fixed(0.002); label("Inactive APC clearance rate constant (1/hr)")
    Tquies_clearance_rate_k      <- fixed(0.002); label("Quiescent T-cell clearance rate constant (1/hr)")
    TEffInact_clearance_rate_k   <- fixed(0.002); label("Inactive effector T-cell clearance rate constant (1/hr)")
    TRegInact_clearance_rate_k   <- fixed(0.002); label("Inactive Treg clearance rate constant (1/hr)")
    TRegAct_clearance_rate_k     <- fixed(0.008); label("Activated Treg clearance rate constant (1/hr)")
    IL10_Treg_prolif_EC50        <- fixed(1); label("IL-10 for half-maximal stimulation of Treg proliferation (normalized AU)")
    IL10_Treg_prolif_Emax        <- fixed(20); label("Maximal stimulation of Treg proliferation by IL-10 (unitless)")
    Treg_prolif_rate_k           <- fixed(1.77972984907436e-05); label("Treg proliferation rate constant (1/hr)")
    InfCK_nTreg_recruit_EC50     <- fixed(1); label("Inflammatory cytokines for half-maximal stimulation of natural Treg recruitment (normalized AU)")
    InfCK_nTreg_recruit_Emax     <- fixed(20); label("Maximal stimulation of natural Treg recruitment by inflammatory cytokines (unitless)")
    TEffCK_nTreg_recruit_EC50    <- fixed(1); label("Effector T-cell cytokines for half-maximal stimulation of natural Treg recruitment (normalized AU)")
    TEffCK_nTreg_recruit_Emax    <- fixed(20); label("Maximal stimulation of natural Treg recruitment by effector T-cell cytokines (unitless)")
    nTreg_recruitment_rate_k     <- fixed(0.000811374504936677); label("Natural Treg recruitment rate constant from blood (cells/mm2/hr per unit source)")
    nTreg_production_rate_k      <- fixed(0.006); label("Circulating natural Treg production rate (AU/hr)")
    TEffCK_EC_act_Emax           <- fixed(10); label("Maximal stimulation of epithelial-cell activation by effector T-cell cytokines (unitless)")
    InfCK_EC_act_Emax            <- fixed(10); label("Maximal stimulation of epithelial-cell activation by inflammatory cytokines (unitless)")
    InfCK_EC_act_EC50            <- fixed(1); label("Inflammatory cytokines for half-maximal stimulation of epithelial-cell activation (normalized AU)")
    TEffCK_EC_act_diff_EC50      <- fixed(1); label("Effector T-cell cytokines for half-maximal stimulation of epithelial-cell activation (normalized AU)")
    APCinfCK_prod_rate_k         <- fixed(0.01); label("Inflammatory cytokine production rate constant per activated APC (AU/cell/hr)")
    Inf_CK_clearance_rate_k      <- fixed(2.77258872223978); label("Inflammatory cytokine clearance rate constant (1/hr; 15-min half-life)")
    Epi_act_apoptosis_rate_k     <- fixed(0.0506686535497036); label("Activated epithelial-cell apoptosis rate constant (1/hr)")
    PI3Ka_ECact_apop_Imax        <- fixed(0.5); label("Maximal inhibition of activated epithelial-cell apoptosis by PI3K-alpha activity (fraction)")
    PI3Ka_ECact_apop_IC50        <- fixed(1); label("PI3K-alpha activity for half-maximal inhibition of activated epithelial-cell apoptosis (activity, unitless)")
    nTreg_circ_clearance_rate_k  <- fixed(0.002); label("Circulating natural Treg clearance rate constant (1/hr)")
    PI3Kd_nTreg_prod_EC50        <- fixed(1); label("PI3K-delta activity for half-maximal stimulation of natural Treg production (blood) (activity, unitless)")
    PI3Kd_nTreg_prod_Emax        <- fixed(2); label("Maximal stimulation of natural Treg production (blood) by PI3K-delta activity (unitless)")
    ECact_ECprolif_Emax          <- fixed(1); label("Maximal stimulation of epithelial-cell formation by activated epithelial cells (unitless)")
    Neutro_clearance_rate_k      <- fixed(0.03); label("Mucosal neutrophil clearance rate constant (1/hr)")
    InfCK_Neutro_recruit_EC50    <- fixed(0.5); label("Inflammatory cytokines for half-maximal stimulation of neutrophil recruitment (normalized AU)")
    InfCK_Neutro_recruit_Emax    <- fixed(80); label("Maximal stimulation of neutrophil recruitment by inflammatory cytokines (unitless)")
    Neutro_recruitment_rate_k    <- fixed(0.0720103769372933); label("Neutrophil recruitment rate constant from blood (cells/mm2/hr per unit source)")
    TEffCK_Neutro_recruit_EC50   <- fixed(1); label("Effector T-cell cytokines for half-maximal stimulation of neutrophil recruitment (normalized AU)")
    TEffCK_Neutro_recruit_Emax   <- fixed(10); label("Maximal stimulation of neutrophil recruitment by effector T-cell cytokines (unitless)")
    Neutro_EC_act_EC50           <- fixed(50); label("Neutrophil density for half-maximal stimulation of epithelial-cell activation (cells/mm2)")
    Neutro_EC_act_Emax           <- fixed(1); label("Maximal stimulation of epithelial-cell activation by neutrophils (unitless)")
    PI3Kg_Neutro_recruit_EC50    <- fixed(1); label("PI3K-gamma activity for half-maximal stimulation of neutrophil recruitment (activity, unitless)")
    PI3Kg_Neutro_recruit_Emax    <- fixed(1); label("Maximal stimulation of neutrophil recruitment by PI3K-gamma activity (unitless)")
    PI3Kg_nTreg_recruit_EC50     <- fixed(1); label("PI3K-gamma activity for half-maximal stimulation of natural Treg recruitment (activity, unitless)")
    PI3Kg_nTreg_recruit_Emax     <- fixed(2); label("Maximal stimulation of natural Treg recruitment by PI3K-gamma activity (unitless)")
    EQ_APC_act_IC50              <- fixed(0.5); label("Epithelial quality for half-maximal inhibition of APC activation (unitless)")
    EQ_APC_act_Imax              <- fixed(1); label("Maximal inhibition of APC activation by epithelial quality (fraction)")
    EQ_APC_act_nH                <- fixed(5); label("Hill coefficient of epithelial-quality inhibition of APC activation (unitless)")
    InfCK_APC_act_EC50           <- fixed(1); label("Inflammatory cytokines for half-maximal stimulation of APC activation (normalized AU)")
    InfCK_APC_act_Emax           <- fixed(20); label("Maximal stimulation of APC activation by inflammatory cytokines (unitless)")
    ECact_ECprolif_EC50          <- fixed(0.1); label("Activated epithelial-cell level for half-maximal stimulation of epithelial-cell formation (AU)")
    APC_TEff_diff_Emax           <- fixed(1); label("Maximal stimulation of effector T-cell differentiation by activated APCs (unitless)")
    APC_Treg_diff_Emax           <- fixed(1); label("Maximal stimulation of Treg differentiation by activated APCs (unitless)")
    APC_TEff_prolif_Emax         <- fixed(1); label("Maximal stimulation of effector T-cell proliferation by activated APCs (unitless)")
    APC_TEff_act_Emax            <- fixed(1); label("Maximal stimulation of effector T-cell activation by activated APCs (unitless)")
    APC_Treg_prolif_Emax         <- fixed(1); label("Maximal stimulation of Treg proliferation by activated APCs (unitless)")
    APC_Treg_act_Emax            <- fixed(1); label("Maximal stimulation of Treg activation by activated APCs (unitless)")
    NeutroInfCK_prod_rate_k      <- fixed(0.005); label("Inflammatory cytokine production rate constant per neutrophil (AU/cell/hr)")
    InfCK_TEff_prolif_EC50       <- fixed(1); label("Inflammatory cytokines for half-maximal stimulation of effector T-cell proliferation (normalized AU)")
    InfCK_TEff_prolif_Emax       <- fixed(20); label("Maximal stimulation of effector T-cell proliferation by inflammatory cytokines (unitless)")
    InfCK_TEff_activation_EC50   <- fixed(0.5); label("Inflammatory cytokines for half-maximal stimulation of effector T-cell activation (normalized AU)")
    InfCK_TEff_activation_Emax   <- fixed(20); label("Maximal stimulation of effector T-cell activation by inflammatory cytokines (unitless)")
    Treg_TEff_prolif_IC50        <- fixed(0.1); label("Activated-Treg to inactive-effector ratio for half-maximal inhibition of effector T-cell proliferation (unitless)")
    Treg_TEff_prolif_Imax        <- fixed(0.7); label("Maximal inhibition of effector T-cell proliferation by Tregs (fraction)")
    PI3Kg_Treg_prolif_EC50       <- fixed(1); label("PI3K-gamma activity for half-maximal stimulation of Treg proliferation (activity, unitless)")
    PI3Kg_Treg_prolif_Emax       <- fixed(1); label("Maximal stimulation of Treg proliferation by PI3K-gamma activity (unitless)")
    PI3Kg_Treg_activation_EC50   <- fixed(1); label("PI3K-gamma activity for half-maximal stimulation of Treg activation (activity, unitless)")
    PI3Kg_Treg_activation_Emax   <- fixed(1); label("Maximal stimulation of Treg activation by PI3K-gamma activity (unitless)")
    Treg_TEff_act_IC50           <- fixed(0.1); label("Activated-Treg to inactive-effector ratio for half-maximal inhibition of effector T-cell activation (unitless)")
    Treg_TEff_act_Imax           <- fixed(0.7); label("Maximal inhibition of effector T-cell activation by Tregs (fraction)")
    InfCK_Neutro_recruit_nH      <- fixed(2); label("Hill coefficient of inflammatory-cytokine stimulation of neutrophil recruitment (unitless)")
    Act_event_clear_rate_k       <- fixed(0.1); label("Clearance rate constant of the acute-injury activating event (1/hr)")
    PI3Kd_TEff_activation_EC50   <- fixed(1); label("PI3K-delta activity for half-maximal stimulation of effector T-cell activation (activity, unitless)")
    PI3Kd_TEff_activation_Emax   <- fixed(0.5); label("Maximal stimulation of effector T-cell activation by PI3K-delta activity (unitless)")
    PI3Kd_Treg_prolif_EC50       <- fixed(1); label("PI3K-delta activity for half-maximal stimulation of Treg proliferation (activity, unitless)")
    PI3Kd_Treg_prolif_Emax       <- fixed(1); label("Maximal stimulation of Treg proliferation by PI3K-delta activity (unitless)")
    PI3Ka_EC_prolif_EC50         <- fixed(1); label("PI3K-alpha inhibition (1 - activity) for half-maximal reduction of epithelial-cell formation (unitless)")
    PI3Ka_EC_prolif_Emax         <- fixed(0); label("Maximal reduction of epithelial-cell formation by PI3K-alpha inhibition (fraction)")
    kEC_activation               <- fixed(5e-04); label("Healthy-to-activated epithelial-cell transition rate constant (1/hr)")
    kEC_formation                <- fixed(0.007); label("Zero-order healthy epithelial-cell formation rate (AU/hr)")
    Tact_GI                      <- fixed(150); label("Activated effector T-cell reference density in the GI risk score (cells/mm2)")
    APCact_GI                    <- fixed(250); label("Activated APC reference density in the GI risk score (cells/mm2)")
    neutroGI                     <- fixed(250); label("Neutrophil reference density in the GI risk score (cells/mm2)")
    threshold                    <- fixed(0.25); label("GI risk score below which the cumulative risk state heals (unitless)")
    k_healing                    <- fixed(0.05); label("Accumulation and healing rate constant of the cumulative risk state RISK (1/hr)")
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)
    fdepot <- exp(lfdepot)
    tlag <- exp(ltlag)
    d1 <- exp(ld1)
    kp_gut <- exp(lkp_gut)
    ic50_pi3ka <- exp(lic50_pi3ka)
    ic50_pi3kb <- exp(lic50_pi3kb)
    ic50_pi3kg <- exp(lic50_pi3kg)
    ic50_pi3kd <- exp(lic50_pi3kd)

    Cc <- 1000 * central / vc
    Cu <- fu * Cc
    Cmucosa <- kp_gut * Cu

    PI3Ka_act <- 1 - imax_pi3ka * Cmucosa / (ic50_pi3ka + Cmucosa)
    PI3Kb_act <- 1 - imax_pi3kb * Cmucosa / (ic50_pi3kb + Cmucosa)
    PI3Kg_act <- 1 - imax_pi3kg * Cmucosa / (ic50_pi3kg + Cmucosa)
    PI3Kd_act <- 1 - imax_pi3kd * Cmucosa / (ic50_pi3kd + Cmucosa)
    PI3Kg_act_circ <- 1 - imax_pi3kg * Cu / (ic50_pi3kg + Cu)
    PI3Kd_act_circ <- 1 - imax_pi3kd * Cu / (ic50_pi3kd + Cu)

    APC_circ <- 1
    T_quies_circ <- 1
    Neutro_circ <- 1

    Epithelial_Quality <- min(1, EChealthy + 0.2 * ECact)
    GI_Risk_Score <- Epi_Qual_Colitis_w * (1 - Epithelial_Quality) +
      (1 - Epi_Qual_Colitis_w) * (T_Eff_act / Tact_GI + APC_act / APCact_GI + Neutrophils / neutroGI) / 3
    APC_Demand <- T_quies / Tquies_per_APC_diff + T_Eff_inact / TEff_per_APC_activation +
      Treg_inact / Treg_per_APC_activation + T_Eff_inact / TEff_per_APC_prolif +
      Treg_inact / Treg_per_APC_prolif
    apc_ratio <- APC_act / APC_Demand

    APC_recruitment <- APC_recruitment_rate_k * APC_circ *
      (1 + TEffCK_APC_recruit_Emax * T_Eff_CKs / (TEffCK_APC_recruit_EC50 + T_Eff_CKs) +
        InfCK_APC_recruit_Emax * Inf_CKs / (InfCK_APC_recruit_EC50 + Inf_CKs)) *
      (1 + PI3Kg_APC_recruit_Emax * PI3Kg_act_circ / (PI3Kg_APC_recruit_EC50 + PI3Kg_act_circ))
    T_cell_recruitment <- Tquies_recruitment_rate_k * T_quies_circ *
      (1 + TEffCK_Tquies_recruit_Emax * T_Eff_CKs / (TEffCK_Tquies_recruit_EC50 + T_Eff_CKs) +
        InfCK_Tquies_recruit_Emax * Inf_CKs / (InfCK_Tquies_recruit_EC50 + Inf_CKs)) *
      (1 + PI3Kg_Tquies_recruit_Emax * PI3Kg_act_circ / (PI3Kg_Tquies_recruit_EC50 + PI3Kg_act_circ))
    APC_activation_rate <- APC_activation_rate_k * APC_inact *
      (1 + InfCK_APC_act_Emax * Inf_CKs / (InfCK_APC_act_EC50 + Inf_CKs)) *
      (1 - Treg_APC_act_Imax * Treg_act / (Treg_APC_act_IC50 + Treg_act)) *
      (1 - EQ_APC_act_Imax * Epithelial_Quality^EQ_APC_act_nH /
        (EQ_APC_act_IC50^EQ_APC_act_nH + Epithelial_Quality^EQ_APC_act_nH))
    IL10_production_rate <- Treg_IL10_prod_rate_k * Treg_act + APCact_IL10_prod_rate_k * APC_act *
      (1 + PI3Kd_APCIL10_prod_rate_Emax * PI3Kd_act / (PI3Kd_APCIL10_prod_rate_EC50 + PI3Kd_act))
    IL12_production_rate <- APCact_IL12_prod_rate_k * APC_act *
      (1 - PI3Kd_APCIL12_prod_rate_Imax * PI3Kd_act / (PI3Kd_APCIL12_prod_rate_IC50 + PI3Kd_act))
    T_Eff_activation_rate <- T_Eff_act_rate_k * T_Eff_inact *
      (1 + APC_TEff_act_Emax * apc_ratio / (1 / TEff_per_APC_activation + apc_ratio)) *
      (1 + TEffCK_TEff_activation_Emax * T_Eff_CKs / (TEffCK_TEff_activation_EC50 + T_Eff_CKs) +
        InfCK_TEff_activation_Emax * Inf_CKs / (InfCK_TEff_activation_EC50 + Inf_CKs)) *
      (1 + PI3Kg_TEff_activation_Emax * PI3Kg_act / (PI3Kg_TEff_activation_EC50 + PI3Kg_act)) *
      (1 + PI3Kd_TEff_activation_Emax * PI3Kd_act / (PI3Kd_TEff_activation_EC50 + PI3Kd_act)) *
      (1 - Treg_TEff_act_Imax * (Treg_act / T_Eff_inact) / (Treg_TEff_act_IC50 + Treg_act / T_Eff_inact))
    Treg_activation_rate <- Treg_act_rate_k * Treg_inact *
      (1 + APC_Treg_act_Emax * apc_ratio / (1 / Treg_per_APC_activation + apc_ratio)) *
      (1 + IL10_Treg_activation_Emax * IL_10 / (IL10_Treg_activation_EC50 + IL_10)) *
      (1 + PI3Kd_Treg_activation_Emax * PI3Kd_act / (PI3Kd_Treg_activation_EC50 + PI3Kd_act)) *
      (1 + PI3Kg_Treg_activation_Emax * PI3Kg_act / (PI3Kg_Treg_activation_EC50 + PI3Kg_act))
    T_Eff_proliferation_rate <- T_Eff_prolif_rate_k * T_Eff_inact *
      (1 + APC_TEff_prolif_Emax * apc_ratio / (1 / TEff_per_APC_prolif + apc_ratio)) *
      (1 + TEffCK_TEff_prolif_Emax * T_Eff_CKs / (TEffCK_TEff_prolif_EC50 + T_Eff_CKs) +
        InfCK_TEff_prolif_Emax * Inf_CKs / (InfCK_TEff_prolif_EC50 + Inf_CKs)) *
      (1 + PI3Kg_TEff_prolif_Emax * PI3Kg_act / (PI3Kg_TEff_prolif_EC50 + PI3Kg_act)) *
      (1 + PI3Kd_TEff_prolif_Emax * PI3Kd_act / (PI3Kd_TEff_prolif_EC50 + PI3Kd_act)) *
      (1 - Treg_TEff_prolif_Imax * (Treg_act / T_Eff_inact) / (Treg_TEff_prolif_IC50 + Treg_act / T_Eff_inact))
    T_Eff_differentiation_rate <- T_Eff_diff_rate_k * T_quies *
      (1 + APC_TEff_diff_Emax * apc_ratio / (1 / Tquies_per_APC_diff + apc_ratio)) *
      (1 + PI3Kg_TEff_diff_Emax * PI3Kg_act / (PI3Kg_TEff_diff_EC50 + PI3Kg_act) +
        PI3Kd_TEff_diff_Emax * PI3Kd_act / (PI3Kd_TEff_diff_EC50 + PI3Kd_act) +
        IL12_TEff_diff_Emax * IL_12 / (IL12_TEff_diff_EC50 + IL_12)) *
      (1 - IL10_TEff_diff_Imax * IL_10 / (IL10_TEff_diff_IC50 + IL_10))
    Treg_differentiation_rate <- Treg_diff_rate_k * T_quies *
      (1 + APC_Treg_diff_Emax * apc_ratio / (1 / Tquies_per_APC_diff + apc_ratio)) *
      (1 + IL10_Treg_diff_Emax * IL_10 / (IL10_Treg_diff_EC50 + IL_10)) *
      (1 + PI3Kg_Treg_diff_Emax * PI3Kg_act / (PI3Kg_Treg_diff_EC50 + PI3Kg_act)) *
      (1 + PI3Kd_Treg_diff_Emax * PI3Kd_act / (PI3Kd_Treg_diff_EC50 + PI3Kd_act)) *
      (1 - IL12_Treg_diff_Imax * IL_12 / (IL12_Treg_diff_IC50 + IL_12))
    Treg_proliferation_rate <- Treg_prolif_rate_k * Treg_inact *
      (1 + APC_Treg_prolif_Emax * apc_ratio / (1 / Treg_per_APC_prolif + apc_ratio)) *
      (1 + IL10_Treg_prolif_Emax * IL_10 / (IL10_Treg_prolif_EC50 + IL_10)) *
      (1 + PI3Kg_Treg_prolif_Emax * PI3Kg_act / (PI3Kg_Treg_prolif_EC50 + PI3Kg_act)) *
      (1 + PI3Kd_Treg_prolif_Emax * PI3Kd_act / (PI3Kd_Treg_prolif_EC50 + PI3Kd_act))
    nTreg_recruitment_rate <- nTreg_recruitment_rate_k * nTreg_circ *
      (1 + TEffCK_nTreg_recruit_Emax * T_Eff_CKs / (TEffCK_nTreg_recruit_EC50 + T_Eff_CKs) +
        InfCK_nTreg_recruit_Emax * Inf_CKs / (InfCK_nTreg_recruit_EC50 + Inf_CKs)) *
      (1 + PI3Kg_nTreg_recruit_Emax * PI3Kg_act_circ / (PI3Kg_nTreg_recruit_EC50 + PI3Kg_act_circ))
    nTreg_production_rate <- nTreg_production_rate_k *
      (1 + PI3Kd_nTreg_prod_Emax * PI3Kd_act_circ / (PI3Kd_nTreg_prod_EC50 + PI3Kd_act_circ))
    Neutro_recruitment_rate <- Neutro_recruitment_rate_k * Neutro_circ *
      (1 + TEffCK_Neutro_recruit_Emax * T_Eff_CKs / (TEffCK_Neutro_recruit_EC50 + T_Eff_CKs) +
        InfCK_Neutro_recruit_Emax * Inf_CKs^InfCK_Neutro_recruit_nH /
          (InfCK_Neutro_recruit_EC50^InfCK_Neutro_recruit_nH + Inf_CKs^InfCK_Neutro_recruit_nH)) *
      (1 + PI3Kg_Neutro_recruit_Emax * PI3Kg_act_circ / (PI3Kg_Neutro_recruit_EC50 + PI3Kg_act_circ))
    Inf_CK_production_rate <- APCinfCK_prod_rate_k * APC_act + NeutroInfCK_prod_rate_k * Neutrophils
    ECactivation <- kEC_activation * EChealthy *
      (1 + TEffCK_EC_act_Emax * T_Eff_CKs / (TEffCK_EC_act_diff_EC50 + T_Eff_CKs) +
        InfCK_EC_act_Emax * Inf_CKs / (InfCK_EC_act_EC50 + Inf_CKs) +
        Neutro_EC_act_Emax * Neutrophils / (Neutro_EC_act_EC50 + Neutrophils) + Activating_Event)
    ECformation <- kEC_formation * (1 + ECact_ECprolif_Emax * ECact / (ECact_ECprolif_EC50 + ECact) *
      (1 - PI3Ka_EC_prolif_Emax * (1 - PI3Ka_act) / (PI3Ka_EC_prolif_EC50 + (1 - PI3Ka_act)) -
        PI3Kb_EC_prolif_Emax * (1 - PI3Kb_act) / (PI3Kb_EC_prolif_EC50 + (1 - PI3Kb_act))))
    EChealthy_apoptosis <- Epi_quies_apoptosis_rate_k * EChealthy *
      (1 - PI3Ka_ECquies_apop_Imax * PI3Ka_act / (PI3Ka_ECquies_apop_IC50 + PI3Ka_act))
    ECact_apoptosis <- Epi_act_apoptosis_rate_k * ECact *
      (1 - PI3Ka_ECact_apop_Imax * PI3Ka_act / (PI3Ka_ECact_apop_IC50 + PI3Ka_act))

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * fdepot * depot - cl / vc * central - (q / vc * central - q / vp * peripheral1)
    d/dt(peripheral1) <- q / vc * central - q / vp * peripheral1
    d/dt(T_Eff_act) <- T_Eff_activation_rate - TEffAct_clearance_rate_k * T_Eff_act
    d/dt(APC_act) <- APC_activation_rate - APCact_clearance_rate_k * APC_act
    d/dt(Treg_act) <- Treg_activation_rate - TRegAct_clearance_rate_k * Treg_act
    d/dt(APC_inact) <- APC_recruitment - APC_activation_rate - APCinact_clearance_rate_k * APC_inact
    d/dt(IL_12) <- IL12_production_rate - IL_12_clearance_rate_k * IL_12
    d/dt(IL_10) <- IL10_production_rate - IL_10_clearance_rate_k * IL_10
    d/dt(T_Eff_inact) <- -T_Eff_activation_rate + T_Eff_proliferation_rate + T_Eff_differentiation_rate -
      TEffInact_clearance_rate_k * T_Eff_inact
    d/dt(Treg_inact) <- -Treg_activation_rate + Treg_differentiation_rate - TRegInact_clearance_rate_k * Treg_inact +
      Treg_proliferation_rate + nTreg_recruitment_rate
    d/dt(T_Eff_CKs) <- TEffCK_prod_rate_k * T_Eff_act - T_eff_CK_clear_rate_k * T_Eff_CKs
    d/dt(T_quies) <- T_cell_recruitment - T_Eff_differentiation_rate - Treg_differentiation_rate -
      Tquies_clearance_rate_k * T_quies
    d/dt(Inf_CKs) <- Inf_CK_production_rate - Inf_CK_clearance_rate_k * Inf_CKs
    d/dt(Activating_Event) <- -Act_event_clear_rate_k * Activating_Event
    d/dt(Neutrophils) <- Neutro_recruitment_rate - Neutro_clearance_rate_k * Neutrophils
    d/dt(PI3K_Inh_Tissue_AUC) <- Cmucosa
    d/dt(EChealthy) <- -ECactivation + ECformation - EChealthy_apoptosis
    d/dt(ECact) <- ECactivation - ECact_apoptosis
    d/dt(RISK) <- GI_Risk_Score * k_healing - k_healing * RISK * (GI_Risk_Score < threshold)
    d/dt(nTreg_circ) <- -nTreg_recruitment_rate + nTreg_production_rate - nTreg_circ_clearance_rate_k * nTreg_circ

    alag(depot) <- tlag
    dur(depot) <- d1
    alag(Activating_Event) <- Epi_act_event_lag
    dur(Activating_Event) <- Epi_act_event_dur

    # Drug-free steady state of the baseline virtual patient. The deposit
    # starts from Table S1 densities with no epithelium (EChealthy =
    # ECact = 0) and runs two 240-day drug-free run-ins before dosing
    # (simulate_Vpop.m); these are the converged values of that run-in
    # (relative change < 1e-14 between 1e5 and 2e5 h). Activating_Event,
    # PI3K_Inh_Tissue_AUC and the PK states start at 0.
    T_Eff_act(0) <- 2.8148915
    APC_act(0) <- 5.9915551
    Treg_act(0) <- 5.7220501
    APC_inact(0) <- 353.61305
    IL_12(0) <- 0.032414951
    IL_10(0) <- 0.096869552
    T_Eff_inact(0) <- 42.864047
    Treg_inact(0) <- 66.839997
    T_Eff_CKs(0) <- 0.034315704
    T_quies(0) <- 417.11661
    Inf_CKs(0) <- 0.032433815
    Neutrophils(0) <- 6.0020155
    EChealthy(0) <- 0.97774269
    ECact(0) <- 0.022552759
    RISK(0) <- 0.021346745
    nTreg_circ(0) <- 2.0982378
  })
}
