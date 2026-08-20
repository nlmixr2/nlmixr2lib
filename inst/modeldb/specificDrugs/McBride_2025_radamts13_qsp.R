McBride_2025_radamts13_qsp <- function() {
  description <- "QSP. Mechanistic ADAMTS13-VWF-platelet model for congenital thrombotic thrombocytopenic purpura (cTTP), driven by a two-compartment PK model of recombinant ADAMTS13 (rADAMTS13, TAK-755). Mass-action binding of elongated (active) and globular VWF to ADAMTS13, rADAMTS13, extracellular hemoglobin and thrombospondin-1; VWF-platelet binding, platelet aggregate formation and ADAMTS13-mediated aggregate lysis; constant platelet synthesis with first-order loss. 26 ODEs. Predicts platelet count, ADAMTS13 activity and active-VWF fraction, and was used for the virtual clinical trial simulations submitted to FDA supporting the cTTP approval (McBride 2025 CPT Pharmacometrics Syst Pharmacol)."
  reference <- "McBride C, Jiang J, Zhang Z, Tolsma J, Patwari P, Mellgard B, Vakilynejad M, Bhattacharya I, Zhu AZX. Quantitative Systems Pharmacology Modeling of Platelet Responses to Recombinant ADAMTS13 in Patients With Congenital Thrombotic Thrombocytopenic Purpura. CPT Pharmacometrics Syst Pharmacol. 2025;14(9):1575-1585. doi:10.1002/psp4.70063"
  vignette <- "McBride_2025_radamts13_qsp"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # The model has no dosing-independent baseline built in: every ODE state
  # starts at the rxode2 default of 0 and relaxes to its own steady state.
  # Initial conditions are deliberately NOT set with `state(0) <-` so that
  # rxSolve(inits = ) keeps working. Run the system with no dosing for
  # >= 2000 h (Methods S1; the shipped Run_ADAMTS13_model.m uses 5000 h)
  # before the first dose, then start dosing from the resulting state
  # vector. The validation vignette demonstrates this lead-in.

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Allometric scaling of the rADAMTS13 PK on a 68.7 kg reference subject (Table S1 'Human PK parameters'): exponent 0.75 on CL and Q, 1.0 on Vc and Vp. Also converts the IU/kg dose to a delivered mg amount in the event table (see vignette).",
      source_name        = "WT"
    )
  )

  # 26 ODE states. The soluble species (VWF forms, ADAMTS13, hemoglobin,
  # TSP-1 and all their complexes) are plasma analytes in nM; the
  # platelet-containing species are whole-blood analytes in nM of platelet
  # equivalents (a platelet count of 1e9/L corresponds to 1/6.02e5 nM).
  # Verified against Methods S1 Eqs 4-25 and ADAMTS13_Model.m (Data S2).
  compartmentData <- list(
    vg                    = list(analyte = "von Willebrand factor, globular monomer unit", units = "nmol", specimen = "plasma", verified = TRUE),
    ve                    = list(analyte = "von Willebrand factor, elongated (active) monomer unit", units = "nmol", specimen = "plasma", verified = TRUE),
    adam                  = list(analyte = "endogenous ADAMTS13", units = "nmol", specimen = "plasma", verified = TRUE),
    hb                    = list(analyte = "extracellular hemoglobin", units = "nmol", specimen = "plasma", verified = TRUE),
    tsp1                  = list(analyte = "thrombospondin-1", units = "nmol", specimen = "plasma", verified = TRUE),
    vg_adam               = list(analyte = "globular VWF:ADAMTS13 complex", units = "nmol", specimen = "plasma", verified = TRUE),
    ve_adam               = list(analyte = "elongated VWF:ADAMTS13 complex", units = "nmol", specimen = "plasma", verified = TRUE),
    vg_hb                 = list(analyte = "globular VWF:hemoglobin complex", units = "nmol", specimen = "plasma", verified = TRUE),
    ve_hb                 = list(analyte = "elongated VWF:hemoglobin complex", units = "nmol", specimen = "plasma", verified = TRUE),
    ve_tsp1               = list(analyte = "elongated VWF:thrombospondin-1 complex", units = "nmol", specimen = "plasma", verified = TRUE),
    vg_frag               = list(analyte = "globular VWF cleavage fragments", units = "nmol", specimen = "plasma", verified = TRUE),
    ve_frag               = list(analyte = "elongated VWF cleavage fragments", units = "nmol", specimen = "plasma", verified = TRUE),
    vg_radam              = list(analyte = "globular VWF:rADAMTS13 complex", units = "nmol", specimen = "plasma", verified = TRUE),
    ve_radam              = list(analyte = "elongated VWF:rADAMTS13 complex", units = "nmol", specimen = "plasma", verified = TRUE),
    platelet              = list(analyte = "free platelets", units = "nmol", specimen = "whole blood", verified = TRUE),
    ve_platelet           = list(analyte = "elongated VWF:platelet complex", units = "nmol", specimen = "whole blood", verified = TRUE),
    ve_adam_platelet      = list(analyte = "elongated VWF:ADAMTS13:platelet complex", units = "nmol", specimen = "whole blood", verified = TRUE),
    ve_tsp1_platelet      = list(analyte = "elongated VWF:thrombospondin-1:platelet complex", units = "nmol", specimen = "whole blood", verified = TRUE),
    ve_hb_platelet        = list(analyte = "elongated VWF:hemoglobin:platelet complex", units = "nmol", specimen = "whole blood", verified = TRUE),
    ve_radam_platelet     = list(analyte = "elongated VWF:rADAMTS13:platelet complex", units = "nmol", specimen = "whole blood", verified = TRUE),
    ve_platelet_coag      = list(analyte = "elongated VWF:multi-platelet aggregate", units = "nmol", specimen = "whole blood", verified = TRUE),
    ve_tsp1_platelet_coag = list(analyte = "elongated VWF:thrombospondin-1:multi-platelet aggregate", units = "nmol", specimen = "whole blood", verified = TRUE),
    ve_hb_platelet_coag   = list(analyte = "elongated VWF:hemoglobin:multi-platelet aggregate", units = "nmol", specimen = "whole blood", verified = TRUE),
    radam                 = list(analyte = "recombinant ADAMTS13 (TAK-755) in the pharmacodynamic compartment", units = "nmol", specimen = "plasma", verified = TRUE),
    central               = list(analyte = "recombinant ADAMTS13 (TAK-755)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1           = list(analyte = "recombinant ADAMTS13 (TAK-755)", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 40,
    n_studies      = 5,
    age_range      = "adults and adolescents with cTTP (Phase 3 NCT03393975 enrolled patients >= 0 years; the QSP virtual population is the non-pregnant adult cTTP phenotype)",
    weight_median  = "68.7 kg (allometric reference weight, Table S1 'Human PK parameters')",
    disease_state  = "Congenital thrombotic thrombocytopenic purpura (cTTP), a severe inherited ADAMTS13 deficiency causing thrombotic microangiopathy",
    dose_range     = "rADAMTS13 40 IU/kg IV Q1W or Q2W prophylaxis; plasma-based therapy (PBT) exposure equivalent to ~10 IU/kg ADAMTS13 (10-20 mL/kg fresh-frozen or solvent/detergent-treated plasma) Q1W or Q2W",
    notes          = paste(
      "Calibration/validation data (Table S2 and Figure 2): Phase 1 cTTP PK (NCT02216084, Scully 2017, N = 15);",
      "two adamts13-knockout mouse studies (Kopic 2016 and in-house #248.220.5098, 5 male + 5 female mice per dose);",
      "human baseline platelet counts and VWF activity from healthy individuals and untreated cTTP patients;",
      "cTTP FFP standard-of-care literature; Phase 2 iTTP (NCT03922308, N = 16 rADAMTS13-treated);",
      "and the Phase 3 cTTP crossover study (NCT03393975, N = 48 enrolled; data from 40 patients used here,",
      "split 25 estimation / 15 validation). n_subjects = 40 refers to the Phase 3 patients that informed the",
      "population parameter variability. Clinical trial simulations used 1000 virtual patients per arm,",
      "phenotype-matched on body weight and baseline ADAMTS13 and platelet levels; the paper states that",
      "sex, age and height were not model inputs."
    ),
    virtual_population = paste(
      "Population variability was estimated by iterative two-stage (ITS) on seven individual parameters",
      "(ksyn_platelet, ve_frac_base, kd_ve_platelet, WT, spike_amount, vwf_pct, spike_day; Figure S2) and",
      "reported only as histograms - the mean vector and covariance matrix are not tabulated, so no IIV is",
      "encoded in this file. Two fully specified virtual patients are shipped in Data S2",
      "(ADAMTS13_Model.m) and are reproduced in the validation vignette."
    )
  )

  # Paper-mechanistic QSP species; none map onto a canonical nlmixr2lib
  # compartment role. Names follow the published Methods S1 / Table S1
  # nomenclature (VE = elongated/active VWF monomer unit, VG = globular
  # VWF monomer unit, PLT = platelet, coag = multi-platelet aggregate).
  paper_specific_compartments <- c(
    "vg", "ve", "adam", "hb", "tsp1",
    "vg_adam", "ve_adam", "vg_hb", "ve_hb", "ve_tsp1",
    "vg_frag", "ve_frag", "vg_radam", "ve_radam",
    "platelet", "ve_platelet", "ve_adam_platelet", "ve_tsp1_platelet",
    "ve_hb_platelet", "ve_radam_platelet",
    "ve_platelet_coag", "ve_tsp1_platelet_coag", "ve_hb_platelet_coag",
    "radam"
  )

  ini({
    # ------------------------------------------------------------------
    # rADAMTS13 (TAK-755) population PK - two compartment, IV infusion.
    # Table S1 "Human PK parameters" reports CL = 0.0398 and Q = 0.456
    # (x (WT/68.7)^0.75); the shipped model code (Data S2,
    # ADAMTS13_Model.m lines 75-78) uses CL = 0.0389 and CLD = 0.0456.
    # The code values are used here because Data S2 is the executable
    # artifact that produced the published figures. Note that the two
    # candidate Q values CANNOT be told apart from the on-disk data: both
    # land within ~10% of the Phase 3 observed mean ADAMTS13 Cmax
    # (Table S2: 101% rADAMTS13, 19% PBT). The ten-fold discrepancy is an
    # unresolved inconsistency in the published sources; see the vignette
    # Errata.
    # ------------------------------------------------------------------
    lcl     <- fixed(log(0.0389)); label("Clearance at the reference body weight (L/h)")                        # Data S2 ADAMTS13_Model.m L77 (Table S1 Human PK: 0.0398)
    lvc     <- fixed(log(2.69));   label("Central volume of distribution at the reference body weight (L)")     # Table S1 Human PK "V1"; Data S2 L75
    lq      <- fixed(log(0.0456)); label("Intercompartmental clearance at the reference body weight (L/h)")     # Data S2 ADAMTS13_Model.m L78 (Table S1 Human PK "Q": 0.456)
    lvp     <- fixed(log(3.71));   label("Peripheral volume of distribution at the reference body weight (L)")  # Table S1 Human PK "V2"; Data S2 L76
    e_wt_cl <- fixed(0.75);        label("Allometric exponent on CL and Q (unitless)")                          # Table S1 Human PK; Data S2 L77-L78
    e_wt_vc <- fixed(1.0);         label("Allometric exponent on Vc and Vp (unitless)")                         # Table S1 Human PK; Data S2 L75-L76
    wt_ref  <- fixed(68.7);        label("Reference body weight for allometric scaling (kg)")                   # Table S1 Human PK
    lfrel   <- fixed(log(1 - 0.0737)); label("Relative ADAMTS13 activity of the dosed product (unitless)")      # Table S1 Human PK "Frel"; rADAMTS13 = (1 - 0.0737), PBT = (1 - 0.390)

    # ------------------------------------------------------------------
    # ADAMTS13 / VWF binding and catalysis (Table S1). VE = elongated
    # (active) VWF monomer unit, VG = globular VWF monomer unit.
    # ------------------------------------------------------------------
    kd_ve_adam     <- fixed(0.0101);  label("Dissociation constant, elongated VWF : (r)ADAMTS13 (nM)")       # Data S2 L97 (Table S1 KD_VE_ADAM: 0.01 nM); calibrated in Step 2
    kd_ve_hb       <- fixed(457.5);   label("Dissociation constant, elongated VWF : hemoglobin (nM)")        # Table S1 KD_VE_Hb
    kd_ve_tsp1     <- fixed(0.10);    label("Dissociation constant, elongated VWF : thrombospondin-1 (nM)")  # Table S1 KD_VE_TSP; Novelli 2013
    kd_ve_platelet <- fixed(4.37);    label("Dissociation constant, elongated VWF : platelet (nM)")          # Table S1 KD_VE_Platelet; calibrated in Step 2
    kd_vg_adam     <- fixed(80);      label("Dissociation constant, globular VWF : (r)ADAMTS13 (nM)")        # Table S1 KD_VG_ADAM; Zanardelli 2009, Feys 2009
    kd_vg_hb       <- fixed(20587.5); label("Dissociation constant, globular VWF : hemoglobin (nM)")         # Table S1 KD_VG_Hb

    kon_ve_adam     <- fixed(527);  label("Association rate constant, elongated VWF : (r)ADAMTS13 (1/nM/h)")      # Table S1 kon_VE_ADAM
    kon_ve_hb       <- fixed(0.36); label("Association rate constant, elongated VWF : hemoglobin (1/nM/h)")       # Table S1 kon_VE_Hb; Schlosshauer 2004
    kon_ve_tsp1     <- fixed(0.36); label("Association rate constant, elongated VWF : thrombospondin-1 (1/nM/h)") # Table S1 kon_VE_TSP
    kon_ve_platelet <- fixed(2.36); label("Association rate constant, elongated VWF : platelet (1/nM/h)")         # Table S1 kon_VE_platelet
    kon_vg_adam     <- fixed(0.36); label("Association rate constant, globular VWF : (r)ADAMTS13 (1/nM/h)")       # Table S1 Kon_VG_ADAM
    kon_vg_hb       <- fixed(0.36); label("Association rate constant, globular VWF : hemoglobin (1/nM/h)")        # Table S1 Kon_VG_Hb (Data S2 L109: 1e-4 * 3600 = 0.36)

    kcat_ve        <- fixed(0.168);  label("Cleavage rate of elongated VWF by (r)ADAMTS13 (1/h)")                    # Table S1 kcat_E
    kcat_vg        <- fixed(0);      label("Cleavage rate of globular VWF by (r)ADAMTS13 (1/h)")                     # Table S1 kcatG; Crawley 2011 (globular VWF is not a substrate)
    kcat_ve_ag     <- fixed(0.0383); label("Cleavage rate of VWF-platelet aggregates by (r)ADAMTS13 (1/nM/h)")       # Table S1 kcatS_ag; the "/hour" in Table S1 omits the /nM (Data S2 L99)
    radam_ag_ec50  <- fixed(0.218);  label("ADAMTS13 activity giving half-maximal aggregate cleavage (IU/mL)")       # Table S1 rADAM_Ag_EC50
    k_ag           <- fixed(3.58e8); label("Platelet aggregate formation rate (1/nM^2/h)")                           # Table S1 k_ag
    n_pl           <- fixed(2);      label("Aggregation-reaction order in the summed VWF:platelet species")          # Table S1 n_pl
    f_cat_platelet <- fixed(1);      label("Ratio of the (r)ADAMTS13 cleavage rate on platelet-bound vs free elongated VWF (unitless)")  # Methods S1 Eq 2d (kcat_VE:PLT = kcat_VE); Data S2 L85 adam_pl_cat_factor

    # ------------------------------------------------------------------
    # Degradation rates (Table S1; each derived from a reported half-life)
    # ------------------------------------------------------------------
    kdeg_adam      <- fixed(log(2) / 60);  label("Degradation rate of endogenous ADAMTS13 (1/h)")          # Table S1 kdeg_ADAM = 0.0116/h, half-life 60 h
    kdeg_vwf       <- fixed(log(2) / 15);  label("Degradation rate of VWF and its complexes (1/h)")        # Table S1 kdeg_VWF = 0.0462/h, half-life 15 h; Denis 2008
    kdeg_vwf_frag  <- fixed(log(2) / 0.5); label("Degradation rate of VWF cleavage fragments (1/h)")       # Data S2 L116, half-life 0.5 h
    kdeg_hb        <- fixed(log(2) / 1);   label("Degradation rate of extracellular hemoglobin (1/h)")     # Table S1 Kdeg_Hb = 0.6931/h, half-life 1 h
    kdeg_tsp1      <- fixed(log(2) / 9);   label("Degradation rate of thrombospondin-1 (1/h)")             # Table S1 kdeg_TSP = 0.0770/h, half-life 9 h; Barclay 2016
    kdeg_platelet  <- fixed(log(2) / 168); label("Degradation rate of platelets (1/h)")                    # Table S1 kdeg_Platelet = 0.0041/h, half-life 7 days; Grossman 1960
    kdeg_ag        <- fixed(log(2) / 24);  label("Degradation rate of platelet aggregates (1/h)")          # Table S1 kdeg_Ag = 0.0289/h, half-life 24 h; Boudjeltia 2002

    # ------------------------------------------------------------------
    # Baselines and steady-state anchors
    # ------------------------------------------------------------------
    ksyn_platelet <- fixed(2.37e-6); label("Platelet synthesis rate (nM/h)")                                           # Table S1 ksyn_Plate_human (constant by design; assumption 4 of the Methods)
    vwf_ss_ref    <- fixed(46.8);    label("Total VWF monomer-unit concentration in a healthy reference subject (nM)") # Data S2 L120 (46.8 nM x 220 kDa = 10.3 ug/mL, the Table S1 healthy VWF_SS)
    vwf_pct       <- fixed(100);     label("Baseline total VWF relative to the healthy reference (%)")                 # Data S2 L62 (individual value); 100% = the healthy reference
    adam_pct      <- fixed(1e-5);    label("Baseline endogenous ADAMTS13 activity in cTTP, relative to 1 IU/mL healthy (unitless)")  # Data S2 L90 (Table S1 ADAMTS13_SS cTTP: 0.005 IU/mL; see vignette Errata)
    hb_ss         <- fixed(89);      label("Baseline extracellular hemoglobin (ug/mL)")                                # Data S2 L86
    tsp1_ss       <- fixed(3200);    label("Baseline thrombospondin-1 in cTTP (ng/mL)")                                # Table S1 TSP1_SS cTTP = 3.2 ug/mL; Liu 2005, Gonzalez 2004, Tan 2008
    adam_act_iu_ug <- fixed(1.5);    label("ADAMTS13 activity per unit mass (IU/ug)")                                  # Table S1 ADAM_Activity_UperUg; Rieger 2006

    # ------------------------------------------------------------------
    # Elongated (active) VWF fraction and the acute-event pulse.
    # A TTP-triggering event is modelled as a pulse-shaped rise in the
    # active VWF fraction (Methods, "The potential physiological factors
    # that trigger a TTP episode were modeled as a pulse-shape rise in
    # the active amount of VWF from baseline").
    # ve_frac_base is the Table S1 population value; spike_amount,
    # spike_day and spike_duration are individual-level event descriptors
    # for which Table S1 reports no population value, so the shipped
    # virtual patient VP #621 values (Data S2 L57-L63, L72) are used.
    # Setting spike_amount == ve_frac_base makes the pulse a no-op.
    # spike_day is on the model's own time axis (t = 0), NOT relative to
    # the first dose, so a simulation that begins with a no-dose lead-in
    # must add the lead-in length to it.
    # ------------------------------------------------------------------
    ve_frac_base   <- fixed(0.0271);      label("Baseline fraction of total VWF in the elongated (active) form (unitless)")  # Table S1 VE_Frac
    spike_amount   <- fixed(0.056882853); label("Elongated VWF fraction during an acute TTP event (unitless)")               # Data S2 L61 (virtual patient #621); not paper-tabulated
    spike_day      <- fixed(56.04229659); label("Start day of the acute TTP event (day)")                                    # Data S2 L63 (virtual patient #621); not paper-tabulated
    spike_duration <- fixed(168);         label("Duration of the acute TTP event (h)")                                       # Data S2 L72
    spike_steep    <- fixed(0.5);         label("Transition half-width of the event pulse (h)")                              # Data S2 L168-L169
    tau_pd         <- fixed(0.1);         label("Transfer time constant from the PK central compartment into the PD ADAMTS13 pool (h)")  # Data S2 L89
  })

  model({
    # ==================================================================
    # 0. Physicochemical constants and unit conversions (Table S1
    #    "Physicochemical parameters"). Molecular weights in Da (g/mol).
    # ==================================================================
    mw_adamts13 <- 190000   # Da; Table S1 MW_ADAM; Levy 2005
    mw_hb       <- 64458    # Da; Table S1 MW_Hb (tetramer); Van Beekvelt 2001
    mw_tsp1     <- 155000   # Da; Table S1 MW_TSP1; Reynolds 2021
    mw_vwf_mono <- 220000   # Da; Table S1 MW_VWF Monomer; Stockschlaeder 2014
    # 1 nM of platelets = 6.02e5 x 1e9/L (Data S2 Run_ADAMTS13_model.m L116:
    # count = state[nM] * 1e-9 mol/L/nM * 6.02e23 /mol * 1e-9 (1e9/L))
    plt_nM_to_1e9L <- 6.02e5

    # ==================================================================
    # 1. rADAMTS13 PK - two-compartment, allometrically scaled
    # ==================================================================
    cl <- exp(lcl) * (WT / wt_ref)^e_wt_cl   # L/h
    vc <- exp(lvc) * (WT / wt_ref)^e_wt_vc   # L
    q  <- exp(lq)  * (WT / wt_ref)^e_wt_cl   # L/h
    vp <- exp(lvp) * (WT / wt_ref)^e_wt_vc   # L

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ==================================================================
    # 2. koff from KD and kon (Methods S1 Eqs 1a-1f)
    # ==================================================================
    koff_ve_adam     <- kd_ve_adam     * kon_ve_adam
    koff_ve_hb       <- kd_ve_hb       * kon_ve_hb
    koff_ve_tsp1     <- kd_ve_tsp1     * kon_ve_tsp1
    koff_ve_platelet <- kd_ve_platelet * kon_ve_platelet
    koff_vg_adam     <- kd_vg_adam     * kon_vg_adam
    koff_vg_hb       <- kd_vg_hb       * kon_vg_hb

    # ==================================================================
    # 3. Steady-state anchors and synthesis rates (Methods S1 Eqs 2a-2c)
    #    ksyn = kdeg * (steady-state concentration in nM)
    # ==================================================================
    hb_input   <- hb_ss   * 1e6 / mw_hb                       # ug/mL -> nM
    tsp1_input <- tsp1_ss * 1000 / mw_tsp1                    # ng/mL -> nM
    vwf_ss     <- vwf_ss_ref * vwf_pct / 100                  # nM
    adam_ss    <- (1 / adam_act_iu_ug) * 1e6 / mw_adamts13 * adam_pct  # nM; 1 IU/mL healthy scaled by adam_pct

    ksyn_hb   <- kdeg_hb   * hb_input       # nM/h
    ksyn_vwf  <- kdeg_vwf  * vwf_ss         # nM/h
    ksyn_adam <- kdeg_adam * adam_ss        # nM/h
    ksyn_tsp1 <- kdeg_tsp1 * tsp1_input     # nM/h

    # ==================================================================
    # 4. Time-varying active-VWF fraction (Data S2 L168-L169).
    #    A symmetric tanh pulse from ve_frac_base up to
    #    spike_amount * (168 / spike_duration), starting at spike_day.
    # ==================================================================
    spike_start <- spike_day * 24
    ve_frac <- ve_frac_base +
      (spike_amount * (168 / spike_duration) - ve_frac_base) *
      0.5 * (tanh((t - spike_start) / spike_steep) -
               tanh((t - spike_start - spike_duration) / spike_steep))

    # ==================================================================
    # 5. Derived quantities (Methods S1 Eqs 3a-3c, 2d)
    # ==================================================================
    ve_plt_sum <- ve_platelet + ve_tsp1_platelet + ve_hb_platelet          # nM
    # Total ADAMTS13 activity (endogenous + recombinant), IU/mL
    adam_iu_ml <- (radam + adam) * mw_adamts13 * 1e-6 * adam_act_iu_ug
    # Saturable enhancement of aggregate lysis by ADAMTS13 (Eq 3c)
    f_a13      <- 1 / (1 + radam_ag_ec50 / adam_iu_ml)
    kcat_ve_platelet <- kcat_ve * f_cat_platelet                            # Eq 2d
    # rADAMTS13 concentration in the central compartment, nM
    cc_central_nM <- (central / vc) / mw_adamts13 * 1e6

    # ==================================================================
    # 6. Reaction rates r1-r61 (Data S2 ADAMTS13_Model.m L182-L276;
    #    equivalent to Methods S1 Eqs 4-25). Numbering follows Data S2.
    # ==================================================================
    r1  <- koff_ve_adam * ve_adam - kon_ve_adam * ve * adam
    r2  <- koff_ve_hb   * ve_hb   - kon_ve_hb   * ve * hb
    r3  <- koff_ve_tsp1 * ve_tsp1 - kon_ve_tsp1 * ve * tsp1
    r4  <- ksyn_vwf * ve_frac
    r5  <- kdeg_vwf * ve
    r6  <- koff_vg_adam * vg_adam - kon_vg_adam * vg * adam
    r7  <- koff_vg_hb   * vg_hb   - kon_vg_hb   * vg * hb
    r8  <- ksyn_vwf * (1 - ve_frac)
    r9  <- kdeg_vwf * vg
    r10 <- kcat_vg * vg_adam
    r11 <- kcat_ve * ve_adam
    r12 <- ksyn_hb
    r13 <- kdeg_hb * hb
    r14 <- kdeg_vwf * (ve_hb + vg_hb)
    r15 <- kdeg_vwf_frag * ve_frag
    r16 <- kdeg_vwf_frag * vg_frag
    r17 <- kdeg_vwf * ve_adam
    r18 <- kdeg_vwf * vg_adam
    r19 <- kdeg_vwf * ve_hb
    r20 <- kdeg_vwf * vg_hb
    r21 <- kdeg_vwf * ve_tsp1
    r22 <- kon_vg_adam * vg * radam - koff_vg_adam * vg_radam
    r23 <- kon_ve_adam * ve * radam - koff_ve_adam * ve_radam
    r24 <- kcat_vg * vg_radam
    r25 <- kcat_ve * ve_radam
    r26 <- kdeg_vwf * vg_radam
    r27 <- kdeg_vwf * ve_radam
    r28 <- ksyn_adam
    r29 <- kdeg_adam * adam
    r30 <- kon_ve_platelet * ve * platelet - koff_ve_platelet * ve_platelet
    r31 <- ksyn_platelet
    r32 <- kdeg_platelet * platelet
    r33 <- kdeg_platelet * ve_platelet
    r34 <- kon_ve_adam * ve_platelet * adam - koff_ve_adam * ve_adam_platelet
    r35 <- kon_ve_platelet * ve_adam * platelet - koff_ve_platelet * ve_adam_platelet
    r36 <- kcat_ve_platelet * ve_adam_platelet
    r37 <- kdeg_platelet * ve_adam_platelet
    r38 <- kon_ve_platelet * ve_hb * platelet - koff_ve_platelet * ve_hb_platelet
    r39 <- kon_ve_hb * ve_platelet * hb - koff_ve_hb * ve_hb_platelet
    r40 <- kdeg_platelet * ve_hb_platelet
    r41 <- kon_ve_platelet * ve_tsp1 * platelet - koff_ve_platelet * ve_tsp1_platelet
    r42 <- kon_ve_tsp1 * ve_platelet * tsp1 - koff_ve_tsp1 * ve_tsp1_platelet
    r43 <- kdeg_platelet * ve_tsp1_platelet
    r44 <- kon_ve_platelet * ve_radam * platelet - koff_ve_platelet * ve_radam_platelet
    r45 <- kon_ve_adam * ve_platelet * radam - koff_ve_adam * ve_radam_platelet
    r46 <- kcat_ve_platelet * ve_radam_platelet
    r47 <- kdeg_platelet * ve_radam_platelet
    r48 <- ksyn_tsp1
    r49 <- kdeg_tsp1 * tsp1
    r50 <- k_ag * ve_platelet      * ve_plt_sum^n_pl
    r51 <- k_ag * ve_tsp1_platelet * ve_plt_sum^n_pl
    r52 <- k_ag * ve_hb_platelet   * ve_plt_sum^n_pl
    r53 <- kdeg_ag * ve_platelet_coag
    r54 <- kdeg_ag * ve_tsp1_platelet_coag
    r55 <- kdeg_ag * ve_hb_platelet_coag
    # Data S2 sets the de-aggregation rates r56-r58 to zero; aggregates
    # resolve only by degradation (r53-r55) or ADAMTS13 lysis (r59-r61).
    r59 <- kcat_ve_ag * (adam + radam) * ve_platelet_coag      * f_a13
    r60 <- kcat_ve_ag * (adam + radam) * ve_tsp1_platelet_coag * f_a13
    r61 <- kcat_ve_ag * (adam + radam) * ve_hb_platelet_coag   * f_a13

    # ==================================================================
    # 7. ODE system
    # ==================================================================
    d/dt(vg)   <- r6 + r7 + r8 - r9 - r22                                  # Methods S1 Eq 7
    d/dt(ve)   <- r1 + r2 + r3 + r4 - r5 - r23 - r30                       # Methods S1 Eq 6
    d/dt(adam) <- r6 + r1 + r10 + r11 + r28 - r29 + r17 + r18 - r34 + r36  # Methods S1 Eq 4
    d/dt(hb)   <- r2 + r7 + r12 - r13 + r14 - r39 + r61                    # Methods S1 Eq 9
    d/dt(tsp1) <- r3 - r42 + r48 - r49 + r60                               # Methods S1 Eq 10

    d/dt(vg_adam)  <- -r6 - r10 - r18                                      # Methods S1 Eq 15
    d/dt(ve_adam)  <- -r1 - r11 - r17 - r35                                # Methods S1 Eq 11
    d/dt(vg_hb)    <- -r7 - r20                                            # Methods S1 Eq 17
    d/dt(ve_hb)    <- -r2 - r19 - r38                                      # Methods S1 Eq 13
    d/dt(ve_tsp1)  <- -r3 - r21 - r41                                      # Methods S1 Eq 14

    d/dt(vg_frag)  <- r10 - r16 + r24                                      # Data S2 L289
    d/dt(ve_frag)  <- r11 - r15 + r25 + r36 + r46 + r59 + r60 + r61        # Data S2 L290

    d/dt(vg_radam) <- r22 - r24 - r26                                      # Methods S1 Eq 16
    d/dt(ve_radam) <- r23 - r25 - r27 - r44                                # Methods S1 Eq 12

    d/dt(platelet)          <- r36 + r46 - r30 - r35 - r38 - r41 - r44 + r31 - r32 + r59 + r60 + r61  # Methods S1 Eq 8
    d/dt(ve_platelet)       <- r30 - r33 - r34 - r39 - r42 - r45 - r50     # Methods S1 Eq 18
    d/dt(ve_adam_platelet)  <- r34 + r35 - r36 - r37                       # Methods S1 Eq 19
    d/dt(ve_tsp1_platelet)  <- r41 + r42 - r43 - r51                       # Methods S1 Eq 22
    d/dt(ve_hb_platelet)    <- r38 + r39 - r40 - r52                       # Methods S1 Eq 21
    d/dt(ve_radam_platelet) <- r44 + r45 - r46 - r47                       # Methods S1 Eq 20

    d/dt(ve_platelet_coag)      <- r50 - r53 - r59                         # Methods S1 Eq 23
    d/dt(ve_tsp1_platelet_coag) <- r51 - r54 - r60                         # Methods S1 Eq 24
    d/dt(ve_hb_platelet_coag)   <- r52 - r55 - r61                         # Methods S1 Eq 25

    # rADAMTS13 in the PD compartment tracks the PK central concentration
    # with time constant tau_pd. Data S2 L316 implements the PK as a
    # forcing function on the PD system ("individual predictions from a
    # two-compartment model represented the PK of rADAMTS13 for input
    # into the model", Figure 1 footnote a); the binding reactions
    # r22/r23/r45 therefore do not deplete this pool. Methods S1 Eq 5
    # prints only the binding terms and omits the PK input; see the
    # vignette Errata.
    d/dt(radam) <- (cc_central_nM - radam) / tau_pd

    # rADAMTS13 PK, amounts in mg (Data S2 L318-L319 in concentration form)
    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(central) <- exp(lfrel)

    # ==================================================================
    # 8. Outputs
    # ==================================================================
    Cc  <- central / vc                       # rADAMTS13 in plasma, ug/mL
    # plt is THE published platelet count: Data S2 Run_ADAMTS13_model.m
    # L116 converts only the free `Platelet` state to 1e9/L, so platelets
    # sequestered in VWF complexes and aggregates are not counted.
    plt <- platelet * plt_nM_to_1e9L          # free platelet count, 1e9/L
    # Diagnostic only: platelets in every species, for the mass-balance
    # check. One platelet per aggregate: Methods S1 Eq 23 converts one
    # VE:PLT species into one aggregate species, and `ve_plt_sum^n_pl` is
    # a nonlinear rate dependency (Table S1 n_pl, "coagulation reaction
    # dependency on sum of VE platelet bound compounds"), not a
    # stoichiometric coefficient. Data S2 L300 counts aggregates the same
    # way in its own mass-conservation accumulator.
    plt_total <- (platelet + ve_platelet + ve_adam_platelet + ve_tsp1_platelet +
                    ve_hb_platelet + ve_radam_platelet +
                    ve_platelet_coag + ve_tsp1_platelet_coag +
                    ve_hb_platelet_coag) * plt_nM_to_1e9L
    adamts13Activity <- adam_iu_ml            # total ADAMTS13 activity, IU/mL
    vwfTotal   <- vg + ve + vg_adam + ve_adam + vg_hb + ve_hb + ve_tsp1 +
      vg_radam + ve_radam + ve_platelet + ve_adam_platelet +
      ve_tsp1_platelet + ve_hb_platelet + ve_radam_platelet +
      ve_platelet_coag + ve_tsp1_platelet_coag + ve_hb_platelet_coag   # nM
    vwfActive  <- ve + ve_adam + ve_hb + ve_tsp1 + ve_radam +
      ve_platelet + ve_adam_platelet + ve_tsp1_platelet +
      ve_hb_platelet + ve_radam_platelet +
      ve_platelet_coag + ve_tsp1_platelet_coag + ve_hb_platelet_coag   # nM
    veFrac <- ve_frac                         # instantaneous elongated-VWF fraction
    # Total VWF on the mass scale Table S1 reports it (VWF_SS, ug/mL).
    vwfTotal_ugmL <- vwfTotal * mw_vwf_mono / 1e6
  })
}
