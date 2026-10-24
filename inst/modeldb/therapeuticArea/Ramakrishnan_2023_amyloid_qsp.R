Ramakrishnan_2023_amyloid_qsp <- function() {
  description <- paste(
    "QSP. Amyloid-beta (A-beta) pathway natural-history model in Alzheimer's",
    "disease (Genentech / Rosa & Co. SimBiology model): APP synthesis and",
    "cleavage by BACE1 to C99 and sAPP-beta, gamma-secretase cleavage of C99",
    "to A-beta40 and A-beta42 monomers, a sequestered A-beta42 pool,",
    "pathway-resolved monomer clearance (astrocyte receptor uptake,",
    "enzymatic degradation, microglial uptake, bulk flow to CSF, active",
    "transport to plasma), aggregation of monomers to oligomers (10-mers),",
    "fibrils and plaque in a single well-mixed brain compartment, and",
    "transport of soluble species to CSF and plasma with peripheral",
    "A-beta production. Fourteen ODE states. Calibrated to an APOE-epsilon4",
    "carrier virtual patient; APOE4_CARRIER = 0 applies the paper's Table 3",
    "mechanistic differences for the APOE-epsilon4 non-carrier virtual",
    "patient. Soluble species start at their steady state and insoluble",
    "A-beta grows linearly. The anti-A-beta antibody layer of the paper",
    "(solanezumab, crenezumab, aducanumab, gantenerumab) is NOT included:",
    "its PK, binding-rate, compartment-volume and microglial-activation",
    "values are not published (Table S3 lists the names only).",
    "Deterministic: no between-subject variability or residual error.",
    sep = " "
  )
  reference <- paste(
    "Ramakrishnan V, Friedrich C, Witt C, Sheehan R, Pryor M, Atwal JK,",
    "Wildsmith K, Kudrycki K, Lee SH, Mazer N, Hofmann C, Fuji RN, Jin JY,",
    "Ramanujan S, Dolton M, Quartino A. Quantitative systems pharmacology",
    "model of the amyloid pathway in Alzheimer's disease: Insights into the",
    "therapeutic mechanisms of clinical candidates. CPT Pharmacometrics Syst",
    "Pharmacol. 2023;12(1):62-73. doi:10.1002/psp4.12876.",
    "ODEs, initial amounts (Table S1), parameter values (Table S2) and rule",
    "definitions (Table S4) from the Supplementary Material; APOE-epsilon4",
    "non-carrier values from main-text Table 3.",
    sep = " "
  )
  vignette <- "Ramakrishnan_2023_amyloid_qsp"

  # None of these states maps onto a canonical nlmixr2lib compartment role.
  # The `_brain` / `_csf` / `_plasma` suffix follows the paper's three
  # physiologic compartments (Table S1 'Compartment' column).
  paper_specific_compartments <- c(
    "app",
    "sappb",
    "c99",
    "ab40_brain",
    "ab42_brain",
    "olig_brain",
    "fibril_brain",
    "plaque_brain",
    "ab40_csf",
    "ab42_csf",
    "olig_csf",
    "ab40_plasma",
    "ab42_plasma",
    "olig_plasma"
  )

  units <- list(
    time = "h",
    dosing = "none",
    concentration = paste(
      "app, sappb and c99 are neuronal concentrations in umol/L (Table S1).",
      "All A-beta states are amounts in pmol (Table S1 reports moles; the",
      "model rescales by 1e12 so the states are numerically well scaled).",
      "Oligomer, fibril and plaque amounts are in oligomer units (1 oligomer",
      "unit = 10 monomers). ratio_ab4240_<cmt> outputs are dimensionless",
      "A-beta42:A-beta40 monomer ratios; insol_brain is fibril + plaque in",
      "pmol oligomer units; pct_insol is its percent change from the",
      "calibrated baseline.",
      sep = " "
    )
  )

  compartmentData <- list(
    app = list(
      analyte = "amyloid precursor protein (APP)",
      units = "umol/L (neuronal concentration)",
      specimen = "tissue",
      verified = TRUE
    ),
    sappb = list(
      analyte = "soluble APP-beta (sAPP-beta)",
      units = "umol/L (neuronal concentration)",
      specimen = "tissue",
      verified = TRUE
    ),
    c99 = list(
      analyte = "C99 fragment of APP",
      units = "umol/L (neuronal concentration)",
      specimen = "tissue",
      verified = TRUE
    ),
    ab40_brain = list(analyte = "A-beta40 monomer", units = "pmol", specimen = "brain ISF", verified = TRUE),
    ab42_brain = list(analyte = "A-beta42 monomer", units = "pmol", specimen = "brain ISF", verified = TRUE),
    olig_brain = list(
      analyte = "A-beta oligomer (A-beta40 + A-beta42 10-mers)",
      units = "pmol oligomer units",
      specimen = "brain ISF",
      verified = TRUE
    ),
    fibril_brain = list(
      analyte = "A-beta fibril (insoluble)",
      units = "pmol oligomer units",
      specimen = "tissue",
      verified = TRUE
    ),
    plaque_brain = list(
      analyte = "A-beta plaque (insoluble)",
      units = "pmol oligomer units",
      specimen = "tissue",
      verified = TRUE
    ),
    ab40_csf = list(analyte = "A-beta40 monomer", units = "pmol", specimen = "CSF", verified = TRUE),
    ab42_csf = list(analyte = "A-beta42 monomer", units = "pmol", specimen = "CSF", verified = TRUE),
    olig_csf = list(analyte = "A-beta oligomer", units = "pmol oligomer units", specimen = "CSF", verified = TRUE),
    ab40_plasma = list(analyte = "A-beta40 monomer", units = "pmol", specimen = "plasma", verified = TRUE),
    ab42_plasma = list(analyte = "A-beta42 monomer", units = "pmol", specimen = "plasma", verified = TRUE),
    olig_plasma = list(
      analyte = "A-beta oligomer",
      units = "pmol oligomer units",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    APOE4_CARRIER = list(
      description = paste(
        "APOE-epsilon4 carrier indicator selecting the virtual patient:",
        "1 = APOE-epsilon4 carrier (the calibrated virtual patient, Table S2",
        "values), 0 = APOE-epsilon4 non-carrier (main-text Table 3 values for",
        "BACE1 Km and Vmax and for the astrocyte-receptor, brain-to-CSF and",
        "brain-to-plasma monomer clearance fractions).",
        sep = " "
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "1 = APOE-epsilon4 carrier (the calibrated virtual patient)",
      notes = paste(
        "The paper contrasts two virtual patients only; it does not separate",
        "heterozygous from homozygous carriers, so the binary carrier",
        "canonical is used. Time-invariant. The default initial amounts are",
        "the carrier's Table S1 steady state; the paper reports no separate",
        "non-carrier initial state. Its Figure 4 shows both virtual patients",
        "progressing linearly from week 0, which happens only once the",
        "non-carrier's soluble and fibril states have relaxed to their own",
        "steady state; the vignette shows how to pre-equilibrate them.",
        sep = " "
      ),
      source_name = "APOE e4 + VP / APOE e4 - VP (Table 3 column headers)"
    )
  )

  population <- list(
    species = "human (in silico virtual patients)",
    n_subjects = NA_integer_,
    n_studies = NA_integer_,
    disease_state = paste(
      "Alzheimer's disease during the ~20-year phase of linear A-beta",
      "aggregation (the upward slope of the amyloid S curve). Two virtual",
      "patients: an APOE-epsilon4 carrier (calibrated) and an APOE-epsilon4",
      "non-carrier (test case).",
      sep = " "
    ),
    dose_range = paste(
      "No dosing in this file. The paper calibrated its antibody layer to",
      "solanezumab, crenezumab, aducanumab and gantenerumab trials; those",
      "parameter values are not published.",
      sep = " "
    ),
    notes = paste(
      "Deterministic typical-value QSP model calibrated to literature and",
      "internal Genentech data rather than fitted to a population. Table 1",
      "lists the CSF A-beta42 and amyloid-PET SUVR values (normal,",
      "prodromal, AD) used to guide virtual patient development. Initial",
      "amounts of brain A-beta species derive from McDonald 2015 and Roberts",
      "2017; CSF A-beta42 from Genentech study ABE4869g; plasma A-beta from",
      "Roher 2009, Ovod 2017 and internal data (Table S1).",
      sep = " "
    )
  )

  ini({
    # ================================================================
    # All values are from the Supplementary Material of Ramakrishnan 2023
    # (Table S1 initial amounts, Table S2 parameters) unless the comment
    # says main-text Table 3. They are mechanistic constants of a
    # deterministic SimBiology model, reported as point values with no
    # uncertainty, so they are kept on the LINEAR scale under snake-case
    # names derived from the paper's own identifiers, and every one is
    # wrapped in fixed(). Table S2 marks some values 'Fitted' (during
    # calibration), but nothing is estimable from a user's data here.
    # Rates are per hour; amounts given in moles in the paper are
    # converted to pmol (x 1e12).
    # ================================================================

    # ---- APP processing (neuronal concentrations, umol/L) ----
    k_app_prod <- fixed(0.399); label("APP production rate APP_production_rate_k (umol/L/h)") # Table S2 'APP_production_rate_k' 0.399, Fitted
    k_app_c99 <- fixed(1); label("APP cleavage rate to C99 APP_to_C99_rate_k (1/h)") # Table S2 'APP_to_C99_rate_k' 1, Normalized
    bace1_emax <- fixed(0.455); label("BACE1 maximum activation of APP cleavage, APOE4 carrier BACE1_APP_toC99_emax (unitless)") # Table S2 'BACE1_APP_toC99_emax' 0.455; Table 3 'BACE1 Vmax' APOE4+ 0.455
    bace1_emax_nc <- fixed(0.253); label("BACE1 maximum activation of APP cleavage, APOE4 non-carrier (unitless)") # Table 3 'BACE1 Vmax' APOE4- VP value 0.253
    bace1_ec50 <- fixed(7.76); label("APP level for half-maximal BACE1 cleavage, APOE4 carrier BACE1_APP_toC99_ec50 (umol/L)") # Table S2 'BACE1_APP_toC99_ec50' 7.76; Table 3 'BACE1 Km' APOE4+ 7.760
    bace1_ec50_nc <- fixed(8.438); label("APP level for half-maximal BACE1 cleavage, APOE4 non-carrier (umol/L)") # Table 3 'BACE1 Km' APOE4- VP value 8.438
    k_app_sappa <- fixed(0.0142); label("APP cleavage rate to sAPP-alpha APP_to_sAPPalpha_rate_k (1/h)") # Table S2 'APP_to_sAPPalpha_rate_k' 0.0142, Fitted
    k_sappb_cl <- fixed(0.7); label("sAPP-beta clearance rate sAPPbeta_clearance_rate_k (1/h)") # Table S2 'sAPPbeta_clearance_rate_k' 0.7, Assumption
    k_c99_abeta <- fixed(8.3e-08); label("C99 conversion to A-beta and secretion into ISF C99_to_Abeta_rate_k ((L/h)*(mol/umol))") # Table S2 'C99_to_Abeta_rate_k' 8.3E-08, Fitted
    gsec_ab42_emax <- fixed(0.00095); label("Gamma-secretase maximum activation of C99 cleavage to A-beta42 (unitless)") # Table S2 'gammaSec_C99_toAbeta42_emax' 0.00095
    gsec_ab42_ec50 <- fixed(5.16); label("C99 level for half-maximal cleavage to A-beta42 (umol/L)") # Table S2 'gammaSec_C99_toAbeta42_ec50' 5.16
    gsec_ab40_emax <- fixed(0.00166); label("Gamma-secretase maximum activation of C99 cleavage to A-beta40 (unitless)") # Table S2 'gammaSec_C99_toAbeta40_emax' 0.00166
    gsec_ab40_ec50 <- fixed(0.97); label("C99 level for half-maximal cleavage to A-beta40 (umol/L)") # Table S2 'gammaSec_C99_toAbeta40_ec50' 0.97
    k_c99_cl <- fixed(0.0508); label("C99 clearance rate C99_clearance_rate_k (1/h)") # Table S2 'C99_clearance_rate_k' 0.0508, Fitted
    bace1 <- fixed(1); label("BACE1 activity, normalized (unitless)") # Table S1 'BACE1' initial amount 1, normalized
    gsec <- fixed(1); label("Gamma-secretase activity, normalized (unitless)") # Table S1 'gamma_Secretase' initial amount 1, normalized

    # ---- A-beta40 monomer, brain ISF ----
    k_ab40_cl <- fixed(0.0834); label("A-beta40 total clearance rate in brain ISF Abeta40M_clearance_rate_k (1/h)") # Table S2 'Abeta40M_clearance_rate_k' 0.0834 (Bateman 2006)
    f_ab40_astro <- fixed(0.1); label("Fraction of A-beta40 clearance by astrocyte receptor uptake, APOE4 carrier (unitless)") # Table S2 'ReceptorAstro_Abeta40M_clearance_fraction' 0.1; Table 3 APOE4+ 0.100
    f_ab40_astro_nc <- fixed(0.2); label("Fraction of A-beta40 clearance by astrocyte receptor uptake, APOE4 non-carrier (unitless)") # Table 3 'Abeta40 astrocyte receptor clearance fraction' APOE4- 0.200
    f_ab40_degenz <- fixed(0.199); label("Fraction of A-beta40 clearance by enzymatic degradation (unitless)") # Table S2 'AbetaDegEnz_Abeta40M_clearance_fraction' 0.199
    f_ab40_mglia <- fixed(0.05); label("Fraction of A-beta40 clearance by microglial phagocytosis (unitless)") # Table S2 'MicrogliaAct_Abeta40M_clearance_fraction' 0.05
    f_ab40_olig <- fixed(0.001); label("Fraction of A-beta40 monomer clearance incorporated into oligomer (unitless)") # Table S2 'Abeta40M_oligagg_clearance_fraction' 0.001, Fitted
    f_ab40_csf <- fixed(0.301); label("Fraction of A-beta40 clearance by bulk flow to CSF, APOE4 carrier (unitless)") # Table S2 'Abeta40M_BraintoCSF_clearance_fraction' 0.301; Table 3 APOE4+ 0.301
    f_ab40_csf_nc <- fixed(0.602); label("Fraction of A-beta40 clearance by bulk flow to CSF, APOE4 non-carrier (unitless)") # Table 3 'Abeta40 brain to CSF clearance fraction' APOE4- 0.602
    f_ab40_bbb <- fixed(0.349); label("Fraction of A-beta40 clearance by active BBB transport to plasma, APOE4 carrier (unitless)") # Table S2 'Abeta40M_BraintoPlasma_Active_clearance_fraction' 0.349; Table 3 APOE4+ 0.349
    f_ab40_bbb_nc <- fixed(0.698); label("Fraction of A-beta40 clearance by active BBB transport to plasma, APOE4 non-carrier (unitless)") # Table 3 'Abeta40 brain to plasma active clearance fraction' APOE4- 0.698

    # ---- A-beta42 monomer, brain ISF ----
    k_ab42_cl <- fixed(0.0356); label("A-beta42 total clearance rate in brain ISF Abeta42M_clearance_rate_k (1/h)") # Table S2 'Abeta42M_clearance_rate_k' 0.0356, fitted to the A-beta42:A-beta40 ratio
    f_ab42_astro <- fixed(0.117); label("Fraction of A-beta42 clearance by astrocyte receptor uptake, APOE4 carrier (unitless)") # Table S2 'ReceptorAstro_Abeta42M_clearance_fraction' 0.117; Table 3 APOE4+ 0.117
    f_ab42_astro_nc <- fixed(0.235); label("Fraction of A-beta42 clearance by astrocyte receptor uptake, APOE4 non-carrier (unitless)") # Table 3 'Abeta42 astrocyte receptor clearance fraction' APOE4- 0.235
    f_ab42_degenz <- fixed(0.159); label("Fraction of A-beta42 clearance by enzymatic degradation (unitless)") # Table S2 'AbetaDegEnz_Abeta42M_clearance_fraction' 0.159
    f_ab42_mglia <- fixed(0.117); label("Fraction of A-beta42 clearance by microglial phagocytosis (unitless)") # Table S2 'MicrogliaAct_Abeta42M_clearance_fraction' 0.117
    f_ab42_olig <- fixed(0.007); label("Fraction of A-beta42 monomer clearance incorporated into oligomer (unitless)") # Table S2 'Abeta42M_oligagg_clearance_fraction' 0.007, Fitted
    f_ab42_csf <- fixed(0.353); label("Fraction of A-beta42 clearance by bulk flow to CSF, APOE4 carrier (unitless)") # Table S2 'Abeta42M_BraintoCSF_clearance_fraction' 0.353; Table 3 APOE4+ 0.353
    f_ab42_csf_nc <- fixed(0.707); label("Fraction of A-beta42 clearance by bulk flow to CSF, APOE4 non-carrier (unitless)") # Table 3 'Abeta42 brain to CSF clearance fraction' APOE4- 0.707
    f_ab42_bbb <- fixed(0.246); label("Fraction of A-beta42 clearance by active BBB transport to plasma, APOE4 carrier (unitless)") # Table S2 'Abeta42M_BraintoPlasma_Active_clearance_fraction' 0.246; Table 3 APOE4+ 0.246
    f_ab42_bbb_nc <- fixed(0.492); label("Fraction of A-beta42 clearance by active BBB transport to plasma, APOE4 non-carrier (unitless)") # Table 3 'Abeta42 brain to plasma active clearance fraction' APOE4- 0.492
    k_ab42_seq <- fixed(3.57e-11 * 1e12); label("A-beta42 secretion from sequestered brain pools Abeta42M_SequesteredPool_secretion_rate (pmol/h)") # Table S2 'Abeta42M_SequesteredPool_secretion_rate' 3.57E-11 moles/hour
    f_ab42_dilution <- fixed(0.125); label("Brain-to-CSF and brain-to-plasma dilution factor for A-beta42 (unitless)") # constant 0.125 printed in the d(Abeta42_M_CSF)/dt and d(Abeta42_M_Pl)/dt ODEs (Supplemental methods)

    # ---- A-beta oligomer, brain ISF ----
    k_olig_cl <- fixed(0.000557); label("A-beta oligomer total clearance rate in brain ISF AbetaAllO_clearance_rate_k (1/h)") # Table S2 'AbetaAllO_clearance_rate_k' 0.000557
    f_olig_astro <- fixed(0.00553); label("Fraction of oligomer clearance by astrocyte receptor uptake (unitless)") # Table S2 'ReceptorAstro_AbetaAllO_clearance_fraction' 0.00553
    f_olig_degenz <- fixed(0.00734); label("Fraction of oligomer clearance by enzymatic degradation (unitless)") # Table S2 'AbetaDegEnz_AbetaAllO_clearance_fraction' 0.00734
    f_olig_mglia <- fixed(0.00922); label("Fraction of oligomer clearance by microglial phagocytosis (unitless)") # Table S2 'MicrogliaAct_AbetaAllO_clearance_fraction' 0.00922
    f_olig_fib <- fixed(0.962); label("Fraction of oligomer clearance by aggregation into fibrils (unitless)") # Table S2 'AbetaAllO_fibrilgagg_clearance_fraction' 0.962, Assumption
    f_olig_csf <- fixed(0.0111); label("Fraction of oligomer clearance by transport to CSF (unitless)") # Table S2 'AbetaAllO_BraintoCSF_clearance_fraction' 0.0111
    f_olig_bbb <- fixed(0.00515); label("Fraction of oligomer clearance by transport to plasma (unitless)") # Table S2 'AbetaAllO_BraintoPlasma_Active_clearance_fraction' 0.00515

    # ---- A-beta fibril and plaque, brain ----
    k_fib_degrade <- fixed(0); label("Fibril disaggregation rate to oligomer AbetaAllF_degrade_rate_k (1/h)") # Table S2 'AbetaAllF_degrade_rate_k' 0, Assumption
    k_fib_cl <- fixed(0.00004); label("Fibril total clearance rate AbetaAllF_clearance_rate_k (1/h)") # Table S2 'AbetaAllF_clearance_rate_k' 0.00004, Assumption
    f_fib_mglia <- fixed(0.3); label("Fraction of fibril clearance by microglia (unitless)") # Table S2 'MicrogliaAct_AbetaAllF_clearance_fraction' 0.3, Assumption
    f_fib_plaque <- fixed(0.7); label("Fraction of fibril clearance by aggregation into plaque (unitless)") # Table S2 'AbetaAllF_plaqueagg_clearance_fraction' 0.7, Assumption
    k_plaque_cl <- fixed(0); label("Endogenous plaque clearance rate AbetaAllP_clearance_rate_k (1/h)") # Table S2 'AbetaAllP_clearance_rate_k' 0, Assumption

    # ---- CSF clearance ----
    k_ab40_csf_cl <- fixed(0.125); label("A-beta40 clearance rate in CSF Abeta40M_CSF_clearance_rate_k (1/h)") # Table S2 'Abeta40M_CSF_clearance_rate_k' 0.125
    k_ab42_csf_cl <- fixed(0.125); label("A-beta42 clearance rate in CSF Abeta42M_CSF_clearance_rate_k (1/h)") # Table S2 'Abeta42M_CSF_clearance_rate_k' 0.125
    k_olig_csf_cl <- fixed(0.125); label("Oligomer clearance rate in CSF AbetaAllO_CSF_clearance_rate_k (1/h)") # Table S2 'AbetaAllO_CSF_clearance_rate_k' 0.125

    # ---- Plasma production and clearance ----
    k_ab40_periph <- fixed(5.72e-11 * 1e12); label("Peripheral A-beta40 production rate Abeta40M_periph_production_rate_k (pmol/h)") # Table S2 'Abeta40M_periph_production_rate_k' 5.72E-11 moles/hour, Fitted
    k_ab42_periph <- fixed(3.75e-12 * 1e12); label("Peripheral A-beta42 production rate Abeta42M_periph_production_rate_k (pmol/h)") # Table S2 'Abeta42M_periph_production_rate_k' 3.75E-12 moles/hour, Fitted
    k_ab40_pl_cl <- fixed(0.231); label("A-beta40 clearance rate in plasma Abeta40M_Pl_clearance_rate_k (1/h)") # Table S2 'Abeta40M_Pl_clearance_rate_k' 0.231
    k_ab42_pl_cl <- fixed(0.231); label("A-beta42 clearance rate in plasma Abeta42M_Pl_clearance_rate_k (1/h)") # Table S2 'Abeta42M_Pl_clearance_rate_k' 0.231
    k_olig_pl_cl <- fixed(0.231); label("Oligomer clearance rate in plasma AbetaAllO_Pl_clearance_rate_k (1/h)") # Table S2 'AbetaAllO_Pl_clearance_rate_k' 0.231

    # ---- Initial amounts (Table S1) ----
    bl_app <- fixed(10); label("Initial APP concentration in neurons (umol/L)") # Table S1 'APP' 10 uM
    bl_sappb <- fixed(0.37); label("Initial sAPP-beta concentration in neurons (umol/L)") # Table S1 'sAPPbeta' 0.37 uM
    bl_c99 <- fixed(5.05); label("Initial C99 concentration in neurons (umol/L)") # Table S1 'C99' 5.05 uM
    bl_ab40_brain <- fixed(1.37e-09 * 1e12); label("Initial A-beta40 monomer amount in brain (pmol)") # Table S1 'Abeta40_M' 1.37E-09 moles
    bl_ab42_brain <- fixed(1.98e-09 * 1e12); label("Initial A-beta42 monomer amount in brain (pmol)") # Table S1 'Abeta42_M' 1.98E-09 moles
    bl_olig_brain <- fixed(1.09e-09 * 1e12); label("Initial A-beta oligomer amount in brain (pmol oligomer units)") # Table S1 'AbetaAll_O' 1.09E-09 moles
    bl_fibril_brain <- fixed(1.46e-08 * 1e12); label("Initial A-beta fibril amount in brain (pmol oligomer units)") # Table S1 'AbetaAll_F' 1.46E-08 moles
    bl_plaque_brain <- fixed(5.89e-08 * 1e12); label("Initial A-beta plaque amount in brain (pmol oligomer units)") # Table S1 'AbetaAll_P' 5.89E-08 moles
    bl_ab40_csf <- fixed(1.38e-10 * 1e12); label("Initial A-beta40 monomer amount in CSF (pmol)") # Table S1 'Abeta40_M_CSF' 1.38E-10 moles
    bl_ab42_csf <- fixed(1.24e-11 * 1e12); label("Initial A-beta42 monomer amount in CSF (pmol)") # Table S1 'Abeta42_M_CSF' 1.24E-11 moles
    bl_olig_csf <- fixed(2.69e-14 * 1e12); label("Initial A-beta oligomer amount in CSF (pmol oligomer units)") # Table S1 'AbetaAll_O_CSF' 2.69E-14 moles
    bl_ab40_plasma <- fixed(4.95e-10 * 1e12); label("Initial A-beta40 monomer amount in plasma (pmol)") # Table S1 'Abeta40_M_Pl' 4.95E-10 moles
    bl_ab42_plasma <- fixed(3.23e-11 * 1e12); label("Initial A-beta42 monomer amount in plasma (pmol)") # Table S1 'Abeta42_M_Pl' 3.23E-11 moles
    bl_olig_plasma <- fixed(2.81e-14 * 1e12); label("Initial A-beta oligomer amount in plasma (pmol oligomer units)") # Table S1 'AbetaAll_O_Pl' 2.81E-14 moles
  })
  model({
    # ---- APOE-epsilon4 virtual patient (main-text Table 3) ----
    bace1_emax_i <- bace1_emax * APOE4_CARRIER + bace1_emax_nc * (1 - APOE4_CARRIER)
    bace1_ec50_i <- bace1_ec50 * APOE4_CARRIER + bace1_ec50_nc * (1 - APOE4_CARRIER)
    f_ab40_astro_i <- f_ab40_astro * APOE4_CARRIER + f_ab40_astro_nc * (1 - APOE4_CARRIER)
    f_ab40_csf_i <- f_ab40_csf * APOE4_CARRIER + f_ab40_csf_nc * (1 - APOE4_CARRIER)
    f_ab40_bbb_i <- f_ab40_bbb * APOE4_CARRIER + f_ab40_bbb_nc * (1 - APOE4_CARRIER)
    f_ab42_astro_i <- f_ab42_astro * APOE4_CARRIER + f_ab42_astro_nc * (1 - APOE4_CARRIER)
    f_ab42_csf_i <- f_ab42_csf * APOE4_CARRIER + f_ab42_csf_nc * (1 - APOE4_CARRIER)
    f_ab42_bbb_i <- f_ab42_bbb * APOE4_CARRIER + f_ab42_bbb_nc * (1 - APOE4_CARRIER)

    # ---- Rules (Table S4), untreated system ----
    # Microglia_Activation is 0 without antibody: every term it multiplies
    # in the printed ODEs and rules drops out, and the soluble species are
    # at steady state at the Table S1 amounts only with it at 0.
    k_ab40_olig <- f_ab40_olig * k_ab40_cl # Abeta40M_oligagg_rate_k_Rule
    k_ab42_olig <- f_ab42_olig * k_ab42_cl # Abeta42M_oligagg_rate_k_Rule
    k_olig_clr <- k_olig_cl * (f_olig_astro + f_olig_degenz + f_olig_mglia) # ABetaAllO_clearance_rate_k_Rule
    k_olig_fib <- f_olig_fib * k_olig_cl # ABetaAllO_fibrilgagg_rate_k_Rule
    k_olig_bbb <- f_olig_bbb * k_olig_cl # ABetaAllO_BraintoPlasma_rate_k_Rule
    k_olig_csf <- f_olig_csf * k_olig_cl # ABetaAllO_BraintoCSF_rate_k_Rule
    k_fib_clr <- k_fib_cl * f_fib_mglia # ABetaAllF_clearance_rate_k_Rule
    k_fib_plaque <- f_fib_plaque * k_fib_cl # ABetaAllF_plaqueagg_rate_k_Rule
    k_plaque_clr <- k_plaque_cl # ABetaAllP_clearance_rate_k_Rule

    # ---- Fluxes ----
    flux_app_c99 <- k_app_c99 * app * (bace1_emax_i * bace1 / (bace1_ec50_i + app))
    # C99 -> A-beta flux in mol/h (k_c99_abeta carries mol/umol). As printed,
    # the same mol/h term is subtracted from the umol/L C99 state; at the
    # baseline it is ~1.5e-10 against a C99 turnover of ~0.26 umol/L/h.
    flux_c99_ab40 <- k_c99_abeta * c99 * gsec_ab40_emax * gsec / (gsec_ab40_ec50 + c99)
    flux_c99_ab42 <- k_c99_abeta * c99 * gsec_ab42_emax * gsec / (gsec_ab42_ec50 + c99)

    # ---- APP processing ----
    d/dt(app) <- k_app_prod - flux_app_c99 - k_app_sappa * app
    d/dt(sappb) <- flux_app_c99 - k_sappb_cl * sappb
    d/dt(c99) <- -flux_c99_ab42 - flux_c99_ab40 + flux_app_c99 - k_c99_cl * c99

    # ---- Brain monomers, oligomers, fibrils, plaque ----
    # The factor 10 converts monomers lost to oligomer units gained.
    d/dt(ab40_brain) <- 1e12 * flux_c99_ab40 -
      f_ab40_csf_i * k_ab40_cl * ab40_brain -
      f_ab40_bbb_i * k_ab40_cl * ab40_brain -
      10 * k_ab40_olig * ab40_brain -
      k_ab40_cl * ab40_brain * (f_ab40_astro_i + f_ab40_degenz + f_ab40_mglia)
    d/dt(ab42_brain) <- 1e12 * flux_c99_ab42 -
      f_ab42_csf_i * k_ab42_cl * ab42_brain -
      f_ab42_bbb_i * k_ab42_cl * ab42_brain -
      10 * k_ab42_olig * ab42_brain -
      k_ab42_cl * ab42_brain * (f_ab42_astro_i + f_ab42_degenz + f_ab42_mglia) +
      k_ab42_seq
    d/dt(olig_brain) <- -k_olig_clr * olig_brain -
      (k_olig_fib * olig_brain - k_fib_degrade * fibril_brain) +
      k_ab40_olig * ab40_brain + k_ab42_olig * ab42_brain -
      k_olig_bbb * olig_brain - k_olig_csf * olig_brain
    d/dt(fibril_brain) <- -k_fib_clr * fibril_brain - k_fib_plaque * fibril_brain +
      (k_olig_fib * olig_brain - k_fib_degrade * fibril_brain)
    d/dt(plaque_brain) <- k_fib_plaque * fibril_brain - k_plaque_clr * plaque_brain

    # ---- CSF ----
    # The CSF clearance term appears twice in each printed CSF ODE while
    # only one copy reaches plasma; the doubled loss is what balances the
    # Table S1 CSF amounts at steady state, so it is kept as printed.
    d/dt(ab40_csf) <- f_ab40_csf_i * k_ab40_cl * ab40_brain -
      2 * k_ab40_csf_cl * ab40_csf
    d/dt(ab42_csf) <- f_ab42_dilution * f_ab42_csf_i * k_ab42_cl * ab42_brain -
      2 * k_ab42_csf_cl * ab42_csf
    d/dt(olig_csf) <- k_olig_csf * olig_brain - 2 * k_olig_csf_cl * olig_csf

    # ---- Plasma ----
    d/dt(ab40_plasma) <- f_ab40_bbb_i * k_ab40_cl * ab40_brain +
      k_ab40_csf_cl * ab40_csf + k_ab40_periph - k_ab40_pl_cl * ab40_plasma
    d/dt(ab42_plasma) <- f_ab42_dilution * f_ab42_bbb_i * k_ab42_cl * ab42_brain +
      k_ab42_csf_cl * ab42_csf - k_ab42_pl_cl * ab42_plasma + k_ab42_periph
    d/dt(olig_plasma) <- k_olig_bbb * olig_brain + k_olig_csf_cl * olig_csf -
      k_olig_pl_cl * olig_plasma

    app(0) <- bl_app
    sappb(0) <- bl_sappb
    c99(0) <- bl_c99
    ab40_brain(0) <- bl_ab40_brain
    ab42_brain(0) <- bl_ab42_brain
    olig_brain(0) <- bl_olig_brain
    fibril_brain(0) <- bl_fibril_brain
    plaque_brain(0) <- bl_plaque_brain
    ab40_csf(0) <- bl_ab40_csf
    ab42_csf(0) <- bl_ab42_csf
    olig_csf(0) <- bl_olig_csf
    ab40_plasma(0) <- bl_ab40_plasma
    ab42_plasma(0) <- bl_ab42_plasma
    olig_plasma(0) <- bl_olig_plasma

    # ---- Outputs ----
    ratio_ab4240_brain <- ab42_brain / ab40_brain
    ratio_ab4240_csf <- ab42_csf / ab40_csf
    ratio_ab4240_plasma <- ab42_plasma / ab40_plasma
    # The paper takes amyloid-PET SUVR as linearly proportional to insoluble
    # (fibril + plaque) A-beta; the proportionality constant is not
    # published, so the SUVR proxy is the percent change from baseline.
    insol_brain <- fibril_brain + plaque_brain
    pct_insol <- 100 * (insol_brain / (bl_fibril_brain + bl_plaque_brain) - 1)
  })
}
