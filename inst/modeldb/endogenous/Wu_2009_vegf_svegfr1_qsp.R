Wu_2009_vegf_svegfr1_qsp <- function() {
  description <- paste(
    "QSP. Whole-body three-compartment (blood, calf muscle, rest of body)",
    "model of endogenous VEGF and soluble VEGF receptor-1 (sVEGFR1, sFlt-1)",
    "distribution in a healthy 70 kg human (Wu 2009): 61 mass-action ODEs",
    "(28 species per tissue, 5 in blood). Tissue species are VEGF121 and",
    "VEGF165 secreted by myocytes (1:10), sVEGFR1 secreted abluminally by",
    "endothelial cells, VEGF-sVEGFR1 complexes, heparan-sulfate binding sites",
    "in the extracellular matrix and the endothelial and parenchymal basement",
    "membranes, and abluminal endothelial VEGFR1, VEGFR2 and neuropilin-1",
    "(NRP1) with all their ligand-bound and NRP1-coupled complexes; receptor",
    "insertion balances internalization so total surface densities are",
    "constant. Soluble species move between interstitium and plasma by",
    "bidirectional vascular permeability and lymphatic drainage, and are",
    "cleared from plasma first order. sVEGFR1 acts only as a ligand trap",
    "(no heterodimerization with surface receptors). Deterministic, no IIV or",
    "residual error. States are concentrations per cm^3 of whole tissue or",
    "whole blood (pM); the reported outputs are free concentrations in",
    "available interstitial fluid or plasma (pM). Time is in seconds. The",
    "model starts with all receptors and matrix sites free and no ligand,",
    "and reaches the paper's healthy control steady state after about 90",
    "days of simulated time; to perturb from that steady state, keep the",
    "equilibration period in the same solve and change parameters",
    "(e.g. lymph flow) afterwards through time-varying data columns.",
    sep = " "
  )
  reference <- paste(
    "Wu FTH, Stefanini MO, Mac Gabhann F, Popel AS.",
    "A compartment model of VEGF distribution in humans in the presence of",
    "soluble VEGF receptor-1 acting as a ligand trap.",
    "PLoS One. 2009;4(4):e5108.",
    "doi:10.1371/journal.pone.0005108.",
    sep = " "
  )
  vignette <- "Wu_2009_vegf_svegfr1_qsp"
  units <- list(time = "s", dosing = "pmol", concentration = "pM")

  # One state per molecular species per compartment. Suffix _norm = "normal"
  # (rest of body), _calf = healthy calf, _blood = blood. Species names:
  # v121 / v165 = free VEGF121 / VEGF165; sr1 = free sVEGFR1; r1 / r2 / n1 =
  # free VEGFR1 / VEGFR2 / NRP1; ecm / ebm / pbm = free matrix binding sites;
  # underscores join the members of a complex (v121_r1_n1 = VEGF121-VEGFR1-
  # NRP1, r2_v165_n1 = VEGFR2-VEGF165-NRP1). Every state is a concentration
  # per cm^3 of whole tissue (tissue states) or whole blood (blood states),
  # in pM (fmol/cm^3); divide by kav_<compartment> for the concentration in
  # available interstitial fluid or plasma.
  compartmentData <- local({
    analyteMap <- c(
      r1 = "free VEGFR1 (abluminal endothelial surface)",
      r2 = "free VEGFR2 (abluminal endothelial surface)",
      n1 = "free neuropilin-1 (abluminal endothelial surface)",
      v121_r1 = "VEGF121-VEGFR1 complex",
      v121_r2 = "VEGF121-VEGFR2 complex",
      v121_n1 = "VEGF121-NRP1 complex",
      v165_r1 = "VEGF165-VEGFR1 complex",
      v165_r2 = "VEGF165-VEGFR2 complex",
      v165_n1 = "VEGF165-NRP1 complex",
      r1_n1 = "VEGFR1-NRP1 complex",
      sr1_n1 = "sVEGFR1-NRP1 complex",
      v121_r1_n1 = "VEGF121-VEGFR1-NRP1 complex",
      v121_sr1_n1 = "VEGF121-sVEGFR1-NRP1 complex",
      r2_v165_n1 = "VEGFR2-VEGF165-NRP1 complex",
      ecm = "free VEGF/sVEGFR1-binding sites, extracellular matrix",
      ebm = "free VEGF/sVEGFR1-binding sites, endothelial basement membrane",
      pbm = "free VEGF/sVEGFR1-binding sites, parenchymal basement membrane",
      v165_ecm = "VEGF165 bound to extracellular-matrix sites",
      v165_ebm = "VEGF165 bound to endothelial-basement-membrane sites",
      v165_pbm = "VEGF165 bound to parenchymal-basement-membrane sites",
      sr1_ecm = "sVEGFR1 bound to extracellular-matrix sites",
      sr1_ebm = "sVEGFR1 bound to endothelial-basement-membrane sites",
      sr1_pbm = "sVEGFR1 bound to parenchymal-basement-membrane sites",
      v121 = "free VEGF121",
      v165 = "free VEGF165",
      v121_sr1 = "VEGF121-sVEGFR1 complex",
      v165_sr1 = "VEGF165-sVEGFR1 complex",
      sr1 = "free sVEGFR1"
    )
    bloodSpecies <- c("v121", "v165", "v121_sr1", "v165_sr1", "sr1")
    nm <- c(
      paste0(names(analyteMap), "_norm"),
      paste0(names(analyteMap), "_calf"),
      paste0(bloodSpecies, "_blood")
    )
    species <- sub("_(norm|calf|blood)$", "", nm)
    where <- c(
      norm = "rest-of-body tissue",
      calf = "calf muscle tissue",
      blood = "blood"
    )[sub("^.*_(norm|calf|blood)$", "\\1", nm)]
    Map(
      list,
      analyte = stats::setNames(paste0(analyteMap[species], ", ", where), nm),
      units = ifelse(grepl("_blood$", nm), "pM (per L of whole blood)", "pM (per L of whole tissue)"),
      specimen = ifelse(grepl("_blood$", nm), "plasma", "tissue"),
      verified = TRUE
    )
  })
  paper_specific_compartments <- names(compartmentData)

  covariateData <- list()

  population <- list(
    species = "human",
    disease_state = paste(
      "Healthy adult (virtual typical 70 kg man; resting, supine, awake",
      "control). No patient data were fitted: geometry, kinetics, transport",
      "and receptor densities come from the literature, and the VEGF and",
      "sVEGFR1 secretion rates were calibrated so that plasma free VEGF is",
      "1.5 pM and plasma free sVEGFR1 is 100 pM (healthy controls of",
      "peripheral-arterial-disease studies).",
      sep = " "
    ),
    weight_range = "70 kg reference man (65 kg non-blood tissue, 5 L blood)",
    notes = paste(
      "Calf compartment = unilateral gastrocnemius and soleus (868 cm^3);",
      "rest of body = 60,453 cm^3 treated as homogeneous skeletal muscle",
      "parameterised from vastus lateralis histology (Wu 2009 Table 1).",
      "Calibration targets: plasma free VEGF 1-2 pM (32-113 pg/mL), plasma",
      "free sVEGFR1 21-22 ng/mL (about 100 pM) (Tables 9-10).",
      sep = " "
    )
  )

  ini({
    # All values are literature-derived or calibrated (Wu 2009 Methods
    # section 3, Tables 1-11); none was estimated from individual data, so
    # every parameter is fixed. Concentration-type values are converted from
    # the paper's mol/cm^3 to pM (fmol/cm^3; factor 1e15) and second-order
    # rate constants from (mol/cm^3)^-1 s^-1 to 1/(pM*s) (factor 1e-15),
    # written as '<paper value> * <factor>' so the source number stays
    # visible. Tissue-specific in vivo values are Table 3's 'Simulation
    # Control Values (In Vivo Conversions)'.

    # --- Geometry (Table 1; Section 3.1) ---------------------------------
    vol_norm <- fixed(60453); label("Volume of the rest-of-body tissue compartment U_N (cm^3)") # Table 1 'Compartment Volume' Normal = 60453 cm3
    vol_calf <- fixed(868); label("Volume of the calf tissue compartment U_D (cm^3)") # Table 1 'Compartment Volume' Healthy Calf = 868 cm3
    vol_blood <- fixed(5000); label("Blood volume U_B (cm^3)") # Section 3.1iii: blood volume 5 L
    kav_norm <- fixed(0.03); label("Available interstitial fluid volume fraction, rest of body KAV_N (fraction)") # Table 1 'Available IF Volume Fraction' Normal = 3%
    kav_calf <- fixed(0.11); label("Available interstitial fluid volume fraction, calf KAV_D (fraction)") # Table 1 'Available IF Volume Fraction' Healthy Calf = 11%
    kav_blood <- fixed(0.6); label("Plasma volume fraction of blood KAV_B (fraction)") # Section 3.1iii: 60% of blood volume is plasma
    esav_norm <- fixed(73); label("Abluminal endothelial surface area per tissue volume, rest of body ESAV_N (cm^2/cm^3)") # Table 1 'ESAV (abluminal)' Normal = 73 cm2/cm3 tissue
    esav_calf <- fixed(44); label("Abluminal endothelial surface area per tissue volume, calf ESAV_D (cm^2/cm^3)") # Table 1 'ESAV (abluminal)' Healthy Calf = 44 cm2/cm3 tissue

    # --- Binding and coupling kinetics (Table 3) -------------------------
    kon_v_r1_norm <- fixed(1e12 * 1e-15); label("VEGF-VEGFR1 association, rest of body (1/(pM*s))") # Table 3 'VEGF binding to VEGFR1' kon Normal = 10^12 (mol/cm3 tissue)^-1 s^-1
    kon_v_r1_calf <- fixed(2.73e11 * 1e-15); label("VEGF-VEGFR1 association, calf (1/(pM*s))") # Table 3 'VEGF binding to VEGFR1' kon Healthy Calf = 2.73*10^11
    koff_v_r1 <- fixed(1e-3); label("VEGF-VEGFR1 dissociation (1/s)") # Table 3 'VEGF binding to VEGFR1' koff = 10^-3 s^-1
    kon_v_r2_norm <- fixed(3.33e11 * 1e-15); label("VEGF-VEGFR2 association, rest of body (1/(pM*s))") # Table 3 'VEGF binding to VEGFR2' kon Normal = 3.33*10^11
    kon_v_r2_calf <- fixed(9.09e10 * 1e-15); label("VEGF-VEGFR2 association, calf (1/(pM*s))") # Table 3 'VEGF binding to VEGFR2' kon Healthy Calf = 9.09*10^10
    koff_v_r2 <- fixed(1e-3); label("VEGF-VEGFR2 dissociation (1/s)") # Table 3 'VEGF binding to VEGFR2' koff = 10^-3 s^-1
    kon_v121_n1_norm <- fixed(0); label("VEGF121-NRP1 association, rest of body; off in the base case (1/(pM*s))") # Table 3 'VEGF121 binding to NRP1' kon = 0; Section 2.1ii 'turned off (kon = 0) in the base case'
    kon_v121_n1_calf <- fixed(0); label("VEGF121-NRP1 association, calf; off in the base case (1/(pM*s))") # Table 3 'VEGF121 binding to NRP1' kon = 0
    koff_v121_n1 <- fixed(0); label("VEGF121-NRP1 dissociation; off in the base case (1/s)") # Table 3 'VEGF121 binding to NRP1' koff = 0
    kon_v165_n1_norm <- fixed(1.04e11 * 1e-15); label("VEGF165-NRP1 association, rest of body (1/(pM*s))") # Table 3 'VEGF165 binding to NRP1' kon Normal = 1.04*10^11
    kon_v165_n1_calf <- fixed(2.84e10 * 1e-15); label("VEGF165-NRP1 association, calf (1/(pM*s))") # Table 3 'VEGF165 binding to NRP1' kon Healthy Calf = 2.84*10^10
    koff_v165_n1 <- fixed(1e-3); label("VEGF165-NRP1 dissociation (1/s)") # Table 3 'VEGF165 binding to NRP1' koff = 10^-3 s^-1
    kon_v165_m_norm <- fixed(1.4e10 * 1e-15); label("VEGF165-matrix site association, rest of body (1/(pM*s))") # Table 3 'VEGF165 binding to GAG' kon Normal = 1.4*10^10
    kon_v165_m_calf <- fixed(3.82e9 * 1e-15); label("VEGF165-matrix site association, calf (1/(pM*s))") # Table 3 'VEGF165 binding to GAG' kon Healthy Calf = 3.82*10^9
    koff_v165_m <- fixed(1e-2); label("VEGF165-matrix site dissociation (1/s)") # Table 3 'VEGF165 binding to GAG' koff = 10^-2 s^-1
    kc_r1_n1_norm <- fixed(1.37e12 * 1e-15); label("VEGFR1-NRP1 coupling, rest of body (1/(pM*s))") # Table 3 'VEGFR1 coupling to NRP1' kc Normal = 1.37*10^12
    kc_r1_n1_calf <- fixed(2.27e12 * 1e-15); label("VEGFR1-NRP1 coupling, calf (1/(pM*s))") # Table 3 'VEGFR1 coupling to NRP1' kc Healthy Calf = 2.27*10^12
    kdissoc_r1_n1 <- fixed(1e-2); label("VEGFR1-NRP1 decoupling (1/s)") # Table 3 'VEGFR1 coupling to NRP1' kdissoc = 10^-2 s^-1
    kc_n1_vr2_norm <- fixed(4.25e11 * 1e-15); label("VEGF165-VEGFR2 coupling to NRP1, rest of body (1/(pM*s))") # Table 3 'VEGF.VEGFR2 coupling to NRP1' kc Normal = 4.25*10^11
    kc_n1_vr2_calf <- fixed(7.05e11 * 1e-15); label("VEGF165-VEGFR2 coupling to NRP1, calf (1/(pM*s))") # Table 3 'VEGF.VEGFR2 coupling to NRP1' kc Healthy Calf = 7.05*10^11
    koff_n1_vr2 <- fixed(1e-3); label("NRP1 release from VEGFR2-VEGF165-NRP1 (1/s)") # Table 3 'VEGF.VEGFR2 coupling to NRP1' koff = 10^-3 s^-1
    kc_r2_vn1_norm <- fixed(1.37e12 * 1e-15); label("VEGF165-NRP1 coupling to VEGFR2, rest of body (1/(pM*s))") # Table 3 'VEGF.NRP1 coupling to VEGFR2' kc Normal = 1.37*10^12
    kc_r2_vn1_calf <- fixed(2.27e12 * 1e-15); label("VEGF165-NRP1 coupling to VEGFR2, calf (1/(pM*s))") # Table 3 'VEGF.NRP1 coupling to VEGFR2' kc Healthy Calf = 2.27*10^12
    koff_r2_vn1 <- fixed(1e-3); label("VEGFR2 release from VEGFR2-VEGF165-NRP1 (1/s)") # Table 3 'VEGF.NRP1 coupling to VEGFR2' koff = 10^-3 s^-1
    kon_v_sr1_norm <- fixed(1e12 * 1e-15); label("VEGF-sVEGFR1 association, rest of body (1/(pM*s))") # Table 3 'VEGF binding to sVEGFR1' kon Normal = 10^12
    kon_v_sr1_calf <- fixed(2.73e11 * 1e-15); label("VEGF-sVEGFR1 association, calf (1/(pM*s))") # Table 3 'VEGF binding to sVEGFR1' kon Healthy Calf = 2.73*10^11
    kon_v_sr1_blood <- fixed(5e10 * 1e-15); label("VEGF-sVEGFR1 association, blood (1/(pM*s))") # Table 3 'VEGF binding to sVEGFR1' kon Blood = 5*10^10 (mol/cm3 blood)^-1 s^-1
    koff_v_sr1 <- fixed(1e-3); label("VEGF-sVEGFR1 dissociation (1/s)") # Table 3 'VEGF binding to sVEGFR1' koff = 10^-3 s^-1
    kon_sr1_m_norm <- fixed(1.4e10 * 1e-15); label("sVEGFR1-matrix site association, rest of body (1/(pM*s))") # Table 3 'sVEGFR1 coupling to GAG' kon Normal = 1.4*10^10
    kon_sr1_m_calf <- fixed(3.82e9 * 1e-15); label("sVEGFR1-matrix site association, calf (1/(pM*s))") # Table 3 'sVEGFR1 coupling to GAG' kon Healthy Calf = 3.82*10^9
    koff_sr1_m <- fixed(1e-2); label("sVEGFR1-matrix site dissociation (1/s)") # Table 3 'sVEGFR1 coupling to GAG' koff = 10^-2 s^-1
    kon_sr1_n1_norm <- fixed(1.85e11 * 1e-15); label("sVEGFR1-NRP1 association, rest of body (1/(pM*s))") # Table 3 'sVEGFR1 binding to NRP1' kon Normal = 1.85*10^11
    kon_sr1_n1_calf <- fixed(5.05e10 * 1e-15); label("sVEGFR1-NRP1 association, calf (1/(pM*s))") # Table 3 'sVEGFR1 binding to NRP1' kon Healthy Calf = 5.05*10^10
    koff_sr1_n1 <- fixed(1e-2); label("sVEGFR1-NRP1 dissociation (1/s)") # Table 3 'sVEGFR1 binding to NRP1' koff = 10^-2 s^-1
    kint_free <- fixed(2.8e-4); label("Internalization rate of free surface receptors (1/s)") # Table 3 'Receptor Internalization' kint,R = 2.8*10^-4 s^-1
    kint_bound <- fixed(2.8e-4); label("Internalization rate of receptor complexes (1/s)") # Table 3 'Receptor Internalization' kint,C = 2.8*10^-4 s^-1

    # --- Total surface receptor densities (Table 8: 10,000 per EC) ------
    r1tot_norm <- fixed(1.21e-13 * 1e15); label("Total VEGFR1 density, rest of body (pM)") # Table 8 VEGFR1 Simulation Control 1.21*10^-13 moles/cm3 Normal tissue
    r1tot_calf <- fixed(7.31e-14 * 1e15); label("Total VEGFR1 density, calf (pM)") # Table 8 VEGFR1 Simulation Control 7.31*10^-14 moles/cm3 Healthy Calf tissue
    r2tot_norm <- fixed(1.21e-13 * 1e15); label("Total VEGFR2 density, rest of body (pM)") # Table 8 VEGFR2 Simulation Control 1.21*10^-13 moles/cm3 Normal tissue
    r2tot_calf <- fixed(7.31e-14 * 1e15); label("Total VEGFR2 density, calf (pM)") # Table 8 VEGFR2 Simulation Control 7.31*10^-14 moles/cm3 Healthy Calf tissue
    n1tot_norm <- fixed(1.21e-13 * 1e15); label("Total NRP1 density, rest of body (pM)") # Table 8 NRP1 Simulation Control 1.21*10^-13 moles/cm3 Normal tissue
    n1tot_calf <- fixed(7.31e-14 * 1e15); label("Total NRP1 density, calf (pM)") # Table 8 NRP1 Simulation Control 7.31*10^-14 moles/cm3 Healthy Calf tissue

    # --- Total interstitial matrix binding-site densities (Table 7) ----
    ecmtot_norm <- fixed(2.15e-11 * 1e15); label("Total ECM binding sites, rest of body (pM)") # Table 7 [VEGF-binding HS-sites]ECM 2.15*10^-11 moles/cm3 Normal tissue
    ecmtot_calf <- fixed(8.24e-11 * 1e15); label("Total ECM binding sites, calf (pM)") # Table 7 [VEGF-binding HS-sites]ECM 8.24*10^-11 moles/cm3 Healthy Calf tissue
    ebmtot_norm <- fixed(2e-12 * 1e15); label("Total EBM binding sites, rest of body (pM)") # Table 7 [VEGF-binding HS-sites]EBM 2*10^-12 moles/cm3 Normal tissue
    ebmtot_calf <- fixed(4e-12 * 1e15); label("Total EBM binding sites, calf (pM)") # Table 7 [VEGF-binding HS-sites]EBM 4*10^-12 moles/cm3 Healthy Calf tissue
    pbmtot_norm <- fixed(2e-11 * 1e15); label("Total PBM binding sites, rest of body (pM)") # Table 7 [VEGF-binding HS-sites]PBM 2*10^-11 moles/cm3 Normal tissue
    pbmtot_calf <- fixed(4.8e-11 * 1e15); label("Total PBM binding sites, calf (pM)") # Table 7 [VEGF-binding HS-sites]PBM 4.8*10^-11 moles/cm3 Healthy Calf tissue

    # --- Secretion rates (Table 11; calibrated, Results section 1.1) ---
    ksyn_v121_norm <- fixed(8.44e-19 * 1e15); label("VEGF121 secretion rate, rest of body (pM/s)") # Table 11 qV121 Normal = 8.44*10^-19 mole/(cm3 tissue)/s (0.0175 molecule/MD/s)
    ksyn_v121_calf <- fixed(5.21e-19 * 1e15); label("VEGF121 secretion rate, calf (pM/s)") # Table 11 qV121 Healthy Calf = 5.21*10^-19 mole/(cm3 tissue)/s (0.0105 molecule/MD/s)
    ksyn_v165_norm <- fixed(8.44e-18 * 1e15); label("VEGF165 secretion rate, rest of body (pM/s)") # Table 11 qV165 Normal = 8.44*10^-18 mole/(cm3 tissue)/s (0.175 molecule/MD/s)
    ksyn_v165_calf <- fixed(5.21e-18 * 1e15); label("VEGF165 secretion rate, calf (pM/s)") # Table 11 qV165 Healthy Calf 0.105 molecule/MD/s x 4.96*10^-17 (Table 11 footnote) = 5.21*10^-18; the table's printed 2.73*10^-18 is a misprint (see vignette)
    ksyn_sr1_norm <- fixed(1.30e-19 * 1e15); label("sVEGFR1 secretion rate, rest of body (pM/s)") # Table 11 qsR1 Normal = 1.30*10^-19 mole/(cm3 tissue)/s (0.0107 molecule/EC/s)
    ksyn_sr1_calf <- fixed(1.53e-19 * 1e15); label("sVEGFR1 secretion rate, calf (pM/s)") # Table 11 qsR1 Healthy Calf = 1.53*10^-19 mole/(cm3 tissue)/s (0.021 molecule/EC/s)

    # --- Transport (Tables 4-6; supine awake control) ------------------
    kp_v <- fixed(4.29e-8); label("Basal vascular permeability of free VEGF (cm/s)") # Table 4 VEGF kP = 4.29*10^-8 cm/s
    kp_sr1 <- fixed(1.86e-8); label("Basal vascular permeability of free sVEGFR1 (cm/s)") # Table 4 sVEGFR1 kP = 1.86*10^-8 cm/s
    kp_v_sr1 <- fixed(1.86e-8); label("Basal vascular permeability of the sVEGFR1-VEGF complex (cm/s)") # Table 4 sVEGFR1-VEGF kP = 1.86*10^-8 cm/s
    kcl_v <- fixed(1.08e-3); label("Plasma clearance rate constant of free VEGF (1/s)") # Table 4 VEGF kC = 1.08*10^-3 s^-1
    kcl_sr1 <- fixed(5e-6); label("Plasma clearance rate constant of free sVEGFR1 (1/s)") # Table 4 sVEGFR1 kC = 5*10^-6 s^-1
    kcl_v_sr1 <- fixed(5e-6); label("Plasma clearance rate constant of the sVEGFR1-VEGF complex (1/s)") # Table 4 sVEGFR1-VEGF kC = 5*10^-6 s^-1
    klymph_norm <- fixed(0.14); label("Lymph flow from the rest-of-body interstitium kL_N (cm^3/s)") # Table 6 'Supine; Day/Awake (Simulation Control) kL' Normal = 0.14 cm3/s
    klymph_calf <- fixed(0.0026); label("Lymph flow from the calf interstitium kL_D (cm^3/s)") # Table 6 'Supine; Day/Awake (Simulation Control) kL' Calf = 0.0026 cm3/s
    gamma_norm <- fixed(1); label("Endothelial surface-area recruitment factor, rest of body (unitless)") # Table 5 Supine (control) Normal = 1
    gamma_calf <- fixed(1); label("Endothelial surface-area recruitment factor, calf (unitless)") # Table 5 Supine (control) Healthy Calf = 1
  })

  model({
    # Wu 2009 Methods section 2 (Eqs. 7-39), written as reaction fluxes
    # j_<reaction>_<compartment> (net forward rate, pM/s on the whole-tissue
    # or whole-blood basis) so that each flux appears with the opposite sign
    # in the states it connects. Names in the j_ fluxes: j_a_b = a + b <->
    # a_b; j_ab_c and j_a_bc distinguish the two routes to a triple complex
    # (e.g. j_v121r1_n1 = VEGF121-VEGFR1 + NRP1, j_v121_r1n1 = VEGF121 +
    # VEGFR1-NRP1). Equal forward and reverse vascular permeabilities
    # (Table 4 'bidirectional') let each permeability term be written once.

    # ---- "Normal" compartment: rest of the body (suffix _norm) ---------------------------------------
    # Interstitial-matrix binding (Eqs. 7-15). One affinity per tissue for
    # all three matrix regions (ECM, EBM, PBM); sVEGFR1 competes for the
    # same sites as VEGF165 (Section 3.4i).
    j_v165_ecm_norm <- kon_v165_m_norm * v165_norm * ecm_norm - koff_v165_m * v165_ecm_norm
    j_v165_ebm_norm <- kon_v165_m_norm * v165_norm * ebm_norm - koff_v165_m * v165_ebm_norm
    j_v165_pbm_norm <- kon_v165_m_norm * v165_norm * pbm_norm - koff_v165_m * v165_pbm_norm
    j_sr1_ecm_norm <- kon_sr1_m_norm * sr1_norm * ecm_norm - koff_sr1_m * sr1_ecm_norm
    j_sr1_ebm_norm <- kon_sr1_m_norm * sr1_norm * ebm_norm - koff_sr1_m * sr1_ebm_norm
    j_sr1_pbm_norm <- kon_sr1_m_norm * sr1_norm * pbm_norm - koff_sr1_m * sr1_pbm_norm

    # Abluminal endothelial-surface binding and coupling (Eqs. 16-29).
    j_v121_r1_norm <- kon_v_r1_norm * v121_norm * r1_norm - koff_v_r1 * v121_r1_norm
    j_v165_r1_norm <- kon_v_r1_norm * v165_norm * r1_norm - koff_v_r1 * v165_r1_norm
    j_v121_r2_norm <- kon_v_r2_norm * v121_norm * r2_norm - koff_v_r2 * v121_r2_norm
    j_v165_r2_norm <- kon_v_r2_norm * v165_norm * r2_norm - koff_v_r2 * v165_r2_norm
    j_v121_n1_norm <- kon_v121_n1_norm * v121_norm * n1_norm - koff_v121_n1 * v121_n1_norm
    j_v165_n1_norm <- kon_v165_n1_norm * v165_norm * n1_norm - koff_v165_n1 * v165_n1_norm
    j_r1_n1_norm <- kc_r1_n1_norm * r1_norm * n1_norm - kdissoc_r1_n1 * r1_n1_norm
    j_v121r1_n1_norm <- kc_r1_n1_norm * v121_r1_norm * n1_norm - kdissoc_r1_n1 * v121_r1_n1_norm
    j_v121_r1n1_norm <- kon_v_r1_norm * v121_norm * r1_n1_norm - koff_v_r1 * v121_r1_n1_norm
    j_v165r2_n1_norm <- kc_n1_vr2_norm * n1_norm * v165_r2_norm - koff_n1_vr2 * r2_v165_n1_norm
    j_r2_v165n1_norm <- kc_r2_vn1_norm * r2_norm * v165_n1_norm - koff_r2_vn1 * r2_v165_n1_norm
    j_sr1_n1_norm <- kon_sr1_n1_norm * sr1_norm * n1_norm - koff_sr1_n1 * sr1_n1_norm
    j_v121sr1_n1_norm <- kon_sr1_n1_norm * v121_sr1_norm * n1_norm - koff_sr1_n1 * v121_sr1_n1_norm
    j_v121_sr1n1_norm <- kon_v_sr1_norm * v121_norm * sr1_n1_norm - koff_v_sr1 * v121_sr1_n1_norm

    # Interstitial VEGF-sVEGFR1 trapping (Eqs. 30-34).
    j_v121_sr1_norm <- kon_v_sr1_norm * v121_norm * sr1_norm - koff_v_sr1 * v121_sr1_norm
    j_v165_sr1_norm <- kon_v_sr1_norm * v165_norm * sr1_norm - koff_v_sr1 * v165_sr1_norm

    # Transport rate constants on the tissue side (1/s). Lymph: kL / U_j
    # applied to the available-fluid concentration [X]_j / KAV_j. Vascular
    # permeability: gamma_j * S_jB / U_j = gamma_j * ESAV_j, because
    # S_jB = ESAV_j * U_j.
    klym_norm <- klymph_norm / vol_norm
    kperm_norm <- gamma_norm * esav_norm

    # Surface species (Eqs. 16-29). Insertion of free receptors
    # s_R = kint * [R_total] holds the total surface density constant.
    d/dt(r1_norm) <- kint_free * r1tot_norm - kint_free * r1_norm - j_v121_r1_norm - j_v165_r1_norm - j_r1_n1_norm
    d/dt(r2_norm) <- kint_free * r2tot_norm - kint_free * r2_norm - j_v121_r2_norm - j_v165_r2_norm - j_r2_v165n1_norm
    d/dt(n1_norm) <- kint_free * n1tot_norm - kint_free * n1_norm - j_v121_n1_norm - j_v165_n1_norm - j_v165r2_n1_norm -
      j_r1_n1_norm - j_v121r1_n1_norm - j_sr1_n1_norm - j_v121sr1_n1_norm
    d/dt(v121_r1_norm) <- -kint_bound * v121_r1_norm + j_v121_r1_norm - j_v121r1_n1_norm
    d/dt(v121_r2_norm) <- -kint_bound * v121_r2_norm + j_v121_r2_norm
    d/dt(v121_n1_norm) <- -kint_bound * v121_n1_norm + j_v121_n1_norm
    d/dt(v165_r1_norm) <- -kint_bound * v165_r1_norm + j_v165_r1_norm
    d/dt(v165_r2_norm) <- -kint_bound * v165_r2_norm + j_v165_r2_norm - j_v165r2_n1_norm
    d/dt(v165_n1_norm) <- -kint_bound * v165_n1_norm + j_v165_n1_norm - j_r2_v165n1_norm
    d/dt(r1_n1_norm) <- -kint_bound * r1_n1_norm + j_r1_n1_norm - j_v121_r1n1_norm
    d/dt(sr1_n1_norm) <- -kint_bound * sr1_n1_norm + j_sr1_n1_norm - j_v121_sr1n1_norm
    d/dt(v121_r1_n1_norm) <- -kint_bound * v121_r1_n1_norm + j_v121r1_n1_norm + j_v121_r1n1_norm
    d/dt(v121_sr1_n1_norm) <- -kint_bound * v121_sr1_n1_norm + j_v121sr1_n1_norm + j_v121_sr1n1_norm
    d/dt(r2_v165_n1_norm) <- -kint_bound * r2_v165_n1_norm + j_v165r2_n1_norm + j_r2_v165n1_norm

    # Interstitial-matrix species (Eqs. 7-15).
    d/dt(ecm_norm) <- -j_v165_ecm_norm - j_sr1_ecm_norm
    d/dt(ebm_norm) <- -j_v165_ebm_norm - j_sr1_ebm_norm
    d/dt(pbm_norm) <- -j_v165_pbm_norm - j_sr1_pbm_norm
    d/dt(v165_ecm_norm) <- j_v165_ecm_norm
    d/dt(v165_ebm_norm) <- j_v165_ebm_norm
    d/dt(v165_pbm_norm) <- j_v165_pbm_norm
    d/dt(sr1_ecm_norm) <- j_sr1_ecm_norm
    d/dt(sr1_ebm_norm) <- j_sr1_ebm_norm
    d/dt(sr1_pbm_norm) <- j_sr1_pbm_norm

    # Interstitial-fluid soluble species (Eqs. 30-34).
    d/dt(v121_norm) <- ksyn_v121_norm - klym_norm * v121_norm / kav_norm +
      kperm_norm * kp_v * (v121_blood / kav_blood - v121_norm / kav_norm) -
      j_v121_n1_norm - j_v121_r2_norm - j_v121_r1_norm - j_v121_r1n1_norm - j_v121_sr1_norm - j_v121_sr1n1_norm
    d/dt(v165_norm) <- ksyn_v165_norm - klym_norm * v165_norm / kav_norm +
      kperm_norm * kp_v * (v165_blood / kav_blood - v165_norm / kav_norm) -
      j_v165_ebm_norm - j_v165_ecm_norm - j_v165_pbm_norm - j_v165_n1_norm - j_v165_r2_norm - j_v165_r1_norm - j_v165_sr1_norm
    d/dt(v121_sr1_norm) <- -klym_norm * v121_sr1_norm / kav_norm +
      kperm_norm * kp_v_sr1 * (v121_sr1_blood / kav_blood - v121_sr1_norm / kav_norm) +
      j_v121_sr1_norm - j_v121sr1_n1_norm
    d/dt(v165_sr1_norm) <- -klym_norm * v165_sr1_norm / kav_norm +
      kperm_norm * kp_v_sr1 * (v165_sr1_blood / kav_blood - v165_sr1_norm / kav_norm) +
      j_v165_sr1_norm
    d/dt(sr1_norm) <- ksyn_sr1_norm - klym_norm * sr1_norm / kav_norm +
      kperm_norm * kp_sr1 * (sr1_blood / kav_blood - sr1_norm / kav_norm) -
      j_sr1_ebm_norm - j_sr1_ecm_norm - j_sr1_pbm_norm - j_sr1_n1_norm - j_v121_sr1_norm - j_v165_sr1_norm

    # ---- "Calf" compartment: unilateral gastrocnemius + soleus (suffix _calf) ---------------------------------------
    # Interstitial-matrix binding (Eqs. 7-15). One affinity per tissue for
    # all three matrix regions (ECM, EBM, PBM); sVEGFR1 competes for the
    # same sites as VEGF165 (Section 3.4i).
    j_v165_ecm_calf <- kon_v165_m_calf * v165_calf * ecm_calf - koff_v165_m * v165_ecm_calf
    j_v165_ebm_calf <- kon_v165_m_calf * v165_calf * ebm_calf - koff_v165_m * v165_ebm_calf
    j_v165_pbm_calf <- kon_v165_m_calf * v165_calf * pbm_calf - koff_v165_m * v165_pbm_calf
    j_sr1_ecm_calf <- kon_sr1_m_calf * sr1_calf * ecm_calf - koff_sr1_m * sr1_ecm_calf
    j_sr1_ebm_calf <- kon_sr1_m_calf * sr1_calf * ebm_calf - koff_sr1_m * sr1_ebm_calf
    j_sr1_pbm_calf <- kon_sr1_m_calf * sr1_calf * pbm_calf - koff_sr1_m * sr1_pbm_calf

    # Abluminal endothelial-surface binding and coupling (Eqs. 16-29).
    j_v121_r1_calf <- kon_v_r1_calf * v121_calf * r1_calf - koff_v_r1 * v121_r1_calf
    j_v165_r1_calf <- kon_v_r1_calf * v165_calf * r1_calf - koff_v_r1 * v165_r1_calf
    j_v121_r2_calf <- kon_v_r2_calf * v121_calf * r2_calf - koff_v_r2 * v121_r2_calf
    j_v165_r2_calf <- kon_v_r2_calf * v165_calf * r2_calf - koff_v_r2 * v165_r2_calf
    j_v121_n1_calf <- kon_v121_n1_calf * v121_calf * n1_calf - koff_v121_n1 * v121_n1_calf
    j_v165_n1_calf <- kon_v165_n1_calf * v165_calf * n1_calf - koff_v165_n1 * v165_n1_calf
    j_r1_n1_calf <- kc_r1_n1_calf * r1_calf * n1_calf - kdissoc_r1_n1 * r1_n1_calf
    j_v121r1_n1_calf <- kc_r1_n1_calf * v121_r1_calf * n1_calf - kdissoc_r1_n1 * v121_r1_n1_calf
    j_v121_r1n1_calf <- kon_v_r1_calf * v121_calf * r1_n1_calf - koff_v_r1 * v121_r1_n1_calf
    j_v165r2_n1_calf <- kc_n1_vr2_calf * n1_calf * v165_r2_calf - koff_n1_vr2 * r2_v165_n1_calf
    j_r2_v165n1_calf <- kc_r2_vn1_calf * r2_calf * v165_n1_calf - koff_r2_vn1 * r2_v165_n1_calf
    j_sr1_n1_calf <- kon_sr1_n1_calf * sr1_calf * n1_calf - koff_sr1_n1 * sr1_n1_calf
    j_v121sr1_n1_calf <- kon_sr1_n1_calf * v121_sr1_calf * n1_calf - koff_sr1_n1 * v121_sr1_n1_calf
    j_v121_sr1n1_calf <- kon_v_sr1_calf * v121_calf * sr1_n1_calf - koff_v_sr1 * v121_sr1_n1_calf

    # Interstitial VEGF-sVEGFR1 trapping (Eqs. 30-34).
    j_v121_sr1_calf <- kon_v_sr1_calf * v121_calf * sr1_calf - koff_v_sr1 * v121_sr1_calf
    j_v165_sr1_calf <- kon_v_sr1_calf * v165_calf * sr1_calf - koff_v_sr1 * v165_sr1_calf

    # Transport rate constants on the tissue side (1/s). Lymph: kL / U_j
    # applied to the available-fluid concentration [X]_j / KAV_j. Vascular
    # permeability: gamma_j * S_jB / U_j = gamma_j * ESAV_j, because
    # S_jB = ESAV_j * U_j.
    klym_calf <- klymph_calf / vol_calf
    kperm_calf <- gamma_calf * esav_calf

    # Surface species (Eqs. 16-29). Insertion of free receptors
    # s_R = kint * [R_total] holds the total surface density constant.
    d/dt(r1_calf) <- kint_free * r1tot_calf - kint_free * r1_calf - j_v121_r1_calf - j_v165_r1_calf - j_r1_n1_calf
    d/dt(r2_calf) <- kint_free * r2tot_calf - kint_free * r2_calf - j_v121_r2_calf - j_v165_r2_calf - j_r2_v165n1_calf
    d/dt(n1_calf) <- kint_free * n1tot_calf - kint_free * n1_calf - j_v121_n1_calf - j_v165_n1_calf - j_v165r2_n1_calf -
      j_r1_n1_calf - j_v121r1_n1_calf - j_sr1_n1_calf - j_v121sr1_n1_calf
    d/dt(v121_r1_calf) <- -kint_bound * v121_r1_calf + j_v121_r1_calf - j_v121r1_n1_calf
    d/dt(v121_r2_calf) <- -kint_bound * v121_r2_calf + j_v121_r2_calf
    d/dt(v121_n1_calf) <- -kint_bound * v121_n1_calf + j_v121_n1_calf
    d/dt(v165_r1_calf) <- -kint_bound * v165_r1_calf + j_v165_r1_calf
    d/dt(v165_r2_calf) <- -kint_bound * v165_r2_calf + j_v165_r2_calf - j_v165r2_n1_calf
    d/dt(v165_n1_calf) <- -kint_bound * v165_n1_calf + j_v165_n1_calf - j_r2_v165n1_calf
    d/dt(r1_n1_calf) <- -kint_bound * r1_n1_calf + j_r1_n1_calf - j_v121_r1n1_calf
    d/dt(sr1_n1_calf) <- -kint_bound * sr1_n1_calf + j_sr1_n1_calf - j_v121_sr1n1_calf
    d/dt(v121_r1_n1_calf) <- -kint_bound * v121_r1_n1_calf + j_v121r1_n1_calf + j_v121_r1n1_calf
    d/dt(v121_sr1_n1_calf) <- -kint_bound * v121_sr1_n1_calf + j_v121sr1_n1_calf + j_v121_sr1n1_calf
    d/dt(r2_v165_n1_calf) <- -kint_bound * r2_v165_n1_calf + j_v165r2_n1_calf + j_r2_v165n1_calf

    # Interstitial-matrix species (Eqs. 7-15).
    d/dt(ecm_calf) <- -j_v165_ecm_calf - j_sr1_ecm_calf
    d/dt(ebm_calf) <- -j_v165_ebm_calf - j_sr1_ebm_calf
    d/dt(pbm_calf) <- -j_v165_pbm_calf - j_sr1_pbm_calf
    d/dt(v165_ecm_calf) <- j_v165_ecm_calf
    d/dt(v165_ebm_calf) <- j_v165_ebm_calf
    d/dt(v165_pbm_calf) <- j_v165_pbm_calf
    d/dt(sr1_ecm_calf) <- j_sr1_ecm_calf
    d/dt(sr1_ebm_calf) <- j_sr1_ebm_calf
    d/dt(sr1_pbm_calf) <- j_sr1_pbm_calf

    # Interstitial-fluid soluble species (Eqs. 30-34).
    d/dt(v121_calf) <- ksyn_v121_calf - klym_calf * v121_calf / kav_calf +
      kperm_calf * kp_v * (v121_blood / kav_blood - v121_calf / kav_calf) -
      j_v121_n1_calf - j_v121_r2_calf - j_v121_r1_calf - j_v121_r1n1_calf - j_v121_sr1_calf - j_v121_sr1n1_calf
    d/dt(v165_calf) <- ksyn_v165_calf - klym_calf * v165_calf / kav_calf +
      kperm_calf * kp_v * (v165_blood / kav_blood - v165_calf / kav_calf) -
      j_v165_ebm_calf - j_v165_ecm_calf - j_v165_pbm_calf - j_v165_n1_calf - j_v165_r2_calf - j_v165_r1_calf - j_v165_sr1_calf
    d/dt(v121_sr1_calf) <- -klym_calf * v121_sr1_calf / kav_calf +
      kperm_calf * kp_v_sr1 * (v121_sr1_blood / kav_blood - v121_sr1_calf / kav_calf) +
      j_v121_sr1_calf - j_v121sr1_n1_calf
    d/dt(v165_sr1_calf) <- -klym_calf * v165_sr1_calf / kav_calf +
      kperm_calf * kp_v_sr1 * (v165_sr1_blood / kav_blood - v165_sr1_calf / kav_calf) +
      j_v165_sr1_calf
    d/dt(sr1_calf) <- ksyn_sr1_calf - klym_calf * sr1_calf / kav_calf +
      kperm_calf * kp_sr1 * (sr1_blood / kav_blood - sr1_calf / kav_calf) -
      j_sr1_ebm_calf - j_sr1_ecm_calf - j_sr1_pbm_calf - j_sr1_n1_calf - j_v121_sr1_calf - j_v165_sr1_calf

    # ---- Blood (Eqs. 35-39) ---------------------------------------------
    # Lymph and permeability flows leave a tissue of volume U_j and enter
    # blood of volume U_B, so the blood-side rate constants carry
    # U_j / U_B (lymph: kL_j / U_B; permeability: gamma_j * ESAV_j * U_j / U_B).
    jb_v121_sr1 <- kon_v_sr1_blood * v121_blood * sr1_blood - koff_v_sr1 * v121_sr1_blood
    jb_v165_sr1 <- kon_v_sr1_blood * v165_blood * sr1_blood - koff_v_sr1 * v165_sr1_blood
    blym_norm <- klymph_norm / vol_blood
    blym_calf <- klymph_calf / vol_blood
    bperm_norm <- gamma_norm * esav_norm * vol_norm / vol_blood
    bperm_calf <- gamma_calf * esav_calf * vol_calf / vol_blood

    d/dt(v121_blood) <- -kcl_v * v121_blood +
      blym_norm * v121_norm / kav_norm + blym_calf * v121_calf / kav_calf +
      bperm_norm * kp_v * (v121_norm / kav_norm - v121_blood / kav_blood) +
      bperm_calf * kp_v * (v121_calf / kav_calf - v121_blood / kav_blood) -
      jb_v121_sr1
    d/dt(v165_blood) <- -kcl_v * v165_blood +
      blym_norm * v165_norm / kav_norm + blym_calf * v165_calf / kav_calf +
      bperm_norm * kp_v * (v165_norm / kav_norm - v165_blood / kav_blood) +
      bperm_calf * kp_v * (v165_calf / kav_calf - v165_blood / kav_blood) -
      jb_v165_sr1
    d/dt(v121_sr1_blood) <- -kcl_v_sr1 * v121_sr1_blood +
      blym_norm * v121_sr1_norm / kav_norm + blym_calf * v121_sr1_calf / kav_calf +
      bperm_norm * kp_v_sr1 * (v121_sr1_norm / kav_norm - v121_sr1_blood / kav_blood) +
      bperm_calf * kp_v_sr1 * (v121_sr1_calf / kav_calf - v121_sr1_blood / kav_blood) +
      jb_v121_sr1
    d/dt(v165_sr1_blood) <- -kcl_v_sr1 * v165_sr1_blood +
      blym_norm * v165_sr1_norm / kav_norm + blym_calf * v165_sr1_calf / kav_calf +
      bperm_norm * kp_v_sr1 * (v165_sr1_norm / kav_norm - v165_sr1_blood / kav_blood) +
      bperm_calf * kp_v_sr1 * (v165_sr1_calf / kav_calf - v165_sr1_blood / kav_blood) +
      jb_v165_sr1
    d/dt(sr1_blood) <- -kcl_sr1 * sr1_blood +
      blym_norm * sr1_norm / kav_norm + blym_calf * sr1_calf / kav_calf +
      bperm_norm * kp_sr1 * (sr1_norm / kav_norm - sr1_blood / kav_blood) +
      bperm_calf * kp_sr1 * (sr1_calf / kav_calf - sr1_blood / kav_blood) -
      jb_v121_sr1 - jb_v165_sr1

    # ---- Initial state ---------------------------------------------------
    # The paper solves the system as an initial-value problem to steady
    # state without stating the starting point. Here every receptor and
    # matrix site starts free and every ligand-containing species starts at
    # zero; the healthy steady state is reached (all outputs within 0.1%)
    # after about 90 days of simulated time (see the vignette). These
    # assignments take precedence over an `inits` argument to rxSolve().
    r1_norm(0) <- r1tot_norm
    r2_norm(0) <- r2tot_norm
    n1_norm(0) <- n1tot_norm
    ecm_norm(0) <- ecmtot_norm
    ebm_norm(0) <- ebmtot_norm
    pbm_norm(0) <- pbmtot_norm
    r1_calf(0) <- r1tot_calf
    r2_calf(0) <- r2tot_calf
    n1_calf(0) <- n1tot_calf
    ecm_calf(0) <- ecmtot_calf
    ebm_calf(0) <- ebmtot_calf
    pbm_calf(0) <- pbmtot_calf

    # ---- Outputs: free concentrations in available fluid (pM) -----------
    # [X]_pl = [X]_B / KAV_B and [X]_IS,j = [X]_j / KAV_j (Eqs. 1 and 6).
    # VEGF outputs sum the two isoforms, as the paper reports them.
    vegf_free_plasma <- (v121_blood + v165_blood) / kav_blood
    svegfr1_free_plasma <- sr1_blood / kav_blood
    vegf_svegfr1_plasma <- (v121_sr1_blood + v165_sr1_blood) / kav_blood
    vegf_free_is_norm <- (v121_norm + v165_norm) / kav_norm
    vegf_free_is_calf <- (v121_calf + v165_calf) / kav_calf
    svegfr1_free_is_norm <- sr1_norm / kav_norm
    svegfr1_free_is_calf <- sr1_calf / kav_calf
  })
}
