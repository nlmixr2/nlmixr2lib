Kolodkin_2013_cortisol_qsp <- function() {
  description <- paste(
    "QSP. Kolodkin 2013 human-hepatocyte gene-regulatory network of the",
    "stress hormone cortisol and its two nuclear receptors: the high-affinity",
    "glucocorticoid receptor (GR, Kd 10 nM) and the low-affinity pregnane",
    "X receptor (PXR, Kd 10 uM). Liganded GR represses GR transcription",
    "(feedback), induces PXR transcription (feed-forward) and induces the",
    "pharmacodynamic marker tyrosine aminotransferase (TAT); liganded PXR",
    "induces CYP3A4, which metabolises cortisol (the pharmacokinetic",
    "response). An optional pool of integrated low-affinity endogenous PXR",
    "ligands (L2, 100 nM by default) and a dexamethasone branch are included.",
    "Free blood cortisol exchanges with the cell and binds corticosteroid-",
    "binding globulin and albumin; a stress episode is a cortisol dose into",
    "the cort_added depot. 36 ODEs, mass-action and competitive",
    "Michaelis-Menten kinetics, time in minutes, concentrations in nM.",
    "Deterministic: no IIV and no residual error.",
    sep = " "
  )
  reference <- paste(
    "Kolodkin A, Sahin N, Phillips A, Hood SR, Bruggeman FJ, Westerhoff HV,",
    "Plant N (2013). Optimization of stress response through the nuclear",
    "receptor-mediated cortisol signalling network. Nat Commun 4:1792.",
    "doi:10.1038/ncomms2799. PMCID: PMC3644104.",
    "Rate laws and parameter values from Supplementary Tables S1 and S5,",
    "balance equations from Supplementary Tables S2 and S6, initial",
    "conditions from Supplementary Tables S4 and S8; implementation details",
    "cross-checked against the authors' SBML deposit, BioModels",
    "BIOMD0000000576.",
    sep = " "
  )
  vignette <- "Kolodkin_2013_cortisol_qsp"
  dosing <- c("cort_added", "dex_out")

  # Every state is a species of the published reaction network (gene,
  # mRNA, protein, ligand-receptor or receptor-gene complex, plasma binding
  # protein); none maps onto a canonical PK compartment role. Names are
  # snake_case forms of the Supplementary Table S2 / S6 species
  # (PXRproteinCortCYPgene -> cypgene_pxrcort, CortOUT -> cort_blood, ...).
  paper_specific_compartments <- c(
    "cort",
    "cort_added",
    "cort_blood",
    "cbg",
    "cbg_cort",
    "alb",
    "alb_cort",
    "ligand2",
    "dex",
    "dex_out",
    "gr_gene",
    "gr_mrna",
    "gr_prot",
    "gr_cort",
    "gr_dex",
    "grgene_grcort",
    "grgene_grdex",
    "pxr_gene",
    "pxr_mrna",
    "pxr_prot",
    "pxr_cort",
    "pxr_l2",
    "pxr_dex",
    "pxrgene_grcort",
    "pxrgene_grdex",
    "cyp_gene",
    "cyp_mrna",
    "cyp_prot",
    "cypgene_pxrcort",
    "cypgene_pxrl2",
    "cypgene_pxrdex",
    "tat_gene",
    "tat_mrna",
    "tat_prot",
    "tatgene_grcort",
    "tatgene_grdex"
  )

  units <- list(
    time = "min",
    dosing = "nM (a stress episode is a bolus of cortisol concentration into cort_added; dexamethasone exposure is a bolus into dex_out, which then holds that medium concentration)",
    concentration = "nM"
  )

  covariateData <- list()

  compartmentData <- list(
    cort = list(analyte = "free cortisol", units = "nM", specimen = "tissue", verified = TRUE),
    cort_added = list(
      analyte = "cortisol released into blood by a stress episode (input depot)",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    cort_blood = list(analyte = "free cortisol", units = "nM", specimen = "plasma", verified = TRUE),
    cbg = list(analyte = "free corticosteroid-binding globulin", units = "nM", specimen = "plasma", verified = TRUE),
    cbg_cort = list(
      analyte = "cortisol bound to corticosteroid-binding globulin",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    alb = list(analyte = "free albumin", units = "nM", specimen = "plasma", verified = TRUE),
    alb_cort = list(analyte = "cortisol bound to albumin", units = "nM", specimen = "plasma", verified = TRUE),
    ligand2 = list(
      analyte = "free pool of integrated low-affinity endogenous PXR ligands (L2)",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    dex = list(analyte = "free dexamethasone", units = "nM", specimen = "tissue", verified = TRUE),
    dex_out = list(
      analyte = "dexamethasone in the culture medium",
      units = "nM",
      specimen = "not applicable",
      verified = TRUE
    ),
    gr_gene = list(analyte = "free GR gene response elements", units = "nM", specimen = "tissue", verified = TRUE),
    gr_mrna = list(analyte = "GR mRNA", units = "nM", specimen = "tissue", verified = TRUE),
    gr_prot = list(analyte = "free GR protein", units = "nM", specimen = "tissue", verified = TRUE),
    gr_cort = list(analyte = "GR-cortisol complex (GRL)", units = "nM", specimen = "tissue", verified = TRUE),
    gr_dex = list(analyte = "GR-dexamethasone complex", units = "nM", specimen = "tissue", verified = TRUE),
    grgene_grcort = list(analyte = "GRL bound to the GR gene", units = "nM", specimen = "tissue", verified = TRUE),
    grgene_grdex = list(
      analyte = "GR-dexamethasone bound to the GR gene",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    pxr_gene = list(analyte = "free PXR gene response elements", units = "nM", specimen = "tissue", verified = TRUE),
    pxr_mrna = list(analyte = "PXR mRNA", units = "nM", specimen = "tissue", verified = TRUE),
    pxr_prot = list(analyte = "free PXR protein", units = "nM", specimen = "tissue", verified = TRUE),
    pxr_cort = list(analyte = "PXR-cortisol complex (PXRL)", units = "nM", specimen = "tissue", verified = TRUE),
    pxr_l2 = list(analyte = "PXR-L2 complex (PXRL2)", units = "nM", specimen = "tissue", verified = TRUE),
    pxr_dex = list(analyte = "PXR-dexamethasone complex", units = "nM", specimen = "tissue", verified = TRUE),
    pxrgene_grcort = list(analyte = "GRL bound to the PXR gene", units = "nM", specimen = "tissue", verified = TRUE),
    pxrgene_grdex = list(
      analyte = "GR-dexamethasone bound to the PXR gene",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    cyp_gene = list(
      analyte = "free CYP3A4 gene response elements",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    cyp_mrna = list(analyte = "CYP3A4 mRNA", units = "nM", specimen = "tissue", verified = TRUE),
    cyp_prot = list(analyte = "CYP3A4 protein", units = "nM", specimen = "tissue", verified = TRUE),
    cypgene_pxrcort = list(
      analyte = "PXRL bound to the CYP3A4 gene",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    cypgene_pxrl2 = list(
      analyte = "PXRL2 bound to the CYP3A4 gene",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    cypgene_pxrdex = list(
      analyte = "PXR-dexamethasone bound to the CYP3A4 gene",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    tat_gene = list(analyte = "free TAT gene response elements", units = "nM", specimen = "tissue", verified = TRUE),
    tat_mrna = list(analyte = "TAT mRNA", units = "nM", specimen = "tissue", verified = TRUE),
    tat_prot = list(analyte = "TAT protein", units = "nM", specimen = "tissue", verified = TRUE),
    tatgene_grcort = list(analyte = "GRL bound to the TAT gene", units = "nM", specimen = "tissue", verified = TRUE),
    tatgene_grdex = list(
      analyte = "GR-dexamethasone bound to the TAT gene",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = NA_integer_,
    disease_state = paste(
      "Healthy human liver (hepatocyte) under normal physiology and under",
      "repeated acute stress episodes. No patients were studied.",
      sep = " "
    ),
    dose_range = paste(
      "Calibration: primary human hepatocytes exposed for 48 h to vehicle or",
      "0.1-50000 nM cortisol or dexamethasone, dosed daily (Figure 1b-c,",
      "Supplementary Figure S2). Simulation: intracellular cortisol clamped at",
      "10 nM (Figures 2-3), or stress episodes of about 130 nM cortisol added",
      "to blood, each doubling free blood cortisol, once per week, day or hour",
      "for 21 days (Figures 4-5).",
      sep = " "
    ),
    notes = paste(
      "A mechanistic systems-biology model. Rate constants come from the",
      "literature where available (Supplementary Table S1 references) and the",
      "remainder (mRNA and protein synthesis rates, liganded-receptor",
      "dissociation from response elements) were hand-fitted to transcript",
      "levels measured by TaqMan qPCR in primary human hepatocytes (cells",
      "donated by GlaxoSmithKline, quadruplicate cultures) within",
      "biologically realistic bounds; no estimation uncertainty is reported.",
      "The 1-ligand model of Supplementary Tables S1-S4 is this model with",
      "bl_ligand2 = 0.",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # Every value is a fixed constant of the published network: Kolodkin
    # 2013 took rate constants from the literature or hand-fitted them to
    # hepatocyte transcript data, and reports no estimation uncertainty.
    # Values are Supplementary Table S5 (the cortisol + L2 model), which
    # repeats Table S1 apart from switching on the L2 reactions. Time is in
    # MINUTES and concentrations in nM. Reaction numbers (v1, v2, ...) are
    # those of Supplementary Table S1.
    # =====================================================================

    # --- Basal (unliganded-gene) transcription, mRNA decay, translation and
    #     protein decay.
    ktx_cyp <- fixed(0.00321) ; label("CYP3A4 mRNA synthesis from the free gene (1/min)")  # Table S1/S5 v1 k1 = 0.00321 min-1
    kdeg_cyp_mrna <- fixed(0.04) ; label("CYP3A4 mRNA degradation (1/min)")  # Table S1/S5 v2 k2 = 0.04 min-1
    ktl_cyp <- fixed(2.5) ; label("CYP3A4 protein synthesis per mRNA (1/min)")  # Table S1/S5 v3 k3 = 2.5 min-1
    kdeg_cyp <- fixed(0.0015) ; label("CYP3A4 protein degradation (1/min)")
    # ^^ Table S1/S5 v4 print k4 = 0.0005 min-1; the authors' SBML deposit
    #    (BIOMD0000000576, reaction re4) uses 0.0015. The deposit value is
    #    used because it reproduces the paper's own numbers exactly: the
    #    cortisol-free baseline CYP3A4 protein is 111.0 nM (1-ligand) and
    #    1628.3 nM (with L2) with 0.0015, against the Results' '1,628 nM
    #    versus 111 nM'; 0.0005 gives 333 and 4885 nM. See vignette Errata.
    ktx_pxr <- fixed(5.52e-05) ; label("PXR mRNA synthesis from the free gene (1/min)")  # Table S1/S5 v5 k5 = 5.52e-05 min-1
    kdeg_pxr_mrna <- fixed(0.006) ; label("PXR mRNA degradation (1/min)")  # Table S1/S5 v6 k6 = 0.006 min-1
    ktl_pxr <- fixed(10) ; label("PXR protein synthesis per mRNA (1/min)")  # Table S1/S5 v7 k7 = 10 min-1
    kdeg_pxr <- fixed(0.003) ; label("Free PXR protein degradation (1/min)")  # Table S1/S5 v8 k8 = 0.003 min-1
    ktx_gr <- fixed(3.2e-06) ; label("GR mRNA synthesis from the free gene (1/min)")  # Table S1/S5 v9 k9 = 3.2e-06 min-1
    kdeg_gr_mrna <- fixed(0.003) ; label("GR mRNA degradation (1/min)")  # Table S1/S5 v10 k10 = 0.003 min-1
    ktl_gr <- fixed(19.98) ; label("GR protein synthesis per mRNA (1/min)")  # Table S1/S5 v11 k11 = 19.98 min-1
    kdeg_gr <- fixed(0.001) ; label("Free GR protein degradation (1/min)")  # Table S1/S5 v12 k12 = 0.001 min-1
    ktx_tat <- fixed(0.000855) ; label("TAT mRNA synthesis from the free gene (1/min)")  # Table S1/S5 v13 k13 = 0.000855 min-1
    kdeg_tat_mrna <- fixed(0.064) ; label("TAT mRNA degradation (1/min)")  # Table S1/S5 v14 k14 = 0.064 min-1
    ktl_tat <- fixed(0.5) ; label("TAT protein synthesis per mRNA (1/min)")  # Table S1/S5 v43 k43 = 0.5 min-1
    kdeg_tat <- fixed(0.012) ; label("TAT protein degradation (1/min)")  # Table S1/S5 v44 k44 = 0.012 min-1

    # --- Ligand-receptor binding (diffusion-limited association) and the
    #     degradation of each liganded receptor.
    kon_gr_cort <- fixed(60) ; label("Cortisol-GR association (1/(nM*min))")  # Table S1/S5 v15 kf15 = 60 nM-1 min-1 (diffusion limited)
    koff_gr_cort <- fixed(600) ; label("Cortisol-GR dissociation (1/min)")  # Table S1/S5 v15 kb15 = 600 min-1 (Kd = 10 nM)
    kdeg_gr_cort <- fixed(0.001) ; label("GR-cortisol complex degradation (1/min)")  # Table S1/S5 v16 k16 = 0.001 min-1
    kon_pxr_cort <- fixed(60) ; label("Cortisol-PXR association (1/(nM*min))")  # Table S1/S5 v17 kf17 = 60 nM-1 min-1
    koff_pxr_cort <- fixed(600000) ; label("Cortisol-PXR dissociation (1/min)")  # Table S1/S5 v17 kb17 = 600000 min-1 (Kd = 10 uM)
    kdeg_pxr_cort <- fixed(0.0015) ; label("PXR-cortisol complex degradation (1/min)")  # Table S1/S5 v18 k18 = 0.0015 min-1
    kon_pxr_l2 <- fixed(60) ; label("L2-PXR association (1/(nM*min))")  # Table S5 V_L2_PXR_binding kf = 60 nM-1 min-1 (Table S1: 0, L2 absent)
    koff_pxr_l2 <- fixed(600000) ; label("L2-PXR dissociation (1/min)")  # Table S5 V_L2_PXR_binding kb = 600000 min-1 (Kd = 10 uM)
    kdeg_pxr_l2 <- fixed(0.0015) ; label("PXR-L2 complex degradation (1/min)")  # Table S5 V_L2_PXR_deg k = 0.0015 min-1
    kon_gr_dex <- fixed(60) ; label("Dexamethasone-GR association (1/(nM*min))")  # Table S1/S5 v28 kf28 = 60 nM-1 min-1
    koff_gr_dex <- fixed(60) ; label("Dexamethasone-GR dissociation (1/min)")  # Table S1/S5 v28 kb28 = 60 min-1 (Kd about 1 nM)
    kdeg_gr_dex <- fixed(0.001) ; label("GR-dexamethasone complex degradation (1/min)")  # Table S1/S5 v29 k29 = 0.001 min-1
    kon_pxr_dex <- fixed(60) ; label("Dexamethasone-PXR association (1/(nM*min))")  # Table S1/S5 v30 kf30 = 60 nM-1 min-1
    koff_pxr_dex <- fixed(60000) ; label("Dexamethasone-PXR dissociation (1/min)")  # Table S1/S5 v30 kb30 = 60000 min-1 (Kd about 1000 nM)
    kdeg_pxr_dex <- fixed(0.0015) ; label("PXR-dexamethasone complex degradation (1/min)")  # Table S1/S5 v31 k31 = 0.0015 min-1

    # --- Liganded receptors binding their response elements, and
    #     transcription from the occupied gene. One association constant
    #     (60 nM-1 min-1, diffusion limited) is shared by every gene-binding
    #     reaction in the tables; it is split per target gene here so that a
    #     single regulatory arm can be switched off (the paper's 'stubborn'
    #     network variants) without touching the others. The occupied-gene
    #     dissociation and transcription constants are likewise shared by
    #     the cortisol, L2 and dexamethasone forms of each receptor.
    kon_cypgene <- fixed(60) ; label("Liganded PXR binding to the CYP3A4 gene (1/(nM*min))")  # Table S1/S5 v20 kf20, v33 kf33, V_CYPmRNA_PXR_L2_binding kf = 60 nM-1 min-1
    koff_cypgene <- fixed(1) ; label("Liganded PXR release from the CYP3A4 gene (1/min)")  # Table S1/S5 v20 kb20, v33 kb33, V_CYPmRNA_PXR_L2_binding kb = 1 min-1
    ktx_cyp_act <- fixed(0.05) ; label("CYP3A4 mRNA synthesis from the PXR-occupied gene (1/min)")  # Table S1/S5 v21 k21, v34 k34, V_CYPmRNA_synt_PXR_L2 k = 0.05 min-1
    kon_pxrgene <- fixed(60) ; label("Liganded GR binding to the PXR gene (1/(nM*min))")  # Table S1/S5 v22 kf22, v35 kf35 = 60 nM-1 min-1
    koff_pxrgene <- fixed(200) ; label("Liganded GR release from the PXR gene (1/min)")  # Table S1/S5 v22 kb22, v35 kb35 = 200 min-1
    ktx_pxr_act <- fixed(0.00011) ; label("PXR mRNA synthesis from the GR-occupied gene (1/min)")  # Table S1/S5 v23 k23, v36 k36 = 0.00011 min-1
    kon_grgene <- fixed(60) ; label("Liganded GR binding to the GR gene (1/(nM*min))")  # Table S1/S5 v24 kf24, v37 kf37 = 60 nM-1 min-1
    koff_grgene <- fixed(60) ; label("Liganded GR release from the GR gene (1/min)")  # Table S1/S5 v24 kb24, v37 kb37 = 60 min-1
    ktx_gr_act <- fixed(1.2e-06) ; label("GR mRNA synthesis from the GR-occupied gene (1/min)")  # Table S1/S5 v25 k25, v38 k38 = 1.2e-06 min-1
    kon_tatgene <- fixed(60) ; label("Liganded GR binding to the TAT gene (1/(nM*min))")  # Table S1/S5 v26 kf26, v39 kf39 = 60 nM-1 min-1
    koff_tatgene <- fixed(300) ; label("Liganded GR release from the TAT gene (1/min)")  # Table S1/S5 v26 kb26, v39 kb39 = 300 min-1
    ktx_tat_act <- fixed(0.005) ; label("TAT mRNA synthesis from the GR-occupied gene (1/min)")  # Table S1/S5 v27 k27, v40 k40 = 0.005 min-1

    # --- CYP3A4-mediated metabolism, with cortisol, L2 and dexamethasone as
    #     mutually competitive substrates.
    kcat_cort <- fixed(0.083) ; label("CYP3A4 turnover number for cortisol (1/min)")  # Table S1/S5 v19 kcat = 0.083; BIOMD0000000576 re19 Vm = 0.083
    km_cort <- fixed(15000) ; label("Michaelis constant of CYP3A4 for cortisol (nM)")  # Table S1/S5 v19 KMC = 15000 nM
    km_l2 <- fixed(15000) ; label("Competitive-inhibition constant of L2 on CYP3A4 (nM)")  # BIOMD0000000576 re19 Kms2 and re32 Kms2 = 15000; Table S1 v19 lists the KML term but no value ('Since L2=0, competition neglected')
    kcat_dex <- fixed(0.00425) ; label("CYP3A4 turnover number for dexamethasone (1/min)")  # Table S1/S5 v32 kcat = VmDEX / 1200 nM = 5.1 / 1200 = 0.00425
    km_dex <- fixed(23000) ; label("Michaelis constant of CYP3A4 for dexamethasone (nM)")  # Table S1/S5 v32 KMDEX = 23000 nM

    # --- Cortisol supply, membrane transport and blood binding.
    cortisone <- fixed(24) ; label("Intracellular cortisone, an infinite reservoir for cortisol synthesis (nM)")
    # ^^ Supplementary Table S4/S8 print Cortisone(0) = 23 nM; the SBML
    #    deposit sets the constant species Cortisone = 24 nM. The deposit
    #    value is used: it gives the free blood cortisol of about 15 nM that
    #    Supplementary Figure S3 quotes ('doubling free cortisol
    #    concentration from 15-30nM'). See vignette Errata.
    kf_cort_syn <- fixed(0.016) ; label("Cortisol synthesis from cortisone (1/min)")  # Table S1/S5 v42 kf42 = 0.016 min-1 (fitted to give a realistic cortisol concentration)
    kb_cort_syn <- fixed(0.016) ; label("Cortisol conversion back to cortisone (1/min)")  # Table S1/S5 v42 kb42 = 0.016 min-1
    kin_cort <- fixed(1000) ; label("Free cortisol transport, blood to cell (1/min)")  # Table S1/S5 V_cortisolTransport kctf = 1000 min-1
    kout_cort <- fixed(1000) ; label("Free cortisol transport, cell to blood (1/min)")  # Table S1/S5 V_cortisolTransport kctb = 1000 min-1
    kdist_cort <- fixed(1000) ; label("Distribution of added cortisol into free blood cortisol (1/min)")  # Table S1/S5 V_cort_distribution k = 1000 min-1
    kon_cbg <- fixed(60) ; label("Cortisol-CBG association (1/(nM*min))")  # Table S1/S5 V_cort_CBG kcortCBG_1 = 60 nM-1 min-1
    koff_cbg <- fixed(270) ; label("Cortisol-CBG dissociation (1/min)")  # Table S1/S5 V_cort_CBG kcortCBG_2 = 270 min-1
    kon_alb <- fixed(60) ; label("Cortisol-albumin association (1/(nM*min))")  # Table S1/S5 V_cort_albumin kcortALB_1 = 60 nM-1 min-1
    koff_alb <- fixed(900000) ; label("Cortisol-albumin dissociation (1/min)")  # Table S1/S5 V_cort_albumin kcortALB_2 = 900000 min-1
    vratio_blood <- fixed(5) ; label("Blood volume relative to the hepatocyte volume (unitless)")  # BIOMD0000000576 compartment 'blood' size = 5 (cell compartment 'default' size = 1); not in the supplementary tables
    kin_dex <- fixed(100) ; label("Dexamethasone transport, medium to cell (1/min)")  # Table S1/S5 v41 kf41 = 100 min-1 (high, not to limit uptake)
    kout_dex <- fixed(100) ; label("Dexamethasone transport, cell to medium (1/min)")  # Table S1/S5 v41 kb41 = 100 min-1

    # --- Initial pool of integrated endogenous PXR ligands. 100 nM is the
    #     cortisol + L2 model of Supplementary Tables S5-S8 and Figures 3-5;
    #     set to 0 for the cortisol-only model of Tables S1-S4 and Figure 2.
    bl_ligand2 <- fixed(100) ; label("Initial free L2 ligand pool (nM)")  # Table S8 Ligand2(0) = 100 nM; Results 'This pool was set to a concentration of 100 nM, with a Kd=10 uM'
  })

  model({
    # ===================================================================
    # 1. Rates (nM/min). Each line is the Supplementary Table S1 / S5 rate
    #    law of the named reaction.
    # ===================================================================
    v1 <- ktx_cyp * cyp_gene
    v2 <- kdeg_cyp_mrna * cyp_mrna
    v3 <- ktl_cyp * cyp_mrna
    v4 <- kdeg_cyp * cyp_prot
    v5 <- ktx_pxr * pxr_gene
    v6 <- kdeg_pxr_mrna * pxr_mrna
    v7 <- ktl_pxr * pxr_mrna
    v8 <- kdeg_pxr * pxr_prot
    v9 <- ktx_gr * gr_gene
    v10 <- kdeg_gr_mrna * gr_mrna
    v11 <- ktl_gr * gr_mrna
    v12 <- kdeg_gr * gr_prot
    v13 <- ktx_tat * tat_gene
    v14 <- kdeg_tat_mrna * tat_mrna
    v15 <- kon_gr_cort * cort * gr_prot - koff_gr_cort * gr_cort
    v16 <- kdeg_gr_cort * gr_cort
    v17 <- kon_pxr_cort * pxr_prot * cort - koff_pxr_cort * pxr_cort
    v18 <- kdeg_pxr_cort * pxr_cort

    # Competitive Michaelis-Menten metabolism by CYP3A4 (v19 cortisol, v32
    # dexamethasone). The denominators follow the deposit's LigandDegrOld
    # function: all three CYP3A4 substrates (cortisol, L2, dexamethasone)
    # compete. L2 itself is not metabolised.
    v19 <- kcat_cort * cyp_prot * (cort / km_cort) /
      (1 + cort / km_cort + ligand2 / km_l2 + dex / km_dex)
    v32 <- kcat_dex * cyp_prot * (dex / km_dex) /
      (1 + dex / km_dex + ligand2 / km_l2 + cort / km_cort)

    v20 <- kon_cypgene * cyp_gene * pxr_cort - koff_cypgene * cypgene_pxrcort
    v21 <- ktx_cyp_act * cypgene_pxrcort
    v22 <- kon_pxrgene * pxr_gene * gr_cort - koff_pxrgene * pxrgene_grcort
    v23 <- ktx_pxr_act * pxrgene_grcort
    v24 <- kon_grgene * gr_gene * gr_cort - koff_grgene * grgene_grcort
    v25 <- ktx_gr_act * grgene_grcort
    v26 <- kon_tatgene * tat_gene * gr_cort - koff_tatgene * tatgene_grcort
    v27 <- ktx_tat_act * tatgene_grcort
    v28 <- kon_gr_dex * gr_prot * dex - koff_gr_dex * gr_dex
    v29 <- kdeg_gr_dex * gr_dex
    v30 <- kon_pxr_dex * pxr_prot * dex - koff_pxr_dex * pxr_dex
    v31 <- kdeg_pxr_dex * pxr_dex
    v33 <- kon_cypgene * cyp_gene * pxr_dex - koff_cypgene * cypgene_pxrdex
    v34 <- ktx_cyp_act * cypgene_pxrdex
    v35 <- kon_pxrgene * pxr_gene * gr_dex - koff_pxrgene * pxrgene_grdex
    v36 <- ktx_pxr_act * pxrgene_grdex
    v37 <- kon_grgene * gr_gene * gr_dex - koff_grgene * grgene_grdex
    v38 <- ktx_gr_act * grgene_grdex
    v39 <- kon_tatgene * tat_gene * gr_dex - koff_tatgene * tatgene_grdex
    v40 <- ktx_tat_act * tatgene_grdex
    v41 <- kin_dex * dex_out - kout_dex * dex
    v42 <- kf_cort_syn * cortisone - kb_cort_syn * cort
    v43 <- ktl_tat * tat_mrna
    v44 <- kdeg_tat * tat_prot

    v_l2_bind <- kon_pxr_l2 * pxr_prot * ligand2 - koff_pxr_l2 * pxr_l2
    v_l2_deg <- kdeg_pxr_l2 * pxr_l2
    v_l2_cypgene <- kon_cypgene * cyp_gene * pxr_l2 - koff_cypgene * cypgene_pxrl2
    v_l2_cyp_mrna <- ktx_cyp_act * cypgene_pxrl2

    # Blood side (Supplementary Table S1 V_cortisolTransport,
    # V_cort_distribution, V_cort_CBG, V_cort_albumin). The transport flux
    # is per unit cell volume, so the blood concentration changes by
    # 1/vratio_blood of it; the binding and distribution rates are already
    # per unit blood volume (deposit reactions cortisolTransport,
    # Cortisol_CBG, Cort_Alb, cort_distribution).
    v_transport <- kin_cort * cort_blood - kout_cort * cort
    v_dist <- kdist_cort * cort_added
    v_cbg <- kon_cbg * cort_blood * cbg - koff_cbg * cbg_cort
    v_alb <- kon_alb * alb * cort_blood - koff_alb * alb_cort

    # ===================================================================
    # 2. Balance equations, Supplementary Table S6, with the corrections
    #    listed in the vignette Errata (all four agree with the deposit):
    #    d/dt(gr_prot) carries v12 once, not twice; CYP3A4 mRNA is made by
    #    the L2 transcription flux, not the L2 gene-binding flux; the
    #    L2-occupied CYP3A4 gene gains the binding flux only; the free
    #    CYP3A4 gene loses the L2 binding flux.
    # ===================================================================
    d/dt(cort) <- -v15 + v16 - v17 + v18 - v19 + v42 + v_transport
    d/dt(cort_added) <- -v_dist
    d/dt(cort_blood) <- v_dist - v_cbg - v_alb - v_transport / vratio_blood
    d/dt(cbg) <- -v_cbg
    d/dt(cbg_cort) <- v_cbg
    d/dt(alb) <- -v_alb
    d/dt(alb_cort) <- v_alb
    d/dt(ligand2) <- -v_l2_bind + v_l2_deg
    d/dt(dex) <- -v28 + v29 - v30 + v31 - v32 + v41
    # Culture-medium dexamethasone is a constant boundary species in the
    # deposit (DEXout, constant = TRUE): a dose sets the medium level and it
    # is not depleted by uptake.
    d/dt(dex_out) <- 0

    d/dt(gr_gene) <- -v24 - v37
    d/dt(gr_mrna) <- v9 - v10 + v25 + v38
    d/dt(gr_prot) <- v11 - v12 - v15 - v28
    d/dt(gr_cort) <- v15 - v16 - v22 - v24 - v26
    d/dt(gr_dex) <- v28 - v29 - v35 - v37 - v39
    d/dt(grgene_grcort) <- v24
    d/dt(grgene_grdex) <- v37

    d/dt(pxr_gene) <- -v22 - v35
    d/dt(pxr_mrna) <- v5 - v6 + v23 + v36
    d/dt(pxr_prot) <- v7 - v8 - v17 - v30 - v_l2_bind
    d/dt(pxr_cort) <- v17 - v18 - v20
    d/dt(pxr_l2) <- v_l2_bind - v_l2_deg - v_l2_cypgene
    d/dt(pxr_dex) <- v30 - v31 - v33
    d/dt(pxrgene_grcort) <- v22
    d/dt(pxrgene_grdex) <- v35

    d/dt(cyp_gene) <- -v20 - v33 - v_l2_cypgene
    d/dt(cyp_mrna) <- v1 - v2 + v21 + v34 + v_l2_cyp_mrna
    d/dt(cyp_prot) <- v3 - v4
    d/dt(cypgene_pxrcort) <- v20
    d/dt(cypgene_pxrl2) <- v_l2_cypgene
    d/dt(cypgene_pxrdex) <- v33

    d/dt(tat_gene) <- -v26 - v39
    d/dt(tat_mrna) <- v13 - v14 + v27 + v40
    d/dt(tat_prot) <- v43 - v44
    d/dt(tatgene_grcort) <- v26
    d/dt(tatgene_grdex) <- v39

    # ===================================================================
    # 3. Initial conditions, Supplementary Table S8 (states not listed
    #    start at rxode2's default of 0). These are the authors' starting
    #    point, not a steady state; the published simulations first run the
    #    network to steady state (see the vignette).
    # ===================================================================
    cort(0) <- 1.14
    cbg(0) <- 550
    alb(0) <- 60000
    ligand2(0) <- bl_ligand2
    cyp_gene(0) <- 0.82976
    cyp_prot(0) <- 104
    cyp_mrna(0) <- 0.075
    cypgene_pxrcort(0) <- 0.00024
    gr_gene(0) <- 0.5
    gr_prot(0) <- 47.24
    gr_mrna(0) <- 0.0008
    grgene_grcort(0) <- 0.33
    gr_cort(0) <- 32.36
    pxr_gene(0) <- 0.78
    pxr_prot(0) <- 99.9
    pxr_mrna(0) <- 0.007
    pxrgene_grcort(0) <- 0.05
    pxr_cort(0) <- 0.01
    tat_gene(0) <- 0.81
    tat_mrna(0) <- 0.4
    tat_prot(0) <- 0.462028
    tatgene_grcort(0) <- 0.02

    # ===================================================================
    # 4. Reported outputs. Total receptor protein (free plus every liganded
    #    and gene-bound form) is what the figures plot as 'GR' and 'PXR'
    #    (deposit assignment rules GRprotein and PXRprotein).
    # ===================================================================
    gr_total <- gr_prot + gr_cort + gr_dex + grgene_grcort + grgene_grdex +
      pxrgene_grcort + pxrgene_grdex + tatgene_grcort + tatgene_grdex
    pxr_total <- pxr_prot + pxr_cort + pxr_l2 + pxr_dex + cypgene_pxrcort +
      cypgene_pxrl2 + cypgene_pxrdex
    pxr_active <- pxr_cort + pxr_l2
  })
}
