Desai_2024_crisprCas9_nhp <- function() {
  description <- paste(
    "Preclinical (cynomolgus monkey). QSP. Translational quantitative",
    "systems pharmacology platform for in vivo CRISPR-Cas9 gene editing",
    "delivered by lipid nanoparticle (LNP), covering the whole-body",
    "disposition of all three components of the therapy - the LNP",
    "delivery vehicle, the single guide RNA (sgRNA) and the Cas9",
    "messenger RNA (mRNA) - plus intracellular Cas9 translation and",
    "ribonucleoprotein (RNP) assembly (Desai 2024). The PK layer is a",
    "minimal-PBPK structure with plasma, lymph, a three-region liver",
    "(vascular / endosomal / interstitial-cellular), kidney and a",
    "lumped remainder tissue; the kidney compartment is present for LNP",
    "and sgRNA but absent for mRNA, which is too large to be renally",
    "cleared. LNP additionally undergoes opsonin bio-corona formation",
    "(kass / kdis), phagocytosis into the mononuclear phagocyte system",
    "(kint) and LDL-receptor-mediated endocytosis (kon,LNP / koff,LNP)",
    "with receptor recycling. Two pharmacodynamic case studies are",
    "attached to the liver RNP concentration: NTLA-2001 (transthyretin",
    "amyloidosis) via a sigmoidal indirect-response model inhibiting",
    "TTR production, and VERVE-101 (LDL-cholesterol lowering) via a",
    "Friberg-type transit / feedback model on serum PCSK9 feeding a",
    "precursor-dependent LDL-cholesterol model. Estimated against",
    "cynomolgus monkey LNP plasma PK (1-3 mg/kg total RNA), serum TTR",
    "(1.5-6 mg/kg) and serum PCSK9 / LDL-cholesterol (0.75-1.5 mg/kg),",
    "all by 2 h IV infusion. Sibling extractions:",
    "Desai_2024_crisprCas9_mouse and Desai_2024_crisprCas9_human.",
    sep = " "
  )
  reference <- paste(
    "Desai DA, Schmidt S, Cristofoletti R (2024).",
    "A quantitative systems pharmacology (QSP) platform for preclinical",
    "to clinical translation of in-vivo CRISPR-Cas therapy.",
    "Frontiers in Pharmacology 15:1454785.",
    "doi:10.3389/fphar.2024.1454785. PMCID: PMC11449743.",
    "Model equations from Supplementary Information S1; physiological",
    "volumes and flows from Supplementary Tables S1 and S2;",
    "drug-specific parameters from Table 2.",
    sep = " "
  )
  vignette <- "Desai_2024_crisprCas9"

  # Every state below is a paper-mechanistic QSP species x physiological
  # sub-compartment. Naming scheme is <molecule>_<region> for the three
  # disposition species, plus the intracellular / receptor / PD states.
  paper_specific_compartments <- c(
    "lnp_plasma", "lnp_lymph", "lnp_livervas", "lnp_mps", "lnp_opsonin",
    "lnp_liverendo", "lnp_ldlr", "lnp_liverinter", "lnp_kidney", "lnp_rem",
    "ldlr",
    "sgrna_plasma", "sgrna_lymph", "sgrna_livervas", "sgrna_liverendo",
    "sgrna_liverinter", "sgrna_kidney", "sgrna_rem", "sgrna_cell",
    "mrna_plasma", "mrna_lymph", "mrna_livervas", "mrna_liverendo",
    "mrna_liverinter", "mrna_rem",
    "cas9", "rnp",
    "ttr", "pcsk9_prol", "pcsk9_transit1", "pcsk9_circ", "ldlc"
  )

  units <- list(
    time = "h",
    dosing = paste(
      "ug. Three simultaneous IV inputs are required: the LNP dose into",
      "lnp_plasma, the sgRNA dose into sgrna_plasma and the Cas9 mRNA",
      "dose into mrna_plasma. Desai 2024 assumed immediate release of",
      "the transgene product after administration (Model assumptions,",
      "point 2), so the RNA cargo is dosed directly into plasma",
      "alongside the vehicle. Total RNA splits 33.3% sgRNA / 66.7%",
      "mRNA and the paired LNP dose is given in Table 1 (e.g. 1 mg/kg",
      "total RNA in a 5 kg NHP = 1665 ug sgRNA + 3335 ug mRNA with",
      "18.5 mg/kg = 92500 ug LNP).",
      sep = " "
    ),
    concentration = paste(
      "ug/mL for all LNP, sgRNA, mRNA, Cas9 and RNP concentrations.",
      "The three PD outputs (TTR, PCSK9, LDLC) are percent of baseline.",
      sep = " "
    )
  )

  covariateData <- list()

  compartmentData <- list(
    lnp_plasma = list(analyte = "lipid nanoparticle", units = "ug", specimen = "plasma", verified = TRUE),
    lnp_lymph = list(analyte = "lipid nanoparticle", units = "ug", specimen = "lymph", verified = TRUE),
    lnp_livervas = list(analyte = "lipid nanoparticle", units = "ug", specimen = "plasma", verified = TRUE),
    lnp_mps = list(analyte = "lipid nanoparticle", units = "ug", specimen = "tissue", verified = TRUE),
    lnp_opsonin = list(analyte = "opsonin-bound lipid nanoparticle", units = "ug", specimen = "plasma", verified = TRUE),
    lnp_liverendo = list(analyte = "lipid nanoparticle", units = "ug", specimen = "endosome", verified = TRUE),
    lnp_ldlr = list(analyte = "LNP-LDL receptor complex", units = "ug", specimen = "endosome", verified = TRUE),
    lnp_liverinter = list(analyte = "lipid nanoparticle", units = "ug", specimen = "tissue", verified = TRUE),
    lnp_kidney = list(analyte = "lipid nanoparticle", units = "ug", specimen = "tissue", verified = TRUE),
    lnp_rem = list(analyte = "lipid nanoparticle", units = "ug", specimen = "tissue", verified = TRUE),
    ldlr = list(analyte = "LDL receptor", units = "ug", specimen = "endosome", verified = TRUE),
    sgrna_plasma = list(analyte = "single guide RNA", units = "ug", specimen = "plasma", verified = TRUE),
    sgrna_lymph = list(analyte = "single guide RNA", units = "ug", specimen = "lymph", verified = TRUE),
    sgrna_livervas = list(analyte = "single guide RNA", units = "ug", specimen = "plasma", verified = TRUE),
    sgrna_liverendo = list(analyte = "single guide RNA", units = "ug", specimen = "endosome", verified = TRUE),
    sgrna_liverinter = list(analyte = "single guide RNA", units = "ug", specimen = "tissue", verified = TRUE),
    sgrna_kidney = list(analyte = "single guide RNA", units = "ug", specimen = "tissue", verified = TRUE),
    sgrna_rem = list(analyte = "single guide RNA", units = "ug", specimen = "tissue", verified = TRUE),
    sgrna_cell = list(analyte = "single guide RNA", units = "ug", specimen = "tissue", verified = TRUE),
    mrna_plasma = list(analyte = "Cas9 messenger RNA", units = "ug", specimen = "plasma", verified = TRUE),
    mrna_lymph = list(analyte = "Cas9 messenger RNA", units = "ug", specimen = "lymph", verified = TRUE),
    mrna_livervas = list(analyte = "Cas9 messenger RNA", units = "ug", specimen = "plasma", verified = TRUE),
    mrna_liverendo = list(analyte = "Cas9 messenger RNA", units = "ug", specimen = "endosome", verified = TRUE),
    mrna_liverinter = list(analyte = "Cas9 messenger RNA", units = "ug", specimen = "tissue", verified = TRUE),
    mrna_rem = list(analyte = "Cas9 messenger RNA", units = "ug", specimen = "tissue", verified = TRUE),
    cas9 = list(analyte = "Cas9 protein", units = "ug", specimen = "tissue", verified = TRUE),
    rnp = list(analyte = "Cas9-sgRNA ribonucleoprotein complex", units = "ug", specimen = "tissue", verified = TRUE),
    ttr = list(analyte = "transthyretin", units = "% of baseline", specimen = "serum", verified = TRUE),
    pcsk9_prol = list(analyte = "PCSK9 proliferating precursor", units = "% of baseline", specimen = "not applicable", verified = TRUE),
    pcsk9_transit1 = list(analyte = "PCSK9 transit pool", units = "% of baseline", specimen = "not applicable", verified = TRUE),
    pcsk9_circ = list(analyte = "PCSK9", units = "% of baseline", specimen = "serum", verified = TRUE),
    ldlc = list(analyte = "LDL cholesterol", units = "% of baseline", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "cynomolgus monkey (non-human primate)",
    n_subjects = NA_integer_,
    n_studies = 2L,
    weight_median = "5 kg",
    disease_state = paste(
      "Healthy cynomolgus monkeys used as the preclinical model for",
      "transthyretin amyloidosis (NTLA-2001) and for LDL-cholesterol",
      "lowering (VERVE-101).",
      sep = " "
    ),
    dose_range = paste(
      "LNP plasma PK: 1, 2 and 3 mg/kg total RNA (18.5, 36.7 and 55.5",
      "mg/kg LNP). Serum TTR: 1.5, 3 and 6 mg/kg total RNA (27.75,",
      "68.82 and 137.64 mg/kg LNP). Serum PCSK9 / LDL-cholesterol:",
      "0.75 and 1.5 mg/kg total RNA (17.2 and 27.75 mg/kg LNP).",
      "All by 2 h IV infusion.",
      sep = " "
    ),
    notes = paste(
      "Desai 2024 Table 1. No individual-level cohort is described: the",
      "analysis was performed on mean profiles digitised from Gillmore",
      "2021 (NTLA-2001 LNP PK and serum TTR) and Lee 2023 (VERVE-101",
      "serum PCSK9 and LDL cholesterol) with WebPlotDigitizer, so",
      "n_subjects is not reported. The reference body weight for the",
      "NHP physiology in Supplementary Tables S1 and S2 is 5 kg.",
      sep = " "
    )
  )

  ini({
    # ---------------------------------------------------------------
    # Estimated drug-specific attributes -- Table 2, NHP column
    # ---------------------------------------------------------------
    lkin_endo <- log(0.039); label("Rate of endocytosis for LNP, sgRNA and mRNA (1/h)")   # Table 2, NHP (43% RSE)
    lkout_exo <- log(2690); label("Rate of exocytosis for LNP, sgRNA and mRNA (1/h)")     # Table 2, NHP (65.4% RSE)
    lkdeg_dr <- log(2.04); label("Rate of degradation of the LNP-LDL receptor complex (1/h)")  # Table 2, NHP (12% RSE)
    lldltot <- log(539); label("Total LDL receptor concentration (ug/mL)")                # Table 2, NHP (66.2% RSE)
    lkdis <- log(8.64); label("Rate of dissociation from opsonins (1/h)")                 # Table 2, NHP (39% RSE)
    lkass <- log(5.06); label("Rate of association to opsonins (1/h)")                    # Table 2, NHP (14.5% RSE)
    lkrelease <- log(0.634); label("Rate of release of sgRNA and mRNA from LNP (1/h)")    # Table 2, NHP (5.31% RSE)

    # ---------------------------------------------------------------
    # Drug-specific attributes held constant -- Table 2, NHP column
    # ---------------------------------------------------------------
    kdeg_lnp <- fixed(1.6486); label("Rate of degradation of unbound LNP (1/h); literature value from Miyazawa 2024")  # Table 2, NHP
    lkel <- fixed(log(0.231)); label("Rate of elimination of the LDL receptor (1/h); literature value from Harwood 1997")    # Table 2, NHP
    kdeg_sgrna <- fixed(2.01); label("Rate of degradation of sgRNA (1/h); taken from the step-1 mechanistic model")    # Table 2, NHP footnote c
    kdeg_mrna <- fixed(0.1232); label("Rate of degradation of mRNA (1/h); literature value from Miyazawa 2024")        # Table 2, NHP
    koff_rnp <- fixed(0.00188); label("Rate of dissociation of the RNP complex (1/h); taken from the step-1 mechanistic model")  # Table 2, NHP footnote c
    kd <- fixed(0.49); label("Equilibrium dissociation constant for Cas9-sgRNA (nM); literature value from Sternberg 2014")      # Table 2
    kint <- fixed(0.9063); label("Rate of phagocytosis into the mononuclear phagocyte system (1/h); literature value from Miyazawa 2024")  # Table 2
    kon_lnp <- fixed(0.18); label("Rate of association of LNP with the LDL receptor (mL/ug/h); literature value from Harwood 1997")       # Table 2
    koff_lnp <- fixed(33.12); label("Rate of dissociation of the LNP-LDL receptor complex (1/h); literature value from Harwood 1997")     # Table 2
    fu_lnp <- fixed(0.002); label("Plasma free fraction of LNP (unitless); literature value from Mager 2012")   # Table 2
    fu_rna <- fixed(0.15); label("Plasma free fraction of sgRNA and mRNA (unitless); literature value from Ayyar 2021")  # Table 2
    ktrans <- fixed(0.36); label("Rate of translation from mRNA to Cas9 (1/h); taken from the step-1 mechanistic model")  # Table 2, NHP footnote c

    # ---------------------------------------------------------------
    # TTR pharmacodynamics (NTLA-2001) -- Table 2, NHP column
    # ---------------------------------------------------------------
    rbase_ttr <- fixed(100); label("Baseline serum TTR (% of baseline)")             # Table 2, NHP
    lkout_ttr <- log(0.493); label("First-order degradation rate of TTR protein (1/day)")  # Table 2, NHP (15% RSE)
    limax_ttr <- log(0.961); label("Maximum inhibition of TTR production (unitless)")      # Table 2, NHP (0.625% RSE)
    lic50_ttr <- log(4.77); label("Liver RNP concentration at 50% inhibition of TTR production (ug/mL)")  # Table 2, NHP (11.5% RSE)
    lhill_ttr <- log(0.31); label("Hill coefficient of the TTR inhibition function (unitless)")           # Table 2, NHP gamma (13.3% RSE)

    # ---------------------------------------------------------------
    # PCSK9 and LDL-cholesterol pharmacodynamics (VERVE-101)
    # -- Table 2, NHP column
    # ---------------------------------------------------------------
    rbase_pcsk9 <- fixed(100); label("Baseline serum PCSK9 (% of baseline)")   # Table 2, NHP
    lmtt <- log(14.5); label("Mean transit time for PCSK9 (day)")              # Table 2, NHP (0.563% RSE)
    limax_pcsk9 <- log(0.771); label("Maximum inhibition of PCSK9 precursor proliferation (unitless)")  # Table 2, NHP (14.2% RSE)
    lic50_pcsk9 <- log(21.5); label("Liver RNP concentration at 50% inhibition of PCSK9 (ug/mL)")       # Table 2, NHP (0.994% RSE)
    gamma_pcsk9 <- fixed(1.1); label("Feedback exponent on circulating PCSK9 (unitless)")               # Table 2, NHP
    rbase_ldlc <- fixed(100); label("Baseline serum LDL cholesterol (% of baseline)")                   # Suppl Info S1, LDL equation initial condition
    lgamma_ldlc <- log(0.672); label("Precursor exponent linking PCSK9 to LDL-cholesterol synthesis (unitless)")  # Table 2, NHP (18.5% RSE)
    lkdeg_ldlc <- log(4.66); label("Rate of degradation of LDL cholesterol (1/day)")                    # Table 2 continued, NHP (67.4% RSE)

    # ---------------------------------------------------------------
    # Residual error -- proportional per output (Methods, Software).
    # Desai 2024 states a proportional error model was used for every
    # output but does not report any sigma_slope value.
    # ---------------------------------------------------------------
    propSd <- fixed(0); label("Proportional residual SD for LNP plasma concentration (fraction); not published")     # Methods, Software
    propSd_Cc_sgrna <- fixed(0); label("Proportional residual SD for sgRNA plasma concentration (fraction); not published")  # Methods, Software
    propSd_Cc_mrna <- fixed(0); label("Proportional residual SD for mRNA plasma concentration (fraction); not published")    # Methods, Software
    propSd_TTR <- fixed(0); label("Proportional residual SD for serum TTR (fraction); not published")                # Methods, Software
    propSd_PCSK9 <- fixed(0); label("Proportional residual SD for serum PCSK9 (fraction); not published")            # Methods, Software
    propSd_LDLC <- fixed(0); label("Proportional residual SD for serum LDL cholesterol (fraction); not published")   # Methods, Software
  })

  model({
    # ===============================================================
    # 0. Species physiology -- held as structural constants of this
    #    species-specific file rather than as estimable parameters.
    #    Cynomolgus monkey, 5 kg reference body weight.
    # ===============================================================
    # Physiological volumes (mL) -- Supplementary Table S1, NHP column
    vplasma <- 187      # Plasma volume in circulation (mL)
    vlymph <- 25.1      # Volume of lymph (mL)
    vlivervas <- 15.9   # Plasma volume in liver vasculature (mL)
    vliverendo <- 0.934 # Endosomal volume of liver (mL)
    vliverinter <- 157.4 # Cellular and interstitial volume of liver (mL)
    vmps <- 12          # Volume of the mononuclear phagocyte system (mL)
    vkidney <- 25.9     # Plasma, interstitial and cellular volume of kidney (mL)
    vrem <- 1087        # Extracellular volume of remainder tissue (mL)

    # Physiological flows (mL/h) -- Supplementary Table S2, NHP column
    ql <- 1251          # Plasma flow to liver (mL/h)
    ll <- 2.502         # Lymph flow to liver (mL/h)
    qt <- 10358         # Plasma flow to remainder tissue (mL/h)
    lt <- 20.716        # Lymph flow to remainder tissue (mL/h)
    qk <- 3237          # Plasma flow to kidney (mL/h)
    lk <- 6.464         # Lymph flow to kidney (mL/h)
    gfr <- 1138         # Glomerular filtration rate (mL/h)

    # Reflection coefficients (unitless) -- Supplementary Table S2.
    # sigma_v gates the transcapillary convective flux out of the
    # vascular space and is the row Desai 2024 labels "Reflection
    # coefficient (pinocytosis)"; the Methods define sigma_V as the
    # resistance of the vascular endothelium on the paracellular
    # macropinocytosis route. sigma_l gates the interstitium-to-lymph
    # convective return and is the row labelled "Reflection
    # coefficient"; 0.2 is the lymphatic reflection coefficient of the
    # Shah & Betts 2012 platform these values are fixed from.
    sigma_v <- 0.9
    sigma_l <- 0.2

    # ===============================================================
    # 1. Individual parameters
    # ===============================================================
    kin_endo <- exp(lkin_endo)
    kout_exo <- exp(lkout_exo)
    kdeg_dr <- exp(lkdeg_dr)
    ldltot <- exp(lldltot)
    kdis <- exp(lkdis)
    kass <- exp(lkass)
    krelease <- exp(lkrelease)
    kel <- exp(lkel)

    # ===============================================================
    # 2. Derived drug-specific attributes -- Suppl Info S1,
    #    "Drug Specific attributes"
    # ===============================================================
    # CL_in / CL_out are rate constants multiplied by the physiological
    # volume of the region they act on. Every CL_out in the published
    # ODEs (CL_out,DR and CL_out,liverendo) is defined as
    # k_out,exo * V_liverendo.
    clin_livervas <- kin_endo * vlivervas
    clin_liverinter <- kin_endo * vliverinter
    clin_dr <- kin_endo * vliverendo
    clout_endo <- kout_exo * vliverendo

    # Renal clearance of the free fraction, CL_R = fu * GFR.
    clr_lnp <- fu_lnp * gfr
    clr_rna <- fu_rna * gfr

    # LDL-receptor synthesis balances elimination at baseline,
    # k_syn = LDL_tot * k_el.
    ksyn_ldlr <- ldltot * kel

    # Cas9-sgRNA association rate, kon = koff / KD (Table 2). KD is
    # tabulated in nM while every concentration in the model is in
    # ug/mL; the ratio is evaluated numerically as printed.
    kon_rnp <- koff_rnp / kd

    # ===============================================================
    # 3. Concentrations (ug/mL) from the amount states
    # ===============================================================
    c_lnp_plasma <- lnp_plasma / vplasma
    c_lnp_lymph <- lnp_lymph / vlymph
    c_lnp_livervas <- lnp_livervas / vlivervas
    c_lnp_mps <- lnp_mps / vmps
    c_lnp_liverendo <- lnp_liverendo / vliverendo
    c_lnp_ldlr <- lnp_ldlr / vliverendo
    c_lnp_liverinter <- lnp_liverinter / vliverinter
    c_lnp_kidney <- lnp_kidney / vkidney
    c_lnp_rem <- lnp_rem / vrem
    c_ldlr <- ldlr / vliverendo

    c_sgrna_plasma <- sgrna_plasma / vplasma
    c_sgrna_lymph <- sgrna_lymph / vlymph
    c_sgrna_livervas <- sgrna_livervas / vlivervas
    c_sgrna_liverendo <- sgrna_liverendo / vliverendo
    c_sgrna_liverinter <- sgrna_liverinter / vliverinter
    c_sgrna_kidney <- sgrna_kidney / vkidney
    c_sgrna_rem <- sgrna_rem / vrem
    c_sgrna_cell <- sgrna_cell / vliverinter

    c_mrna_plasma <- mrna_plasma / vplasma
    c_mrna_lymph <- mrna_lymph / vlymph
    c_mrna_livervas <- mrna_livervas / vlivervas
    c_mrna_liverendo <- mrna_liverendo / vliverendo
    c_mrna_liverinter <- mrna_liverinter / vliverinter
    c_mrna_rem <- mrna_rem / vrem

    c_cas9 <- cas9 / vliverinter
    c_rnp <- rnp / vliverinter

    # ===============================================================
    # 4. Lipid nanoparticle disposition -- Suppl Info S1, "Lipid
    #    Nanoparticle (LNP)". States are amounts; each published
    #    dC/dt equation is the amount equation divided by the region
    #    volume, so the amount ODE is the printed numerator.
    # ===============================================================
    d/dt(lnp_plasma) <-
      (ql - ll) * c_lnp_livervas + (qt - lt) * c_lnp_rem +
      (qk - lk) * c_lnp_kidney + (ll + lt + lk) * c_lnp_lymph -
      ql * c_lnp_plasma - qt * c_lnp_plasma - qk * c_lnp_plasma

    d/dt(lnp_lymph) <-
      ll * (1 - sigma_l) * c_lnp_liverinter +
      lt * (1 - sigma_l) * c_lnp_rem +
      lk * (1 - sigma_l) * c_lnp_kidney -
      (ll + lk + lt) * c_lnp_lymph

    # The opsonin (bio-corona) compartment is carried as an amount.
    # Desai 2024 writes it as a concentration divided by V_opsonins,
    # a volume that is not reported in Supplementary Table S1; because
    # the opsonin term enters the liver-vascular equation as
    # k_dis * V_opsonins * C_Opsonins, V_opsonins cancels exactly and
    # the amount form below is identical to the printed system.
    d/dt(lnp_livervas) <-
      ql * c_lnp_plasma - (ql - ll) * c_lnp_livervas -
      clin_livervas * c_lnp_livervas -
      ll * (1 - sigma_v) * c_lnp_livervas -
      kint * vlivervas * c_lnp_livervas -
      kass * vlivervas * c_lnp_livervas +
      kdis * lnp_opsonin +
      clout_endo * c_lnp_ldlr

    d/dt(lnp_mps) <-
      kint * vlivervas * c_lnp_livervas - kdeg_lnp * vmps * c_lnp_mps

    d/dt(lnp_opsonin) <-
      kass * vlivervas * c_lnp_livervas - kdis * lnp_opsonin

    d/dt(lnp_liverendo) <-
      clin_livervas * c_lnp_livervas -
      kdeg_lnp * vliverendo * c_lnp_liverendo +
      koff_lnp * vliverendo * c_lnp_ldlr -
      kon_lnp * c_lnp_liverendo * c_ldlr * vliverendo +
      clin_liverinter * c_lnp_liverinter

    d/dt(ldlr) <-
      koff_lnp * vliverendo * c_lnp_ldlr -
      kon_lnp * c_lnp_liverendo * c_ldlr * vliverendo -
      kel * vliverendo * c_ldlr +
      ksyn_ldlr * vliverendo +
      clin_dr * c_lnp_ldlr

    d/dt(lnp_ldlr) <-
      kon_lnp * c_lnp_liverendo * c_ldlr * vliverendo -
      koff_lnp * vliverendo * c_lnp_ldlr -
      clin_dr * c_lnp_ldlr -
      kdeg_dr * vliverendo * c_lnp_ldlr

    # The k_release term is subtracted TWICE from the interstitial LNP
    # pool: Desai 2024 lumps the release rate of sgRNA and of Cas9 mRNA
    # into a single k_release, so one LNP-release flux feeds the sgRNA
    # interstitial pool and a second feeds the mRNA interstitial pool.
    d/dt(lnp_liverinter) <-
      clout_endo * c_lnp_ldlr -
      clin_liverinter * c_lnp_liverinter +
      ll * (1 - sigma_v) * c_lnp_livervas -
      krelease * vliverinter * c_lnp_liverinter -
      krelease * vliverinter * c_lnp_liverinter -
      ll * (1 - sigma_l) * c_lnp_liverinter

    d/dt(lnp_kidney) <-
      qk * c_lnp_plasma - (qk - lk) * c_lnp_kidney -
      clr_lnp * c_lnp_kidney - lk * (1 - sigma_l) * c_lnp_kidney

    d/dt(lnp_rem) <-
      qt * c_lnp_plasma - (qt - lt) * c_lnp_rem -
      kdeg_lnp * vrem * c_lnp_rem - lt * (1 - sigma_l) * c_lnp_rem

    # ===============================================================
    # 5. Single guide RNA disposition -- Suppl Info S1, "Single guide
    #    RNA (sgRNA)"
    # ===============================================================
    d/dt(sgrna_plasma) <-
      (ql - ll) * c_sgrna_livervas + (qt - lt) * c_sgrna_rem +
      (qk - lk) * c_sgrna_kidney + (ll + lt + lk) * c_sgrna_lymph -
      ql * c_sgrna_plasma - qt * c_sgrna_plasma - qk * c_sgrna_plasma

    d/dt(sgrna_lymph) <-
      ll * (1 - sigma_l) * c_sgrna_liverinter +
      lt * (1 - sigma_l) * c_sgrna_rem -
      (ll + lt + lk) * c_sgrna_lymph +
      lk * (1 - sigma_l) * c_sgrna_kidney

    d/dt(sgrna_livervas) <-
      ql * c_sgrna_plasma - (ql - ll) * c_sgrna_livervas -
      clin_livervas * c_sgrna_livervas +
      clout_endo * c_sgrna_liverendo -
      ll * (1 - sigma_v) * c_sgrna_livervas

    # Two exocytosis fluxes leave the endosomal pool, one returning to
    # the liver vasculature and one delivering to the interstitium;
    # both are printed in Suppl Info S1 and both are matched by gains
    # in the receiving equations.
    d/dt(sgrna_liverendo) <-
      clin_livervas * c_sgrna_livervas -
      clout_endo * c_sgrna_liverendo -
      kdeg_sgrna * vliverendo * c_sgrna_liverendo +
      clin_liverinter * c_sgrna_liverinter -
      clout_endo * c_sgrna_liverendo

    d/dt(sgrna_liverinter) <-
      clout_endo * c_sgrna_liverendo -
      clin_liverinter * c_sgrna_liverinter +
      ll * (1 - sigma_v) * c_sgrna_livervas +
      krelease * vliverinter * c_lnp_liverinter -
      ll * (1 - sigma_l) * c_sgrna_liverinter -
      krelease * vliverinter * c_sgrna_liverinter

    d/dt(sgrna_kidney) <-
      qk * c_sgrna_plasma - (qk - lk) * c_sgrna_kidney -
      clr_rna * c_sgrna_kidney - lk * (1 - sigma_l) * c_sgrna_kidney

    d/dt(sgrna_rem) <-
      qt * c_sgrna_plasma - (qt - lt) * c_sgrna_rem -
      kdeg_sgrna * vrem * c_sgrna_rem -
      lt * (1 - sigma_l) * c_sgrna_rem

    # ===============================================================
    # 6. Cas9 messenger RNA disposition -- Suppl Info S1, "Messenger
    #    RNA (mRNA) and Cas protein". No kidney compartment: mRNA is
    #    too large to undergo renal elimination (Methods, Model
    #    structure and workflow).
    # ===============================================================
    d/dt(mrna_plasma) <-
      (ql - ll) * c_mrna_livervas + (qt - lt) * c_mrna_rem +
      (ll + lt) * c_mrna_lymph -
      ql * c_mrna_plasma - qt * c_mrna_plasma

    d/dt(mrna_lymph) <-
      ll * (1 - sigma_l) * c_mrna_liverinter +
      lt * (1 - sigma_l) * c_mrna_rem -
      (ll + lt) * c_mrna_lymph

    d/dt(mrna_livervas) <-
      ql * c_mrna_plasma - (ql - ll) * c_mrna_livervas -
      clin_livervas * c_mrna_livervas +
      clout_endo * c_mrna_liverendo -
      ll * (1 - sigma_v) * c_mrna_livervas

    d/dt(mrna_liverendo) <-
      clin_livervas * c_mrna_livervas -
      clout_endo * c_mrna_liverendo -
      kdeg_mrna * vliverendo * c_mrna_liverendo +
      clin_liverinter * c_mrna_liverinter -
      clout_endo * c_mrna_liverendo

    d/dt(mrna_liverinter) <-
      clout_endo * c_mrna_liverendo -
      clin_liverinter * c_mrna_liverinter +
      ll * (1 - sigma_v) * c_mrna_livervas +
      krelease * vliverinter * c_lnp_liverinter -
      ll * (1 - sigma_l) * c_mrna_liverinter -
      ktrans * vliverinter * c_mrna_liverinter

    d/dt(mrna_rem) <-
      qt * c_mrna_plasma - (qt - lt) * c_mrna_rem -
      kdeg_mrna * vrem * c_mrna_rem -
      lt * (1 - sigma_l) * c_mrna_rem

    # ===============================================================
    # 7. Intracellular ribonucleoprotein assembly -- Suppl Info S1,
    #    "Ribonucleoprotein Complex (RNP)". Degradation of the RNP
    #    complex is assumed negligible (Model assumptions, point 4).
    # ===============================================================
    d/dt(sgrna_cell) <-
      krelease * vliverinter * c_sgrna_liverinter -
      kon_rnp * c_sgrna_cell * c_cas9 * vliverinter +
      koff_rnp * vliverinter * c_rnp

    d/dt(cas9) <-
      ktrans * vliverinter * c_mrna_liverinter +
      koff_rnp * vliverinter * c_rnp -
      kon_rnp * c_sgrna_cell * c_cas9 * vliverinter

    d/dt(rnp) <-
      kon_rnp * c_sgrna_cell * c_cas9 * vliverinter -
      koff_rnp * vliverinter * c_rnp

    # ===============================================================
    # 8. Pharmacodynamics
    # ===============================================================
    # All PD rate constants are reported per day in Table 2 while the
    # model time unit is hours, so each is divided by 24.
    kout_ttr <- exp(lkout_ttr) / 24
    imax_ttr <- exp(limax_ttr)
    ic50_ttr <- exp(lic50_ttr)
    hill_ttr <- exp(lhill_ttr)
    kin_ttr <- kout_ttr * rbase_ttr

    mtt <- exp(lmtt) * 24
    ktr <- 2 / mtt
    imax_pcsk9 <- exp(limax_pcsk9)
    ic50_pcsk9 <- exp(lic50_pcsk9)
    gamma_ldlc <- exp(lgamma_ldlc)
    kdeg_ldlc <- exp(lkdeg_ldlc) / 24
    ksyn_ldlc <- kdeg_ldlc * rbase_ldlc

    # TTR: indirect-response model I (inhibition of production) driven
    # by the liver RNP concentration -- Suppl Info S1, "Reduction of
    # TTR proteins". Table 2 reports an estimated gamma coefficient for
    # the NHP TTR model (0.31, 13.3% RSE) that the printed equation
    # omits; it is applied here in the canonical sigmoidal Imax form,
    # which collapses to the printed equation when hill = 1.
    inh_ttr <- imax_ttr * c_rnp^hill_ttr / (ic50_ttr^hill_ttr + c_rnp^hill_ttr)
    d/dt(ttr) <- kin_ttr * (1 - inh_ttr) - kout_ttr * ttr
    ttr(0) <- rbase_ttr

    # PCSK9: Friberg-type transit / feedback model -- Suppl Info S1,
    # "Reduction of PCSK9 and LDL cholesterol". k_prol and k_circ are
    # not tabulated; both are pinned to k_tr by the stationarity of the
    # printed initial conditions (Prol0 = T1_0 = Circ0 = 100).
    inh_pcsk9 <- imax_pcsk9 * c_rnp / (ic50_pcsk9 + c_rnp)
    d/dt(pcsk9_prol) <-
      ktr * pcsk9_prol * (1 - inh_pcsk9) *
      (rbase_pcsk9 / pcsk9_circ)^gamma_pcsk9 -
      ktr * pcsk9_prol
    d/dt(pcsk9_transit1) <- ktr * pcsk9_prol - ktr * pcsk9_transit1
    d/dt(pcsk9_circ) <- ktr * pcsk9_transit1 - ktr * pcsk9_circ
    pcsk9_prol(0) <- rbase_pcsk9
    pcsk9_transit1(0) <- rbase_pcsk9
    pcsk9_circ(0) <- rbase_pcsk9

    # LDL cholesterol: precursor-dependent model driven by circulating
    # PCSK9. The synthesis rate is not tabulated; it is pinned by the
    # printed initial condition LDL0 = 100 at Circ = Circ0.
    d/dt(ldlc) <-
      ksyn_ldlc * (pcsk9_circ / rbase_pcsk9)^gamma_ldlc -
      kdeg_ldlc * ldlc
    ldlc(0) <- rbase_ldlc

    # ===============================================================
    # 9. Initial condition for the free LDL receptor pool
    # ===============================================================
    ldlr(0) <- ldltot * vliverendo

    # ===============================================================
    # 10. Observations
    # ===============================================================
    Cc <- c_lnp_plasma
    Cc_sgrna <- c_sgrna_plasma
    Cc_mrna <- c_mrna_plasma
    TTR <- ttr
    PCSK9 <- pcsk9_circ
    LDLC <- ldlc

    Cc ~ prop(propSd)
    Cc_sgrna ~ prop(propSd_Cc_sgrna)
    Cc_mrna ~ prop(propSd_Cc_mrna)
    TTR ~ prop(propSd_TTR)
    PCSK9 ~ prop(propSd_PCSK9)
    LDLC ~ prop(propSd_LDLC)
  })
}
