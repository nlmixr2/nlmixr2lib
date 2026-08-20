Desai_2024_crisprCas9_human <- function() {
  description <- paste(
    "QSP. Translational quantitative systems pharmacology platform",
    "for in vivo CRISPR-Cas9 gene editing",
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
    "with receptor recycling. This is the clinical arm of the platform",
    "(step 4 of the paper's four-step workflow), scaled to human",
    "physiology and re-estimated against NTLA-2001 LNP plasma PK and",
    "serum transthyretin from the first-in-human study (0.1-1 mg/kg",
    "total RNA by 2 h IV infusion). The pharmacodynamic layer is an",
    "indirect-response model in which the liver ribonucleoprotein",
    "concentration inhibits TTR production; unlike the NHP arm the",
    "human TTR model carries no Hill coefficient. Sibling extractions:",
    "Desai_2024_crisprCas9_mouse and Desai_2024_crisprCas9_nhp.",
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
    "drug-specific parameters from Table 2. Observed human LNP plasma",
    "PK is from Abdelhady et al. (2023) and observed serum TTR from",
    "Gane et al. (2022) J Hepatology 77:S58-S59;",
    "doi:10.1016/s0168-8278(22)00520-7.",
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
    "cas9", "rnp", "ttr"
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
      "mRNA and the paired LNP dose is given in Table 1 (e.g. 0.3",
      "mg/kg total RNA in a 71 kg adult = 7093 ug sgRNA + 14207 ug",
      "mRNA with 5.55 mg/kg = 394050 ug LNP).",
      sep = " "
    ),
    concentration = paste(
      "ug/mL for all LNP, sgRNA, mRNA, Cas9 and RNP concentrations.",
      "The PD output (TTR) is percent of baseline.",
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
    ttr = list(analyte = "transthyretin", units = "% of baseline", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 2L,
    weight_median = "71 kg",
    disease_state = "Adults with hereditary transthyretin amyloidosis treated with NTLA-2001.",
    dose_range = paste(
      "0.1, 0.3, 0.7 and 1 mg/kg total RNA (1.85, 5.55, 12.9 and 18.5",
      "mg/kg LNP) by 2 h IV infusion.",
      sep = " "
    ),
    notes = paste(
      "Desai 2024 Table 1. No individual-level cohort is described: the",
      "analysis was performed on mean LNP plasma profiles digitised",
      "from Abdelhady 2023 (sampled 0-50 h) and mean serum TTR",
      "reduction digitised from Gane 2022 (sampled at 7, 14 and 28",
      "days) with WebPlotDigitizer, so n_subjects is not reported. The",
      "reference body weight for the human physiology in Supplementary",
      "Tables S1 and S2 is 71 kg.",
      sep = " "
    )
  )

  ini({
    # ---------------------------------------------------------------
    # Estimated drug-specific attributes -- Table 2, Human column
    # ---------------------------------------------------------------
    lkin_endo <- log(0.007); label("Rate of endocytosis for LNP, sgRNA and mRNA (1/h)")   # Table 2, human (1.44% RSE)
    lkout_exo <- log(775); label("Rate of exocytosis for LNP, sgRNA and mRNA (1/h)")      # Table 2, human (0.602% RSE)
    lkdeg_dr <- log(5.15); label("Rate of degradation of the LNP-LDL receptor complex (1/h)")  # Table 2, human (6.47% RSE)
    lldltot <- log(84.5); label("Total LDL receptor concentration (ug/mL)")               # Table 2, human (0.508% RSE)
    lkdis <- log(0.186); label("Rate of dissociation from opsonins (1/h)")                # Table 2, human (27.2% RSE)
    lkass <- log(56.8); label("Rate of association to opsonins (1/h)")                    # Table 2, human (16% RSE)
    lkdeg_lnp <- log(0.101); label("Rate of degradation of unbound LNP (1/h)")            # Table 2, human (35.9% RSE)
    lkel <- log(0.009); label("Rate of elimination of the LDL receptor (1/h)")            # Table 2, human (44.7% RSE)

    # ---------------------------------------------------------------
    # Drug-specific attributes held constant -- Table 2, Human column
    # ---------------------------------------------------------------
    krelease <- fixed(0.0056); label("Rate of release of sgRNA and mRNA from LNP (1/h); literature value from Miyazawa 2024")  # Table 2, human
    kdeg_sgrna <- fixed(2.01); label("Rate of degradation of sgRNA (1/h); taken from the step-1 mechanistic model")    # Table 2, human footnote c
    kdeg_mrna <- fixed(0.1232); label("Rate of degradation of mRNA (1/h); literature value from Miyazawa 2024")        # Table 2, human
    koff_rnp <- fixed(0.00188); label("Rate of dissociation of the RNP complex (1/h); taken from the step-1 mechanistic model")  # Table 2, human footnote c
    kd <- fixed(0.49); label("Equilibrium dissociation constant for Cas9-sgRNA (nM); literature value from Sternberg 2014")      # Table 2
    kint <- fixed(0.9063); label("Rate of phagocytosis into the mononuclear phagocyte system (1/h); literature value from Miyazawa 2024")  # Table 2
    kon_lnp <- fixed(0.18); label("Rate of association of LNP with the LDL receptor (mL/ug/h); literature value from Harwood 1997")       # Table 2
    koff_lnp <- fixed(33.12); label("Rate of dissociation of the LNP-LDL receptor complex (1/h); literature value from Harwood 1997")     # Table 2
    fu_lnp <- fixed(0.002); label("Plasma free fraction of LNP (unitless); literature value from Mager 2012")   # Table 2
    fu_rna <- fixed(0.15); label("Plasma free fraction of sgRNA and mRNA (unitless); literature value from Ayyar 2021")  # Table 2
    ktrans <- fixed(0.36); label("Rate of translation from mRNA to Cas9 (1/h); taken from the step-1 mechanistic model")  # Table 2, human footnote c

    # ---------------------------------------------------------------
    # TTR pharmacodynamics (NTLA-2001) -- Table 2, Human column.
    # Table 2 reports no gamma coefficient for the human TTR model, so
    # the inhibition function is the plain Imax form printed in
    # Supplementary Information S1.
    # ---------------------------------------------------------------
    rbase_ttr <- fixed(100); label("Baseline serum TTR (% of baseline)")             # Table 2, human
    lkout_ttr <- log(0.247); label("First-order degradation rate of TTR protein (1/day)")  # Table 2, human (17.5% RSE)
    limax_ttr <- log(0.959); label("Maximum inhibition of TTR production (unitless)")      # Table 2, human (6.55% RSE)
    lic50_ttr <- log(0.3); label("Liver RNP concentration at 50% inhibition of TTR production (ug/mL)")  # Table 2, human (0.02% RSE)

    # ---------------------------------------------------------------
    # Residual error -- proportional per output (Methods, Software).
    # Desai 2024 states a proportional error model was used for every
    # output but does not report any sigma_slope value.
    # ---------------------------------------------------------------
    propSd <- fixed(0); label("Proportional residual SD for LNP plasma concentration (fraction); not published")     # Methods, Software
    propSd_Cc_sgrna <- fixed(0); label("Proportional residual SD for sgRNA plasma concentration (fraction); not published")  # Methods, Software
    propSd_Cc_mrna <- fixed(0); label("Proportional residual SD for mRNA plasma concentration (fraction); not published")    # Methods, Software
    propSd_TTR <- fixed(0); label("Proportional residual SD for serum TTR (fraction); not published")                # Methods, Software
  })

  model({
    # ===============================================================
    # 0. Species physiology -- held as structural constants of this
    #    species-specific file rather than as estimable parameters.
    #    Human, 71 kg reference body weight.
    # ===============================================================
    # Physiological volumes (mL) -- Supplementary Table S1, human column
    vplasma <- 3126     # Plasma volume in circulation (mL)
    vlymph <- 274       # Volume of lymph (mL)
    vlivervas <- 183    # Plasma volume in liver vasculature (mL)
    vliverendo <- 10.7  # Endosomal volume of liver (mL)
    vliverinter <- 429  # Cellular and interstitial volume of liver (mL)
    vmps <- 137.1       # Volume of the mononuclear phagocyte system (mL)
    vkidney <- 332      # Plasma, interstitial and cellular volume of kidney (mL)
    vrem <- 12394       # Extracellular volume of remainder tissue (mL)

    # Physiological flows (mL/h) -- Supplementary Table S2, human column
    ql <- 13210         # Plasma flow to liver (mL/h)
    ll <- 26.42         # Lymph flow to liver (mL/h)
    qt <- 97677         # Plasma flow to remainder tissue (mL/h)
    lt <- 195.354       # Lymph flow to remainder tissue (mL/h)
    qk <- 36402         # Plasma flow to kidney (mL/h)
    lk <- 72.804        # Lymph flow to kidney (mL/h)
    gfr <- 7200         # Glomerular filtration rate (mL/h)

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
    kdeg_lnp <- exp(lkdeg_lnp)
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
    kin_ttr <- kout_ttr * rbase_ttr

    # TTR: indirect-response model I (inhibition of production) driven
    # by the liver RNP concentration -- Suppl Info S1, "Reduction of
    # TTR proteins", exactly as printed. Table 2 reports no gamma
    # coefficient for the human TTR model.
    inh_ttr <- imax_ttr * c_rnp / (ic50_ttr + c_rnp)
    d/dt(ttr) <- kin_ttr * (1 - inh_ttr) - kout_ttr * ttr
    ttr(0) <- rbase_ttr

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

    Cc ~ prop(propSd)
    Cc_sgrna ~ prop(propSd_Cc_sgrna)
    Cc_mrna ~ prop(propSd_Cc_mrna)
    TTR ~ prop(propSd_TTR)
  })
}
