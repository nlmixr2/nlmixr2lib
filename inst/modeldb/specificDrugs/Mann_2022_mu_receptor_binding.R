Mann_2022_mu_receptor_binding <- function() {
  description <- paste(
    "QSP. Competitive mu-opioid receptor binding kinetics layer of the",
    "Mann 2022 translational opioid-overdose model. Tracks the fraction",
    "of receptors bound by an opioid agonist (RL_op) and by an opioid",
    "antagonist (RL_antag) under simultaneous exposure to both ligands.",
    "All 12 ligands characterised by Mann 2022 Supplement 1 Table S2 are",
    "carried inline as fixed parameters (Kon in pM^-n s^-1, Koff in",
    "s^-1, slope n unitless); the OPIOID_ID and ANTAGONIST_ID integer",
    "covariates select which ligand occupies each binding slot at",
    "simulation time, so the same compiled model can simulate any",
    "agonist-antagonist pair from the Table-S2 panel without re-",
    "instantiation. Ligand effect-site concentrations enter the model",
    "as the time-varying covariate columns L_OPIOID_pM and",
    "L_ANTAGONIST_pM, typically piped from the PK layer in a composed",
    "chain (e.g., Mann_2022_fentanyl_iv or Mann_2022_carfentanil_iv for",
    "the opioid slot; Laffont_2024_naloxone or Laffont_2024_nalmefene",
    "for the antagonist slot)."
  )
  reference <- paste(
    "Mann J, Samieegohar M, Chaturbedi A, Zirkle J, Han X, Ahmadi SF,",
    "Eshleman A, Janowsky A, Wolfrum K, Swanson T, Bloom S, Dahan A,",
    "Olofsen E, Florian J, Strauss DG, Li Z.",
    "Development of a Translational Model to Assess the Impact of Opioid",
    "Overdose and Naloxone Dosing on Respiratory Depression and",
    "Cardiac Arrest. Clin Pharmacol Ther. 2022;112(5):1020-1032.",
    "doi:10.1002/cpt.2696. PMID 35766413.",
    "Binding parameter table (Table S2) is the FDA Division of Applied",
    "Regulatory Science's published fit of association/dissociation",
    "assays at 12 ligands against transfected human mu-opioid receptor",
    "in rat C6 glioma cells (C6-hMOR). Equation 1 (single-ligand) and",
    "the implicit multi-ligand competitive extension are stated in",
    "Supplement 1 Receptor Binding Component section.",
    "Nalmefene (13th ligand) Kon, Koff, n values are from Laffont CM,",
    "Purohit P, Delcamp N, Gonzalez-Garcia I, Skolnick P. Comparison of",
    "intranasal naloxone and intranasal nalmefene in a translational",
    "model assessing the impact of synthetic opioid overdose on",
    "respiratory depression and cardiac arrest. Front Psychiatry.",
    "2024;15:1399803. doi:10.3389/fpsyt.2024.1399803, Supplementary",
    "Table S3 (scaling approach: Laffont 2024 Methods 2.3 multiplied",
    "the Mann 2022 naloxone Kon and Koff by the nalmefene/naloxone",
    "ratios reported in Cassel et al. 2005, [3H]Alvimopan binding to",
    "the mu-opioid receptor: comparative binding kinetics of opioid",
    "antagonists. Eur J Pharmacol. 2005;520(1-3):29-36.",
    "doi:10.1016/j.ejphar.2005.08.008; the steepness parameter n was",
    "assumed identical to naloxone)."
  )
  vignette <- "Laffont_2025_opioid_overdose_reversal_simulation"
  units <- list(
    time          = "min",
    dosing        = "(not applicable; ligand concentrations are time-varying covariates)",
    concentration = "fraction (RL_op and RL_antag are fractions of total mu-opioid receptor pool)"
  )

  covariateData <- list(
    OPIOID_ID = list(
      description        = "Integer 1..13 selecting which Table-S2 (Mann 2022) or Table-S3 (Laffont 2024) ligand fills the opioid agonist slot at runtime",
      units              = "(categorical)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Mapping (Mann 2022 Supplement 1 Table S2 row order, alphabetical",
        "after first letter, preserved verbatim in this model file):",
        "1 = alfentanil; 2 = buprenorphine; 3 = butyryl fentanyl;",
        "4 = carfentanil; 5 = fluorobutyryl fentanyl;",
        "6 = fentanyl; 7 = fluoroisobutyryl fentanyl;",
        "8 = furanyl fentanyl; 9 = isobutyryl fentanyl;",
        "10 = naloxone; 11 = remifentanil; 12 = sufentanil;",
        "13 = nalmefene (added by task 131 from Laffont 2024",
        "Supplementary Table S3 scaling approach over Cassel 2005).",
        "Used in model() to gate which Kon_<ligand>, Koff_<ligand>,",
        "n_<ligand> values are routed into the opioid binding ODE."
      ),
      source_name        = "(none; new canonical covariate registered for this model)"
    ),
    ANTAGONIST_ID = list(
      description        = "Integer 1..13 selecting which Table-S2 (Mann 2022) or Table-S3 (Laffont 2024) ligand fills the opioid antagonist slot at runtime",
      units              = "(categorical)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Same integer-to-ligand mapping as OPIOID_ID. In the Mann 2022",
        "context only naloxone (10) is used as the antagonist; in the",
        "Laffont 2024 / Laffont 2025 expansion, nalmefene (13) is the",
        "second antagonist option, parameterised from Laffont 2024",
        "Supplementary Table S3. Buprenorphine (2) is a partial agonist",
        "with the slowest dissociation kinetics in Table S2; routing it",
        "into the antagonist slot models buprenorphine occupying the",
        "receptor without producing the full ventilatory-depression",
        "effect coupled to the agonist slot."
      ),
      source_name        = "(none; new canonical covariate registered for this model)"
    ),
    L_OPIOID_pM = list(
      description        = "Time-varying opioid agonist effect-site concentration (pM) feeding the opioid binding ODE",
      units              = "pM",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "In the Mann 2022 chain, the upstream IV-opioid PK layer",
        "(Mann_2022_fentanyl_iv or Mann_2022_carfentanil_iv) exposes its",
        "effect-site concentration Ce_pM as this covariate. The pM",
        "scale matches the Table S2 Kon units (pM^-n s^-1). Standalone",
        "use: supply L_OPIOID_pM as a time-varying data column on the",
        "subject's records."
      ),
      source_name        = "(none; new canonical covariate registered for this model)"
    ),
    L_ANTAGONIST_pM = list(
      description        = "Time-varying opioid antagonist effect-site concentration (pM) feeding the antagonist binding ODE",
      units              = "pM",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Antagonist analogue of L_OPIOID_pM. In the Mann 2022 +",
        "Laffont 2024 chain, the upstream naloxone or nalmefene PK",
        "layer exposes its effect-site concentration in pM as this",
        "covariate."
      ),
      source_name        = "(none; new canonical covariate registered for this model)"
    )
  )

  population <- list(
    species        = "in vitro (rat C6 glioma cells stably transfected with human mu-opioid receptor cDNA, C6-hMOR cell line) for binding-parameter fitting; the receptor-binding layer is then deployed in the Mann 2022 human-overdose simulation chain.",
    n_subjects     = NA_integer_,
    n_studies      = 1L,
    age_range      = NA_character_,
    weight_range   = NA_character_,
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = "In vitro receptor binding assay; not a clinical population.",
    dose_range     = paste(
      "Ligand concentrations at the binding-assay scale ranged from",
      "0.01 to 81.7 nM across ligands (per-ligand ranges in Mann 2022",
      "Supplement 1 Experimental Procedures). In the integrated overdose",
      "chain, ligand concentrations are driven by the upstream PK layer."
    ),
    regions        = NA_character_,
    notes          = paste(
      "Parameter values are best-fit point estimates from association",
      "and dissociation experiments simultaneously fit to Equation 1 of",
      "Supplement 1: dRL/dt = Kon * L^n * R - Koff * RL, with",
      "R = 1 - RL (fraction). For multi-ligand competition, the same",
      "equation is applied per ligand against the shared free-receptor",
      "fraction R_free = 1 - RL_op - RL_antag (Mann 2022 Results,",
      "Independent validation of the receptor binding component).",
      "Free-ligand approximation: the Mann 2022 binding assays observed",
      "<= 5% ligand loss to binding even at the lowest concentrations,",
      "so L is treated as an exogenous input rather than a state.",
      "95% confidence intervals for Kon/Koff/n are tabulated in Table S2",
      "but not encoded here -- the registry carries point estimates",
      "only (the binding model is a typical-value mechanism without",
      "between-subject variability)."
    ),
    n_ligands      = 13L
  )

  ini({
    # 12-ligand Table S2 binding-kinetic parameter set. Mann 2022
    # reports Kon in pM^-n s^-1 and Koff in s^-1; the values are
    # preserved here exactly as published so each ini() line is a
    # one-to-one round-trip of Table S2. The per-second to per-minute
    # conversion happens inside model() at the dispatch step.

    # 1 - alfentanil (Table S2)
    kon_alfentanil  <- fixed(3.32e-07)
    label("Alfentanil association rate Kon (pM^-n s^-1, FIXED, Table S2)")  # Table S2 row 1
    koff_alfentanil <- fixed(1.41e-02)
    label("Alfentanil dissociation rate Koff (s^-1, FIXED, Table S2)")      # Table S2 row 1
    n_alfentanil    <- fixed(1.0)
    label("Alfentanil binding slope n (unitless, FIXED, Table S2)")         # Table S2 row 1

    # 2 - buprenorphine (Table S2)
    kon_buprenorphine  <- fixed(4.83e-05)
    label("Buprenorphine association rate Kon (pM^-n s^-1, FIXED, Table S2)")  # Table S2 row 2
    koff_buprenorphine <- fixed(1.37e-04)
    label("Buprenorphine dissociation rate Koff (s^-1, FIXED, Table S2)")      # Table S2 row 2; slowest dissociation in panel (t1/2 = 84.4 min)
    n_buprenorphine    <- fixed(0.5)
    label("Buprenorphine binding slope n (unitless, FIXED, Table S2)")         # Table S2 row 2

    # 3 - butyryl fentanyl (Table S2)
    kon_butyrylfentanyl  <- fixed(1.14e-05)
    label("Butyryl-fentanyl association rate Kon (pM^-n s^-1, FIXED, Table S2)")  # Table S2 row 3
    koff_butyrylfentanyl <- fixed(5.32e-03)
    label("Butyryl-fentanyl dissociation rate Koff (s^-1, FIXED, Table S2)")      # Table S2 row 3
    n_butyrylfentanyl    <- fixed(0.99)
    label("Butyryl-fentanyl binding slope n (unitless, FIXED, Table S2)")         # Table S2 row 3

    # 4 - carfentanil (Table S2)
    kon_carfentanil  <- fixed(9.95e-06)
    label("Carfentanil association rate Kon (pM^-n s^-1, FIXED, Table S2)")  # Table S2 row 4
    koff_carfentanil <- fixed(2.47e-04)
    label("Carfentanil dissociation rate Koff (s^-1, FIXED, Table S2)")      # Table S2 row 4; second-slowest dissociation in panel (t1/2 = 46.6 min) -- driver of difficult-to-reverse carfentanil overdose
    n_carfentanil    <- fixed(1.023)
    label("Carfentanil binding slope n (unitless, FIXED, Table S2)")         # Table S2 row 4

    # 5 - fluorobutyryl fentanyl (Table S2)
    kon_fluorobutyrylfentanyl  <- fixed(3.96e-04)
    label("Fluorobutyryl-fentanyl association rate Kon (pM^-n s^-1, FIXED, Table S2)")  # Table S2 row 5
    koff_fluorobutyrylfentanyl <- fixed(1.54e-02)
    label("Fluorobutyryl-fentanyl dissociation rate Koff (s^-1, FIXED, Table S2)")      # Table S2 row 5
    n_fluorobutyrylfentanyl    <- fixed(0.5)
    label("Fluorobutyryl-fentanyl binding slope n (unitless, FIXED, Table S2)")         # Table S2 row 5

    # 6 - fentanyl (Table S2)
    kon_fentanyl  <- fixed(3.08e-05)
    label("Fentanyl association rate Kon (pM^-n s^-1, FIXED, Table S2)")  # Table S2 row 6
    koff_fentanyl <- fixed(4.33e-03)
    label("Fentanyl dissociation rate Koff (s^-1, FIXED, Table S2)")      # Table S2 row 6; dissociation t1/2 = 2.66 min
    n_fentanyl    <- fixed(0.844)
    label("Fentanyl binding slope n (unitless, FIXED, Table S2)")         # Table S2 row 6

    # 7 - fluoroisobutyryl fentanyl (Table S2)
    kon_fluoroisobutyrylfentanyl  <- fixed(1.77e-04)
    label("Fluoroisobutyryl-fentanyl association rate Kon (pM^-n s^-1, FIXED, Table S2)")  # Table S2 row 7
    koff_fluoroisobutyrylfentanyl <- fixed(6.29e-03)
    label("Fluoroisobutyryl-fentanyl dissociation rate Koff (s^-1, FIXED, Table S2)")      # Table S2 row 7
    n_fluoroisobutyrylfentanyl    <- fixed(0.5)
    label("Fluoroisobutyryl-fentanyl binding slope n (unitless, FIXED, Table S2)")         # Table S2 row 7

    # 8 - furanyl fentanyl (Table S2)
    kon_furanylfentanyl  <- fixed(3.42e-05)
    label("Furanyl-fentanyl association rate Kon (pM^-n s^-1, FIXED, Table S2)")  # Table S2 row 8
    koff_furanylfentanyl <- fixed(5.43e-03)
    label("Furanyl-fentanyl dissociation rate Koff (s^-1, FIXED, Table S2)")      # Table S2 row 8
    n_furanylfentanyl    <- fixed(1.01)
    label("Furanyl-fentanyl binding slope n (unitless, FIXED, Table S2)")         # Table S2 row 8

    # 9 - isobutyryl fentanyl (Table S2)
    kon_isobutyrylfentanyl  <- fixed(1.27e-04)
    label("Isobutyryl-fentanyl association rate Kon (pM^-n s^-1, FIXED, Table S2)")  # Table S2 row 9
    koff_isobutyrylfentanyl <- fixed(3.11e-02)
    label("Isobutyryl-fentanyl dissociation rate Koff (s^-1, FIXED, Table S2)")      # Table S2 row 9
    n_isobutyrylfentanyl    <- fixed(0.69)
    label("Isobutyryl-fentanyl binding slope n (unitless, FIXED, Table S2)")         # Table S2 row 9

    # 10 - naloxone (Table S2)
    kon_naloxone  <- fixed(1.67e-04)
    label("Naloxone association rate Kon (pM^-n s^-1, FIXED, Table S2)")  # Table S2 row 10; primary antagonist in Mann 2022 overdose-rescue chain
    koff_naloxone <- fixed(3.96e-02)
    label("Naloxone dissociation rate Koff (s^-1, FIXED, Table S2)")      # Table S2 row 10; dissociation t1/2 = 0.29 min
    n_naloxone    <- fixed(0.86)
    label("Naloxone binding slope n (unitless, FIXED, Table S2)")         # Table S2 row 10

    # 11 - remifentanil (Table S2)
    kon_remifentanil  <- fixed(8.08e-06)
    label("Remifentanil association rate Kon (pM^-n s^-1, FIXED, Table S2)")  # Table S2 row 11
    koff_remifentanil <- fixed(2.08e-03)
    label("Remifentanil dissociation rate Koff (s^-1, FIXED, Table S2)")      # Table S2 row 11
    n_remifentanil    <- fixed(0.70)
    label("Remifentanil binding slope n (unitless, FIXED, Table S2)")         # Table S2 row 11

    # 12 - sufentanil (Table S2)
    kon_sufentanil  <- fixed(1.16e-05)
    label("Sufentanil association rate Kon (pM^-n s^-1, FIXED, Table S2)")  # Table S2 row 12
    koff_sufentanil <- fixed(1.07e-03)
    label("Sufentanil dissociation rate Koff (s^-1, FIXED, Table S2)")      # Table S2 row 12
    n_sufentanil    <- fixed(1.02)
    label("Sufentanil binding slope n (unitless, FIXED, Table S2)")         # Table S2 row 12

    # 13 - nalmefene (Laffont 2024 Supplement Table S3; added by task 131
    # since Mann 2022 Table S2 did not include nalmefene). Values derived
    # by Laffont 2024 Methods 2.3 scaling approach: nalmefene-to-naloxone
    # Kon and Koff ratios from Cassel 2005 multiplied by the Mann 2022
    # naloxone Kon and Koff to keep the binding panel internally
    # consistent. n was assumed identical to naloxone (0.86).
    kon_nalmefene  <- fixed(2.06e-04)
    label("Nalmefene association rate Kon (pM^-n s^-1, FIXED, Laffont 2024 Supp Table S3)")  # Laffont 2024 Supp Table S3: nalmefene kon = 2.06e-04 pM^-n s^-1 (scaling from Cassel 2005 over Mann 2022 naloxone)
    koff_nalmefene <- fixed(1.35e-02)
    label("Nalmefene dissociation rate Koff (s^-1, FIXED, Laffont 2024 Supp Table S3)")      # Laffont 2024 Supp Table S3: nalmefene koff = 1.35e-02 s^-1; ~2.9x slower than naloxone (3.96e-02) -> dissociation t1/2 = 0.86 min vs naloxone 0.29 min
    n_nalmefene    <- fixed(0.86)
    label("Nalmefene binding slope n (unitless, FIXED, assumed identical to naloxone)")     # Laffont 2024 Supp Table S3 / Methods 2.3: n for nalmefene assumed equal to naloxone (0.86)
  })

  model({
    # Per-second to per-minute conversion (Table S2 reports rates per
    # second; the binding ODE is integrated on a per-minute time scale to
    # match the Mann 2022 PK layers that compose with this model).
    sec_to_min <- 60.0

    # OPIOID slot: dispatch by integer OPIOID_ID covariate (1..12 per
    # Table S2 row order documented in covariateData$OPIOID_ID$notes).
    # When OPIOID_ID is out of range (e.g., 0 or NA-encoded), every
    # indicator evaluates to FALSE = 0 and the slot rate is zero (no
    # binding from that slot), so the model degrades safely.
    Kon_op_per_s <-
      (OPIOID_ID ==  1) * kon_alfentanil               +
      (OPIOID_ID ==  2) * kon_buprenorphine            +
      (OPIOID_ID ==  3) * kon_butyrylfentanyl          +
      (OPIOID_ID ==  4) * kon_carfentanil              +
      (OPIOID_ID ==  5) * kon_fluorobutyrylfentanyl    +
      (OPIOID_ID ==  6) * kon_fentanyl                 +
      (OPIOID_ID ==  7) * kon_fluoroisobutyrylfentanyl +
      (OPIOID_ID ==  8) * kon_furanylfentanyl          +
      (OPIOID_ID ==  9) * kon_isobutyrylfentanyl       +
      (OPIOID_ID == 10) * kon_naloxone                 +
      (OPIOID_ID == 11) * kon_remifentanil             +
      (OPIOID_ID == 12) * kon_sufentanil               +
      (OPIOID_ID == 13) * kon_nalmefene
    Koff_op_per_s <-
      (OPIOID_ID ==  1) * koff_alfentanil               +
      (OPIOID_ID ==  2) * koff_buprenorphine            +
      (OPIOID_ID ==  3) * koff_butyrylfentanyl          +
      (OPIOID_ID ==  4) * koff_carfentanil              +
      (OPIOID_ID ==  5) * koff_fluorobutyrylfentanyl    +
      (OPIOID_ID ==  6) * koff_fentanyl                 +
      (OPIOID_ID ==  7) * koff_fluoroisobutyrylfentanyl +
      (OPIOID_ID ==  8) * koff_furanylfentanyl          +
      (OPIOID_ID ==  9) * koff_isobutyrylfentanyl       +
      (OPIOID_ID == 10) * koff_naloxone                 +
      (OPIOID_ID == 11) * koff_remifentanil             +
      (OPIOID_ID == 12) * koff_sufentanil               +
      (OPIOID_ID == 13) * koff_nalmefene
    n_op <-
      (OPIOID_ID ==  1) * n_alfentanil               +
      (OPIOID_ID ==  2) * n_buprenorphine            +
      (OPIOID_ID ==  3) * n_butyrylfentanyl          +
      (OPIOID_ID ==  4) * n_carfentanil              +
      (OPIOID_ID ==  5) * n_fluorobutyrylfentanyl    +
      (OPIOID_ID ==  6) * n_fentanyl                 +
      (OPIOID_ID ==  7) * n_fluoroisobutyrylfentanyl +
      (OPIOID_ID ==  8) * n_furanylfentanyl          +
      (OPIOID_ID ==  9) * n_isobutyrylfentanyl       +
      (OPIOID_ID == 10) * n_naloxone                 +
      (OPIOID_ID == 11) * n_remifentanil             +
      (OPIOID_ID == 12) * n_sufentanil               +
      (OPIOID_ID == 13) * n_nalmefene

    # ANTAGONIST slot: same dispatch by ANTAGONIST_ID. Both slots may be
    # set to the same integer if a single-ligand simulation is desired,
    # but the per-slot ligand-concentration covariates L_OPIOID_pM and
    # L_ANTAGONIST_pM must be set consistently.
    Kon_an_per_s <-
      (ANTAGONIST_ID ==  1) * kon_alfentanil               +
      (ANTAGONIST_ID ==  2) * kon_buprenorphine            +
      (ANTAGONIST_ID ==  3) * kon_butyrylfentanyl          +
      (ANTAGONIST_ID ==  4) * kon_carfentanil              +
      (ANTAGONIST_ID ==  5) * kon_fluorobutyrylfentanyl    +
      (ANTAGONIST_ID ==  6) * kon_fentanyl                 +
      (ANTAGONIST_ID ==  7) * kon_fluoroisobutyrylfentanyl +
      (ANTAGONIST_ID ==  8) * kon_furanylfentanyl          +
      (ANTAGONIST_ID ==  9) * kon_isobutyrylfentanyl       +
      (ANTAGONIST_ID == 10) * kon_naloxone                 +
      (ANTAGONIST_ID == 11) * kon_remifentanil             +
      (ANTAGONIST_ID == 12) * kon_sufentanil               +
      (ANTAGONIST_ID == 13) * kon_nalmefene
    Koff_an_per_s <-
      (ANTAGONIST_ID ==  1) * koff_alfentanil               +
      (ANTAGONIST_ID ==  2) * koff_buprenorphine            +
      (ANTAGONIST_ID ==  3) * koff_butyrylfentanyl          +
      (ANTAGONIST_ID ==  4) * koff_carfentanil              +
      (ANTAGONIST_ID ==  5) * koff_fluorobutyrylfentanyl    +
      (ANTAGONIST_ID ==  6) * koff_fentanyl                 +
      (ANTAGONIST_ID ==  7) * koff_fluoroisobutyrylfentanyl +
      (ANTAGONIST_ID ==  8) * koff_furanylfentanyl          +
      (ANTAGONIST_ID ==  9) * koff_isobutyrylfentanyl       +
      (ANTAGONIST_ID == 10) * koff_naloxone                 +
      (ANTAGONIST_ID == 11) * koff_remifentanil             +
      (ANTAGONIST_ID == 12) * koff_sufentanil               +
      (ANTAGONIST_ID == 13) * koff_nalmefene
    n_an <-
      (ANTAGONIST_ID ==  1) * n_alfentanil               +
      (ANTAGONIST_ID ==  2) * n_buprenorphine            +
      (ANTAGONIST_ID ==  3) * n_butyrylfentanyl          +
      (ANTAGONIST_ID ==  4) * n_carfentanil              +
      (ANTAGONIST_ID ==  5) * n_fluorobutyrylfentanyl    +
      (ANTAGONIST_ID ==  6) * n_fentanyl                 +
      (ANTAGONIST_ID ==  7) * n_fluoroisobutyrylfentanyl +
      (ANTAGONIST_ID ==  8) * n_furanylfentanyl          +
      (ANTAGONIST_ID ==  9) * n_isobutyrylfentanyl       +
      (ANTAGONIST_ID == 10) * n_naloxone                 +
      (ANTAGONIST_ID == 11) * n_remifentanil             +
      (ANTAGONIST_ID == 12) * n_sufentanil               +
      (ANTAGONIST_ID == 13) * n_nalmefene

    # Convert per-second to per-minute time scale (matches the per-minute
    # time scale used by the Mann 2022 PK layers Mann_2022_fentanyl_iv /
    # Mann_2022_carfentanil_iv that compose with this binding model).
    Kon_op  <- Kon_op_per_s  * sec_to_min
    Koff_op <- Koff_op_per_s * sec_to_min
    Kon_an  <- Kon_an_per_s  * sec_to_min
    Koff_an <- Koff_an_per_s * sec_to_min

    # Free-receptor fraction R_free = 1 - RL_op - RL_antag. Equation 1 of
    # Mann 2022 Supplement 1, multi-ligand competitive extension stated
    # implicitly in Results (Independent validation of the receptor
    # binding component): both ligands compete for the same shared pool.
    R_free <- 1 - RL_op - RL_antag

    # Multi-ligand competitive binding ODEs (Mann 2022 Supplement 1
    # Equation 1, per slot). L_OPIOID_pM and L_ANTAGONIST_pM enter as
    # time-varying covariates from the upstream PK layer (or from the
    # subject's data records in standalone use).
    d/dt(RL_op)    <- Kon_op * (L_OPIOID_pM    )^n_op - Koff_op * RL_op +
                      -Kon_op * (L_OPIOID_pM   )^n_op * (RL_op + RL_antag)
    d/dt(RL_antag) <- Kon_an * (L_ANTAGONIST_pM)^n_an - Koff_an * RL_antag +
                      -Kon_an * (L_ANTAGONIST_pM)^n_an * (RL_op + RL_antag)
    # The two expressions above are the algebraic expansion of
    # Kon * L^n * R_free - Koff * RL with R_free = 1 - RL_op - RL_antag;
    # they are written in expanded form because rxode2 derivative-rhs
    # syntax composes more cleanly without referencing R_free directly
    # in the d/dt() rhs (R_free is still defined above for output use).

    # Outputs (named per Mann 2022 paper convention):
    CAR <- RL_op                  # fraction of mu receptors bound by the opioid agonist
    CBR <- RL_antag               # fraction of mu receptors bound by the antagonist
  })
}
