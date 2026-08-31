Echterhof_2026_phage_rat_pbpk <- function() {
  description <- paste(
    "Preclinical (rat).",
    "PBPK (whole-body, permeability-limited, Monolix 2024R1 SAEM, naive-pooled)",
    "for the biodistribution of intravenously administered lytic bacteriophage.",
    "INTERSPECIES-SCALING PROJECTION, not a fit: every estimated parameter is",
    "carried over unchanged from the mouse model",
    "(modellib('Echterhof_2026_phage_mouse_pbpk')), which was fitted to",
    "125I-Sulfo-SHPP-labelled Luz24, PAML-31-1 and OMKO1 biodistribution in",
    "CD-1 mice, and only the anatomy and the tissue-resident phagocyte",
    "densities are replaced with the rat column of Table S3. The paper used this",
    "configuration to externally qualify the model against digitised Klebsiella",
    "pneumoniae phage data in rats, and reports a significant overprediction at",
    "the highest 10^11 PFU dose while the 10^9 and 10^6 PFU groups were more",
    "accurately predicted -- the authors attribute this to model misspecification",
    "of the saturable clearance processes. Twelve organs (lung, muscle, bone, liver,",
    "stomach, spleen, small intestine, large intestine, brain, skin, kidney and",
    "a carcass mass-balance remainder) plus venous and arterial blood. Each",
    "organ carries a permeability-limited vascular / interstitial pair and a",
    "saturable mononuclear-phagocyte-system (MPS) pool whose capacity is set by",
    "the tissue-resident phagocyte count; the brain carries no MPS pool (the",
    "supplement's brain equations omit it). Elimination is by MPS degradation",
    "and by a hypothetical first-order 'active' clearance at epithelial",
    "surfaces, which drains the kidney to urine and the liver, stomach and",
    "intestines to a gut lumen (the measured 'stomach contents'). The MPS",
    "capacity is expressed per 100,000 phagocytes in PFU, so the administered",
    "dose DOSE_PHAGE_PFU sets the saturation scale of the %ID state system.",
    "IMPORTANT: the fitted observable was 125I scintillation counts, i.e.",
    "labelled phage PLUS free label, and was co-generated with a twin free-125I",
    "PBPK model that the publication describes only in prose. That twin model is",
    "NOT reproducible from any published source and is NOT implemented here, so",
    "this file reproduces the phage-disposition system alone. See the vignette",
    "Errata.",
    sep = " "
  )
  reference <- paste(
    "Echterhof A, Dharmaraj T, Blankenberg P, Targ B, Nguyen TD, Bollyky PL,",
    "Smith NM, Blankenberg FG (2026).",
    "Whole-body distribution of three Pseudomonas phages characterized by a",
    "translational physiologically based pharmacokinetic model.",
    "Antimicrob Agents Chemother 70(1):e01506-25.",
    "doi:10.1128/aac.01506-25. PMCID PMC12777566.",
    "Structural parameter estimates: Table 1 (NOT the abstract or the Results",
    "narrative, which print a different, stale parameter set -- see the vignette",
    "Errata).",
    "Per-organ differential equations and the anatomical parameter table:",
    "Supplemental material AAC01506-25-S0001.docx, 'Equations' section and",
    "Table S3.",
    "MPS carrying capacity derivation: Table S2, from Bichet et al.",
    "No erratum or corrigendum was located for this article.",
    sep = " "
  )
  vignette <- "Echterhof_2026_phage_pbpk"

  # The intravenous bolus enters the venous blood pool, not a state named depot
  # or central, so buildModelDb()'s two-name heuristic cannot infer it and the
  # registry would otherwise report no dosing compartment.
  dosing <- "venous"

  units <- list(
    time          = "h",
    dosing        = "%ID",
    concentration = "%ID/g",
    amount        = "%ID",
    weight        = "kg"
  )

  # `gut_lumen` here is the "GI Elimination" box of Figure 1 -- the luminal
  # reservoir fed by active clearance from liver, stomach, small intestine and
  # large intestine and drained by first-order gastric transit. It is the state
  # the authors measured as "Stomach Contents" (Fig. S11), reported in %ID
  # rather than %ID/g. The supplement prints no equation for it; it is
  # reconstructed by analogy with the printed urine equation (see Errata).
  # Its `specimen` is recorded as "faeces" because that is the closest member
  # of the package specimen vocabulary for gastrointestinal luminal content
  # destined for excretion; the paper measured it in the stomach, not the
  # faeces, so read the label as "GI lumen".
  compartmentData <- list(
    venous                = list(analyte = "125I-labelled phage", units = "%ID", specimen = "plasma", verified = TRUE),
    int_venous            = list(analyte = "125I-labelled phage", units = "%ID", specimen = "endosome", verified = TRUE),
    arterial              = list(analyte = "125I-labelled phage", units = "%ID", specimen = "plasma", verified = TRUE),
    vp_lung               = list(analyte = "125I-labelled phage", units = "%ID", specimen = "plasma", verified = TRUE),
    is_lung               = list(analyte = "125I-labelled phage", units = "%ID", specimen = "tissue", verified = TRUE),
    int_lung              = list(analyte = "125I-labelled phage", units = "%ID", specimen = "endosome", verified = TRUE),
    vp_muscle             = list(analyte = "125I-labelled phage", units = "%ID", specimen = "plasma", verified = TRUE),
    is_muscle             = list(analyte = "125I-labelled phage", units = "%ID", specimen = "tissue", verified = TRUE),
    int_muscle            = list(analyte = "125I-labelled phage", units = "%ID", specimen = "endosome", verified = TRUE),
    vp_bone               = list(analyte = "125I-labelled phage", units = "%ID", specimen = "plasma", verified = TRUE),
    is_bone               = list(analyte = "125I-labelled phage", units = "%ID", specimen = "tissue", verified = TRUE),
    int_bone              = list(analyte = "125I-labelled phage", units = "%ID", specimen = "endosome", verified = TRUE),
    vp_liver              = list(analyte = "125I-labelled phage", units = "%ID", specimen = "plasma", verified = TRUE),
    is_liver              = list(analyte = "125I-labelled phage", units = "%ID", specimen = "tissue", verified = TRUE),
    int_liver             = list(analyte = "125I-labelled phage", units = "%ID", specimen = "endosome", verified = TRUE),
    vp_stomach            = list(analyte = "125I-labelled phage", units = "%ID", specimen = "plasma", verified = TRUE),
    is_stomach            = list(analyte = "125I-labelled phage", units = "%ID", specimen = "tissue", verified = TRUE),
    int_stomach           = list(analyte = "125I-labelled phage", units = "%ID", specimen = "endosome", verified = TRUE),
    vp_spleen             = list(analyte = "125I-labelled phage", units = "%ID", specimen = "plasma", verified = TRUE),
    is_spleen             = list(analyte = "125I-labelled phage", units = "%ID", specimen = "tissue", verified = TRUE),
    int_spleen            = list(analyte = "125I-labelled phage", units = "%ID", specimen = "endosome", verified = TRUE),
    vp_small_intestine    = list(analyte = "125I-labelled phage", units = "%ID", specimen = "plasma", verified = TRUE),
    is_small_intestine    = list(analyte = "125I-labelled phage", units = "%ID", specimen = "tissue", verified = TRUE),
    int_small_intestine   = list(analyte = "125I-labelled phage", units = "%ID", specimen = "endosome", verified = TRUE),
    vp_large_intestine    = list(analyte = "125I-labelled phage", units = "%ID", specimen = "plasma", verified = TRUE),
    is_large_intestine    = list(analyte = "125I-labelled phage", units = "%ID", specimen = "tissue", verified = TRUE),
    int_large_intestine   = list(analyte = "125I-labelled phage", units = "%ID", specimen = "endosome", verified = TRUE),
    vp_brain              = list(analyte = "125I-labelled phage", units = "%ID", specimen = "plasma", verified = TRUE),
    is_brain              = list(analyte = "125I-labelled phage", units = "%ID", specimen = "tissue", verified = TRUE),
    vp_skin               = list(analyte = "125I-labelled phage", units = "%ID", specimen = "plasma", verified = TRUE),
    is_skin               = list(analyte = "125I-labelled phage", units = "%ID", specimen = "tissue", verified = TRUE),
    int_skin              = list(analyte = "125I-labelled phage", units = "%ID", specimen = "endosome", verified = TRUE),
    vp_remainder          = list(analyte = "125I-labelled phage", units = "%ID", specimen = "plasma", verified = TRUE),
    is_remainder          = list(analyte = "125I-labelled phage", units = "%ID", specimen = "tissue", verified = TRUE),
    int_remainder         = list(analyte = "125I-labelled phage", units = "%ID", specimen = "endosome", verified = TRUE),
    vp_kidney             = list(analyte = "125I-labelled phage", units = "%ID", specimen = "plasma", verified = TRUE),
    is_kidney             = list(analyte = "125I-labelled phage", units = "%ID", specimen = "tissue", verified = TRUE),
    int_kidney            = list(analyte = "125I-labelled phage", units = "%ID", specimen = "endosome", verified = TRUE),
    urine                 = list(analyte = "125I-labelled phage", units = "%ID", specimen = "urine", verified = TRUE),
    gut_lumen             = list(analyte = "125I-labelled phage", units = "%ID", specimen = "faeces", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Scales cardiac output allometrically (CO = CO_C * WT^0.75), every organ",
        "and blood volume linearly (V = fv * WT), and the active surface",
        "clearance linearly (CL = CL_C * WT, the L/h/kg unit of Table 1).",
        "Table S3 tabulates 0.25 kg as the generic rat; no rat weights are",
        "reported because no rats were studied here (the rat arm is a",
        "simulation against digitised literature data). Time-invariant."
      ),
      source_name        = "BW"
    ),
    DOSE_PHAGE_PFU = list(
      description        = "Administered phage dose, in plaque-forming units",
      units              = "PFU",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Does NOT set the mass of the dose event -- the state system is in %ID,",
        "so every dose record carries amt = 100. DOSE_PHAGE_PFU enters ONLY the",
        "MPS carrying capacity, converting the literature capacity",
        "10^3.81 PFU per 100,000 phagocytes (Table 1 A_MPS; derived in Table S2)",
        "into the %ID units the state system uses:",
        "A_RES [%ID/1e5 cells] = 10^3.81 * 100 / DOSE_PHAGE_PFU.",
        "This is the transformation the paper describes in reverse under",
        "'Interspecies scaling': 'To simulate PFU per milliliter, doses were",
        "parameterized in PFU, and the UNTRANSFORMED macrophage capacity term",
        "was utilized (log10 PFU/100,000 cells)'. A dedicated DOSE_ column is",
        "required rather than the generic DOSE canonical because amt is in %ID",
        "while the capacity formula is calibrated in PFU; carrying both in one",
        "column would silently mix units (same rationale as DOSE_MTX_MGM2).",
        "Table S1 mean doses: 10^10.7 PFU (Luz24), 10^10.1 PFU (OMKO1),",
        "10^11.9 PFU (PAML-31-1)."
      ),
      source_name        = "Mean Phage Dose (Log10 PFU)"
    ),
    STUDY_PAML31 = list(
      description        = "1 = PAML-31-1 phage study arm; 0 = Luz24 or OMKO1 arm",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Luz24 and OMKO1 arms, which share one CL_Active)",
      notes              = paste(
        "Selects the covariate-adjusted active clearance CL_Active,PAML of",
        "Table 1 in place of the population CL_Active. The paper identified this",
        "as a categorical covariate after a single expectation step with",
        "inter-animal variability fixed to omega = 1, detecting a correlation",
        ">0.3 between CL_Active and 'experiments conducted using the PAML-31",
        "phage'; each phage was run as its own study (Table S1 is captioned",
        "'Phage dose and mouse weights for each phage study'), which is why this",
        "is registered on the STUDY_<id> family rather than as a phage-strain",
        "column. NOTE: the model applies the two Table 1 CL values DIRECTLY and",
        "never evaluates theta_PAML31 = -1.94, which reconciles with neither",
        "reported clearance under any standard transform -- see the vignette",
        "Errata."
      ),
      source_name        = "PAML-31"
    )
  )

  population <- list(
    species        = "rat",
    n_subjects     = NA_integer_,
    n_studies      = 1L,
    age_range      = NA_character_,
    weight_range   = NA_character_,
    weight_median  = "0.25 kg (the generic rat tabulated in Table S3)",
    sex_female_pct = NA_real_,
    disease_state  = "Healthy, uninfected rats (external-qualification simulation)",
    dose_range     = paste(
      "10^6, 10^9 and 10^11 PFU intravenously, administered as a rapid 5 min",
      "push (Methods 'External validation')."
    ),
    regions        = NA_character_,
    notes          = paste(
      "NO rat animals were studied here. The rat arm is a simulation: the",
      "estimated parameters come entirely from the CD-1 mouse fit and the rat",
      "enters only through the Table S3 anatomy and phagocyte densities. The",
      "comparison data are DIGITISED from the published Klebsiella pneumoniae",
      "phage study of Li et al. using WebPlotDigitizer (Methods 'External",
      "validation'), so no subject count, strain, age or sex split is",
      "recoverable. The mouse fit was naive-pooled with no random effects, so",
      "no IIV is encoded here either. Note the phage genus differs between the",
      "calibration (Pseudomonas) and the qualification (Klebsiella) data, which",
      "the authors flag as a deliberate stress test of extrapolation potential",
      "between phage genera."
    )
  )

  ini({
    # ==================================================================
    # Rat anatomy and physiology -- Supplemental Table S3, "Rat" column.
    # Every value is a literature-fixed physiologic constant, not an estimate.
    # ==================================================================
    # Table S3 labels CO as "L/h/kg". Read literally, a 70 kg human would have a
    # cardiac output of 19 L/min (3x physiologic) and a 20 g mouse 5.5 mL/min
    # (about a third of physiologic). Read as the allometric coefficient of
    # CO = CO_C * BW^0.75, all three tabulated species land on textbook values
    # simultaneously -- mouse 14.6 mL/min, rat 88 mL/min, human 6.66 L/min. The
    # allometric reading is adopted; see the vignette Errata.
    CO_C <- fixed(15) ; label("Cardiac output coefficient (L/h/kg^0.75)")                      # Table S3 CO, Rat

    # Fractional cardiac output by organ (unitless). The lung sits in series and
    # takes the whole cardiac output; the carcass takes the balance of the
    # eleven tabulated systemic organs, so no F_Q_Crs row is needed.
    fq_brain            <- fixed(0.02)    ; label("Fraction of cardiac output to brain (unitless)")            # Table S3 F Q_Brn, Rat
    fq_liver            <- fixed(0.021)    ; label("Fraction of cardiac output to liver, hepatic artery (unitless)")  # Table S3 F Q_Lvr, Rat
    fq_spleen           <- fixed(0.0085)    ; label("Fraction of cardiac output to spleen (unitless)")           # Table S3 F Q_Spn, Rat
    fq_small_intestine  <- fixed(0.104)    ; label("Fraction of cardiac output to small intestine (unitless)")  # Table S3 F Q_SI, Rat
    fq_large_intestine  <- fixed(0.036)   ; label("Fraction of cardiac output to large intestine (unitless)")  # Table S3 F Q_LI, Rat
    fq_kidney           <- fixed(0.141)    ; label("Fraction of cardiac output to kidney (unitless)")           # Table S3 F Q_Kid, Rat
    fq_muscle           <- fixed(0.278)    ; label("Fraction of cardiac output to muscle (unitless)")           # Table S3 F Q_Msc, Rat
    fq_skin             <- fixed(0.058)    ; label("Fraction of cardiac output to skin (unitless)")             # Table S3 F Q_Ski, Rat
    fq_bone             <- fixed(0.122)   ; label("Fraction of cardiac output to bone (unitless)")             # Table S3 F Q_Bon, Rat
    fq_stomach          <- fixed(0.013) ; label("Fraction of cardiac output to stomach (unitless)")          # Table S3 F Q_Sto, Rat

    # Organ volume as a fraction of body weight (unitless). Carcass closes the
    # balance against 1, so no F_V_Crs row is needed.
    fv_lung             <- fixed(0.005)   ; label("Lung volume as a fraction of body weight (unitless)")            # Table S3 F V_Lun, Rat
    fv_brain            <- fixed(0.00057)   ; label("Brain volume as a fraction of body weight (unitless)")           # Table S3 F V_Brn, Rat
    fv_liver            <- fixed(0.0366)   ; label("Liver volume as a fraction of body weight (unitless)")           # Table S3 F V_Lvr, Rat
    fv_spleen           <- fixed(0.002)    ; label("Spleen volume as a fraction of body weight (unitless)")          # Table S3 F V_Spn, Rat (row captioned "Fraction of BW as Lung" -- a copy-paste label error; the row is spleen, see Errata)
    fv_small_intestine  <- fixed(0.014)   ; label("Small-intestine volume as a fraction of body weight (unitless)") # Table S3 F V_SI, Rat
    fv_large_intestine  <- fixed(0.0084)   ; label("Large-intestine volume as a fraction of body weight (unitless)") # Table S3 F V_LI, Rat
    fv_kidney           <- fixed(0.0073)   ; label("Kidney volume as a fraction of body weight (unitless)")          # Table S3 F V_Kid, Rat
    fv_muscle           <- fixed(0.4043)    ; label("Muscle volume as a fraction of body weight (unitless)")          # Table S3 F V_Msc, Rat
    fv_skin             <- fixed(0.1903)   ; label("Skin volume as a fraction of body weight (unitless)")            # Table S3 F V_Ski, Rat
    fv_bone             <- fixed(0.02)     ; label("Bone volume as a fraction of body weight (unitless)")            # Table S3 F V_Mar, Rat
    fv_stomach          <- fixed(0.0046)    ; label("Stomach volume as a fraction of body weight (unitless)")         # Table S3 Fr_V_Sto, Rat
    fv_blood            <- fixed(0.074) ; label("Blood volume as a fraction of body weight (unitless)")           # Table S3 F V_Pla, Rat

    # Vascular (capillary blood) fraction of each organ's volume (unitless).
    bv_lung             <- fixed(0.36)   ; label("Vascular fraction of lung volume (unitless)")             # Table S3 F V_Lun_Bld, Rat
    bv_brain            <- fixed(0.03)  ; label("Vascular fraction of brain volume (unitless)")            # Table S3 F V_Brn_Bld, Rat
    bv_liver            <- fixed(0.21)  ; label("Vascular fraction of liver volume (unitless)")            # Table S3 F V_Lvr_Bld, Rat
    bv_spleen           <- fixed(0.22)  ; label("Vascular fraction of spleen volume (unitless)")           # Table S3 F V_Spn_Bld, Rat
    bv_small_intestine  <- fixed(0.3)  ; label("Vascular fraction of small-intestine volume (unitless)")  # Table S3 F V_SI_Bld, Rat
    bv_large_intestine  <- fixed(0.3)  ; label("Vascular fraction of large-intestine volume (unitless)")  # Table S3 F V_LI_Bld, Rat
    bv_kidney           <- fixed(0.16)  ; label("Vascular fraction of kidney volume (unitless)")           # Table S3 F V_Kid_Bld, Rat
    bv_muscle           <- fixed(0.04)  ; label("Vascular fraction of muscle volume (unitless)")           # Table S3 F V_Msc_Bld, Rat
    bv_skin             <- fixed(0.02)  ; label("Vascular fraction of skin volume (unitless)")             # Table S3 F V_Ski_Bld, Rat
    bv_remainder        <- fixed(0.04)  ; label("Vascular fraction of carcass volume (unitless)")          # Table S3 F V_Crs_Bld, Rat
    bv_bone             <- fixed(0.04)   ; label("Vascular fraction of bone volume (unitless)")             # Table S3 F V_Bon_Bld, Rat
    bv_stomach          <- fixed(0.3)  ; label("Vascular fraction of stomach volume (unitless)")          # Table S3 F V_Sto_Bld, Rat

    # Tissue-resident phagocyte density (cells/g). Sets each organ's MPS
    # carrying capacity.
    m_blood             <- fixed(14300)  ; label("Phagocyte density in blood (cells/g)")             # Table S3 M Pla, Rat
    m_lung              <- fixed(9.21e6) ; label("Phagocyte density in lung (cells/g)")              # Table S3 M Lun, Rat
    # Table S3 also tabulates M Brn (mouse 1.8e5 cells/g), but the supplement's
    # brain block prints no MPS uptake, release or degradation terms, so the
    # brain has no MPS pool in this model and its phagocyte density is unused.
    # It is deliberately not declared here rather than declared and ignored.
    m_liver             <- fixed(2.70e7) ; label("Phagocyte density in liver (cells/g)")             # Table S3 M Lvr, Rat
    m_spleen            <- fixed(2.28e8) ; label("Phagocyte density in spleen (cells/g)")            # Table S3 M Spn, Rat
    m_small_intestine   <- fixed(60000)  ; label("Phagocyte density in small intestine (cells/g)")   # Table S3 M SI, Rat
    m_large_intestine   <- fixed(15000)  ; label("Phagocyte density in large intestine (cells/g)")   # Table S3 M LI, Rat
    m_kidney            <- fixed(390000)  ; label("Phagocyte density in kidney (cells/g)")            # Table S3 M Kid, Rat
    m_remainder         <- fixed(6.35e6)  ; label("Phagocyte density in carcass (cells/g)")           # Table S3 M Crs, Rat
    m_muscle            <- fixed(40000)  ; label("Phagocyte density in muscle (cells/g)")            # Table S3 M Msc, Rat
    m_skin              <- fixed(40000)  ; label("Phagocyte density in skin (cells/g)")              # Table S3 M Ski, Rat
    m_bone              <- fixed(1.49e7)  ; label("Phagocyte density in bone (cells/g)")              # Table S3 M Bon, Rat
    m_stomach           <- fixed(6.35e6)  ; label("Phagocyte density in stomach (cells/g)")           # Table S3 M Sto, Rat

    # Venous share of total blood volume. NOT reported by Echterhof 2026; 0.8 is
    # the convention of the nanoparticle-PBPK framework lineage this model is
    # built on (Carlander 2016 / Deng 2019, Table S3 citations 6-7; the same 0.8
    # appears in Cheng 2020 and in Kutumova_2024_albuminNanoparticles_mouse_pbpk.R).
    # The printed venous equation is not mass-conserving and references an
    # arterial concentration C_a(t) that is never given an equation, so the
    # two-pool venous -> lung -> arterial -> organs -> venous topology is taken
    # from Figure 1. Changing the split, or collapsing to a single blood pool,
    # moves every reported output by less than 3 percent. See the Errata.
    f_venous <- fixed(0.8) ; label("Venous share of total blood volume (unitless)")  # Not reported by Echterhof 2026; framework convention, see Errata

    # ==================================================================
    # Estimated structural parameters -- Table 1.
    # The abstract and the Results narrative print DIFFERENT values for six of
    # these (K_p,LunKid 0.000138, K_p,LvrSpnSto 0.627, K_p,Crs 0.220, P_S
    # 0.0227, CL_Active 0.0145, T_MPS,deg 0.0301) and the abstract additionally
    # swaps which organ group gets the small partition coefficient. Table 1 is
    # the final set: it reproduces the paper's own reported blood terminal
    # half-life and liver mean residence time, and the narrative set does not.
    # See the vignette Errata for the discriminating simulation.
    # ==================================================================
    kp_lunkid     <- 0.395     ; label("Partition coefficient, high-flow organs: lung and kidney (unitless)")   # Table 1 K p,LunKid, 0.395 (RSE 5.15%, 0.358-0.437)
    kp_lvrspnsto  <- 2.12e-4   ; label("Partition coefficient, liver, spleen and stomach (unitless)")           # Table 1 K p,LvrSpnSto, 2.12e-4 (RSE 7.17%, 1.84e-4 to 2.44e-4)
    kp_crs        <- 0.934     ; label("Partition coefficient, rest of body (unitless)")                        # Table 1 K p,Crs, 0.934 (RSE 6.17%, 0.828-1.05)
    ps            <- 4.27e-2   ; label("Whole-body permeability-surface coefficient, as a fraction of organ blood flow (unitless)")  # Table 1 P s, 4.27e-2 (RSE 14.3%, 0.0324-0.0562)

    a_mps_log10   <- fixed(3.81)  ; label("MPS carrying capacity, log10(PFU per 100,000 phagocytes)")           # Table 1 A MPS, 3.81, Fixed from Bichet et al.; derived in Table S2 as log10(6.59e3) = 3.81913
    t_mps_up      <- fixed(0.001) ; label("MPS uptake half-life (h)")                                           # Table 1 T MPS,up, 0.001, Fixed for identifiability (endosomal degradation assumed rate-limiting)
    t_mps_deg     <- 4.37e-2      ; label("MPS degradation half-life (h)")                                      # Table 1 T MPS,deg, 4.37e-2 (RSE 12.2%, 3.45e-2 to 5.54e-2)

    clactive_c      <- 1.29e-2 ; label("Active surface clearance, Luz24 and OMKO1 (L/h/kg)")                    # Table 1 CL Active, 1.29e-2 (RSE 5.10%, 1.16e-2 to 1.42e-2)
    clactive_paml_c <- 1.34e-3 ; label("Active surface clearance, PAML-31-1 (L/h/kg)")                          # Table 1 CL Active,PAML (covariate-adjusted row), 1.34e-3 (RSE 0.205%, 0.0013-0.0014)

    t_sto_out     <- 3.28 ; label("Gastric transit half-life out of the gut lumen (h)")                         # Table 1 T Sto,out, 3.28 (RSE 29.8%, 1.89-5.68)
    t_u_out       <- 1.73 ; label("Urine transit half-life (h)")                                                # Table 1 T U,out, 1.73 (RSE 31.3%, 0.971-3.07)

    # ==================================================================
    # Residual error -- Table 1 "Error model parameters". Reported as
    # "residual variance" in %ID/g, but the values are on the scale of the
    # observations themselves (blood peaks near 15 %ID/g against sigma_Pla
    # 4.28), so they are additive standard deviations, which is also Monolix's
    # default constant-error parameterisation. No IIV: the fit was naive-pooled
    # with no random effects (Methods), so every sigma is the only stochastic
    # term in the model.
    # ==================================================================
    addSd                     <- 4.28  ; label("Additive residual SD, blood (%ID/g)")                    # Table 1 sigma Pla, 4.28E+00 (RSE 6.00%)
    addSd_Clung               <- 0.721 ; label("Additive residual SD, lung (%ID/g)")                      # Table 1 sigma Lun, 7.21E-01 (RSE 7.50%)
    addSd_Cliver              <- 1.23  ; label("Additive residual SD, liver (%ID/g)")                     # Table 1 sigma Lvr, 1.23E+00 (RSE 7.50%)
    addSd_Cspleen             <- 3.26  ; label("Additive residual SD, spleen (%ID/g)")                    # Table 1 sigma Spn, 3.26E+00 (RSE 7.50%)
    addSd_Csmall_intestine    <- 1.18  ; label("Additive residual SD, small intestine (%ID/g)")           # Table 1 sigma Smi, 1.18E+00 (RSE 7.50%)
    addSd_Clarge_intestine    <- 1.39  ; label("Additive residual SD, large intestine (%ID/g)")           # Table 1 sigma Lgi, 1.39E+00 (RSE 7.50%)
    addSd_Cmuscle             <- 0.902 ; label("Additive residual SD, muscle (%ID/g)")                    # Table 1 sigma Msc, 9.02E-01 (RSE 7.50%)
    addSd_Ckidney             <- 0.899 ; label("Additive residual SD, kidney (%ID/g)")                    # Table 1 sigma Kid, 8.99E-01 (RSE 7.50%)
    addSd_Cstomach            <- 4.15  ; label("Additive residual SD, stomach (%ID/g)")                   # Table 1 sigma Sto, 4.15E+00 (RSE 9.21%)
    addSd_Cbone               <- 1.14  ; label("Additive residual SD, bone (%ID/g)")                      # Table 1 sigma Bon, 1.14E+00 (RSE 7.50%)
    addSd_gut_lumen           <- 6.53  ; label("Additive residual SD, stomach contents (%ID)")            # Table 1 sigma Stc, 6.53E+00 (RSE 8.51%)
    addSd_urine               <- 6.12  ; label("Additive residual SD, urine (%ID)")                       # Table 1 sigma Urn, 6.12E+00 (RSE 9.53%)
  })

  model({
    # ================================================================
    # 1. Blood flows (mL/h). Volumes are handled in mL and amounts in %ID, so
    #    concentrations come out directly as %ID/g (tissue density 1 g/mL).
    # ================================================================
    q_co <- CO_C * WT^0.75 * 1000

    q_brain            <- fq_brain            * q_co
    q_liver_art        <- fq_liver            * q_co
    q_spleen           <- fq_spleen           * q_co
    q_small_intestine  <- fq_small_intestine  * q_co
    q_large_intestine  <- fq_large_intestine  * q_co
    q_kidney           <- fq_kidney           * q_co
    q_muscle           <- fq_muscle           * q_co
    q_skin             <- fq_skin             * q_co
    q_bone             <- fq_bone             * q_co
    q_stomach          <- fq_stomach          * q_co

    # Hepatic portal vein: spleen + stomach + small and large intestine drain
    # into the liver rather than back to venous blood (supplement, Liver block:
    # Q_Hpv = Q_Spn + Q_Sto + Q_SI + Q_LI).
    q_hpv <- q_spleen + q_stomach + q_small_intestine + q_large_intestine

    # The carcass takes the balance of cardiac output; the lung sits in series
    # and takes all of it.
    q_remainder <- q_co - q_brain - q_liver_art - q_kidney - q_muscle -
      q_skin - q_bone - q_hpv
    q_lung <- q_co

    # ================================================================
    # 2. Organ volumes (mL). Carcass closes the volume balance.
    # ================================================================
    v_lung             <- fv_lung             * WT * 1000
    v_brain            <- fv_brain            * WT * 1000
    v_liver            <- fv_liver            * WT * 1000
    v_spleen           <- fv_spleen           * WT * 1000
    v_small_intestine  <- fv_small_intestine  * WT * 1000
    v_large_intestine  <- fv_large_intestine  * WT * 1000
    v_kidney           <- fv_kidney           * WT * 1000
    v_muscle           <- fv_muscle           * WT * 1000
    v_skin             <- fv_skin             * WT * 1000
    v_bone             <- fv_bone             * WT * 1000
    v_stomach          <- fv_stomach          * WT * 1000
    v_blood            <- fv_blood            * WT * 1000

    v_remainder <- WT * 1000 - v_lung - v_brain - v_liver - v_spleen -
      v_small_intestine - v_large_intestine - v_kidney - v_muscle -
      v_skin - v_bone - v_stomach - v_blood

    # Vascular / interstitial split of each organ
    vv_lung            <- v_lung            * bv_lung
    vv_brain           <- v_brain           * bv_brain
    vv_liver           <- v_liver           * bv_liver
    vv_spleen          <- v_spleen          * bv_spleen
    vv_small_intestine <- v_small_intestine * bv_small_intestine
    vv_large_intestine <- v_large_intestine * bv_large_intestine
    vv_kidney          <- v_kidney          * bv_kidney
    vv_muscle          <- v_muscle          * bv_muscle
    vv_skin            <- v_skin            * bv_skin
    vv_bone            <- v_bone            * bv_bone
    vv_stomach         <- v_stomach         * bv_stomach
    vv_remainder       <- v_remainder       * bv_remainder

    vi_lung            <- v_lung            - vv_lung
    vi_brain           <- v_brain           - vv_brain
    vi_liver           <- v_liver           - vv_liver
    vi_spleen          <- v_spleen          - vv_spleen
    vi_small_intestine <- v_small_intestine - vv_small_intestine
    vi_large_intestine <- v_large_intestine - vv_large_intestine
    vi_kidney          <- v_kidney          - vv_kidney
    vi_muscle          <- v_muscle          - vv_muscle
    vi_skin            <- v_skin            - vv_skin
    vi_bone            <- v_bone            - vv_bone
    vi_stomach         <- v_stomach         - vv_stomach
    vi_remainder       <- v_remainder       - vv_remainder

    v_venous   <- f_venous * v_blood
    v_arterial <- v_blood - v_venous

    # ================================================================
    # 3. Rate constants. Every "T" parameter in Table 1 is a HALF-life
    #    (Table 1 definitions: "Uptake half-life", "Gastric transit half-life",
    #    "Urine transit half-life"; Results calls T_MPS,deg "a degradation
    #    half-life"), so k = ln(2) / T throughout.
    # ================================================================
    k_up_max <- log(2) / t_mps_up
    # Supplement, note under the Venous equation: "in all cases k_rel = k_up,max,
    # due to limited ability to differentiate the difference between
    # degradation-limited and cell entry-limited kinetics."
    k_rel        <- k_up_max
    k_deg        <- log(2) / t_mps_deg
    k_gi_out     <- log(2) / t_sto_out
    k_urine_out  <- log(2) / t_u_out

    # Active surface clearance (mL/h); Table 1 unit is L/h/kg.
    clactive <- (clactive_c * (1 - STUDY_PAML31) +
                   clactive_paml_c * STUDY_PAML31) * WT * 1000

    # ================================================================
    # 4. MPS carrying capacity (%ID), supplement A_i,max = M_i * V_i * A_RES/1e5.
    #    A_RES is tabulated in PFU per 100,000 phagocytes and is converted to
    #    the %ID scale of the state system by the administered dose.
    # ================================================================
    a_res_unit <- 10^a_mps_log10 * 100 / DOSE_PHAGE_PFU

    amax_venous            <- m_blood            * v_venous           * a_res_unit / 1e5
    amax_lung              <- m_lung             * v_lung             * a_res_unit / 1e5
    amax_liver             <- m_liver            * v_liver            * a_res_unit / 1e5
    amax_spleen            <- m_spleen           * v_spleen           * a_res_unit / 1e5
    amax_small_intestine   <- m_small_intestine  * v_small_intestine  * a_res_unit / 1e5
    amax_large_intestine   <- m_large_intestine  * v_large_intestine  * a_res_unit / 1e5
    amax_kidney            <- m_kidney           * v_kidney           * a_res_unit / 1e5
    amax_muscle            <- m_muscle           * v_muscle           * a_res_unit / 1e5
    amax_skin              <- m_skin             * v_skin             * a_res_unit / 1e5
    amax_bone              <- m_bone             * v_bone             * a_res_unit / 1e5
    amax_stomach           <- m_stomach          * v_stomach          * a_res_unit / 1e5
    amax_remainder         <- m_remainder        * v_remainder        * a_res_unit / 1e5

    # ================================================================
    # 5. Concentrations (%ID/g)
    # ================================================================
    c_venous   <- venous   / v_venous
    c_arterial <- arterial / v_arterial

    cv_lung            <- vp_lung            / vv_lung
    cv_brain           <- vp_brain           / vv_brain
    cv_liver           <- vp_liver           / vv_liver
    cv_spleen          <- vp_spleen          / vv_spleen
    cv_small_intestine <- vp_small_intestine / vv_small_intestine
    cv_large_intestine <- vp_large_intestine / vv_large_intestine
    cv_kidney          <- vp_kidney          / vv_kidney
    cv_muscle          <- vp_muscle          / vv_muscle
    cv_skin            <- vp_skin            / vv_skin
    cv_bone            <- vp_bone            / vv_bone
    cv_stomach         <- vp_stomach         / vv_stomach
    cv_remainder       <- vp_remainder       / vv_remainder

    ci_lung            <- is_lung            / vi_lung
    ci_brain           <- is_brain           / vi_brain
    ci_liver           <- is_liver           / vi_liver
    ci_spleen          <- is_spleen          / vi_spleen
    ci_small_intestine <- is_small_intestine / vi_small_intestine
    ci_large_intestine <- is_large_intestine / vi_large_intestine
    ci_kidney          <- is_kidney          / vi_kidney
    ci_muscle          <- is_muscle          / vi_muscle
    ci_skin            <- is_skin            / vi_skin
    ci_bone            <- is_bone            / vi_bone
    ci_stomach         <- is_stomach         / vi_stomach
    ci_remainder       <- is_remainder       / vi_remainder

    # ================================================================
    # 6. Transcapillary diffusion, P_S * Q_i * (C_i,V - C_i,int / K_p)
    # ================================================================
    diff_lung            <- ps * q_lung            * (cv_lung            - ci_lung            / kp_lunkid)
    diff_brain           <- ps * q_brain           * (cv_brain           - ci_brain           / kp_crs)
    diff_liver           <- ps * q_liver_art       * (cv_liver           - ci_liver           / kp_lvrspnsto)
    diff_spleen          <- ps * q_spleen          * (cv_spleen          - ci_spleen          / kp_lvrspnsto)
    diff_small_intestine <- ps * q_small_intestine * (cv_small_intestine - ci_small_intestine / kp_crs)
    diff_large_intestine <- ps * q_large_intestine * (cv_large_intestine - ci_large_intestine / kp_crs)
    diff_kidney          <- ps * q_kidney          * (cv_kidney          - ci_kidney          / kp_lunkid)
    diff_muscle          <- ps * q_muscle          * (cv_muscle          - ci_muscle          / kp_crs)
    diff_skin            <- ps * q_skin            * (cv_skin            - ci_skin            / kp_crs)
    diff_bone            <- ps * q_bone            * (cv_bone            - ci_bone            / kp_crs)
    diff_stomach         <- ps * q_stomach         * (cv_stomach         - ci_stomach         / kp_lvrspnsto)
    diff_remainder       <- ps * q_remainder       * (cv_remainder       - ci_remainder       / kp_crs)

    # ================================================================
    # 7. Saturable MPS uptake, k_i,up(t) = k_up,max * (1 - A_i,res / A_i,max),
    #    exactly as printed. Do NOT wrap this in max(0, .): the relation is
    #    self-correcting (an overshoot above capacity simply reverses the flux,
    #    and the pools never in fact go negative), whereas the clamp puts a
    #    derivative kink precisely at the operating point A_i,res ~ A_i,max,
    #    which stalls the implicit solver -- the PAML-31-1 arm fails to
    #    integrate at all with the clamp in place and solves instantly without.
    # ================================================================
    up_venous            <- k_up_max * (1 - int_venous            / amax_venous)            * venous
    up_lung              <- k_up_max * (1 - int_lung              / amax_lung)              * is_lung
    up_liver             <- k_up_max * (1 - int_liver             / amax_liver)             * is_liver
    up_spleen            <- k_up_max * (1 - int_spleen            / amax_spleen)            * is_spleen
    up_small_intestine   <- k_up_max * (1 - int_small_intestine   / amax_small_intestine)   * is_small_intestine
    up_large_intestine   <- k_up_max * (1 - int_large_intestine   / amax_large_intestine)   * is_large_intestine
    up_kidney            <- k_up_max * (1 - int_kidney            / amax_kidney)            * is_kidney
    up_muscle            <- k_up_max * (1 - int_muscle            / amax_muscle)            * is_muscle
    up_skin              <- k_up_max * (1 - int_skin              / amax_skin)              * is_skin
    up_bone              <- k_up_max * (1 - int_bone              / amax_bone)              * is_bone
    up_stomach           <- k_up_max * (1 - int_stomach           / amax_stomach)           * is_stomach
    up_remainder         <- k_up_max * (1 - int_remainder         / amax_remainder)         * is_remainder

    # ================================================================
    # 8. ODE system (40 states, all amounts in %ID)
    # ================================================================
    # --- Blood. The intravenous bolus lands in venous blood. The printed
    #     venous equation is not mass-conserving (it debits Q_Lun * C_Lun,v
    #     rather than Q_Lun * C_V) and the printed organ equations are driven by
    #     an arterial concentration C_a(t) that the supplement never defines, so
    #     the venous -> lung -> arterial -> organs -> venous topology is taken
    #     from Figure 1. See the vignette Errata.
    d/dt(venous) <- q_brain * cv_brain +
      (q_liver_art + q_hpv) * cv_liver +
      q_kidney * cv_kidney + q_muscle * cv_muscle + q_skin * cv_skin +
      q_bone * cv_bone + q_remainder * cv_remainder -
      q_lung * c_venous - up_venous + k_rel * int_venous
    d/dt(int_venous) <- up_venous - k_rel * int_venous - k_deg * int_venous

    d/dt(arterial) <- q_lung * cv_lung - q_co * c_arterial

    # --- Lung: in series, fed by venous blood, draining to arterial.
    d/dt(vp_lung)  <- q_lung * c_venous - q_lung * cv_lung - diff_lung
    d/dt(is_lung)  <- diff_lung - up_lung + k_rel * int_lung
    d/dt(int_lung) <- up_lung - k_rel * int_lung - k_deg * int_lung

    # --- Liver: dual input from the hepatic artery and the portal vein, plus
    #     active surface (hepatobiliary) clearance.
    d/dt(vp_liver) <- q_liver_art * c_arterial +
      q_spleen * cv_spleen + q_stomach * cv_stomach +
      q_small_intestine * cv_small_intestine +
      q_large_intestine * cv_large_intestine -
      (q_liver_art + q_hpv) * cv_liver - diff_liver - clactive * cv_liver
    d/dt(is_liver)  <- diff_liver - up_liver + k_rel * int_liver
    d/dt(int_liver) <- up_liver - k_rel * int_liver - k_deg * int_liver

    # --- Stomach, small and large intestine: drain to the liver by the portal
    #     vein and carry active surface clearance.
    d/dt(vp_stomach)  <- q_stomach * (c_arterial - cv_stomach) - diff_stomach -
      clactive * cv_stomach
    d/dt(is_stomach)  <- diff_stomach - up_stomach + k_rel * int_stomach
    d/dt(int_stomach) <- up_stomach - k_rel * int_stomach - k_deg * int_stomach

    d/dt(vp_small_intestine) <- q_small_intestine * (c_arterial - cv_small_intestine) -
      diff_small_intestine - clactive * cv_small_intestine
    d/dt(is_small_intestine) <- diff_small_intestine - up_small_intestine +
      k_rel * int_small_intestine
    d/dt(int_small_intestine) <- up_small_intestine -
      k_rel * int_small_intestine - k_deg * int_small_intestine

    d/dt(vp_large_intestine) <- q_large_intestine * (c_arterial - cv_large_intestine) -
      diff_large_intestine - clactive * cv_large_intestine
    d/dt(is_large_intestine) <- diff_large_intestine - up_large_intestine +
      k_rel * int_large_intestine
    d/dt(int_large_intestine) <- up_large_intestine -
      k_rel * int_large_intestine - k_deg * int_large_intestine

    # --- Spleen: drains to the liver, no active clearance.
    d/dt(vp_spleen)  <- q_spleen * (c_arterial - cv_spleen) - diff_spleen
    d/dt(is_spleen)  <- diff_spleen - up_spleen + k_rel * int_spleen
    d/dt(int_spleen) <- up_spleen - k_rel * int_spleen - k_deg * int_spleen

    # --- Kidney: active clearance drains to urine.
    d/dt(vp_kidney)  <- q_kidney * (c_arterial - cv_kidney) - diff_kidney -
      clactive * cv_kidney
    d/dt(is_kidney)  <- diff_kidney - up_kidney + k_rel * int_kidney
    d/dt(int_kidney) <- up_kidney - k_rel * int_kidney - k_deg * int_kidney

    # --- Brain: the supplement's brain block prints ONLY the vascular and
    #     interstitial equations, with no MPS uptake, release or degradation
    #     terms, even though Figure 1 draws an MPS box for the brain. The
    #     printed equations are followed. See the vignette Errata.
    d/dt(vp_brain) <- q_brain * (c_arterial - cv_brain) - diff_brain
    d/dt(is_brain) <- diff_brain

    # --- Muscle, skin, bone and the carcass remainder.
    d/dt(vp_muscle)  <- q_muscle * (c_arterial - cv_muscle) - diff_muscle
    d/dt(is_muscle)  <- diff_muscle - up_muscle + k_rel * int_muscle
    d/dt(int_muscle) <- up_muscle - k_rel * int_muscle - k_deg * int_muscle

    d/dt(vp_skin)  <- q_skin * (c_arterial - cv_skin) - diff_skin
    d/dt(is_skin)  <- diff_skin - up_skin + k_rel * int_skin
    d/dt(int_skin) <- up_skin - k_rel * int_skin - k_deg * int_skin

    d/dt(vp_bone)  <- q_bone * (c_arterial - cv_bone) - diff_bone
    d/dt(is_bone)  <- diff_bone - up_bone + k_rel * int_bone
    d/dt(int_bone) <- up_bone - k_rel * int_bone - k_deg * int_bone

    d/dt(vp_remainder)  <- q_remainder * (c_arterial - cv_remainder) - diff_remainder
    d/dt(is_remainder)  <- diff_remainder - up_remainder + k_rel * int_remainder
    d/dt(int_remainder) <- up_remainder - k_rel * int_remainder - k_deg * int_remainder

    # --- Excretion. The urine equation is printed in the supplement's Kidney
    #     block; the GI-elimination equation is NOT printed anywhere and is
    #     reconstructed by analogy with it, collecting the active clearance of
    #     the four organs Figure 1 draws feeding the "GI Elimination" box and
    #     draining at the gastric transit half-life. See the vignette Errata.
    d/dt(urine) <- clactive * cv_kidney - k_urine_out * urine
    d/dt(gut_lumen) <- clactive * (cv_liver + cv_stomach +
                                     cv_small_intestine + cv_large_intestine) -
      k_gi_out * gut_lumen

    # ================================================================
    # 9. Observations. Organs were excised whole, washed and counted, so the
    #    measured %ID/g is the total organ content -- vascular plus
    #    interstitial plus the MPS pool -- over the whole organ weight.
    #    Stomach contents and urine were reported as %ID, i.e. as amounts.
    # ================================================================
    Cc                 <- c_venous
    Clung              <- (vp_lung            + is_lung            + int_lung)            / v_lung
    Cliver             <- (vp_liver           + is_liver           + int_liver)           / v_liver
    Cspleen            <- (vp_spleen          + is_spleen          + int_spleen)          / v_spleen
    Csmall_intestine   <- (vp_small_intestine + is_small_intestine + int_small_intestine) / v_small_intestine
    Clarge_intestine   <- (vp_large_intestine + is_large_intestine + int_large_intestine) / v_large_intestine
    Cmuscle            <- (vp_muscle          + is_muscle          + int_muscle)          / v_muscle
    Ckidney            <- (vp_kidney          + is_kidney          + int_kidney)          / v_kidney
    Cstomach           <- (vp_stomach         + is_stomach         + int_stomach)         / v_stomach
    Cbone              <- (vp_bone            + is_bone            + int_bone)            / v_bone
    # Hypothesis-generating organs -- no data were available to estimate their
    # parameters (Methods: skin and brain "were included for
    # hypothesis-generating purposes"), so they carry no residual-error model.
    Cskin              <- (vp_skin      + is_skin      + int_skin)      / v_skin
    Cbrain             <- (vp_brain     + is_brain)                     / v_brain
    Cremainder         <- (vp_remainder + is_remainder + int_remainder) / v_remainder

    # Total 125I-phage in the system, for the vignette mass-balance check: set
    # clactive and t_mps_deg elimination to zero and this must hold at 100 %ID.
    amt_total <- venous + arterial + int_venous +
      vp_lung + is_lung + int_lung +
      vp_liver + is_liver + int_liver +
      vp_spleen + is_spleen + int_spleen +
      vp_small_intestine + is_small_intestine + int_small_intestine +
      vp_large_intestine + is_large_intestine + int_large_intestine +
      vp_kidney + is_kidney + int_kidney +
      vp_muscle + is_muscle + int_muscle +
      vp_skin + is_skin + int_skin +
      vp_bone + is_bone + int_bone +
      vp_stomach + is_stomach + int_stomach +
      vp_brain + is_brain +
      vp_remainder + is_remainder + int_remainder +
      urine + gut_lumen

    Cc                ~ add(addSd)
    Clung             ~ add(addSd_Clung)
    Cliver            ~ add(addSd_Cliver)
    Cspleen           ~ add(addSd_Cspleen)
    Csmall_intestine  ~ add(addSd_Csmall_intestine)
    Clarge_intestine  ~ add(addSd_Clarge_intestine)
    Cmuscle           ~ add(addSd_Cmuscle)
    Ckidney           ~ add(addSd_Ckidney)
    Cstomach          ~ add(addSd_Cstomach)
    Cbone             ~ add(addSd_Cbone)
    gut_lumen         ~ add(addSd_gut_lumen)
    urine             ~ add(addSd_urine)
  })
}
