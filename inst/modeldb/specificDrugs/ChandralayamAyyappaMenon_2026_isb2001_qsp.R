ChandralayamAyyappaMenon_2026_isb2001_qsp <- function() {
  description <- "QSP. Human quantitative systems pharmacology model of ISB 2001, a CD38 x BCMA x CD3 trispecific T cell engager, used to select the first-in-human dose and predict the clinical efficacy dose range in relapsed/refractory multiple myeloma (TRIgnite-1, NCT05862012). A minimal PBPK disposition model (plasma, tight tissue, lymph, leaky tissue split into bone marrow and non-marrow leaky tissue) is coupled to a full sequential target-engagement network in plasma and bone marrow: mass-action binding to CD3 on CD4+ and CD8+ T cells, membrane BCMA and CD38 on plasma cells, CD38 on monocytes/NK cells/neutrophils, and soluble BCMA and CD38, forming dimers, trimers and CD38-BCMA tetramers. The trimers and tetramers that crosslink a T cell to a tumour cell are summed as active species (ACT) and normalised per tumour cell (nACT), the exposure metric the paper used for dose prediction and clinical validation."
  reference <- paste(
    "Chandralayam Ayyappa Menon V, Matsuura T, Holkova B, Gudi GS, Drake A,",
    "Pihlgren M, van der Graaf PH, Sunitha GN, Garton A, Perro M, Konto C, Pacaud L.",
    "Clinical validation of a QSP model for ISB 2001, a trispecific T cell engager",
    "to support optimal FIH study design in RRMM patients.",
    "Clin Pharmacol Ther. 2026;120(2):452-464. doi:10.1002/cpt.70319.",
    "The minimal-PBPK backbone follows Cao & Jusko (2014) / Shah & Betts (2012)",
    "(references 15-17 of the source); the CD3-bispecific target-engagement and",
    "trimer framework follows Betts et al. (2019) AAPS J 21:66 (reference 8 of the",
    "source), see modellib('Betts_2019_pf_06671008_qsp').",
    sep = " "
  )
  vignette <- "ChandralayamAyyappaMenon_2026_isb2001"

  # Time is DAYS throughout, because every mPBPK flow (L, L1, L2), the clearance
  # and the SC absorption rate constant are reported per day (Table S2, Table 1).
  # The binding constants are reported per SECOND and the internalisation rate
  # constants per HOUR; both are kept in ini() exactly as printed and converted
  # to a per-day basis inside model(), so each ini() value can be read straight
  # off the source table.
  #
  # Every drug state holds a CONCENTRATION in nM, exactly as in the Simbiology
  # listing of Data S1, EXCEPT `depot`, which holds the subcutaneous dose as an
  # AMOUNT in nmol (Data S1 `Dose_sc`; its flux enters plasma as Ka*BA*Dose_sc
  # divided by the plasma volume). Doses are therefore given in nmol.
  #
  # The source never reports the molecular weight of ISB 2001, so a ug/kg dose
  # cannot be converted to nmol from on-disk information alone. The vignette
  # derives an apparent MW of ~1.9e5 g/mol from the paper's own paired reporting
  # of the same in-vitro potencies in nM (Table 2) and ng/mL (Results, "MABEL
  # dose") and states that derivation explicitly; it is NOT used here.
  units <- list(time = "day", dosing = "nmol", concentration = "nM")

  # The five mPBPK compartments use the Cao 2013 mAb mPBPK naming exception
  # (`plasma` / `tight` / `leaky` / `lymph`), `bonemarrow` follows
  # Poels_2025_elranatamab_qsp, and `depot` is the canonical SC dosing state.
  # Every remaining state is a drug-target complex of the sequential
  # target-engagement network and is declared paper-specific below: these are
  # combinatorial binding species of one particular trispecific molecule and do
  # not generalise to other models.
  paper_specific_compartments <- c(
    "drug_cd3cd4", "drug_cd3cd8", "drug_bcma", "drug_cd38", "drug_cd38_other",
    "trimer_cd3cd4_bcma", "trimer_cd3cd8_bcma", "trimer_cd3cd4_cd38",
    "trimer_cd3cd8_cd38", "trimer_bcma_cd38", "trimer_cd3cd4_cd38other",
    "trimer_cd3cd8_cd38other", "tetramer_cd3cd4", "tetramer_cd3cd8",
    "complex_sbcma", "complex_scd38", "complex_cd3cd4_sbcma",
    "complex_cd3cd8_sbcma", "complex_sbcma_cd38", "complex_cd3cd4_scd38",
    "complex_cd3cd8_scd38", "complex_sbcma_scd38",
    "complex_cd3cd4_sbcma_cd38", "complex_cd3cd8_sbcma_cd38",
    "drug_cd3cd4_bonemarrow", "drug_cd3cd8_bonemarrow", "drug_bcma_bonemarrow",
    "drug_cd38_bonemarrow", "trimer_cd3cd4_bcma_bonemarrow",
    "trimer_cd3cd8_bcma_bonemarrow", "trimer_cd3cd4_cd38_bonemarrow",
    "trimer_cd3cd8_cd38_bonemarrow", "trimer_bcma_cd38_bonemarrow",
    "tetramer_cd3cd4_bonemarrow", "tetramer_cd3cd8_bonemarrow",
    "complex_sbcma_bonemarrow", "complex_scd38_bonemarrow",
    "complex_cd3cd4_sbcma_bonemarrow", "complex_cd3cd8_sbcma_bonemarrow",
    "complex_sbcma_cd38_bonemarrow", "complex_cd3cd4_scd38_bonemarrow",
    "complex_cd3cd8_scd38_bonemarrow", "complex_sbcma_scd38_bonemarrow",
    "complex_cd3cd4_sbcma_cd38_bonemarrow", "complex_cd3cd8_sbcma_cd38_bonemarrow"
  )

  # Issue #482 state map. The Data S1 symbol for each state is given in the
  # analyte string so the ODE listing can be audited line by line. `verified`
  # is TRUE because every entry was checked against the "Biological species in
  # the model" table of Data S1.
  compartmentData <- list(
    depot       = list(analyte = "ISB 2001 (Data S1 Dose_sc)", units = "nmol", specimen = "administration site", verified = TRUE),
    plasma      = list(analyte = "ISB 2001, free (Data S1 Drug)", units = "nM", specimen = "plasma", verified = TRUE),
    tight       = list(analyte = "ISB 2001, free (Data S1 Drug_t)", units = "nM", specimen = "tissue", verified = TRUE),
    leaky       = list(analyte = "ISB 2001, free (Data S1 Drug_L)", units = "nM", specimen = "tissue", verified = TRUE),
    lymph       = list(analyte = "ISB 2001, free (Data S1 Drug_ly)", units = "nM", specimen = "lymph", verified = TRUE),
    bonemarrow  = list(analyte = "ISB 2001, free (Data S1 Drug_bm)", units = "nM", specimen = "tumor", verified = TRUE),

    drug_cd3cd4     = list(analyte = "ISB 2001-CD3(CD4+ T cell) dimer (Data S1 CX1)", units = "nM", specimen = "plasma", verified = TRUE),
    drug_cd3cd8     = list(analyte = "ISB 2001-CD3(CD8+ T cell) dimer (Data S1 CX2)", units = "nM", specimen = "plasma", verified = TRUE),
    drug_bcma       = list(analyte = "ISB 2001-membrane BCMA dimer (Data S1 CX3)", units = "nM", specimen = "plasma", verified = TRUE),
    drug_cd38       = list(analyte = "ISB 2001-membrane CD38 dimer, plasma cell (Data S1 CX4)", units = "nM", specimen = "plasma", verified = TRUE),
    drug_cd38_other = list(analyte = "ISB 2001-CD38 dimer on monocyte/NK/neutrophil (Data S1 CX5)", units = "nM", specimen = "plasma", verified = TRUE),

    trimer_cd3cd4_bcma      = list(analyte = "ISB 2001-CD3(CD4+)-BCMA trimer (Data S1 TX1)", units = "nM", specimen = "plasma", verified = TRUE),
    trimer_cd3cd8_bcma      = list(analyte = "ISB 2001-CD3(CD8+)-BCMA trimer (Data S1 TX2)", units = "nM", specimen = "plasma", verified = TRUE),
    trimer_cd3cd4_cd38      = list(analyte = "ISB 2001-CD3(CD4+)-CD38 trimer (Data S1 TX3)", units = "nM", specimen = "plasma", verified = TRUE),
    trimer_cd3cd8_cd38      = list(analyte = "ISB 2001-CD3(CD8+)-CD38 trimer (Data S1 TX4)", units = "nM", specimen = "plasma", verified = TRUE),
    trimer_bcma_cd38        = list(analyte = "ISB 2001-BCMA-CD38 trimer, both on tumour cell (Data S1 TX5)", units = "nM", specimen = "plasma", verified = TRUE),
    trimer_cd3cd4_cd38other = list(analyte = "ISB 2001-CD3(CD4+)-CD38 trimer, off-tumour CD38 (Data S1 TX6)", units = "nM", specimen = "plasma", verified = TRUE),
    trimer_cd3cd8_cd38other = list(analyte = "ISB 2001-CD3(CD8+)-CD38 trimer, off-tumour CD38 (Data S1 TX7)", units = "nM", specimen = "plasma", verified = TRUE),
    tetramer_cd3cd4         = list(analyte = "ISB 2001-CD3(CD4+)-BCMA-CD38 tetramer (Data S1 TeX1)", units = "nM", specimen = "plasma", verified = TRUE),
    tetramer_cd3cd8         = list(analyte = "ISB 2001-CD3(CD8+)-BCMA-CD38 tetramer (Data S1 TeX2)", units = "nM", specimen = "plasma", verified = TRUE),

    complex_sbcma             = list(analyte = "ISB 2001-soluble BCMA dimer (Data S1 sCX1)", units = "nM", specimen = "plasma", verified = TRUE),
    complex_scd38             = list(analyte = "ISB 2001-soluble CD38 dimer (Data S1 sCX2)", units = "nM", specimen = "plasma", verified = TRUE),
    complex_cd3cd4_sbcma      = list(analyte = "ISB 2001-CD3(CD4+)-soluble BCMA (Data S1 sTX1)", units = "nM", specimen = "plasma", verified = TRUE),
    complex_cd3cd8_sbcma      = list(analyte = "ISB 2001-CD3(CD8+)-soluble BCMA (Data S1 sTX2)", units = "nM", specimen = "plasma", verified = TRUE),
    complex_sbcma_cd38        = list(analyte = "ISB 2001-soluble BCMA-membrane CD38 (Data S1 sTX3)", units = "nM", specimen = "plasma", verified = TRUE),
    complex_cd3cd4_scd38      = list(analyte = "ISB 2001-CD3(CD4+)-soluble CD38 (Data S1 sTX4)", units = "nM", specimen = "plasma", verified = TRUE),
    complex_cd3cd8_scd38      = list(analyte = "ISB 2001-CD3(CD8+)-soluble CD38 (Data S1 sTX5)", units = "nM", specimen = "plasma", verified = TRUE),
    complex_sbcma_scd38       = list(analyte = "ISB 2001-soluble BCMA-soluble CD38 (Data S1 sTX6)", units = "nM", specimen = "plasma", verified = TRUE),
    complex_cd3cd4_sbcma_cd38 = list(analyte = "ISB 2001-CD3(CD4+)-soluble BCMA-membrane CD38 (Data S1 sTeX1)", units = "nM", specimen = "plasma", verified = TRUE),
    complex_cd3cd8_sbcma_cd38 = list(analyte = "ISB 2001-CD3(CD8+)-soluble BCMA-membrane CD38 (Data S1 sTeX2)", units = "nM", specimen = "plasma", verified = TRUE),

    drug_cd3cd4_bonemarrow = list(analyte = "ISB 2001-CD3(CD4+ T cell) dimer (Data S1 CX1_bm)", units = "nM", specimen = "tumor", verified = TRUE),
    drug_cd3cd8_bonemarrow = list(analyte = "ISB 2001-CD3(CD8+ T cell) dimer (Data S1 CX2_bm)", units = "nM", specimen = "tumor", verified = TRUE),
    drug_bcma_bonemarrow   = list(analyte = "ISB 2001-membrane BCMA dimer (Data S1 CX3_bm)", units = "nM", specimen = "tumor", verified = TRUE),
    drug_cd38_bonemarrow   = list(analyte = "ISB 2001-membrane CD38 dimer (Data S1 CX4_bm)", units = "nM", specimen = "tumor", verified = TRUE),

    trimer_cd3cd4_bcma_bonemarrow = list(analyte = "ISB 2001-CD3(CD4+)-BCMA trimer (Data S1 TX1_bm)", units = "nM", specimen = "tumor", verified = TRUE),
    trimer_cd3cd8_bcma_bonemarrow = list(analyte = "ISB 2001-CD3(CD8+)-BCMA trimer (Data S1 TX2_bm)", units = "nM", specimen = "tumor", verified = TRUE),
    trimer_cd3cd4_cd38_bonemarrow = list(analyte = "ISB 2001-CD3(CD4+)-CD38 trimer (Data S1 TX3_bm)", units = "nM", specimen = "tumor", verified = TRUE),
    trimer_cd3cd8_cd38_bonemarrow = list(analyte = "ISB 2001-CD3(CD8+)-CD38 trimer (Data S1 TX4_bm)", units = "nM", specimen = "tumor", verified = TRUE),
    trimer_bcma_cd38_bonemarrow   = list(analyte = "ISB 2001-BCMA-CD38 trimer, both on tumour cell (Data S1 TX5_bm)", units = "nM", specimen = "tumor", verified = TRUE),
    tetramer_cd3cd4_bonemarrow    = list(analyte = "ISB 2001-CD3(CD4+)-BCMA-CD38 tetramer (Data S1 TeX1_bm)", units = "nM", specimen = "tumor", verified = TRUE),
    tetramer_cd3cd8_bonemarrow    = list(analyte = "ISB 2001-CD3(CD8+)-BCMA-CD38 tetramer (Data S1 TeX2_bm)", units = "nM", specimen = "tumor", verified = TRUE),

    complex_sbcma_bonemarrow             = list(analyte = "ISB 2001-soluble BCMA dimer (Data S1 sCX1_bm)", units = "nM", specimen = "tumor", verified = TRUE),
    complex_scd38_bonemarrow             = list(analyte = "ISB 2001-soluble CD38 dimer (Data S1 sCX2_bm)", units = "nM", specimen = "tumor", verified = TRUE),
    complex_cd3cd4_sbcma_bonemarrow      = list(analyte = "ISB 2001-CD3(CD4+)-soluble BCMA (Data S1 sTX1_bm)", units = "nM", specimen = "tumor", verified = TRUE),
    complex_cd3cd8_sbcma_bonemarrow      = list(analyte = "ISB 2001-CD3(CD8+)-soluble BCMA (Data S1 sTX2_bm)", units = "nM", specimen = "tumor", verified = TRUE),
    complex_sbcma_cd38_bonemarrow        = list(analyte = "ISB 2001-soluble BCMA-membrane CD38 (Data S1 sTX3_bm)", units = "nM", specimen = "tumor", verified = TRUE),
    complex_cd3cd4_scd38_bonemarrow      = list(analyte = "ISB 2001-CD3(CD4+)-soluble CD38 (Data S1 sTX4_bm)", units = "nM", specimen = "tumor", verified = TRUE),
    complex_cd3cd8_scd38_bonemarrow      = list(analyte = "ISB 2001-CD3(CD8+)-soluble CD38 (Data S1 sTX5_bm)", units = "nM", specimen = "tumor", verified = TRUE),
    complex_sbcma_scd38_bonemarrow       = list(analyte = "ISB 2001-soluble BCMA-soluble CD38 (Data S1 sTX6_bm)", units = "nM", specimen = "tumor", verified = TRUE),
    complex_cd3cd4_sbcma_cd38_bonemarrow = list(analyte = "ISB 2001-CD3(CD4+)-soluble BCMA-membrane CD38 (Data S1 sTeX1_bm)", units = "nM", specimen = "tumor", verified = TRUE),
    complex_cd3cd8_sbcma_cd38_bonemarrow = list(analyte = "ISB 2001-CD3(CD8+)-soluble BCMA-membrane CD38 (Data S1 sTeX2_bm)", units = "nM", specimen = "tumor", verified = TRUE)
  )

  # The published QSP model carries no covariate columns: it is a deterministic
  # typical-patient model. The virtual-patient simulation of the source
  # ("Virtual patient simulation", Methods) varies body weight, SC
  # bioavailability and the sBCMA / CD38 / CD3 densities directly, which the
  # vignette does by perturbing the corresponding ini() parameters rather than
  # by carrying data columns.

  population <- list(
    species        = "human",
    n_subjects     = 30,
    n_studies      = 1,
    disease_state  = "Relapsed and/or refractory multiple myeloma (RRMM)",
    dose_range     = paste(
      "Subcutaneous. Fractionated step-up in cycle 1: priming dose on C1D1 (10% of the target dose),",
      "step-up dose on C1D4 (30% of the target dose), then the target dose weekly from C1D8.",
      "Dose levels DL1-DL8; C1D1 doses 0.5-15 ug/kg, with the C1D1 dose fixed at 15 ug/kg from DL4 onwards.",
      "Selected FIH target dose 5.0 ug/kg (DL1); first strong clinical responses at the 50 ug/kg target dose (DL3).",
      sep = " "
    ),
    regions        = "United States, Australia, India",
    notes          = paste(
      "Clinical validation cohort is the 30 patients dosed at DL1-DL8 in the phase 1 TRIgnite-1 study",
      "(NCT05862012) as of 10 February 2025 (Table 3, Figure S7A); DL1-DL3 were single-patient cohorts.",
      "The model itself is deterministic and was parameterised entirely from in-vitro binding data,",
      "measured receptor densities and literature physiology (Table 1, Table S2) rather than fitted to",
      "these patients; the patient data are the external validation set (Figure 4).",
      "The virtual patient population used for the concentration-profile comparison is n = 100 with",
      "body weight 40-100 kg (normal), SC bioavailability 0.5-1.0, and 30% CV on the typical sBCMA,",
      "CD38 and CD3 densities.",
      sep = " "
    )
  )

  ini({
    # ================================================================
    # Binding kinetics. ISB 2001 affinities and dissociation rate
    # constants were MEASURED (Table S2); the association rate constants
    # are back-calculated as kon = koff / kd by the Data S1 initial
    # assignment rule, reproducing the Table 1 "Calculated" kon values
    # (CD3 1.80e-3, BCMA 2.40e-4, CD38 1.00e-3 1/(nM*s)).
    # All are fixed: none was estimated.
    # ================================================================
    kd_cd3    <- fixed(41);      label("ISB 2001 binding affinity to CD3 (nM)")                      # Table S2, Target Engagement
    koff_cd3  <- fixed(7.20e-2); label("ISB 2001 dissociation rate constant from CD3 (1/s)")         # Table S2, Target Engagement
    kd_bcma   <- fixed(1.9);     label("ISB 2001 binding affinity to BCMA (nM)")                     # Table S2, Target Engagement
    koff_bcma <- fixed(4.60e-4); label("ISB 2001 dissociation rate constant from BCMA (1/s)")        # Table S2, Target Engagement
    kd_cd38   <- fixed(2.70);    label("ISB 2001 binding affinity to CD38 (nM)")                     # Table S2, Target Engagement
    koff_cd38 <- fixed(2.80e-3); label("ISB 2001 dissociation rate constant from CD38 (1/s)")        # Table S2, Target Engagement

    # Avidity factor for the second tumour-cell arm engaging after the
    # first. Data S1 "Constant parameters": KAI, fixed to 1 (no avidity).
    kai <- fixed(1); label("Avidity factor for second tumour-antigen arm (unitless)")                # Data S1, KAI

    # ================================================================
    # Receptor internalisation. Only NON-crosslinking complexes are
    # internalised; every CD3-engaged trimer / tetramer is fixed to zero
    # (Data S1 footnote b, and Methods "Only non-crosslinking complexes
    # were assumed to be internalized").
    # ================================================================
    kint_cd3  <- fixed(0.173); label("CD3 internalisation rate constant (1/h)")                      # Table S2; from internalisation t1/2 ~ 4 h
    kint_bcma <- fixed(0.35);  label("BCMA internalisation rate constant (1/h)")                     # Table S2; from internalisation t1/2 ~ 2 h
    kint_cd38 <- fixed(0.35);  label("CD38 internalisation rate constant (1/h)")                     # Table S2; from internalisation t1/2 ~ 2 h
    kint_crosslinked <- fixed(0); label("Internalisation rate constant of CD3-crosslinked complexes (1/h)")  # Data S1 footnote b: Kint6-11, kint16-17 fixed to 0

    # ================================================================
    # Minimal PBPK disposition (human). Physiological volumes and lymph
    # flows are literature IgG values; the leaky-tissue compartment of the
    # published mPBPK is split into bone marrow and non-marrow leaky
    # tissue so that bone-marrow nACT can be simulated.
    # ================================================================
    lvp            <- fixed(log(2.81)); label("Plasma volume (L)")                                        # Table S2, Human QSP model
    vtight        <- fixed(8.1);  label("Interstitial fluid volume of tight tissue (L)")            # Table S2, Human QSP model
    vlymph        <- fixed(5.2);  label("Lymph volume (L)")                                         # Table S2, Human QSP model
    vleaky_wo_bm  <- fixed(4);    label("Interstitial fluid volume of leaky tissue excluding bone marrow (L)")  # Table S2, Human QSP model
    vbonemarrow   <- fixed(0.37); label("Interstitial fluid volume of bone marrow (L)")             # Table S2, Human QSP model

    lflow  <- fixed(2.9);  label("Total lymph flow (L/day)")                                        # Table S2, Human QSP model: L
    lflow1 <- fixed(0.96); label("Lymph flow for tight tissue (L/day)")                             # Table S2, Human QSP model: L1
    lflow2 <- fixed(1.95); label("Lymph flow for leaky tissue (L/day)")                             # Table S2, Human QSP model: L2

    sigma_tight <- fixed(0.945); label("Vascular reflection coefficient, tight tissue (unitless)")   # Table S2, Human QSP model: sigma1
    sigma_leaky <- fixed(0.697); label("Vascular reflection coefficient, leaky tissue (unitless)")   # Table S2, Human QSP model: sigma2
    sigma_l     <- fixed(0.2);   label("Lymphatic capillary reflection coefficient (unitless)")      # Table S2, Human QSP model: sigmaL

    # Non-specific clearance, allometrically scaled to a 70 kg human from
    # the hFcRn Tg32 mouse estimate (0.272 mL/day, Table 1). Bioavailability
    # and absorption rate constant after SC injection are literature values
    # for ISB 2001 (source reference 36), not fitted here.
    lcl <- fixed(log(0.253)); label("Non-specific plasma clearance (L/day)")                               # Table 1, Human QSP model; scaled from hFcRn mouse
    lka <- fixed(log(0.5));   label("Subcutaneous absorption rate constant (1/day)")                       # Table S2, Human QSP model: KA
    ba <- fixed(0.75);  label("Subcutaneous bioavailability (unitless)")                             # Table S2, Human QSP model: BA

    # ================================================================
    # Target densities (molecules per cell) and soluble target
    # concentrations. All literature or measured; none estimated.
    # ================================================================
    cd3_per_cell            <- fixed(100000); label("CD3 density on CD4+/CD8+ T cells (molecule/cell)")   # Table S2, Human QSP model
    bcma_per_cell           <- fixed(2000);   label("BCMA density on plasma cell (molecule/cell)")        # Table S2, Human QSP model
    cd38_per_cell           <- fixed(120486); label("CD38 density on plasma cell (molecule/cell)")        # Table S2, Human QSP model
    cd38_per_monocyte       <- fixed(12000);  label("CD38 density on monocyte (molecule/cell)")           # Table S2, Human QSP model
    cd38_per_nk             <- fixed(22000);  label("CD38 density on NK cell (molecule/cell)")            # Table S2, Human QSP model
    cd38_per_neutrophil     <- fixed(1800);   label("CD38 density on neutrophil (molecule/cell)")         # Table S2, Human QSP model

    sbcma_plasma     <- fixed(15);   label("Soluble BCMA concentration in plasma (nM)")               # Table S2, Human QSP model
    scd38_plasma     <- fixed(0.01); label("Soluble CD38 concentration in plasma (nM)")               # Table S2, Human QSP model
    sbcma_bonemarrow <- fixed(15);   label("Soluble BCMA concentration in bone marrow (nM)")          # Table S2: sBCMA in plasma and BM compartment
    scd38_bonemarrow <- fixed(0.01); label("Soluble CD38 concentration in bone marrow (nM)")          # Table S2: sCD38 in plasma and BM compartment

    # ================================================================
    # Cell counts. These are TOTAL counts in the compartment, converted to
    # a molar target concentration in model() by multiplying by the target
    # density, dividing by Avogadro's number and by the compartment volume.
    # ================================================================
    n_cd4_plasma        <- fixed(5.25e9);  label("CD4+ T cell count in the plasma compartment (cells)")     # Table S2: Central CD4+ T-cell
    n_cd8_plasma        <- fixed(2.69e9);  label("CD8+ T cell count in the plasma compartment (cells)")     # Table S2: Central CD8+ T-cell
    n_nk_plasma         <- fixed(1.08e9);  label("NK cell count in the plasma compartment (cells)")         # Table S2: Central NK cell
    n_monocyte_plasma   <- fixed(2.62e9);  label("Monocyte count in the plasma compartment (cells)")        # Table S2: Central Monocyte
    n_neutrophil_plasma <- fixed(2.02e10); label("Neutrophil count in the plasma compartment (cells)")      # Table S2: Central Neutrophil

    # Plasma cells (the CD38/BCMA-expressing tumour cells) reside in the
    # bone marrow; Data S1 sets the plasma-compartment count TC to 0.
    n_tumor_plasma     <- fixed(0);       label("Plasma (myeloma) cell count in the plasma compartment (cells)")  # Data S1 background parameters: TC = 0
    n_tumor_bonemarrow <- fixed(6.55e9);  label("Plasma (myeloma) cell count in bone marrow (cells)")       # Table S2: Plasma cell; Data S1 TC_bm
    n_cd4_bonemarrow   <- fixed(9.00e8);  label("CD4+ T cell count in bone marrow (cells)")                 # Table S2: BM CD4+ T-cell
    n_cd8_bonemarrow   <- fixed(9.00e8);  label("CD8+ T cell count in bone marrow (cells)")                 # Table S2: BM CD8+ T-cell

    # ================================================================
    # Residual error. The source is a deterministic QSP simulation model:
    # it reports no residual-error model and no between-subject variance,
    # so the proportional residual error is carried at zero rather than
    # invented. See the vignette "Assumptions and deviations" section.
    # ================================================================
    propSd <- fixed(0); label("Proportional residual error (fraction; ZERO - not reported in source)")
  })

  model({
    # ================================================================
    # 1. Unit harmonisation.
    #    Binding constants are reported per second and internalisation
    #    rate constants per hour; the rest of the model runs in days.
    # ================================================================
    # Back-transform the log-scale fixed effects (parameter naming convention).
    vp <- exp(lvp)
    cl <- exp(lcl)
    ka <- exp(lka)

    kon_cd3_d  <- koff_cd3 / kd_cd3   * 86400   # 1/(nM*day); Data S1: kon_CD3 = koff_CD3/Kd_CD3
    kon_bcma_d <- koff_bcma / kd_bcma * 86400   # 1/(nM*day); Data S1: kon_BCMA = koff_BCMA/Kd_BCMA
    kon_cd38_d <- koff_cd38 / kd_cd38 * 86400   # 1/(nM*day); Data S1: kon_CD38 = koff_CD38/Kd_CD38
    koff_cd3_d  <- koff_cd3  * 86400            # 1/day
    koff_bcma_d <- koff_bcma * 86400            # 1/day
    koff_cd38_d <- koff_cd38 * 86400            # 1/day

    kint_cd3_d  <- kint_cd3  * 24               # 1/day
    kint_bcma_d <- kint_bcma * 24               # 1/day
    kint_cd38_d <- kint_cd38 * 24               # 1/day
    kint_x_d    <- kint_crosslinked * 24        # 1/day (zero: CD3-crosslinked complexes are not internalised)

    kel <- cl / vp                              # 1/day; Data S1 initial assignment: kel = CL/plasma

    # Fractional split of the leaky-tissue lymph flow between bone marrow
    # and the remaining leaky tissue, exactly as written in Data S1.
    f_leaky <- vleaky_wo_bm / (vleaky_wo_bm + vbonemarrow)
    f_bm    <- vbonemarrow  / (vleaky_wo_bm + vbonemarrow)

    # ================================================================
    # 2. Total target concentrations (nM).
    #    Data S1 variable assignment rules, e.g.
    #      CD3CD4tot = CD4tot*CD3percell/AC/plasma
    #    The 1e9 converts mol/L to nmol/L (nM); Simbiology carried units
    #    implicitly, so the factor is not written in the source listing.
    # ================================================================
    avogadro <- 6.02e23                                                    # 1/mol; Data S1 background parameter AC

    cd3cd4_tot <- n_cd4_plasma * cd3_per_cell  / avogadro / vp * 1e9
    cd3cd8_tot <- n_cd8_plasma * cd3_per_cell  / avogadro / vp * 1e9
    bcma_tot   <- n_tumor_plasma * bcma_per_cell / avogadro / vp * 1e9
    cd38_tot   <- n_tumor_plasma * cd38_per_cell / avogadro / vp * 1e9
    cd38_other_tot <- (cd38_per_monocyte * n_monocyte_plasma +
                         cd38_per_nk * n_nk_plasma +
                         cd38_per_neutrophil * n_neutrophil_plasma) / avogadro / vp * 1e9

    cd3cd4_tot_bm <- n_cd4_bonemarrow * cd3_per_cell / avogadro / vbonemarrow * 1e9
    cd3cd8_tot_bm <- n_cd8_bonemarrow * cd3_per_cell / avogadro / vbonemarrow * 1e9
    bcma_tot_bm   <- n_tumor_bonemarrow * bcma_per_cell / avogadro / vbonemarrow * 1e9
    cd38_tot_bm   <- n_tumor_bonemarrow * cd38_per_cell / avogadro / vbonemarrow * 1e9

    # ================================================================
    # 3. Free (unbound) target concentrations, by conservation.
    #    Data S1 "Variable condition assignment" block, transcribed term
    #    for term.
    # ================================================================
    u_cd3cd4 <- cd3cd4_tot - drug_cd3cd4 - trimer_cd3cd4_bcma - trimer_cd3cd4_cd38 -
      tetramer_cd3cd4 - complex_cd3cd4_sbcma - complex_cd3cd4_sbcma_cd38 -
      complex_cd3cd4_scd38 - trimer_cd3cd4_cd38other
    u_cd3cd8 <- cd3cd8_tot - drug_cd3cd8 - trimer_cd3cd8_bcma - trimer_cd3cd8_cd38 -
      tetramer_cd3cd8 - complex_cd3cd8_sbcma - complex_cd3cd8_sbcma_cd38 -
      complex_cd3cd8_scd38 - trimer_cd3cd8_cd38other
    u_bcma <- bcma_tot - drug_bcma - trimer_cd3cd4_bcma - trimer_cd3cd8_bcma -
      trimer_bcma_cd38 - tetramer_cd3cd4 - tetramer_cd3cd8
    u_cd38 <- cd38_tot - drug_cd38 - trimer_cd3cd4_cd38 - trimer_cd3cd8_cd38 -
      trimer_bcma_cd38 - tetramer_cd3cd4 - tetramer_cd3cd8 - complex_sbcma_cd38 -
      complex_cd3cd4_sbcma_cd38 - complex_cd3cd8_sbcma_cd38
    u_sbcma <- sbcma_plasma - (complex_sbcma + complex_cd3cd4_sbcma + complex_cd3cd8_sbcma +
                                 complex_sbcma_cd38 + complex_sbcma_scd38 +
                                 complex_cd3cd4_sbcma_cd38 + complex_cd3cd8_sbcma_cd38)
    u_scd38 <- scd38_plasma - (complex_scd38 + complex_cd3cd4_scd38 +
                                 complex_cd3cd8_scd38 + complex_sbcma_scd38)
    u_cd38_other <- cd38_other_tot - (drug_cd38_other + trimer_cd3cd4_cd38other +
                                        trimer_cd3cd8_cd38other)

    u_cd3cd4_bm <- cd3cd4_tot_bm - drug_cd3cd4_bonemarrow - trimer_cd3cd4_bcma_bonemarrow -
      trimer_cd3cd4_cd38_bonemarrow - tetramer_cd3cd4_bonemarrow -
      complex_cd3cd4_sbcma_bonemarrow - complex_cd3cd4_sbcma_cd38_bonemarrow -
      complex_cd3cd4_scd38_bonemarrow
    u_cd3cd8_bm <- cd3cd8_tot_bm - drug_cd3cd8_bonemarrow - trimer_cd3cd8_bcma_bonemarrow -
      trimer_cd3cd8_cd38_bonemarrow - tetramer_cd3cd8_bonemarrow -
      complex_cd3cd8_sbcma_bonemarrow - complex_cd3cd8_sbcma_cd38_bonemarrow -
      complex_cd3cd8_scd38_bonemarrow
    u_bcma_bm <- bcma_tot_bm - drug_bcma_bonemarrow - trimer_cd3cd4_bcma_bonemarrow -
      trimer_cd3cd8_bcma_bonemarrow - trimer_bcma_cd38_bonemarrow -
      tetramer_cd3cd4_bonemarrow - tetramer_cd3cd8_bonemarrow
    u_cd38_bm <- cd38_tot_bm - drug_cd38_bonemarrow - trimer_cd3cd4_cd38_bonemarrow -
      trimer_cd3cd8_cd38_bonemarrow - trimer_bcma_cd38_bonemarrow -
      tetramer_cd3cd4_bonemarrow - tetramer_cd3cd8_bonemarrow -
      complex_sbcma_cd38_bonemarrow - complex_cd3cd4_sbcma_cd38_bonemarrow -
      complex_cd3cd8_sbcma_cd38_bonemarrow
    u_sbcma_bm <- sbcma_bonemarrow - (complex_sbcma_bonemarrow +
                                        complex_cd3cd4_sbcma_bonemarrow +
                                        complex_cd3cd8_sbcma_bonemarrow +
                                        complex_sbcma_cd38_bonemarrow +
                                        complex_cd3cd4_sbcma_cd38_bonemarrow +
                                        complex_cd3cd8_sbcma_cd38_bonemarrow)
    # Data S1 writes the bone-marrow sCD38 balance against the PLASMA
    # soluble CD38 pool (usCD38) for the sTX6_bm and sCX1_bm/sCX2_bm terms;
    # the free bone-marrow sCD38 itself is not assigned there. It is
    # carried here as the bone-marrow pool less its own bound species,
    # which is the only mass-conserving reading.
    u_scd38_bm <- scd38_bonemarrow - (complex_scd38_bonemarrow +
                                        complex_cd3cd4_scd38_bonemarrow +
                                        complex_cd3cd8_scd38_bonemarrow +
                                        complex_sbcma_scd38_bonemarrow)

    # ================================================================
    # 4. ODE system - PLASMA compartment.
    #    Data S1 divides every flux by the compartment volume and
    #    multiplies the binding fluxes back by it, so the binding terms
    #    are in concentration units directly and only the transport and
    #    absorption terms carry a 1/V.
    # ================================================================
    d/dt(plasma) <-
      -(kon_cd3_d * plasma * u_cd3cd4 - koff_cd3_d * drug_cd3cd4) -
      (kon_cd3_d * plasma * u_cd3cd8 - koff_cd3_d * drug_cd3cd8) -
      (kon_bcma_d * plasma * u_bcma - koff_bcma_d * drug_bcma) -
      (kon_cd38_d * plasma * u_cd38 - koff_cd38_d * drug_cd38) -
      (kon_bcma_d * plasma * u_sbcma - koff_bcma_d * complex_sbcma) -
      (kon_cd38_d * plasma * u_scd38 - koff_cd38_d * complex_scd38) -
      (kon_cd38_d * plasma * u_cd38_other - koff_cd38_d * drug_cd38_other) -
      kel * plasma +
      (-(1 - sigma_leaky) * lflow2 * f_leaky * plasma -
         (1 - sigma_leaky) * lflow2 * f_bm * plasma -
         (1 - sigma_tight) * lflow1 * plasma +
         lflow * lymph +
         ka * ba * depot) / vp

    d/dt(drug_cd3cd4) <-
      (kon_cd3_d * plasma * u_cd3cd4 - koff_cd3_d * drug_cd3cd4) -
      kint_cd3_d * drug_cd3cd4 -
      (kon_bcma_d * drug_cd3cd4 * u_bcma - koff_bcma_d * trimer_cd3cd4_bcma) -
      (kon_cd38_d * drug_cd3cd4 * u_cd38 - koff_cd38_d * trimer_cd3cd4_cd38) -
      (kon_bcma_d * drug_cd3cd4 * u_sbcma - koff_bcma_d * complex_cd3cd4_sbcma) -
      (kon_cd38_d * drug_cd3cd4 * u_scd38 - koff_cd38_d * complex_cd3cd4_scd38) -
      (kon_cd38_d * drug_cd3cd4 * u_cd38_other - koff_cd38_d * trimer_cd3cd4_cd38other)

    d/dt(drug_cd3cd8) <-
      (kon_cd3_d * plasma * u_cd3cd8 - koff_cd3_d * drug_cd3cd8) -
      kint_cd3_d * drug_cd3cd8 -
      (kon_cd38_d * u_cd38 * drug_cd3cd8 - koff_cd38_d * trimer_cd3cd8_cd38) -
      (kon_bcma_d * drug_cd3cd8 * u_bcma - koff_bcma_d * trimer_cd3cd8_bcma) -
      (kon_bcma_d * drug_cd3cd8 * u_sbcma - koff_bcma_d * complex_cd3cd8_sbcma) -
      (kon_cd38_d * drug_cd3cd8 * u_scd38 - koff_cd38_d * complex_cd3cd8_scd38) -
      (kon_cd38_d * drug_cd3cd8 * u_cd38_other - koff_cd38_d * trimer_cd3cd8_cd38other)

    d/dt(drug_bcma) <-
      (kon_bcma_d * plasma * u_bcma - koff_bcma_d * drug_bcma) -
      kint_bcma_d * drug_bcma -
      (kon_cd3_d * drug_bcma * u_cd3cd4 - koff_cd3_d * trimer_cd3cd4_bcma) -
      (kon_cd3_d * u_cd3cd8 * drug_bcma - koff_cd3_d * trimer_cd3cd8_bcma) -
      (kai * kon_cd38_d * u_cd38 * drug_bcma - koff_cd38_d * trimer_bcma_cd38) -
      (kon_bcma_d * drug_bcma * u_sbcma - koff_bcma_d * complex_sbcma_cd38)

    d/dt(drug_cd38) <-
      (kon_cd38_d * plasma * u_cd38 - koff_cd38_d * drug_cd38) -
      kint_cd38_d * drug_cd38 -
      (kon_cd3_d * drug_cd38 * u_cd3cd4 - koff_cd3_d * trimer_cd3cd4_cd38) -
      (kon_cd3_d * drug_cd38 * u_cd3cd8 - koff_cd3_d * trimer_cd3cd8_cd38) -
      (kai * kon_bcma_d * drug_cd38 * u_bcma - koff_bcma_d * trimer_bcma_cd38)

    d/dt(trimer_cd3cd8_cd38) <-
      (kon_cd38_d * u_cd38 * drug_cd3cd8 - koff_cd38_d * trimer_cd3cd8_cd38) +
      (kon_cd3_d * drug_cd38 * u_cd3cd8 - koff_cd3_d * trimer_cd3cd8_cd38) -
      (kai * kon_bcma_d * trimer_cd3cd8_cd38 * u_bcma - koff_bcma_d * tetramer_cd3cd8) -
      kint_x_d * trimer_cd3cd8_cd38

    d/dt(trimer_cd3cd8_bcma) <-
      (kon_bcma_d * drug_cd3cd8 * u_bcma - koff_bcma_d * trimer_cd3cd8_bcma) +
      (kon_cd3_d * u_cd3cd8 * drug_bcma - koff_cd3_d * trimer_cd3cd8_bcma) -
      2 * (kai * kon_cd38_d * trimer_cd3cd8_bcma * u_cd38 - koff_cd38_d * tetramer_cd3cd8) -
      kint_x_d * trimer_cd3cd8_bcma

    d/dt(trimer_cd3cd4_bcma) <-
      (kon_bcma_d * drug_cd3cd4 * u_bcma - koff_bcma_d * trimer_cd3cd4_bcma) +
      (kon_cd3_d * drug_bcma * u_cd3cd4 - koff_cd3_d * trimer_cd3cd4_bcma) -
      (kai * kon_cd38_d * trimer_cd3cd4_bcma * u_cd38 - koff_cd38_d * tetramer_cd3cd4) -
      kint_x_d * trimer_cd3cd4_bcma

    d/dt(trimer_cd3cd4_cd38) <-
      (kon_cd38_d * drug_cd3cd4 * u_cd38 - koff_cd38_d * trimer_cd3cd4_cd38) +
      (kon_cd3_d * drug_cd38 * u_cd3cd4 - koff_cd3_d * trimer_cd3cd4_cd38) -
      (kai * kon_bcma_d * trimer_cd3cd4_cd38 * u_bcma - koff_bcma_d * tetramer_cd3cd4) -
      kint_x_d * trimer_cd3cd4_cd38

    d/dt(tetramer_cd3cd8) <-
      (kai * kon_cd38_d * trimer_cd3cd8_bcma * u_cd38 - koff_cd38_d * tetramer_cd3cd8) +
      (kai * kon_bcma_d * u_bcma * trimer_cd3cd8_cd38 - koff_bcma_d * tetramer_cd3cd8) +
      (kon_cd3_d * trimer_bcma_cd38 * u_cd3cd8 - koff_cd3_d * tetramer_cd3cd8) -
      kint_x_d * tetramer_cd3cd8

    d/dt(tetramer_cd3cd4) <-
      (kai * kon_cd38_d * trimer_cd3cd4_bcma * u_cd38 - koff_cd38_d * tetramer_cd3cd4) +
      (kai * kon_bcma_d * trimer_cd3cd4_cd38 * u_bcma - koff_bcma_d * tetramer_cd3cd4) +
      (kon_cd3_d * trimer_bcma_cd38 * u_cd3cd4 - koff_cd3_d * tetramer_cd3cd4) -
      kint_x_d * tetramer_cd3cd4

    d/dt(trimer_bcma_cd38) <-
      (kai * kon_bcma_d * drug_cd38 * u_bcma - koff_bcma_d * trimer_bcma_cd38) +
      (kai * kon_cd38_d * u_cd38 * drug_bcma - koff_cd38_d * trimer_bcma_cd38) -
      (kon_cd3_d * trimer_bcma_cd38 * u_cd3cd4 - koff_cd3_d * tetramer_cd3cd4) -
      (kon_cd3_d * trimer_bcma_cd38 * u_cd3cd8 - koff_cd3_d * tetramer_cd3cd8) -
      kint_bcma_d * trimer_bcma_cd38

    d/dt(complex_sbcma) <-
      (kon_bcma_d * plasma * u_sbcma - koff_bcma_d * complex_sbcma) -
      (kon_cd3_d * complex_sbcma * u_cd3cd4 - koff_cd3_d * complex_cd3cd4_sbcma) -
      (kon_cd3_d * complex_sbcma * u_cd3cd8 - koff_cd3_d * complex_cd3cd8_sbcma) -
      (kon_cd38_d * complex_sbcma * u_cd38 - koff_cd38_d * complex_sbcma_cd38) -
      (kon_cd38_d * complex_sbcma * u_scd38 - koff_cd38_d * complex_sbcma_scd38) -
      kel * complex_sbcma

    d/dt(complex_cd3cd4_sbcma) <-
      (kon_cd3_d * complex_sbcma * u_cd3cd4 - koff_cd3_d * complex_cd3cd4_sbcma) +
      (kon_bcma_d * drug_cd3cd4 * u_sbcma - koff_bcma_d * complex_cd3cd4_sbcma) -
      (kon_cd38_d * complex_cd3cd4_sbcma * u_cd38 - koff_cd38_d * complex_cd3cd4_sbcma_cd38) -
      kint_cd3_d * complex_cd3cd4_sbcma

    d/dt(complex_cd3cd8_sbcma) <-
      (kon_bcma_d * drug_cd3cd8 * u_sbcma - koff_bcma_d * complex_cd3cd8_sbcma) +
      (kon_cd3_d * complex_sbcma * u_cd3cd8 - koff_cd3_d * complex_cd3cd8_sbcma) -
      (kon_cd38_d * complex_cd3cd8_sbcma * u_cd38 - koff_cd38_d * complex_cd3cd8_sbcma_cd38) -
      kint_cd3_d * complex_cd3cd8_sbcma

    d/dt(complex_sbcma_cd38) <-
      (kon_bcma_d * drug_bcma * u_sbcma - koff_bcma_d * complex_sbcma_cd38) +
      (kon_cd38_d * complex_sbcma * u_cd38 - koff_cd38_d * complex_sbcma_cd38) -
      (kon_cd3_d * complex_sbcma_cd38 * u_cd3cd4 - koff_cd3_d * complex_cd3cd4_sbcma_cd38) -
      (kon_cd3_d * complex_sbcma_cd38 * u_cd3cd8 - koff_cd3_d * complex_cd3cd8_sbcma_cd38) -
      kint_cd38_d * complex_sbcma_cd38

    d/dt(complex_cd3cd4_sbcma_cd38) <-
      (kon_cd38_d * complex_cd3cd4_sbcma * u_cd38 - koff_cd38_d * complex_cd3cd4_sbcma_cd38) +
      (kon_cd3_d * complex_sbcma_cd38 * u_cd3cd4 - koff_cd3_d * complex_cd3cd4_sbcma_cd38) -
      kint_x_d * complex_cd3cd4_sbcma_cd38

    d/dt(complex_cd3cd8_sbcma_cd38) <-
      (kon_cd38_d * complex_cd3cd8_sbcma * u_cd38 - koff_cd38_d * complex_cd3cd8_sbcma_cd38) +
      (kon_cd3_d * complex_sbcma_cd38 * u_cd3cd8 - koff_cd3_d * complex_cd3cd8_sbcma_cd38) -
      kint_x_d * complex_cd3cd8_sbcma_cd38

    d/dt(complex_scd38) <-
      (kon_cd38_d * plasma * u_scd38 - koff_cd38_d * complex_scd38) -
      (kon_cd3_d * complex_scd38 * u_cd3cd4 - koff_cd3_d * complex_cd3cd4_scd38) -
      (kon_cd3_d * complex_scd38 * u_cd3cd8 - koff_cd3_d * complex_cd3cd8_scd38) -
      (kon_bcma_d * complex_scd38 * u_sbcma - koff_bcma_d * complex_sbcma_scd38) -
      kel * complex_scd38

    d/dt(drug_cd38_other) <-
      -kint_cd38_d * drug_cd38_other +
      (kon_cd38_d * plasma * u_cd38_other - koff_cd38_d * drug_cd38_other) -
      (kon_cd3_d * drug_cd38_other * u_cd3cd4 - koff_cd3_d * trimer_cd3cd4_cd38other) -
      (kon_cd3_d * drug_cd38_other * u_cd3cd8 - koff_cd3_d * trimer_cd3cd8_cd38other)

    d/dt(complex_cd3cd4_scd38) <-
      (kon_cd3_d * complex_scd38 * u_cd3cd4 - koff_cd3_d * complex_cd3cd4_scd38) +
      (kon_cd38_d * drug_cd3cd4 * u_scd38 - koff_cd38_d * complex_cd3cd4_scd38) -
      kint_cd3_d * complex_cd3cd4_scd38

    d/dt(complex_cd3cd8_scd38) <-
      (kon_cd38_d * drug_cd3cd8 * u_scd38 - koff_cd38_d * complex_cd3cd8_scd38) +
      (kon_cd3_d * complex_scd38 * u_cd3cd8 - koff_cd3_d * complex_cd3cd8_scd38) -
      kint_cd3_d * complex_cd3cd8_scd38

    d/dt(complex_sbcma_scd38) <-
      (kon_cd38_d * complex_sbcma * u_scd38 - koff_cd38_d * complex_sbcma_scd38) +
      (kon_bcma_d * complex_scd38 * u_sbcma - koff_bcma_d * complex_sbcma_scd38) -
      kel * complex_sbcma_scd38

    d/dt(trimer_cd3cd4_cd38other) <-
      (kon_cd38_d * drug_cd3cd4 * u_cd38_other - koff_cd38_d * trimer_cd3cd4_cd38other) +
      (kon_cd3_d * drug_cd38_other * u_cd3cd4 - koff_cd3_d * trimer_cd3cd4_cd38other) -
      kint_x_d * trimer_cd3cd4_cd38other

    d/dt(trimer_cd3cd8_cd38other) <-
      (kon_cd38_d * drug_cd3cd8 * u_cd38_other - koff_cd38_d * trimer_cd3cd8_cd38other) +
      (kon_cd3_d * drug_cd38_other * u_cd3cd8 - koff_cd3_d * trimer_cd3cd8_cd38other) -
      kint_x_d * trimer_cd3cd8_cd38other

    # ================================================================
    # 5. ODE system - BONE MARROW compartment.
    #    Identical target-engagement network, except that CD38 on
    #    monocytes / NK cells / neutrophils is modelled in plasma only.
    # ================================================================
    d/dt(bonemarrow) <-
      -(kon_cd3_d * bonemarrow * u_cd3cd4_bm - koff_cd3_d * drug_cd3cd4_bonemarrow) -
      (kon_cd3_d * bonemarrow * u_cd3cd8_bm - koff_cd3_d * drug_cd3cd8_bonemarrow) -
      (kon_bcma_d * bonemarrow * u_bcma_bm - koff_bcma_d * drug_bcma_bonemarrow) -
      (kon_cd38_d * bonemarrow * u_cd38_bm - koff_cd38_d * drug_cd38_bonemarrow) -
      (kon_bcma_d * bonemarrow * u_sbcma_bm - koff_bcma_d * complex_sbcma_bonemarrow) -
      (kon_cd38_d * bonemarrow * u_scd38_bm - koff_cd38_d * complex_scd38_bonemarrow) +
      ((1 - sigma_leaky) * lflow2 * f_bm * plasma -
         lflow2 * (1 - sigma_l) * f_bm * bonemarrow) / vbonemarrow

    d/dt(drug_cd3cd4_bonemarrow) <-
      (kon_cd3_d * bonemarrow * u_cd3cd4_bm - koff_cd3_d * drug_cd3cd4_bonemarrow) -
      kint_cd3_d * drug_cd3cd4_bonemarrow -
      (kon_bcma_d * drug_cd3cd4_bonemarrow * u_bcma_bm - koff_bcma_d * trimer_cd3cd4_bcma_bonemarrow) -
      (kon_cd38_d * drug_cd3cd4_bonemarrow * u_cd38_bm - koff_cd38_d * trimer_cd3cd4_cd38_bonemarrow) -
      (kon_bcma_d * drug_cd3cd4_bonemarrow * u_sbcma_bm - koff_bcma_d * complex_cd3cd4_sbcma_bonemarrow) -
      (kon_cd38_d * drug_cd3cd4_bonemarrow * u_scd38_bm - koff_cd38_d * complex_cd3cd4_scd38_bonemarrow)

    d/dt(drug_cd3cd8_bonemarrow) <-
      (kon_cd3_d * bonemarrow * u_cd3cd8_bm - koff_cd3_d * drug_cd3cd8_bonemarrow) -
      kint_cd3_d * drug_cd3cd8_bonemarrow -
      (kon_cd38_d * u_cd38_bm * drug_cd3cd8_bonemarrow - koff_cd38_d * trimer_cd3cd8_cd38_bonemarrow) -
      (kon_bcma_d * drug_cd3cd8_bonemarrow * u_bcma_bm - koff_bcma_d * trimer_cd3cd8_bcma_bonemarrow) -
      (kon_bcma_d * drug_cd3cd8_bonemarrow * u_sbcma_bm - koff_bcma_d * complex_cd3cd8_sbcma_bonemarrow) -
      (kon_cd38_d * drug_cd3cd8_bonemarrow * u_scd38_bm - koff_cd38_d * complex_cd3cd8_scd38_bonemarrow)

    d/dt(drug_bcma_bonemarrow) <-
      (kon_bcma_d * bonemarrow * u_bcma_bm - koff_bcma_d * drug_bcma_bonemarrow) -
      kint_bcma_d * drug_bcma_bonemarrow -
      (kon_cd3_d * drug_bcma_bonemarrow * u_cd3cd4_bm - koff_cd3_d * trimer_cd3cd4_bcma_bonemarrow) -
      (kon_cd3_d * u_cd3cd8_bm * drug_bcma_bonemarrow - koff_cd3_d * trimer_cd3cd8_bcma_bonemarrow) -
      (kai * kon_cd38_d * u_cd38_bm * drug_bcma_bonemarrow - koff_cd38_d * trimer_bcma_cd38_bonemarrow) -
      (kon_bcma_d * drug_bcma_bonemarrow * u_sbcma_bm - koff_bcma_d * complex_sbcma_cd38_bonemarrow)

    d/dt(drug_cd38_bonemarrow) <-
      (kon_cd38_d * bonemarrow * u_cd38_bm - koff_cd38_d * drug_cd38_bonemarrow) -
      kint_cd38_d * drug_cd38_bonemarrow -
      (kon_cd3_d * drug_cd38_bonemarrow * u_cd3cd4_bm - koff_cd3_d * trimer_cd3cd4_cd38_bonemarrow) -
      (kon_cd3_d * drug_cd38_bonemarrow * u_cd3cd8_bm - koff_cd3_d * trimer_cd3cd8_cd38_bonemarrow) -
      (kai * kon_bcma_d * drug_cd38_bonemarrow * u_bcma_bm - koff_bcma_d * trimer_bcma_cd38_bonemarrow)

    d/dt(trimer_cd3cd8_cd38_bonemarrow) <-
      (kon_cd38_d * u_cd38_bm * drug_cd3cd8_bonemarrow - koff_cd38_d * trimer_cd3cd8_cd38_bonemarrow) +
      (kon_cd3_d * drug_cd38_bonemarrow * u_cd3cd8_bm - koff_cd3_d * trimer_cd3cd8_cd38_bonemarrow) -
      (kai * kon_bcma_d * trimer_cd3cd8_cd38_bonemarrow * u_bcma_bm - koff_bcma_d * tetramer_cd3cd8_bonemarrow) -
      kint_x_d * trimer_cd3cd8_cd38_bonemarrow

    d/dt(trimer_cd3cd8_bcma_bonemarrow) <-
      (kon_bcma_d * drug_cd3cd8_bonemarrow * u_bcma_bm - koff_bcma_d * trimer_cd3cd8_bcma_bonemarrow) +
      (kon_cd3_d * u_cd3cd8_bm * drug_bcma_bonemarrow - koff_cd3_d * trimer_cd3cd8_bcma_bonemarrow) -
      2 * (kai * kon_cd38_d * trimer_cd3cd8_bcma_bonemarrow * u_cd38_bm - koff_cd38_d * tetramer_cd3cd8_bonemarrow) -
      kint_x_d * trimer_cd3cd8_bcma_bonemarrow

    d/dt(trimer_cd3cd4_bcma_bonemarrow) <-
      (kon_bcma_d * drug_cd3cd4_bonemarrow * u_bcma_bm - koff_bcma_d * trimer_cd3cd4_bcma_bonemarrow) +
      (kon_cd3_d * drug_bcma_bonemarrow * u_cd3cd4_bm - koff_cd3_d * trimer_cd3cd4_bcma_bonemarrow) -
      (kai * kon_cd38_d * trimer_cd3cd4_bcma_bonemarrow * u_cd38_bm - koff_cd38_d * tetramer_cd3cd4_bonemarrow) -
      kint_x_d * trimer_cd3cd4_bcma_bonemarrow

    d/dt(trimer_cd3cd4_cd38_bonemarrow) <-
      (kon_cd38_d * drug_cd3cd4_bonemarrow * u_cd38_bm - koff_cd38_d * trimer_cd3cd4_cd38_bonemarrow) +
      (kon_cd3_d * drug_cd38_bonemarrow * u_cd3cd4_bm - koff_cd3_d * trimer_cd3cd4_cd38_bonemarrow) -
      (kai * kon_bcma_d * trimer_cd3cd4_cd38_bonemarrow * u_bcma_bm - koff_bcma_d * tetramer_cd3cd4_bonemarrow) -
      kint_x_d * trimer_cd3cd4_cd38_bonemarrow

    d/dt(tetramer_cd3cd8_bonemarrow) <-
      (kai * kon_cd38_d * trimer_cd3cd8_bcma_bonemarrow * u_cd38_bm - koff_cd38_d * tetramer_cd3cd8_bonemarrow) +
      (kai * kon_bcma_d * u_bcma_bm * trimer_cd3cd8_cd38_bonemarrow - koff_bcma_d * tetramer_cd3cd8_bonemarrow) +
      (kon_cd3_d * trimer_bcma_cd38_bonemarrow * u_cd3cd8_bm - koff_cd3_d * tetramer_cd3cd8_bonemarrow) -
      kint_x_d * tetramer_cd3cd8_bonemarrow

    d/dt(tetramer_cd3cd4_bonemarrow) <-
      (kai * kon_cd38_d * trimer_cd3cd4_bcma_bonemarrow * u_cd38_bm - koff_cd38_d * tetramer_cd3cd4_bonemarrow) +
      (kai * kon_bcma_d * trimer_cd3cd4_cd38_bonemarrow * u_bcma_bm - koff_bcma_d * tetramer_cd3cd4_bonemarrow) +
      (kon_cd3_d * trimer_bcma_cd38_bonemarrow * u_cd3cd4_bm - koff_cd3_d * tetramer_cd3cd4_bonemarrow) -
      kint_x_d * tetramer_cd3cd4_bonemarrow

    d/dt(trimer_bcma_cd38_bonemarrow) <-
      (kai * kon_bcma_d * drug_cd38_bonemarrow * u_bcma_bm - koff_bcma_d * trimer_bcma_cd38_bonemarrow) +
      (kai * kon_cd38_d * u_cd38_bm * drug_bcma_bonemarrow - koff_cd38_d * trimer_bcma_cd38_bonemarrow) -
      (kon_cd3_d * trimer_bcma_cd38_bonemarrow * u_cd3cd4_bm - koff_cd3_d * tetramer_cd3cd4_bonemarrow) -
      (kon_cd3_d * trimer_bcma_cd38_bonemarrow * u_cd3cd8_bm - koff_cd3_d * tetramer_cd3cd8_bonemarrow) -
      kint_bcma_d * trimer_bcma_cd38_bonemarrow

    d/dt(complex_sbcma_bonemarrow) <-
      (kon_bcma_d * bonemarrow * u_sbcma_bm - koff_bcma_d * complex_sbcma_bonemarrow) -
      (kon_cd3_d * complex_sbcma_bonemarrow * u_cd3cd4_bm - koff_cd3_d * complex_cd3cd4_sbcma_bonemarrow) -
      (kon_cd3_d * complex_sbcma_bonemarrow * u_cd3cd8_bm - koff_cd3_d * complex_cd3cd8_sbcma_bonemarrow) -
      (kon_cd38_d * complex_sbcma_bonemarrow * u_cd38_bm - koff_cd38_d * complex_sbcma_cd38_bonemarrow) -
      (kon_cd38_d * complex_sbcma_bonemarrow * u_scd38_bm - koff_cd38_d * complex_sbcma_scd38_bonemarrow) -
      kel * complex_sbcma_bonemarrow

    d/dt(complex_cd3cd4_sbcma_bonemarrow) <-
      (kon_cd3_d * complex_sbcma_bonemarrow * u_cd3cd4_bm - koff_cd3_d * complex_cd3cd4_sbcma_bonemarrow) +
      (kon_bcma_d * drug_cd3cd4_bonemarrow * u_sbcma_bm - koff_bcma_d * complex_cd3cd4_sbcma_bonemarrow) -
      (kon_cd38_d * complex_cd3cd4_sbcma_bonemarrow * u_cd38_bm - koff_cd38_d * complex_cd3cd4_sbcma_cd38_bonemarrow) -
      kint_cd3_d * complex_cd3cd4_sbcma_bonemarrow

    d/dt(complex_cd3cd8_sbcma_bonemarrow) <-
      (kon_bcma_d * drug_cd3cd8_bonemarrow * u_sbcma_bm - koff_bcma_d * complex_cd3cd8_sbcma_bonemarrow) +
      (kon_cd3_d * complex_sbcma_bonemarrow * u_cd3cd8_bm - koff_cd3_d * complex_cd3cd8_sbcma_bonemarrow) -
      (kon_cd38_d * complex_cd3cd8_sbcma_bonemarrow * u_cd38_bm - koff_cd38_d * complex_cd3cd8_sbcma_cd38_bonemarrow) -
      kint_cd3_d * complex_cd3cd8_sbcma_bonemarrow

    d/dt(complex_sbcma_cd38_bonemarrow) <-
      (kon_bcma_d * drug_bcma_bonemarrow * u_sbcma_bm - koff_bcma_d * complex_sbcma_cd38_bonemarrow) +
      (kon_cd38_d * complex_sbcma_bonemarrow * u_cd38_bm - koff_cd38_d * complex_sbcma_cd38_bonemarrow) -
      (kon_cd3_d * complex_sbcma_cd38_bonemarrow * u_cd3cd4_bm - koff_cd3_d * complex_cd3cd4_sbcma_cd38_bonemarrow) -
      (kon_cd3_d * complex_sbcma_cd38_bonemarrow * u_cd3cd8_bm - koff_cd3_d * complex_cd3cd8_sbcma_cd38_bonemarrow) -
      kint_cd38_d * complex_sbcma_cd38_bonemarrow

    d/dt(complex_cd3cd4_sbcma_cd38_bonemarrow) <-
      (kon_cd38_d * complex_cd3cd4_sbcma_bonemarrow * u_cd38_bm - koff_cd38_d * complex_cd3cd4_sbcma_cd38_bonemarrow) +
      (kon_cd3_d * complex_sbcma_cd38_bonemarrow * u_cd3cd4_bm - koff_cd3_d * complex_cd3cd4_sbcma_cd38_bonemarrow) -
      kint_x_d * complex_cd3cd4_sbcma_cd38_bonemarrow

    d/dt(complex_cd3cd8_sbcma_cd38_bonemarrow) <-
      (kon_cd38_d * complex_cd3cd8_sbcma_bonemarrow * u_cd38_bm - koff_cd38_d * complex_cd3cd8_sbcma_cd38_bonemarrow) +
      (kon_cd3_d * complex_sbcma_cd38_bonemarrow * u_cd3cd8_bm - koff_cd3_d * complex_cd3cd8_sbcma_cd38_bonemarrow) -
      kint_x_d * complex_cd3cd8_sbcma_cd38_bonemarrow

    d/dt(complex_scd38_bonemarrow) <-
      (kon_cd38_d * bonemarrow * u_scd38_bm - koff_cd38_d * complex_scd38_bonemarrow) -
      (kon_cd3_d * complex_scd38_bonemarrow * u_cd3cd4_bm - koff_cd3_d * complex_cd3cd4_scd38_bonemarrow) -
      (kon_cd3_d * complex_scd38_bonemarrow * u_cd3cd8_bm - koff_cd3_d * complex_cd3cd8_scd38_bonemarrow) -
      (kon_bcma_d * complex_scd38_bonemarrow * u_sbcma_bm - koff_bcma_d * complex_sbcma_scd38_bonemarrow) -
      kel * complex_scd38_bonemarrow

    d/dt(complex_cd3cd4_scd38_bonemarrow) <-
      (kon_cd38_d * drug_cd3cd4_bonemarrow * u_scd38_bm - koff_cd38_d * complex_cd3cd4_scd38_bonemarrow) +
      (kon_cd3_d * complex_scd38_bonemarrow * u_cd3cd4_bm - koff_cd3_d * complex_cd3cd4_scd38_bonemarrow) -
      kint_cd3_d * complex_cd3cd4_scd38_bonemarrow

    d/dt(complex_cd3cd8_scd38_bonemarrow) <-
      (kon_cd38_d * drug_cd3cd8_bonemarrow * u_scd38_bm - koff_cd38_d * complex_cd3cd8_scd38_bonemarrow) +
      (kon_cd3_d * complex_scd38_bonemarrow * u_cd3cd8_bm - koff_cd3_d * complex_cd3cd8_scd38_bonemarrow) -
      kint_cd3_d * complex_cd3cd8_scd38_bonemarrow

    d/dt(complex_sbcma_scd38_bonemarrow) <-
      (kon_cd38_d * complex_sbcma_bonemarrow * u_scd38_bm - koff_cd38_d * complex_sbcma_scd38_bonemarrow) +
      (kon_bcma_d * complex_scd38_bonemarrow * u_sbcma_bm - koff_bcma_d * complex_sbcma_scd38_bonemarrow) -
      kel * complex_sbcma_scd38_bonemarrow

    # ================================================================
    # 6. ODE system - remaining mPBPK compartments and the SC depot.
    # ================================================================
    d/dt(tight) <- (-(1 - sigma_l) * lflow1 * tight + (1 - sigma_tight) * lflow1 * plasma) / vtight

    d/dt(lymph) <- (lflow2 * (1 - sigma_l) * f_leaky * leaky +
                      (1 - sigma_l) * lflow1 * tight +
                      lflow2 * (1 - sigma_l) * f_bm * bonemarrow -
                      lflow * lymph) / vlymph

    d/dt(leaky) <- (-lflow2 * (1 - sigma_l) * f_leaky * leaky +
                      (1 - sigma_leaky) * lflow2 * f_leaky * plasma) / vleaky_wo_bm

    # The SC depot loses drug both to systemic absorption (Ka*BA) and to
    # presystemic loss (Ka*(1-BA)), so it empties at Ka regardless of BA.
    d/dt(depot) <- -ka * (1 - ba) * depot - ka * ba * depot

    # ================================================================
    # 7. Active species and the nACT exposure metric.
    #    Data S1 assignment rules: the trimers and tetramers that
    #    crosslink a T cell to a tumour cell are summed per T-cell subset
    #    and normalised per tumour cell. The 1e-9 converts nM to mol/L
    #    before multiplying by Avogadro's number.
    # ================================================================
    act_cd4_bonemarrow <- trimer_cd3cd4_bcma_bonemarrow + trimer_cd3cd4_cd38_bonemarrow +
      tetramer_cd3cd4_bonemarrow + complex_cd3cd4_sbcma_cd38_bonemarrow
    act_cd8_bonemarrow <- trimer_cd3cd8_bcma_bonemarrow + trimer_cd3cd8_cd38_bonemarrow +
      tetramer_cd3cd8_bonemarrow + complex_cd3cd8_sbcma_cd38_bonemarrow

    # Molecules of active species per bone-marrow tumour (plasma) cell.
    # This is the metric the paper used for FIH and efficacy dose
    # prediction: nACT EC50 = 3.1 and EC90 = 28 molecule/cell from the
    # in-vivo mouse PKPD model (Table 2).
    nact <- (act_cd4_bonemarrow + act_cd8_bonemarrow) * 1e-9 * vbonemarrow *
      avogadro / n_tumor_bonemarrow
    # Per-T-cell normalisations (Data S1 ACTperCD4_bm / ACTperCD8_bm).
    nact_cd4 <- act_cd4_bonemarrow * 1e-9 * vbonemarrow * avogadro / n_cd4_bonemarrow
    nact_cd8 <- act_cd8_bonemarrow * 1e-9 * vbonemarrow * avogadro / n_cd8_bonemarrow

    # ================================================================
    # 8. Observations.
    #    The clinical assay quantifies FREE ISB 2001 in serum, which is
    #    the free plasma state (Methods, "Free ISB 2001 concentrations in
    #    the serum samples were quantified using a validated
    #    electrochemiluminescence-based immunoassay").
    # ================================================================
    Cc          <- plasma                      # free ISB 2001 in serum (nM)
    Cbonemarrow <- bonemarrow                  # free ISB 2001 in bone-marrow interstitial fluid (nM)

    Cc ~ prop(propSd)
  })
}
