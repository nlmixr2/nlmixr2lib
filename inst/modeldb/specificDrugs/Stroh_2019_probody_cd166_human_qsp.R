Stroh_2019_probody_cd166_human_qsp <- function() {
  description <- paste(
    "QSP. Human projection of the CytomX quantitative systems pharmacology model",
    "for a PROBODY therapeutic (Pb-Tx) directed against CD166 (ALCAM). A Pb-Tx is",
    "a monoclonal antibody whose two paratopes each carry a peptide mask tethered",
    "by a tumour-protease-cleavable substrate. The model resolves every",
    "combination of the two arms being mask-closed (c), mask-open / 'breathing'",
    "(o) or cleaved (m), free or CD166-bound, in three compartments (plasma,",
    "peripheral tissue, tumour): 15 Pb-Tx species plus free CD166 per compartment,",
    "54 states in all. Reversible breathing is parameterised by the fold-masking",
    "ratio Kmask = kclose / kopen; substrate proteolysis is the pseudo-first-order",
    "kcleave, raised in tumour by fcleave_tumor to represent elevated",
    "tumour-associated protease activity. Full (non-approximated) TMDD: monovalent",
    "and bivalent CD166 binding with receptor turnover and complex endocytosis in",
    "peripheral tissue and tumour, first-order elimination of free species from",
    "plasma and periphery only. The human parameter set is allometrically scaled",
    "and assumption-driven (Table 2) rather than fitted; the cynomolgus-monkey",
    "calibration is the sibling model Stroh_2019_probody_cd166_monkey_qsp.",
    "Deterministic mechanism model: no IIV and no residual-error model, matching",
    "the source.",
    sep = " "
  )
  reference <- paste(
    "Stroh M, Sagert J, Burke JM, Apgar JF, Lin L, Millard BL, Kavanaugh WM.",
    "Quantitative systems pharmacology model of a masked, tumor-activated antibody.",
    "CPT Pharmacometrics Syst Pharmacol. 2019;8(9):676-684. doi:10.1002/psp4.12448.",
    "Structural ODE system (the k / A1 / A2 matrices of the KroneckerBio model)",
    "from Supporting Information 'Model code' (PSP4-8-676-s003).",
    sep = " "
  )
  vignette <- "Stroh_2019_probody_cd166_qsp"

  units <- list(time = "s", dosing = "nmol", concentration = "nM")

  covariateData <- list()

  # Every state is a paper-mechanistic Pb-Tx species or a CD166 pool; none maps
  # onto a canonical PK compartment name.
  paper_specific_compartment_pattern <- "^(pb_|target_)"

  compartmentData <- list(
    target_central = list(analyte = "CD166 (ALCAM)", units = "nmol", specimen = "plasma", verified = TRUE),
    pb_c_c_central = list(
      analyte = "anti-CD166 PROBODY therapeutic",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    pb_c_o_central = list(
      analyte = "anti-CD166 PROBODY therapeutic",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    pb_o_o_central = list(
      analyte = "anti-CD166 PROBODY therapeutic",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    pb_c_m_central = list(
      analyte = "anti-CD166 PROBODY therapeutic",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    pb_o_m_central = list(
      analyte = "anti-CD166 PROBODY therapeutic",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    pb_m_m_central = list(
      analyte = "anti-CD166 PROBODY therapeutic",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    pb_c_oR_central = list(
      analyte = "anti-CD166 PROBODY therapeutic-CD166 complex",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    pb_o_oR_central = list(
      analyte = "anti-CD166 PROBODY therapeutic-CD166 complex",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    pb_oR_oR_central = list(
      analyte = "anti-CD166 PROBODY therapeutic-CD166 complex",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    pb_c_mR_central = list(
      analyte = "anti-CD166 PROBODY therapeutic-CD166 complex",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    pb_o_mR_central = list(
      analyte = "anti-CD166 PROBODY therapeutic-CD166 complex",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    pb_oR_m_central = list(
      analyte = "anti-CD166 PROBODY therapeutic-CD166 complex",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    pb_oR_mR_central = list(
      analyte = "anti-CD166 PROBODY therapeutic-CD166 complex",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    pb_m_mR_central = list(
      analyte = "anti-CD166 PROBODY therapeutic-CD166 complex",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    pb_mR_mR_central = list(
      analyte = "anti-CD166 PROBODY therapeutic-CD166 complex",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    pb_el_central = list(
      analyte = "anti-CD166 PROBODY therapeutic (cumulative eliminated)",
      units = "nmol",
      specimen = "not applicable",
      verified = TRUE
    ),
    target_peripheral = list(analyte = "CD166 (ALCAM)", units = "nmol", specimen = "tissue", verified = TRUE),
    pb_c_c_peripheral = list(
      analyte = "anti-CD166 PROBODY therapeutic",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    pb_c_o_peripheral = list(
      analyte = "anti-CD166 PROBODY therapeutic",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    pb_o_o_peripheral = list(
      analyte = "anti-CD166 PROBODY therapeutic",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    pb_c_m_peripheral = list(
      analyte = "anti-CD166 PROBODY therapeutic",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    pb_o_m_peripheral = list(
      analyte = "anti-CD166 PROBODY therapeutic",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    pb_m_m_peripheral = list(
      analyte = "anti-CD166 PROBODY therapeutic",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    pb_c_oR_peripheral = list(
      analyte = "anti-CD166 PROBODY therapeutic-CD166 complex",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    pb_o_oR_peripheral = list(
      analyte = "anti-CD166 PROBODY therapeutic-CD166 complex",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    pb_oR_oR_peripheral = list(
      analyte = "anti-CD166 PROBODY therapeutic-CD166 complex",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    pb_c_mR_peripheral = list(
      analyte = "anti-CD166 PROBODY therapeutic-CD166 complex",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    pb_o_mR_peripheral = list(
      analyte = "anti-CD166 PROBODY therapeutic-CD166 complex",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    pb_oR_m_peripheral = list(
      analyte = "anti-CD166 PROBODY therapeutic-CD166 complex",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    pb_oR_mR_peripheral = list(
      analyte = "anti-CD166 PROBODY therapeutic-CD166 complex",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    pb_m_mR_peripheral = list(
      analyte = "anti-CD166 PROBODY therapeutic-CD166 complex",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    pb_mR_mR_peripheral = list(
      analyte = "anti-CD166 PROBODY therapeutic-CD166 complex",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    pb_el_peripheral = list(
      analyte = "anti-CD166 PROBODY therapeutic (cumulative eliminated)",
      units = "nmol",
      specimen = "not applicable",
      verified = TRUE
    ),
    pb_tmdd_o_peripheral = list(
      analyte = "anti-CD166 PROBODY therapeutic (cumulative receptor-mediated uptake via a masked/open arm)",
      units = "nmol",
      specimen = "not applicable",
      verified = TRUE
    ),
    pb_tmdd_m_peripheral = list(
      analyte = "anti-CD166 PROBODY therapeutic (cumulative receptor-mediated uptake via a cleaved arm)",
      units = "nmol",
      specimen = "not applicable",
      verified = TRUE
    ),
    target_tumor = list(analyte = "CD166 (ALCAM)", units = "nmol", specimen = "tumor", verified = TRUE),
    pb_c_c_tumor = list(
      analyte = "anti-CD166 PROBODY therapeutic",
      units = "nmol",
      specimen = "tumor",
      verified = TRUE
    ),
    pb_c_o_tumor = list(
      analyte = "anti-CD166 PROBODY therapeutic",
      units = "nmol",
      specimen = "tumor",
      verified = TRUE
    ),
    pb_o_o_tumor = list(
      analyte = "anti-CD166 PROBODY therapeutic",
      units = "nmol",
      specimen = "tumor",
      verified = TRUE
    ),
    pb_c_m_tumor = list(
      analyte = "anti-CD166 PROBODY therapeutic",
      units = "nmol",
      specimen = "tumor",
      verified = TRUE
    ),
    pb_o_m_tumor = list(
      analyte = "anti-CD166 PROBODY therapeutic",
      units = "nmol",
      specimen = "tumor",
      verified = TRUE
    ),
    pb_m_m_tumor = list(
      analyte = "anti-CD166 PROBODY therapeutic",
      units = "nmol",
      specimen = "tumor",
      verified = TRUE
    ),
    pb_c_oR_tumor = list(
      analyte = "anti-CD166 PROBODY therapeutic-CD166 complex",
      units = "nmol",
      specimen = "tumor",
      verified = TRUE
    ),
    pb_o_oR_tumor = list(
      analyte = "anti-CD166 PROBODY therapeutic-CD166 complex",
      units = "nmol",
      specimen = "tumor",
      verified = TRUE
    ),
    pb_oR_oR_tumor = list(
      analyte = "anti-CD166 PROBODY therapeutic-CD166 complex",
      units = "nmol",
      specimen = "tumor",
      verified = TRUE
    ),
    pb_c_mR_tumor = list(
      analyte = "anti-CD166 PROBODY therapeutic-CD166 complex",
      units = "nmol",
      specimen = "tumor",
      verified = TRUE
    ),
    pb_o_mR_tumor = list(
      analyte = "anti-CD166 PROBODY therapeutic-CD166 complex",
      units = "nmol",
      specimen = "tumor",
      verified = TRUE
    ),
    pb_oR_m_tumor = list(
      analyte = "anti-CD166 PROBODY therapeutic-CD166 complex",
      units = "nmol",
      specimen = "tumor",
      verified = TRUE
    ),
    pb_oR_mR_tumor = list(
      analyte = "anti-CD166 PROBODY therapeutic-CD166 complex",
      units = "nmol",
      specimen = "tumor",
      verified = TRUE
    ),
    pb_m_mR_tumor = list(
      analyte = "anti-CD166 PROBODY therapeutic-CD166 complex",
      units = "nmol",
      specimen = "tumor",
      verified = TRUE
    ),
    pb_mR_mR_tumor = list(
      analyte = "anti-CD166 PROBODY therapeutic-CD166 complex",
      units = "nmol",
      specimen = "tumor",
      verified = TRUE
    ),
    pb_tmdd_o_tumor = list(
      analyte = "anti-CD166 PROBODY therapeutic (cumulative receptor-mediated uptake via a masked/open arm)",
      units = "nmol",
      specimen = "not applicable",
      verified = TRUE
    ),
    pb_tmdd_m_tumor = list(
      analyte = "anti-CD166 PROBODY therapeutic (cumulative receptor-mediated uptake via a cleaved arm)",
      units = "nmol",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human (in silico projection; no clinical data were fitted)",
    n_subjects = 1L,
    disease_state = "cancer patient carrying a solid tumour; CD166 is expressed on both tumour and healthy tissue",
    dose_range = paste(
      "A single 4.5 mg/kg dose for the plasma / peripheral / tumour projections of",
      "Figure 4 and multiple 3 mg/kg doses for the intact-versus-cleaved",
      "projection of Figure 5. Doses are entered into this model in nmol; see the",
      "vignette Errata for why mg/kg cannot be converted from on-disk sources.",
      sep = " "
    ),
    notes = paste(
      "The human model is a forward projection, not a fit: it carries the",
      "cynomolgus-monkey elementary rate constants unchanged (kendo, kon1, kon2,",
      "koff1, Kmask, kcleave), allometrically scales kel, k12 and k21, and assumes",
      "a 0.01 L breast tumour and a 70 kg body weight (Table 2). No human subjects",
      "contributed data, so n_subjects records the single typical individual the",
      "projection describes.",
      sep = " "
    )
  )

  ini({
    # ---- physiologic volumes (Table 2) ----
    lvc <- fixed(log(2.6)); label("Plasma volume V1 (L)")  # Table 2, 'Plasma volume, V 1 (L)' = 2.6, from Davies & Morris 1993
    lvtumor <- fixed(log(0.01)); label("Tumour volume V3 (L)")  # Table 2, 'Tumor volume, V 3 (L)' = 0.01, based on a breast tumour

    # ---- disposition (Table 2) ----
    lkel <- fixed(log(6.3e-7)); label("First-order elimination rate constant of free Pb-Tx from plasma and periphery (1/s)")  # Table 2, allometric scaling of the monkey value: kel * (BW_human/BW_monkey)^0.85 * (V_monkey/V_human)
    lk12 <- fixed(log(4.8e-6)); label("Plasma to peripheral transport rate constant (1/s)")  # Table 2, allometric scaling: k12 * (BW_human/BW_monkey)^(-0.25)
    lk21 <- fixed(log(4.4e-6)); label("Peripheral to plasma transport rate constant (1/s)")  # Table 2, allometric scaling: k21 * (BW_human/BW_monkey)^(-0.25)
    lk13 <- fixed(log(1.9e-8)); label("Plasma to tumour transport rate constant (1/s)")  # Table 2, derived as k13 = Q * p / (p + V1/V3) with Q = 1e-5 1/s and p = 0.5
    lk31 <- fixed(log(1.0e-5)); label("Tumour to plasma transport rate constant (1/s)")  # Table 2, derived as k31 = Q / (1 + p * V3/V1) with Q = 1e-5 1/s and p = 0.5

    # ---- CD166 turnover and binding (Table 2) ----
    lkint <- fixed(log(1.0e-4)); label("CD166 endocytosis rate constant, free receptor and complex (1/s)")  # Table 2, 'Target endocytosis rate constant, k endo' = 1.0e-4, same as monkey
    ksynr_central <- fixed(0); label("CD166 synthesis rate in plasma (nmol/s)")  # Figure 2b draws no CD166 in the plasma compartment; the deposited code carries the plasma binding machinery but never synthesises plasma receptor, so it is inert. Not log-transformed because log(0) is undefined.
    lksynr_peripheral <- fixed(log(2.3e-3)); label("CD166 synthesis rate in peripheral tissue (nmol/s)")  # Table 2, 'Target synthesis rate, k synR (nmol/second) peripheral' = 2.3e-3, from ksynR = kendo * R_T
    lksynr_tumor <- fixed(log(2.7e-4)); label("CD166 synthesis rate in tumour (nmol/s)")  # Table 2, 'Target synthesis rate, k synR (nmol/second) tumor' = 2.7e-4, from ksynR = kendo * R_T
    lkon1 <- fixed(log(1e-3)); label("Monovalent CD166 association rate constant (1/(nM*s))")  # Table 2, 'Forward binding rate, k on1' = 1e-3, assumed the same as monkey
    lkon2 <- fixed(log(1e-3)); label("Second-arm (bivalent) CD166 association rate constant (1/(nM*s))")  # Table 2, 'Forward binding rate, k on2' = 1e-3, assumed the same as monkey
    lkoff1 <- fixed(log(2e-3)); label("CD166 dissociation rate constant per bound arm (1/s)")  # Table 2, 'Reverse binding rate constant, k off1' = 2e-3, assumed the same as monkey

    # ---- PROBODY mask and substrate ----
    lkopen <- fixed(log(1.16e-2)); label("Mask opening ('breathing') rate constant per closed arm (1/s)")  # not in Table 1 or Table 2: taken from the Supporting Information 'Model code' A1 matrix (entry (3,2) = 2.32e-2 = 2 * kopen). Corroborated by Ippolito 2024 Table S2, which cites Stroh 2019 for 'Probody unmasking rate' k_o = 0.0116 1/second.
    lkmask <- fixed(log(220)); label("Fold-masking Kmask = kclose / kopen for mask M1 (unitless)")  # Table 2, 'Fold-masking, K mask' = 220 (M1); the weaker mask M2 is 57. Assumed the same as monkey.
    lkcleave <- fixed(log(3e-7)); label("Substrate cleavage rate constant per uncleaved arm in plasma and periphery (1/s)")  # Table 2, 'Rate constant for cleavage reaction, k cleave' < 3e-7. Results: the marginal probability gave no lower bound, and 3e-7 'was carried forward as the likely maximum value'.
    lfcleave_tumor <- fixed(log(10)); label("Fold increase in kcleave inside the tumour (unitless)")  # Discussion: 'Under the scenario of a 10-fold increased k cleave in the tumor relative to the periphery, Figure 4d ...'. The paper scans this factor; 10 is the scenario Figure 4d,e reports.
  })

  model({
    # Pb-Tx species are named exactly as in the Supporting Information 'Model code'
    # state list. Each Pb-Tx carries two arms; an arm is in one of three states:
    #   c  = mask closed   -- masked, cannot bind CD166
    #   o  = mask open     -- 'breathing' open, binding-competent, mask still attached
    #   m  = mAb arm       -- substrate cleaved, mask shed, binding-competent
    # A trailing R marks an arm bound to CD166 (oR, mR). Arms are unordered, so the
    # six free species are c_c, c_o, o_o, c_m, o_m, m_m and the nine bound species are
    # c_oR, o_oR, oR_oR, c_mR, o_mR, oR_m, oR_mR, m_mR, mR_mR (Figure 2a).

    vc <- exp(lvc)
    vtumor <- exp(lvtumor)
    kel <- exp(lkel)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    k13 <- exp(lk13)
    k31 <- exp(lk31)
    kint <- exp(lkint)
    ksynr_peripheral <- exp(lksynr_peripheral)
    ksynr_tumor <- exp(lksynr_tumor)
    kon1 <- exp(lkon1)
    kon2 <- exp(lkon2)
    koff1 <- exp(lkoff1)
    kopen <- exp(lkopen)
    kmask <- exp(lkmask)
    kclose <- kmask * kopen
    kcleave <- exp(lkcleave)
    kcleave_tumor <- exp(lfcleave_tumor) * kcleave

    # CD166 sits at its synthesis / endocytosis steady state before the first dose.
    target_central(0) <- ksynr_central / kint
    target_peripheral(0) <- ksynr_peripheral / kint
    target_tumor(0) <- ksynr_tumor / kint

    # ---- plasma compartment (deposited-code states 1-17) ----
    d/dt(target_central) <-
      ksynr_central - kint * target_central +
      koff1 * pb_c_oR_central + koff1 * pb_o_oR_central + koff1 * pb_c_mR_central + koff1 * pb_o_mR_central + koff1 * pb_oR_m_central + koff1 * pb_m_mR_central +
      2 * koff1 * pb_oR_oR_central + 2 * koff1 * pb_oR_mR_central + 2 * koff1 * pb_mR_mR_central -
      kon1 * target_central * (pb_c_o_central + 2 * pb_o_o_central + pb_c_m_central + 2 * pb_o_m_central + 2 * pb_m_m_central) -
      kon2 * target_central * (pb_o_oR_central + pb_o_mR_central + pb_oR_m_central + pb_m_mR_central)
    d/dt(pb_c_c_central) <-
      kclose * pb_c_o_central + k21 * pb_c_c_peripheral + k31 * pb_c_c_tumor -
      (2 * kopen + 2 * kcleave + kel + k12 + k13) * pb_c_c_central
    d/dt(pb_c_o_central) <-
      2 * kopen * pb_c_c_central + 2 * kclose * pb_o_o_central + koff1 * pb_c_oR_central + k21 * pb_c_o_peripheral + k31 * pb_c_o_tumor -
      (kclose + kopen + 2 * kcleave + kon1 * target_central + kel + k12 + k13) * pb_c_o_central
    d/dt(pb_o_o_central) <-
      kopen * pb_c_o_central + koff1 * pb_o_oR_central + k21 * pb_o_o_peripheral + k31 * pb_o_o_tumor -
      (2 * kclose + 2 * kcleave + 2 * kon1 * target_central + kel + k12 + k13) * pb_o_o_central
    d/dt(pb_c_m_central) <-
      2 * kcleave * pb_c_c_central + kcleave * pb_c_o_central + kclose * pb_o_m_central + koff1 * pb_c_mR_central + k21 * pb_c_m_peripheral + k31 * pb_c_m_tumor -
      (kopen + kcleave + kon1 * target_central + kel + k12 + k13) * pb_c_m_central
    d/dt(pb_o_m_central) <-
      kcleave * pb_c_o_central + 2 * kcleave * pb_o_o_central + kopen * pb_c_m_central + koff1 * (pb_o_mR_central + pb_oR_m_central) + k21 * pb_o_m_peripheral + k31 * pb_o_m_tumor -
      (kclose + kcleave + 2 * kon1 * target_central + kel + k12 + k13) * pb_o_m_central
    d/dt(pb_m_m_central) <-
      kcleave * (pb_c_m_central + pb_o_m_central) + koff1 * pb_m_mR_central + k21 * pb_m_m_peripheral + k31 * pb_m_m_tumor -
      (2 * kon1 * target_central + kel + k12 + k13) * pb_m_m_central
    d/dt(pb_c_oR_central) <-
      kon1 * target_central * pb_c_o_central + kclose * pb_o_oR_central -
      (koff1 + kopen + 2 * kcleave + kint) * pb_c_oR_central
    d/dt(pb_o_oR_central) <-
      kopen * pb_c_oR_central + 2 * kon1 * target_central * pb_o_o_central + 2 * koff1 * pb_oR_oR_central -
      (koff1 + kclose + 2 * kcleave + kon2 * target_central + kint) * pb_o_oR_central
    d/dt(pb_oR_oR_central) <-
      kon2 * target_central * pb_o_oR_central -
      (2 * koff1 + 2 * kcleave + kint) * pb_oR_oR_central
    d/dt(pb_c_mR_central) <-
      kon1 * target_central * pb_c_m_central + kcleave * pb_c_oR_central + kclose * pb_o_mR_central -
      (koff1 + kopen + kcleave + kint) * pb_c_mR_central
    d/dt(pb_o_mR_central) <-
      kon1 * target_central * pb_o_m_central + kcleave * pb_o_oR_central + kopen * pb_c_mR_central + koff1 * pb_oR_mR_central -
      (koff1 + kclose + kcleave + kon2 * target_central + kint) * pb_o_mR_central
    d/dt(pb_oR_m_central) <-
      kon1 * target_central * pb_o_m_central + kcleave * (pb_c_oR_central + pb_o_oR_central) + koff1 * pb_oR_mR_central -
      (koff1 + kcleave + kon2 * target_central + kint) * pb_oR_m_central
    d/dt(pb_oR_mR_central) <-
      kon2 * target_central * (pb_o_mR_central + pb_oR_m_central) + 2 * kcleave * pb_oR_oR_central -
      (2 * koff1 + kcleave + kint) * pb_oR_mR_central
    d/dt(pb_m_mR_central) <-
      2 * kon1 * target_central * pb_m_m_central + kcleave * (pb_c_mR_central + pb_o_mR_central + pb_oR_m_central) + 2 * koff1 * pb_mR_mR_central -
      (koff1 + kon2 * target_central + kint) * pb_m_mR_central
    d/dt(pb_mR_mR_central) <-
      kon2 * target_central * pb_m_mR_central + kcleave * pb_oR_mR_central -
      (2 * koff1 + kint) * pb_mR_mR_central
    d/dt(pb_el_central) <- kel * (pb_c_c_central + pb_c_o_central + pb_o_o_central + pb_c_m_central + pb_o_m_central + pb_m_m_central)

    # ---- peripheral-tissue compartment (deposited-code states 18-36) ----
    d/dt(target_peripheral) <-
      ksynr_peripheral - kint * target_peripheral +
      koff1 * pb_c_oR_peripheral + koff1 * pb_o_oR_peripheral + koff1 * pb_c_mR_peripheral + koff1 * pb_o_mR_peripheral + koff1 * pb_oR_m_peripheral + koff1 * pb_m_mR_peripheral +
      2 * koff1 * pb_oR_oR_peripheral + 2 * koff1 * pb_oR_mR_peripheral + 2 * koff1 * pb_mR_mR_peripheral -
      kon1 * target_peripheral * (pb_c_o_peripheral + 2 * pb_o_o_peripheral + pb_c_m_peripheral + 2 * pb_o_m_peripheral + 2 * pb_m_m_peripheral) -
      kon2 * target_peripheral * (pb_o_oR_peripheral + pb_o_mR_peripheral + pb_oR_m_peripheral + pb_m_mR_peripheral)
    d/dt(pb_c_c_peripheral) <-
      kclose * pb_c_o_peripheral + k12 * pb_c_c_central -
      (2 * kopen + 2 * kcleave + kel + k21) * pb_c_c_peripheral
    d/dt(pb_c_o_peripheral) <-
      2 * kopen * pb_c_c_peripheral + 2 * kclose * pb_o_o_peripheral + koff1 * pb_c_oR_peripheral + k12 * pb_c_o_central -
      (kclose + kopen + 2 * kcleave + kon1 * target_peripheral + kel + k21) * pb_c_o_peripheral
    d/dt(pb_o_o_peripheral) <-
      kopen * pb_c_o_peripheral + koff1 * pb_o_oR_peripheral + k12 * pb_o_o_central -
      (2 * kclose + 2 * kcleave + 2 * kon1 * target_peripheral + kel + k21) * pb_o_o_peripheral
    d/dt(pb_c_m_peripheral) <-
      2 * kcleave * pb_c_c_peripheral + kcleave * pb_c_o_peripheral + kclose * pb_o_m_peripheral + koff1 * pb_c_mR_peripheral + k12 * pb_c_m_central -
      (kopen + kcleave + kon1 * target_peripheral + kel + k21) * pb_c_m_peripheral
    d/dt(pb_o_m_peripheral) <-
      kcleave * pb_c_o_peripheral + 2 * kcleave * pb_o_o_peripheral + kopen * pb_c_m_peripheral + koff1 * (pb_o_mR_peripheral + pb_oR_m_peripheral) + k12 * pb_o_m_central -
      (kclose + kcleave + 2 * kon1 * target_peripheral + kel + k21) * pb_o_m_peripheral
    d/dt(pb_m_m_peripheral) <-
      kcleave * (pb_c_m_peripheral + pb_o_m_peripheral) + koff1 * pb_m_mR_peripheral + k12 * pb_m_m_central -
      (2 * kon1 * target_peripheral + kel + k21) * pb_m_m_peripheral
    d/dt(pb_c_oR_peripheral) <-
      kon1 * target_peripheral * pb_c_o_peripheral + kclose * pb_o_oR_peripheral -
      (koff1 + kopen + 2 * kcleave + kint) * pb_c_oR_peripheral
    d/dt(pb_o_oR_peripheral) <-
      kopen * pb_c_oR_peripheral + 2 * kon1 * target_peripheral * pb_o_o_peripheral + 2 * koff1 * pb_oR_oR_peripheral -
      (koff1 + kclose + 2 * kcleave + kon2 * target_peripheral + kint) * pb_o_oR_peripheral
    d/dt(pb_oR_oR_peripheral) <-
      kon2 * target_peripheral * pb_o_oR_peripheral -
      (2 * koff1 + 2 * kcleave + kint) * pb_oR_oR_peripheral
    d/dt(pb_c_mR_peripheral) <-
      kon1 * target_peripheral * pb_c_m_peripheral + kcleave * pb_c_oR_peripheral + kclose * pb_o_mR_peripheral -
      (koff1 + kopen + kcleave + kint) * pb_c_mR_peripheral
    d/dt(pb_o_mR_peripheral) <-
      kon1 * target_peripheral * pb_o_m_peripheral + kcleave * pb_o_oR_peripheral + kopen * pb_c_mR_peripheral + koff1 * pb_oR_mR_peripheral -
      (koff1 + kclose + kcleave + kon2 * target_peripheral + kint) * pb_o_mR_peripheral
    d/dt(pb_oR_m_peripheral) <-
      kon1 * target_peripheral * pb_o_m_peripheral + kcleave * (pb_c_oR_peripheral + pb_o_oR_peripheral) + koff1 * pb_oR_mR_peripheral -
      (koff1 + kcleave + kon2 * target_peripheral + kint) * pb_oR_m_peripheral
    d/dt(pb_oR_mR_peripheral) <-
      kon2 * target_peripheral * (pb_o_mR_peripheral + pb_oR_m_peripheral) + 2 * kcleave * pb_oR_oR_peripheral -
      (2 * koff1 + kcleave + kint) * pb_oR_mR_peripheral
    d/dt(pb_m_mR_peripheral) <-
      2 * kon1 * target_peripheral * pb_m_m_peripheral + kcleave * (pb_c_mR_peripheral + pb_o_mR_peripheral + pb_oR_m_peripheral) + 2 * koff1 * pb_mR_mR_peripheral -
      (koff1 + kon2 * target_peripheral + kint) * pb_m_mR_peripheral
    d/dt(pb_mR_mR_peripheral) <-
      kon2 * target_peripheral * pb_m_mR_peripheral + kcleave * pb_oR_mR_peripheral -
      (2 * koff1 + kint) * pb_mR_mR_peripheral
    d/dt(pb_tmdd_m_peripheral) <- kint * (pb_c_mR_peripheral + pb_o_mR_peripheral + pb_oR_mR_peripheral + pb_m_mR_peripheral + pb_mR_mR_peripheral)
    d/dt(pb_tmdd_o_peripheral) <- kint * (pb_c_oR_peripheral + pb_o_oR_peripheral + pb_oR_oR_peripheral + pb_oR_m_peripheral)
    d/dt(pb_el_peripheral) <- kel * (pb_c_c_peripheral + pb_c_o_peripheral + pb_o_o_peripheral + pb_c_m_peripheral + pb_o_m_peripheral + pb_m_m_peripheral)

    # ---- tumour compartment (deposited-code states 37-54) ----
    d/dt(target_tumor) <-
      ksynr_tumor - kint * target_tumor +
      koff1 * pb_c_oR_tumor + koff1 * pb_o_oR_tumor + koff1 * pb_c_mR_tumor + koff1 * pb_o_mR_tumor + koff1 * pb_oR_m_tumor + koff1 * pb_m_mR_tumor +
      2 * koff1 * pb_oR_oR_tumor + 2 * koff1 * pb_oR_mR_tumor + 2 * koff1 * pb_mR_mR_tumor -
      kon1 * target_tumor * (pb_c_o_tumor + 2 * pb_o_o_tumor + pb_c_m_tumor + 2 * pb_o_m_tumor + 2 * pb_m_m_tumor) -
      kon2 * target_tumor * (pb_o_oR_tumor + pb_o_mR_tumor + pb_oR_m_tumor + pb_m_mR_tumor)
    d/dt(pb_c_c_tumor) <-
      kclose * pb_c_o_tumor + k13 * pb_c_c_central -
      (2 * kopen + 2 * kcleave_tumor + k31) * pb_c_c_tumor
    d/dt(pb_c_o_tumor) <-
      2 * kopen * pb_c_c_tumor + 2 * kclose * pb_o_o_tumor + koff1 * pb_c_oR_tumor + k13 * pb_c_o_central -
      (kclose + kopen + 2 * kcleave_tumor + kon1 * target_tumor + k31) * pb_c_o_tumor
    d/dt(pb_o_o_tumor) <-
      kopen * pb_c_o_tumor + koff1 * pb_o_oR_tumor + k13 * pb_o_o_central -
      (2 * kclose + 2 * kcleave_tumor + 2 * kon1 * target_tumor + k31) * pb_o_o_tumor
    d/dt(pb_c_m_tumor) <-
      2 * kcleave_tumor * pb_c_c_tumor + kcleave_tumor * pb_c_o_tumor + kclose * pb_o_m_tumor + koff1 * pb_c_mR_tumor + k13 * pb_c_m_central -
      (kopen + kcleave_tumor + kon1 * target_tumor + k31) * pb_c_m_tumor
    d/dt(pb_o_m_tumor) <-
      kcleave_tumor * pb_c_o_tumor + 2 * kcleave_tumor * pb_o_o_tumor + kopen * pb_c_m_tumor + koff1 * (pb_o_mR_tumor + pb_oR_m_tumor) + k13 * pb_o_m_central -
      (kclose + kcleave_tumor + 2 * kon1 * target_tumor + k31) * pb_o_m_tumor
    d/dt(pb_m_m_tumor) <-
      kcleave_tumor * (pb_c_m_tumor + pb_o_m_tumor) + koff1 * pb_m_mR_tumor + k13 * pb_m_m_central -
      (2 * kon1 * target_tumor + k31) * pb_m_m_tumor
    d/dt(pb_c_oR_tumor) <-
      kon1 * target_tumor * pb_c_o_tumor + kclose * pb_o_oR_tumor -
      (koff1 + kopen + 2 * kcleave_tumor + kint) * pb_c_oR_tumor
    d/dt(pb_o_oR_tumor) <-
      kopen * pb_c_oR_tumor + 2 * kon1 * target_tumor * pb_o_o_tumor + 2 * koff1 * pb_oR_oR_tumor -
      (koff1 + kclose + 2 * kcleave_tumor + kon2 * target_tumor + kint) * pb_o_oR_tumor
    d/dt(pb_oR_oR_tumor) <-
      kon2 * target_tumor * pb_o_oR_tumor -
      (2 * koff1 + 2 * kcleave_tumor + kint) * pb_oR_oR_tumor
    d/dt(pb_c_mR_tumor) <-
      kon1 * target_tumor * pb_c_m_tumor + kcleave_tumor * pb_c_oR_tumor + kclose * pb_o_mR_tumor -
      (koff1 + kopen + kcleave_tumor + kint) * pb_c_mR_tumor
    d/dt(pb_o_mR_tumor) <-
      kon1 * target_tumor * pb_o_m_tumor + kcleave_tumor * pb_o_oR_tumor + kopen * pb_c_mR_tumor + koff1 * pb_oR_mR_tumor -
      (koff1 + kclose + kcleave_tumor + kon2 * target_tumor + kint) * pb_o_mR_tumor
    d/dt(pb_oR_m_tumor) <-
      kon1 * target_tumor * pb_o_m_tumor + kcleave_tumor * (pb_c_oR_tumor + pb_o_oR_tumor) + koff1 * pb_oR_mR_tumor -
      (koff1 + kcleave_tumor + kon2 * target_tumor + kint) * pb_oR_m_tumor
    d/dt(pb_oR_mR_tumor) <-
      kon2 * target_tumor * (pb_o_mR_tumor + pb_oR_m_tumor) + 2 * kcleave_tumor * pb_oR_oR_tumor -
      (2 * koff1 + kcleave_tumor + kint) * pb_oR_mR_tumor
    d/dt(pb_m_mR_tumor) <-
      2 * kon1 * target_tumor * pb_m_m_tumor + kcleave_tumor * (pb_c_mR_tumor + pb_o_mR_tumor + pb_oR_m_tumor) + 2 * koff1 * pb_mR_mR_tumor -
      (koff1 + kon2 * target_tumor + kint) * pb_m_mR_tumor
    d/dt(pb_mR_mR_tumor) <-
      kon2 * target_tumor * pb_m_mR_tumor + kcleave_tumor * pb_oR_mR_tumor -
      (2 * koff1 + kint) * pb_mR_mR_tumor
    d/dt(pb_tmdd_m_tumor) <- kint * (pb_c_mR_tumor + pb_o_mR_tumor + pb_oR_mR_tumor + pb_m_mR_tumor + pb_mR_mR_tumor)
    d/dt(pb_tmdd_o_tumor) <- kint * (pb_c_oR_tumor + pb_o_oR_tumor + pb_oR_oR_tumor + pb_oR_m_tumor)

    # ---- observations (Figure 3 plots 'Total Drug (nM)'; Figure 5 splits it) ----
    Cc <- (pb_c_c_central + pb_c_o_central + pb_o_o_central + pb_c_m_central + pb_o_m_central + pb_m_m_central + pb_c_oR_central + pb_o_oR_central + pb_oR_oR_central + pb_c_mR_central + pb_o_mR_central + pb_oR_m_central + pb_oR_mR_central + pb_m_mR_central + pb_mR_mR_central) / vc
    pbIntactPlasma <- (pb_c_c_central + pb_c_o_central + pb_o_o_central + pb_c_oR_central + pb_o_oR_central + pb_oR_oR_central) / vc
    pbCleavedPlasma <- Cc - pbIntactPlasma
    pbPeripheral <- pb_c_c_peripheral + pb_c_o_peripheral + pb_o_o_peripheral + pb_c_m_peripheral + pb_o_m_peripheral + pb_m_m_peripheral + pb_c_oR_peripheral + pb_o_oR_peripheral + pb_oR_oR_peripheral + pb_c_mR_peripheral + pb_o_mR_peripheral + pb_oR_m_peripheral + pb_oR_mR_peripheral + pb_m_mR_peripheral + pb_mR_mR_peripheral
    pbTumor <- (pb_c_c_tumor + pb_c_o_tumor + pb_o_o_tumor + pb_c_m_tumor + pb_o_m_tumor + pb_m_m_tumor + pb_c_oR_tumor + pb_o_oR_tumor + pb_oR_oR_tumor + pb_c_mR_tumor + pb_o_mR_tumor + pb_oR_m_tumor + pb_oR_mR_tumor + pb_m_mR_tumor + pb_mR_mR_tumor) / vtumor
    uptakePeripheral <- pb_tmdd_m_peripheral + pb_tmdd_o_peripheral
    uptakeTumor <- pb_tmdd_m_tumor + pb_tmdd_o_tumor
  })
}
