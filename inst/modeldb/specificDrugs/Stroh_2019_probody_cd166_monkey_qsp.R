Stroh_2019_probody_cd166_monkey_qsp <- function() {
  description <- paste(
    "QSP. Preclinical (cynomolgus monkey). CytomX quantitative systems",
    "pharmacology model for a PROBODY therapeutic (Pb-Tx) directed against CD166",
    "(ALCAM), calibrated to plasma PK of one parental antibody and five Pb-Tx of",
    "differing mask strength and substrate cleavability. A Pb-Tx is a monoclonal",
    "antibody whose two paratopes each carry a peptide mask tethered by a",
    "protease-cleavable substrate. The model resolves every combination of the two",
    "arms being mask-closed (c), mask-open / 'breathing' (o) or cleaved (m), free",
    "or CD166-bound, in two compartments (plasma, peripheral tissue): 15 Pb-Tx",
    "species plus free CD166 per compartment, 36 states in all. Reversible",
    "breathing is parameterised by the fold-masking ratio Kmask = kclose / kopen;",
    "substrate proteolysis is the pseudo-first-order kcleave. Full",
    "(non-approximated) TMDD: monovalent and bivalent CD166 binding with receptor",
    "turnover and complex endocytosis in peripheral tissue, first-order",
    "elimination of free species from plasma and periphery. Non-tumour-bearing",
    "monkeys were studied, so this model carries no tumour compartment (Table 1",
    "reports no tumour volume, perfusion or partition parameters); the",
    "tumour-bearing human projection is the sibling model",
    "Stroh_2019_probody_cd166_human_qsp. Deterministic mechanism model: no IIV and",
    "no residual-error model, matching the source.",
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
    )
  )

  population <- list(
    species = "cynomolgus monkey (Macaca fascicularis)",
    n_subjects = 20L,
    n_studies = 1L,
    weight_median = "2.6 kg",
    disease_state = "experimentally naive, non-tumour-bearing",
    dose_range = "3, 5 and 10 mg/kg intravenous slow bolus, as a single dose or as two doses three weeks apart",
    regions = "Charles River Laboratories",
    notes = paste(
      "Methods report n = 2 per dosing group and Figure 3 shows ten",
      "molecule-by-dose groups (four at 3 mg/kg, three at 5 mg/kg, three at",
      "10 mg/kg), so n_subjects = 20 assumes no animal was re-used across groups;",
      "the paper does not state a total. Six molecules were fitted",
      "simultaneously: the parental antibody mAb(0,0) and five Pb-Tx built from it",
      "with masks M1 / M2 and substrates S1 / S2 (Pb-Tx M1,S1; M2,S1; M1,S2; M2,S2)",
      "plus Pb-Tx M1,0, which has no protease-activatable substrate. Physiologic",
      "parameters and the parameters describing shared molecular features were",
      "held common across all six; only Kmask and kcleave distinguish them. Whole",
      "blood was sampled for up to 21 days post dose and Pb-Tx and mAb were assayed",
      "by sandwich colorimetric ELISA.",
      sep = " "
    )
  )

  ini({
    # ---- physiologic volumes (Table 1) ----
    lvc <- fixed(log(0.1)); label("Plasma volume V1 (L)")  # Table 1, 'Plasma and peripheral volumes, V 1 and V 2 , respectively (L)' = 0.1, from Davies & Morris 1993
    lvp <- fixed(log(0.1)); label("Peripheral volume V2 (L)")  # Table 1, same row as V1 = 0.1

    # ---- disposition (Table 1) ----
    lkel <- log(1.0e-6); label("First-order elimination rate constant of free Pb-Tx (1/s)")  # Table 1, 'First-order elimination rate, k el' = 1.0e-6, 'From fit to beta-phase PK data'
    lk12 <- log(1.1e-5); label("Plasma to peripheral transport rate constant (1/s)")  # Table 1, 'Intercompartment transport rate constant, k 12' = 1.1e-5, 'From fit to alpha-phase PK'
    lk21 <- log(1.0e-5); label("Peripheral to plasma transport rate constant (1/s)")  # Table 1, 'Intercompartment transport rate constant, k 21' = 1.0e-5, same fit

    # ---- CD166 turnover and binding (Table 1) ----
    lkint <- fixed(log(1.0e-4)); label("CD166 endocytosis rate constant, free receptor and complex (1/s)")  # Table 1, 'Target endocytosis rate constant, k endo' = 1.0e-4, 'Typical value for receptor turnover rate' (Lauffenburger & Linderman 1993)
    ksynr_central <- fixed(0); label("CD166 synthesis rate in plasma (nmol/s)")  # Results: 'model provisions were included for binding to target within the peripheral and tumor compartments'; Figure 2b draws no CD166 in plasma. The binding machinery is carried in plasma but is inert. Not log-transformed because log(0) is undefined.
    lksynr_peripheral <- log(9.0e-5); label("CD166 synthesis rate in peripheral tissue (nmol/s)")  # Table 1, 'Target synthesis rate, k synR (nmol/second)' = 9.0e-5, 'estimated directly by fitting the model to PK data'; corresponds to ~1e4 receptors/cell
    lkon1 <- fixed(log(1e-3)); label("Monovalent CD166 association rate constant (1/(nM*s))")  # Table 1, 'Forward binding rate, k on1' = 1e-3, 'Set equal to a standard value from the literature' (Schlosshauer & Baker 2004)
    lkon2 <- fixed(log(1e-3)); label("Second-arm (bivalent) CD166 association rate constant (1/(nM*s))")  # Table 1, 'Forward binding rate, k on2' = 1e-3; Results: 'k on2 was estimated under the assumption that avidity is not influential', i.e. set equal to kon1
    lkoff1 <- fixed(log(2e-3)); label("CD166 dissociation rate constant per bound arm (1/s)")  # Table 1, 'Reverse binding rate constant, k off1' = 2e-3, derived from kon1 and the apparent affinity Kapp = 1.0 nM. See vignette Errata: the printed relation 'k off1 = k on1 * K app' omits the factor 2 that the tabulated value and the deposited A1 matrix both carry.

    # ---- PROBODY mask and substrate ----
    lkopen <- fixed(log(1.16e-2)); label("Mask opening ('breathing') rate constant per closed arm (1/s)")  # not in Table 1: taken from the Supporting Information 'Model code' A1 matrix (entry (3,2) = 2.32e-2 = 2 * kopen). Corroborated by Ippolito 2024 Table S2, which cites Stroh 2019 for 'Probody unmasking rate' k_o = 0.0116 1/second.
    lkmask <- log(220); label("Fold-masking Kmask = kclose / kopen for mask M1 (unitless)")  # Table 1, 'Fold-masking, K mask' = 220 (M1) and 57 (M2), 'From fit of PK using parametric scan'
    lkcleave <- fixed(log(3e-7)); label("Substrate cleavage rate constant per uncleaved arm (1/s)")  # Table 1, 'Rate constant for cleavage reaction, k cleave' < 3e-7. Results: no lower bound could be distinguished from zero and 3e-7 'was carried forward as the likely maximum value'.
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
    vp <- exp(lvp)
    kel <- exp(lkel)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    kint <- exp(lkint)
    ksynr_peripheral <- exp(lksynr_peripheral)
    kon1 <- exp(lkon1)
    kon2 <- exp(lkon2)
    koff1 <- exp(lkoff1)
    kopen <- exp(lkopen)
    kmask <- exp(lkmask)
    kclose <- kmask * kopen
    kcleave <- exp(lkcleave)

    # CD166 sits at its synthesis / endocytosis steady state before the first dose.
    target_central(0) <- ksynr_central / kint
    target_peripheral(0) <- ksynr_peripheral / kint

    # ---- plasma compartment (deposited-code states 1-17) ----
    d/dt(target_central) <-
      ksynr_central - kint * target_central +
      koff1 * pb_c_oR_central + koff1 * pb_o_oR_central + koff1 * pb_c_mR_central + koff1 * pb_o_mR_central + koff1 * pb_oR_m_central + koff1 * pb_m_mR_central +
      2 * koff1 * pb_oR_oR_central + 2 * koff1 * pb_oR_mR_central + 2 * koff1 * pb_mR_mR_central -
      kon1 * target_central * (pb_c_o_central + 2 * pb_o_o_central + pb_c_m_central + 2 * pb_o_m_central + 2 * pb_m_m_central) -
      kon2 * target_central * (pb_o_oR_central + pb_o_mR_central + pb_oR_m_central + pb_m_mR_central)
    d/dt(pb_c_c_central) <-
      kclose * pb_c_o_central + k21 * pb_c_c_peripheral -
      (2 * kopen + 2 * kcleave + kel + k12) * pb_c_c_central
    d/dt(pb_c_o_central) <-
      2 * kopen * pb_c_c_central + 2 * kclose * pb_o_o_central + koff1 * pb_c_oR_central + k21 * pb_c_o_peripheral -
      (kclose + kopen + 2 * kcleave + kon1 * target_central + kel + k12) * pb_c_o_central
    d/dt(pb_o_o_central) <-
      kopen * pb_c_o_central + koff1 * pb_o_oR_central + k21 * pb_o_o_peripheral -
      (2 * kclose + 2 * kcleave + 2 * kon1 * target_central + kel + k12) * pb_o_o_central
    d/dt(pb_c_m_central) <-
      2 * kcleave * pb_c_c_central + kcleave * pb_c_o_central + kclose * pb_o_m_central + koff1 * pb_c_mR_central + k21 * pb_c_m_peripheral -
      (kopen + kcleave + kon1 * target_central + kel + k12) * pb_c_m_central
    d/dt(pb_o_m_central) <-
      kcleave * pb_c_o_central + 2 * kcleave * pb_o_o_central + kopen * pb_c_m_central + koff1 * (pb_o_mR_central + pb_oR_m_central) + k21 * pb_o_m_peripheral -
      (kclose + kcleave + 2 * kon1 * target_central + kel + k12) * pb_o_m_central
    d/dt(pb_m_m_central) <-
      kcleave * (pb_c_m_central + pb_o_m_central) + koff1 * pb_m_mR_central + k21 * pb_m_m_peripheral -
      (2 * kon1 * target_central + kel + k12) * pb_m_m_central
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

    # ---- observations (Figure 3 plots 'Total Drug (nM)' against time in days) ----
    Cc <- (pb_c_c_central + pb_c_o_central + pb_o_o_central + pb_c_m_central + pb_o_m_central + pb_m_m_central + pb_c_oR_central + pb_o_oR_central + pb_oR_oR_central + pb_c_mR_central + pb_o_mR_central + pb_oR_m_central + pb_oR_mR_central + pb_m_mR_central + pb_mR_mR_central) / vc
    pbIntactPlasma <- (pb_c_c_central + pb_c_o_central + pb_o_o_central + pb_c_oR_central + pb_o_oR_central + pb_oR_oR_central) / vc
    pbCleavedPlasma <- Cc - pbIntactPlasma
    pbPeripheral <- (pb_c_c_peripheral + pb_c_o_peripheral + pb_o_o_peripheral + pb_c_m_peripheral + pb_o_m_peripheral + pb_m_m_peripheral + pb_c_oR_peripheral + pb_o_oR_peripheral + pb_oR_oR_peripheral + pb_c_mR_peripheral + pb_o_mR_peripheral + pb_oR_m_peripheral + pb_oR_mR_peripheral + pb_m_mR_peripheral + pb_mR_mR_peripheral) / vp
    uptakePeripheral <- pb_tmdd_m_peripheral + pb_tmdd_o_peripheral
  })
}
