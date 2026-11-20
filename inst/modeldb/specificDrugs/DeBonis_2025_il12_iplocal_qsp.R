DeBonis_2025_il12_iplocal_qsp <- function() {
  description <- paste(
    "QSP. Local-production model of continuous intraperitoneal (i.p.)",
    "interleukin-12 (IL-12) delivery, as from an implanted cell or gene",
    "therapy (DeBonis 2025, Eqs. S31-S39). Structurally the i.p. form of",
    "the reduced-bioavailability model: IL-12 produced at a zero-order rate",
    "kpr(t) in the peritoneal cavity (depot_ip) drains by lymph flow through",
    "the lymphatic system (lnode) into blood (central), binding the T-cell",
    "IL-12 receptor in the lymphatic system, where STAT4 activation by the",
    "complex up-regulates receptor production. Predictive parameter set",
    "only (Table S11): the i.p. transport rates come from the Lenzi 2002",
    "fit and the receptor parameters from the Motzer 1998 fit. Deliver",
    "kpr(t) as infusion records into depot_ip. Deterministic, no IIV, no",
    "residual error."
  )
  reference <- paste(
    "DeBonis J, Veiseh O, Igoshin OA. Uncovering the interleukin-12",
    "pharmacokinetic desensitization mechanism and its consequences with",
    "mathematical modeling. CPT Pharmacometrics Syst Pharmacol.",
    "2025;14(2):217-228. doi:10.1002/psp4.13258",
    sep = " "
  )
  vignette <- "DeBonis_2025_il12_qsp"
  units <- list(time = "day", dosing = "pmol", concentration = "pg/mL")
  dosing <- "depot_ip"

  compartmentData <- list(
    depot_ip = list(analyte = "IL-12", units = "pmol", specimen = "administration site", verified = TRUE),
    lnode = list(analyte = "IL-12 (free)", units = "pmol", specimen = "lymph", verified = TRUE),
    target = list(
      analyte = "IL-12 receptor (free, on lymphatic T cells)",
      units = "pmol",
      specimen = "lymph",
      verified = TRUE
    ),
    complex = list(analyte = "IL-12:IL-12 receptor complex", units = "pmol", specimen = "lymph", verified = TRUE),
    central = list(analyte = "IL-12", units = "pmol", specimen = "serum", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_studies = 2L,
    disease_state = paste(
      "Simulation scenario for local, continuous i.p. IL-12 production; no",
      "clinical data were fitted with this parameter set"
    ),
    dose_range = paste(
      "Zero-order i.p. production profiles (constant, step-wise and",
      "linearly increasing) reaching the same final rate (Figure 5a)"
    ),
    notes = paste(
      "Every value in Table S11 is carried from another fit: J_LP and J_BL",
      "from the Lenzi 2002 i.p. fit (Table S10, DeBonis_2025_il12_ip_qsp)",
      "and the remaining parameters from the Motzer 1998 subcutaneous fit",
      "(Table S7, DeBonis_2025_il12_qsp)."
    )
  )

  ini({
    lka_ip <- fixed(log(0.94))
    label("First-order IL-12 transport rate from the peritoneal cavity to the lymphatic system JLP (1/day)") # Table S11 'J_LP' (from Table S10)
    lk_lnode_central <- fixed(log(43))
    label("First-order IL-12 transport rate from the lymphatic system to blood JBL (1/day)") # Table S11 'J_BL' (from Table S10)
    lkel <- fixed(log(1.6))
    label("First-order systemic IL-12 clearance rate constant from blood kCB (1/day)") # Table S11 'k_CB' (from Table S7, Motzer)
    lkon <- fixed(log(200))
    label("IL-12:IL-12 receptor association rate constant kon (mL/pmol/day)") # Table S11 'k_on' (from Table S7, Motzer)
    lkint <- fixed(log(0.56))
    label("IL-12:IL-12 receptor complex internalisation rate constant kint (1/day)") # Table S11 'k_int' (from Table S7, Motzer)
    lkdeg <- fixed(log(0.46))
    label("Unbound IL-12 receptor internalisation rate constant kCR (1/day)") # Table S11 'k_CR' (from Table S7, Motzer)
    lntcell_lnode <- fixed(log(3.8e10))
    label("Number of IL-12R+ T cells in the lymphatic system NL (cells)") # Table S11 'N_L' (from Table S7, Motzer)
    lr0_cell <- fixed(log(5.6e-9))
    label("Basal IL-12 receptor per T cell R0 (pmol/cell)") # Table S11 'R_0' (from Table S7, Motzer)
    lkupreg <- fixed(log(5.0e-6))
    label("STAT4-driven receptor up-regulation production coefficient alpha (pmol/cell/day)") # Table S11 'alpha' (from Table S7, Motzer)
    lic_pstat <- fixed(log(1.4e6))
    label("Complexes per T cell for half-maximal STAT4 activation ICpSTAT (complexes/cell)") # Table S11 'IC_pstat' (from Table S7, Motzer)
    lfdepot <- fixed(log(1))
    label("Absolute bioavailability F from the peritoneal cavity (fraction)") # Table S3 'F': 1 for the RB model; Eq. S39
    lv_ip <- fixed(log(1020))
    label("Peritoneal cavity volume VP (mL)") # Not printed in the paper; carried from the i.p.-injection model, authors' deposited code build_RB_Lenzi_model.m ('parameters.Vs = 1020')
    lv_lnode <- fixed(log(33))
    label("Volume of the lymphatic system IL-12 transits VL (mL)") # Table S3 'V_L' = 33 (Varkhede 2018)
    lvc <- fixed(log(5100))
    label("Blood volume VB (mL)") # Table S3 'V_B' = 5100, 'Average human blood volume'
  })

  model({
    ka_ip <- exp(lka_ip)
    k_lnode_central <- exp(lk_lnode_central)
    kel <- exp(lkel)
    kon <- exp(lkon)
    kint <- exp(lkint)
    kdeg <- exp(lkdeg)
    ntcell_lnode <- exp(lntcell_lnode)
    r0_cell <- exp(lr0_cell)
    kupreg <- exp(lkupreg)
    ic_pstat <- exp(lic_pstat)
    fdepot <- exp(lfdepot)
    v_ip <- exp(lv_ip)
    v_lnode <- exp(lv_lnode)
    vc <- exp(lvc)

    # Peritoneal degradation that sets the bioavailability F (Eq. S39, whose
    # J_LS is the i.p. transport rate J_LP); zero when F = 1.
    kdg <- (1 - fdepot) / fdepot * ka_ip
    # Basal receptor production per cell (Eq. S38).
    krec_cell <- kdeg * r0_cell

    # Complexes per T cell (Eq. S36; 6.023e11 molecules per pmol) and the
    # quasi-steady-state STAT4 activation fraction (Eq. S37).
    complex_cell <- 6.023e11 * complex / ntcell_lnode
    pstat <- complex_cell / (complex_cell + ic_pstat)

    # Binding flux in the lymphatic system, pmol/day.
    bind <- kon * lnode * target / v_lnode

    # Eqs. S31-S35 written for amounts (pmol). The zero-order production
    # kpr(t) of Eq. S31 (pmol/day) enters depot_ip as infusion records.
    d/dt(depot_ip) <- -ka_ip * depot_ip - kdg * depot_ip
    d/dt(lnode) <- ka_ip * depot_ip - bind - k_lnode_central * lnode
    d/dt(target) <- ntcell_lnode * krec_cell + ntcell_lnode * kupreg * pstat - bind - kdeg * target
    d/dt(complex) <- bind - kint * complex
    d/dt(central) <- k_lnode_central * lnode - kel * central

    target(0) <- r0_cell * ntcell_lnode

    # Serum and peritoneal (local) IL-12 in pg/mL (70000 pg/pmol).
    Cc <- central / vc * 70000
    Cip <- depot_ip / v_ip * 70000
  })
}
