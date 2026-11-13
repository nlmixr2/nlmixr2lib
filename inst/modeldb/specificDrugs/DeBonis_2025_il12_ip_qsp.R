DeBonis_2025_il12_ip_qsp <- function() {
  description <- paste(
    "QSP. Reduced-bioavailability model of recombinant human interleukin-12",
    "(IL-12) adapted to a single intraperitoneal (i.p.) injection (DeBonis",
    "2025, Eqs. S21-S30), fitted to i.p. and serum IL-12 in patients with",
    "Mullerian, gastrointestinal or mesothelioma peritoneal malignancies",
    "(Lenzi 2002). IL-12 drains by lymph flow from the peritoneal cavity",
    "(depot_ip) through the lymphatic system (lnode) into blood (central);",
    "in the lymphatic system it binds the T-cell IL-12 receptor, and STAT4",
    "activation by the complex up-regulates receptor production. Only a",
    "fraction beta of the nominal dose appears in the peritoneal cavity.",
    "Deterministic mechanism model with no IIV and no residual error",
    "(Table S10)."
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
    n_studies = 1L,
    disease_state = paste(
      "Mullerian carcinoma, gastrointestinal primary malignancy or",
      "mesothelioma confined to the peritoneal cavity (Lenzi 2002 phase I",
      "intraperitoneal rhIL-12 trial)"
    ),
    dose_range = "Single 0.3 ug/kg i.p. rhIL-12 dose (21 ug, 300 pmol, at an assumed 70 kg)",
    weight_assumed = "70 kg (ug/kg dose converted to pmol)",
    notes = paste(
      "Fitted to one mean i.p. and one mean serum IL-12 time course from",
      "Lenzi 2002 (Figure S11). The supplement does not give the number of",
      "patients in that profile. The paper refits only this one dataset",
      "with the i.p. model; the local-production model",
      "(DeBonis_2025_il12_iplocal_qsp) combines its transport parameters",
      "with the Motzer receptor parameters."
    )
  )

  ini({
    lka_ip <- log(0.94)
    label("First-order IL-12 transport rate from the peritoneal cavity to the lymphatic system JLP (1/day)") # Table S10 row printed 'J_LS' (the i.p.-to-lymph rate, J_LP in Eq. S22)
    lk_lnode_central <- log(43)
    label("First-order IL-12 transport rate from the lymphatic system to blood JBL (1/day)") # Table S10 'J_BL'
    lkel <- log(2.0)
    label("First-order systemic IL-12 clearance rate constant from blood kCB (1/day)") # Table S10 'k_CB'
    lkon <- log(5.4)
    label("IL-12:IL-12 receptor association rate constant kon (mL/pmol/day)") # Table S10 'k_on'
    lkint <- log(0.46)
    label("IL-12:IL-12 receptor complex internalisation rate constant kint (1/day)") # Table S10 'k_int'
    lkdeg <- log(0.46)
    label("Unbound IL-12 receptor internalisation rate constant kCR (1/day)") # Table S10 'k_CR'
    lntcell_lnode <- log(1.0e9)
    label("Number of IL-12R+ T cells in the lymphatic system NL (cells)") # Table S10 'N_L'
    lr0_cell <- log(1.7e-9)
    label("Basal IL-12 receptor per T cell R0 (pmol/cell)") # Table S10 'R_0'
    lkupreg <- log(1.0e-5)
    label("STAT4-driven receptor up-regulation production coefficient alpha (pmol/cell/day)") # Table S10 'alpha'
    lic_pstat <- log(9.1e2)
    label("Complexes per T cell for half-maximal STAT4 activation ICpSTAT (complexes/cell)") # Table S10 'IC_pstat'
    lfdose_ip <- log(5.3e-3)
    label("Fraction of the nominal i.p. dose present in the peritoneal cavity at injection beta (fraction)") # Table S10 'beta'; Eq. S30
    lfdepot <- fixed(log(1))
    label("Absolute bioavailability F from the peritoneal cavity (fraction)") # Table S3 'F': 1 for the RB model; Eq. S29
    lv_ip <- fixed(log(1020))
    label("Peritoneal cavity volume including the injection volume VP (mL)") # Not printed in the paper; authors' deposited code build_RB_Lenzi_model.m ('parameters.Vs = 1020', i.p. volume + injection volume)
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
    fdose_ip <- exp(lfdose_ip)
    fdepot <- exp(lfdepot)
    v_ip <- exp(lv_ip)
    v_lnode <- exp(lv_lnode)
    vc <- exp(lvc)

    # Peritoneal degradation that sets the bioavailability F (Eq. S29, whose
    # J_LS is the i.p. transport rate J_LP); zero when F = 1.
    kdg <- (1 - fdepot) / fdepot * ka_ip
    # Basal receptor production per cell (Eq. S28).
    krec_cell <- kdeg * r0_cell

    # Complexes per T cell (Eq. S26; 6.023e11 molecules per pmol) and the
    # quasi-steady-state STAT4 activation fraction (Eq. S27).
    complex_cell <- 6.023e11 * complex / ntcell_lnode
    pstat <- complex_cell / (complex_cell + ic_pstat)

    # Binding flux in the lymphatic system, pmol/day.
    bind <- kon * lnode * target / v_lnode

    # Eqs. S21-S25 written for amounts (pmol).
    d/dt(depot_ip) <- -ka_ip * depot_ip - kdg * depot_ip
    d/dt(lnode) <- ka_ip * depot_ip - bind - k_lnode_central * lnode
    d/dt(target) <- ntcell_lnode * krec_cell + ntcell_lnode * kupreg * pstat - bind - kdeg * target
    d/dt(complex) <- bind - kint * complex
    d/dt(central) <- k_lnode_central * lnode - kel * central

    # Eq. S30: IL12_P,0 = beta * D / V_P.
    f(depot_ip) <- fdose_ip
    target(0) <- r0_cell * ntcell_lnode

    # Serum and peritoneal IL-12 in pg/mL (70000 pg/pmol).
    Cc <- central / vc * 70000
    Cip <- depot_ip / v_ip * 70000
  })
}
