DeBonis_2025_il12_qsp <- function() {
  description <- paste(
    "QSP. Reduced-bioavailability model of subcutaneous recombinant human",
    "interleukin-12 (IL-12) pharmacokinetic desensitization in adults with",
    "advanced renal cell carcinoma (DeBonis 2025, the paper's preferred",
    "mechanism). IL-12 moves by lymph flow from the subcutaneous depot",
    "through the lymphatic system (lnode) into blood (central). In the",
    "lymphatic system it binds the IL-12 receptor on T cells irreversibly;",
    "the complex is internalised and, through quasi-steady-state STAT4",
    "activation, up-regulates receptor production, so repeated doses are",
    "sequestered during lymphatic transit and less IL-12 reaches the blood.",
    "Deterministic mechanism model with no IIV and no residual error;",
    "parameters are the global-optimisation fit to the Motzer 1998 phase I",
    "data (Table S7), the set the paper uses for all later analyses."
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

  compartmentData <- list(
    depot = list(analyte = "IL-12", units = "pmol", specimen = "administration site", verified = TRUE),
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
    n_subjects = 51L,
    n_studies = 1L,
    disease_state = "Advanced renal cell carcinoma (phase I rhIL-12 trial)",
    dose_range = paste(
      "0.1-1.5 ug/kg rhIL-12 SC on days 1, 8 and 15 of 28-day cycles: a fixed",
      "cycle at the same dose and an escalating cycle (0.1, 0.5, then the",
      "target dose) (Motzer 1998; Table S1)"
    ),
    weight_assumed = "70 kg for every subject (ug/kg doses converted to pmol)",
    notes = paste(
      "Fitted to mean serum IL-12 time courses digitised from Motzer 1998",
      "(51 patients enrolled; per-cohort counts in Table S1: fixed cycle",
      "15 at 0.5 and 6 at 1.0 ug/kg, escalating cycle 3, 3, 6 and 12 at",
      "target doses 0.5, 1.0, 1.25 and 1.5 ug/kg), not to individual",
      "data. The same structure was refitted separately to",
      "Rakhit 1999 (n = 24) and Portielje 1999 (n = 25 analysed) data",
      "(Table S7); those sets are shown in the vignette."
    )
  )

  ini({
    lka <- log(1.6)
    label("First-order IL-12 transport rate from the SC space to the lymphatic system JLS (1/day)") # Table S7 'J_LS', Motzer column
    lk_lnode_central <- log(70)
    label("First-order IL-12 transport rate from the lymphatic system to blood JBL (1/day)") # Table S7 'J_BL', Motzer column
    lkel <- log(1.6)
    label("First-order systemic IL-12 clearance rate constant from blood kCB (1/day)") # Table S7 'k_CB', Motzer column
    lkon <- log(200)
    label("IL-12:IL-12 receptor association rate constant kon (mL/pmol/day)") # Table S7 'k_on', Motzer column
    lkint <- log(0.56)
    label("IL-12:IL-12 receptor complex internalisation rate constant kint (1/day)") # Table S7 'k_int', Motzer column
    lkdeg <- log(0.46)
    label("Unbound IL-12 receptor internalisation rate constant kCR (1/day)") # Table S7 'k_CR', Motzer column
    lntcell_lnode <- log(3.8e10)
    label("Number of IL-12R+ T cells in the lymphatic system NL (cells)") # Table S7 'N_L', Motzer column
    lr0_cell <- log(5.6e-9)
    label("Basal IL-12 receptor per T cell R0 (pmol/cell)") # Table S7 'R_0', Motzer column
    lkupreg <- log(5.0e-6)
    label("STAT4-driven receptor up-regulation production coefficient alpha (pmol/cell/day)") # Table S7 'alpha', Motzer column
    lic_pstat <- log(1.4e6)
    label("Complexes per T cell for half-maximal STAT4 activation ICpSTAT (complexes/cell)") # Table S7 'IC_pstat', Motzer column
    lfdepot <- fixed(log(1))
    label("Absolute SC bioavailability F (fraction)") # Table S3 'F': 1 for the RB model, 'Assumed to be one to isolate bioavailability mechanism'
    lv_lnode <- fixed(log(33))
    label("Volume of the lymphatic system IL-12 transits VL (mL)") # Table S3 'V_L' = 33 (Varkhede 2018)
    lvc <- fixed(log(5100))
    label("Blood volume VB (mL)") # Table S3 'V_B' = 5100, 'Average human blood volume'
  })

  model({
    ka <- exp(lka)
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
    v_lnode <- exp(lv_lnode)
    vc <- exp(lvc)

    # SC degradation that sets the bioavailability F (Eq. 4 / Eq. S19); zero
    # when F = 1, as in this model.
    kdg <- (1 - fdepot) / fdepot * ka
    # Basal receptor production per cell, constrained so the receptor pool
    # sits at R0 * NL before any IL-12 (Eq. 11 / Eq. S18).
    krec_cell <- kdeg * r0_cell

    # Complexes per T cell (Eq. S16; 6.023e11 molecules per pmol) and the
    # quasi-steady-state STAT4 activation fraction (Eq. S17).
    complex_cell <- 6.023e11 * complex / ntcell_lnode
    pstat <- complex_cell / (complex_cell + ic_pstat)

    # Binding flux in the lymphatic system, pmol/day: kon * IL12_L * R_L *
    # V_L with IL12_L = lnode / V_L and R_L = target / V_L.
    bind <- kon * lnode * target / v_lnode

    # Eqs. S11-S15 written for amounts (pmol) rather than concentrations;
    # the paper's concentration form is each state divided by its volume.
    d/dt(depot) <- -ka * depot - kdg * depot
    d/dt(lnode) <- ka * depot - bind - k_lnode_central * lnode
    d/dt(target) <- ntcell_lnode * krec_cell + ntcell_lnode * kupreg * pstat - bind - kdeg * target
    d/dt(complex) <- bind - kint * complex
    d/dt(central) <- k_lnode_central * lnode - kel * central

    # Eq. 10 / Eq. S10 baseline: R_SS,L = R0 * NL / VL, i.e. R0 * NL pmol.
    target(0) <- r0_cell * ntcell_lnode

    # Serum IL-12 in pg/mL: pmol/mL times the 70 kDa molecular weight
    # (70000 pg/pmol).
    Cc <- central / vc * 70000
  })
}
