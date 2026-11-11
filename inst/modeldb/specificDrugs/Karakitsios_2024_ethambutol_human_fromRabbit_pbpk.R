Karakitsios_2024_ethambutol_human_fromRabbit_pbpk <- function() {
  description <- paste(
    "PBPK (middle-out, permeability-limited lung; Monolix Suite 2021R1).",
    "Human extrapolation, from rabbits, of ethambutol distribution into healthy lung,",
    "including the intracellular (alveolar-cell) concentration, in a typical healthy adult.",
    "A two-compartment serum model with lagged first-order absorption, fitted by the authors",
    "to the Peloquin 1999 mean data (Table S11), drives a lung split into pulmonary blood (in",
    "instantaneous equilibrium with extracellular water), extracellular water and",
    "intracellular water, with permeability-limited, pH-dependent transfer of the unbound",
    "non-ionised drug across the lung-cell membrane. Cell permeability is carried over",
    "unchanged from the rabbit fit and the intracellular unbound fraction is rescaled through",
    "a species-constant acidic-phospholipid affinity. The lung states hold CONCENTRATIONS, not",
    "amounts. Typical-value model (no between-subject variability).",
    sep = " "
  )
  reference <- paste(
    "Karakitsios E, Dokoumetzidis A.",
    "Extrapolation of lung pharmacokinetics of antitubercular drugs from preclinical species",
    "to humans using PBPK modelling. J Antimicrob Chemother. 2024;79(6):1362.",
    "doi:10.1093/jac/dkae109.",
    "Serum data from Peloquin CA, Bulpitt AE, Jaresko GS, et al. Pharmacokinetics of",
    "ethambutol under fasting conditions, with food, and with antacids. Antimicrob Agents",
    "Chemother. 1999;43:568-572. doi:10.1128/AAC.43.3.568.",
    sep = " "
  )
  vignette <- "Karakitsios_2024_antitubercular_lung_pbpk"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    depot = list(analyte = "ethambutol", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ethambutol", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "ethambutol", units = "mg", specimen = "serum", verified = TRUE),
    lung_ew = list(analyte = "ethambutol", units = "ug/mL", specimen = "tissue", verified = TRUE),
    lung_iw = list(analyte = "ethambutol", units = "ug/mL", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 14L,
    n_studies = 1L,
    sex_female_pct = NA_real_,
    weight_median = "80.7 kg (reference human of Table S4)",
    disease_state = paste(
      "Healthy fasted adults for the serum model (Peloquin 1999, n = 14; Table S11); the lung",
      "prediction was compared with alveolar-cell concentrations of healthy individuals after",
      "15 mg/kg (Supplementary Section 5; Figure S12).",
      sep = " "
    ),
    dose_range = "25 mg/kg oral for the serum fit (Table S11); 15 mg/kg oral for the lung extrapolation (Table 1).",
    notes = paste(
      "The serum model was fitted to digitised MEAN concentrations, so it carries no",
      "between-subject variability (Table S11). Ethambutol is absent from Table 4: its",
      "extrapolation is shown graphically only (Figure S12).",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # 1. Serum layer: Table S11, the authors' human two-compartment model.
    # =====================================================================
    lalag <- log(0.38)   ; label("Absorption lag time Tlag (h)")                   # Table S11 Tlag = 0.38 h (RSE 8.25%)
    lka   <- log(0.51)   ; label("Absorption rate constant ka (1/h)")              # Table S11 Ka = 0.51 1/h (RSE 23.7%)
    lcl   <- log(59.52)  ; label("Apparent clearance CL/F (L/h)")                  # Table S11 Cl = 59.52 L/h (RSE 24.5%)
    lvc   <- log(214.62) ; label("Apparent central volume V1/F (L)")               # Table S11 V1 = 214.62 L (RSE 23.6%)
    lq    <- log(44.22)  ; label("Apparent inter-compartmental clearance Q/F (L/h)")  # Table S11 Q = 44.22 L/h (RSE 31.2%)
    lvp   <- log(739.1)  ; label("Apparent peripheral volume V2/F (L)")            # Table S11 V2 = 739.1 L (RSE 83.7%)

    # =====================================================================
    # 2. Drug-related parameters, independent of species (Table S5).
    # =====================================================================
    pka1 <- fixed(9.55) ; label("First basic pKa (unitless)")    # Table S5 ETH pKa1 = 9.55 (base)
    pka2 <- fixed(6.5)  ; label("Second basic pKa (unitless)")   # Table S5 ETH pKa2 = 6.5 (base)

    # =====================================================================
    # 3. Species-related lung physiology for humans (Table S4).
    # =====================================================================
    qlung     <- fixed(356 * 1000) ; label("Lung blood flow rate Qlung (mL/h)")             # Table S4 human blood flow rate = 356 L/hr
    v_pb      <- fixed(101.195)    ; label("Volume of pulmonary blood VPB (mL)")            # Table S4 human volume of pulmonary blood = 101.195 mL
    v_ew      <- fixed(89.161)     ; label("Volume of lung extracellular water VEW (mL)")   # Table S4 human extracellular volume of lung tissue = 89.161 mL
    v_iw      <- fixed(356.644)    ; label("Volume of lung intracellular water VIW (mL)")   # Table S4 human intracellular volume of lung tissue = 356.644 mL
    sa_cells  <- fixed(96)         ; label("Surface area of lung cells SAcells (m^2)")      # Table S4 human surface area of lung cells = 96 m^2
    ph_plasma <- fixed(7.4)        ; label("Plasma pH (unitless)")                          # Table S4 human pH plasma = 7.4
    ph_ew     <- fixed(7.4)        ; label("Lung extracellular-water pH (unitless)")        # Table S4 human pH EW-Lung = 7.4
    ph_iw     <- fixed(6.69)       ; label("Lung intracellular-water pH (unitless)")        # Table S4 human pH IW-Lung = 6.69

    # =====================================================================
    # 4. Drug- and species-related parameters (Tables S8 and S12).
    # =====================================================================
    bp        <- fixed(1.04)  ; label("Blood-to-plasma concentration ratio B/P (unitless)")  # Table S8 ETH blood to plasma ratio = 1.04
    fu_plasma <- fixed(0.853) ; label("Fraction unbound in plasma fupls (unitless)")         # Table S8 ETH fu plasma = 0.853
    fu_ew     <- fixed(1)     ; label("Fraction unbound in lung extracellular water fuEW")   # Table S8 / Table S12 ETH fu EW = 1 (assumed)
    fu_iw     <- fixed(0.25)  ; label("Fraction unbound in lung intracellular water fuIW")   # Table S8 / Table S12 ETH fu IW = 0.25 (from rabbits, equations 5-6)

    # =====================================================================
    # 5. Cell permeability estimated in rabbits and held constant across
    #    species (Table 3, Table S8).
    # =====================================================================
    lperm_cells <- fixed(log(0.023 * 10)) ; label("Permeability of lung cells PermCells (cm/h, log scale)")   # Table S8 ETH 0.023 dm/h (from rabbits) = 0.23 cm/h

    # =====================================================================
    # 6. Serum residual error (Table S11, constant error model).
    # =====================================================================
    addSd <- 0.16 ; label("Serum additive residual error (ug/mL)")   # Table S11 a = 0.16 mg/L (RSE 14.1%)
  })

  model({
    # 1. Serum layer (Table S11): two compartments, lagged first-order
    #    absorption.
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl / vc * central -
      q / vc * central + q / vp * peripheral1
    d/dt(peripheral1) <- q / vc * central - q / vp * peripheral1
    alag(depot) <- exp(lalag)

    # 2. Ionisation (diprotic base; product form reproduces Table S8 fni).
    x_plasma <- (1 + 10^(pka1 - ph_plasma)) * (1 + 10^(pka2 - ph_plasma))
    x_ew     <- (1 + 10^(pka1 - ph_ew))     * (1 + 10^(pka2 - ph_ew))
    x_iw     <- (1 + 10^(pka1 - ph_iw))     * (1 + 10^(pka2 - ph_iw))
    fni_plasma <- 1 / x_plasma
    fni_ew     <- 1 / x_ew
    fni_iw     <- 1 / x_iw

    perm_cells <- exp(lperm_cells)

    # 3. Lung constants (instantaneous PB-EW equilibrium; PS in mL/h).
    kewb     <- fu_plasma * fni_plasma / (fu_ew * fni_ew * bp)
    ps_cells <- perm_cells * sa_cells * 1e4
    v_ew_eff <- v_ew + v_pb / kewb

    cblood <- central / vc * bp
    jcells <- ps_cells * (lung_ew * fni_ew - fu_iw * lung_iw * fni_iw)

    # 4. Healthy lung, equations (1) and (2). lung_ew holds the UNBOUND EW
    #    concentration; lung_iw holds the TOTAL IW concentration.
    d/dt(lung_ew) <- fu_ew / v_ew_eff *
      (qlung * cblood - qlung * lung_ew / (fu_ew * kewb) - jcells)
    d/dt(lung_iw) <- jcells / v_iw

    # 5. Observations. The intracellular (alveolar-cell) concentration is
    #    the quantity compared with human data (Figure S12).
    Cc <- central / vc

    c_lung_ew <- lung_ew / fu_ew
    c_lung_pb <- c_lung_ew / kewb
    c_lung <- (c_lung_pb * v_pb + c_lung_ew * v_ew + lung_iw * v_iw) /
      (v_pb + v_ew + v_iw)
    c_lung_iw <- lung_iw

    cu_plasma <- Cc * fu_plasma
    cu_ew     <- lung_ew

    Cc ~ add(addSd)
  })
}
