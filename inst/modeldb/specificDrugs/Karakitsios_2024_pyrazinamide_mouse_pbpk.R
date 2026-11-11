Karakitsios_2024_pyrazinamide_mouse_pbpk <- function() {
  description <- paste(
    "PBPK (middle-out, permeability-limited lung; Monolix Suite 2021R1).",
    "Preclinical (mouse, BALB/c). Pyrazinamide distribution into healthy/uninvolved lung and",
    "into the cellular lesions of Mycobacterium tuberculosis-infected BALB/c mice after oral",
    "dosing. The published Irwin 2016 one-compartment plasma model with first-order",
    "absorption drives a lung split into pulmonary blood (in instantaneous equilibrium with",
    "extracellular water), extracellular water and intracellular water, with",
    "permeability-limited transfer of the unbound non-ionised drug across the lung-cell",
    "membrane; extracellular water also feeds a cellular lesion by permeability-limited",
    "transfer. Pyrazinamide is a weak base, so binding is to extracellular albumin and the",
    "fraction unbound in extracellular water is the estimated quantity. The lung and lesion",
    "states hold CONCENTRATIONS, not amounts.",
    sep = " "
  )
  reference <- paste(
    "Karakitsios E, Dokoumetzidis A.",
    "Extrapolation of lung pharmacokinetics of antitubercular drugs from preclinical species",
    "to humans using PBPK modelling. J Antimicrob Chemother. 2024;79(6):1362.",
    "doi:10.1093/jac/dkae109.",
    "The plasma layer is the BALB/c model of Irwin SM, Prideaux B, Lyon ER, et al.",
    "Bedaquiline and pyrazinamide treatment responses are affected by pulmonary lesion",
    "heterogeneity in Mycobacterium tuberculosis infected C3HeB/FeJ mice. ACS Infect Dis.",
    "2016;2:251-267. doi:10.1021/acsinfecdis.5b00127.",
    sep = " "
  )
  vignette <- "Karakitsios_2024_antitubercular_lung_pbpk"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    depot = list(analyte = "pyrazinamide", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "pyrazinamide", units = "mg", specimen = "plasma", verified = TRUE),
    lung_ew = list(analyte = "pyrazinamide", units = "ug/mL", specimen = "tissue", verified = TRUE),
    lung_iw = list(analyte = "pyrazinamide", units = "ug/mL", specimen = "tissue", verified = TRUE),
    lesion = list(analyte = "pyrazinamide", units = "ug/mL", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "mouse (BALB/c)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    sex_female_pct = NA_real_,
    weight_median = "0.025 kg (reference mouse of Table S2)",
    disease_state = paste(
      "Mycobacterium tuberculosis-infected BALB/c mice, which develop cellular (not caseous)",
      "lesions (Irwin 2016; Supplementary Section 5).",
      sep = " "
    ),
    dose_range = "150 mg/kg oral, single dose: the human-equivalent dose in mice (Table 1; Figure S7).",
    notes = paste(
      "Mean lung and lesion concentrations were fitted. The Irwin 2016 model describes PLASMA",
      "concentrations, which are converted to blood with the blood-to-plasma ratio before",
      "entering the lung (Supplementary Section 1).",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # 1. Plasma layer: Irwin 2016 Supplemental Table 1 BALB/c pyrazinamide,
    #    one compartment with first-order absorption (population geometric
    #    means), reported per kilogram and scaled to the reference mouse.
    #    Apparent values; held constant here.
    # =====================================================================
    bw_ref <- fixed(0.025)      ; label("Reference body weight of a typical mouse (kg)")   # Table S2 mouse body weight = 0.025 kg
    lka    <- fixed(log(18.48)) ; label("Absorption rate constant ka (1/h)")               # Irwin 2016 Supplemental Table 1 BALB/c PZA ka GM = 18.48 1/h
    lcl    <- fixed(log(0.46))  ; label("Plasma clearance CL (L/h/kg)")                    # Irwin 2016 Supplemental Table 1 BALB/c PZA CL GM = 0.46 L/h/kg
    lvc    <- fixed(log(0.96))  ; label("Volume of distribution V (L/kg)")                 # Irwin 2016 Supplemental Table 1 BALB/c PZA V GM = 0.96 L/kg

    # =====================================================================
    # 2. Drug-related parameters, independent of species (Table S5).
    # =====================================================================
    pka <- fixed(3.53) ; label("Basic pKa (unitless)")   # Table S5 PZA pKa = 3.53 (base)

    # =====================================================================
    # 3. Species-related lung physiology for the mouse (Table S2).
    # =====================================================================
    qlung     <- fixed(14.01 * 60) ; label("Lung blood flow rate Qlung (mL/h)")             # Table S2 mouse blood flow rate = 14.01 mL/min
    v_pb      <- fixed(0.03135)    ; label("Volume of pulmonary blood VPB (mL)")            # Table S2 mouse volume of pulmonary blood = 0.03135 mL
    v_ew      <- fixed(0.019155)   ; label("Volume of lung extracellular water VEW (mL)")   # Table S2 mouse extracellular volume of lung tissue = 0.019155 mL
    v_iw      <- fixed(0.099495)   ; label("Volume of lung intracellular water VIW (mL)")   # Table S2 mouse intracellular volume of lung tissue = 0.099495 mL
    sa_cells  <- fixed(0.05)       ; label("Surface area of lung cells SAcells (m^2)")      # Table S2 mouse surface area of lung cells = 0.05 m^2
    ph_plasma <- fixed(7.391)      ; label("Plasma pH (unitless)")                          # Table S2 mouse pH plasma = 7.391
    ph_ew     <- fixed(7.4)        ; label("Lung extracellular-water pH (unitless)")        # Table S2 mouse pH EW-Lung = 7.4
    ph_iw     <- fixed(7.36)       ; label("Lung intracellular-water pH (unitless)")        # Table S2 mouse pH IW-Lung = 7.36

    # =====================================================================
    # 4. Drug- and species-related parameters (Table S6). fuIW follows
    #    from an acidic-phospholipid affinity of 0 (Table 3 footnote c).
    # =====================================================================
    bp        <- fixed(0.79)  ; label("Blood-to-plasma concentration ratio B/P (unitless)")  # Table S6 PZA blood to plasma ratio = 0.79
    fu_plasma <- fixed(0.591) ; label("Fraction unbound in plasma fupls (unitless)")         # Table S6 PZA fu plasma = 0.591
    fu_iw     <- fixed(0.95)  ; label("Fraction unbound in lung intracellular water fuIW")   # Table 3 / Table S6 PZA mice fuIW = 0.95 (predicted)

    # =====================================================================
    # 5. Lung and lesion parameters estimated in BALB/c mice (Table 3).
    #    The lesion rate is printed in L/h but enters equation (10) as 1/h.
    # =====================================================================
    lperm_cells   <- log(0.000021 * 10) ; label("Permeability of lung cells PermCells (cm/h, log scale)")       # Table S6 PZA 0.000021 dm/h = 2.1E-4 cm/h (Table 3, RSE 63.9%)
    lfu_ew        <- log(0.25)          ; label("Fraction unbound in lung extracellular water fuEW (log scale)")  # Table 3 PZA mice fuEW = 0.25 (RSE 6.78%)
    lrate_lesion  <- log(72.5)          ; label("Rate of distribution into the cellular lesion (1/h, log scale)") # Table 3 PZA mice = 72.5 (RSE 54.4%)
    lfnifu_lesion <- log(0.53)          ; label("Product fni x fu in the cellular lesion (unitless, log scale)")   # Table 3 PZA mice = 0.53 (RSE 3.04%)
  })

  model({
    # 1. Plasma layer (Irwin 2016), per kg times the reference body weight.
    ka <- exp(lka)
    cl <- exp(lcl) * bw_ref
    vc <- exp(lvc) * bw_ref

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl / vc * central

    # 2. Ionisation (monoprotic base).
    x_plasma <- 1 + 10^(pka - ph_plasma)
    x_ew     <- 1 + 10^(pka - ph_ew)
    x_iw     <- 1 + 10^(pka - ph_iw)
    fni_plasma <- 1 / x_plasma
    fni_ew     <- 1 / x_ew
    fni_iw     <- 1 / x_iw

    perm_cells   <- exp(lperm_cells)
    fu_ew        <- exp(lfu_ew)
    rate_lesion  <- exp(lrate_lesion)
    fnifu_lesion <- exp(lfnifu_lesion)

    # 3. Lung constants (instantaneous PB-EW equilibrium; PS in mL/h).
    kewb     <- fu_plasma * fni_plasma / (fu_ew * fni_ew * bp)
    ps_cells <- perm_cells * sa_cells * 1e4
    v_ew_eff <- v_ew + v_pb / kewb

    cblood <- central / vc * bp
    jcells <- ps_cells * (lung_ew * fni_ew - fu_iw * lung_iw * fni_iw)

    # 4. Healthy/uninvolved lung, equations (1) and (2). lung_ew holds the
    #    UNBOUND EW concentration; lung_iw holds the TOTAL IW concentration.
    d/dt(lung_ew) <- fu_ew / v_ew_eff *
      (qlung * cblood - qlung * lung_ew / (fu_ew * kewb) - jcells)
    d/dt(lung_iw) <- jcells / v_iw

    # 5. Cellular lesion, equation (10).
    d/dt(lesion) <- rate_lesion * (lung_ew * fni_ew - fnifu_lesion * lesion)

    # 6. Observations: volume-weighted uninvolved-lung homogenate (Table S1).
    Cc <- central / vc

    c_lung_ew <- lung_ew / fu_ew
    c_lung_pb <- c_lung_ew / kewb
    c_lung <- (c_lung_pb * v_pb + c_lung_ew * v_ew + lung_iw * v_iw) /
      (v_pb + v_ew + v_iw)
    c_lesion <- lesion
  })
}
