Karakitsios_2024_isoniazid_human_fromMouse_pbpk <- function() {
  description <- paste(
    "PBPK (middle-out, permeability-limited lung; Monolix Suite 2021R1).",
    "Human extrapolation, from mice, of isoniazid distribution into healthy/uninvolved lung in",
    "a typical adult with pulmonary tuberculosis. The published Strydom 2019 two-compartment",
    "plasma model with lagged first-order absorption drives a lung split into pulmonary blood",
    "(in instantaneous equilibrium with extracellular water), extracellular water and",
    "intracellular water, with permeability-limited transfer of the unbound non-ionised drug",
    "across the lung-cell membrane. Cell permeability is carried over unchanged from the mouse",
    "fit; the human extracellular unbound fraction is predicted from plasma protein binding",
    "and the lung albumin level. The lung states hold CONCENTRATIONS, not amounts.",
    "Typical-value model (no between-subject variability).",
    sep = " "
  )
  reference <- paste(
    "Karakitsios E, Dokoumetzidis A.",
    "Extrapolation of lung pharmacokinetics of antitubercular drugs from preclinical species",
    "to humans using PBPK modelling. J Antimicrob Chemother. 2024;79(6):1362.",
    "doi:10.1093/jac/dkae109.",
    "The plasma layer is the isoniazid model of Strydom N, Gupta SV, Fox WS, et al.",
    "Tuberculosis drugs' distribution and emergence of resistance in patient's lung lesions:",
    "a mechanistic model and tool for regimen and dose optimization. PLoS Med.",
    "2019;16:e1002773. doi:10.1371/journal.pmed.1002773.",
    sep = " "
  )
  vignette <- "Karakitsios_2024_antitubercular_lung_pbpk"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")
  dosing <- c("depot", "central")

  compartmentData <- list(
    depot = list(analyte = "isoniazid", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "isoniazid", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "isoniazid", units = "mg", specimen = "plasma", verified = TRUE),
    lung_ew = list(analyte = "isoniazid", units = "ug/mL", specimen = "tissue", verified = TRUE),
    lung_iw = list(analyte = "isoniazid", units = "ug/mL", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 15L,
    n_studies = 1L,
    age_range = "23-58 years (Strydom 2019 Table 1)",
    sex_female_pct = 33.3,
    weight_median = "80.7 kg (reference human of Table S4)",
    disease_state = paste(
      "Adults with pulmonary tuberculosis who had failed treatment and underwent lung",
      "resection surgery (Strydom 2019, South Korea). NAT2 phenotypes: 11 fast, 2",
      "intermediate and 2 slow acetylators (Strydom 2019 Table 1; Supplementary Section 5).",
      sep = " "
    ),
    dose_range = "300 mg oral, single dose (Table 1, Table 4).",
    regions = "South Korea (Strydom 2019 clinical study NCT00816426).",
    notes = paste(
      "The extrapolated profiles represent a typical mean patient (Methods). Strydom 2019",
      "estimated separate clearances for slow, intermediate and fast acetylators; the typical",
      "mean patient uses their cohort-weighted mean (see the clearance line in ini()).",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # 1. Plasma layer: Strydom 2019 Table 2 isoniazid, two compartments
    #    with lagged first-order absorption. Apparent values, held constant.
    # =====================================================================
    lka   <- fixed(log(0.738)) ; label("Absorption rate constant ka (1/h)")             # Strydom 2019 Table 2 INH rate of absorption = 0.738 1/h
    lalag <- fixed(log(0.130)) ; label("Absorption lag time (h)")                       # Strydom 2019 Table 2 INH lag time = 0.130 h
    # Cohort-weighted acetylator clearance: 11 fast (38.1), 2 intermediate
    # (24.0) and 2 slow (9.41 L/h) of 15 patients = 32.4 L/h. This value
    # reproduces the Table 4 'Observed' isoniazid lung AUC, Cmax and Tmax.
    lcl   <- fixed(log((11 * 38.1 + 2 * 24.0 + 2 * 9.41) / 15)) ; label("Apparent clearance CL/F (L/h)")   # Strydom 2019 Table 2 INH CL s 9.41, i 24.0, f 38.1 L/h; Table 1 phenotypes
    lvc   <- fixed(log(48.9))  ; label("Apparent central volume V/F (L)")              # Strydom 2019 Table 2 INH central volume = 48.9 L
    lq    <- fixed(log(6.43))  ; label("Apparent inter-compartmental clearance Q/F (L/h)")  # Strydom 2019 Table 2 INH intercomp clearance = 6.43 L/h
    lvp   <- fixed(log(40.4))  ; label("Apparent peripheral volume Vp/F (L)")           # Strydom 2019 Table 2 INH peripheral volume = 40.4 L

    # =====================================================================
    # 2. Drug-related parameters, independent of species (Table S5).
    # =====================================================================
    pka <- fixed(3.27) ; label("Basic pKa (unitless)")   # Table S5 INH pKa = 3.27 (base)

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
    bp        <- fixed(0.89)  ; label("Blood-to-plasma concentration ratio B/P (unitless)")  # Table S8 INH blood to plasma ratio = 0.89
    fu_plasma <- fixed(0.516) ; label("Fraction unbound in plasma fupls (unitless)")         # Table S8 INH fu plasma = 0.516
    fu_ew     <- fixed(0.64)  ; label("Fraction unbound in lung extracellular water fuEW")   # Table S8 / Table S12 INH fu EW = 0.64 (equations 8-9)
    fu_iw     <- fixed(0.984) ; label("Fraction unbound in lung intracellular water fuIW")   # Table S8 / Table S12 INH fu IW = 0.984 (acidic-phospholipid affinity 0)

    # =====================================================================
    # 5. Cell permeability estimated in mice and held constant across
    #    species (Table 3, Table S8).
    # =====================================================================
    lperm_cells <- fixed(log(0.000037 * 10)) ; label("Permeability of lung cells PermCells (cm/h, log scale)")   # Table S8 INH 0.000037 dm/h = 3.7E-4 cm/h
  })

  model({
    # 1. Plasma layer (Strydom 2019): two compartments, lagged first-order
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

    # 2. Ionisation (monoprotic base).
    x_plasma <- 1 + 10^(pka - ph_plasma)
    x_ew     <- 1 + 10^(pka - ph_ew)
    x_iw     <- 1 + 10^(pka - ph_iw)
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

    # 4. Healthy/uninvolved lung, equations (1) and (2). lung_ew holds the
    #    UNBOUND EW concentration; lung_iw holds the TOTAL IW concentration.
    d/dt(lung_ew) <- fu_ew / v_ew_eff *
      (qlung * cblood - qlung * lung_ew / (fu_ew * kewb) - jcells)
    d/dt(lung_iw) <- jcells / v_iw

    # 5. Observations: total (homogenate) and unbound concentrations.
    Cc <- central / vc

    c_lung_ew <- lung_ew / fu_ew
    c_lung_pb <- c_lung_ew / kewb
    c_lung <- (c_lung_pb * v_pb + c_lung_ew * v_ew + lung_iw * v_iw) /
      (v_pb + v_ew + v_iw)

    cu_plasma <- Cc * fu_plasma
    cu_ew     <- lung_ew
  })
}
