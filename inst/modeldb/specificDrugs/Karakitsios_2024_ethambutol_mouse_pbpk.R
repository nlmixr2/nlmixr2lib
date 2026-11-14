Karakitsios_2024_ethambutol_mouse_pbpk <- function() {
  description <- paste(
    "PBPK (middle-out, permeability-limited lung; Monolix Suite 2021R1).",
    "Preclinical (mouse). Ethambutol distribution into healthy/uninvolved lung in mice after",
    "oral dosing. The published Chen 2016 two-compartment blood model with first-order",
    "absorption and an absorption lag drives a lung split into pulmonary blood (in",
    "instantaneous equilibrium with extracellular water), extracellular water and",
    "intracellular water, with permeability-limited, pH-dependent transfer of the unbound",
    "non-ionised drug across the lung-cell membrane. The lung states hold CONCENTRATIONS,",
    "not amounts.",
    sep = " "
  )
  reference <- paste(
    "Karakitsios E, Dokoumetzidis A.",
    "Extrapolation of lung pharmacokinetics of antitubercular drugs from preclinical species",
    "to humans using PBPK modelling. J Antimicrob Chemother. 2024;79(6):1362.",
    "doi:10.1093/jac/dkae109.",
    "The blood layer is the mouse model of Chen C, Ortega F, Alameda L, et al. Population",
    "pharmacokinetics, optimized design and sample size determination for rifampicin,",
    "isoniazid, ethambutol and pyrazinamide in the mouse. Eur J Pharm Sci. 2016;93:319-333.",
    "doi:10.1016/j.ejps.2016.07.017.",
    sep = " "
  )
  vignette <- "Karakitsios_2024_antitubercular_lung_pbpk"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")
  dosing <- c("depot", "central")

  compartmentData <- list(
    depot = list(analyte = "ethambutol", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ethambutol", units = "mg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "ethambutol", units = "mg", specimen = "whole blood", verified = TRUE),
    lung_ew = list(analyte = "ethambutol", units = "ug/mL", specimen = "tissue", verified = TRUE),
    lung_iw = list(analyte = "ethambutol", units = "ug/mL", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "mouse",
    n_subjects = NA_integer_,
    n_studies = 2L,
    sex_female_pct = NA_real_,
    weight_median = "0.025 kg (reference mouse of Table S2)",
    disease_state = "Uninvolved (healthy) lung (Muliaditan 2022 lung data; Table 1).",
    dose_range = "100 mg/kg oral, single dose: the human-equivalent dose in mice (Table 1; Figure S6).",
    notes = paste(
      "Mean lung concentrations were fitted. The Chen 2016 model describes BLOOD",
      "concentrations, so Cc here is a blood concentration and enters the lung directly",
      "(Supplementary Section 1). Chen 2016 estimated no between-animal variability for",
      "ethambutol because of its sparse one-sample-per-mouse design.",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # 1. Blood layer: Chen 2016 Table 5 ethambutol (EMB), two compartments
    #    with first-order absorption and a lag, reported per kilogram and
    #    scaled to the reference mouse. Inherited and held constant here.
    # =====================================================================
    bw_ref  <- fixed(0.025)            ; label("Reference body weight of a typical mouse (kg)")   # Table S2 mouse body weight = 0.025 kg
    lka     <- fixed(log(0.869))       ; label("Absorption rate constant ka (1/h)")               # Chen 2016 Table 5 EMB ka = 0.869 1/h
    lcl     <- fixed(log(2560 / 1000)) ; label("Blood clearance CL (L/h/kg)")                     # Chen 2016 Table 5 EMB CL = 2560 mL/h/kg
    lvc     <- fixed(log(2180 / 1000)) ; label("Central volume V (L/kg)")                         # Chen 2016 Table 5 EMB V = 2180 mL/kg
    lq      <- fixed(log(1760 / 1000)) ; label("Inter-compartmental clearance Q (L/h/kg)")        # Chen 2016 Table 5 EMB Q = 1760 mL/h/kg
    lvp     <- fixed(log(4910 / 1000)) ; label("Peripheral volume V2 (L/kg)")                     # Chen 2016 Table 5 EMB V2 = 4910 mL/kg
    lfdepot <- fixed(log(0.64))        ; label("Oral bioavailability F (fraction)")               # Chen 2016 Table 5 EMB F = 64.0%
    lalag   <- fixed(log(0.0577))      ; label("Absorption lag time Tlag (h)")                    # Chen 2016 Table 5 EMB Tlag = 0.0577 h (3.5 min, Results 3.3)

    # =====================================================================
    # 2. Drug-related parameters, independent of species (Table S5).
    #    Ethambutol is a diprotic base.
    # =====================================================================
    pka1 <- fixed(9.55) ; label("First basic pKa (unitless)")    # Table S5 ETH pKa1 = 9.55 (base)
    pka2 <- fixed(6.5)  ; label("Second basic pKa (unitless)")   # Table S5 ETH pKa2 = 6.5 (base)

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
    # 4. Drug- and species-related parameters (Table S6).
    # =====================================================================
    bp        <- fixed(1.08)  ; label("Blood-to-plasma concentration ratio B/P (unitless)")  # Table S6 ETH blood to plasma ratio = 1.08
    fu_plasma <- fixed(0.826) ; label("Fraction unbound in plasma fupls (unitless)")         # Table S6 ETH fu plasma = 0.826
    fu_ew     <- fixed(1)     ; label("Fraction unbound in lung extracellular water fuEW")   # Table S6 ETH fu EW = 1; Table 3 footnote d

    # =====================================================================
    # 5. Lung parameters estimated in mice (Table 3).
    # =====================================================================
    lperm_cells <- log(0.034 * 10) ; label("Permeability of lung cells PermCells (cm/h, log scale)")       # Table S6 ETH 0.034 dm/h = 0.34 cm/h (Table 3, RSE 13.6%)
    lfu_iw      <- log(0.056)      ; label("Fraction unbound in lung intracellular water fuIW (log scale)")  # Table 3 ETH mice fuIW = 0.056 (RSE 12.4%)
  })

  model({
    # 1. Blood layer (Chen 2016), per kg times the reference body weight.
    ka <- exp(lka)
    cl <- exp(lcl) * bw_ref
    vc <- exp(lvc) * bw_ref
    q  <- exp(lq) * bw_ref
    vp <- exp(lvp) * bw_ref

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl / vc * central -
      q / vc * central + q / vp * peripheral1
    d/dt(peripheral1) <- q / vc * central - q / vp * peripheral1
    f(depot) <- exp(lfdepot)
    alag(depot) <- exp(lalag)

    # 2. Ionisation (diprotic base; product form reproduces Table S6 fni).
    x_plasma <- (1 + 10^(pka1 - ph_plasma)) * (1 + 10^(pka2 - ph_plasma))
    x_ew     <- (1 + 10^(pka1 - ph_ew))     * (1 + 10^(pka2 - ph_ew))
    x_iw     <- (1 + 10^(pka1 - ph_iw))     * (1 + 10^(pka2 - ph_iw))
    fni_plasma <- 1 / x_plasma
    fni_ew     <- 1 / x_ew
    fni_iw     <- 1 / x_iw

    perm_cells <- exp(lperm_cells)
    fu_iw      <- exp(lfu_iw)

    # 3. Lung constants (instantaneous PB-EW equilibrium; PS in mL/h).
    kewb     <- fu_plasma * fni_plasma / (fu_ew * fni_ew * bp)
    ps_cells <- perm_cells * sa_cells * 1e4
    v_ew_eff <- v_ew + v_pb / kewb

    # Chen 2016 modelled blood, so the central concentration IS Cblood.
    cblood <- central / vc
    jcells <- ps_cells * (lung_ew * fni_ew - fu_iw * lung_iw * fni_iw)

    # 4. Healthy/uninvolved lung, equations (1) and (2). lung_ew holds the
    #    UNBOUND EW concentration; lung_iw holds the TOTAL IW concentration.
    d/dt(lung_ew) <- fu_ew / v_ew_eff *
      (qlung * cblood - qlung * lung_ew / (fu_ew * kewb) - jcells)
    d/dt(lung_iw) <- jcells / v_iw

    # 5. Observations. Cc is the BLOOD concentration of the Chen 2016 model.
    Cc <- central / vc

    c_lung_ew <- lung_ew / fu_ew
    c_lung_pb <- c_lung_ew / kewb
    c_lung <- (c_lung_pb * v_pb + c_lung_ew * v_ew + lung_iw * v_iw) /
      (v_pb + v_ew + v_iw)
    c_lung_iw <- lung_iw
  })
}
