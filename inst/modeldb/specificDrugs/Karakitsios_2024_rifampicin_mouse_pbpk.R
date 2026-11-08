Karakitsios_2024_rifampicin_mouse_pbpk <- function() {
  description <- paste(
    "PBPK (middle-out, permeability-limited lung; Monolix Suite 2021R1).",
    "Preclinical (mouse). Rifampicin distribution into healthy/uninvolved lung in mice after",
    "oral dosing. The published Chen 2016 one-compartment blood model with first-order",
    "absorption drives a lung split into pulmonary blood (in instantaneous equilibrium with",
    "extracellular water), extracellular water and intracellular water, with",
    "permeability-limited, pH-dependent transfer of the unbound non-ionised drug across the",
    "lung-cell membrane. The lung states hold CONCENTRATIONS, not amounts.",
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

  compartmentData <- list(
    depot = list(analyte = "rifampicin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "rifampicin", units = "mg", specimen = "whole blood", verified = TRUE),
    lung_ew = list(analyte = "rifampicin", units = "ug/mL", specimen = "tissue", verified = TRUE),
    lung_iw = list(analyte = "rifampicin", units = "ug/mL", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "mouse",
    n_subjects = NA_integer_,
    n_studies = 2L,
    sex_female_pct = NA_real_,
    weight_median = "0.025 kg (reference mouse of Table S2)",
    disease_state = "Uninvolved (healthy) lung (Muliaditan 2022 lung data; Table 1).",
    dose_range = "10 mg/kg oral, single dose: the human-equivalent dose in mice (Table 1; Figure S2).",
    notes = paste(
      "Mean lung concentrations were fitted. The Chen 2016 model describes BLOOD",
      "concentrations, so Cc here is a blood concentration and enters the lung directly",
      "(Supplementary Section 1). Chen 2016 estimated a larger volume at the lowest dose",
      "(1 mg/kg) and a higher clearance on day 3 (auto-induction); the single-dose 10 mg/kg",
      "values are used.",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # 1. Blood layer: Chen 2016 Table 5 rifampicin, one compartment with
    #    first-order absorption, reported per kilogram and scaled to the
    #    reference mouse. Inherited and held constant here.
    # =====================================================================
    bw_ref  <- fixed(0.025)          ; label("Reference body weight of a typical mouse (kg)")   # Table S2 mouse body weight = 0.025 kg
    lka     <- fixed(log(1.02))      ; label("Absorption rate constant ka (1/h)")               # Chen 2016 Table 5 RIF ka = 1.02 1/h
    lcl     <- fixed(log(79.3 / 1000)) ; label("Blood clearance CL (L/h/kg)")                   # Chen 2016 Table 5 RIF CL = 79.3 mL/h/kg (days 1 and 2)
    lvc     <- fixed(log(1250 / 1000)) ; label("Volume of distribution V (L/kg)")               # Chen 2016 Table 5 RIF V = 1250 mL/kg (doses above the lowest)
    lfdepot <- fixed(log(0.656))     ; label("Oral bioavailability F (fraction)")               # Chen 2016 Table 5 RIF F = 65.6%

    # =====================================================================
    # 2. Drug-related parameters, independent of species (Table S5).
    # =====================================================================
    pka_acid <- fixed(1.7) ; label("Acidic pKa (unitless)")   # Table S5 RIF pKa1 = 1.7
    pka_base <- fixed(7.9) ; label("Basic pKa (unitless)")    # Table S5 RIF pKa2 = 7.9

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
    bp        <- fixed(0.58)  ; label("Blood-to-plasma concentration ratio B/P (unitless)")  # Table S6 RIF blood to plasma ratio = 0.58
    fu_plasma <- fixed(0.029) ; label("Fraction unbound in plasma fupls (unitless)")         # Table S6 RIF fu plasma = 0.029
    fu_ew     <- fixed(1)     ; label("Fraction unbound in lung extracellular water fuEW")   # Table S6 RIF fu EW = 1; Table 3 footnote d

    # =====================================================================
    # 5. Lung parameters estimated in mice (Table 3).
    # =====================================================================
    lperm_cells <- log(144942.33 * 10) ; label("Permeability of lung cells PermCells (cm/h, log scale)")       # Table S6 RIF 144942.33 dm/h = 1.45E+6 cm/h (Table 3, RSE 61.0%)
    lfu_iw      <- log(0.0084)         ; label("Fraction unbound in lung intracellular water fuIW (log scale)")  # Table 3 RIF mice fuIW = 8.4E-3 (RSE 9.7%)
  })

  model({
    # 1. Blood layer (Chen 2016), per kg times the reference body weight.
    ka <- exp(lka)
    cl <- exp(lcl) * bw_ref
    vc <- exp(lvc) * bw_ref

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl / vc * central
    f(depot) <- exp(lfdepot)

    # 2. Ionisation (ampholyte; product form reproduces Table S6 fni).
    x_plasma <- (1 + 10^(ph_plasma - pka_acid)) * (1 + 10^(pka_base - ph_plasma))
    x_ew     <- (1 + 10^(ph_ew - pka_acid))     * (1 + 10^(pka_base - ph_ew))
    x_iw     <- (1 + 10^(ph_iw - pka_acid))     * (1 + 10^(pka_base - ph_iw))
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
  })
}
