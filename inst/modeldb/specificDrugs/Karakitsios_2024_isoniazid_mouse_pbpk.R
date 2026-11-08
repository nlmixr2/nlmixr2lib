Karakitsios_2024_isoniazid_mouse_pbpk <- function() {
  description <- paste(
    "PBPK (middle-out, permeability-limited lung; Monolix Suite 2021R1).",
    "Preclinical (mouse). Isoniazid distribution into healthy/uninvolved lung in mice after",
    "oral dosing. The published Chen 2016 one-compartment blood model with first-order",
    "absorption and dose-dependent clearance drives a lung split into pulmonary blood (in",
    "instantaneous equilibrium with extracellular water), extracellular water and",
    "intracellular water, with permeability-limited transfer of the unbound non-ionised drug",
    "across the lung-cell membrane. Isoniazid is a weak base, so binding is to extracellular",
    "albumin and the fraction unbound in extracellular water is the estimated quantity. The",
    "lung states hold CONCENTRATIONS, not amounts.",
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

  covariateData <- list(
    DOSE_ISONIAZID_MGKG = list(
      description = "Administered oral isoniazid dose level per kilogram of body weight",
      units = "mg/kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives the Chen 2016 dose-dependent blood clearance",
        "CL = 855 * (1 - 0.752 * (DOSE - 0.2) / (17.5 + (DOSE - 0.2))) mL/h/kg, typical clearance",
        "at 0.2 mg/kg. The lung model was fitted at 5 mg/kg (Supplementary Section 5), which",
        "gives CL = 716.6 mL/h/kg. Supply the same value as the dose given (mg/kg).",
        sep = " "
      ),
      source_name = "dose"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "isoniazid", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "isoniazid", units = "mg", specimen = "whole blood", verified = TRUE),
    lung_ew = list(analyte = "isoniazid", units = "ug/mL", specimen = "tissue", verified = TRUE),
    lung_iw = list(analyte = "isoniazid", units = "ug/mL", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "mouse",
    n_subjects = NA_integer_,
    n_studies = 2L,
    sex_female_pct = NA_real_,
    weight_median = "0.025 kg (reference mouse of Table S2)",
    disease_state = "Uninvolved (healthy) lung (Muliaditan 2022 lung data; Table 1).",
    dose_range = paste(
      "5 mg/kg oral, single dose (Table 1; Figure S8). Chosen over the 25 mg/kg lung data as",
      "the closer human-equivalent dose for a cohort dominated by fast acetylators",
      "(Supplementary Section 5).",
      sep = " "
    ),
    notes = paste(
      "Mean lung concentrations were fitted. The Chen 2016 model describes BLOOD",
      "concentrations, so Cc here is a blood concentration and enters the lung directly",
      "(Supplementary Section 1). Chen 2016 also estimated a smaller volume at the lowest",
      "dose (0.2 mg/kg, 437 mL/kg); the volume for higher doses is used.",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # 1. Blood layer: Chen 2016 Table 5 isoniazid, one compartment with
    #    first-order absorption and dose-dependent clearance, reported per
    #    kilogram and scaled to the reference mouse. Held constant here.
    # =====================================================================
    bw_ref  <- fixed(0.025)           ; label("Reference body weight of a typical mouse (kg)")         # Table S2 mouse body weight = 0.025 kg
    lka     <- fixed(log(12.6))       ; label("Absorption rate constant ka (1/h)")                     # Chen 2016 Table 5 INH ka = 12.6 1/h
    lcl     <- fixed(log(855 / 1000)) ; label("Blood clearance at a 0.2 mg/kg dose TVCL (L/h/kg)")     # Chen 2016 Table 5 INH CL = 855 mL/h/kg (typical clearance at 0.2 mg/kg, Results 3.2)
    clmax   <- fixed(0.752)           ; label("Maximum fractional decrease in clearance CLMAX (fraction)")  # Chen 2016 Table 5 INH CLMAX = 75.2%
    cl50    <- fixed(17.5)            ; label("Dose giving half of CLMAX CL50 (mg/kg)")                # Chen 2016 Table 5 INH CL50 = 17,500 ug/kg
    lvc     <- fixed(log(989 / 1000)) ; label("Volume of distribution V (L/kg)")                       # Chen 2016 Table 5 INH V = 989 mL/kg (doses above the lowest)
    lfdepot <- fixed(log(0.843))      ; label("Oral bioavailability F (fraction)")                     # Chen 2016 Table 5 INH F = 84.3%

    # =====================================================================
    # 2. Drug-related parameters, independent of species (Table S5).
    # =====================================================================
    pka <- fixed(3.27) ; label("Basic pKa (unitless)")   # Table S5 INH pKa = 3.27 (base)

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
    bp        <- fixed(0.84)  ; label("Blood-to-plasma concentration ratio B/P (unitless)")  # Table S6 INH blood to plasma ratio = 0.84
    fu_plasma <- fixed(0.396) ; label("Fraction unbound in plasma fupls (unitless)")         # Table S6 INH fu plasma = 0.396
    fu_iw     <- fixed(0.97)  ; label("Fraction unbound in lung intracellular water fuIW")   # Table 3 / Table S6 INH mice fuIW = 0.97 (predicted)

    # =====================================================================
    # 5. Lung parameters estimated in mice (Table 3).
    # =====================================================================
    lperm_cells <- log(0.000037 * 10) ; label("Permeability of lung cells PermCells (cm/h, log scale)")       # Table S6 INH 0.000037 dm/h = 3.7E-4 cm/h (Table 3, RSE 19.6%)
    lfu_ew      <- log(0.17)          ; label("Fraction unbound in lung extracellular water fuEW (log scale)")  # Table 3 INH mice fuEW = 0.17 (RSE 3.83%)
  })

  model({
    # 1. Blood layer (Chen 2016): dose-dependent clearance, per kg times
    #    the reference body weight.
    ka <- exp(lka)
    cl <- exp(lcl) * (1 - clmax * (DOSE_ISONIAZID_MGKG - 0.2) /
      (cl50 + (DOSE_ISONIAZID_MGKG - 0.2))) * bw_ref
    vc <- exp(lvc) * bw_ref

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl / vc * central
    f(depot) <- exp(lfdepot)

    # 2. Ionisation (monoprotic base).
    x_plasma <- 1 + 10^(pka - ph_plasma)
    x_ew     <- 1 + 10^(pka - ph_ew)
    x_iw     <- 1 + 10^(pka - ph_iw)
    fni_plasma <- 1 / x_plasma
    fni_ew     <- 1 / x_ew
    fni_iw     <- 1 / x_iw

    perm_cells <- exp(lperm_cells)
    fu_ew      <- exp(lfu_ew)

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
