Karakitsios_2024_pyrazinamide_human_fromMouse_pbpk <- function() {
  description <- paste(
    "PBPK (middle-out, permeability-limited lung; Monolix Suite 2021R1).",
    "Human extrapolation, from BALB/c mice, of pyrazinamide distribution into healthy/uninvolved",
    "lung and into the cellular lesions of a typical adult with pulmonary tuberculosis. The",
    "published Strydom 2019 one-compartment plasma model with first-order absorption drives a",
    "lung split into pulmonary blood (in instantaneous equilibrium with extracellular water),",
    "extracellular water and intracellular water, with permeability-limited transfer of the",
    "unbound non-ionised drug across the lung-cell membrane; extracellular water also feeds a",
    "cellular lesion. Cell permeability and the lesion rate are carried over unchanged from",
    "the mouse fit; the human extracellular unbound fraction is predicted from plasma protein",
    "binding and the lung albumin level. The lung and lesion states hold CONCENTRATIONS, not",
    "amounts. Typical-value model (no between-subject variability).",
    sep = " "
  )
  reference <- paste(
    "Karakitsios E, Dokoumetzidis A.",
    "Extrapolation of lung pharmacokinetics of antitubercular drugs from preclinical species",
    "to humans using PBPK modelling. J Antimicrob Chemother. 2024;79(6):1362.",
    "doi:10.1093/jac/dkae109.",
    "The plasma layer is the pyrazinamide model of Strydom N, Gupta SV, Fox WS, et al.",
    "Tuberculosis drugs' distribution and emergence of resistance in patient's lung lesions:",
    "a mechanistic model and tool for regimen and dose optimization. PLoS Med.",
    "2019;16:e1002773. doi:10.1371/journal.pmed.1002773.",
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
    species = "human",
    n_subjects = 15L,
    n_studies = 1L,
    age_range = "23-58 years (Strydom 2019 Table 1)",
    sex_female_pct = 33.3,
    weight_median = "80.7 kg (reference human of Table S4)",
    disease_state = paste(
      "Adults with pulmonary tuberculosis who had failed treatment and underwent lung",
      "resection surgery (Strydom 2019, South Korea).",
      sep = " "
    ),
    dose_range = "1500 mg oral, single dose (Table 1, Table 4); multiple daily 1500 mg doses in Figure 4.",
    regions = "South Korea (Strydom 2019 clinical study NCT00816426).",
    notes = paste(
      "The extrapolated profiles represent a typical mean patient (Methods). Because BALB/c",
      "mice develop only cellular lesions, the cellular lesion was compared with the Strydom",
      "2019 empirical small-cellular-nodule model (Supplementary Section 5; Table 4).",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # 1. Plasma layer: Strydom 2019 Table 2 pyrazinamide, one compartment
    #    with first-order absorption. Apparent values, held constant.
    # =====================================================================
    lka <- fixed(log(0.554)) ; label("Absorption rate constant ka (1/h)")   # Strydom 2019 Table 2 PZA rate of absorption = 0.554 1/h
    lcl <- fixed(log(2.05))  ; label("Apparent clearance CL/F (L/h)")       # Strydom 2019 Table 2 PZA clearance = 2.05 L/h
    lvc <- fixed(log(30.0))  ; label("Apparent central volume V/F (L)")     # Strydom 2019 Table 2 PZA central volume = 30.0 L

    # =====================================================================
    # 2. Drug-related parameters, independent of species (Table S5).
    # =====================================================================
    pka <- fixed(3.53) ; label("Basic pKa (unitless)")   # Table S5 PZA pKa = 3.53 (base)

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
    ph_lesion <- fixed(5.84)       ; label("Cellular-lesion pH (unitless)")                 # Supplementary Section 1: mean pH 5.84 in human cellular lesions

    # =====================================================================
    # 4. Drug- and species-related parameters (Tables S8 and S12).
    # =====================================================================
    bp        <- fixed(1.12)  ; label("Blood-to-plasma concentration ratio B/P (unitless)")  # Table S8 PZA blood to plasma ratio = 1.12
    fu_plasma <- fixed(0.719) ; label("Fraction unbound in plasma fupls (unitless)")         # Table S8 PZA fu plasma = 0.719
    fu_ew     <- fixed(0.81)  ; label("Fraction unbound in lung extracellular water fuEW")   # Table S8 / Table S12 PZA fu EW = 0.81 (equations 8-9)
    fu_iw     <- fixed(0.98)  ; label("Fraction unbound in lung intracellular water fuIW")   # Table S8 / Table S12 PZA fu IW = 0.98 (acidic-phospholipid affinity 0)
    fu_lesion <- fixed(0.991) ; label("Fraction unbound in the cellular lesion fuLesion")    # Table S12 PZA fu cellular lesion = 0.991 (equation 14)

    # =====================================================================
    # 5. Distribution parameters estimated in mice and held constant
    #    across species (Table 3, Table S8). The lesion rate enters as 1/h.
    # =====================================================================
    lperm_cells  <- fixed(log(0.000021 * 10)) ; label("Permeability of lung cells PermCells (cm/h, log scale)")       # Table S8 PZA 0.000021 dm/h = 2.1E-4 cm/h
    lrate_lesion <- fixed(log(72.5))          ; label("Rate of distribution into the cellular lesion (1/h, log scale)") # Table 3 PZA mice = 72.5, kept constant across species
  })

  model({
    # 1. Plasma layer (Strydom 2019): one compartment, first-order absorption.
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)

    d/dt(depot)   <- -ka * depot
    d/dt(central) <- ka * depot - cl / vc * central

    # 2. Ionisation (monoprotic base).
    x_plasma <- 1 + 10^(pka - ph_plasma)
    x_ew     <- 1 + 10^(pka - ph_ew)
    x_iw     <- 1 + 10^(pka - ph_iw)
    x_lesion <- 1 + 10^(pka - ph_lesion)
    fni_plasma <- 1 / x_plasma
    fni_ew     <- 1 / x_ew
    fni_iw     <- 1 / x_iw
    fni_lesion <- 1 / x_lesion

    perm_cells  <- exp(lperm_cells)
    rate_lesion <- exp(lrate_lesion)

    # Human lesion product fni x fu: Table S12 fu at the lesion pH.
    fnifu_lesion <- fni_lesion * fu_lesion

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

    # 6. Observations: total (homogenate) and unbound concentrations.
    Cc <- central / vc

    c_lung_ew <- lung_ew / fu_ew
    c_lung_pb <- c_lung_ew / kewb
    c_lung <- (c_lung_pb * v_pb + c_lung_ew * v_ew + lung_iw * v_iw) /
      (v_pb + v_ew + v_iw)
    c_lesion <- lesion

    cu_plasma <- Cc * fu_plasma
    cu_ew     <- lung_ew
    cu_lesion <- lesion * fu_lesion
  })
}
