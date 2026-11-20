Karakitsios_2024_rifampicin_human_fromRabbit_pbpk <- function() {
  description <- paste(
    "PBPK (middle-out, permeability-limited lung; Monolix Suite 2021R1).",
    "Human extrapolation, from rabbits, of rifampicin distribution into healthy/uninvolved",
    "lung and into tuberculosis-infected lung tissue (cellular lesion and cavity caseum) in",
    "a typical adult with pulmonary tuberculosis. The published Strydom 2019 one-compartment",
    "plasma model with transit-compartment oral absorption drives a lung split into pulmonary",
    "blood (in instantaneous equilibrium with extracellular water), extracellular water and",
    "intracellular water, with permeability-limited, pH-dependent transfer of the unbound",
    "non-ionised drug across the lung-cell membrane. Extracellular water feeds a cellular",
    "lesion by permeability-limited transfer, and the lesion feeds the cavity caseum by passive",
    "diffusion. Cell permeability, lesion and caseum rates and the caseum partition ratio are",
    "carried over unchanged from the rabbit fit; human physiology, protein binding and",
    "ionisation replace the rabbit values. The lung, lesion and caseum states hold",
    "CONCENTRATIONS, not amounts. Typical-value model (no between-subject variability).",
    sep = " "
  )
  reference <- paste(
    "Karakitsios E, Dokoumetzidis A.",
    "Extrapolation of lung pharmacokinetics of antitubercular drugs from preclinical species",
    "to humans using PBPK modelling. J Antimicrob Chemother. 2024;79(6):1362.",
    "doi:10.1093/jac/dkae109.",
    "The plasma layer is the rifampicin model of Strydom N, Gupta SV, Fox WS, et al.",
    "Tuberculosis drugs' distribution and emergence of resistance in patient's lung lesions:",
    "a mechanistic model and tool for regimen and dose optimization. PLoS Med.",
    "2019;16:e1002773. doi:10.1371/journal.pmed.1002773.",
    sep = " "
  )
  vignette <- "Karakitsios_2024_antitubercular_lung_pbpk"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    depot = list(analyte = "rifampicin", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "rifampicin", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "rifampicin", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "rifampicin", units = "mg", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "rifampicin", units = "mg", specimen = "administration site", verified = TRUE),
    transit5 = list(analyte = "rifampicin", units = "mg", specimen = "administration site", verified = TRUE),
    transit6 = list(analyte = "rifampicin", units = "mg", specimen = "administration site", verified = TRUE),
    transit7 = list(analyte = "rifampicin", units = "mg", specimen = "administration site", verified = TRUE),
    transit8 = list(analyte = "rifampicin", units = "mg", specimen = "administration site", verified = TRUE),
    transit9 = list(analyte = "rifampicin", units = "mg", specimen = "administration site", verified = TRUE),
    transit10 = list(analyte = "rifampicin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "rifampicin", units = "mg", specimen = "plasma", verified = TRUE),
    lung_ew = list(analyte = "rifampicin", units = "ug/mL", specimen = "tissue", verified = TRUE),
    lung_iw = list(analyte = "rifampicin", units = "ug/mL", specimen = "tissue", verified = TRUE),
    lesion = list(analyte = "rifampicin", units = "ug/mL", specimen = "tissue", verified = TRUE),
    caseum = list(analyte = "rifampicin", units = "ug/mL", specimen = "tissue", verified = TRUE)
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
      "resection surgery (Strydom 2019, South Korea), providing plasma and lesion-specific",
      "drug concentrations.",
      sep = " "
    ),
    dose_range = "600 mg oral, single dose (Table 1, Table 4); multiple daily 600 mg doses in Figure 4.",
    regions = "South Korea (Strydom 2019 clinical study NCT00816426).",
    notes = paste(
      "The extrapolated profiles represent a typical mean patient (Methods). Human lung data",
      "were not available as raw data, so the authors compared their predictions with",
      "simulations of the Strydom 2019 empirical lung and lesion models (Table 4 'Observed').",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # 1. Plasma layer: Strydom 2019 Table 2 rifampicin, one compartment
    #    with first-order absorption from a transit chain (Savic approach).
    #    Apparent values, inherited and held constant.
    # =====================================================================
    lka <- fixed(log(1.55)) ; label("Absorption rate constant ka (1/h)")   # Strydom 2019 Table 2 RIF rate of absorption = 1.55 1/h
    lcl <- fixed(log(5.72)) ; label("Apparent clearance CL/F (L/h)")       # Strydom 2019 Table 2 RIF clearance = 5.72 L/h
    lvc <- fixed(log(52.3)) ; label("Apparent central volume V/F (L)")     # Strydom 2019 Table 2 RIF central volume = 52.3 L
    # Strydom 2019 does not report its transit-chain parameters. MTT was
    # back-solved from the Table 4 'Observed' Tmax of all four rifampicin
    # tissues; the number of transits is not identifiable and is set to 10.
    lmtt <- fixed(log(1.4)) ; label("Mean transit time MTT (h)")             # back-solved from Table 4 'Observed' Tmax (not reported by Strydom 2019)
    ntr  <- fixed(10)       ; label("Number of transit compartments (count)") # not identifiable from Table 4; set to 10

    # =====================================================================
    # 2. Drug-related parameters, independent of species (Table S5).
    # =====================================================================
    pka_acid <- fixed(1.7) ; label("Acidic pKa (unitless)")   # Table S5 RIF pKa1 = 1.7
    pka_base <- fixed(7.9) ; label("Basic pKa (unitless)")    # Table S5 RIF pKa2 = 7.9

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
    bp        <- fixed(0.59)  ; label("Blood-to-plasma concentration ratio B/P (unitless)")  # Table S8 RIF blood to plasma ratio = 0.59
    fu_plasma <- fixed(0.078) ; label("Fraction unbound in plasma fupls (unitless)")         # Table S8 RIF fu plasma = 0.078
    fu_ew     <- fixed(1)     ; label("Fraction unbound in lung extracellular water fuEW")   # Table S8 / Table S12 RIF fu EW = 1 (assumed)
    fu_iw     <- fixed(0.058) ; label("Fraction unbound in lung intracellular water fuIW")   # Table S8 / Table S12 RIF fu IW = 0.058 (from rabbits, equations 5-6)
    fu_lesion <- fixed(0.174) ; label("Fraction unbound in the cellular lesion fuLesion")    # Table S12 RIF fu cellular lesion = 0.174 (equation 14)
    fu_caseum <- fixed(0.073) ; label("Fraction unbound in caseum fuCaseum")                 # Sarathy 2016 Supplementary Table 1 RIF caseum fu = 7.3% (Supplementary Section 1)

    # =====================================================================
    # 5. Distribution parameters estimated in rabbits and held constant
    #    across species (Table 3, Table S8). Rates enter as 1/h.
    # =====================================================================
    lperm_cells  <- fixed(log(3216.69 * 10)) ; label("Permeability of lung cells PermCells (cm/h, log scale)")       # Table S8 RIF 3216.69 dm/h (from rabbits) = 3.22E+4 cm/h
    lrate_lesion <- fixed(log(1.58e7))       ; label("Rate of distribution into the cellular lesion (1/h, log scale)") # Table 3 RIF rabbits = 1.58E+7, kept constant across species
    lrate_caseum <- fixed(log(0.22))         ; label("Rate of distribution into the caseum (1/h, log scale)")          # Table 3 RIF rabbits = 0.22, kept constant across species
    lpc_caseum   <- fixed(log(1.5))          ; label("Caseum-to-lesion partition ratio PC (unitless, log scale)")      # Table 3 RIF rabbits = 1.5, kept constant across species
  })

  model({
    # 1. Plasma layer (Strydom 2019): depot -> 10 transit compartments at
    #    ktr = ntr / MTT -> central at ka.
    ka  <- exp(lka)
    cl  <- exp(lcl)
    vc  <- exp(lvc)
    ktr <- ntr / exp(lmtt)

    d/dt(depot)     <- -ktr * depot
    d/dt(transit1)  <- ktr * depot - ktr * transit1
    d/dt(transit2)  <- ktr * transit1 - ktr * transit2
    d/dt(transit3)  <- ktr * transit2 - ktr * transit3
    d/dt(transit4)  <- ktr * transit3 - ktr * transit4
    d/dt(transit5)  <- ktr * transit4 - ktr * transit5
    d/dt(transit6)  <- ktr * transit5 - ktr * transit6
    d/dt(transit7)  <- ktr * transit6 - ktr * transit7
    d/dt(transit8)  <- ktr * transit7 - ktr * transit8
    d/dt(transit9)  <- ktr * transit8 - ktr * transit9
    d/dt(transit10) <- ktr * transit9 - ka * transit10
    d/dt(central)   <- ka * transit10 - cl / vc * central

    # 2. Ionisation (ampholyte; product form reproduces Table S8 fni).
    x_plasma <- (1 + 10^(ph_plasma - pka_acid)) * (1 + 10^(pka_base - ph_plasma))
    x_ew     <- (1 + 10^(ph_ew - pka_acid))     * (1 + 10^(pka_base - ph_ew))
    x_iw     <- (1 + 10^(ph_iw - pka_acid))     * (1 + 10^(pka_base - ph_iw))
    x_lesion <- (1 + 10^(ph_lesion - pka_acid)) * (1 + 10^(pka_base - ph_lesion))
    fni_plasma <- 1 / x_plasma
    fni_ew     <- 1 / x_ew
    fni_iw     <- 1 / x_iw
    fni_lesion <- 1 / x_lesion

    perm_cells  <- exp(lperm_cells)
    rate_lesion <- exp(lrate_lesion)
    rate_caseum <- exp(lrate_caseum)
    pc_caseum   <- exp(lpc_caseum)

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

    # 5. Infected lung, equations (10) and (11).
    d/dt(lesion) <- rate_lesion * (lung_ew * fni_ew - fnifu_lesion * lesion)
    d/dt(caseum) <- rate_caseum * (pc_caseum * lesion - caseum)

    # 6. Observations: total (homogenate) and unbound concentrations.
    Cc <- central / vc

    c_lung_ew <- lung_ew / fu_ew
    c_lung_pb <- c_lung_ew / kewb
    c_lung <- (c_lung_pb * v_pb + c_lung_ew * v_ew + lung_iw * v_iw) /
      (v_pb + v_ew + v_iw)
    c_lesion <- lesion
    c_caseum <- caseum

    cu_plasma <- Cc * fu_plasma
    cu_ew     <- lung_ew
    cu_lesion <- lesion * fu_lesion
    cu_caseum <- caseum * fu_caseum
  })
}
