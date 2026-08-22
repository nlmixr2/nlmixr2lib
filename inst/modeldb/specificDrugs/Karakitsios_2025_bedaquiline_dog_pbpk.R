Karakitsios_2025_bedaquiline_dog_pbpk <- function() {
  description <- paste(
    "PBPK (middle-out, permeability-limited lung; Monolix Suite 2023R1).",
    "Preclinical (beagle dog). Bedaquiline distribution into healthy/uninvolved lung tissue in",
    "beagle dogs. An empirical three-compartment oral plasma model drives a permeability-limited",
    "lung model split into pulmonary blood, extracellular water, intracellular water and",
    "lysosomes, with pH-dependent ion trapping in the lysosomes. The lung states hold",
    "CONCENTRATIONS, not amounts. The dog plasma model was fitted to dog data; the lung",
    "parameters were estimated in mice and carried over unchanged, with the lysosome:cytosol",
    "partition rescaled to the lower canine intracellular pH (6.69 in humans, 7.04 in dogs, 7.36",
    "in rodents), which lowers the lysosomal trapping gradient relative to rodents. Dogs were a",
    "validation species for the healthy-lung extrapolation -- no infected-lung data exist in",
    "dogs, so this model has no lesion or caseum compartments.",
    sep = " "
  )
  reference <- paste(
    "Karakitsios E, Della Pasqua O, Dokoumetzidis A.",
    "Extrapolation of lung pharmacokinetics of bedaquiline across species using",
    "physiologically-based pharmacokinetic modelling.",
    "Br J Clin Pharmacol. 2025;91(11):3167-3178.",
    "doi:10.1002/bcp.70163.",
    sep = " "
  )
  vignette <- "Karakitsios_2025_bedaquiline_lung_pbpk"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    depot       = list(analyte = "bedaquiline", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "bedaquiline", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "bedaquiline", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "bedaquiline", units = "mg", specimen = "plasma", verified = TRUE),
    lung_ew     = list(analyte = "bedaquiline", units = "ug/mL", specimen = "tissue", verified = TRUE),
    lung_iw     = list(analyte = "bedaquiline", units = "ug/mL", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Not a model covariate. Table S7 reports the dog plasma model per kilogram, and this file",
        "multiplies those values by the fixed 10 kg Simcyp reference dog of Table S3 footnote a,",
        "because the Table S3 lung blood flow, sub-volumes and cell surface area are fixed",
        "per-species quantities that are not weight-scalable within a species. Dose in mg for a",
        "10 kg animal (e.g. 10 mg/kg = 100 mg).",
        sep = " "
      )
    )
  )

  population <- list(
    species        = "beagle dog",
    n_subjects     = NA_integer_,
    n_studies      = 1L,
    sex_female_pct = NA_real_,
    weight_median  = "10 kg (Simcyp reference dog; Table S3 footnote a)",
    disease_state  = "Healthy (uninfected) dogs; uninvolved lung tissue only.",
    dose_range     = paste(
      "Oral. 2.5 mg/kg/day and 10 mg/kg/day multiple doses (male and female beagles).",
      "Janssen R&D data, Table S1. Only single-timepoint lung tissue homogenate concentrations",
      "were available in dogs.",
      sep = " "
    ),
    notes          = paste(
      "Mean concentrations were modelled rather than individual data because of the small number",
      "of animals per sampling time (Methods 2.1). Dogs served as a validation species for the",
      "mouse-to-larger-species extrapolation of healthy-lung disposition; Figure 3 shows the",
      "observed-versus-predicted comparison after multiple doses, with all predictions within",
      "2-fold of observed and no systematic bias. High doses and intravenous data were excluded",
      "to avoid phospholipidosis, which in dogs made tissue concentrations rise",
      "disproportionately to plasma exposure (Discussion).",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # 1. Empirical plasma model (Table S7), reported per kilogram and scaled
    #    to the fixed reference dog weight below.
    # =====================================================================
    bw_ref  <- fixed(10)    ; label("Reference body weight of a typical beagle dog (kg)")        # Table S3 footnote a: Simcyp dog weight = 10 kg
    lka     <- log(0.48)    ; label("Absorption rate constant ka (1/h)")                         # Table S7 ka = 0.48 (RSE 21.5%)
    lcl     <- log(0.085)   ; label("Clearance CL (L/h/kg)")                                     # Table S7 Cl = 0.085 L/h/kg (RSE 9.29%)
    lvc     <- log(1.86)    ; label("Central volume V1 (L/kg)")                                  # Table S7 V1 = 1.86 L/kg (RSE 24.2%)
    lq      <- log(0.27)    ; label("Inter-compartmental clearance Q2 (L/h/kg)")                 # Table S7 Q2 = 0.27 L/h/kg (RSE 7.33%)
    lvp     <- log(52.99)   ; label("First peripheral volume V2 (L/kg)")                         # Table S7 V2 = 52.99 L/kg (RSE 21.4%)
    lq2     <- log(0.57)    ; label("Inter-compartmental clearance Q3 (L/h/kg)")                 # Table S7 Q3 = 0.57 L/h/kg (RSE 32.1%)
    lvp2    <- log(1.87)    ; label("Second peripheral volume V3 (L/kg)")                        # Table S7 V3 = 1.87 L/kg (RSE 16.5%)
    lfdepot <- log(0.53)    ; label("Oral bioavailability F (unitless)")                         # Table S7 F = 0.53 (RSE 22.2%)

    # =====================================================================
    # 2. Drug-related parameters, independent of species (Table S2).
    # =====================================================================
    pka1 <- fixed(8.5) ; label("First (strong) basic pKa (unitless)")   # Table S2 pKa1 = 8.5 (base), Janssen R&D data
    pka2 <- fixed(1.1) ; label("Second (weak) basic pKa (unitless)")    # Table S2 pKa2 = 1.1 (base), ACD predicted

    # =====================================================================
    # 3. Species-related lung physiology for the dog (Table S3).
    # =====================================================================
    qlung       <- fixed(122940) ; label("Lung blood flow rate Qlung (mL/h)")                    # Table S3 dog blood flow rate = 122940 mL/h
    sa_cells    <- fixed(506000) ; label("Surface area of lung cells SAcells (cm^2)")            # Table S3 dog surface area lung cells = 506000 cm^2
    v_pb        <- fixed(12.54)  ; label("Volume of pulmonary blood VPB (mL)")                   # Table S3 dog volume of pulmonary blood = 12.54 mL
    v_ew        <- fixed(14.98)  ; label("Volume of lung extracellular water VEW (mL)")          # Table S3 dog volume of EW = 14.98 mL
    v_lysosome  <- fixed(0.8)    ; label("Volume of lung lysosomes VLYS (mL)")                   # Table S3 dog volume of lysosomes = 0.8 mL
    v_iw        <- fixed(51.7)   ; label("Volume of lung intracellular water VIW (mL)")          # Table S3 dog volume of IW = 51.7 mL
    ph_plasma   <- fixed(7.4)    ; label("Plasma pH (unitless)")                                 # Table S3 dog pH_plasma = 7.4
    ph_ew       <- fixed(7.4)    ; label("Lung extracellular-water pH (unitless)")               # Table S3 dog pH_EW = 7.4
    ph_iw       <- fixed(7.04)   ; label("Lung intracellular-water pH (unitless)")               # Table S3 dog pH_IW = 7.04

    # =====================================================================
    # 4. Drug- AND species-related parameters (Table S4).
    # =====================================================================
    bp        <- fixed(1.0)    ; label("Blood-to-plasma concentration ratio B/P (unitless)")     # Table S4 dog blood to plasma ratio = 1.0
    fu_plasma <- fixed(6.9e-4) ; label("Fraction unbound in plasma fupls (unitless)")            # Table S4 dog fupls = 6.9E-4
    fu_ew     <- fixed(1.1e-3) ; label("Fraction unbound in lung extracellular water fuEW")      # Table S4 dog fuEW = 1.1E-3
    fu_iw     <- fixed(1.0)    ; label("Fraction unbound in lung intracellular water fuIW")      # Supporting Information: "a fixed value of 1 was used for fuIW"

    # =====================================================================
    # 5. Lung parameters carried over from the mouse fit (Table S8) and
    #    rescaled to the canine pH_IW by equation (S9) in model().
    # =====================================================================
    lperm_cells   <- fixed(log(55.7))       ; label("Permeability of lung cells PermCells (cm/h, log scale)")  # Table S8 Permcells = 5.57E+1 cm/h, estimated in mice and held constant across species
    lfulys_fnilys <- fixed(log(1.9411e-8))  ; label("Product fuLYS*fniLYS in lysosomes (unitless, log scale)") # derived from Table S8 KpLYS:IW/fuIW = 3.48E+6 in mice via equations (S8)-(S9)

    # =====================================================================
    # 6. Random effects (Table S7).
    # =====================================================================
    etalfdepot ~ 1.7689  # Table S7 omega_F = 1.33 (Monolix SD of the log-scale random effect); variance = 1.33^2
    etalvc     ~ 0.0961  # Table S7 omega_V1 = 0.31; variance = 0.31^2
    etalvp     ~ 0.3721  # Table S7 omega_V2 = 0.61; variance = 0.61^2 (the paper's V2 is this file's vp)
    etalq2     ~ 0.3136  # Table S7 omega_Q3 = 0.56; variance = 0.56^2 (the paper's Q3 is this file's q2)

    # =====================================================================
    # 7. Residual error (Table S7) -- combined additive and proportional. The
    #    lung layer was extrapolated rather than fitted in dogs, so no dog
    #    lung residual error is reported.
    # =====================================================================
    addSd  <- 0.74 ; label("Plasma additive residual error (ug/mL)")           # Table S7 a = 0.74 (RSE 35.7%)
    propSd <- 0.14 ; label("Plasma proportional residual error (fraction)")    # Table S7 b = 0.14 (RSE 6.76%)
  })

  model({
    # 1. Individual plasma parameters, scaled to the reference dog weight.
    ka  <- exp(lka)
    cl  <- exp(lcl)          * bw_ref
    vc  <- exp(lvc + etalvc) * bw_ref
    q   <- exp(lq)           * bw_ref
    vp  <- exp(lvp + etalvp) * bw_ref
    q2  <- exp(lq2 + etalq2) * bw_ref
    vp2 <- exp(lvp2)         * bw_ref

    perm_cells <- exp(lperm_cells)

    # 2. Ionisation (diprotic base product form; see vignette Errata).
    x_plasma <- (1 + 10^(pka1 - ph_plasma)) * (1 + 10^(pka2 - ph_plasma))
    x_ew     <- (1 + 10^(pka1 - ph_ew))     * (1 + 10^(pka2 - ph_ew))
    x_iw     <- (1 + 10^(pka1 - ph_iw))     * (1 + 10^(pka2 - ph_iw))

    fi_plasma <- 1 / x_plasma
    fi_ew     <- 1 / x_ew
    fi_iw     <- 1 / x_iw

    kp_lysosome <- fi_iw / exp(lfulys_fnilys)   # equation (S9)

    # 3. Derived lung constants.
    kewb     <- (fu_plasma / fu_ew) * (x_ew / x_plasma)
    ps       <- perm_cells * sa_cells
    v_ew_eff <- (v_ew + v_pb / kewb) / fu_ew
    v_iw_eff <- v_iw + kp_lysosome * v_lysosome

    cblood <- (central / vc) * bp
    jperm  <- ps * (lung_ew * fi_ew - fu_iw * lung_iw * fi_iw)

    # 4. Plasma ODEs.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - cl * central / vc -
                          q * central / vc + q * peripheral1 / vp -
                          q2 * central / vc + q2 * peripheral2 / vp2
    d/dt(peripheral1) <-  q * central / vc - q * peripheral1 / vp
    d/dt(peripheral2) <-  q2 * central / vc - q2 * peripheral2 / vp2

    f(depot) <- exp(lfdepot + etalfdepot)

    # 5. Healthy/uninvolved lung -- equations (S1) and (S4).
    d/dt(lung_ew) <- (qlung * cblood - (qlung / fu_ew) * (lung_ew / kewb) - jperm) / v_ew_eff
    d/dt(lung_iw) <- jperm / v_iw_eff

    # 6. Observations.
    Cc <- central / vc

    c_lung_ew  <- lung_ew / fu_ew
    c_lung_pb  <- c_lung_ew / kewb
    c_lysosome <- kp_lysosome * lung_iw
    c_lung <- (c_lung_pb * v_pb + c_lung_ew * v_ew +
               lung_iw * v_iw + c_lysosome * v_lysosome) /
              (v_pb + v_ew + v_iw + v_lysosome)
    cu_plasma <- Cc * fu_plasma

    Cc ~ add(addSd) + prop(propSd)
  })
}
