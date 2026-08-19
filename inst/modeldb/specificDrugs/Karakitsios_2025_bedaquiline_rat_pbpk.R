Karakitsios_2025_bedaquiline_rat_pbpk <- function() {
  description <- paste(
    "PBPK (middle-out, permeability-limited lung; Monolix Suite 2023R1).",
    "Preclinical (rat). Bedaquiline distribution into healthy/uninvolved lung tissue in",
    "Sprague-Dawley rats. An empirical three-compartment oral plasma model drives a",
    "permeability-limited lung model split into pulmonary blood, extracellular water,",
    "intracellular water and lysosomes, with pH-dependent ion trapping in the lysosomes. The",
    "lung states hold CONCENTRATIONS, not amounts. The rat plasma model was fitted to rat data;",
    "the lung parameters were estimated in mice and carried over unchanged, with the",
    "lysosome:cytosol partition rescaled to the rat intracellular pH. Rats were the validation",
    "species for the healthy-lung extrapolation -- no infected-lung data exist in rats, so this",
    "model has no lesion or caseum compartments.",
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
        "Not a model covariate. Table S6 reports the rat plasma model per kilogram, and this file",
        "multiplies those values by the fixed 0.25 kg Simcyp reference rat of Table S3 footnote a,",
        "because the Table S3 lung blood flow, sub-volumes and cell surface area are fixed",
        "per-species quantities that are not weight-scalable within a species. Dose in mg for a",
        "0.25 kg animal (e.g. 20 mg/kg = 5 mg).",
        sep = " "
      )
    )
  )

  population <- list(
    species        = "rat (Sprague-Dawley)",
    n_subjects     = NA_integer_,
    n_studies      = 1L,
    sex_female_pct = NA_real_,
    weight_median  = "0.25 kg (Simcyp reference rat; Table S3 footnote a)",
    disease_state  = "Healthy (uninfected) rats; uninvolved lung tissue only.",
    dose_range     = paste(
      "Oral. 20 mg/kg single dose (male Sprague-Dawley); 6 mg/kg/day and 20 mg/kg/day multiple",
      "doses (male and female Sprague-Dawley). Janssen R&D data, Table S1.",
      sep = " "
    ),
    notes          = paste(
      "Mean concentrations were modelled rather than individual data because of the small number",
      "of animals per sampling time (Methods 2.1). Rats served as a validation species for the",
      "mouse-to-larger-species extrapolation of healthy-lung disposition: Figure S4 compares",
      "observed and predicted lung AUC0-24 (340.67 vs 502.35 ug*h/g), Cmax (22.26 vs 35.22 ug/g)",
      "and Tmax (8 vs 6.75 h) after a single 20 mg/kg oral dose, and Figure 3 shows the",
      "multiple-dose observed-versus-predicted comparison. All predictions fell within 2-fold of",
      "observed. High doses and intravenous data were excluded to avoid phospholipidosis",
      "(Discussion).",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # 1. Empirical plasma model (Table S6), reported per kilogram and scaled
    #    to the fixed reference rat weight below.
    # =====================================================================
    bw_ref  <- fixed(0.25)  ; label("Reference body weight of a typical rat (kg)")               # Table S3 footnote a: Simcyp rat weight = 0.25 kg
    lka     <- log(0.2)     ; label("Absorption rate constant ka (1/h)")                         # Table S6 ka = 0.2 (RSE 33.7%)
    lcl     <- log(0.51)    ; label("Clearance CL (L/h/kg)")                                     # Table S6 Cl = 0.51 L/h/kg (RSE 21.1%)
    lvc     <- log(1.36)    ; label("Central volume V1 (L/kg)")                                  # Table S6 V1 = 1.36 L/kg (RSE 58.0%)
    lq      <- log(5.84)    ; label("Inter-compartmental clearance Q2 (L/h/kg)")                 # Table S6 Q2 = 5.84 L/h/kg (RSE 48.5%)
    lvp     <- log(3.01)    ; label("First peripheral volume V2 (L/kg)")                         # Table S6 V2 = 3.01 L/kg (RSE 25.6%)
    lq2     <- log(0.052)   ; label("Inter-compartmental clearance Q3 (L/h/kg)")                 # Table S6 Q3 = 0.052 L/h/kg (RSE 53.5%)
    lvp2    <- log(19.86)   ; label("Second peripheral volume V3 (L/kg)")                        # Table S6 V3 = 19.86 L/kg (RSE 46.0%)
    lfdepot <- log(0.65)    ; label("Oral bioavailability F (unitless)")                         # Table S6 F = 0.65 (RSE 21.5%)

    # =====================================================================
    # 2. Drug-related parameters, independent of species (Table S2).
    # =====================================================================
    pka1 <- fixed(8.5) ; label("First (strong) basic pKa (unitless)")   # Table S2 pKa1 = 8.5 (base), Janssen R&D data
    pka2 <- fixed(1.1) ; label("Second (weak) basic pKa (unitless)")    # Table S2 pKa2 = 1.1 (base), ACD predicted

    # =====================================================================
    # 3. Species-related lung physiology for the rat (Table S3).
    # =====================================================================
    qlung       <- fixed(4800)   ; label("Lung blood flow rate Qlung (mL/h)")                    # Table S3 rat blood flow rate = 4800 mL/h
    sa_cells    <- fixed(3900)   ; label("Surface area of lung cells SAcells (cm^2)")            # Table S3 rat surface area lung cells = 3900 cm^2
    v_pb        <- fixed(0.314)  ; label("Volume of pulmonary blood VPB (mL)")                   # Table S3 rat volume of pulmonary blood = 0.314 mL
    v_ew        <- fixed(0.103)  ; label("Volume of lung extracellular water VEW (mL)")          # Table S3 rat volume of EW = 0.103 mL
    v_lysosome  <- fixed(0.0124) ; label("Volume of lung lysosomes VLYS (mL)")                   # Table S3 rat volume of lysosomes = 0.0124 mL
    v_iw        <- fixed(0.811)  ; label("Volume of lung intracellular water VIW (mL)")          # Table S3 rat volume of IW = 0.811 mL
    ph_plasma   <- fixed(7.4)    ; label("Plasma pH (unitless)")                                 # Table S3 rat pH_plasma = 7.4
    ph_ew       <- fixed(7.4)    ; label("Lung extracellular-water pH (unitless)")               # Table S3 rat pH_EW = 7.4
    ph_iw       <- fixed(7.36)   ; label("Lung intracellular-water pH (unitless)")               # Table S3 rat pH_IW = 7.36

    # =====================================================================
    # 4. Drug- AND species-related parameters (Table S4).
    # =====================================================================
    bp        <- fixed(1.0)    ; label("Blood-to-plasma concentration ratio B/P (unitless)")     # Table S4 rat blood to plasma ratio = 1.0
    fu_plasma <- fixed(7.2e-4) ; label("Fraction unbound in plasma fupls (unitless)")            # Table S4 rat fupls = 7.2E-4
    fu_ew     <- fixed(1.1e-3) ; label("Fraction unbound in lung extracellular water fuEW")      # Table S4 rat fuEW = 1.1E-3
    fu_iw     <- fixed(1.0)    ; label("Fraction unbound in lung intracellular water fuIW")      # Supporting Information: "a fixed value of 1 was used for fuIW"

    # =====================================================================
    # 5. Lung parameters carried over from the mouse fit (Table S8).
    #    Methods 2.2: "The cellular permeability value (Permcells) was kept
    #    constant across species considering each species' surface areas",
    #    so PermCells is FIXED here rather than estimated. The lysosomal
    #    product fuLYS*fniLYS is likewise held constant across species and
    #    equation (S9) rescales KpLYS:IW/fuIW to the rat pH_IW in model().
    #    Because rat pH_IW (7.36) equals mouse pH_IW, the rat value recovers
    #    the mouse estimate 3.48e6 exactly.
    # =====================================================================
    lperm_cells   <- fixed(log(55.7))       ; label("Permeability of lung cells PermCells (cm/h, log scale)")  # Table S8 Permcells = 5.57E+1 cm/h, estimated in mice and held constant across species
    lfulys_fnilys <- fixed(log(1.9411e-8))  ; label("Product fuLYS*fniLYS in lysosomes (unitless, log scale)") # derived from Table S8 KpLYS:IW/fuIW = 3.48E+6 in mice via equations (S8)-(S9)

    # =====================================================================
    # 6. Random effects (Table S6).
    # =====================================================================
    etalka  ~ 0.3844   # Table S6 omega_ka = 0.62 (Monolix SD of the log-scale random effect); variance = 0.62^2
    etalcl  ~ 0.09     # Table S6 omega_Cl = 0.3; variance = 0.3^2
    etalq2  ~ 0.6724   # Table S6 omega_Q3 = 0.82; variance = 0.82^2 (the paper's Q3 is this file's q2)

    # =====================================================================
    # 7. Residual error (Table S6). The lung layer was extrapolated rather
    #    than fitted in rats, so no rat lung residual error is reported.
    # =====================================================================
    propSd <- 0.22 ; label("Plasma proportional residual error (fraction)")  # Table S6 b = 0.22 (RSE 9.98%)
  })

  model({
    # 1. Individual plasma parameters, scaled to the reference rat weight.
    ka  <- exp(lka + etalka)
    cl  <- exp(lcl + etalcl) * bw_ref
    vc  <- exp(lvc)          * bw_ref
    q   <- exp(lq)           * bw_ref
    vp  <- exp(lvp)          * bw_ref
    q2  <- exp(lq2 + etalq2) * bw_ref
    vp2 <- exp(lvp2)         * bw_ref

    perm_cells <- exp(lperm_cells)

    # 2. Ionisation. Bedaquiline is a diprotic BASE (Table S2, Figure S8), so
    #    the four-term Henderson-Hasselbalch factor printed in the Supporting
    #    Information is applied in its product form
    #    (1+10^(pKa1-pH))*(1+10^(pKa2-pH)). See vignette Errata.
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

    f(depot) <- exp(lfdepot)

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

    Cc ~ prop(propSd)
  })
}
