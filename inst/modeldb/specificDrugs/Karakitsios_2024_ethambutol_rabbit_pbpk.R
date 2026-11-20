Karakitsios_2024_ethambutol_rabbit_pbpk <- function() {
  description <- paste(
    "PBPK (middle-out, permeability-limited lung; Monolix Suite 2021R1).",
    "Preclinical (rabbit). Ethambutol distribution into healthy/uninvolved lung in",
    "Mycobacterium tuberculosis-infected rabbits after oral dosing. A two-compartment plasma",
    "model with first-order absorption, fitted by the authors to the Zimmerman 2017 rabbit",
    "plasma data (Table S10), drives a lung split into pulmonary blood (in instantaneous",
    "equilibrium with extracellular water), extracellular water and intracellular water, with",
    "permeability-limited, pH-dependent transfer of the unbound non-ionised drug across the",
    "lung-cell membrane. The lung states hold CONCENTRATIONS, not amounts.",
    sep = " "
  )
  reference <- paste(
    "Karakitsios E, Dokoumetzidis A.",
    "Extrapolation of lung pharmacokinetics of antitubercular drugs from preclinical species",
    "to humans using PBPK modelling. J Antimicrob Chemother. 2024;79(6):1362.",
    "doi:10.1093/jac/dkae109.",
    "Plasma and lung data from Zimmerman M, Lestner J, Prideaux B, et al. Ethambutol",
    "partitioning in tuberculous pulmonary lesions explains its clinical efficacy.",
    "Antimicrob Agents Chemother. 2017;61:e00924-17. doi:10.1128/AAC.00924-17.",
    sep = " "
  )
  vignette <- "Karakitsios_2024_antitubercular_lung_pbpk"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")
  dosing <- c("depot", "central")

  compartmentData <- list(
    depot = list(analyte = "ethambutol", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ethambutol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ethambutol", units = "mg", specimen = "plasma", verified = TRUE),
    lung_ew = list(analyte = "ethambutol", units = "ug/mL", specimen = "tissue", verified = TRUE),
    lung_iw = list(analyte = "ethambutol", units = "ug/mL", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "rabbit (New Zealand White)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    sex_female_pct = NA_real_,
    weight_median = "2.5 kg (reference rabbit of Table S3)",
    disease_state = "Mycobacterium tuberculosis-infected rabbits (Zimmerman 2017); uninvolved lung only.",
    dose_range = "100 mg/kg oral, single and multiple doses (Table 1; Figure S4).",
    notes = paste(
      "The plasma model (Table S10) was fitted by the authors in Monolix 2021R1 to the",
      "Zimmerman 2017 plasma data after single and multiple oral 100 mg/kg doses; it has",
      "between-animal variability and a proportional residual error. Its clearances and",
      "volumes are printed in mL/h and mL without a per-kg qualifier and are read here as",
      "absolute values for the 2.5 kg reference rabbit (see vignette).",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # 1. Plasma layer: Table S10, the authors' rabbit two-compartment
    #    model (apparent, F-relative values).
    # =====================================================================
    lka <- log(0.76)         ; label("Absorption rate constant ka (1/h)")           # Table S10 Ka = 0.76 1/h (RSE 38.4%)
    lcl <- log(23503.5 / 1000) ; label("Plasma clearance CL (L/h)")                 # Table S10 Cl = 23503.5 mL/h (RSE 11.7%)
    lvc <- log(18280.47 / 1000) ; label("Central volume V1 (L)")                    # Table S10 V1 = 18280.47 mL (RSE 56.0%)
    lq  <- log(7643.31 / 1000) ; label("Inter-compartmental clearance Q (L/h)")     # Table S10 Q = 7643.31 mL/h (RSE 52.2%)
    lvp <- log(66960.95 / 1000) ; label("Peripheral volume V2 (L)")                 # Table S10 V2 = 66960.95 mL (RSE 12.0%)

    # =====================================================================
    # 2. Drug-related parameters, independent of species (Table S5).
    #    Ethambutol is a diprotic base.
    # =====================================================================
    pka1 <- fixed(9.55) ; label("First basic pKa (unitless)")    # Table S5 ETH pKa1 = 9.55 (base)
    pka2 <- fixed(6.5)  ; label("Second basic pKa (unitless)")   # Table S5 ETH pKa2 = 6.5 (base)

    # =====================================================================
    # 3. Species-related lung physiology for the rabbit (Table S3).
    # =====================================================================
    qlung     <- fixed(5.5211 * 3600) ; label("Lung blood flow rate Qlung (mL/h)")              # Table S3 rabbit blood flow rate = 5.5211 mL/sec
    v_pb      <- fixed(3.135)         ; label("Volume of pulmonary blood VPB (mL)")             # Table S3 rabbit volume of pulmonary blood = 3.135 mL
    v_ew      <- fixed(2.577)         ; label("Volume of lung extracellular water VEW (mL)")    # Table S3 rabbit extracellular volume of lung tissue = 2.577 mL
    v_iw      <- fixed(11.288)        ; label("Volume of lung intracellular water VIW (mL)")    # Table S3 rabbit intracellular volume of lung tissue = 11.288 mL
    sa_cells  <- fixed(3.36)          ; label("Surface area of lung cells SAcells (m^2)")       # Table S3 rabbit surface area of lung cells = 3.36 m^2
    ph_plasma <- fixed(7.389)         ; label("Plasma pH (unitless)")                           # Table S3 rabbit pH plasma = 7.389
    ph_ew     <- fixed(7.4)           ; label("Lung extracellular-water pH (unitless)")         # Table S3 rabbit pH EW-Lung = 7.4
    ph_iw     <- fixed(7.14)          ; label("Lung intracellular-water pH (unitless)")         # Table S3 rabbit pH IW-Lung = 7.14

    # =====================================================================
    # 4. Drug- and species-related parameters (Table S7).
    # =====================================================================
    bp        <- fixed(0.92872) ; label("Blood-to-plasma concentration ratio B/P (unitless)")  # Table S7 ETH blood to plasma ratio = 0.92872
    fu_plasma <- fixed(0.802)   ; label("Fraction unbound in plasma fupls (unitless)")         # Table S7 ETH fu plasma = 0.802 (Sarathy 2016)
    fu_ew     <- fixed(1)       ; label("Fraction unbound in lung extracellular water fuEW")   # Table S7 ETH fu EW = 1; Table 3 footnote d

    # =====================================================================
    # 5. Lung parameters estimated in rabbits (Table 3).
    # =====================================================================
    lperm_cells <- log(0.023 * 10) ; label("Permeability of lung cells PermCells (cm/h, log scale)")       # Table S7 ETH 0.023 dm/h = 0.23 cm/h (Table 3, RSE 18.8%)
    lfu_iw      <- log(0.08)       ; label("Fraction unbound in lung intracellular water fuIW (log scale)")  # Table 3 ETH rabbits fuIW = 0.08 (RSE 10.5%)

    # =====================================================================
    # 6. Plasma between-animal variability and residual error (Table S10).
    #    Monolix omegas are SDs of the log-scale random effect.
    # =====================================================================
    etalcl ~ 0.1156   # Table S10 omega_Cl = 0.34; variance = 0.34^2
    etalvc ~ 0.0225   # Table S10 omega_V1 = 0.15; variance = 0.15^2
    etalq  ~ 1.5129   # Table S10 omega_Q = 1.23; variance = 1.23^2

    propSd <- 0.18 ; label("Plasma proportional residual error (fraction)")  # Table S10 b = 0.18 (RSE 16.4%)
  })

  model({
    # 1. Plasma layer (Table S10): two compartments, first-order absorption.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq + etalq)
    vp <- exp(lvp)

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl / vc * central -
      q / vc * central + q / vp * peripheral1
    d/dt(peripheral1) <- q / vc * central - q / vp * peripheral1

    # 2. Ionisation (diprotic base; product form reproduces Table S7 fni).
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

    cblood <- central / vc * bp
    jcells <- ps_cells * (lung_ew * fni_ew - fu_iw * lung_iw * fni_iw)

    # 4. Healthy/uninvolved lung, equations (1) and (2). lung_ew holds the
    #    UNBOUND EW concentration; lung_iw holds the TOTAL IW concentration.
    d/dt(lung_ew) <- fu_ew / v_ew_eff *
      (qlung * cblood - qlung * lung_ew / (fu_ew * kewb) - jcells)
    d/dt(lung_iw) <- jcells / v_iw

    # 5. Observations: volume-weighted uninvolved-lung homogenate (Table S1).
    Cc <- central / vc

    c_lung_ew <- lung_ew / fu_ew
    c_lung_pb <- c_lung_ew / kewb
    c_lung <- (c_lung_pb * v_pb + c_lung_ew * v_ew + lung_iw * v_iw) /
      (v_pb + v_ew + v_iw)
    c_lung_iw <- lung_iw

    Cc ~ prop(propSd)
  })
}
