Karakitsios_2024_rifampicin_rabbit_pbpk <- function() {
  description <- paste(
    "PBPK (middle-out, permeability-limited lung; Monolix Suite 2021R1).",
    "Preclinical (rabbit). Rifampicin distribution into healthy/uninvolved lung and into",
    "tuberculosis-infected lung tissue (cellular lesion and cavity caseum) in Mycobacterium",
    "tuberculosis-infected New Zealand White rabbits after intravenous infusion. The published",
    "Rifat 2018 one-compartment plasma model drives a lung split into pulmonary blood (in",
    "instantaneous equilibrium with extracellular water), extracellular water and",
    "intracellular water, with permeability-limited, pH-dependent transfer of the unbound",
    "non-ionised drug across the lung-cell membrane. Extracellular water feeds a cellular",
    "lesion by permeability-limited transfer, and the lesion feeds the cavity caseum by",
    "passive diffusion. The lung, lesion and caseum states hold CONCENTRATIONS, not amounts.",
    "This is the species in which the rifampicin lung, lesion and caseum parameters used for",
    "the rabbit-to-human extrapolation were estimated.",
    sep = " "
  )
  reference <- paste(
    "Karakitsios E, Dokoumetzidis A.",
    "Extrapolation of lung pharmacokinetics of antitubercular drugs from preclinical species",
    "to humans using PBPK modelling. J Antimicrob Chemother. 2024;79(6):1362.",
    "doi:10.1093/jac/dkae109.",
    "The plasma layer is the rabbit model of Rifat D, Prideaux B, Savic RM, et al.",
    "Pharmacokinetics of rifapentine and rifampin in a rabbit model of tuberculosis and",
    "correlation with clinical trial data. Sci Transl Med. 2018;10:eaai7786.",
    "doi:10.1126/scitranslmed.aai7786.",
    sep = " "
  )
  vignette <- "Karakitsios_2024_antitubercular_lung_pbpk"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    central = list(analyte = "rifampicin", units = "mg", specimen = "plasma", verified = TRUE),
    lung_ew = list(analyte = "rifampicin", units = "ug/mL", specimen = "tissue", verified = TRUE),
    lung_iw = list(analyte = "rifampicin", units = "ug/mL", specimen = "tissue", verified = TRUE),
    lesion = list(analyte = "rifampicin", units = "ug/mL", specimen = "tissue", verified = TRUE),
    caseum = list(analyte = "rifampicin", units = "ug/mL", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "rabbit (New Zealand White)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    sex_female_pct = NA_real_,
    weight_median = "2.5 kg (reference rabbit of Table S3)",
    disease_state = paste(
      "Mycobacterium tuberculosis-infected rabbits with cavitary pulmonary tuberculosis",
      "(Rifat 2018). Lung samples were dissected into uninvolved lung, cellular lesion, cavity",
      "wall and cavity caseum.",
      sep = " "
    ),
    dose_range = paste(
      "10 mg/kg intravenous infusion over about 10 minutes, single and multiple doses",
      "(Table 1; Rifat 2018).",
      sep = " "
    ),
    notes = paste(
      "Mean (digitised) concentrations were fitted rather than individual data. The cellular",
      "lesion was fitted to the pooled cellular-lesion and cavity-wall data (five data points",
      "after single and multiple dosing; Methods). Lung, lesion and caseum PK profiles were",
      "fitted simultaneously, so the caseum prediction depends on the lesion and lung fits",
      "(Supplementary Figure S1).",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # 1. Plasma layer: Rifat 2018 Table 2 rabbit one-compartment model,
    #    reported per kilogram and scaled to the reference rabbit weight.
    #    Inherited from the upstream publication and held fixed here.
    # =====================================================================
    bw_ref <- fixed(2.5)        ; label("Reference body weight of a typical rabbit (kg)")   # Table S3 rabbit body weight = 2.5 kg
    lcl    <- fixed(log(0.175)) ; label("Plasma clearance CL (L/h/kg)")                     # Rifat 2018 Table 2: CL = 0.175 L/h per kg
    lvc    <- fixed(log(0.54))  ; label("Plasma volume of distribution V (L/kg)")           # Rifat 2018 Table 2: V = 0.54 L/kg

    # =====================================================================
    # 2. Drug-related parameters, independent of species (Table S5).
    #    Rifampicin is an ampholyte.
    # =====================================================================
    pka_acid <- fixed(1.7) ; label("Acidic pKa (unitless)")   # Table S5 RIF pKa1 = 1.7
    pka_base <- fixed(7.9) ; label("Basic pKa (unitless)")    # Table S5 RIF pKa2 = 7.9

    # =====================================================================
    # 3. Species-related lung physiology for the rabbit (Table S3).
    # =====================================================================
    qlung    <- fixed(5.5211 * 3600) ; label("Lung blood flow rate Qlung (mL/h)")             # Table S3 rabbit blood flow rate = 5.5211 mL/sec
    v_pb     <- fixed(3.135)         ; label("Volume of pulmonary blood VPB (mL)")            # Table S3 rabbit volume of pulmonary blood = 3.135 mL
    v_ew     <- fixed(2.577)         ; label("Volume of lung extracellular water VEW (mL)")   # Table S3 rabbit extracellular volume of lung tissue = 2.577 mL
    v_iw     <- fixed(11.288)        ; label("Volume of lung intracellular water VIW (mL)")   # Table S3 rabbit intracellular volume of lung tissue = 11.288 mL
    sa_cells <- fixed(3.36)          ; label("Surface area of lung cells SAcells (m^2)")      # Table S3 rabbit surface area of lung cells = 3.36 m^2
    ph_plasma <- fixed(7.389)        ; label("Plasma pH (unitless)")                          # Table S3 rabbit pH plasma = 7.389
    ph_ew     <- fixed(7.4)          ; label("Lung extracellular-water pH (unitless)")        # Table S3 rabbit pH EW-Lung = 7.4
    ph_iw     <- fixed(7.14)         ; label("Lung intracellular-water pH (unitless)")        # Table S3 rabbit pH IW-Lung = 7.14

    # =====================================================================
    # 4. Drug- and species-related parameters (Table S7).
    # =====================================================================
    bp        <- fixed(0.66178) ; label("Blood-to-plasma concentration ratio B/P (unitless)")  # Table S7 RIF blood to plasma ratio = 0.66178
    fu_plasma <- fixed(0.0605)  ; label("Fraction unbound in plasma fupls (unitless)")         # Table S7 RIF fu plasma = 0.0605 (Sarathy 2016)
    fu_ew     <- fixed(1)       ; label("Fraction unbound in lung extracellular water fuEW")   # Table S7 RIF fu EW = 1; Table 3 footnote d: assumed 1 for moderate-to-strong bases

    # =====================================================================
    # 5. Lung, lesion and caseum parameters estimated in rabbits (Table 3).
    #    Table 3 prints the two rates in L/h; equations (10) and (11) need
    #    them in 1/h, which is how they are used here (see vignette).
    # =====================================================================
    lperm_cells   <- log(3216.69 * 10) ; label("Permeability of lung cells PermCells (cm/h, log scale)")       # Table S7 RIF 3216.69 dm/h = 3.22E+4 cm/h (Table 3, RSE 31.1%)
    lfu_iw        <- log(0.015)        ; label("Fraction unbound in lung intracellular water fuIW (log scale)")  # Table 3 RIF rabbits fuIW = 0.015 (RSE 8.6%)
    lrate_lesion  <- log(1.58e7)       ; label("Rate of distribution into the cellular lesion (1/h, log scale)") # Table 3 RIF rabbits = 1.58E+7 (RSE 15.0%)
    lfnifu_lesion <- log(1.6e-8)       ; label("Product fni x fu in the cellular lesion (unitless, log scale)")   # Table 3 RIF rabbits = 1.6E-8 (RSE 12.2%)
    lrate_caseum  <- log(0.22)         ; label("Rate of distribution into the caseum (1/h, log scale)")          # Table 3 RIF rabbits = 0.22 (RSE 23.8%)
    lpc_caseum    <- log(1.5)          ; label("Caseum-to-lesion partition ratio PC (unitless, log scale)")      # Table 3 RIF rabbits ratio of distribution in caseum = 1.5 (RSE 32.9%)
  })

  model({
    # 1. Plasma layer (Rifat 2018): one compartment, intravenous infusion
    #    into central, parameters per kg times the reference body weight.
    cl <- exp(lcl) * bw_ref
    vc <- exp(lvc) * bw_ref

    d/dt(central) <- -cl / vc * central

    # 2. Ionisation. The tabulated non-ionised fractions of Tables S5-S8 are
    #    reproduced by the product form (1 + 10^(pH - pKa_acid)) *
    #    (1 + 10^(pKa_base - pH)) for this ampholyte.
    x_plasma <- (1 + 10^(ph_plasma - pka_acid)) * (1 + 10^(pka_base - ph_plasma))
    x_ew     <- (1 + 10^(ph_ew - pka_acid))     * (1 + 10^(pka_base - ph_ew))
    x_iw     <- (1 + 10^(ph_iw - pka_acid))     * (1 + 10^(pka_base - ph_iw))
    fni_plasma <- 1 / x_plasma
    fni_ew     <- 1 / x_ew
    fni_iw     <- 1 / x_iw

    perm_cells  <- exp(lperm_cells)
    fu_iw       <- exp(lfu_iw)
    rate_lesion <- exp(lrate_lesion)
    fnifu_lesion <- exp(lfnifu_lesion)
    rate_caseum <- exp(lrate_caseum)
    pc_caseum   <- exp(lpc_caseum)

    # 3. Lung constants. K_EW:B is total EW over total blood concentration
    #    when the unbound non-ionised concentrations of plasma and EW are
    #    equal (instantaneous PB-EW equilibrium). PS converts m^2 to cm^2,
    #    giving mL/h.
    kewb     <- fu_plasma * fni_plasma / (fu_ew * fni_ew * bp)
    ps_cells <- perm_cells * sa_cells * 1e4
    v_ew_eff <- v_ew + v_pb / kewb

    cblood <- central / vc * bp
    jcells <- ps_cells * (lung_ew * fni_ew - fu_iw * lung_iw * fni_iw)

    # 4. Healthy/uninvolved lung, equations (1) and (2). lung_ew holds the
    #    UNBOUND EW concentration Cu_EW; lung_iw holds the TOTAL IW
    #    concentration.
    d/dt(lung_ew) <- fu_ew / v_ew_eff *
      (qlung * cblood - qlung * lung_ew / (fu_ew * kewb) - jcells)
    d/dt(lung_iw) <- jcells / v_iw

    # 5. Infected lung, equations (10) and (11).
    d/dt(lesion) <- rate_lesion * (lung_ew * fni_ew - fnifu_lesion * lesion)
    d/dt(caseum) <- rate_caseum * (pc_caseum * lesion - caseum)

    # 6. Observations. The uninvolved-lung homogenate concentration is the
    #    volume-weighted mean of pulmonary blood, EW and IW (Table S1).
    Cc <- central / vc

    c_lung_ew <- lung_ew / fu_ew
    c_lung_pb <- c_lung_ew / kewb
    c_lung <- (c_lung_pb * v_pb + c_lung_ew * v_ew + lung_iw * v_iw) /
      (v_pb + v_ew + v_iw)
    c_lesion <- lesion
    c_caseum <- caseum
  })
}
