Karakitsios_2024_moxifloxacin_rabbit_pbpk <- function() {
  description <- paste(
    "PBPK (middle-out, permeability-limited lung; Monolix Suite 2021R1).",
    "Preclinical (rabbit). Moxifloxacin distribution into healthy/uninvolved lung and into",
    "tuberculosis-infected lung tissue (cellular lesion and caseum) in Mycobacterium",
    "tuberculosis-infected rabbits after oral dosing. The published Pienaar 2017",
    "two-compartment plasma model with first-order absorption drives a lung split into",
    "pulmonary blood, extracellular water and intracellular water. Because moxifloxacin is a",
    "P-glycoprotein substrate, pulmonary blood is a separate state with permeability-limited",
    "capillary transfer to extracellular water instead of the instantaneous equilibrium used",
    "for the other drugs; transfer across the lung-cell membrane is also permeability limited",
    "and pH dependent. Extracellular water feeds a cellular lesion, which feeds the caseum by",
    "passive diffusion. The lung, lesion and caseum states hold CONCENTRATIONS, not amounts.",
    sep = " "
  )
  reference <- paste(
    "Karakitsios E, Dokoumetzidis A.",
    "Extrapolation of lung pharmacokinetics of antitubercular drugs from preclinical species",
    "to humans using PBPK modelling. J Antimicrob Chemother. 2024;79(6):1362.",
    "doi:10.1093/jac/dkae109.",
    "The plasma layer is the rabbit moxifloxacin model of Pienaar E, Sarathy J, Prideaux B,",
    "et al. Comparing efficacies of moxifloxacin, levofloxacin and gatifloxacin in",
    "tuberculosis granulomas using a multi-scale systems pharmacology approach.",
    "PLoS Comput Biol. 2017;13:e1005650. doi:10.1371/journal.pcbi.1005650.",
    sep = " "
  )
  vignette <- "Karakitsios_2024_antitubercular_lung_pbpk"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")
  dosing <- c("depot", "central")

  compartmentData <- list(
    depot = list(analyte = "moxifloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "moxifloxacin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "moxifloxacin", units = "mg", specimen = "plasma", verified = TRUE),
    vp_lung = list(analyte = "moxifloxacin", units = "ug/mL", specimen = "whole blood", verified = TRUE),
    lung_ew = list(analyte = "moxifloxacin", units = "ug/mL", specimen = "tissue", verified = TRUE),
    lung_iw = list(analyte = "moxifloxacin", units = "ug/mL", specimen = "tissue", verified = TRUE),
    lesion = list(analyte = "moxifloxacin", units = "ug/mL", specimen = "tissue", verified = TRUE),
    caseum = list(analyte = "moxifloxacin", units = "ug/mL", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "rabbit (New Zealand White)",
    n_subjects = NA_integer_,
    n_studies = 2L,
    sex_female_pct = NA_real_,
    weight_median = "2.5 kg (reference rabbit of Table S3)",
    disease_state = paste(
      "Mycobacterium tuberculosis-infected rabbits with cellular and caseous pulmonary",
      "granulomas and cavities (Pienaar 2017; Sarathy 2019).",
      sep = " "
    ),
    dose_range = "100 mg/kg oral, single dose (Table 1).",
    notes = paste(
      "Mean (digitised) concentrations were fitted. Cavity-wall data were pooled with",
      "cellular-lesion data (4 data points), and laser-capture-microdissection data for normal",
      "lung, cellular areas and caseum of caseous granulomas were added to enrich the data",
      "(Supplementary Section 5, Figure S3).",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # 1. Plasma layer: Pienaar 2017 Table 3 rabbit two-compartment oral
    #    model, reported per kilogram and scaled to the reference rabbit.
    #    Apparent (F-relative) values; inherited and held constant here.
    # =====================================================================
    bw_ref  <- fixed(2.5)       ; label("Reference body weight of a typical rabbit (kg)")       # Table S3 rabbit body weight = 2.5 kg
    lka     <- fixed(log(0.9))  ; label("Absorption rate constant ka (1/h)")                     # Pienaar 2017 Table 3 rabbit MXF ka = 0.9 1/h
    lcl     <- fixed(log(1.3))  ; label("Plasma clearance CL (L/h/kg)")                          # Pienaar 2017 Table 3 rabbit MXF CL = 1.3 L/h/kg
    lvc     <- fixed(log(5.35)) ; label("Plasma volume of distribution Vp (L/kg)")               # Pienaar 2017 Table 3 rabbit MXF Vp = 5.35 L/kg
    lq      <- fixed(log(2.86)) ; label("Inter-compartmental clearance Q (L/h/kg)")              # Pienaar 2017 Table 3 rabbit MXF Q = 2.86 L/h/kg
    lvp     <- fixed(log(80.8)) ; label("Peripheral volume of distribution Vpe (L/kg)")          # Pienaar 2017 Table 3 rabbit MXF Vpe = 80.8 L/kg

    # =====================================================================
    # 2. Drug-related parameters, independent of species (Table S5).
    # =====================================================================
    pka_acid <- fixed(6.25) ; label("Acidic pKa (unitless)")   # Table S5 MXF pKa1 = 6.25
    pka_base <- fixed(9.29) ; label("Basic pKa (unitless)")    # Table S5 MXF pKa2 = 9.29

    # =====================================================================
    # 3. Species-related lung physiology for the rabbit (Table S3).
    # =====================================================================
    qlung     <- fixed(5.5211 * 3600) ; label("Lung blood flow rate Qlung (mL/h)")              # Table S3 rabbit blood flow rate = 5.5211 mL/sec
    v_pb      <- fixed(3.135)         ; label("Volume of pulmonary blood VPB (mL)")             # Table S3 rabbit volume of pulmonary blood = 3.135 mL
    v_ew      <- fixed(2.577)         ; label("Volume of lung extracellular water VEW (mL)")    # Table S3 rabbit extracellular volume of lung tissue = 2.577 mL
    v_iw      <- fixed(11.288)        ; label("Volume of lung intracellular water VIW (mL)")    # Table S3 rabbit intracellular volume of lung tissue = 11.288 mL
    sa_cells  <- fixed(3.36)          ; label("Surface area of lung cells SAcells (m^2)")       # Table S3 rabbit surface area of lung cells = 3.36 m^2
    sa_capillaries <- fixed(3.10)     ; label("Surface area of lung capillaries SAcap (m^2)")   # Table S3 rabbit surface area of capillaries in lung = 3.10 m^2
    ph_plasma <- fixed(7.389)         ; label("Plasma pH (unitless)")                           # Table S3 rabbit pH plasma = 7.389
    ph_ew     <- fixed(7.4)           ; label("Lung extracellular-water pH (unitless)")         # Table S3 rabbit pH EW-Lung = 7.4
    ph_iw     <- fixed(7.14)          ; label("Lung intracellular-water pH (unitless)")         # Table S3 rabbit pH IW-Lung = 7.14

    # =====================================================================
    # 4. Drug- and species-related parameters (Table S7).
    # =====================================================================
    bp        <- fixed(0.890308) ; label("Blood-to-plasma concentration ratio B/P (unitless)")  # Table S7 MXF blood to plasma ratio = 0.890308
    fu_plasma <- fixed(0.695)    ; label("Fraction unbound in plasma fupls (unitless)")         # Table S7 MXF fu plasma = 0.695 (Sarathy 2016)
    fu_ew     <- fixed(1)        ; label("Fraction unbound in lung extracellular water fuEW")   # Table S7 MXF fu EW = 1; Table 3 footnote d

    # =====================================================================
    # 5. Lung, lesion and caseum parameters estimated in rabbits (Table 3).
    #    The two rates are printed in L/h but enter equations (10)-(11) as
    #    1/h (see vignette).
    # =====================================================================
    lperm_cells       <- log(0.99 * 10) ; label("Permeability of lung cells PermCells (cm/h, log scale)")         # Table S7 MXF 0.99 dm/h = 9.9 cm/h (Table 3, RSE 0.10%)
    lperm_capillaries <- log(0.6 * 10)  ; label("Permeability of lung capillaries PermCap (cm/h, log scale)")      # Table S7 MXF 0.6 dm/h = 6.0 cm/h (Table 3, RSE 19.7%)
    lfu_iw            <- log(0.056)     ; label("Fraction unbound in lung intracellular water fuIW (log scale)")   # Table 3 MXF rabbits fuIW = 0.056 (RSE 20.4%)
    lrate_lesion      <- log(2.19e4)    ; label("Rate of distribution into the cellular lesion (1/h, log scale)")  # Table 3 MXF rabbits = 2.19E+4 (RSE 34.5%)
    lfnifu_lesion     <- log(4.7e-5)    ; label("Product fni x fu in the cellular lesion (unitless, log scale)")    # Table 3 MXF rabbits = 4.7E-5 (RSE 22.7%)
    lrate_caseum      <- log(9.64)      ; label("Rate of distribution into the caseum (1/h, log scale)")           # Table 3 MXF rabbits = 9.64 (RSE 135%)
    lpc_caseum        <- log(0.49)      ; label("Caseum-to-lesion partition ratio PC (unitless, log scale)")       # Table 3 MXF rabbits ratio of distribution in caseum = 0.49 (RSE 14.0%)
  })

  model({
    # 1. Plasma layer (Pienaar 2017): two compartments, first-order
    #    absorption, parameters per kg times the reference body weight.
    ka <- exp(lka)
    cl <- exp(lcl) * bw_ref
    vc <- exp(lvc) * bw_ref
    q  <- exp(lq) * bw_ref
    vp <- exp(lvp) * bw_ref

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl / vc * central -
      q / vc * central + q / vp * peripheral1
    d/dt(peripheral1) <- q / vc * central - q / vp * peripheral1

    # 2. Ionisation (ampholyte; product form reproduces Table S7 fni).
    x_plasma <- (1 + 10^(ph_plasma - pka_acid)) * (1 + 10^(pka_base - ph_plasma))
    x_ew     <- (1 + 10^(ph_ew - pka_acid))     * (1 + 10^(pka_base - ph_ew))
    x_iw     <- (1 + 10^(ph_iw - pka_acid))     * (1 + 10^(pka_base - ph_iw))
    fni_plasma <- 1 / x_plasma
    fni_ew     <- 1 / x_ew
    fni_iw     <- 1 / x_iw

    perm_cells       <- exp(lperm_cells)
    perm_capillaries <- exp(lperm_capillaries)
    fu_iw        <- exp(lfu_iw)
    rate_lesion  <- exp(lrate_lesion)
    fnifu_lesion <- exp(lfnifu_lesion)
    rate_caseum  <- exp(lrate_caseum)
    pc_caseum    <- exp(lpc_caseum)

    # 3. Permeability-surface products (m^2 -> cm^2, giving mL/h).
    ps_cells       <- perm_cells * sa_cells * 1e4
    ps_capillaries <- perm_capillaries * sa_capillaries * 1e4

    cblood <- central / vc * bp

    # Capillary flux PB -> EW of unbound non-ionised drug: PB blood
    # concentration divided by B/P gives the PB plasma concentration.
    jcap   <- ps_capillaries *
      (fu_plasma * fni_plasma * vp_lung / bp - lung_ew * fni_ew)
    jcells <- ps_cells * (lung_ew * fni_ew - fu_iw * lung_iw * fni_iw)

    # 4. Lung, equations (3), (4) and (2). Equations (3)-(4) are encoded
    #    in their mass-conserving form: PB is perfused (arterial in, venous
    #    out at the PB concentration) and the capillary flux leaves PB and
    #    enters EW (see vignette). lung_ew holds the UNBOUND EW
    #    concentration; vp_lung and lung_iw hold TOTAL concentrations.
    d/dt(vp_lung) <- (qlung * (cblood - vp_lung) - jcap) / v_pb
    d/dt(lung_ew) <- fu_ew / v_ew * (jcap - jcells)
    d/dt(lung_iw) <- jcells / v_iw

    # 5. Infected lung, equations (10) and (11).
    d/dt(lesion) <- rate_lesion * (lung_ew * fni_ew - fnifu_lesion * lesion)
    d/dt(caseum) <- rate_caseum * (pc_caseum * lesion - caseum)

    # 6. Observations: volume-weighted uninvolved-lung homogenate (Table S1).
    Cc <- central / vc

    c_lung_ew <- lung_ew / fu_ew
    c_lung <- (vp_lung * v_pb + c_lung_ew * v_ew + lung_iw * v_iw) /
      (v_pb + v_ew + v_iw)
    c_lesion <- lesion
    c_caseum <- caseum
  })
}
