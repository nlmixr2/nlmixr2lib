Karakitsios_2024_moxifloxacin_human_fromRabbit_pbpk <- function() {
  description <- paste(
    "PBPK (middle-out, permeability-limited lung; Monolix Suite 2021R1).",
    "Human extrapolation, from rabbits, of moxifloxacin distribution into healthy/uninvolved",
    "lung and into tuberculosis-infected lung tissue (cellular lesion and cavity caseum) in a",
    "typical adult with pulmonary tuberculosis. The published Strydom 2019 one-compartment",
    "plasma model with first-order absorption drives a lung split into pulmonary blood,",
    "extracellular water and intracellular water. Because moxifloxacin is a P-glycoprotein",
    "substrate, pulmonary blood is a separate state with permeability-limited capillary",
    "transfer to extracellular water; transfer across the lung-cell membrane is also",
    "permeability limited and pH dependent. Extracellular water feeds a cellular lesion, which",
    "feeds the caseum by passive diffusion. Capillary and cell permeabilities, lesion and",
    "caseum rates and the caseum partition ratio are carried over unchanged from the rabbit",
    "fit; human physiology, protein binding and ionisation replace the rabbit values. The",
    "lung, lesion and caseum states hold CONCENTRATIONS, not amounts. Typical-value model (no",
    "between-subject variability).",
    sep = " "
  )
  reference <- paste(
    "Karakitsios E, Dokoumetzidis A.",
    "Extrapolation of lung pharmacokinetics of antitubercular drugs from preclinical species",
    "to humans using PBPK modelling. J Antimicrob Chemother. 2024;79(6):1362.",
    "doi:10.1093/jac/dkae109.",
    "The plasma layer is the moxifloxacin model of Strydom N, Gupta SV, Fox WS, et al.",
    "Tuberculosis drugs' distribution and emergence of resistance in patient's lung lesions:",
    "a mechanistic model and tool for regimen and dose optimization. PLoS Med.",
    "2019;16:e1002773. doi:10.1371/journal.pmed.1002773.",
    sep = " "
  )
  vignette <- "Karakitsios_2024_antitubercular_lung_pbpk"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    depot = list(analyte = "moxifloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "moxifloxacin", units = "mg", specimen = "plasma", verified = TRUE),
    vp_lung = list(analyte = "moxifloxacin", units = "ug/mL", specimen = "whole blood", verified = TRUE),
    lung_ew = list(analyte = "moxifloxacin", units = "ug/mL", specimen = "tissue", verified = TRUE),
    lung_iw = list(analyte = "moxifloxacin", units = "ug/mL", specimen = "tissue", verified = TRUE),
    lesion = list(analyte = "moxifloxacin", units = "ug/mL", specimen = "tissue", verified = TRUE),
    caseum = list(analyte = "moxifloxacin", units = "ug/mL", specimen = "tissue", verified = TRUE)
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
    dose_range = "400 mg oral, single dose (Table 1, Table 4); multiple daily 400 mg doses in Figure 4.",
    regions = "South Korea (Strydom 2019 clinical study NCT00816426).",
    notes = paste(
      "The extrapolated profiles represent a typical mean patient (Methods). Predictions were",
      "compared with the Strydom 2019 empirical uninvolved-lung, cavity-wall and cavity-caseum",
      "models (Table 4 'Observed'); the cavity wall was matched to the model's cellular lesion",
      "(Supplementary Section 5).",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # 1. Plasma layer: Strydom 2019 Table 2 moxifloxacin, one compartment
    #    with first-order absorption. Apparent values, held constant.
    # =====================================================================
    lka <- fixed(log(1.55)) ; label("Absorption rate constant ka (1/h)")   # Strydom 2019 Table 2 MFX rate of absorption = 1.55 1/h
    lcl <- fixed(log(8.93)) ; label("Apparent clearance CL/F (L/h)")       # Strydom 2019 Table 2 MFX clearance = 8.93 L/h
    lvc <- fixed(log(147))  ; label("Apparent central volume V/F (L)")     # Strydom 2019 Table 2 MFX central volume = 147 L

    # =====================================================================
    # 2. Drug-related parameters, independent of species (Table S5).
    # =====================================================================
    pka_acid <- fixed(6.25) ; label("Acidic pKa (unitless)")   # Table S5 MXF pKa1 = 6.25
    pka_base <- fixed(9.29) ; label("Basic pKa (unitless)")    # Table S5 MXF pKa2 = 9.29

    # =====================================================================
    # 3. Species-related lung physiology for humans (Table S4).
    # =====================================================================
    qlung     <- fixed(356 * 1000) ; label("Lung blood flow rate Qlung (mL/h)")             # Table S4 human blood flow rate = 356 L/hr
    v_pb      <- fixed(101.195)    ; label("Volume of pulmonary blood VPB (mL)")            # Table S4 human volume of pulmonary blood = 101.195 mL
    v_ew      <- fixed(89.161)     ; label("Volume of lung extracellular water VEW (mL)")   # Table S4 human extracellular volume of lung tissue = 89.161 mL
    v_iw      <- fixed(356.644)    ; label("Volume of lung intracellular water VIW (mL)")   # Table S4 human intracellular volume of lung tissue = 356.644 mL
    sa_cells  <- fixed(96)         ; label("Surface area of lung cells SAcells (m^2)")      # Table S4 human surface area of lung cells = 96 m^2
    sa_capillaries <- fixed(72.3)  ; label("Surface area of lung capillaries SAcap (m^2)")  # Table S4 human surface area of capillaries in lung = 72.3 m^2
    ph_plasma <- fixed(7.4)        ; label("Plasma pH (unitless)")                          # Table S4 human pH plasma = 7.4
    ph_ew     <- fixed(7.4)        ; label("Lung extracellular-water pH (unitless)")        # Table S4 human pH EW-Lung = 7.4
    ph_iw     <- fixed(6.69)       ; label("Lung intracellular-water pH (unitless)")        # Table S4 human pH IW-Lung = 6.69
    ph_lesion <- fixed(5.84)       ; label("Cellular-lesion pH (unitless)")                 # Supplementary Section 1: mean pH 5.84 in human cellular lesions

    # =====================================================================
    # 4. Drug- and species-related parameters (Tables S8 and S12).
    # =====================================================================
    bp        <- fixed(1.045)  ; label("Blood-to-plasma concentration ratio B/P (unitless)")  # Table S8 MXF blood to plasma ratio = 1.045
    fu_plasma <- fixed(0.7168) ; label("Fraction unbound in plasma fupls (unitless)")         # Table S8 MXF fu plasma = 0.7168
    fu_ew     <- fixed(1)      ; label("Fraction unbound in lung extracellular water fuEW")   # Table S8 / Table S12 MXF fu EW = 1 (assumed)
    fu_iw     <- fixed(0.191)  ; label("Fraction unbound in lung intracellular water fuIW")   # Table S8 / Table S12 MXF fu IW = 0.191 (equations 5-6)
    fu_lesion <- fixed(0.757)  ; label("Fraction unbound in the cellular lesion fuLesion")    # Table S12 MXF fu cellular lesion = 0.757 (equation 14)
    fu_caseum <- fixed(0.168)  ; label("Fraction unbound in caseum fuCaseum")                 # Sarathy 2016 Supplementary Table 1 MXF caseum fu = 16.8% (Supplementary Section 1)

    # =====================================================================
    # 5. Distribution parameters estimated in rabbits and held constant
    #    across species (Table 3, Table S8). Rates enter as 1/h.
    # =====================================================================
    lperm_cells       <- fixed(log(0.99 * 10))  ; label("Permeability of lung cells PermCells (cm/h, log scale)")        # Table S8 MXF 0.99 dm/h = 9.9 cm/h
    lperm_capillaries <- fixed(log(0.604 * 10)) ; label("Permeability of lung capillaries PermCap (cm/h, log scale)")     # Table S8 MXF 0.604 dm/h = 6.04 cm/h
    lrate_lesion      <- fixed(log(2.19e4))     ; label("Rate of distribution into the cellular lesion (1/h, log scale)") # Table 3 MXF rabbits = 2.19E+4, kept constant across species
    lrate_caseum      <- fixed(log(9.64))       ; label("Rate of distribution into the caseum (1/h, log scale)")          # Table 3 MXF rabbits = 9.64, kept constant across species
    lpc_caseum        <- fixed(log(0.49))       ; label("Caseum-to-lesion partition ratio PC (unitless, log scale)")      # Table 3 MXF rabbits = 0.49, kept constant across species
  })

  model({
    # 1. Plasma layer (Strydom 2019): one compartment, first-order absorption.
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)

    d/dt(depot)   <- -ka * depot
    d/dt(central) <- ka * depot - cl / vc * central

    # 2. Ionisation (ampholyte; product form reproduces Table S8 fni).
    x_plasma <- (1 + 10^(ph_plasma - pka_acid)) * (1 + 10^(pka_base - ph_plasma))
    x_ew     <- (1 + 10^(ph_ew - pka_acid))     * (1 + 10^(pka_base - ph_ew))
    x_iw     <- (1 + 10^(ph_iw - pka_acid))     * (1 + 10^(pka_base - ph_iw))
    x_lesion <- (1 + 10^(ph_lesion - pka_acid)) * (1 + 10^(pka_base - ph_lesion))
    fni_plasma <- 1 / x_plasma
    fni_ew     <- 1 / x_ew
    fni_iw     <- 1 / x_iw
    fni_lesion <- 1 / x_lesion

    perm_cells       <- exp(lperm_cells)
    perm_capillaries <- exp(lperm_capillaries)
    rate_lesion <- exp(lrate_lesion)
    rate_caseum <- exp(lrate_caseum)
    pc_caseum   <- exp(lpc_caseum)

    # Human lesion product fni x fu: Table S12 fu at the lesion pH.
    fnifu_lesion <- fni_lesion * fu_lesion

    # 3. Permeability-surface products (m^2 -> cm^2, giving mL/h).
    ps_cells       <- perm_cells * sa_cells * 1e4
    ps_capillaries <- perm_capillaries * sa_capillaries * 1e4

    cblood <- central / vc * bp

    # Capillary flux PB -> EW of unbound non-ionised drug: PB blood
    # concentration divided by B/P gives the PB plasma concentration.
    jcap   <- ps_capillaries *
      (fu_plasma * fni_plasma * vp_lung / bp - lung_ew * fni_ew)
    jcells <- ps_cells * (lung_ew * fni_ew - fu_iw * lung_iw * fni_iw)

    # 4. Lung, equations (3), (4) and (2), in mass-conserving form (see
    #    vignette). lung_ew holds the UNBOUND EW concentration; vp_lung and
    #    lung_iw hold TOTAL concentrations.
    d/dt(vp_lung) <- (qlung * (cblood - vp_lung) - jcap) / v_pb
    d/dt(lung_ew) <- fu_ew / v_ew * (jcap - jcells)
    d/dt(lung_iw) <- jcells / v_iw

    # 5. Infected lung, equations (10) and (11).
    d/dt(lesion) <- rate_lesion * (lung_ew * fni_ew - fnifu_lesion * lesion)
    d/dt(caseum) <- rate_caseum * (pc_caseum * lesion - caseum)

    # 6. Observations: total (homogenate) and unbound concentrations.
    Cc <- central / vc

    c_lung_ew <- lung_ew / fu_ew
    c_lung <- (vp_lung * v_pb + c_lung_ew * v_ew + lung_iw * v_iw) /
      (v_pb + v_ew + v_iw)
    c_lesion <- lesion
    c_caseum <- caseum

    cu_plasma <- Cc * fu_plasma
    cu_ew     <- lung_ew
    cu_lesion <- lesion * fu_lesion
    cu_caseum <- caseum * fu_caseum
  })
}
