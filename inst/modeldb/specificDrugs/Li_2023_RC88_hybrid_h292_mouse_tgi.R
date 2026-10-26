Li_2023_RC88_hybrid_h292_mouse_tgi <- function() {
  description <- paste(
    "Preclinical (mouse, Balb/c nude bearing H292-MSLN-9C8 lung-cancer",
    "xenografts). Hybrid tumor-growth-inhibition model of",
    "the mesothelin-targeting antibody-drug conjugate RC88: exponential,",
    "linear and logistic (carrying-capacity) growth combined, a saturable",
    "(Kmax, KC50) kill term in serum RC88 concentration, and three transit",
    "compartments of dying cells, driven by the fixed mouse two-compartment",
    "PK model with parallel linear and Michaelis-Menten elimination (Li 2023)."
  )
  reference <- paste(
    "Li Q, Wang L, Zhang J, Zhao G, Liu Z, Ma X, Jiang J.",
    "Translation of the efficacy of antibody-drug conjugates from preclinical",
    "to clinical using a semimechanistic PK/PD model: A case study with RC88.",
    "Clin Transl Sci. 2023;16(7):1232-1242. doi:10.1111/cts.13526"
  )
  vignette <- "Li_2023_RC88"
  units <- list(time = "h", dosing = "mg/kg", concentration = "ug/mL")

  covariateData <- list()

  compartmentData <- list(
    central = list(analyte = "RC88 (antibody-drug conjugate)", units = "mg/kg", specimen = "serum", verified = TRUE),
    peripheral1 = list(
      analyte = "RC88 (antibody-drug conjugate)",
      units = "mg/kg",
      specimen = "tissue",
      verified = TRUE
    ),
    cycling_cells = list(analyte = "proliferating tumor cells", units = "mm^3", specimen = "tumor", verified = TRUE),
    damaged_cells1 = list(
      analyte = "dying tumor cells (transit state 1)",
      units = "mm^3",
      specimen = "tumor",
      verified = TRUE
    ),
    damaged_cells2 = list(
      analyte = "dying tumor cells (transit state 2)",
      units = "mm^3",
      specimen = "tumor",
      verified = TRUE
    ),
    damaged_cells3 = list(
      analyte = "dying tumor cells (transit state 3)",
      units = "mm^3",
      specimen = "tumor",
      verified = TRUE
    )
  )

  population <- list(
    species = "mouse (Balb/c nude, H292-MSLN-9C8 xenograft)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    age_range = "4-5 weeks at purchase",
    sex_female_pct = 100,
    disease_state = "subcutaneous H292-MSLN-9C8 (mesothelin-transfected lung cancer) xenograft, staged at about 250 mm^3",
    dose_range = "vehicle, 0.75, 1.5 and 3 mg/kg RC88 IV q1w x 3",
    regions = "China (preclinical study, RemeGen)",
    notes = paste(
      "Table 1 lists 6 animals per group. Tumors were measured twice a week.",
      "The mouse PK parameters were fixed from Li_2023_RC88_mouse."
    )
  )

  ini({
    # Mouse PK, fixed from the PK fit (Li 2023 Table 2 'Mouse'; Results:
    # 'The mouse PK parameters derived above were fixed in the subsequent PD
    # modeling'). See Li_2023_RC88_mouse for the Vmax unit reading.
    lcl <- fixed(log(0.527)); label("Linear clearance CL_ADC (mL/h/kg)") # Table 2: CLADC = 0.527 mL/h/kg
    lq <- fixed(log(1.31)); label("Intercompartmental clearance CLD_ADC (mL/h/kg)") # Table 2: CLDADC = 1.31 mL/h/kg
    lvc <- fixed(log(67.1)); label("Central volume of distribution V1 (mL/kg)") # Table 2: V1 = 67.1 mL/kg
    lvp <- fixed(log(49.1)); label("Peripheral volume of distribution V2 (mL/kg)") # Table 2: V2 = 49.1 mL/kg
    lkm <- fixed(log(5.02)); label("Michaelis-Menten constant Km (ug/mL)") # Table 2: Km = 5.02 ug/mL
    lvmax <- fixed(log(1.73)); label("Maximum Michaelis-Menten elimination rate Vmax (ug/h/kg; Table 2 prints ug/mL/h)") # Table 2: Vmax = 1.73

    # Hybrid model, Table S4 column 'H292-MSLN-9C8'.
    ltumorLinGrowth <- log(7.17); label("Linear tumor growth rate Kgrow (mm^3/h)") # Table S4: Kgrow = 7.17 mm^3/h (-)
    ltumorExpGrowth <- log(0.005); label("Exponential tumor growth rate KgEx (1/h)") # Table S4: KgEx = 0.005 1/h (-)
    ltau <- log(146); label("Transit time of each dying-cell compartment t (h)") # Table S4: t = 146 h (92.7%)
    lkmax <- log(0.012); label("Maximum kill rate Kmax (1/h)") # Table S4: Kmax = 0.012 1/h (86.7%)
    lkc50 <- log(3.20); label("RC88 concentration giving half-maximal kill KC50 (ug/mL)") # Table S4: KC50 = 3.20 ug/mL (129%)
    ltsmax <- log(57365); label("Tumor carrying capacity Emax (mm^3)") # Table S4: Emax = 57365 mm^3 (-)
    psi <- fixed(20); label("Switch sharpness between exponential and linear growth psi (unitless)") # Table S4: psi = 20 (fixed)

    # The measured initial tumor volume of each animal was the initial
    # condition of E1 (Methods); treatment started at about 250 mm^3. Override
    # with the measured volume when available.
    lrbase_tumor <- fixed(log(250)); label("Tumor volume at the start of treatment (mm^3)") # Methods: tumors staged at ~250 mm^3

    # Residual and interindividual variability are not reported. Carried as
    # fixed(0) so the model simulates the typical-value trajectory.
    propSd_tumor_vol <- fixed(0); label("Proportional residual error on tumor volume (fraction; not reported in the source)")
  })
  model({
    cl <- exp(lcl)
    q <- exp(lq)
    vc <- exp(lvc)
    vp <- exp(lvp)
    km <- exp(lkm)
    vmax <- exp(lvmax)
    tumorLinGrowth <- exp(ltumorLinGrowth)
    tumorExpGrowth <- exp(ltumorExpGrowth)
    tau <- exp(ltau)
    kmax <- exp(lkmax)
    kc50 <- exp(lkc50)
    tsmax <- exp(ltsmax)
    rbase_tumor <- exp(lrbase_tumor)

    # Mouse PK: doses mg/kg, volumes mL/kg; x1000 gives ug/mL. The
    # Michaelis-Menten term is an amount rate in ug/h/kg (/1000 -> mg/kg/h).
    Cc <- 1000 * central / vc
    mm_elim <- vmax * Cc / (km + Cc) / 1000

    d/dt(central) <- -(cl / vc) * central - (q / vc) * central + (q / vp) * peripheral1 - mm_elim
    d/dt(peripheral1) <- (q / vc) * central - (q / vp) * peripheral1

    # Li 2023 Equation 2. The per-volume growth rate Kg falls from KgEx to
    # Kgrow/E once E passes Kgrow/KgEx, and to zero as E approaches Emax.
    tumor_vol <- cycling_cells + damaged_cells1 + damaged_cells2 + damaged_cells3
    kg <- tumorExpGrowth * (1 - tumor_vol / tsmax) / (1 + (tumorExpGrowth / tumorLinGrowth * tumor_vol)^psi)^(1 / psi)
    kkill <- kmax * Cc / (kc50 + Cc)

    d/dt(cycling_cells) <- (kg - kkill) * cycling_cells
    d/dt(damaged_cells1) <- kkill * cycling_cells - damaged_cells1 / tau
    d/dt(damaged_cells2) <- damaged_cells1 / tau - damaged_cells2 / tau
    d/dt(damaged_cells3) <- damaged_cells2 / tau - damaged_cells3 / tau

    cycling_cells(0) <- rbase_tumor

    tumor_vol ~ prop(propSd_tumor_vol)
  })
}
