Li_2023_RC88_jumbe_h292_mouse_tgi <- function() {
  description <- paste(
    "Preclinical (mouse, Balb/c nude bearing H292-MSLN-9C8 lung-cancer",
    "xenografts). Jumbe tumor-growth-inhibition model of",
    "the mesothelin-targeting antibody-drug conjugate RC88: surface-limited",
    "(volume^(2/3), spherical-tumor) growth and kill, a saturable (Kmax, KC50)",
    "kill term in serum RC88 concentration, and three transit compartments of",
    "dying cells, driven by the fixed mouse two-compartment PK model with",
    "parallel linear and Michaelis-Menten elimination (Li 2023)."
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

    # Jumbe model, Table S2 column 'H292-MSLN-9C8'.
    lkgrow <- log(0.042); label("Net tumor growth rate Kgrow (1/h as printed; multiplies volume^(2/3))") # Table S2: Kgrow = 0.042 1/h (4.23%)
    # Table S2 lists t as a 'Transit rate constant' with unit '(h)'. The
    # value enters Equation 3 wherever 1/t is printed: only then does the
    # paper's own TSC formula return Table 3's TSC (2.38 ug/mL here).
    ldamageTransit <- log(0.027); label("Transit rate constant of dying cells 1/t (1/h)") # Table S2: t = 0.027 (38.1%)
    lkmax <- log(21.6); label("Maximal tumor killing Kmax (unitless)") # Table S2: Kmax = 21.6 (-)
    lkc50 <- log(30.9); label("RC88 concentration giving half-maximal kill KC50 (ug/mL)") # Table S2: KC50 = 30.9 ug/mL (-)

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
    kgrow <- exp(lkgrow)
    damageTransit <- exp(ldamageTransit)
    kmax <- exp(lkmax)
    kc50 <- exp(lkc50)
    rbase_tumor <- exp(lrbase_tumor)

    # Mouse PK: doses mg/kg, volumes mL/kg; x1000 gives ug/mL. The
    # Michaelis-Menten term is an amount rate in ug/h/kg (/1000 -> mg/kg/h).
    Cc <- 1000 * central / vc
    mm_elim <- vmax * Cc / (km + Cc) / 1000

    d/dt(central) <- -(cl / vc) * central - (q / vc) * central + (q / vp) * peripheral1 - mm_elim
    d/dt(peripheral1) <- (q / vc) * central - (q / vp) * peripheral1

    # Li 2023 Equation 3: growth, kill and transit all act on volume^(2/3)
    # (the surface of a spherical tumor). A volume^(2/3) rate empties a state
    # in finite time, after which the solver can step it slightly below zero
    # where ^(2/3) is undefined; max(., 0) holds an emptied state at zero.
    tumor_vol <- cycling_cells + damaged_cells1 + damaged_cells2 + damaged_cells3
    kkill <- kmax * Cc / (kc50 + Cc)
    surf_cycling <- max(cycling_cells, 0)^(2 / 3)
    surf_damaged1 <- max(damaged_cells1, 0)^(2 / 3)
    surf_damaged2 <- max(damaged_cells2, 0)^(2 / 3)
    surf_damaged3 <- max(damaged_cells3, 0)^(2 / 3)

    d/dt(cycling_cells) <- (kgrow - damageTransit * kkill) * surf_cycling
    d/dt(damaged_cells1) <- kkill * damageTransit * surf_cycling - damageTransit * surf_damaged1
    d/dt(damaged_cells2) <- damageTransit * surf_damaged1 - damageTransit * surf_damaged2
    d/dt(damaged_cells3) <- damageTransit * surf_damaged2 - damageTransit * surf_damaged3

    cycling_cells(0) <- rbase_tumor

    tumor_vol ~ prop(propSd_tumor_vol)
  })
}
