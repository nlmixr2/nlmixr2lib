Oggianu_2023_trazodone_gabapentin_mouse <- function() {
  description <- paste(
    "Preclinical (mouse, male CD-1). Joint PK/PD model of the",
    "antinociceptive interaction of oral trazodone and intraperitoneal",
    "gabapentin in the acetic-acid writhing test. Each drug has",
    "zero-order input and two-compartment disposition parameterised by",
    "micro-constants. The intraperitoneal acetic-acid stimulus is a unit",
    "dose into a stimulus compartment scaled by the stimulus amplitude S;",
    "it passes through nine transit compartments into the observed",
    "writhing compartment, which empties first-order. Each drug removes",
    "stimulus in proportion to its plasma molar concentration, and a",
    "product term on the two concentrations carries the departure from",
    "additivity (synergy). Time is in minutes; the PK parameters keep",
    "their published hour units and are converted in model()."
  )
  reference <- paste(
    "Oggianu L, Ottaviani M, Maddalena F, Leoni M, Pellizzoni G,",
    "Gaudesi E, Piccinini L, Comerio M, Furlan R, Salvatore C, Tongiani S,",
    "Calisti F. PK/PD analysis of trazodone and gabapentin in neuropathic",
    "pain rodent models: Translational PK-PD modeling from nonclinical to",
    "clinical development. Clin Transl Sci. 2023;16(4):606-617.",
    "doi:10.1111/cts.13472. Mouse PK: Table 1; writhing PD: Table 2 and",
    "Table 4B; equations: Supplementary Table S1 and Figure S1."
  )
  vignette <- "Oggianu_2023_trazodone_gabapentin"
  units <- list(
    time = "min",
    dosing = "mg/kg (drugs); unit stimulus (acetic acid)",
    concentration = "ng/mL"
  )

  # The acetic-acid stimulus site and the writhing count are specific to
  # the writhing-test model; transit1-transit9 follow the registered
  # numbered-chain pattern.
  paper_specific_compartments <- c("stimulus", "writhing")

  covariateData <- list()

  compartmentData <- list(
    central = list(analyte = "trazodone", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "trazodone", units = "mg/kg", specimen = "tissue", verified = TRUE),
    central_gabapentin = list(analyte = "gabapentin", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1_gabapentin = list(analyte = "gabapentin", units = "mg/kg", specimen = "tissue", verified = TRUE),
    stimulus = list(
      analyte = "acetic acid stimulus",
      units = "writhing-count equivalents",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit1 = list(
      analyte = "acetic acid stimulus",
      units = "writhing-count equivalents",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit2 = list(
      analyte = "acetic acid stimulus",
      units = "writhing-count equivalents",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit3 = list(
      analyte = "acetic acid stimulus",
      units = "writhing-count equivalents",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit4 = list(
      analyte = "acetic acid stimulus",
      units = "writhing-count equivalents",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit5 = list(
      analyte = "acetic acid stimulus",
      units = "writhing-count equivalents",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit6 = list(
      analyte = "acetic acid stimulus",
      units = "writhing-count equivalents",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit7 = list(
      analyte = "acetic acid stimulus",
      units = "writhing-count equivalents",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit8 = list(
      analyte = "acetic acid stimulus",
      units = "writhing-count equivalents",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit9 = list(
      analyte = "acetic acid stimulus",
      units = "writhing-count equivalents",
      specimen = "not applicable",
      verified = TRUE
    ),
    writhing = list(analyte = "writhing count", units = "writhings", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "mouse (CD-1, male)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    weight_range = "20-25 g",
    sex_female_pct = 0,
    disease_state = "healthy fasted mice; acute visceral nociception induced by intraperitoneal 0.7% acetic acid (16 uL/g)",
    dose_range = paste(
      "PK (satellite groups, 12 mice per group, sparse sampling to 24 h):",
      "trazodone 1, 3 and 10 mg/kg oral; gabapentin 10, 30 and 100 mg/kg",
      "intraperitoneal; plus 1 + 3 and 10 + 30 mg/kg combinations. Writhing",
      "PD: trazodone 0.075-4.5 mg/kg oral, gabapentin 1.5-15 mg/kg",
      "intraperitoneal, and combinations 0.075 + 1.5, 0.15 + 1.5, 0.5 + 5",
      "and 1.5 + 15 mg/kg (trazodone + gabapentin), each given 1 h before",
      "the acetic acid."
    ),
    regions = "Italy (Angelini Pharma)",
    notes = paste(
      "Writhing groups were n = 6-9 mice (Figure S4 legend); the total",
      "number of animals is not reported. Writhings were counted in 5-min",
      "bins for 30 min after the acetic acid. Analysis was sequential in",
      "NONMEM VI (FOCE with interaction): PK and the control-group PD",
      "parameters were estimated first and fixed, then the drug potencies",
      "and the synergy term. The packaged parameter set is the synergic",
      "combination model of Table 4B, the model the authors carried",
      "forward to the human translation."
    )
  )

  ini({
    # --- Trazodone PK in mice (Oggianu 2023 Table 1a; oral) -------------
    lkel <- log(1.25); label("Trazodone elimination rate constant k (1/h)") # Table 1a k = 1.25 1/h
    lvc <- log(5.95); label("Trazodone apparent central volume V1/F (L/kg)") # Table 1a V1/F = 5.95 L/kg
    lk12 <- log(0.982); label("Trazodone central-to-peripheral rate constant k12 (1/h)") # Table 1a k12 = 0.982 1/h
    lk21 <- log(1.51); label("Trazodone peripheral-to-central rate constant k21 (1/h)") # Table 1a k21 = 1.51 1/h
    ld1 <- fixed(log(0.0833)); label("Trazodone zero-order absorption duration D1 (h)") # Table 1a D1 = 0.0833 h, no SE reported

    # --- Gabapentin PK in mice (Oggianu 2023 Table 1b; intraperitoneal) --
    lkel_gabapentin <- log(0.112); label("Gabapentin elimination rate constant k (1/h)") # Table 1b k = 0.112 1/h
    lvc_gabapentin <- log(0.846); label("Gabapentin apparent central volume V1/F (L/kg)") # Table 1b V1/F = 0.846 L/kg
    lk12_gabapentin <- log(0.935); label("Gabapentin central-to-peripheral rate constant k12 (1/h)") # Table 1b k12 = 0.935 1/h
    lk21_gabapentin <- log(0.110); label("Gabapentin peripheral-to-central rate constant k21 (1/h)") # Table 1b k21 = 0.110 1/h
    ld1_gabapentin <- fixed(log(0.167)); label("Gabapentin zero-order absorption duration D1 (h)") # Table 1b D1 = 0.167 h, no SE reported

    # --- Writhing PD, control-group parameters (Oggianu 2023 Table 2) ----
    # Estimated on the vehicle groups; the drug PK/PD fits then held these,
    # and the PK above, at their estimates (sequential analysis).
    lktr <- log(1.36); label("Stimulus transit rate constant ktr (1/min)") # Table 2 ktr = 1.36 1/min
    lstim <- log(26.10); label("Acetic-acid stimulus amplitude S at injection (writhings)") # Table 2 S = 26.10
    lkout <- log(0.104); label("First-order loss rate constant k of the writhing compartment (1/min)") # Table 2 k = 0.104 1/min

    # --- Drug potencies and synergy (Oggianu 2023 Table 4B) --------------
    lslope <- log(10.9); label("Trazodone potency k2A on stimulus loss (L/umol/min)") # Table 4B trazodone k2 = 10.9 L/umol/min
    lslope_gabapentin <- log(0.0490); label("Gabapentin potency k2B on stimulus loss (L/umol/min)") # Table 4B gabapentin k2 = 0.0490 L/umol/min
    synergy <- 16.0; label("Combination index gamma on the product of the two molar concentrations (L^2/umol^2/min)") # Table 4B gamma = 16.0

    # --- IIV (variances; %CV in the tables is sqrt(omega^2)) -------------
    etalvc ~ 0.0526 # Table 1a omega2 V = 0.0526 (22.9 %CV)
    etalkel ~ 0.00899 # Table 1a omega2 k = 0.00899 (9.48 %CV)
    etalk12 ~ 0.577 # Table 1a omega2 k12 = 0.577 (76.0 %CV)
    etalvc_gabapentin ~ 0.00287 # Table 1b omega2 V = 0.00287 (5.36 %CV)
    etalkel_gabapentin ~ 0.337 # Table 1b omega2 k = 0.337 (58.1 %CV)
    etalk12_gabapentin ~ 1.43 # Table 1b omega2 k12 = 1.43 (120 %CV)
    etalstim ~ 0.0231 # Table 2 omega2 S = 0.0231 (15.2 %CV)
    etalkout ~ 0.0349 # Table 2 omega2 k = 0.0349 (18.7 %CV)
    etalslope ~ 1.42 # Table 4B trazodone omega2 k2 = 1.42 (119.2 %CV)
    etalslope_gabapentin ~ 1.78 # Table 4B gabapentin omega2 k2 = 1.78 (133.4 %CV)

    # --- Residual error (variances converted to SDs) ---------------------
    propSd <- sqrt(0.0752); label("Proportional residual error, trazodone plasma concentration (fraction)") # Table 1a sigma2 prop = 0.0752 (27.4 %CV)
    propSd_gabapentin <- sqrt(0.0598); label("Proportional residual error, gabapentin plasma concentration (fraction)") # Table 1b sigma2 prop = 0.0598 (24.5 %CV)
    addSd_writhing <- sqrt(2.64); label("Additive residual error, writhing count (writhings)") # Table 4B sigma2 add = 2.64
  })

  model({
    # Published PK rate constants and durations are in hours; the PD and
    # the model clock are in minutes.
    minPerHour <- 60
    # Molecular weights (g/mol) of the free bases, the analytes measured in
    # plasma. Concentration in ng/mL divided by MW gives umol/L, the unit
    # the potencies k2 are reported against.
    mwTrazodone <- 371.86
    mwGabapentin <- 171.24

    # --- Trazodone PK -----------------------------------------------------
    kel <- exp(lkel + etalkel) / minPerHour
    vc <- exp(lvc + etalvc)
    k12 <- exp(lk12 + etalk12) / minPerHour
    k21 <- exp(lk21) / minPerHour
    d1 <- exp(ld1) * minPerHour

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    dur(central) <- d1

    # --- Gabapentin PK ----------------------------------------------------
    kel_gabapentin <- exp(lkel_gabapentin + etalkel_gabapentin) / minPerHour
    vc_gabapentin <- exp(lvc_gabapentin + etalvc_gabapentin)
    k12_gabapentin <- exp(lk12_gabapentin + etalk12_gabapentin) / minPerHour
    k21_gabapentin <- exp(lk21_gabapentin) / minPerHour
    d1_gabapentin <- exp(ld1_gabapentin) * minPerHour

    d/dt(central_gabapentin) <- -kel_gabapentin * central_gabapentin -
      k12_gabapentin * central_gabapentin + k21_gabapentin * peripheral1_gabapentin
    d/dt(peripheral1_gabapentin) <- k12_gabapentin * central_gabapentin -
      k21_gabapentin * peripheral1_gabapentin
    dur(central_gabapentin) <- d1_gabapentin

    # Doses are mg/kg and volumes L/kg, so amount/volume is mg/L; x1000
    # reports ng/mL.
    Cc <- 1000 * central / vc
    Cc_gabapentin <- 1000 * central_gabapentin / vc_gabapentin
    cumolTrazodone <- Cc / mwTrazodone
    cumolGabapentin <- Cc_gabapentin / mwGabapentin

    # --- Writhing PD (Supplementary Table S1B; Figure S1D) ----------------
    ktr <- exp(lktr)
    stim <- exp(lstim + etalstim)
    kout <- exp(lkout + etalkout)
    slope <- exp(lslope + etalslope)
    slope_gabapentin <- exp(lslope_gabapentin + etalslope_gabapentin)

    # Drug-driven first-order loss of the stimulus: k2A*CA + k2B*CB +
    # gamma*CA*CB (Figure S1D); gamma = 0 gives the additive model.
    kdrug <- slope * cumolTrazodone + slope_gabapentin * cumolGabapentin +
      synergy * cumolTrazodone * cumolGabapentin

    d/dt(stimulus) <- -kdrug * stimulus - ktr * stimulus
    d/dt(transit1) <- ktr * stimulus - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(transit4) <- ktr * transit3 - ktr * transit4
    d/dt(transit5) <- ktr * transit4 - ktr * transit5
    d/dt(transit6) <- ktr * transit5 - ktr * transit6
    d/dt(transit7) <- ktr * transit6 - ktr * transit7
    d/dt(transit8) <- ktr * transit7 - ktr * transit8
    d/dt(transit9) <- ktr * transit8 - ktr * transit9
    d/dt(writhing) <- ktr * transit9 - kout * writhing

    # The acetic acid is dosed as amt = 1 into `stimulus`; the bioavailability
    # turns that unit dose into the subject's stimulus amplitude S, which is
    # the paper's initial condition a0(0) = S.
    f(stimulus) <- stim

    Cc ~ prop(propSd)
    Cc_gabapentin ~ prop(propSd_gabapentin)
    writhing ~ add(addSd_writhing)
  })
}
