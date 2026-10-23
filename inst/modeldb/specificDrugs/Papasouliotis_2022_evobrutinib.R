Papasouliotis_2022_evobrutinib <- function() {
  description <- paste(
    "Two-compartment population PK model with sequential zero-order (into the depot) then",
    "first-order absorption, absorption lag time and first-order elimination for oral evobrutinib",
    "tablets in healthy adults, with a food effect on bioavailability (+49%) and on the zero-order",
    "input duration (+427%), linked to an irreversible-binding turnover model of the fraction of",
    "unoccupied Bruton's tyrosine kinase (BTK) in peripheral blood mononuclear cells"
  )
  reference <- paste(
    "Papasouliotis O, Mitchell D, Girard P, Dyroff M. Population pharmacokinetic and",
    "pharmacodynamic modeling of evobrutinib in healthy adult participants. Clin Transl Sci.",
    "2022;15(12):2899-2908. doi:10.1111/cts.13417"
  )
  vignette <- "Papasouliotis_2022_evobrutinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    FED = list(
      description = "Prandial state at dosing: 1 = dosed with food (low-fat meal), 0 = dosed fasted",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted)",
      notes = paste(
        "Dose-record (time-varying) indicator: the food-effect studies were crossovers, so a",
        "participant contributes both levels. In the source studies the fed condition was a",
        "low-fat meal (evobrutinib given 30 min after starting the meal in Study MS200527_0017",
        "Part B); the paper labels the covariate generically as 'food state (with food or while",
        "fasted)' and 'Food on F1' / 'Food on D1', so the general FED indicator is used rather",
        "than FED_LOWFAT. Acts on relative bioavailability (x 1.49 when fed), on the zero-order",
        "input duration D1 (x (1 + 4.27) when fed) and selects which of the two D1 etas applies."
      ),
      source_name = "food state"
    )
  )

  covariatesDataExcluded <- list(
    RACE_JAPANESE = list(
      description = "Japanese versus non-Japanese ethnicity",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Tested in the stepwise covariate search on the PK and BTK occupancy parameters and not",
        "retained (Results: 'Ethnicity (Japanese or non-Japanese) was not a significant",
        "covariate')."
      )
    ),
    WT = list(
      description = "Body weight (also body mass index and body surface area)",
      units = "kg",
      type = "continuous",
      notes = "Tested in the stepwise covariate search and not retained."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested in the stepwise covariate search and not retained."
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Tested in the stepwise covariate search and not retained."
    ),
    DOSE = list(
      description = "Evobrutinib dose",
      units = "mg",
      type = "continuous",
      notes = "Tested in the stepwise covariate search and not retained (PK linear over 25-200 mg)."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "evobrutinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "evobrutinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "evobrutinib", units = "mg", specimen = "plasma", verified = TRUE),
    target = list(
      analyte = "Bruton's tyrosine kinase, free (unoccupied) fraction relative to predose",
      units = "(fraction)",
      specimen = "blood cell",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 76,
    n_studies = 2,
    age_range = "adults (range not reported)",
    disease_state = "healthy adult volunteers (Japanese and non-Japanese)",
    dose_range = paste(
      "25, 75 or 200 mg once daily for 6 days fasted (Study MS200527_0017 Part A, n = 42);",
      "single 75 mg dose fed (low-fat meal) and fasted in crossover (Study MS200527_0017 Part B,",
      "n = 16; Study MS200527_0019, n = 18)"
    ),
    regions = "not reported",
    race_ethnicity = c(
      Japanese = "29 of the 58 Study MS200527_0017 participants (7 per Part A cohort, 8 in Part B); not reported for Study MS200527_0019"
    ),
    n_observations = "2326 PK observations (675, 29%, below the 0.600 ng/mL LLOQ, handled with M3); 441 BTK occupancy observations from 41 Part A participants",
    formulation = "tablet (25 mg tablets)",
    notes = paste(
      "Papasouliotis 2022 Abstract, Methods 'Study population' and Table 1. The BTK occupancy",
      "model was fitted to Part A only (14, 13 and 14 participants at 25, 75 and 200 mg q.d.).",
      "Baseline demographics (age, weight, sex distribution) are not tabulated in the paper."
    )
  )

  ini({
    # PK structural parameters, fasted tablet reference (Papasouliotis 2022 Table 2)
    lcl <- log(273); label("Apparent clearance CL/F (L/h)") # Table 2 'CL/F (L/h)' = 273 (RSE 3.7%)
    lvc <- log(61.1); label("Apparent central volume V2/F (L)") # Table 2 'V2/F (L)' = 61.1 (RSE 13.8%)
    lq <- log(37); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 2 'Q/F (L/h)' = 37 (RSE 13.3%)
    lvp <- log(446); label("Apparent peripheral volume V3/F (L)") # Table 2 'V3/F (L)' = 446 (RSE 21.7%)
    lka <- log(0.784); label("First-order absorption rate constant ka (1/h)") # Table 2 'Ka (/h)' = 0.784 (RSE 2.5%)
    ltlag <- log(0.235); label("Absorption lag time ALAG (h)") # Table 2 'ALAG (h)' = 0.235 (RSE 4.4%)
    ld1 <- log(0.181); label("Duration of zero-order input into the depot D1, fasted (h)") # Table 2 'Zero-order absorption D1 (h)' = 0.181 (RSE 20.8%)

    # Food effects (FED = 1 with a low-fat meal)
    e_fed_fdepot <- 1.49; label("Multiplicative food effect on relative bioavailability F1 (unitless)") # Table 2 'Food on F1' = 1.49 (RSE 3.8%); Results 'increased by 49% when given with food'
    e_fed_d1 <- 4.27; label("Fractional food effect on D1, D1 x (1 + e_fed_d1) when fed (unitless)") # Table 2 'Food on D1' = 4.27 (RSE 20.4%); footnote 'D1 = 0.181 *(1 + 4.27) = 0.954 h under fed condition'

    # PK between-participant variability (log-normal; CV% = sqrt(exp(omega^2) - 1))
    etalcl + etaltlag ~ c(0.0958, 0.0558, 0.156) # Table 2 'BPV CL/F' 0.0958 (CV 31.7%), 'Cov CL/F + ALAG' 0.0558, 'BPV ALAG' 0.156 (CV 41.1%)
    etalq + etalvp ~ c(1.17, 0.988, 1.86) # Table 2 'BPV Q/F' 1.17 (CV 149%), 'Cov Q/F + V3/F' 0.988, 'BPV V3/F' 1.86 (CV 233%)
    etalvc ~ 0.627 # Table 2 'BPV V2/F' 0.627 (CV 93.3%)
    etald1_fasted ~ 1.39 # Table 2 'BPV D1 Fasted' 1.39 (CV 174%); separate eta per prandial state, same variance
    etald1_fed ~ 1.39 # Table 2 'BPV D1 Fed' 1.39 (CV 174%)

    # PK residual error
    addSd <- 0.214; label("Additive residual error on evobrutinib plasma concentration (ng/mL)") # Table 2 'Additive error' = 0.214 (RSE 20.8%)
    propSd <- 0.49; label("Proportional residual error on evobrutinib plasma concentration (fraction)") # Table 2 'Proportional error' = 0.49 (RSE 3.9%); Results 'proportional residual error estimated to be 49%'

    # BTK occupancy (irreversible binding) model (Papasouliotis 2022 Table 3)
    lkout <- log(0.00437); label("First-order BTK protein elimination rate constant kout (1/h)") # Table 3 'kout (L/h)' = 0.00437 (RSE 7.7%); unit misprinted, recovery half-life ln(2)/0.00437 = 159 h as stated in Results
    lkirrev <- log(0.0135); label("Second-order irreversible BTK binding rate constant kirrev (mL/ng/h)") # Table 3 'kirrev (ml/ng/h)' = 0.0135 (RSE 8.1%)

    etalkout ~ 0.217 # Table 3 'BPV kout' 0.217 (CV 49.2%)
    etalkirrev ~ 0.103 # Table 3 'BPV kirrev' 0.103 (CV 32.9%)

    addSd_target <- 0.00317; label("Additive residual error on the fraction of unoccupied BTK (fraction)") # Table 3 'Additive error' = 0.00317 (RSE 16.7%)
    propSd_target <- 0.43; label("Proportional residual error on the fraction of unoccupied BTK (fraction)") # Table 3 'Proportional error' = 0.43 (RSE 8.4%)
  })

  model({
    # Individual PK parameters
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)
    ka <- exp(lka)
    tlag <- exp(ltlag + etaltlag)

    # Zero-order input duration: one eta per prandial state (Table 2 'BPV D1
    # Fasted' / 'BPV D1 Fed'), and the fed duration scaled by (1 + 4.27)
    # per the Table 2 footnote. FED is strictly 0 or 1.
    d1_fasted <- exp(ld1 + etald1_fasted)
    d1_fed <- exp(ld1 + etald1_fed) * (1 + e_fed_d1)
    d1 <- (1 - FED) * d1_fasted + FED * d1_fed

    # Relative bioavailability: fasted reference F1 = 1, x 1.49 when fed
    fdepot <- e_fed_fdepot^FED

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # BTK turnover: kin = kout * BTK_FU,0 with BTK_FU,0 = 1 (baseline-corrected data)
    kout <- exp(lkout + etalkout)
    kirrev <- exp(lkirrev + etalkirrev)
    kin <- kout

    # Plasma concentration in ng/mL (dose mg / volume L = mg/L = 1000 ng/mL);
    # kirrev (mL/ng/h) x Cc (ng/mL) gives a first-order rate in 1/h.
    Cc <- 1000 * central / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(target) <- kin - kout * target - kirrev * Cc * target

    target(0) <- 1
    f(depot) <- fdepot
    alag(depot) <- tlag
    dur(depot) <- d1

    # BTK occupancy (fraction occupied) for reporting
    btko <- 1 - target

    Cc ~ add(addSd) + prop(propSd)
    target ~ add(addSd_target) + prop(propSd_target)
  })
}
