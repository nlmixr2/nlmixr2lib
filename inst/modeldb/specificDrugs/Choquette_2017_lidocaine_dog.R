Choquette_2017_lidocaine_dog <- function() {
  description <- "Preclinical (beagle dog). One-compartment population PK model for lidocaine hydrochloride injected perineurally (ultrasound/neurostimulation-guided paravertebral brachial plexus block, 6 mg/kg split over the C6, C7 and C8-T1 roots) with or without adrenaline (5 ug/mL after 1:1 dilution) in isoflurane-anaesthetised beagles (Choquette 2017, Model C). Absorption is two parallel, time-constrained zero-order inputs into the central compartment: a fraction f of the dose enters over Dur1 from the dose time, and the remaining 1 - f enters over Dur2 after a delay tau. CL/F and V/F (per kg) are shared by both formulations; f, Dur1, Dur2 and tau are estimated separately for the plain and the adrenalinated formulation and are selected by CONMED_EPI."
  reference <- "Choquette A, Troncy E, Guillot M, Varin F, del Castillo JRE. Pharmacokinetics of Lidocaine Hydrochloride Administered with or without Adrenaline for the Paravertebral Brachial Plexus Block in Dogs. PLoS One. 2017;12(1):e0169745. doi:10.1371/journal.pone.0169745"
  vignette <- "Choquette_2017_lidocaine_dog"
  units <- list(time = "min", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "CL/F and V/F are reported per kg (Table 1: L/min/kg and L/kg) and are scaled linearly by WT; the dose is in mg (6 mg/kg in the study). Cohort mean 11 (SD 1.3) kg.",
      source_name = "body weight (Methods, Animals)"
    ),
    CONMED_EPI = list(
      description = "Adrenaline (epinephrine) co-formulated with the lidocaine injection",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "1 = adrenalinated lidocaine (20 mg/mL lidocaine HCl + 0.01 mg/mL adrenaline, diluted 1:1 with saline, so 10 mg/mL + 5 ug/mL as injected); 0 = plain lidocaine (same dilution, no adrenaline). Per-administration (crossover) indicator, not per-subject. It selects which formulation's set of absorption parameters (f, Dur1, Dur2, tau and their IIV) applies; it is not a multiplicative covariate effect.",
      source_name = "Formulation (supplementary data S1: 'Plain' / 'Adrenalinated'); subscript j = LI / LA in Fig 1"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "lidocaine", units = "mg", specimen = "administration site", verified = TRUE),
    depot2 = list(analyte = "lidocaine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "lidocaine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "beagle dog",
    n_subjects = 8L,
    n_studies = 1L,
    age_range = "2-5 years",
    weight_range = "mean 11 (SD 1.3) kg",
    sex_female_pct = 100,
    disease_state = "healthy; isoflurane general anaesthesia (medetomidine/butorphanol premedication, propofol induction) during block execution, reversed with atipamezole after the third injection",
    dose_range = "lidocaine hydrochloride 6 mg/kg total (2 mg/kg per nerve root, three injections within 4 min) as plain or adrenalinated solution",
    regions = "Canada (Universite de Montreal)",
    notes = "Randomised, blinded crossover over three periods: 19 profiles (10 plain, 9 adrenalinated) from 8 female beagles; 13 plasma samples per profile to 180 min. Population analysis by MLEM in ADAPT 5 (3000 importance samples), log-normal parameters, diagonal covariance (Methods, Pharmacokinetic analysis; Results)."
  )

  ini({
    # Table 1 reports every Model parameter as 'log-normal mean (log-normal SD)':
    # the arithmetic mean m and SD s of a log-normal distribution. The
    # typical value (median) is m / sqrt(1 + (s/m)^2) and the log-scale
    # variance is log(1 + (s/m)^2).

    # Disposition, shared by both formulations (per kg; scaled by WT in model())
    lcl <- log(0.017 / sqrt(1 + (0.005 / 0.017)^2)); label("Apparent clearance CL/F per kg (L/min/kg)") # Table 1: CL/F = 0.017 (0.005) L/min/kg, both formulations
    lvc <- log(1.850 / sqrt(1 + (0.580 / 1.850)^2)); label("Apparent volume V/F per kg (L/kg)") # Table 1: V/F = 1.850 (0.580) L/kg, both formulations

    # Absorption, plain lidocaine (j = LI)
    lfdepot_plain <- log(0.56 / sqrt(1 + (0.22 / 0.56)^2)); label("Dose fraction entering the early zero-order input, plain (fraction)") # Table 1: f = 0.56 (0.22), plain
    ld1_plain <- log(5.77 / sqrt(1 + (3.24 / 5.77)^2)); label("Duration of the early zero-order input Dur1, plain (min)") # Table 1: Dur1 = 5.77 (3.24) min, plain
    ld2_plain <- log(18.10 / sqrt(1 + (4.87 / 18.10)^2)); label("Duration of the late zero-order input Dur2, plain (min)") # Table 1: Dur2 = 18.10 (4.87) min, plain
    ltlag_plain <- log(6.52 / sqrt(1 + (0.02 / 6.52)^2)); label("Delay tau to the onset of the late zero-order input, plain (min)") # Table 1: tau = 6.52 (0.02) min, plain

    # Absorption, adrenalinated lidocaine (j = LA)
    lfdepot_adr <- log(0.54 / sqrt(1 + (0.15 / 0.54)^2)); label("Dose fraction entering the early zero-order input, adrenalinated (fraction)") # Table 1: f = 0.54 (0.15), adrenalinated
    ld1_adr <- log(11.30 / sqrt(1 + (7.05 / 11.30)^2)); label("Duration of the early zero-order input Dur1, adrenalinated (min)") # Table 1: Dur1 = 11.30 (7.05) min, adrenalinated
    ld2_adr <- log(254.00 / sqrt(1 + (129.00 / 254.00)^2)); label("Duration of the late zero-order input Dur2, adrenalinated (min)") # Table 1: Dur2 = 254.00 (129.00) min, adrenalinated
    ltlag_adr <- log(3.87 / sqrt(1 + (0.03 / 3.87)^2)); label("Delay tau to the onset of the late zero-order input, adrenalinated (min)") # Table 1: tau = 3.87 (0.03) min, adrenalinated

    # IIV: diagonal (Methods: 'a diagonal covariance matrix of estimable parameters')
    etalcl ~ log(1 + (0.005 / 0.017)^2) # Table 1: CL/F SD 0.005 on mean 0.017
    etalvc ~ log(1 + (0.580 / 1.850)^2) # Table 1: V/F SD 0.580 on mean 1.850
    etalfdepot_plain ~ log(1 + (0.22 / 0.56)^2) # Table 1: f SD 0.22 on mean 0.56, plain
    etald1_plain ~ log(1 + (3.24 / 5.77)^2) # Table 1: Dur1 SD 3.24 on mean 5.77, plain
    etald2_plain ~ log(1 + (4.87 / 18.10)^2) # Table 1: Dur2 SD 4.87 on mean 18.10, plain
    etaltlag_plain ~ log(1 + (0.02 / 6.52)^2) # Table 1: tau SD 0.02 on mean 6.52, plain
    etalfdepot_adr ~ log(1 + (0.15 / 0.54)^2) # Table 1: f SD 0.15 on mean 0.54, adrenalinated
    etald1_adr ~ log(1 + (7.05 / 11.30)^2) # Table 1: Dur1 SD 7.05 on mean 11.30, adrenalinated
    etald2_adr ~ log(1 + (129.00 / 254.00)^2) # Table 1: Dur2 SD 129.00 on mean 254.00, adrenalinated
    etaltlag_adr ~ log(1 + (0.03 / 3.87)^2) # Table 1: tau SD 0.03 on mean 3.87, adrenalinated

    # Residual error: the ADAPT 5 variance model is not reported anywhere in the source.
    propSd <- fixed(0); label("Proportional residual error (fraction; ZERO - not reported in source)") # not reported in the paper or its supplement
  })
  model({
    # Individual disposition parameters (per-kg values scaled linearly by WT)
    cl <- exp(lcl + etalcl) * WT
    vc <- exp(lvc + etalvc) * WT
    kel <- cl / vc

    # Formulation-specific absorption parameters; CONMED_EPI selects the set
    fdepot_plain <- exp(lfdepot_plain + etalfdepot_plain)
    d1_plain <- exp(ld1_plain + etald1_plain)
    d2_plain <- exp(ld2_plain + etald2_plain)
    tlag_plain <- exp(ltlag_plain + etaltlag_plain)
    fdepot_adr <- exp(lfdepot_adr + etalfdepot_adr)
    d1_adr <- exp(ld1_adr + etald1_adr)
    d2_adr <- exp(ld2_adr + etald2_adr)
    tlag_adr <- exp(ltlag_adr + etaltlag_adr)
    fdepot <- (1 - CONMED_EPI) * fdepot_plain + CONMED_EPI * fdepot_adr
    d1 <- (1 - CONMED_EPI) * d1_plain + CONMED_EPI * d1_adr
    d2 <- (1 - CONMED_EPI) * d2_plain + CONMED_EPI * d2_adr
    tlag <- (1 - CONMED_EPI) * tlag_plain + CONMED_EPI * tlag_adr
    # f is log-normal, so an individual draw can exceed 1 (about 4% of plain
    # and 1% of adrenalinated draws), which would put a negative amount in
    # depot2. Cap it at 1 so that 1 - f stays a valid fraction.
    if (fdepot > 1) fdepot <- 1

    # Zero-order input rates (Fig 1C): k01 = f * Dose / Dur1 out of depot over
    # [0, Dur1] after its dose, k02 = (1 - f) * Dose / Dur2 out of depot2 over
    # [tau, tau + Dur2]. podo() is the administered amount before f() is
    # applied; tad(depot2) counts from the lagged dose time. tad() is NA before
    # the first dose, hence the guards. A new dose resets podo()/tad(), so the
    # rates are exact for single administrations or for administrations whose
    # inputs do not overlap.
    k01 <- 0
    k02 <- 0
    if (!is.na(tad(depot))) {
      if (tad(depot) < d1) k01 <- fdepot * podo(depot) / d1
    }
    if (!is.na(tad(depot2))) {
      if (tad(depot2) < d2) k02 <- (1 - fdepot) * podo(depot2) / d2
    }

    d/dt(depot) <- -k01
    d/dt(depot2) <- -k02
    d/dt(central) <- k01 + k02 - kel * central

    # Each administration is two dose records with the same amt (one into
    # depot, one into depot2); f() splits the dose and alag() delays the late
    # input by tau.
    f(depot) <- fdepot
    f(depot2) <- 1 - fdepot
    alag(depot2) <- tlag

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
