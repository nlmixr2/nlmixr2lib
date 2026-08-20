Khwarg_2024_donepezil_im <- function() {
  description <- "Two-compartment population PK model for long-acting intramuscular donepezil (GB-5001) with three-phase absorption: two lagged parallel first-order depots plus a simultaneous zero-order input, in healthy adult men"
  reference <- paste(
    "Khwarg J, Lee H, Yu KS, Seol E, Chung JY.",
    "Population Pharmacokinetic Modeling and Simulation for Dose Optimization",
    "of GB-5001, a Long-Acting Intramuscular Injection of Donepezil,",
    "in Healthy Participants. Neurol Ther. 2024;13(5):1453-1466.",
    "doi:10.1007/s40120-024-00643-4.",
    "Companion oral model: modellib('Khwarg_2024_donepezil_oral').",
    sep = " "
  )
  vignette <- "Khwarg_2024_donepezil"
  # Each IM administration is entered as three parallel dose records (Fig. 2:
  # depot 4, depot 5, and the zero-order arm into central). Declared explicitly
  # because the registry's default detection only recognises depot and central.
  dosing <- c("depot", "depot2", "central")
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  compartmentData <- list(
    depot = list(analyte = "donepezil", units = "mg", specimen = "administration site", verified = TRUE),
    depot2 = list(analyte = "donepezil", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "donepezil", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "donepezil", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 27,
    n_studies = 1,
    age_range = "18-55 years (eligibility); cohort mean 39.9 (SD 9.8) years",
    weight_range = "cohort mean 78.7 (SD 7.3) kg",
    sex_female_pct = 0,
    race_ethnicity = c(White = 94.4, Other = 5.6),
    disease_state = "healthy",
    dose_range = "70, 140 and 280 mg single intramuscular (gluteal) dose of GB-5001",
    regions = "USA (Advarra IRB, Columbia MD); sponsor G2GBIO, Republic of Korea",
    notes = paste(
      "Healthy non-smoking men, BMI 18.5-30.0 kg/m2, from the intramuscular arms",
      "(cohorts 1-3) of the NCT05525780 phase 1 dose-escalation study. 36 participants",
      "were enrolled across the three IM cohorts and randomised 3:1 to GB-5001 or",
      "placebo, so 27 contributed donepezil concentrations (9 per dose level, Table 1).",
      "The quoted demographics are for all 36 IM-cohort participants including placebo",
      "(mean BMI 25.3 (SD 2.6) kg/m2, height 176.6 (SD 5.3) cm).",
      "The oral and intramuscular arms were fitted in a single NONMEM run but share",
      "no parameters, so they are distributed as two model files.",
      sep = " "
    )
  )

  ini({
    # Absorption -- Khwarg 2024 Table 2, 'Intramuscular' block. The IM dose is split
    # three ways (Fig. 2): fraction F4 into a lagged first-order depot, fraction F5
    # into a second lagged first-order depot, and the remaining 1 - F4 - F5 as a
    # zero-order input delivered directly into the central compartment over D2.
    lka <- log(0.00402); label("Absorption rate constant from the first IM depot (KA4, 1/h)") # Table 2 IM: KA4 = 0.00402 1/h (RSE 6.7%; bootstrap 0.004, 95% CI 0.004-0.005)
    lka2 <- log(0.0134); label("Absorption rate constant from the second IM depot (KA5, 1/h)") # Table 2 IM: KA5 = 0.0134 1/h (RSE 17.2%; bootstrap 0.015, 95% CI 0.011-0.034)
    lfdepot <- log(0.748); label("Fraction of the IM dose entering the first depot (F4, fraction)") # Table 2 IM: F4 = 0.748 (RSE 1.9%; bootstrap 0.748, 95% CI 0.719-0.778)
    lfdepot2 <- log(0.145); label("Fraction of the IM dose entering the second depot (F5, fraction)") # Table 2 IM: F5 = 0.145 (RSE 11.2%; bootstrap 0.144, 95% CI 0.106-0.172)
    ltlag <- log(235); label("Absorption lag time of the first IM depot (ALAG4, h)") # Table 2 IM: ALAG4 = 235 h (RSE 0.5%; bootstrap 234.934, 95% CI 231.951-236.813)
    ltlag2 <- log(645); label("Absorption lag time of the second IM depot (ALAG5, h)") # Table 2 IM: ALAG5 = 645 h (RSE 0.2%; bootstrap 644.194, 95% CI 637.819-709.858)
    ld1 <- log(648); label("Duration of the zero-order absorption input into central (D2, h)") # Table 2 IM: D2 = 648 h (RSE 0.2%; bootstrap 647.833, 95% CI 574.904-820.454)

    # Disposition -- Khwarg 2024 Table 2, 'Intramuscular' block. Estimated separately
    # from the oral arm because of the different distribution and elimination
    # profiles of the modified-release IM formulation (Khwarg 2024 Results,
    # 'Population Pharmacokinetic Analysis').
    lcl <- log(10.3); label("Clearance (CL, L/h)") # Table 2 IM: CL = 10.3 L/h (RSE 6.1%; bootstrap 10.295, 95% CI 9.112-11.595)
    lvc <- log(503); label("Central volume of distribution (V2, L)") # Table 2 IM: V2 = 503 L (RSE 29%; bootstrap 511.553, 95% CI 259.082-834.504)
    lq <- log(185); label("Inter-compartmental clearance (Q, L/h)") # Table 2 IM: Q = 185 L/h (RSE 18.1%; bootstrap 183.88, 95% CI 128.741-250.004)
    lvp <- log(1160); label("Peripheral volume of distribution (V3, L)") # Table 2 IM: V3 = 1160 L (RSE 8.3%; bootstrap 1176.22, 95% CI 1026.616-1357.204)

    # IIV -- Khwarg 2024 Table 2 reports the log-scale variance directly; the
    # parenthesised percentage is the derived CV via footnote b,
    # CV = sqrt(exp(omega^2) - 1) * 100. Values below are the variances.
    # sqrt(exp(0.0936) - 1) = 31.3%; sqrt(exp(1.13) - 1) = 144.8%;
    # correlation 0.113 / sqrt(0.0936 * 1.13) = 0.347.
    etalcl + etalvc ~ c(0.0936, 0.113, 1.13) # Table 2 IM: IIV CL 0.0936 (31.3%), covariance CL vs V2 0.113, IIV V2 1.13 (144.8%)

    # Residual error -- reported on the SD scale: the additive term carries linear
    # concentration units (pg/mL, not pg^2/mL^2), so both terms are standard
    # deviations rather than NONMEM $SIGMA variances.
    propSd <- 0.223; label("Proportional residual error (fraction)") # Table 2 IM: proportional error = 0.223 (RSE 5%; bootstrap 0.219, 95% CI 0.196-0.241)
    addSd <- 0.0235; label("Additive residual error (ug/L)") # Table 2 IM: additive error = 23.5 pg/mL = 0.0235 ug/L (RSE 14.4%; bootstrap 23.353, 95% CI 15.668-30.726)
  })

  model({
    # Individual parameters
    ka <- exp(lka)
    ka2 <- exp(lka2)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)
    fdepot <- exp(lfdepot)
    fdepot2 <- exp(lfdepot2)
    tlag <- exp(ltlag)
    tlag2 <- exp(ltlag2)
    d1 <- exp(ld1)

    # Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Three-phase IM absorption into a two-compartment disposition model
    # (Khwarg 2024 Fig. 2). Each IM administration is encoded by the USER as
    # THREE parallel dose records in the event table, all with the same amt:
    #   (a) cmt = "depot"   -- first-order arm, fraction F4, lag ALAG4
    #   (b) cmt = "depot2"  -- first-order arm, fraction F5, lag ALAG5
    #   (c) cmt = "central", rate = -2 -- zero-order arm, fraction 1 - F4 - F5,
    #       modelled duration D2 supplied by dur(central) below.
    # The f() multipliers perform the dose split, so each record carries the
    # full nominal dose. See the vignette for a worked event table.
    d/dt(depot) <- -ka * depot
    d/dt(depot2) <- -ka2 * depot2
    d/dt(central) <- ka * depot + ka2 * depot2 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot
    alag(depot) <- tlag
    f(depot2) <- fdepot2
    alag(depot2) <- tlag2
    f(central) <- 1 - fdepot - fdepot2
    dur(central) <- d1

    # Plasma donepezil concentration. Dose is in mg and vc in L, so central/vc is
    # mg/L; the factor 1000 converts to the ug/L (= ng/mL) scale used in
    # Khwarg 2024 Table 1 and for the additive residual error.
    Cc <- central / vc * 1000

    Cc ~ add(addSd) + prop(propSd)
  })
}
