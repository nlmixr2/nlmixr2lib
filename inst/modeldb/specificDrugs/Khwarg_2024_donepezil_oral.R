Khwarg_2024_donepezil_oral <- function() {
  description <- "Two-compartment population PK model for oral donepezil (Aricept 10 mg tablet) with first-order absorption and an absorption lag time, in healthy adult men"
  reference <- paste(
    "Khwarg J, Lee H, Yu KS, Seol E, Chung JY.",
    "Population Pharmacokinetic Modeling and Simulation for Dose Optimization",
    "of GB-5001, a Long-Acting Intramuscular Injection of Donepezil,",
    "in Healthy Participants. Neurol Ther. 2024;13(5):1453-1466.",
    "doi:10.1007/s40120-024-00643-4.",
    "Companion intramuscular model: modellib('Khwarg_2024_donepezil_im').",
    sep = " "
  )
  vignette <- "Khwarg_2024_donepezil"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  compartmentData <- list(
    depot = list(analyte = "donepezil", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "donepezil", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "donepezil", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 12,
    n_studies = 1,
    age_range = "18-55 years (eligibility); cohort mean 39.5 (SD 10.2) years",
    weight_range = "cohort mean 82.4 (SD 9.1) kg",
    sex_female_pct = 0,
    race_ethnicity = c(White = 100),
    disease_state = "healthy",
    dose_range = "10 mg single oral dose (Aricept tablet)",
    regions = "USA (Advarra IRB, Columbia MD); sponsor G2GBIO, Republic of Korea",
    notes = paste(
      "Healthy non-smoking men, BMI 18.5-30.0 kg/m2, from the oral arm (cohort 4)",
      "of the NCT05525780 phase 1 dose-escalation study; demographics in Khwarg 2024",
      "'Demographics' (mean BMI 27.1 (SD 2.1) kg/m2, height 174.3 (SD 4.4) cm).",
      "The oral and intramuscular arms were fitted in a single NONMEM run but share",
      "no parameters, so they are distributed as two model files.",
      sep = " "
    )
  )

  ini({
    # Structural parameters -- Khwarg 2024 Table 2, 'Oral' block, 'Final model' estimate column.
    lka <- log(0.203); label("Absorption rate constant from depot (KA, 1/h)") # Table 2 oral: KA = 0.203 1/h (RSE 11.4%; bootstrap 0.21, 95% CI 0.165-0.257)
    lcl <- log(14.3); label("Clearance (CL, L/h)") # Table 2 oral: CL = 14.3 L/h (RSE 6.5%; bootstrap 14.35, 95% CI 12.602-16.282)
    lvc <- log(39.5); label("Central volume of distribution (V2, L)") # Table 2 oral: V2 = 39.5 L (RSE 36.5%; bootstrap 40.737, 95% CI 19.526-73.488)
    lq <- log(84.9); label("Inter-compartmental clearance (Q, L/h)") # Table 2 oral: Q = 84.9 L/h (RSE 9.6%; bootstrap 86.514, 95% CI 70.362-103.327)
    lvp <- log(1080); label("Peripheral volume of distribution (V3, L)") # Table 2 oral: V3 = 1080 L (RSE 5.9%; bootstrap 1073.55, 95% CI 957.33-1202.795)
    ltlag <- log(0.931); label("Absorption lag time (ALAG1, h)") # Table 2 oral: ALAG1 = 0.931 h (RSE 1.7%; bootstrap 0.93, 95% CI 0.887-0.956)
    lfdepot <- fixed(log(1)); label("Oral bioavailability (F1, fraction)") # Table 2 oral: F1 = 1 (fixed), anchored to the high oral bioavailability of Aricept (Khwarg 2024 ref [14])

    # IIV -- Khwarg 2024 Table 2 reports the log-scale variance directly; the
    # parenthesised percentage is the derived CV via footnote b,
    # CV = sqrt(exp(omega^2) - 1) * 100. Values below are the variances.
    # CL / central-V block: sqrt(exp(0.0475) - 1) = 22.1%; sqrt(exp(1.59) - 1) = 197.6%;
    # correlation -0.144 / sqrt(0.0475 * 1.59) = -0.524, matching the reported -52.0%.
    etalcl + etalvc ~ c(0.0475, -0.144, 1.59) # Table 2 oral: IIV CL 0.0475 (22.1%), covariance CL vs V2 -0.144 (-52.0%), IIV V2 1.59 (197.6%)
    etalvp ~ 0.0153 # Table 2 oral: IIV V3 = 0.0153 (12.4%); sqrt(exp(0.0153) - 1) = 12.4%
    etalfdepot ~ 0.0162 # Table 2 oral: IIV F1 = 0.0162 (12.8%); sqrt(exp(0.0162) - 1) = 12.8%

    # Residual error -- reported on the SD scale: the additive term carries linear
    # concentration units (pg/mL, not pg^2/mL^2), so both terms are standard
    # deviations rather than NONMEM $SIGMA variances.
    propSd <- 0.137; label("Proportional residual error (fraction)") # Table 2 oral: proportional error = 0.137 (RSE 11.5%; bootstrap 0.132, 95% CI 0.106-0.162)
    addSd <- 0.138; label("Additive residual error (ug/L)") # Table 2 oral: additive error = 138 pg/mL = 0.138 ug/L (RSE 13%; bootstrap 137.786, 95% CI 85.692-202.99)
  })

  model({
    # Individual parameters
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp + etalvp)
    fdepot <- exp(lfdepot + etalfdepot)
    tlag <- exp(ltlag)

    # Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Two-compartment disposition with first-order absorption from a single
    # oral depot (Khwarg 2024 Fig. 2: Oral Depot (1) -> Central (2) <-> Peripheral (3)).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Bioavailability and absorption lag on the oral depot
    f(depot) <- fdepot
    alag(depot) <- tlag

    # Plasma donepezil concentration. Dose is in mg and vc in L, so central/vc is
    # mg/L; the factor 1000 converts to the ug/L (= ng/mL) scale used in
    # Khwarg 2024 Table 1 and for the additive residual error.
    Cc <- central / vc * 1000

    Cc ~ add(addSd) + prop(propSd)
  })
}
