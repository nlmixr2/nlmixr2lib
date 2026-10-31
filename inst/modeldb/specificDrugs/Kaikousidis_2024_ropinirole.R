Kaikousidis_2024_ropinirole <- function() {
  description <- paste(
    "One-compartment population PK model for ropinirole after a single 2 mg oral",
    "dose of a prolonged-release tablet (Requip) in healthy adult volunteers",
    "(Kaikousidis 2024, final Model 1). Linear elimination and a double",
    "absorption: a fraction F1 of the dose passes through a Savic",
    "transit-compartment chain (mean transit time Mtt1, transit rate constant",
    "Ktr1, non-integer number of transits Mtt1 * Ktr1 - 1) into an absorption",
    "compartment emptied with first-order rate ka1, while the remaining 1 - F1",
    "is absorbed simultaneously with first-order rate ka2 after a lag time",
    "Tlag2 = Mtt1 + diffTlag2 that always exceeds the mean transit time.",
    "Log-normal IIV on all structural parameters and logit-normal IIV on F1;",
    "proportional residual error. Each dose is given as two event records, one",
    "into depot and one into depot2, with the same amount.",
    sep = " "
  )
  reference <- paste(
    "Kaikousidis C, Bies RR, Dokoumetzidis A.",
    "Simulating realistic patient profiles from pharmacokinetic models by a",
    "machine learning postprocessing correction of residual variability.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(9):1476-1484.",
    "doi:10.1002/psp4.13182.",
    sep = " "
  )
  vignette <- "Kaikousidis_2024_ropinirole"
  units <- list(time = "h", dosing = "mg", concentration = "ng/L")

  # The lagged arm needs its own dose record, so each dose enters two
  # compartments with the same amount.
  dosing <- c("depot", "depot2")

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "ropinirole", units = "mg", specimen = "administration site", verified = TRUE),
    depot2 = list(analyte = "ropinirole", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ropinirole", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 30L,
    n_studies = 1L,
    n_observations = "18 samples per subject (dense sampling over 60 h)",
    age_range = "19-44 years",
    age_mean = "32 years (SD 8)",
    weight_range = "54.7-94.8 kg",
    weight_mean = "74.5 kg (SD 10.3)",
    height_mean = "172.3 cm (SD 8.0; range 158-188)",
    bmi_mean = "25.1 kg/m^2 (SD 2.7; range 19.4-29.8)",
    sex_female_pct = 26.7,
    race_ethnicity = c(Black = 56.6, Caucasian = 20, Asian = 10, Other = 13.4),
    disease_state = "Healthy volunteers.",
    dose_range = "Single oral dose of 2 mg ropinirole as a prolonged-release tablet (Requip, GlaxoSmithKline).",
    notes = paste(
      "Demographics from Supplementary Table S1 (sex reported as 73.3% male).",
      "Sampling at 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 24, 36, 48",
      "and 60 h (Supplementary Information, Dataset). Model fitted in Monolix",
      "2023R1 (SAEM). The paper's machine-learning (random forest) correction of",
      "the individual predictions is a post-processing step outside this model."
    )
  )

  ini({
    # Structural parameters: Kaikousidis 2024 Supplementary Table S2 (Monolix
    # population parameter values for Model 1). CL and V are apparent values;
    # the total fraction absorbed over the two arms is 1.
    lka <- log(0.09886); label("Absorption rate constant of the transit arm, ka1 (1/h)") # Table S2: ka1_pop = 0.09886 h-1 (RSE 13.1%)
    lka2 <- log(0.2501); label("Absorption rate constant of the lagged arm, ka2 (1/h)") # Table S2: ka2_pop = 0.2501 h-1 (RSE 31.6%)
    logitfdepot <- logit(0.8071); label("Logit of the fraction of the dose absorbed through the transit arm, F1 (unitless)") # Table S2: F1_pop = 0.8071 (RSE 3.85%)
    lmtt <- log(0.5545); label("Mean transit time of the transit arm, Mtt1 (h)") # Table S2: Mtt1_pop = 0.5545 h (RSE 64.4%)
    lktr <- log(1.8352); label("Transit rate constant of the transit arm, Ktr1 (1/h)") # Table S2: Ktr1_pop = 1.8352 h-1 (RSE 162%)
    ldtlag2 <- log(11.985); label("Increment of the lagged arm's lag time over Mtt1, diffTlag2 (h)") # Table S2: diffTlag2_pop = 11.985 h (RSE 6.23%)
    lvc <- log(311.7171); label("Apparent volume of distribution V (L)") # Table S2: V_pop = 311.7171 L (RSE 16.5%)
    lcl <- log(72.2331); label("Apparent clearance Cl (L/h)") # Table S2: Cl_pop = 72.2331 L/h (RSE 7.95%)

    # IIV: Monolix reports omega as the SD of the random effect on the
    # transformed scale (log for the positive parameters, logit for F1); the
    # variances below are omega^2.
    etalka ~ 0.04309776 # Table S2: omega_ka1 = 0.2076 (SD)
    etalka2 ~ 0.02745649 # Table S2: omega_ka2 = 0.1657 (SD)
    etalogitfdepot ~ 0.25867396 # Table S2: omega_F1 = 0.5086 (SD, logit scale)
    etalmtt ~ 0.295936 # Table S2: omega_Mtt1 = 0.544 (SD)
    etalktr ~ 3.49503025 # Table S2: omega_Ktr1 = 1.8695 (SD)
    etaldtlag2 ~ 0.03825936 # Table S2: omega_diffTlag2 = 0.1956 (SD)
    etalvc ~ 0.04056196 # Table S2: omega_V = 0.2014 (SD)
    etalcl ~ 0.12524521 # Table S2: omega_Cl = 0.3539 (SD)

    propSd <- 0.222; label("Proportional residual error (fraction)") # Table S2: Prop. RUV = 0.222 (RSE 4.16%)
  })

  model({
    ka <- exp(lka + etalka)
    ka2 <- exp(lka2 + etalka2)
    fdepot <- expit(logitfdepot + etalogitfdepot)
    mtt <- exp(lmtt + etalmtt)
    ktr <- exp(lktr + etalktr)
    dtlag2 <- exp(ldtlag2 + etaldtlag2)
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl)

    # Results, Modeling the residual error: Tlag2 = Mtt1 + diffTlag2.
    tlag2 <- mtt + dtlag2
    kel <- cl / vc

    # Monolix transit-compartment absorption (Mtt, Ktr, ka): the number of
    # transit compartments is the non-integer nn = Mtt * Ktr - 1, so
    # Ktr = (nn + 1) / Mtt, and the flow into the absorption compartment is
    # the Savic 2007 gamma density with shape nn + 1 and rate Ktr. With the
    # large IIV on Ktr1, nn falls in (-1, 0) for about half the subjects,
    # where that density is unbounded at the dose time and the ODE solver
    # fails. The model is therefore written in the equivalent integrated
    # form: depot holds all transit-arm drug not yet absorbed (chain plus
    # absorption compartment), the amount still in the chain is
    # F1 * dose * (1 - P(nn + 1, Ktr * tad)) with P the regularized lower
    # incomplete gamma function, and only the remainder is absorbed at ka1.
    # Like the closed-form Savic input, the chain term follows the most
    # recent depot dose only.
    nn <- mtt * ktr - 1
    tad_depot <- tad(depot)
    chain <- depot
    if (tad_depot > 0) {
      chain <- fdepot * podo(depot) * (1 - gammap(nn + 1, ktr * tad_depot))
    }
    abs_cmt <- depot - chain

    d/dt(depot) <- -ka * abs_cmt
    d/dt(depot2) <- -ka2 * depot2
    d/dt(central) <- ka * abs_cmt + ka2 * depot2 - kel * central

    f(depot) <- fdepot
    f(depot2) <- 1 - fdepot
    alag(depot2) <- tlag2

    # Dose in mg, volume in L; the paper reports concentrations in ng/L.
    Cc <- 1e6 * central / vc
    Cc ~ prop(propSd)
  })
}
