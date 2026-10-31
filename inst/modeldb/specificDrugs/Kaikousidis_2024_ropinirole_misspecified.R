Kaikousidis_2024_ropinirole_misspecified <- function() {
  description <- paste(
    "One-compartment population PK model for ropinirole after a single 2 mg oral",
    "dose of a prolonged-release tablet (Requip) in healthy adult volunteers",
    "(Kaikousidis 2024, Model 2). This is the earlier, simpler model the authors",
    "kept as a deliberately misspecified comparator: it does not capture the",
    "second concentration peak, and its proportional residual error (0.33) is",
    "larger than that of the final model (Kaikousidis_2024_ropinirole). Linear",
    "elimination and a double first-order absorption without lag times: a",
    "fraction F1 of the dose is absorbed with rate ka1 and the remaining 1 - F1",
    "simultaneously with rate ka2. Log-normal IIV on ka1, ka2, V and CL and",
    "logit-normal IIV on F1; proportional residual error. Each dose is given as",
    "two event records, one into depot and one into depot2, with the same amount.",
    sep = " "
  )
  reference <- paste(
    "Kaikousidis C, Bies RR, Dokoumetzidis A.",
    "Simulating realistic patient profiles from pharmacokinetic models by a",
    "machine learning postprocessing correction of residual variability.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(9):1476-1487.",
    "doi:10.1002/psp4.13182.",
    sep = " "
  )
  vignette <- "Kaikousidis_2024_ropinirole"
  units <- list(time = "h", dosing = "mg", concentration = "ng/L")

  # The two absorption arms each need a dose record, so each dose enters two
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
      "Same dataset as Kaikousidis_2024_ropinirole. Demographics from",
      "Supplementary Table S1 (sex reported as 73.3% male). Model fitted in",
      "Monolix 2023R1 (SAEM); OFV 5098.051 versus 4831.367 for the final model."
    )
  )

  ini({
    # Structural parameters: Kaikousidis 2024 Supplementary Table S3 (Monolix
    # population parameter values for Model 2). CL and V are apparent values;
    # the total fraction absorbed over the two arms is 1.
    lka <- log(0.1214); label("Absorption rate constant of the first arm, ka1 (1/h)") # Table S3: ka1_pop = 0.1214 h-1 (RSE 5.24%)
    lka2 <- log(0.2614); label("Absorption rate constant of the second arm, ka2 (1/h)") # Table S3: ka2_pop = 0.2614 h-1 (RSE 147%)
    logitfdepot <- logit(0.9987); label("Logit of the fraction of the dose absorbed through the first arm, F1 (unitless)") # Table S3: F1_pop = 0.9987 (RSE 0.244%)
    lvc <- log(653.4056); label("Apparent volume of distribution V (L)") # Table S3: V_pop = 653.4056 L (RSE 5.3%)
    lcl <- log(80.0181); label("Apparent clearance Cl (L/h)") # Table S3: Cl_pop = 80.0181 L/h (RSE 6.55%)

    # IIV: Monolix reports omega as the SD of the random effect on the
    # transformed scale (log for the positive parameters, logit for F1); the
    # variances below are omega^2.
    etalka ~ 0.03017169 # Table S3: omega_ka1 = 0.1737 (SD)
    etalka2 ~ 2.941225 # Table S3: omega_ka2 = 1.715 (SD; RSE not estimable, printed NaN)
    etalogitfdepot ~ 0.36735721 # Table S3: omega_F1 = 0.6061 (SD, logit scale)
    etalvc ~ 0.03258025 # Table S3: omega_V = 0.1805 (SD)
    etalcl ~ 0.11737476 # Table S3: omega_Cl = 0.3426 (SD)

    propSd <- 0.3293; label("Proportional residual error (fraction)") # Table S3: Prop. RUV = 0.3293 (RSE 3.61%)
  })

  model({
    ka <- exp(lka + etalka)
    ka2 <- exp(lka2 + etalka2)
    fdepot <- expit(logitfdepot + etalogitfdepot)
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(depot2) <- -ka2 * depot2
    d/dt(central) <- ka * depot + ka2 * depot2 - kel * central

    f(depot) <- fdepot
    f(depot2) <- 1 - fdepot

    # Dose in mg, volume in L; the paper reports concentrations in ng/L.
    Cc <- 1e6 * central / vc
    Cc ~ prop(propSd)
  })
}
