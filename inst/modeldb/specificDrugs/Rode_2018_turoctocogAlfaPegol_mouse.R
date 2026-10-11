Rode_2018_turoctocogAlfaPegol_mouse <- function() {
  description <- paste(
    "Preclinical (mouse, FVIII knockout). One-compartment population PK",
    "model for turoctocog alfa pegol (N8-GP, glycoPEGylated recombinant",
    "factor VIII) after subcutaneous and intravenous administration in",
    "haemophilia A (F8-KO) mice (Rode 2018). First-order absorption from a",
    "subcutaneous depot with bioavailability 24%, first-order elimination",
    "parameterised by an elimination rate constant, IIV on the volume only",
    "and a proportional residual error. Parameters are per kg of body",
    "weight, so doses are in IU/kg and FVIII activity is in IU/mL.",
    sep = " "
  )
  reference <- paste(
    "Rode F, Almholt K, Petersen M, Kreilgaard M, Kjalke M, Karpf DM,",
    "Groth AV, Johansen PB, Larsen LF, Loftager M, Haaning J. Preclinical",
    "pharmacokinetics and biodistribution of subcutaneously administered",
    "glycoPEGylated recombinant factor VIII (N8-GP) and development of a",
    "human pharmacokinetic prediction model. J Thromb Haemost.",
    "2018;16(6):1141-1152. doi:10.1111/jth.14013.",
    "Parameter estimates from Supporting Information Table S1.",
    sep = " "
  )
  vignette <- "Rode_2018_turoctocogAlfaPegol"
  units <- list(time = "h", dosing = "IU/kg", concentration = "IU/mL")

  # Amounts are per kilogram of body weight because Table S1 reports V/Bw and
  # the doses are in IU/kg, so amount / volume is still IU/mL.
  compartmentData <- list(
    depot = list(
      analyte = "turoctocog alfa pegol",
      units = "IU/kg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "turoctocog alfa pegol",
      units = "IU/kg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "mouse (FVIII exon 16 knockout, B6.129S4-F8tm1Kaz/J)",
    n_subjects = 36L,
    n_studies = 1L,
    age_range = "12-16 weeks",
    disease_state = "haemophilia A model (FVIII knockout); no endogenous FVIII activity",
    dose_range = paste(
      "2500 IU/kg single SC injection in the flank, neck or lower abdomen",
      "(three groups of 9) and 280 IU/kg single IV injection (9 mice)"
    ),
    regions = "Denmark (Novo Nordisk)",
    notes = paste(
      "Rode 2018 Methods, 'PK studies in F8-KO mice' and Data S1 'Animals':",
      "male and female mice, sparse sampling from 1 to 72 h after dosing.",
      "FVIII:C measured with a chromogenic assay (Coatest SP) calibrated",
      "against N8-GP. The three SC injection sites were pooled in one fit",
      "with no site covariate."
    )
  )

  ini({
    lvc <- log(39); label("Central volume of distribution per kg body weight (mL/kg)") # Table S1 'V/Bw' = 39 mL/kg [95% CI 36, 42]
    lkel <- log(0.049); label("Elimination rate constant (1/h)") # Table S1 'Ke' = 0.049 1/h [95% CI 0.043, 0.055]
    lka <- log(0.094); label("First-order SC absorption rate constant (1/h)") # Table S1 'Ka' = 0.094 1/h [95% CI 0.072, 0.12]
    lfdepot <- log(0.24); label("SC bioavailability (fraction)") # Table S1 'F' = 24% [95% CI 21, 27]

    etalvc ~ 0.047268 # Table S1 IIV on V/Bw 22% CV -> log(1 + 0.22^2)

    propSd <- 0.32; label("Proportional residual error (fraction)") # Table S1 proportional error = 32% [95% CI 26, 38]
  })

  model({
    vc <- exp(lvc + etalvc)
    kel <- exp(lkel)
    ka <- exp(lka)
    fdepot <- exp(lfdepot)

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    f(depot) <- fdepot

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
