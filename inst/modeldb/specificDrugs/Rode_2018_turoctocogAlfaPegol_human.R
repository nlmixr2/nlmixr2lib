Rode_2018_turoctocogAlfaPegol_human <- function() {
  description <- paste(
    "One-compartment human PK prediction model for subcutaneous",
    "turoctocog alfa pegol (N8-GP, glycoPEGylated recombinant factor VIII)",
    "in patients with severe haemophilia A (Rode 2018). A hybrid model:",
    "clearance, volume, their IIV and the residual error come from a",
    "one-compartment population PK fit to the intravenous pathfinder 1",
    "trial (26 PK profiles after 25, 50 or 75 IU/kg), while the first-order",
    "subcutaneous absorption rate constant, bioavailability and their IIV",
    "are carried over unchanged from the cynomolgus monkey model of the same",
    "paper (Rode_2018_turoctocogAlfaPegol_cyno). Parameters are per kg of",
    "body weight, so doses are in IU/kg and the observed FVIII activity is",
    "in IU/mL (1 IU/mL = 100% of normal). Used by the authors to predict",
    "steady-state FVIII:C after 12.5 IU/kg subcutaneously once daily versus",
    "50 IU/kg intravenously every fourth day.",
    sep = " "
  )
  reference <- paste(
    "Rode F, Almholt K, Petersen M, Kreilgaard M, Kjalke M, Karpf DM,",
    "Groth AV, Johansen PB, Larsen LF, Loftager M, Haaning J. Preclinical",
    "pharmacokinetics and biodistribution of subcutaneously administered",
    "glycoPEGylated recombinant factor VIII (N8-GP) and development of a",
    "human pharmacokinetic prediction model. J Thromb Haemost.",
    "2018;16(6):1141-1152. doi:10.1111/jth.14013.",
    "Clearance, volume and residual error from Table S3 (pathfinder 1,",
    "NCT01205724; Tiede A et al. J Thromb Haemost 2013;11:670-678,",
    "doi:10.1111/jth.12161, and data on file); ka and F from Table S2.",
    sep = " "
  )
  vignette <- "Rode_2018_turoctocogAlfaPegol"
  units <- list(time = "h", dosing = "IU/kg", concentration = "IU/mL")

  # Amounts are per kilogram of body weight because Table S3 reports V/Bw and
  # Cl/Bw and the paper doses in IU/kg, so amount / volume is still IU/mL.
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
    species = "human",
    n_subjects = 26L,
    n_studies = 1L,
    disease_state = "severe haemophilia A (FVIII < 1%), previously treated patients",
    dose_range = paste(
      "Fitted data: single intravenous doses of 25, 50 or 75 IU/kg",
      "(n = 7, 9 and 10). Prediction scenarios: 12.5 IU/kg subcutaneously",
      "once daily and 50 IU/kg intravenously every 4 days."
    ),
    weight_reference = "70 kg (Table S3 title; parameters are per kg body weight)",
    regions = "multinational (pathfinder 1)",
    notes = paste(
      "Clearance, volume and their IIV are from the pathfinder 1 intravenous",
      "trial (Rode 2018 Methods, 'Human PK prediction model'; NCT01205724).",
      "Rode 2018 does not tabulate the pathfinder 1 demographics. The",
      "subcutaneous absorption rate constant and bioavailability were not",
      "estimated in humans: the authors assumed them equal to the",
      "cynomolgus monkey estimates (Table S2), by analogy with rFVIIa,",
      "glycoPEGylated rFVIIa and rFIX."
    )
  )

  ini({
    # Human elimination parameters: Rode 2018 Table S3 (pathfinder 1 IV fit).
    lcl <- log(1.7); label("Clearance per kg body weight (mL/h/kg)") # Table S3 'Cl/Bw' = 1.7 mL/kg/h [95% CI 1.6, 1.8]
    lvc <- log(45); label("Volume of distribution per kg body weight (mL/kg)") # Table S3 'V/Bw' = 45 mL/kg [95% CI 44, 46]

    # Subcutaneous absorption: carried over from the cynomolgus model,
    # Rode 2018 Table S2 (500 IU/mL, 0.2 mL/kg reference formulation).
    lka <- log(0.33); label("First-order SC absorption rate constant (1/h)") # Table S2 'ka' = 0.33 1/h [95% CI 0.19, 0.47]
    lfdepot <- log(0.26); label("SC bioavailability (fraction)") # Table S2 'F' = 26% [95% CI 21, 31]

    # IIV: percentages read as CV%, omega^2 = log(CV^2 + 1).
    # Table S3: IIV CL 40% -> 0.14842; IIV V 25% -> 0.060625;
    # 'Cl-V covariance' 24% read as the CL-V correlation coefficient 0.24:
    # cov = 0.24 * sqrt(0.14842 * 0.060625) = 0.022780.
    etalcl + etalvc ~ c(0.14842, 0.022780, 0.060625)
    etalka ~ 0.14177 # Table S2 IIV ka 39% -> log(1 + 0.39^2)
    etalfdepot ~ 0.019414 # Table S2 IIV F 14% -> log(1 + 0.14^2)

    propSd <- 0.117; label("Proportional residual error (fraction)") # Table S3 proportional error = 11.7%
  })

  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    ka <- exp(lka + etalka)
    fdepot <- exp(lfdepot + etalfdepot)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    f(depot) <- fdepot

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
