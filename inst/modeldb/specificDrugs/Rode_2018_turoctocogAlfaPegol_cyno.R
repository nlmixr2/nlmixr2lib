Rode_2018_turoctocogAlfaPegol_cyno <- function() {
  description <- paste(
    "Preclinical (cynomolgus monkey). Two-compartment population PK model",
    "for turoctocog alfa pegol (N8-GP, glycoPEGylated recombinant factor",
    "VIII) after subcutaneous and intravenous administration in",
    "cynomolgus monkeys (Rode 2018). First-order absorption from a",
    "subcutaneous depot with bioavailability 26% for the 500 IU/mL",
    "(0.2 mL/kg) formulation, reduced by 32% for the 5000 IU/mL",
    "(0.02 mL/kg) low-volume formulation; the high intravenous dose",
    "(250 IU/kg) lowers the central volume by 49% and clearance by 29%",
    "relative to 100 IU/kg. IIV on the absorption rate constant and on",
    "bioavailability; proportional residual error. Parameters are per kg",
    "of body weight, so doses are in IU/kg and FVIII activity is in IU/mL.",
    "The absorption rate constant and bioavailability were carried into",
    "the human prediction model Rode_2018_turoctocogAlfaPegol_human.",
    sep = " "
  )
  reference <- paste(
    "Rode F, Almholt K, Petersen M, Kreilgaard M, Kjalke M, Karpf DM,",
    "Groth AV, Johansen PB, Larsen LF, Loftager M, Haaning J. Preclinical",
    "pharmacokinetics and biodistribution of subcutaneously administered",
    "glycoPEGylated recombinant factor VIII (N8-GP) and development of a",
    "human pharmacokinetic prediction model. J Thromb Haemost.",
    "2018;16(6):1141-1152. doi:10.1111/jth.14013.",
    "Parameter estimates from Supporting Information Table S2.",
    sep = " "
  )
  vignette <- "Rode_2018_turoctocogAlfaPegol"
  units <- list(time = "h", dosing = "IU/kg", concentration = "IU/mL")

  covariateData <- list(
    FORM_TUROCTOCOGPEGOL_SOLN5000 = list(
      description = "Subcutaneous dosing-solution strength indicator: 1 = 5000 IU/mL (0.02 mL/kg), 0 = 500 IU/mL (0.2 mL/kg)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (500 IU/mL formulation injected at 0.2 mL/kg)",
      notes = paste(
        "Rode 2018 Table S2 'Relative effect of low SC dose volume on F' =",
        "-32%. Both formulations delivered 100 IU/kg, so the contrast is the",
        "injected volume (0.02 vs 0.2 mL/kg). Set to 0 on intravenous records",
        "(it acts only on depot bioavailability)."
      ),
      source_name = "SC dose volume / formulation (500 vs 5000 IU/mL)"
    ),
    DOSE_HIGH = list(
      description = "High intravenous dose indicator: 1 = 250 IU/kg IV, 0 = 100 IU/kg IV or any SC dose",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (100 IU/kg IV and the 100 IU/kg SC doses)",
      notes = paste(
        "Rode 2018 Table S2 'Relative effect of high IV dose' on V (-49%) and",
        "Cl (-29%). Only 250 IU/kg intravenously maps to 1. Time-fixed per",
        "animal: each monkey received a single route and dose level."
      ),
      source_name = "IV dose (100 vs 250 IU/kg)"
    )
  )

  # Amounts are per kilogram of body weight because Table S2 reports V/Bw,
  # Cl/Bw, V2/Bw and Cl2/Bw and the doses are in IU/kg.
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
    ),
    peripheral1 = list(
      analyte = "turoctocog alfa pegol",
      units = "IU/kg",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "monkey (cynomolgus, Macaca fascicularis)",
    n_subjects = 11L,
    n_studies = 1L,
    disease_state = "healthy (endogenous FVIII present but not measured by the N8-GP-specific assay)",
    dose_range = paste(
      "100 IU/kg SC once daily on days 1, 2 and 3 as a 500 IU/mL (n = 3) or",
      "5000 IU/mL (n = 3) formulation; single IV dose of 100 IU/kg (n = 2)",
      "or 250 IU/kg (n = 3)"
    ),
    regions = "monkeys from Bioculture (Mauritius)",
    notes = paste(
      "Rode 2018 Methods, 'PK studies in cynomolgus monkeys': SC injection",
      "in the lateral upper hind leg, IV via the saphenous vein, sampling",
      "over 96 h after the first dose. FVIII:C measured with an",
      "anti-PEG-capture chromogenic assay specific for N8-GP, so endogenous",
      "monkey FVIII does not contribute and no baseline term is needed.",
      "Fitted in Phoenix NLME 7.0 (FOCE-ELS)."
    )
  )

  ini({
    lvc <- log(61); label("Central volume of distribution per kg body weight (mL/kg)") # Table S2 'V/Bw' = 61 mL/kg [95% CI 55, 68]
    lcl <- log(3.0); label("Clearance per kg body weight (mL/h/kg)") # Table S2 'Cl/Bw' = 3.0 mL/kg/h [95% CI 2.7, 3.2]
    lvp <- log(15); label("Peripheral volume of distribution per kg body weight (mL/kg)") # Table S2 'V2/Bw' = 15 mL/kg [95% CI 5.9, 25]
    lq <- log(1.2); label("Intercompartmental clearance per kg body weight (mL/h/kg)") # Table S2 'Cl2/Bw' = 1.2 mL/kg/h [95% CI 0.32, 2.0]
    lka <- log(0.33); label("First-order SC absorption rate constant (1/h)") # Table S2 'ka' = 0.33 1/h [95% CI 0.19, 0.47]
    lfdepot <- log(0.26); label("SC bioavailability, 500 IU/mL formulation (fraction)") # Table S2 'F' = 26% [95% CI 21, 31]

    e_form_turoctocogpegol_soln5000_fdepot <- -0.32; label("Fractional change in F for the 5000 IU/mL low-volume formulation (unitless)") # Table S2 'Relative effect of low SC dose volume on F' = -32% [95% CI -51, -13]
    e_dose_high_vc <- -0.49; label("Fractional change in V for the 250 IU/kg IV dose (unitless)") # Table S2 'Relative effect of high IV dose on V' = -49% [95% CI -64, -34]
    e_dose_high_cl <- -0.29; label("Fractional change in CL for the 250 IU/kg IV dose (unitless)") # Table S2 'Relative effect of high IV dose on Cl' = -29% [95% CI -42, -16]

    etalka ~ 0.14177 # Table S2 IIV ka 39% CV -> log(1 + 0.39^2)
    etalfdepot ~ 0.019414 # Table S2 IIV F 14% CV -> log(1 + 0.14^2)

    propSd <- 0.39; label("Proportional residual error (fraction)") # Table S2 proportional residual error = 39% [95% CI 34, 45]
  })

  model({
    vc <- exp(lvc) * (1 + e_dose_high_vc * DOSE_HIGH)
    cl <- exp(lcl) * (1 + e_dose_high_cl * DOSE_HIGH)
    vp <- exp(lvp)
    q <- exp(lq)
    ka <- exp(lka + etalka)
    fdepot <- exp(lfdepot + etalfdepot) *
      (1 + e_form_turoctocogpegol_soln5000_fdepot * FORM_TUROCTOCOGPEGOL_SOLN5000)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
