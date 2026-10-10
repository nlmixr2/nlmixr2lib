Duong_2022_gentamicin_bos <- function() {
  description <- "One-compartment population PK model of intravenous gentamicin in critically ill adult ICU patients (Quebec), the Bos 2019 structure re-estimated by Duong 2022: clearance is a linear function of Cockcroft-Gault creatinine clearance centred at 92 mL/min, volume has no covariate, with IIV on CL and V and combined additive plus proportional residual error (magnitudes carried from the original Bos model because the re-estimated values are not reported)."
  reference <- paste(
    "Duong A, Simard C, Williamson D, Marsot A.",
    "Model Re-Estimation: An Alternative for Poor Predictive Performance during External Evaluations?",
    "Example of Gentamicin in Critically Ill Patients.",
    "Pharmaceutics 2022;14(7):1426. doi:10.3390/pharmaceutics14071426.",
    "Structural model from Bos JC, Prins JM, Misticio MC, Nunguiane G, Lang CN, Beirao JC, Mathot RAA, van Hest RM.",
    "Population pharmacokinetics with Monte Carlo simulations of gentamicin in a population of severely ill adult patients from sub-Saharan Africa.",
    "Antimicrob Agents Chemother 2019;63(4):e02328-18. doi:10.1128/AAC.02328-18.",
    sep = " "
  )
  vignette <- "Duong_2022_gentamicin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance by the Cockcroft-Gault equation (not BSA-normalised)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Duong 2022 Methods Eq. 2: CrCl = (140 - Age) x Body weight (kg) x 1.23",
        "x 0.85 (if female) / Scr, with Scr in umol/L; raw mL/min, NOT",
        "BSA-normalised. Enters CL as 1 + e_crcl_cl x (CRCL - 92), centred at",
        "92 mL/min, the combined-cohort mean of 92.2 mL/min (Table 1)."
      ),
      source_name = "CrCl (CLCG in the Methods)"
    )
  )

  compartmentData <- list(
    central = list(analyte = "gentamicin", units = "mg", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 87L,
    n_studies = 2L,
    age_mean = "59.4 +/- 17.9 years",
    weight_mean = "80.0 +/- 21.5 kg",
    sex_female_pct = 37.9,
    disease_state = paste(
      "Critically ill adult ICU patients receiving intravenous gentamicin:",
      "mostly sepsis at Hopital du Sacre-Coeur de Montreal (HSCM, n = 39)",
      "and mostly endocarditis at the Institut universitaire de cardiologie",
      "et pneumologie de Quebec (IUCPQ, n = 48)."
    ),
    dose_range = "Total daily dose 2.4 +/- 1.1 mg/kg (HSCM 2.9 +/- 0.9, IUCPQ 2.0 +/- 0.7); routine TDM sampling.",
    renal_function = "Serum creatinine 96.9 +/- 66.0 umol/L; Cockcroft-Gault CrCl 92.2 +/- 48.9 mL/min; MDRD eGFR 80.9 +/- 31.9 mL/min (combined cohort).",
    regions = "Canada (Quebec): HSCM 2009-2019 and IUCPQ 2014-2020, retrospective chart review.",
    notes = paste(
      "Demographics from Duong 2022 Table 1 (combined column; 54 M / 33 F).",
      "The Bos 2019 model (developed in non-ICU severely ill sub-Saharan",
      "African adults) was refit in NONMEM 7.5 to the combined HSCM + IUCPQ",
      "data after it failed external evaluation (population MDPE -44.0%,",
      "MADPE 47.8%); the re-estimated model reached MDPE 2.00% and MADPE",
      "30.1% (Table 2)."
    )
  )

  ini({
    # Re-estimated typical values, Duong 2022 Table S2 (Bos et al. row).
    lcl <- log(3.44)
    label("Clearance at CrCl 92 mL/min (L/h)") # Table S2: thetaCL = 3.44
    lvc <- log(22.4)
    label("Volume of distribution (L)") # Table S2: thetaV = 22.4
    e_crcl_cl <- 0.00925
    label("Linear effect of CrCl on CL (per mL/min)") # Table S2 equation CL = 3.44 x (1 + 0.00925 x (CrCl - 92))

    # IIV, Duong 2022 Table S2 as CV%; omega^2 = log(1 + CV^2).
    etalcl ~ 0.070869 # Table S2: IIV CL = 27.1% -> log(1 + 0.271^2)
    etalvc ~ 0.1443 # Table S2: IIV V = 39.4% -> log(1 + 0.394^2)

    # Residual error: Table S2 gives the error model (Mixed) but no
    # re-estimated magnitudes, so the original Bos values used for the
    # external evaluation (Table S1) are carried forward.
    propSd <- 0.32
    label("Proportional residual error (fraction)") # Table S1 (original Bos model): Proportional 32%
    addSd <- 0.056
    label("Additive residual error (mg/L)") # Table S1 (original Bos model): Additive 0.056 mg/L
  })

  model({
    cl <- exp(lcl + etalcl) * (1 + e_crcl_cl * (CRCL - 92))
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
