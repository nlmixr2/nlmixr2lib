Duong_2022_gentamicin_rea <- function() {
  description <- "One-compartment population PK model of intravenous gentamicin in critically ill adult ICU patients (Quebec), the Rea 2008 structure re-estimated by Duong 2022: clearance is a sigmoid (Hill exponent 1.2) function of MDRD-estimated GFR, central volume scales linearly with body weight, with IIV on maximal clearance, on the GFR at half-maximal clearance and on volume, and combined additive plus proportional residual error. This was the best performing of the four re-estimated models and the one used for the dosing nomogram."
  reference <- paste(
    "Duong A, Simard C, Williamson D, Marsot A.",
    "Model Re-Estimation: An Alternative for Poor Predictive Performance during External Evaluations?",
    "Example of Gentamicin in Critically Ill Patients.",
    "Pharmaceutics 2022;14(7):1426. doi:10.3390/pharmaceutics14071426.",
    "Structural model from Rea RS, Capitano B, Bies R, Bigos KL, Smith R, Lee H.",
    "Suboptimal aminoglycoside dosing in critically ill patients.",
    "Ther Drug Monit 2008;30(6):674-681. doi:10.1097/FTD.0b013e31818b6b2f.",
    sep = " "
  )
  vignette <- "Duong_2022_gentamicin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description = "Estimated glomerular filtration rate (four-variable MDRD equation)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Duong 2022 Methods Eq. 1: eGFR = 186.3 x (Scr/88.4)^-1.154 x Age^-0.203",
        "x 1.212 (if black) x 0.742 (if female), with Scr in umol/L. The paper",
        "labels the result mL/min; the MDRD equation natively returns",
        "mL/min/1.73 m^2 and no BSA de-normalisation is described, so supply",
        "the MDRD value as computed. Enters CL as the sigmoid",
        "CRCL^1.2 / (crcl50^1.2 + CRCL^1.2) (Duong 2022 Table S2). Combined",
        "cohort mean 80.9 +/- 31.9 mL/min (Table 1)."
      ),
      source_name = "eGFR (GFR in Table S2)"
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear scaling of V with reference weight 70 kg (Duong 2022 Table S2: V1 = thetaV x (BW/70)).",
      source_name = "BW"
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
    albumin = "29.0 +/- 5.6 g/L (IUCPQ only; not recorded at HSCM).",
    regions = "Canada (Quebec): HSCM 2009-2019 and IUCPQ 2014-2020, retrospective chart review.",
    notes = paste(
      "Demographics from Duong 2022 Table 1 (combined column; 54 M / 33 F).",
      "The Rea 2008 model was refit in NONMEM 7.5 to the combined HSCM + IUCPQ",
      "data after it failed external evaluation (population MDPE 44.2%,",
      "MADPE 54.1%); the re-estimated model reached MDPE 2.14% and MADPE",
      "28.1% (Table 2)."
    )
  )

  ini({
    # Re-estimated typical values, Duong 2022 Table S2 (Rea et al. row).
    lclmax <- log(9.31)
    label("Maximal gentamicin clearance at saturating eGFR (L/h)") # Table S2: thetaCL,a = 9.31
    lcrcl50 <- log(129)
    label("eGFR at half-maximal clearance (mL/min)") # Table S2: thetaCL,b = 129
    lhill <- fixed(log(1.2))
    label("Hill exponent of the eGFR-clearance sigmoid (unitless)") # Table S2 equation GFR^1.2, unchanged from the original Rea model in Table S1
    lvc <- log(21.7)
    label("Volume of distribution at 70 kg (L)") # Table S2: thetaV = 21.7
    e_wt_vc <- fixed(1)
    label("Exponent of body weight on V (unitless)") # Table S2 equation V1 = thetaV x (BW/70)

    # IIV, Duong 2022 Table S2 as CV%; omega^2 = log(1 + CV^2).
    etalclmax ~ 0.12766 # Table S2: IIV CL,a = 36.9% -> log(1 + 0.369^2)
    etalcrcl50 ~ 0.03294 # Table S2: IIV CL,b = 18.3% -> log(1 + 0.183^2)
    etalvc ~ 0.079683 # Table S2: IIV V = 28.8% -> log(1 + 0.288^2)

    # Residual error: Table S2 gives the error model (Mixed) but not its
    # re-estimated magnitudes; the bootstrap means of the re-estimated model
    # (Table S4) are the only reported values.
    propSd <- 0.344
    label("Proportional residual error (fraction)") # Table S4 bootstrap mean: Proportional 34.4%
    addSd <- 0.279
    label("Additive residual error (mg/L)") # Table S4 bootstrap mean: Additive 0.279 mg/L
  })

  model({
    clmax <- exp(lclmax + etalclmax)
    crcl50 <- exp(lcrcl50 + etalcrcl50)
    hill <- exp(lhill)
    cl <- clmax * CRCL^hill / (crcl50^hill + CRCL^hill)
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc

    kel <- cl / vc

    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
