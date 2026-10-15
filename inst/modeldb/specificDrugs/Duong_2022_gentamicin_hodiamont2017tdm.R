Duong_2022_gentamicin_hodiamont2017tdm <- function() {
  description <- "Two-compartment population PK model of intravenous gentamicin in critically ill adult ICU patients (Quebec), the Hodiamont 2017 (Ther Drug Monit) structure as specified and re-estimated by Duong 2022: CL and Q scale allometrically (exponent 0.75) and V1 linearly with total body weight, V2 has no covariate, with IIV on CL and V1 and combined additive plus proportional residual error (magnitudes carried from the original model because the re-estimated values are not reported)."
  reference <- paste(
    "Duong A, Simard C, Williamson D, Marsot A.",
    "Model Re-Estimation: An Alternative for Poor Predictive Performance during External Evaluations?",
    "Example of Gentamicin in Critically Ill Patients.",
    "Pharmaceutics 2022;14(7):1426. doi:10.3390/pharmaceutics14071426.",
    "Structural model from Hodiamont CJ, Janssen JM, de Jong MD, Mathot RA, Juffermans NP, van Hest RM.",
    "Therapeutic Drug Monitoring of Gentamicin Peak Concentrations in Critically Ill Patients.",
    "Ther Drug Monit 2017;39(5):522-530. doi:10.1097/FTD.0000000000000432.",
    sep = " "
  )
  vignette <- "Duong_2022_gentamicin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Reference 70 kg; exponent 0.75 on CL and Q, linear on V1 (Duong 2022",
        "Tables S1 and S2). The original Hodiamont 2017 TDM model retained no",
        "body-weight covariate (see Hodiamont_2017_gentamicin); the scaling",
        "is part of the model as Duong 2022 specified and refit it."
      ),
      source_name = "BW"
    )
  )

  compartmentData <- list(
    central = list(analyte = "gentamicin", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "gentamicin", units = "mg", specimen = "serum", verified = TRUE)
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
    regions = "Canada (Quebec): HSCM 2009-2019 and IUCPQ 2014-2020, retrospective chart review.",
    notes = paste(
      "Demographics from Duong 2022 Table 1 (combined column; 54 M / 33 F).",
      "The Hodiamont TDM model was refit in NONMEM 7.5 to the combined",
      "HSCM + IUCPQ data after it failed external evaluation (population",
      "MDPE -31.7%, MADPE 48.8%); the re-estimated model reached MDPE 6.03%",
      "and MADPE 39.2% (Table 2)."
    )
  )

  ini({
    # Re-estimated typical values, Duong 2022 Table S2 (Hodiamont et al. [16] row).
    lcl <- log(1.63)
    label("Clearance at 70 kg (L/h)") # Table S2: thetaCL = 1.63
    lvc <- log(8.67)
    label("Central volume at 70 kg (L)") # Table S2: thetaV1 = 8.67
    lq <- log(0.943)
    label("Intercompartmental clearance at 70 kg (L/h)") # Table S2: thetaQ = 0.943
    lvp <- log(6.78)
    label("Peripheral volume (L)") # Table S2: thetaV2 = 6.78
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of body weight on CL (unitless)") # Table S2 equation CL = 1.63 x (BW/70)^0.75
    e_wt_q <- fixed(0.75)
    label("Allometric exponent of body weight on Q (unitless)") # Table S2 equation Q = 0.943 x (BW/70)^0.75
    e_wt_vc <- fixed(1)
    label("Exponent of body weight on V1 (unitless)") # Table S2 equation V1 = 8.67 x (BW/70)

    # IIV, Duong 2022 Table S2 as CV%; omega^2 = log(1 + CV^2). Duong 2022
    # reports no CL-V1 correlation, so the block is diagonal.
    etalcl ~ 0.2626 # Table S2: IIV CL = 54.8% -> log(1 + 0.548^2)
    etalvc ~ 0.18967 # Table S2: IIV V1 = 45.7% -> log(1 + 0.457^2)

    # Residual error: Table S2 gives the error model (Mixed) but no
    # re-estimated magnitudes, so the original values used for the external
    # evaluation (Table S1) are carried forward.
    propSd <- 0.194
    label("Proportional residual error (fraction)") # Table S1 (original model): Proportional 19.4%
    addSd <- 0.13
    label("Additive residual error (mg/L)") # Table S1 (original model): Additive 0.13 mg/L
  })

  model({
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q <- exp(lq) * (WT / 70)^e_wt_q
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
