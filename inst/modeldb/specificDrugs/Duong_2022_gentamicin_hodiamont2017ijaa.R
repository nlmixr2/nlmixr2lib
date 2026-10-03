Duong_2022_gentamicin_hodiamont2017ijaa <- function() {
  description <- "Two-compartment population PK model of intravenous gentamicin in critically ill adult ICU patients (Quebec), the Hodiamont 2017 (Int J Antimicrob Agents) structure re-estimated by Duong 2022: clearance scales allometrically with ideal body weight, central volume scales linearly with ideal body weight and with a power function of serum albumin, with IIV on CL and V1 and a proportional residual error whose magnitude is reported in neither the original nor the re-estimated model (fixed to 0)."
  reference <- paste(
    "Duong A, Simard C, Williamson D, Marsot A.",
    "Model Re-Estimation: An Alternative for Poor Predictive Performance during External Evaluations?",
    "Example of Gentamicin in Critically Ill Patients.",
    "Pharmaceutics 2022;14(7):1426. doi:10.3390/pharmaceutics14071426.",
    "Structural model from Hodiamont CJ, Juffermans NP, Bouman CS, de Jong MD, Mathot RA, van Hest RM.",
    "Determinants of gentamicin concentrations in critically ill patients: a population pharmacokinetic analysis.",
    "Int J Antimicrob Agents 2017;49(2):204-211. doi:10.1016/j.ijantimicag.2016.10.022.",
    sep = " "
  )
  vignette <- "Duong_2022_gentamicin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    IBW = list(
      description = "Ideal body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Reference 70 kg; allometric exponent 0.75 on CL and linear on V1",
        "(Duong 2022 Table S2). Duong 2022 does not state the IBW formula, and",
        "height is not among the data it extracted from the medical records",
        "(Methods 2.1), so its fit probably held IBW at the typical value",
        "(the stated rule for covariates missing from the dataset, Methods 2.3)."
      ),
      source_name = "IBW"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Reference 22 g/L; power exponent -0.833 on V1 (Duong 2022 Table S2,",
        "ALBM). Recorded only at IUCPQ (29.0 +/- 5.6 g/L, Table 1); HSCM",
        "patients were assigned the typical value per Methods 2.3."
      ),
      source_name = "ALBM"
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
    albumin = "29.0 +/- 5.6 g/L (IUCPQ only; not recorded at HSCM).",
    regions = "Canada (Quebec): HSCM 2009-2019 and IUCPQ 2014-2020, retrospective chart review.",
    notes = paste(
      "Demographics from Duong 2022 Table 1 (combined column; 54 M / 33 F).",
      "The Hodiamont IJAA model (critically ill patients on or off CVVH) was",
      "refit in NONMEM 7.5 to the combined HSCM + IUCPQ data after it failed",
      "external evaluation (population MDPE 66.1%, MADPE 69.9%); the",
      "re-estimated model reached MDPE 2.20% and MADPE 36.9% (Table 2)."
    )
  )

  ini({
    # Re-estimated typical values, Duong 2022 Table S2 (Hodiamont et al. [15] row).
    lcl <- log(2.12)
    label("Clearance at IBW 70 kg (L/h)") # Table S2: thetaCL = 2.12 (also Table S3); the Table S2 equation column prints 2.11
    lvc <- log(23.9)
    label("Central volume at IBW 70 kg and albumin 22 g/L (L)") # Table S2: thetaV1 = 23.9 (also Table S3); the Table S2 equation column still prints the original 21.2
    lq <- log(1.95)
    label("Intercompartmental clearance (L/h)") # Table S2: thetaQ = 1.95
    lvp <- log(18.1)
    label("Peripheral volume (L)") # Table S2: thetaV2 = 18.1
    e_ibw_cl <- fixed(0.75)
    label("Allometric exponent of IBW on CL (unitless)") # Table S2 equation CL = thetaCL x (IBW/70)^0.75
    e_ibw_vc <- fixed(1)
    label("Exponent of IBW on V1 (unitless)") # Table S2 equation V1 = thetaV1 x (IBW/70)
    e_alb_vc <- fixed(-0.833)
    label("Power exponent of albumin on V1 (unitless)") # Table S2 equation (ALBM/22)^-0.833, unchanged from the original model in Table S1

    # IIV, Duong 2022 Table S2 as CV%; omega^2 = log(1 + CV^2).
    etalcl ~ 0.21598 # Table S2: IIV CL = 49.1% -> log(1 + 0.491^2)
    etalvc ~ 0.18967 # Table S2: IIV V1 = 45.7% -> log(1 + 0.457^2)

    # Residual error: Table S1 and Table S2 both give the error model as
    # Proportional but neither reports its magnitude.
    propSd <- fixed(0)
    label("Proportional residual error (fraction; magnitude not reported, set to 0)") # Tables S1 and S2: proportional error, magnitude blank
  })

  model({
    cl <- exp(lcl + etalcl) * (IBW / 70)^e_ibw_cl
    vc <- exp(lvc + etalvc) * (IBW / 70)^e_ibw_vc * (ALB / 22)^e_alb_vc
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
