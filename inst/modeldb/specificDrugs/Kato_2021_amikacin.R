Kato_2021_amikacin <- function() {
  description <- "One-compartment IV population PK model for amikacin in hospitalized elderly Japanese patients (aged 70 years and over) treated for serious infections, with clearance directly proportional to Cockcroft-Gault creatinine clearance. Estimated with the Pmetrics non-parametric adaptive grid (NPAG); the between-subject variability is a log-normal approximation of the published non-parametric distribution (Kato 2021, n = 15 patients, 33 samples)."
  reference <- "Kato H, Parker SL, Roberts JA, Hagihara M, Asai N, Yamagishi Y, Paterson DL, Mikamo H. Population Pharmacokinetics Analysis of Amikacin Initial Dosing Regimen in Elderly Patients. Antibiotics (Basel). 2021;10(2):100. doi:10.3390/antibiotics10020100"
  vignette <- "Kato_2021_amikacin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "amikacin", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Source name CCr. Estimated by the Cockcroft-Gault equation in raw mL/min, NOT BSA-normalized (Kato 2021 Section 4.2 and Table 1 footnote). Time-fixed per subject: the covariate was taken at the initiation of treatment (Section 4.4). Enters clearance as the through-origin linear ratio CL = CLs * (CCr / 52.9) (Section 2.2 final-model equation), so clearance tends to zero as CRCL tends to zero; must be strictly positive. The 52.9 mL/min reference is the cohort MEAN (Table 1: 52.9 +/- 22.8); the Section 2.2 prose says CCr was normalized to 'the median value of the study population of 52.1 mL/min' (the Table 1 median). The printed final-model equation is carried here. Observed range 10.9-94.9 mL/min; the paper's dosing simulations span 10-90 mL/min.",
      source_name = "CCr"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 15L,
    n_studies = 1L,
    n_samples = 33L,
    age_range = "71-95 years",
    age_median = "80.0 years (mean 80.6 +/- 7.3)",
    weight_range = "32.5-67.3 kg",
    weight_median = "42.6 kg (mean 44.8 +/- 8.9)",
    sex_female_pct = 60,
    race_ethnicity = "Japanese (single Japanese university hospital); not tabulated",
    disease_state = "Elderly inpatients (>= 70 years) treated with amikacin for >= 3 days: pneumonia 5, bacteremia 5, urinary tract infection 1, urinary tract infection plus pneumonia 1, pneumonia plus bacteremia 1, peritoneum inflammation plus bacteremia 1, febrile neutropenia 1. Pseudomonas aeruginosa was the most common isolate (6 of 11 culture-positive patients). Baseline albumin 2.5 +/- 0.4 g/dL, serum creatinine 0.84 +/- 0.77 mg/dL (median 0.59).",
    renal_function = "Cockcroft-Gault creatinine clearance mean 52.9 +/- 22.8 mL/min, median 52.1 (range 10.9-94.9); patients on intermittent or continuous renal replacement therapy at the onset of amikacin therapy were excluded.",
    dose_range = "Amikacin 200-1000 mg/day IV (mean 440 +/- 226 mg/day; 9.6 +/- 3.6 mg/kg/day) by 0.5-1.0 h infusion (median 0.5 h); treatment duration 3-20 days.",
    regions = "Japan (Aichi Medical University Hospital)",
    notes = "Single-centre retrospective therapeutic drug monitoring study, September 2009 - February 2015 (Kato 2021 Section 4.1, Table 1). Two or three serum samples per patient: a trough within 30 min before a dose and a peak about 1 h after the start of infusion (30 min after the end of a 30 min infusion). Amikacin assayed by fluorescence polarization immunoassay (limit of detection 0.8 mg/L, calibration range 0.8-40 mg/L, intra- and inter-assay CV within 6%)."
  )

  ini({
    # Structural parameters: Kato 2021 Table 2, 'Mean' column of the Pmetrics
    # NPAG non-parametric population distribution. The abstract and Section
    # 2.2 quote these means as the population PK parameter estimates; Table 2
    # also reports medians, noted per line.
    lcl <- log(2.25)
    label("Clearance CLs at CRCL = 52.9 mL/min (L/h)")
    # Kato 2021 Table 2: CL mean 2.25, SD 0.78, median 2.19, SE 0.24, CV 34.6%
    lvc <- log(18.0)
    label("Volume of distribution Vd (L)")
    # Kato 2021 Table 2: V mean 18.0, SD 3.4, median 17.1, SE 1.0, CV 18.9%

    # Between-subject variability. NPAG estimates a discrete non-parametric
    # distribution, not an omega matrix; Table 2 summarises it per parameter by
    # SD and CV% (SD / mean: 0.78 / 2.25 = 34.7%, 3.4 / 18.0 = 18.9%, so the
    # CV% column is a between-subject statistic). It is carried here as a
    # log-normal random effect via omega^2 = log(CV^2 + 1), a parametric
    # APPROXIMATION of the non-parametric distribution. No correlation between
    # CL and V is reported, so the two etas are independent.
    # Table 2's 'Var' column (2.19, 17.1) repeats the median column and is not
    # SD^2 (0.61, 11.6); it is not used.
    etalcl ~ 0.113075 # Kato 2021 Table 2 CL CV 34.6%: log(0.346^2 + 1)
    etalvc ~ 0.035098 # Kato 2021 Table 2 V CV 18.9%: log(0.189^2 + 1)

    # Residual error. The Pmetrics assay-error polynomial (C0-C3) and the
    # gamma / lambda noise terms are not reported anywhere in the paper, so
    # both terms are carried as fixed(0) rather than invented -- see the
    # vignette 'Assumptions and deviations'.
    propSd <- fixed(0)
    label("Proportional residual SD (fraction; 0 -- not reported in the source)")
    addSd <- fixed(0)
    label("Additive residual SD (mg/L; 0 -- not reported in the source)")
  })

  model({
    # Reference creatinine clearance (mL/min) from the Section 2.2 final-model
    # equation 'CL = CLs x (CCr/52.9)'. 52.9 is the Table 1 cohort mean; the
    # Section 2.2 prose names the median 52.1 instead, and the printed equation
    # is carried here.
    crcl_ref <- 52.9

    # Individual parameters. Clearance is directly proportional to CRCL (no
    # exponent and no non-renal intercept in the printed equation).
    cl <- exp(lcl + etalcl) * (CRCL / crcl_ref)
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    # One-compartment model with first-order elimination; amikacin is given
    # by IV infusion into the central compartment (Kato 2021 Sections 2.2,
    # 4.3).
    d/dt(central) <- -kel * central

    # Dose in mg, volume in L -> Cc in mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
