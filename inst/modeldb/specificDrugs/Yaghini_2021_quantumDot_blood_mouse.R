Yaghini_2021_quantumDot_blood_mouse <- function() {
  description <- paste(
    "Preclinical (mouse, female Balb/c). Two-compartment population PK model with",
    "first-order elimination for indium-based, cadmium-free PEGylated quantum dots",
    "(bio CFQD, InP core) after a single IV bolus, fitted to serum indium",
    "concentrations over the first 24 h. No between-animal variability (destructive",
    "and sparse sampling). The printed parameters do not reproduce the paper's own",
    "Figure 5 profile; see the vignette Errata. Companion organ-disposition model:",
    "modellib('Yaghini_2021_quantumDot_organ_mouse')."
  )
  reference <- paste(
    "Yaghini E, Tacconi E, Pilling A, Rahman P, Broughton J, Naasani I, Keshtgar MRS,",
    "MacRobert A, Della Pasqua O. (2021). Population pharmacokinetic modelling of",
    "indium-based quantum dot nanoparticles: preclinical in vivo studies.",
    "Eur J Pharm Sci 157:105639. doi:10.1016/j.ejps.2020.105639."
  )
  vignette <- "Yaghini_2021_quantumDot"
  units <- list(time = "h", dosing = "ug", concentration = "ug/L")

  covariateData <- list()

  compartmentData <- list(
    central = list(analyte = "indium (quantum dot)", units = "ug", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "mouse (female Balb/c)",
    n_subjects = 9L,
    n_studies = 1L,
    sex_female_pct = 100,
    disease_state = "healthy",
    dose_range = "single IV bolus of 20 mg/kg PEGylated bio CFQD quantum dots (200 uL in PBS via the tail vein)",
    notes = paste(
      "The analysis dataset is not tabulated. Figure 5 shows 21 serum indium",
      "observations at 7 times (5, 20, 40, 60 and 120 min from one set of 3 serially",
      "sampled mice; 4 h and 24 h from two biodistribution groups of 3), so 9 mice is",
      "inferred from the figure. Indium was quantified by ICP-MS in serum; the indium",
      "mass fraction of the quantum dot is not reported, so the dose units that make",
      "the model dimensionally consistent are not stated (see vignette Errata).",
      "Body weight is not reported."
    )
  )

  ini({
    lcl <- log(0.0416); label("Clearance from the central compartment (L/h)") # Table 1: CL = 0.0416 L/h (12% CV)
    lvc <- log(0.0027); label("Central volume of distribution (L)") # Table 1: V1 = 0.0027 L (9% CV)
    lq <- log(0.498); label("Intercompartmental clearance (L/h)") # Table 1: Q = 0.498 L/h (21% CV)
    lvp <- log(0.068); label("Peripheral volume of distribution (L)") # Table 1: V2 = 0.068 L (17% CV)

    # Table 1 prints the residual rows without a scale. Methods eq (2) defines
    # epsilon as having variance sigma^2, so the printed values are read as
    # NONMEM $SIGMA variances and converted to SDs here. The additive row is
    # printed in ng/L; the model works in ug/L, hence the /1000.
    propSd <- 0.2619160; label("Proportional residual error (fraction)") # Table 1: proportional error 0.0686 (8% CV); sqrt(0.0686)
    addSd <- 0.0002236068; label("Additive residual error (ug/L)") # Table 1: additive error 0.05 ng/L (38%), 12% CV; sqrt(0.05)/1000
  })
  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)

    Cc <- linCmt()
    Cc ~ add(addSd) + prop(propSd)
  })
}
