Tosca_2023_men1611_mouse <- function() {
  description <- paste(
    "Preclinical (mouse; tumor-free BALB-nu/nu and CD-1 mice). Two-compartment",
    "PK model with first-order oral absorption for the PI3K inhibitor MEN1611",
    "after IV and oral administration, with separate oral bioavailabilities for",
    "the MCT suspension and the DMSO solution (fixed to 1). Naive-pooled fit to",
    "mean concentration-time profiles with a combined additive + proportional",
    "residual error and no between-animal variability."
  )
  reference <- paste(
    "Tosca EM, Borella E, Piana C, Bouchene S, Merlino G, Fiascarelli A,",
    "Mazzei P, Magni P. Model-based prediction of",
    "effective target exposure for MEN1611 in combination with trastuzumab in",
    "HER2-positive advanced or metastatic breast cancer patients.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(11):1626-1639.",
    "doi:10.1002/psp4.12910. Parameter values from Supplementary Table S4.",
    sep = " "
  )
  vignette <- "Tosca_2023_men1611_trastuzumab"

  units <- list(time = "h", dosing = "mg/kg", concentration = "ng/mL")

  covariateData <- list(
    FORM_SOLUTION = list(
      description = "MEN1611 oral formulation indicator: 1 = solution in DMSO (MEN1611-PK study 1), 0 = suspension in MCT (MEN1611-PK study 2)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (MCT suspension)",
      notes = "Selects the oral bioavailability: fixed to 1 for the DMSO solution, estimated at 0.559 for the MCT suspension. IV doses go straight to central and are unaffected.",
      source_name = "F_solution / F_suspension (Table S4)"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "MEN1611", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "MEN1611", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "MEN1611", units = "mg/kg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "mouse (tumor-free BALB-nu/nu and CD-1)",
    n_subjects = 54L,
    n_studies = 2L,
    disease_state = "tumor-free mice (MEN1611 PK in tumor-bearing mice was judged poorly informative and not used)",
    dose_range = "Study 1: 2 mg/kg IV and 0.5, 1, 2 or 5 mg/kg oral (DMSO solution), single dose. Study 2: 1.25, 6.25 or 12.5 mg/kg oral (MCT suspension) once daily for 5 days.",
    regions = "preclinical",
    notes = paste(
      "MEN1611-PK study 1: BALB-nu/nu mice, 2 animals for IV and 4 animals per",
      "oral group; blood sampled over 24 h. MEN1611-PK study 2: CD-1 mice, 12",
      "animals per group, sampled over 24 h on days 1 and 5. Mean profiles",
      "stratified by study, dose and day were fitted by a naive-pooled approach",
      "in Monolix 2016R1 (Tosca 2023 Methods)."
    )
  )

  ini({
    # Tosca 2023 Supplementary Table S4, 'Estimates of MEN1611 PK model
    # parameters in TF mice'. The table prints CL and Q as 'L/h' beside
    # volumes in 'L/kg'; they are per-kg clearances (L/h/kg) -- see the
    # vignette, where this reading reproduces the paper's own MEN1611 peak
    # concentrations (Table S9) to 0.1%.
    lka <- log(1.84); label("Absorption rate constant (1/h)") # Table S4: Kabs = 1.84 1/h (8% CV)
    lcl <- log(1); label("Clearance (L/h/kg)") # Table S4: CL = 1 (4% CV)
    lvc <- log(0.997); label("Central volume of distribution (L/kg)") # Table S4: V1 = 0.997 L/kg (9% CV)
    lq <- log(1.49); label("Intercompartmental clearance (L/h/kg)") # Table S4: Q = 1.49 (13% CV)
    lvp <- log(7.6); label("Peripheral volume of distribution (L/kg)") # Table S4: V2 = 7.6 L/kg (12% CV)
    lfdepot <- log(0.559); label("Oral bioavailability of the MCT suspension (fraction)") # Table S4: F suspension = 0.559 (9% CV)
    lfdepot_sol <- fixed(log(1)); label("Oral bioavailability of the DMSO solution (fraction)") # Table S4: F solution = 1 (fixed)

    # Combined residual error y = f + (a + b*f)*eps (Tosca 2023 Methods,
    # 'Data analysis and software').
    addSd <- 1.66; label("Additive residual error (ng/mL)") # Table S4: a = 1.66 (5% CV)
    propSd <- 0.476; label("Proportional residual error (fraction)") # Table S4: b = 0.476 (9% CV)
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- exp(lfdepot * (1 - FORM_SOLUTION) + lfdepot_sol * FORM_SOLUTION)

    # central is in mg/kg and vc in L/kg, so central / vc is mg/L; x 1000 gives ng/mL.
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
