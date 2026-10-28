Sten_2023_lumasiran_healthy <- function() {
  description <- paste(
    "One-compartment PK model with first-order absorption for the GalNAc-conjugated siRNA",
    "lumasiran after a single subcutaneous dose in healthy volunteers, fitted to mean plasma",
    "concentrations up to 48 h post-dose (Sten 2023, Supplementary Table S4). Naive-pooled",
    "fit of published mean profiles in Phoenix WinNonlin: no between-subject variability, and",
    "CL and V are apparent (CL/F, V/F) because only subcutaneous data were available.",
    "One of ten independent per-compound / per-population one-compartment fits in the paper.",
    sep = " "
  )
  reference <- paste(
    "Sten S, Cardilin T, Antonsson M, Gennemark P. Plasma Pharmacokinetics of",
    "N-Acetylgalactosamine-Conjugated Small-Interfering Ribonucleic Acids",
    "(GalNAc-Conjugated siRNAs). Clin Pharmacokinet. 2023;62(12):1661-1672.",
    "doi:10.1007/s40262-023-01314-7",
    sep = " "
  )
  vignette <- "Sten_2023_galnacSirna"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "lumasiran", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "lumasiran", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = 24,
    n_studies = 1,
    disease_state = "Healthy volunteers",
    dose_range = "21.6, 72, 216 and 432 mg (0.3, 1, 3 and 6 mg/kg at 72 kg) single SC dose",
    regions = "Not reported by Sten 2023 (literature data from the lumasiran Phase 1/2 study)",
    notes = paste(
      "Mean plasma concentration-time data collected by Sten 2023 from Frishberg et al. 2021 (doi:10.2215/CJN.14730920), Part A single ascending dose",
      "(ESM spreadsheet, sheet Time_course_data). n_subjects is the sum of the per-dose-arm",
      "group sizes listed in the ESM AUC_Cmax sheet for the cohorts behind these profiles.",
      "The fit is naive pooled on group means, so the model describes the group-mean profile",
      "rather than individuals.",
      sep = " "
    )
  )

  ini({
    lka <- log(0.0612); label("Absorption rate constant (1/h)") # Supplementary Table S4, 1st lumasiran column: Ka = 0.0612 (SE 0.011) 1/h
    lcl <- log(25.1); label("Apparent clearance CL/F (L/h)") # Supplementary Table S4, 1st lumasiran column: CL/F = 25.1 (SE 3.1) L/h
    lvc <- log(9.79); label("Apparent volume of distribution V/F (L)") # Supplementary Table S4, 1st lumasiran column: V/F = 9.79 (SE 4.2) L
    propSd <- fixed(0); label("Proportional residual error (fraction); not reported by Sten 2023")
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # mg/L x 1000 = ng/mL
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
