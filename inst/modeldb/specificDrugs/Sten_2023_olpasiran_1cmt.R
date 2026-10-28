Sten_2023_olpasiran_1cmt <- function() {
  description <- paste(
    "One-compartment PK model with first-order absorption for the GalNAc-conjugated siRNA",
    "olpasiran after a single subcutaneous dose in adults with elevated lipoprotein(a), fitted to mean plasma",
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
    depot = list(analyte = "olpasiran", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "olpasiran", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = 48,
    n_studies = 1,
    disease_state = "Elevated lipoprotein(a) concentration",
    dose_range = "3, 9, 30, 75 and 225 mg single SC dose",
    regions = "Not reported by Sten 2023 (literature data from the olpasiran Phase 1 study)",
    notes = paste(
      "Mean plasma concentration-time data collected by Sten 2023 from Koren et al. 2022 (doi:10.1038/s41591-021-01634-w), single ascending dose",
      "(ESM spreadsheet, sheet Time_course_data). n_subjects is the sum of the per-dose-arm",
      "group sizes listed in the ESM AUC_Cmax sheet for the cohorts behind these profiles.",
      "The fit is naive pooled on group means, so the model describes the group-mean profile",
      "rather than individuals.",
      sep = " "
    )
  )

  ini({
    lka <- log(0.0922); label("Absorption rate constant (1/h)") # Supplementary Table S4, olpasiran column: Ka = 0.0922 (SE 0.0064) 1/h
    lcl <- log(31.4); label("Apparent clearance CL/F (L/h)") # Supplementary Table S4, olpasiran column: CL/F = 31.4 (SE 2.2) L/h
    lvc <- log(11.8); label("Apparent volume of distribution V/F (L)") # Supplementary Table S4, olpasiran column: V/F = 11.8 (SE 2.4) L
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
