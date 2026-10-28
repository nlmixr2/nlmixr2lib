Sten_2023_olpasiran_3cmt <- function() {
  description <- paste(
    "Three-compartment PK model with first-order absorption for the GalNAc-conjugated siRNA",
    "olpasiran after a single subcutaneous dose in adults with elevated lipoprotein(a),",
    "fitted to mean plasma concentrations up to 12 weeks post-dose (Sten 2023, Table 3,",
    "olpasiran-only fit). Sten 2023 used it to predict the 24-48 h profiles of inclisiran,",
    "givosiran, lumasiran and vutrisiran (ESM Fig. S1).",
    "Naive-pooled fit of published mean profiles in Phoenix WinNonlin: no between-subject",
    "variability, and clearances and volumes are apparent (divided by F) because only",
    "subcutaneous data were available.",
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
    central = list(analyte = "olpasiran", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "olpasiran", units = "mg", specimen = "not applicable", verified = TRUE),
    peripheral2 = list(analyte = "olpasiran", units = "mg", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = 48,
    n_studies = 1,
    disease_state = "Adults with elevated lipoprotein(a) concentration",
    dose_range = "3, 9, 30, 75 and 225 mg single SC dose",
    regions = "Not reported by Sten 2023 (literature data from the olpasiran Phase 1 study)",
    notes = paste(
      "Mean plasma concentration-time data (to 2016 h) collected by Sten 2023 from Koren et al. 2022",
      "(doi:10.1038/s41591-021-01634-w) single ascending dose study (ESM spreadsheet, sheet",
      "Time_course_data). n_subjects is the sum of the per-dose-arm group sizes in the ESM",
      "AUC_Cmax sheet for the seven olpasiran profiles. The fit is naive pooled on group means,",
      "so the model describes the group-mean profile rather than individuals.",
      sep = " "
    )
  )

  ini({
    lka <- log(0.635); label("Absorption rate constant (1/h)") # Table 3, column 'Model fitted to olpasiran data only': Ka = 0.635 (SE 0.21) 1/h
    lcl <- log(27.7); label("Apparent clearance CL1/F (L/h)") # Table 3, column 'Model fitted to olpasiran data only': CL1/F = 27.7 (SE 2.4) L/h
    lvc <- log(27.9); label("Apparent central volume V1/F (L)") # Table 3, column 'Model fitted to olpasiran data only': V1/F = 27.9 (SE 5.6) L
    lq <- log(2.62); label("Apparent intercompartmental clearance to the first peripheral compartment CL2/F (L/h)") # Table 3, column 'Model fitted to olpasiran data only': CL2/F = 2.62 (SE 1.4) L/h
    lvp <- log(1470); label("Apparent first peripheral volume V2/F (L)") # Table 3, column 'Model fitted to olpasiran data only': V2/F = 1470 (SE 370) L
    lq2 <- log(914); label("Apparent intercompartmental clearance to the second peripheral compartment CL3/F (L/h)") # Table 3, column 'Model fitted to olpasiran data only': CL3/F = 914 (SE 560) L/h
    lvp2 <- log(296); label("Apparent second peripheral volume V3/F (L)") # Table 3, column 'Model fitted to olpasiran data only': V3/F = 296 (SE 39) L
    propSd <- fixed(0); label("Proportional residual error (fraction); not reported by Sten 2023")
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)
    q2 <- exp(lq2)
    vp2 <- exp(lvp2)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # mg/L x 1000 = ng/mL
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
