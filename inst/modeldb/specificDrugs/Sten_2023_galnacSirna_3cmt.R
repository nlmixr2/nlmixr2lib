Sten_2023_galnacSirna_3cmt <- function() {
  description <- paste(
    "Three-compartment PK model with first-order absorption describing the plasma",
    "pharmacokinetics of GalNAc-conjugated siRNAs as a drug class after a single",
    "subcutaneous dose, fitted jointly to mean plasma concentrations of inclisiran,",
    "givosiran, lumasiran, vutrisiran (to 24-48 h) and olpasiran (to 12 weeks) with one",
    "shared parameter set (Sten 2023, Table 3, joint fit; Fig. 5). Intended for",
    "class-level prediction of a new GalNAc-siRNA, with dose in mg of siRNA.",
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
    depot = list(analyte = "GalNAc-conjugated siRNA", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "GalNAc-conjugated siRNA", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "GalNAc-conjugated siRNA", units = "mg", specimen = "not applicable", verified = TRUE),
    peripheral2 = list(analyte = "GalNAc-conjugated siRNA", units = "mg", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 5,
    disease_state = paste(
      "Pooled literature cohorts: healthy volunteers (vutrisiran, lumasiran), acute intermittent",
      "porphyria (givosiran), primary hyperoxaluria type 1 (lumasiran), adults with normal to",
      "severely impaired renal function (inclisiran) and adults with elevated lipoprotein(a)",
      "(olpasiran)",
      sep = " "
    ),
    dose_range = "5-432 mg single SC dose (olpasiran 3-225 mg)",
    regions = "Not reported by Sten 2023 (literature data)",
    notes = paste(
      "Fitted to the 34 mean concentration-time profiles of the ESM spreadsheet (sheet",
      "Time_course_data): givosiran (Agarwal 2020; FDA NDA 212194), inclisiran (Wright 2020),",
      "lumasiran (Frishberg 2021), vutrisiran (Habtemariam 2021) and olpasiran (Koren 2022).",
      "Each profile is a group mean of 3-12 subjects, so the subject count of the pooled fit is",
      "not meaningful. The fit is naive pooled with one parameter set for all five compounds;",
      "it carries no compound or population covariate.",
      sep = " "
    )
  )

  ini({
    lka <- log(0.157); label("Absorption rate constant (1/h)") # Table 3, column 'Model simultaneously fitted to all data': Ka = 0.157 (SE 0.031) 1/h
    lcl <- log(23.8); label("Apparent clearance CL1/F (L/h)") # Table 3, column 'Model simultaneously fitted to all data': CL1/F = 23.8 (SE 0.77) L/h
    lvc <- log(18.1); label("Apparent central volume V1/F (L)") # Table 3, column 'Model simultaneously fitted to all data': V1/F = 18.1 (SE 4.1) L
    lq <- log(3.02); label("Apparent intercompartmental clearance to the first peripheral compartment CL2/F (L/h)") # Table 3, column 'Model simultaneously fitted to all data': CL2/F = 3.02 (SE 0.37) L/h
    lvp <- log(1340); label("Apparent first peripheral volume V2/F (L)") # Table 3, column 'Model simultaneously fitted to all data': V2/F = 1340 (SE 180) L
    lq2 <- log(45.4); label("Apparent intercompartmental clearance to the second peripheral compartment CL3/F (L/h)") # Table 3, column 'Model simultaneously fitted to all data': CL3/F = 45.4 (SE 15) L/h
    lvp2 <- log(100); label("Apparent second peripheral volume V3/F (L)") # Table 3, column 'Model simultaneously fitted to all data': V3/F = 100 (SE 30) L
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
