Widmer_2006_imatinib <- function() {
  description <- paste0(
    "One-compartment population PK model with first-order absorption and ",
    "first-order elimination for oral imatinib in adults with ",
    "gastrointestinal stromal tumor or chronic myeloid leukemia (Widmer ",
    "2006). This is the DEMOGRAPHIC BASE MODEL of the source paper: no ",
    "covariate is retained, so CL/F and Vc/F are typical values only. ",
    "Inter-individual variability is estimated on CL/F and Vc/F ",
    "independently; residual error is proportional. TRANSCRIBED FROM A ",
    "SECONDARY SOURCE: the parameter values come from Table 1 of Yang ",
    "2025, an external evaluation of 15 published imatinib population PK ",
    "models, not from the primary publication. Re-extract from Widmer 2006 ",
    "when that paper is obtained; in particular, the primary develops an ",
    "alpha-1 acid glycoprotein covariate model that this base model omits."
  )
  reference <- paste0(
    "Widmer N, Decosterd LA, Csajka C, Leyvraz S, Duchosal MA, Rosselet A, ",
    "Rochat B, Eap CB, Henry H, Biollaz J, Buclin T. Population ",
    "pharmacokinetics of imatinib and the role of alpha-acid glycoprotein. ",
    "Br J Clin Pharmacol. 2006;62(1):97-112. ",
    "doi:10.1111/j.1365-2125.2006.02719.x. ",
    "PARAMETER SOURCE (secondary): Yang T, Rasmussen ASB, Weimann A, ",
    "Thastrup M, Rank CU, Als-Nielsen B, Malmros J, Wik HS, Lohi O, ",
    "Overgaard U, Johannsdottir IMR, Vaitkeviciene G, Dalhoff K, ",
    "Schmiegelow K, Lund TM. Published population pharmacokinetic models ",
    "of imatinib perform poorly on TDM data from pediatric patients. ",
    "Target Oncol. 2025;20(5):871-886. Table 1, row 'Widmer et al. (2006)' ",
    "and Table 1 footnote c. doi:10.1007/s11523-025-01172-2."
  )
  vignette <- "Yang_2025_imatinib_external_evaluation"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot   = list(analyte = "imatinib", units = "mg", specimen = "administration site", verified = FALSE),
    central = list(analyte = "imatinib", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    AAG = list(
      description = "Alpha-1 acid glycoprotein concentration",
      units       = "g/L",
      type        = "continuous",
      notes       = paste0(
        "The title of the primary publication is 'Population ",
        "pharmacokinetics of imatinib and the role of alpha-acid ",
        "glycoprotein', so AAG is the central scientific subject of Widmer ",
        "2006 and the primary certainly reports an AAG model. It is NOT in ",
        "this model file because Yang 2025 Table 1 footnote c states that ",
        "the external evaluation 'used the demographic base model for ",
        "imatinib in Widmer's paper', and Yang 2025 Table 1 correspondingly ",
        "lists this row's covariates as 'None'. The secondary source ",
        "therefore transcribes only the covariate-free base model, and no ",
        "AAG coefficient is available to encode. This is the single largest ",
        "known gap between this transcription and the primary; re-extract ",
        "from Widmer 2006 to recover the AAG relationship. Compare ",
        "Petain_2008_imatinib.R and DiPaolo_2014_imatinib.R, which do carry ",
        "quantified AAG effects."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 59L,
    n_studies      = 1L,
    n_observations = "321 imatinib plasma concentrations (Yang 2025 Table 1)",
    age_range      = "20-79 years",
    disease_state  = paste0(
      "Adults with gastrointestinal stromal tumor (GIST) or chronic ",
      "myeloid leukemia (CML), plus one patient with acute lymphoblastic ",
      "leukemia (Yang 2025 Table 1)"
    ),
    dose_range     = "Oral imatinib 150-800 mg total daily dose",
    regions        = "Switzerland",
    bioanalytical  = "HPLC, limit of quantification 50 ng/mL (Yang 2025 Table 1)",
    notes          = paste0(
      "Demographic detail beyond the row above (weight range, sex split, ",
      "race) is not reported by the secondary source and must be read from ",
      "the primary publication."
    )
  )

  ini({
    # ----- Structural parameters (Yang 2025 Table 1, Widmer row) -----
    lka <- log(0.61); label("First-order absorption rate constant ka (1/h)")  # Yang 2025 Table 1: Ka = 0.61
    lcl <- log(14.3); label("Apparent oral clearance CL/F (L/h)")  # Yang 2025 Table 1: CL/F = 14.3
    lvc <- log(347); label("Apparent central volume of distribution Vc/F (L)")  # Yang 2025 Table 1: Vc/F = 347

    # ----- Inter-individual variability -----
    # Yang 2025 Table 1 reports 'CV%(CL): 36%' and 'CV%(Vc): 63%' with no
    # covariance term, so the two etas are carried as independent diagonal
    # elements. The tabulated CV% is taken as omega (the log-scale standard
    # deviation), so the variance is (CV/100)^2 -- the same convention used
    # throughout this transcription and in the shipped
    # Jiang_2023_imatinib.R.
    etalcl ~ 0.1296  # Yang 2025 Table 1: CV%(CL) 36% -> omega^2 = 0.36^2
    etalvc ~ 0.3969  # Yang 2025 Table 1: CV%(Vc) 63% -> omega^2 = 0.63^2

    # ----- Residual unexplained variability -----
    propSd <- 0.31; label("Proportional residual error (fraction)")  # Yang 2025 Table 1: Prop 31%
  })

  model({
    # ----- 1. Individual parameters -----
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)

    # ----- 2. Micro-constants -----
    kel <- cl / vc

    # ----- 3. ODE system -----
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ----- 4. Observation and error -----
    # central is mg and vc is L, so central/vc is mg/L; x 1000 gives ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
