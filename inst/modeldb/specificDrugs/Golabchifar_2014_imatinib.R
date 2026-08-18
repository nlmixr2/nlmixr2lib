Golabchifar_2014_imatinib <- function() {
  description <- paste0(
    "One-compartment population PK model with zero-order absorption of ",
    "duration D1 preceded by an absorption lag time, and first-order ",
    "elimination, for oral imatinib in Iranian adults with chronic-phase ",
    "chronic myeloid leukemia (Golabchifar 2014). No covariate is ",
    "retained. This is the only one of the 15 models evaluated by Yang ",
    "2025 that includes an absorption lag time, and the one with by far ",
    "the smallest residual error (8.1% proportional). Inter-individual ",
    "variability is estimated independently on all four structural ",
    "parameters. TRANSCRIBED FROM A SECONDARY SOURCE: the parameter values ",
    "come from Table 1 of Yang 2025, an external evaluation of 15 ",
    "published imatinib population PK models, not from the primary ",
    "publication. Re-extract from Golabchifar 2014 when that paper is ",
    "obtained."
  )
  reference <- paste0(
    "Golabchifar AA, Rezaee S, Ghavamzadeh A, Alimoghaddam K, Dinan NM, ",
    "Rouini MR. Population pharmacokinetics of imatinib in Iranian ",
    "patients with chronic-phase chronic myeloid leukemia. Cancer ",
    "Chemother Pharmacol. 2014;74(1):85-93. ",
    "doi:10.1007/s00280-014-2455-3. ",
    "PARAMETER SOURCE (secondary): Yang T, Rasmussen ASB, Weimann A, ",
    "Thastrup M, Rank CU, Als-Nielsen B, Malmros J, Wik HS, Lohi O, ",
    "Overgaard U, Johannsdottir IMR, Vaitkeviciene G, Dalhoff K, ",
    "Schmiegelow K, Lund TM. Published population pharmacokinetic models ",
    "of imatinib perform poorly on TDM data from pediatric patients. ",
    "Target Oncol. 2025;20(5):871-886. Table 1, row 'Golabchifar et al. ",
    "(2014)'. doi:10.1007/s11523-025-01172-2."
  )
  vignette <- "Yang_2025_imatinib_external_evaluation"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    central = list(analyte = "imatinib", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 61L,
    n_studies      = 1L,
    n_observations = "533 imatinib plasma concentrations (Yang 2025 Table 1)",
    age_range      = "17-67 years",
    disease_state  = "Adults with chronic-phase chronic myeloid leukemia (CML)",
    dose_range     = "Oral imatinib 300-800 mg total daily dose",
    regions        = "Iran",
    bioanalytical  = "HPLC, limit of quantification 62.5 ng/mL (Yang 2025 Table 1)",
    notes          = paste0(
      "With standard allometric scaling applied by Yang 2025, this model ",
      "was one of the two best overall performers on their external ",
      "dataset, meeting the bias criterion (median prediction error ",
      "-7.68%) with the lowest median absolute prediction error of any ",
      "scaled model except Widmer (39.58%) and the highest fraction of ",
      "predictions within +/- 20% (0.32) -- see Yang 2025 Table 3 and the ",
      "Abstract. The allometric variant is Yang 2025's own modification ",
      "and is NOT encoded here; this file is the model as published. ",
      "Demographic detail beyond the row above (weight range, sex split) ",
      "is not reported by the secondary source and must be read from the ",
      "primary publication."
    )
  )

  ini({
    # ----- Structural parameters (Yang 2025 Table 1, Golabchifar row) ----
    lcl <- log(10.8); label("Apparent oral clearance CL/F (L/h)")  # Yang 2025 Table 1: CL/F = 10.8
    lvc <- log(278); label("Apparent central volume of distribution Vc/F (L)")  # Yang 2025 Table 1: Vc/F = 278
    ld1 <- log(1.43); label("Zero-order absorption duration D1 (h)")  # Yang 2025 Table 1: D1 = 1.43
    ltlag <- log(0.197); label("Absorption lag time Tlag (h)")  # Yang 2025 Table 1: Tlag = 0.197

    # ----- Inter-individual variability -----
    # Yang 2025 Table 1 reports 'CV%(CL): 30.2%', 'CV%(Vc): 53.5%',
    # 'CV%(D1): 36%' and 'CV%(Tlag): 68.3%'. No covariances are tabulated,
    # so the four etas are carried as independent diagonal elements. The
    # tabulated CV% is taken as omega (the log-scale standard deviation),
    # so the variance is (CV/100)^2 -- the convention used throughout this
    # transcription and in the shipped Jiang_2023_imatinib.R. This is the
    # only one of the 15 models with IIV on all four of its structural
    # parameters.
    etalcl ~ 0.091204  # Yang 2025 Table 1: CV%(CL) 30.2% -> omega^2 = 0.302^2
    etalvc ~ 0.286225  # Yang 2025 Table 1: CV%(Vc) 53.5% -> omega^2 = 0.535^2
    etald1 ~ 0.1296  # Yang 2025 Table 1: CV%(D1) 36% -> omega^2 = 0.36^2
    etaltlag ~ 0.466489  # Yang 2025 Table 1: CV%(Tlag) 68.3% -> omega^2 = 0.683^2

    # ----- Residual unexplained variability -----
    # 8.1% is by a wide margin the smallest residual error among the 15
    # models evaluated by Yang 2025 (the next smallest is He 2023 at
    # 18.5%), which is consistent with a densely sampled single-centre
    # study rather than routine therapeutic drug monitoring.
    propSd <- 0.081; label("Proportional residual error (fraction)")  # Yang 2025 Table 1: Prop 8.1%
  })

  model({
    # ----- 1. Individual parameters -----
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    d1 <- exp(ld1 + etald1)
    tlag <- exp(ltlag + etaltlag)

    # ----- 2. Micro-constants -----
    kel <- cl / vc

    # ----- 3. ODE system -----
    # One compartment, zero-order input starting after the lag time. Dose
    # records must carry rate = -2 so rxode2 uses the modelled duration
    # dur(central) = d1. Mean absorption time is tlag + d1 / 2 = 0.912 h.
    d/dt(central) <- -kel * central
    dur(central) <- d1
    alag(central) <- tlag

    # ----- 4. Observation and error -----
    # central is mg and vc is L, so central/vc is mg/L; x 1000 gives ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
