Shriyan_2022_imatinib <- function() {
  description <- paste0(
    "One-compartment population PK model with zero-order absorption of ",
    "duration D1 into the central compartment and first-order elimination ",
    "for oral imatinib in Indian adults with chronic myeloid leukemia ",
    "(Shriyan 2022). No covariate is retained in the final model even ",
    "though the source paper is an ADME-gene pharmacogenetic study; ",
    "CL/F and Vc/F are typical values only. Inter-individual variability ",
    "is estimated on CL/F only; residual error is proportional. ",
    "TRANSCRIBED FROM A SECONDARY SOURCE: the parameter values come from ",
    "Table 1 of Yang 2025, an external evaluation of 15 published imatinib ",
    "population PK models, not from the primary publication. Re-extract ",
    "from Shriyan 2022 when that paper is obtained."
  )
  reference <- paste0(
    "Shriyan B, Mehta P, Patil A, Jadhav S, Kumar S, Puri AS, Bagal B, ",
    "Sengar M, Khattry N, Gota V. Role of ADME gene polymorphisms on ",
    "imatinib disposition: results from a population pharmacokinetic study ",
    "in chronic myeloid leukaemia. Eur J Clin Pharmacol. ",
    "2022;78(8):1321-1330. doi:10.1007/s00228-022-03334-x. ",
    "PARAMETER SOURCE (secondary): Yang T, Rasmussen ASB, Weimann A, ",
    "Thastrup M, Rank CU, Als-Nielsen B, Malmros J, Wik HS, Lohi O, ",
    "Overgaard U, Johannsdottir IMR, Vaitkeviciene G, Dalhoff K, ",
    "Schmiegelow K, Lund TM. Published population pharmacokinetic models ",
    "of imatinib perform poorly on TDM data from pediatric patients. ",
    "Target Oncol. 2025;20(5):871-886. Table 1, row 'Shriyan et al. ",
    "(2022)'. doi:10.1007/s11523-025-01172-2."
  )
  vignette <- "Yang_2025_imatinib_external_evaluation"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    central = list(analyte = "imatinib", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    WBC = list(
      description = "White blood cell count",
      units       = "10^9 cells/L",
      type        = "continuous",
      notes       = paste0(
        "Not retained in this model (Yang 2025 Table 1 lists this row's ",
        "covariates as 'None'). Recorded here because the shipped ",
        "Jiang_2023_imatinib.R Discussion cites Shriyan by name on exactly ",
        "this point: a prior CML study found a white-blood-cell-count ",
        "effect on imatinib clearance, but 'no such correlation was found ",
        "in the population described in this study, which is consistent ",
        "with the results obtained by Delbaldo et al and Shriyan et al'. ",
        "Contrast Schmidli_2005_imatinib.R and Demetri_2009_imatinib.R, ",
        "which do retain a WBC effect."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 49L,
    n_studies      = 1L,
    n_observations = "221 imatinib plasma concentrations (Yang 2025 Table 1)",
    age_range      = "18 years and older",
    disease_state  = "Indian adults with chronic myeloid leukemia (CML)",
    dose_range     = "Oral imatinib 400-800 mg total daily dose",
    regions        = "India",
    bioanalytical  = "HPLC; limit of quantification not reported in Yang 2025 Table 1",
    notes          = paste0(
      "With standard allometric scaling applied by Yang 2025, this model ",
      "produced the SMALLEST bias of any of the 30 model variants they ",
      "evaluated: a median prediction error of 1.24% (Yang 2025 Abstract ",
      "and Table 3). Its precision was nonetheless poor (median absolute ",
      "prediction error 42.71%), which is the paper's central finding -- ",
      "acceptable average bias does not imply usable individual ",
      "predictions. The allometric variant is Yang 2025's own modification ",
      "and is NOT encoded here; this file is the model as published. ",
      "Demographic detail beyond the row above (weight range, sex split) ",
      "is not reported by the secondary source and must be read from the ",
      "primary publication."
    )
  )

  ini({
    # ----- Structural parameters (Yang 2025 Table 1, Shriyan row) -----
    lcl <- log(10.2); label("Apparent oral clearance CL/F (L/h)")  # Yang 2025 Table 1: CL/F = 10.2
    lvc <- log(389); label("Apparent central volume of distribution Vc/F (L)")  # Yang 2025 Table 1: Vc/F = 389
    ld1 <- log(2.42); label("Zero-order absorption duration D1 (h)")  # Yang 2025 Table 1: D1 = 2.42

    # ----- Inter-individual variability -----
    # Yang 2025 Table 1 reports only 'CV%(CL): 27.8%' for this row: no
    # random effect on Vc/F or D1 is tabulated. The tabulated CV% is taken
    # as omega (the log-scale standard deviation), so the variance is
    # (CV/100)^2 -- the convention used throughout this transcription and
    # in the shipped Jiang_2023_imatinib.R.
    etalcl ~ 0.077284  # Yang 2025 Table 1: CV%(CL) 27.8% -> omega^2 = 0.278^2

    # ----- Residual unexplained variability -----
    propSd <- 0.21; label("Proportional residual error (fraction)")  # Yang 2025 Table 1: Prop 21%
  })

  model({
    # ----- 1. Individual parameters -----
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    d1 <- exp(ld1)

    # ----- 2. Micro-constants -----
    kel <- cl / vc

    # ----- 3. ODE system -----
    # One compartment, zero-order input. Dose records must carry rate = -2
    # so rxode2 uses the modelled duration dur(central) = d1.
    d/dt(central) <- -kel * central
    dur(central) <- d1

    # ----- 4. Observation and error -----
    # central is mg and vc is L, so central/vc is mg/L; x 1000 gives ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
