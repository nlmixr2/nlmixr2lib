Wang_2019_imatinib <- function() {
  description <- paste0(
    "One-compartment population PK model with first-order absorption and ",
    "first-order elimination for oral imatinib 400 mg once daily in ",
    "Chinese adults with chronic myeloid leukemia (Wang 2019). CL/F ",
    "carries a power-form total body weight effect referenced to 70 kg ",
    "with an exponent of 0.228, well below the allometric 0.75. ",
    "Inter-individual variability is estimated on CL/F only; residual ",
    "error is combined proportional plus a small additive term. ",
    "TRANSCRIBED FROM A SECONDARY SOURCE: the parameter values come from ",
    "Table 1 of Yang 2025, an external evaluation of 15 published imatinib ",
    "population PK models, not from the primary publication. Re-extract ",
    "from Wang 2019 when that paper is obtained; the primary is a ",
    "pharmacogenetic analysis, so it may report genotype effects that the ",
    "secondary source's covariate cell ('TBW') does not carry."
  )
  reference <- paste0(
    "Wang Q, Jiang ZP, Yu EQ, Zeng J, Zhu Y, Cai HL, Zhang M, Zhang BK, ",
    "Xiang DX. Population pharmacokinetic and pharmacogenetics of imatinib ",
    "in Chinese patients with chronic myeloid leukemia. Pharmacogenomics. ",
    "2019;20(4):251-260. doi:10.2217/pgs-2018-0173. ",
    "PARAMETER SOURCE (secondary): Yang T, Rasmussen ASB, Weimann A, ",
    "Thastrup M, Rank CU, Als-Nielsen B, Malmros J, Wik HS, Lohi O, ",
    "Overgaard U, Johannsdottir IMR, Vaitkeviciene G, Dalhoff K, ",
    "Schmiegelow K, Lund TM. Published population pharmacokinetic models ",
    "of imatinib perform poorly on TDM data from pediatric patients. ",
    "Target Oncol. 2025;20(5):871-886. Table 1, row 'Wang et al. (2019)'. ",
    "doi:10.1007/s11523-025-01172-2."
  )
  vignette <- "Yang_2025_imatinib_external_evaluation"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot   = list(analyte = "imatinib", units = "mg", specimen = "administration site", verified = FALSE),
    central = list(analyte = "imatinib", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Enters CL/F as the power function (WT/70)^0.228 (Yang 2025 Table ",
        "1). The reference 70 kg is the centring constant printed inside ",
        "the covariate term. The exponent is well below the allometric ",
        "0.75, so the model's weight effect is much shallower than ",
        "allometric theory predicts: doubling weight raises CL/F by only ",
        "2^0.228 - 1 = 17%, against the 68% an allometric exponent would ",
        "give. That shallowness is why Yang 2025 nonetheless evaluated a ",
        "variant of this model in which STANDARD ALLOMETRIC SCALING ",
        "replaced the published weight term when predicting in a cohort ",
        "containing children; that variant is Yang 2025's own modification ",
        "and is NOT encoded here. In its original form this model was one ",
        "of only four whose median prediction error met Yang 2025's ",
        "+/- 15% bias criterion (-12.65%, Table 3)."
      ),
      source_name        = "TBW"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 170L,
    n_studies      = 1L,
    n_observations = "229 imatinib plasma concentrations (Yang 2025 Table 1)",
    age_range      = "16-82 years",
    disease_state  = "Chinese adults with chronic myeloid leukemia (CML)",
    dose_range     = "Oral imatinib 400 mg total daily dose",
    regions        = "China",
    bioanalytical  = "UPLC-MS/MS, limit of quantification 2.6 ng/mL (Yang 2025 Table 1)",
    notes          = paste0(
      "Sparse sampling: 229 samples from 170 patients, under 1.5 per ",
      "subject, which is consistent with the absence of any random effect ",
      "on Vc/F. Demographic detail beyond the row above (weight range, ",
      "sex split) is not reported by the secondary source and must be read ",
      "from the primary publication."
    )
  )

  ini({
    # ----- Structural parameters (Yang 2025 Table 1, Wang row) -----
    # The typical value is for a 70 kg subject, at which the weight factor
    # equals 1.
    lka <- log(0.329); label("First-order absorption rate constant ka (1/h)")  # Yang 2025 Table 1: Ka = 0.329
    lcl <- log(9.25); label("Apparent oral clearance CL/F at WT = 70 kg (L/h)")  # Yang 2025 Table 1: CL/F = 9.25 x (TBW/70)^0.228
    lvc <- log(222); label("Apparent central volume of distribution Vc/F (L)")  # Yang 2025 Table 1: Vc/F = 222

    # ----- Covariate effect on CL/F -----
    e_wt_cl <- 0.228; label("Power exponent of (WT / 70) on CL/F (unitless)")  # Yang 2025 Table 1: (TBW/70)^0.228

    # ----- Inter-individual variability -----
    # Yang 2025 Table 1 reports only 'CV%(CL): 18.2%' for this row: no
    # random effect on Vc/F or ka is tabulated. The tabulated CV% is taken
    # as omega (the log-scale standard deviation), so the variance is
    # (CV/100)^2 -- the convention used throughout this transcription and
    # in the shipped Jiang_2023_imatinib.R.
    etalcl ~ 0.033124  # Yang 2025 Table 1: CV%(CL) 18.2% -> omega^2 = 0.182^2

    # ----- Residual unexplained variability -----
    # The additive term is tabulated in ng/mL, matching the units of Cc, so
    # no conversion is needed. It is small relative to the assay limit of
    # quantification (2.6 ng/mL) and negligible against the 40.6%
    # proportional term, but is retained so the published error model is
    # reproduced exactly.
    propSd <- 0.406; label("Proportional residual error (fraction)")  # Yang 2025 Table 1: Prop 40.6%
    addSd <- 3.137; label("Additive residual error (ng/mL)")  # Yang 2025 Table 1: Add 3.137 ng/mL
  })

  model({
    # ----- 1. Individual parameters -----
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc)

    # ----- 2. Micro-constants -----
    kel <- cl / vc

    # ----- 3. ODE system -----
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ----- 4. Observation and error -----
    # central is mg and vc is L, so central/vc is mg/L; x 1000 gives ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
