He_2023_imatinib <- function() {
  description <- paste0(
    "One-compartment population PK model with first-order absorption and ",
    "first-order elimination for oral imatinib in Chinese adults with ",
    "chronic myeloid leukemia in a real-world setting (He 2023). CL/F ",
    "carries a power-form hemoglobin effect referenced to 13 g/dL and a ",
    "one-sided renal-function effect that reduces clearance only in ",
    "subjects whose CKD-EPI estimated glomerular filtration rate falls ",
    "below 85 mL/min/1.73 m^2. Inter-individual variability is estimated ",
    "on CL/F only; residual error is proportional. TRANSCRIBED FROM A ",
    "SECONDARY SOURCE: the parameter values come from Table 1 of Yang ",
    "2025, an external evaluation of 15 published imatinib population PK ",
    "models, not from the primary publication. Re-extract from He 2023 ",
    "when that paper is obtained."
  )
  reference <- paste0(
    "He S, Shao Q, Zhao J, Bian J, Zhao Y, Hao X, Li Y, Wang L, Cui C, ",
    "Chen J. Population pharmacokinetics and pharmacogenetics analyses of ",
    "imatinib in Chinese patients with chronic myeloid leukemia in a ",
    "real-world situation. Cancer Chemother Pharmacol. ",
    "2023;92(5):399-410. doi:10.1007/s00280-023-04581-0. ",
    "PARAMETER SOURCE (secondary): Yang T, Rasmussen ASB, Weimann A, ",
    "Thastrup M, Rank CU, Als-Nielsen B, Malmros J, Wik HS, Lohi O, ",
    "Overgaard U, Johannsdottir IMR, Vaitkeviciene G, Dalhoff K, ",
    "Schmiegelow K, Lund TM. Published population pharmacokinetic models ",
    "of imatinib perform poorly on TDM data from pediatric patients. ",
    "Target Oncol. 2025;20(5):871-886. Table 1, row 'He et al. (2023)' ",
    "and Table 1 footnote h. doi:10.1007/s11523-025-01172-2."
  )
  vignette <- "Yang_2025_imatinib_external_evaluation"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot   = list(analyte = "imatinib", units = "mg", specimen = "administration site", verified = FALSE),
    central = list(analyte = "imatinib", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    HGB = list(
      description        = "Hemoglobin concentration",
      units              = "g/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Enters CL/F as the power function (HGB/13)^0.671 (Yang 2025 Table ",
        "1). Yang 2025 Table 1 abbreviation list gives the unit ",
        "explicitly: 'HG hemoglobin (g/dL)'. The reference 13 g/dL is the ",
        "centring constant printed inside the covariate term, matching the ",
        "reference used by Schmidli_2005_imatinib.R. Data assemblers ",
        "holding hemoglobin in g/L must divide by 10 before using this ",
        "model; the register permits either unit for HGB and requires the ",
        "per-model unit to be declared here, which it is. The positive ",
        "exponent means anemic patients have lower apparent clearance and ",
        "therefore higher exposure on the same dose -- the same direction ",
        "as the red-blood-cell-count effect in the shipped ",
        "Jiang_2023_imatinib.R ((RBC/3.7)^0.49), which is the expected ",
        "relationship given that imatinib partitions extensively into ",
        "erythrocytes."
      ),
      source_name        = "HG"
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate (CKD-EPI equation)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Enters CL/F as a ONE-SIDED (hinged) power term. Yang 2025 Table 1 ",
        "footnote h gives it as 'When eGFR < 85 the theta_eGFR = ",
        "(eGFR/85)^0.25'. The footnote is truncated in the published ",
        "table -- no 'else' branch is printed -- so the value above the ",
        "hinge is taken to be 1, which is the only reading that makes the ",
        "function continuous at eGFR = 85 and is the standard form for a ",
        "renal-impairment-only covariate. That assumption is recorded in ",
        "the vignette's Assumptions and deviations section and should be ",
        "confirmed against the primary. Below the hinge the exponent 0.25 ",
        "makes the effect shallow: a patient at eGFR 40 has a factor of ",
        "(40/85)^0.25 = 0.827, i.e. 17% lower apparent clearance. ",
        "Yang 2025 Table 1 abbreviation list specifies 'eGFR estimated ",
        "glomerular filtration rate calculated by CKD-EPI equation ",
        "(mL/min/1.73 m^2)', which matches the canonical CRCL unit exactly ",
        "-- the register pools creatinine-based estimates and ",
        "tracer-measured GFR onto CRCL because they play the same ",
        "operational role -- so no conversion is required. Note that a ",
        "renal covariate on imatinib clearance is mechanistically ",
        "surprising, since imatinib and its metabolites are barely ",
        "excreted renally (a point the shipped Jiang_2023_imatinib.R ",
        "records from its own primary); the effect is most likely a marker ",
        "of general physiological reserve rather than a renal-elimination ",
        "pathway."
      ),
      source_name        = "eGFR"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 230L,
    n_studies      = 1L,
    n_observations = "424 imatinib plasma concentrations (Yang 2025 Table 1)",
    age_range      = "14-87 years",
    disease_state  = "Chinese adults with chronic myeloid leukemia (CML) in a real-world setting",
    dose_range     = "Oral imatinib 200, 300 or 400 mg total daily dose",
    regions        = "China",
    bioanalytical  = "UPLC-MS, limit of quantification 5 ng/mL (Yang 2025 Table 1)",
    notes          = paste0(
      "The second-largest cohort among the 15 models evaluated by Yang ",
      "2025 after Gotta 2014. Its original (unscaled) form was one of only ",
      "four meeting Yang 2025's bias criterion of a median prediction ",
      "error within +/- 15% (6.42%, Table 3), but it had the WORST ",
      "precision of all 15 original models (median absolute prediction ",
      "error 58.4%). Note the age range starts at 14 years, so the cohort ",
      "includes a small number of adolescents. Demographic detail beyond ",
      "the row above (weight range, sex split) is not reported by the ",
      "secondary source and must be read from the primary publication."
    )
  )

  ini({
    # ----- Structural parameters (Yang 2025 Table 1, He row) -----
    # Typical values are for the REFERENCE subject: HGB = 13 g/dL and
    # eGFR at or above 85 mL/min/1.73 m^2, at which both covariate factors
    # equal 1.
    lka <- log(0.329); label("First-order absorption rate constant ka (1/h)")  # Yang 2025 Table 1: Ka = 0.329
    lcl <- log(7.6); label("Apparent oral clearance CL/F at the reference subject (L/h)")  # Yang 2025 Table 1: CL/F = 7.6 x (HG/13)^0.671 x theta_eGFR
    lvc <- log(270); label("Apparent central volume of distribution Vc/F (L)")  # Yang 2025 Table 1: Vc/F = 270

    # ----- Covariate effects on CL/F -----
    e_hgb_cl <- 0.671; label("Power exponent of (HGB / 13) on CL/F (unitless)")  # Yang 2025 Table 1: (HG/13)^0.671
    e_crcl_cl <- 0.25; label("Power exponent of (CRCL / 85) on CL/F, applied only below 85 mL/min/1.73 m^2 (unitless)")  # Yang 2025 Table 1 footnote h: theta_eGFR = (eGFR/85)^0.25 when eGFR < 85

    # ----- Inter-individual variability -----
    # Yang 2025 Table 1 reports only 'CV%(CL): 23.6%' for this row: no
    # random effect on Vc/F or ka is tabulated. The tabulated CV% is taken
    # as omega (the log-scale standard deviation), so the variance is
    # (CV/100)^2 -- the convention used throughout this transcription and
    # in the shipped Jiang_2023_imatinib.R.
    etalcl ~ 0.055696  # Yang 2025 Table 1: CV%(CL) 23.6% -> omega^2 = 0.236^2

    # ----- Residual unexplained variability -----
    propSd <- 0.185; label("Proportional residual error (fraction)")  # Yang 2025 Table 1: Prop 18.5%
  })

  model({
    # ----- 1. Covariate factors on CL/F -----
    hgb_factor <- (HGB / 13)^e_hgb_cl

    # One-sided renal effect: the power term applies only below the hinge
    # at 85 mL/min/1.73 m^2 and collapses to exactly 1 at or above it.
    # Raising the ratio to (exponent x indicator) is equivalent to the
    # piecewise definition and is continuous at the hinge, where the ratio
    # is 1 and the term equals 1 under either branch.
    crcl_low <- (CRCL < 85)
    crcl_factor <- (CRCL / 85)^(e_crcl_cl * crcl_low)

    # ----- 2. Individual parameters -----
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * hgb_factor * crcl_factor
    vc <- exp(lvc)

    # ----- 3. Micro-constants -----
    kel <- cl / vc

    # ----- 4. ODE system -----
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ----- 5. Observation and error -----
    # central is mg and vc is L, so central/vc is mg/L; x 1000 gives ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
