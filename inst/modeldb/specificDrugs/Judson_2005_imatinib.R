Judson_2005_imatinib <- function() {
  description <- paste0(
    "One-compartment population PK model with zero-order absorption of ",
    "duration D1 into the central compartment and first-order elimination ",
    "for oral imatinib in adults with gastrointestinal stromal tumor or ",
    "soft-tissue sarcoma (Judson 2005, EORTC Soft Tissue and Bone Sarcoma ",
    "Group). Apparent oral clearance CL/F carries a power-form total body ",
    "weight effect referenced to 67.5 kg whose exponent (0.002) is ",
    "numerically negligible. Inter-individual variability is estimated on ",
    "CL/F and Vc/F as a correlated 2x2 block; residual error is combined ",
    "proportional plus additive. TRANSCRIBED FROM A SECONDARY SOURCE: the ",
    "parameter values come from Table 1 of Yang 2025, an external ",
    "evaluation of 15 published imatinib population PK models, not from ",
    "the primary publication. Re-extract from Judson 2005 when that paper ",
    "is obtained."
  )
  reference <- paste0(
    "Judson I, Ma P, Peng B, Verweij J, Racine A, di Paola ED, van ",
    "Glabbeke M, Dimitrijevic S, Scurr M, Dumez H, van Oosterom A. ",
    "Imatinib pharmacokinetics in patients with gastrointestinal stromal ",
    "tumour: a retrospective population pharmacokinetic study over time. ",
    "EORTC Soft Tissue and Bone Sarcoma Group. Cancer Chemother Pharmacol. ",
    "2005;55(4):379-386. doi:10.1007/s00280-004-0876-0. ",
    "PARAMETER SOURCE (secondary): Yang T, Rasmussen ASB, Weimann A, ",
    "Thastrup M, Rank CU, Als-Nielsen B, Malmros J, Wik HS, Lohi O, ",
    "Overgaard U, Johannsdottir IMR, Vaitkeviciene G, Dalhoff K, ",
    "Schmiegelow K, Lund TM. Published population pharmacokinetic models ",
    "of imatinib perform poorly on TDM data from pediatric patients. ",
    "Target Oncol. 2025;20(5):871-886. Table 1, row 'Judson et al. ",
    "(2005)'. doi:10.1007/s11523-025-01172-2."
  )
  vignette <- "Yang_2025_imatinib_external_evaluation"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    central = list(analyte = "imatinib", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Enters CL/F as the power function (WT / 67.5)^0.002 (Yang 2025 ",
        "Table 1, 'CL/F: 10.6 x (TBW/67.5)^0.002'). The reference value ",
        "67.5 kg is the centring constant printed inside the covariate ",
        "term. The exponent is essentially zero, so the weight effect on ",
        "CL/F is negligible over any realistic weight range: a 2-fold ",
        "change in weight moves CL/F by 2^0.002 - 1 = 0.14%. It is ",
        "retained here so the published equation is reproduced exactly ",
        "rather than silently simplified. Yang 2025 nonetheless applied ",
        "STANDARD ALLOMETRIC SCALING to this model when evaluating it in a ",
        "cohort containing children, because the published weight effect ",
        "cannot describe pediatric clearance; that allometric variant is ",
        "Yang 2025's own modification and is NOT encoded here."
      ),
      source_name        = "TBW"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 43L,
    n_studies      = 1L,
    n_observations = "517 imatinib plasma concentrations (Yang 2025 Table 1)",
    age_range      = "not reported in Yang 2025 Table 1",
    disease_state  = "Adults with gastrointestinal stromal tumor (GIST) or soft-tissue sarcoma (STS)",
    dose_range     = "Oral imatinib 400-1000 mg total daily dose",
    regions        = "Europe and USA",
    bioanalytical  = "LC-MS/MS, limit of quantification 4 ng/mL (Yang 2025 Table 1)",
    notes          = paste0(
      "Yang 2025 Table 1 footnote a: 'The population PK model was ",
      "described differently on day 1, day 29, and the extension phase. ",
      "The model on the extension phase was used in this current ",
      "validation study.' The parameters encoded here are therefore the ",
      "EXTENSION-PHASE model, not the day-1 or day-29 model. ",
      "Demographic detail beyond the row above (weight range, sex split, ",
      "race) is not reported by the secondary source and must be read from ",
      "the primary publication."
    )
  )

  ini({
    # ----- Structural parameters (Yang 2025 Table 1, Judson row) -----
    # Typical apparent oral clearance at the reference weight of 67.5 kg.
    lcl <- log(10.6); label("Apparent oral clearance CL/F at WT = 67.5 kg (L/h)")  # Yang 2025 Table 1: CL/F = 10.6 x (TBW/67.5)^0.002

    lvc <- log(182); label("Apparent central volume of distribution Vc/F (L)")  # Yang 2025 Table 1: Vc/F = 182

    # Zero-order absorption duration. The imatinib dose is delivered into
    # the central compartment at a constant rate over d1 hours; there is no
    # separate depot state in a zero-order absorption model.
    ld1 <- log(1.51); label("Zero-order absorption duration D1 (h)")  # Yang 2025 Table 1: D1 = 1.51

    # ----- Covariate effect on CL/F -----
    e_wt_cl <- 0.002; label("Power exponent of (WT / 67.5) on CL/F (unitless)")  # Yang 2025 Table 1: (TBW/67.5)^0.002

    # ----- Inter-individual variability -----
    # Yang 2025 Table 1 reports this model's BSV on the VARIANCE scale --
    # 'Var(eta_CL): 0.305', 'Var(eta_Vc): 0.34', 'Cov(eta_CL, eta_Vc):
    # 0.237' -- unlike most rows of that table, which report CV%. The
    # variance reading is not a guess: it is forced by the covariance.
    # Treating 0.305 and 0.34 as variances gives a correlation of
    # 0.237 / sqrt(0.305 * 0.34) = 0.736, which is admissible; treating them
    # as standard deviations would give 0.237 / sqrt(0.093 * 0.116) = 2.29,
    # which is impossible for a correlation. nlmixr2 takes variances and
    # covariances in this block, so the tabulated numbers are used directly.
    etalcl + etalvc ~ c(0.305,
                        0.237, 0.34)  # Yang 2025 Table 1: Var(eta_CL) 0.305, Cov 0.237, Var(eta_Vc) 0.34

    # ----- Residual unexplained variability -----
    # Combined proportional plus additive. Cc is reported in ng/mL, and the
    # additive term is tabulated in ng/mL, so no unit conversion is needed.
    propSd <- 0.34; label("Proportional residual error (fraction)")  # Yang 2025 Table 1: Prop 34%
    addSd <- 273; label("Additive residual error (ng/mL)")  # Yang 2025 Table 1: Add 273 ng/mL
  })

  model({
    # ----- 1. Individual parameters -----
    cl <- exp(lcl + etalcl) * (WT / 67.5)^e_wt_cl
    vc <- exp(lvc + etalvc)
    d1 <- exp(ld1)

    # ----- 2. Micro-constants -----
    kel <- cl / vc

    # ----- 3. ODE system -----
    # One compartment, zero-order input. Dose records must carry rate = -2
    # so rxode2 uses the modelled duration dur(central) = d1; without it the
    # dose collapses to an instantaneous bolus.
    d/dt(central) <- -kel * central
    dur(central) <- d1

    # ----- 4. Observation and error -----
    # `central` is an amount in mg and `vc` is in L, so central/vc is mg/L;
    # the factor 1000 converts to ng/mL, the unit in which imatinib TDM
    # targets (1000-3000 ng/mL trough) and this model's additive residual
    # error are expressed.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
