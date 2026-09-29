Cheng_2019_lobaplatin <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenous lobaplatin in",
    "elderly (>= 65 years) Chinese patients with small cell lung cancer",
    "(Cheng 2019; n = 113 patients across 7 centres, 678 plasma",
    "concentrations). Linear disposition with a small elimination clearance",
    "(0.478 L/h) and a long terminal phase. The central volume is an",
    "ADDITIVE linear function of body surface area (51.4 L, plus 47.7 L per",
    "m^2 above the 1.675 m^2 cohort median), and the inter-compartmental",
    "clearance is 31% higher in patients with creatinine clearance >= 80",
    "mL/min. Exponential inter-individual variability on all four",
    "disposition parameters and a proportional residual error. IMPORTANT:",
    "Cheng 2019 Table 3 labels the two clearances CL1 (0.478 L/h) and CL2",
    "(12.1 L/h) without saying which is elimination and which is",
    "distribution, and the Discussion reads CL2 as the renally relevant",
    "clearance. The paper's own published AUC simulations and its visual",
    "predictive check both identify CL1 as the ELIMINATION clearance and",
    "CL2 as the INTER-COMPARTMENTAL clearance; that assignment is what is",
    "encoded here and is gated in the validation vignette."
  )
  reference <- paste(
    "Cheng Y, Wu L, Liu X, Zhao Y, Liu C, Chen Q, Sun T, Zheng Q (2019).",
    "Population pharmacokinetics and individualized lobaplatin regimen for",
    "the treatment of Chinese small cell lung cancer in the elderly.",
    "Medicine (Baltimore) 98(3):e14136. doi:10.1097/MD.0000000000014136.",
    sep = " "
  )
  vignette <- "Cheng_2019_lobaplatin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Cheng 2019 Methods 2.3/2.4: the 4 mL plasma fraction was
  # assayed for lobaplatin by LC-MS/MS (diphenhydramine internal standard); the
  # separate 1 mL whole-blood aliquot for total platinum is NOT the analyte
  # modelled here. Figure 3 axis: 'Concentration, mg/L'.
  compartmentData <- list(
    central = list(analyte = "lobaplatin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "lobaplatin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The only size descriptor retained in the final model. Cheng 2019",
        "Results 3.5 'Fixed effect screening': 'Since height, weight, BMI,",
        "and BSA are correlated (collinearity), with CL and V theoretically",
        "affected by weight, BSA was preferred for covariate screening'.",
        "The RETAINED form is ADDITIVE LINEAR, not the allometric /",
        "exponential form that the screening paragraph mentions: Results",
        "3.5 and 3.6 both print 'V1 = 51.4 + (BSA-1.675) x 47.7' and both",
        "gloss it as 'V1 was increased by 4.77 L for each 0.1 BSA increment",
        "compared with 1.675', which is self-consistent only for the",
        "additive reading (0.1 x 47.7 = 4.77). The generic linear template",
        "in Methods 2.6.3, 'PTV = theta1 x [1 + theta2 x (COV - COVmedian)]',",
        "would instead give 51.4 x (1 + 47.7 x (BSA - 1.675)), i.e. 296 L at",
        "BSA 1.775, so the printed final equation governs. dOFV -8.825.",
        "Centering value 1.675 m^2 is the cohort median; the cohort range is",
        "1.24 to 2.09 m^2 (Results 3.8.1), which bounds the validated",
        "extrapolation range. Note the additive form makes Vc non-positive",
        "below BSA = 1.675 - 51.4/47.7 = 0.597 m^2."
      ),
      source_name = "BSA"
    ),
    CRCL = list(
      description = "Creatinine clearance, Cockcroft-Gault (raw, NOT BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters the model only through the binary threshold CRCL >= 80",
        "mL/min (the paper's covariate CCRG, the enrolment grouping),",
        "derived inside model(). Cheng 2019 inclusion criteria: 'Ccr >= 60",
        "ml/min, obtained by Cockcroft-Gault formula', so the column is raw",
        "mL/min with no BSA normalization, matching the raw-mL/min variant",
        "of the canonical CRCL column (Delattre_2010_amikacin.R,",
        "Chen_2023_nemonoxacin.R precedents). Group allocation: group A =",
        "Ccr >= 80 mL/min (n = 51), group B = 60 <= Ccr < 80 mL/min",
        "(n = 49); cohort Ccr 80.3 mL/min (Table 1). Results 3.5: 'patient",
        "grouping according to the Ccr also significantly affected CL2",
        "(dOFV -7.999), with CL2 values of 12.10 L/h and 15.85 L/h for Ccr",
        "< 80 ml/min and Ccr >= 80 ml/min, respectively, indicating a 31%",
        "increase'. No patient with Ccr < 60 mL/min was enrolled (Discussion",
        "'The limitations of this study'), so the model is not validated",
        "below that."
      ),
      source_name = "Ccr (grouped as CCRG)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 113,
    n_studies = 1,
    age_range = ">= 65 years (inclusion criterion)",
    age_median = "68 years (range 66-71, Table 1)",
    sex_female_pct = 25,
    race_ethnicity = c(Asian = 100),
    disease_state = "small cell lung cancer (30 limited stage, 70 extensive stage)",
    dose_range = paste(
      "30 mg/m^2 (creatinine clearance >= 80 mL/min) or 20 mg/m^2",
      "(60 <= creatinine clearance < 80 mL/min) as a 2 h intravenous",
      "infusion, up to 4 cycles; 56 patients at each dose level"
    ),
    regions = "China (7 centres)",
    renal_function = "creatinine clearance >= 60 mL/min (Cockcroft-Gault) required for enrolment",
    co_medication = "etoposide (60-100 mg/m^2) or irinotecan (120-200 mg/m^2) in all regimens",
    notes = paste(
      "113 patients contributed to the PK analysis: 100 sparsely sampled in",
      "the population study (4 random time points in cycle 1, 1 sample at",
      "4 h in later cycles) plus 13 richly sampled at Hunan Cancer Hospital",
      "(10 time points). 678 plasma concentrations, 17 below the limit of",
      "quantitation. Baseline demographics are Cheng 2019 Table 1; PK cohort",
      "counts are Results 3.4. The paper does not print the infusion",
      "duration or the cycle length; both are recovered from the figures and",
      "are recorded in the validation vignette's Assumptions section (2 h",
      "infusion from the Figure 2 CWRES-vs-TIME axis, which shifts every",
      "protocol sampling time by exactly +2 h; ~21-day cycles from the",
      "~1500 h span of the 4-cycle panel of Figure 3)."
    )
  )

  ini({
    # Cheng 2019 Table 3 'Final parameter estimation of the Pop PK model',
    # column 'Estimates (SE%)'. Bootstrap medians and 95% CIs from the same
    # table (1000 replicates, 781 successes) are quoted in the comments.
    #
    # CL1 / CL2 ASSIGNMENT. Table 3 prints 'CL1, L/h = 0.478' and
    # 'CL2, L/h = 12.1' and never states which is elimination and which is
    # distribution; the Discussion treats CL2 as the renally relevant
    # clearance. Three independent results in the paper itself show the
    # opposite, i.e. CL1 = elimination CL and CL2 = Q:
    #   (1) The three published AUC-ratio simulations (Results 3.8.2-3.8.4)
    #       are reproduced to within 1 percentage point by CL1 = CL over the
    #       0-24 h sampling window (equal 30 mg/m^2 dosing: group B +8.9% vs
    #       the published +8%; protocol 20 vs 30 mg/m^2: group A +37.8% vs
    #       the published 39%; adjusted 27 vs 30 mg/m^2: group B -2.0% vs the
    #       published -3%). With CL2 = CL the same three come out +29.6%,
    #       +15.7% and +16.6%.
    #   (2) The 0-26 h panel of Figure 3 has a median of ~0.13 mg/L at 26 h.
    #       CL1 = CL predicts 0.18 mg/L; CL2 = CL predicts 0.0025 mg/L.
    #   (3) The 4-cycle panel of Figure 3 shows 0.1-0.3 mg/L persisting to
    #       1500 h with visible cycle-to-cycle accumulation, which requires
    #       the long terminal phase that only CL1 = CL produces.
    # The vignette gates all three.
    lvc <- log(51.4)
    label("Central volume V1 at the median BSA of 1.675 m^2 (L)") # Table 3 'V1, L' = 51.4 (SE 5.1%); bootstrap 51.5 (46.3-57.2)
    lvp <- log(202.0)
    label("Peripheral volume V2 (L)") # Table 3 'V2, L' = 202.0 (SE 5.8%); bootstrap 202.8 (177.6-230.4)
    lcl <- log(0.478)
    label("Elimination clearance, the paper's CL1 (L/h)") # Table 3 'CL1, L/h' = 0.478 (SE 8.3%); bootstrap 0.478 (0.365-0.602)
    lq <- log(12.1)
    label("Inter-compartmental clearance at Ccr < 80 mL/min, the paper's CL2 (L/h)") # Table 3 'CL2, L/h' = 12.1 (SE 5.9%); bootstrap 12.2 (10.5-14.2)

    # Covariate effects.
    # BSA on Vc is an ADDITIVE slope in L per m^2, NOT a power exponent and
    # NOT the multiplicative linear template of Methods 2.6.3 -- see the
    # covariateData[['BSA']] notes for why the printed final equation
    # 'V1 = 51.4 + (BSA-1.675) x 47.7' governs.
    e_bsa_vc <- 47.7
    label("Additive slope of BSA on Vc, centred at 1.675 m^2 (L per m^2)") # Table 3 'Theta of the BSA on V1' = 47.7 (SE 11.5%); bootstrap 47.1 (8.8-78.5); Results 3.5/3.6 'V1 = 51.4 + (BSA-1.675) x 47.7'
    # CCRG on Q is reported as the ratio itself, not as the theta2 of the
    # Methods 2.6.3 piece-wise template 'PTV = theta1 x (1 + theta2)':
    # 12.1 x 1.31 = 15.85 L/h reproduces the printed group-A value exactly,
    # whereas 12.1 x (1 + 1.31) = 27.95 L/h does not.
    e_crcl_q <- 1.31
    label("Multiplicative factor on Q for creatinine clearance >= 80 mL/min (power-form base)") # Table 3 'Theta of CCRG on CL 2' = 1.31 (SE 8.5%); bootstrap 1.30 (1.05-1.59); Results 3.5 '12.10 L/h and 15.85 L/h'

    # Inter-individual variability. Methods 2.6.2 / Results 3.6 state the
    # exponential model 'Pi = PTV x exp(eta_i)', so the Table 3 'Inter-
    # individual variability ... %' rows are read as log-normal CV% and
    # converted exactly with omega^2 = log(1 + CV^2). The alternative
    # reading, that the printed number is the raw NONMEM variance x 100
    # (omega^2 = 0.514 etc.), is excluded by the 0-26 h panel of Figure 3:
    # simulating the study design gives a mean log-SD of 0.60 under that
    # reading against 0.38 measured off the published percentile lines,
    # while the CV% reading gives 0.44. See the vignette Errata.
    # Source-table text is quoted with SINGLE quotes: rxode2 promotes a
    # trailing comment on an eta line into label(), and a double quote there
    # breaks the re-parse.
    etalvc ~ 0.234455 # Table 3 'Inter-individual variability / V 1 , %' = 51.4 (SE 9.5%); bootstrap 50.5 (41.2-60.1); log(1 + 0.514^2)
    etalvp ~ 0.298684 # Table 3 'Inter-individual variability / V 2 , %' = 59.0 (SE 9.3%); bootstrap 58.4 (48.2-69.4); log(1 + 0.590^2)
    etalcl ~ 0.480067 # Table 3 'Inter-individual variability / CL 1 , %' = 78.5 (SE 23.9%); bootstrap 76.1 (36.2-117.3); log(1 + 0.785^2)
    etalq ~ 0.144306 # Table 3 'Inter-individual variability / CL 2 , %' = 39.4 (SE 14.3%); bootstrap 38.6 (27.5-48.3); log(1 + 0.394^2)

    # Residual error. Results 3.6 prints the retained form as
    # 'Cobs,ij = Cpred,ij x (1 + eps_ij,1)', i.e. proportional only, and
    # Table 3 has a single 'Residual variability / Proportional error, %'
    # row = 3.21. That row is read as the NONMEM sigma^2 x 100, giving
    # sqrt(0.0321) = 17.9%. Taken at face value it would be a 3.21%
    # proportional error, which is below the +/-15% precision that the CFDA
    # bioanalytical guideline cited in Methods 2.4 permits for the assay
    # itself, and is contradicted by the paper's own IPRED-vs-OBS panel
    # (Figure 2, upper left), whose scatter about the identity line is
    # roughly 13-18% at every prediction level. See the vignette Errata.
    propSd <- sqrt(0.0321)
    label("Proportional residual error (fraction)") # Table 3 'Proportional error, %' = 3.21 (SE 14.2%); bootstrap 3.17 (2.39-4.18) read as sigma^2 x 100
  })

  model({
    # The paper's CCRG covariate is the enrolment grouping: 1 for group A
    # (Ccr >= 80 mL/min), 0 for group B (60 <= Ccr < 80 mL/min).
    crclHigh <- (CRCL >= 80)

    # Individual disposition parameters. The BSA term is additive and sits
    # INSIDE the exponential IIV, matching the NONMEM pattern
    # 'TVV1 = THETA(1) + (BSA - 1.675) x THETA(5); V1 = TVV1 x EXP(ETA(1))'
    # that Results 3.5/3.6's printed equation and Methods 2.6.2's
    # 'Pi = PTV x exp(eta_i)' together imply.
    vc <- (exp(lvc) + e_bsa_vc * (BSA - 1.675)) * exp(etalvc)
    vp <- exp(lvp + etalvp)
    cl <- exp(lcl + etalcl)
    q <- exp(lq + etalq) * e_crcl_q^crclHigh

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Two-compartment disposition; lobaplatin is given as an intravenous
    # infusion, so there is no absorption compartment and no bioavailability
    # term. The infusion duration is supplied by the event table.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Plasma lobaplatin concentration in mg/L (dose in mg, volumes in L).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
