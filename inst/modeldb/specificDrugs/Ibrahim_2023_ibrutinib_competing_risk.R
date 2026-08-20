Ibrahim_2023_ibrutinib_competing_risk <- function() {
  description <- paste(
    "Competing-risk multistate model for study dropout and death in patients with chronic lymphocytic leukemia",
    "treated with ibrutinib (phase Ib/II PCYC-1102). Three ODE states track the marginal probability of being alive",
    "and on study (s_alive, the only transient state), of having dropped out (s_dropout), and of having died",
    "(s_death); both terminal states are absorbing. The transition rate to dropout falls as the model-predicted",
    "leukocyte count rises (hazard ratio 4.92 for every 10-unit decrease in leukocyte count -- patients whose",
    "predicted count fell below normal, potentially reflecting neutropenia, were the most likely to leave the",
    "study), while the transition rate to death rises with the model-predicted lymph-node burden (hazard ratio 1.35",
    "for every 10-unit increase in SPD) and with the deletion(17p) chromosomal abnormality (hazard ratio 4.16).",
    "The two time-varying predictors are supplied as covariates and are intended to be the model-predicted leukocyte",
    "count and SPD from the sister efficacy model, held constant between visits as in the source analysis. There are",
    "no drug-dosing events in this model and ibrutinib exposure does not enter it directly. Sister model files from",
    "the same paper: modellib('Ibrahim_2023_ibrutinib_leukocyte_spd'), modellib('Ibrahim_2023_ibrutinib_sbp'),",
    "modellib('Ibrahim_2023_ibrutinib_dbp')."
  )
  reference <- paste(
    "Ibrahim EIK, Karlsson MO, Friberg LE.",
    "Assessment of ibrutinib scheduling on leukocyte, lymph node size and blood pressure dynamics in chronic",
    "lymphocytic leukemia through pharmacokinetic-pharmacodynamic models.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(9):1305-1318.",
    "doi:10.1002/psp4.13010.",
    "Open Access under CC BY-NC.",
    "Structural equations from Appendix S2 equations S1-S3; parameter values and covariate forms from Table 2",
    "(competing risk model rows) and Table 2 footnotes e and f.",
    sep = " "
  )
  vignette <- "Ibrahim_2023_ibrutinib"
  paper_specific_compartments <- c("s_alive", "s_dropout", "s_death")
  units <- list(
    time          = "month",
    dosing        = "n/a (no drug-dosing events; the model propagates state-occupancy probabilities)",
    concentration = "probability (all three states are occupancy probabilities, not drug concentrations)"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    s_alive   = list(analyte = "probability of being alive and on study", units = NA_character_, specimen = "not applicable", verified = FALSE),
    s_dropout = list(analyte = "probability of having dropped out", units = NA_character_, specimen = "not applicable", verified = FALSE),
    s_death   = list(analyte = "probability of having died", units = NA_character_, specimen = "not applicable", verified = FALSE)
  )

  covariateData <- list(
    WBC = list(
      description        = "Model-predicted total leukocyte count from the preceding visit.",
      units              = "10^9 cells/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying, held constant between two visits (Ibrahim 2023 Patients and Methods 'Covariate analysis':",
        "'the model-predicted metrics for leukocyte count, SPD, sBP, and dBP dynamics from preceding visits were",
        "explored in the competing risk model. These metrics were assessed as time-varying covariates that remained",
        "constant between two visits'). The intended source is the leukocyte output of the sister model",
        "modellib('Ibrahim_2023_ibrutinib_leukocyte_spd'), which is already in 10^9 cells/L -- no conversion needed.",
        "Enters the dropout hazard in power form lambda12 = 0.00908 * (WBC / 12)^-0.89 (Ibrahim 2023 Table 2",
        "footnote e), i.e. reference value 12 x10^9 cells/L. Check against the reported effect size: a 10-unit",
        "decrease from the reference, 12 -> 2, gives (2/12)^-0.89 = 4.93, matching the reported hazard ratio of 4.92."
      ),
      source_name        = "Leukocyte"
    ),
    TUMSZ = list(
      description        = "Model-predicted lymph-node burden from the preceding visit, as the sum of the products of perpendicular diameters (SPD).",
      units              = "mm^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "UNIT WARNING: Ibrahim 2023 reports SPD in cm^2 throughout, but the canonical unit for an SPPD-type TUMSZ",
        "column in this package is mm^2, so the reference value below is the source's 14 cm^2 expressed as",
        "1400 mm^2. Because the effect is a power form the conversion is numerically invariant -- (TUMSZ/1400) with",
        "TUMSZ in mm^2 equals (SPD/14) with SPD in cm^2 -- but the sister model",
        "modellib('Ibrahim_2023_ibrutinib_leukocyte_spd') emits its tumorSpd output in cm^2, so multiply that output",
        "by 100 before supplying it here. Time-varying, held constant between two visits (Ibrahim 2023 Patients and",
        "Methods 'Covariate analysis'). Enters the death hazard in power form",
        "lambda13 = 0.00275 * (TUMSZ / 1400)^0.563 * exp(1.42 * TUM_17P_DEL) (Ibrahim 2023 Table 2 footnote f).",
        "Check against the reported effect size: a 10-unit increase from the reference, 14 -> 24 cm^2, gives",
        "(24/14)^0.563 = 1.35, matching the reported hazard ratio of 1.35."
      ),
      source_name        = "SPD"
    ),
    TUM_17P_DEL = list(
      description        = "Deletion(17p) chromosomal abnormality in the CLL clone: 1 = del(17p) present, 0 = absent.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no del(17p))",
      notes              = paste(
        "Time-fixed per subject; a baseline cytogenetic call. Enters the death hazard multiplicatively as",
        "exp(1.42 * TUM_17P_DEL), a hazard ratio of exp(1.42) = 4.14 (Ibrahim 2023 reports 4.16 from the unrounded",
        "coefficient; Table 2 'Coefficient of deletion (17p) on lambda13' = 1.42, 95% CI 0.153-2.7). Deletion(17p)",
        "removes the TP53 locus and defines the highest-risk CLL subgroup (10-year overall survival 29%;",
        "Ibrahim 2023 Discussion). Note that despite the Table 2 row being labelled 'on lambda12' in the printed",
        "table, footnote f places the del(17p) term in the lambda13 (death) equation, and the narrative states",
        "'Patients carrying deletion (17p) chromosomal abnormality had a statistically significant higher",
        "probability of DEATH'; the footnote equation is used here (see the vignette Assumptions and deviations)."
      ),
      source_name        = "del(17p)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 120L,
    n_studies      = 1L,
    age_range      = "mean 62.4 (SD 9.9) years",
    weight_range   = "mean 82.3 (SD 17) kg",
    sex_female_pct = 24.2,
    race_ethnicity = NULL,
    disease_state  = "Chronic lymphocytic leukemia (CLL); 20.8% treatment-naive, 79.2% relapsed/refractory",
    dose_range     = "ibrutinib 420 mg once daily (n = 94) or 840 mg once daily (n = 38) in PCYC-1102",
    regions        = "United States (PCYC-1102, phase Ib/II)",
    notes          = paste(
      "Baseline demographics from Ibrahim 2023 Supplementary Table S1. Patients were followed for a maximum of",
      "2.4 years (median 1.7 years). During the study 11 patients died and 22 dropped out before the end of the",
      "study period, for disease progression (n = 5), adverse events (n = 5) and other reasons (n = 12)",
      "(Ibrahim 2023 Appendix S1 section 1). Unlike the three PK-PD sister models, which were fitted with FOCEI in",
      "nlmixr 2.0.6, the competing risk analysis was performed with the msm R package (version 1.6.9)",
      "(Ibrahim 2023 Patients and Methods 'Model development and evaluation'), and Table 2 reports no",
      "between-subject variability for it -- hence this model has no eta parameters."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Ibrahim 2023 Table 2, 'Competing risk model' rows, with the covariate
    # forms given verbatim by Table 2 footnotes e and f:
    #   lambda12 = 0.00908 * e^(-0.89 * LN(Leukocyte / 12))
    #   lambda13 = 0.00275 * e^( 0.563 * LN(SPD / 14) + 1.42 * del(17p) )
    # Both baseline transition rate constants are per MONTH, so this model
    # runs on a month time axis (units$time above) while the three sister
    # PK-PD models run on a day axis.
    # ------------------------------------------------------------------
    llambda12 <- log(0.00908)
    label("Log baseline transition rate constant from alive to dropout (1/month)")  # Table 2 lambda12 = 0.00908 (95% CI 0.0053-0.0155)
    llambda13 <- log(0.00275)
    label("Log baseline transition rate constant from alive to death (1/month)")  # Table 2 lambda13 = 0.00275 (95% CI 0.00102-0.00743)

    e_wbc_lambda12 <- -0.89
    label("Power coefficient of the model-predicted leukocyte count (referenced to 12 x10^9 cells/L) on the dropout transition rate")  # Table 2 'Coefficient of past model-predicted leukocyte count on lambda12' = -0.89 (95% CI -1.47 to -0.304); footnote e
    e_tumsz_lambda13 <- 0.563
    label("Power coefficient of the model-predicted SPD (referenced to 1400 mm^2 = 14 cm^2) on the death transition rate")  # Table 2 'Coefficient of past model-predicted SPD on lambda13' = 0.563 (95% CI 0.077-1.05); footnote f
    e_del17p_lambda13 <- 1.42
    label("Log-hazard coefficient for deletion(17p) on the death transition rate; HR = exp(1.42) = 4.14")  # Table 2 'Coefficient of deletion (17p)' = 1.42 (95% CI 0.153-2.7); footnote f places it on lambda13

    # ------------------------------------------------------------------
    # No IIV. Table 2 reports no between-subject variability for any
    # competing-risk parameter (the msm multistate fit is a population
    # hazard model), so no eta* parameters are added.
    #
    # Residual error placeholder. The source analysis maximises the
    # multistate event likelihood in msm and has no observation-error
    # model. This nlmixr2 translation is intended for forward simulation:
    # the state-occupancy probabilities are exposed as outputs and a tiny
    # additive residual is attached to the overall-survival output so the
    # nlmixr2 likelihood machinery accepts the model. NOT from the source
    # -- see the vignette Assumptions and deviations section.
    # ------------------------------------------------------------------
    addSd <- 0.001
    label("Placeholder additive residual error on the overall-survival probability output (unitless); not from the source")
  })

  model({
    # ------------------------------------------------------------------
    # 1. Transition rate constants (Ibrahim 2023 Table 2 footnotes e, f).
    #    Both covariate effects are power functions of the model-predicted
    #    biomarker referenced to the value printed in the footnote, which
    #    is exactly exp(coef * log(x / ref)).
    # ------------------------------------------------------------------
    lambda12 <- exp(llambda12) * exp(e_wbc_lambda12 * log(WBC / 12))
    lambda13 <- exp(llambda13) *
      exp(e_tumsz_lambda13 * log(TUMSZ / 1400) + e_del17p_lambda13 * TUM_17P_DEL)

    # ------------------------------------------------------------------
    # 2. Multistate ODE system (Ibrahim 2023 Appendix S2, equations
    #    S1-S3). s_alive is the single transient state; s_dropout and
    #    s_death are absorbing. All patients start alive and on study.
    # ------------------------------------------------------------------
    d/dt(s_alive)   <- -lambda12 * s_alive - lambda13 * s_alive  # eq. S1
    s_alive(0)      <- 1

    d/dt(s_dropout) <- lambda12 * s_alive                        # eq. S2
    s_dropout(0)    <- 0

    d/dt(s_death)   <- lambda13 * s_alive                        # eq. S3
    s_death(0)      <- 0

    # ------------------------------------------------------------------
    # 3. Outputs. Overall survival is the complement of the cumulative
    #    death probability; dropout is a competing event that removes a
    #    patient from observation without them having died.
    # ------------------------------------------------------------------
    probAlive   <- s_alive     # alive AND still on study
    probDropout <- s_dropout   # cumulative marginal dropout probability
    probDeath   <- s_death     # cumulative marginal death probability

    # `sur` is the registered canonical survival-probability output; here it
    # is overall survival, i.e. the complement of the cumulative death
    # probability (a patient who has dropped out is still alive).
    sur <- 1 - s_death

    sur ~ add(addSd)
  })
}
