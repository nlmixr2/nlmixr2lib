Takechi_2025_nemolizumab_mbma_iga <- function() {
  description <- paste0(
    "MBMA. Longitudinal logistic model-based meta-analysis of the ",
    "Investigator's Global Assessment (IGA) success rate - the proportion of ",
    "patients achieving at least a 2-point decrease in IGA - in prurigo ",
    "nodularis, comparing nemolizumab with dupilumab against a shared placebo ",
    "response. Fitted to 13 study arms from six randomized controlled trials ",
    "(1170 participants) identified by a PRISMA literature search: the ",
    "Japanese phase II/III study of Yokozeki et al., the phase II study of ",
    "Staender et al., OLYMPIA 1 and OLYMPIA 2 for nemolizumab, and PRIME and ",
    "PRIME2 for dupilumab. On the logit scale the treatment effect is an ",
    "intercept plus a mono-exponential placebo rise shared by every arm plus ",
    "a mono-exponential drug rise whose maximum and onset rate are estimated ",
    "separately for each drug; the responder probability is the logistic ",
    "transform of that sum. Both maxima are constants because most trials ",
    "used the approved dose and no dose-response was detectable, so the ",
    "DOSE_NEMOLIZUMAB_MG and DOSE_DUPILUMAB_MG columns act as arm selectors. ",
    "No covariate was retained for this endpoint. Nemolizumab has both the ",
    "larger maximum (1.71 vs 1.49) and the faster onset (0.416 vs 0.134 per ",
    "week), which is the source's basis for its rapid-onset claim. ",
    "Variability is BETWEEN-STUDY only (one eta added to the logit-scale ",
    "effect); the simulation scope is the study-arm success rate, NOT an ",
    "individual patient trajectory. Sister models from the same paper: ",
    "modellib('Takechi_2025_nemolizumab') (popPK) and ",
    "modellib('Takechi_2025_nemolizumab_ppnrs') (population PD for weekly ",
    "average Peak Pruritus NRS)."
  )

  reference <- paste(
    "Takechi T, Shimizu J, Kabashima K, Ieiri I.",
    "Quantitative evaluation of nemolizumab pharmacokinetics and efficacy in",
    "prurigo nodularis: a population pharmacokinetics and model-based",
    "meta-analysis approach.",
    "Dermatol Ther (Heidelb). 2025;15(12):3615-3632.",
    "doi:10.1007/s13555-025-01554-4.",
    "Structural model: Methods 'Model-Based Meta-analysis', the displayed",
    "logit / EFF / f0 / fdrug equation group and the displayed residual",
    "equation.",
    "Parameter values: Table 4, 'IGA' column.",
    "Included trials and arm sizes: Table 3.",
    sep = " "
  )

  vignette <- "Takechi_2025_nemolizumab"

  # The single random effect is a BETWEEN-STUDY effect added directly to the
  # logit-scale treatment effect EFF, not an IIV on a named fixed-effect
  # parameter, so it cannot follow the popPK eta<transformed-param-name>
  # convention. Declared paper-specific per the Goteti_2024_SLE_mbma
  # precedent.
  paper_specific_etas <- c("eta_study_eff")

  units <- list(
    time          = "week",
    dosing        = paste0(
      "mg/administration (nemolizumab 30 mg or 60 mg Q4W, or 0.5 mg/kg Q4W in ",
      "SPR.115828; dupilumab 300 mg Q2W. Supplied through the ",
      "DOSE_NEMOLIZUMAB_MG and DOSE_DUPILUMAB_MG covariate columns and NOT as ",
      "rxode2 dose events; no dose-response was identifiable, so any non-zero ",
      "value selects that drug's constant maximum effect)"
    ),
    concentration = paste0(
      "fraction/arm (probability that a patient in the study arm achieves a ",
      "decrease of at least 2 points in the Investigator's Global Assessment; ",
      "this output is NOT a drug concentration - the slash satisfies ",
      "checkModelConventions parsing)"
    )
  )

  covariateData <- list(
    DOSE_NEMOLIZUMAB_MG = list(
      description        = "Assigned subcutaneous nemolizumab dose for the study arm; 0 in placebo and dupilumab arms.",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Arm selector only: consumed as (DOSE_NEMOLIZUMAB_MG > 0). Takechi ",
        "2025 Results 'MBMA': 'As most trials were conducted using the ",
        "approved dose, and no clear dose-response relationship was observed, ",
        "the maximum effect was modeled as a constant.' Two of the four ",
        "nemolizumab arms in Table 3 are not clean flat milligram doses - ",
        "OLYMPIA 1 and OLYMPIA 2 report a pooled '30 mg, 60 mg, Q4W' arm, and ",
        "SPR.115828 dosed 0.5 mg/kg Q4W - so only presence is well defined ",
        "here. Setting both this column and DOSE_DUPILUMAB_MG non-zero is ",
        "outside the source's calibration and would make the model ADD the ",
        "two drug effects; no trial in the meta-analysis combined them."
      ),
      source_name        = "Drug / Dose (Takechi 2025 Table 3, Nemolizumab rows)"
    ),
    DOSE_DUPILUMAB_MG = list(
      description        = "Assigned subcutaneous dupilumab dose for the study arm; 0 in placebo and nemolizumab arms.",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Arm selector only: consumed as (DOSE_DUPILUMAB_MG > 0). Every ",
        "included dupilumab arm (PRIME and PRIME2) used 300 mg Q2W, so this ",
        "column takes only the values 0 and 300 across the meta-analysis and ",
        "no dose-response is identifiable. NOTE the different dosing interval ",
        "from nemolizumab (Q2W vs Q4W): the two dose columns are not a common ",
        "per-administration metric and must not be compared numerically. The ",
        "dupilumab trials also report the Worst Itch NRS and IGA for PN ",
        "(IGA PN-S) rather than PP-NRS and IGA; the source treats these as ",
        "conceptually equivalent while cautioning against over-interpretation ",
        "(Discussion, final limitation)."
      ),
      source_name        = "Drug / Dose (Takechi 2025 Table 3, Dupilumab rows)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1170L,
    n_studies      = 6L,
    n_arms         = 13L,
    age_range      = "study-arm mean ages 46.7-59.7 years",
    weight_range   = "study-arm mean weights 60.4-87.1 kg",
    sex_female_pct = 51.9,
    disease_state  = paste0(
      "Moderate-to-severe prurigo nodularis. Across arms the baseline PP-NRS ",
      "mean ranged from 8.3 to 8.6 and the proportion of patients with an IGA ",
      "score of 4 (severe) ranged from 28.0% to 54.5% (Table 3). Takechi 2025 ",
      "notes that these baseline characteristics were similar across studies, ",
      "'suggesting that population heterogeneity had a limited impact on the ",
      "results'."
    ),
    dose_range     = paste0(
      "Nemolizumab 30 mg and 60 mg Q4W (Yokozeki et al., OLYMPIA 1, OLYMPIA ",
      "2) and 0.5 mg/kg Q4W (SPR.115828); dupilumab 300 mg Q2W (PRIME, ",
      "PRIME2); matched placebo in every trial. Trial durations 12-24 weeks."
    ),
    regions        = paste0(
      "Japan (Yokozeki et al. / jRCT2011200017) and multinational North ",
      "America / Europe (SPR.115828, OLYMPIA 1, OLYMPIA 2, PRIME, PRIME2)."
    ),
    notes          = paste0(
      "MBMA AT THE STUDY-ARM LEVEL: each modelled data point is one trial ",
      "arm's reported IGA success rate at one time point, digitised from ",
      "published figures with WebPlotDigitizer 4.5 where not tabulated. The ",
      "model is intended for simulating study-arm success rates and is NOT ",
      "suitable for individual-subject simulation. Sex percentages: the ",
      "sex_female_pct above is 100 minus the participant-weighted mean of the ",
      "per-arm 'Male (%)' column of Table 3 and is approximate. Arm sizes ",
      "(Table 3) are 76/77/76, 36/34, 96/190, 91/183, 76/75 and 82/78, ",
      "totalling 1170; they matter because the residual is scaled by the ",
      "inverse square root of the arm size - see the addSd note in ini(). ",
      "Selection: 104 records identified, 20 duplicates removed, 84 screened, ",
      "20 full texts assessed, 6 trials included (Fig. S7 PRISMA diagram). ",
      "The screened-but-not-retained covariates for this endpoint were age, ",
      "sex, race, body weight, baseline PP-NRS and IGA scores; none improved ",
      "the IGA model. In particular the proportion of patients classified as ",
      "moderate by IGA, which IS retained on the nemolizumab maximum effect ",
      "in the companion PP-NRS meta-analysis, was NOT retained here (Table 4 ",
      "prints a dash in the IGA column of that row)."
    )
  )

  ini({
    # ========================================================================
    # All values are Takechi 2025 Table 4, 'IGA' column.
    #
    # STRUCTURE (Methods 'Model-Based Meta-analysis'):
    #   Pr(IGA) = 1 / (1 + exp(-EFF))
    #   EFF     = f0 + fdrug + eta_study
    #   f0      = A + Pmax  * (1 - exp(-Kp * time))
    #   fdrug   = Edrug * (1 - exp(-Kd * time))
    # with Edrug and Kd estimated separately per drug. Time is in WEEKS: the
    # rate constants are printed as "(/week)" in Table 4 and the source's
    # Figs. 3 and 4 use a 0-24 week axis.
    #
    # VARIANCE SCALE. Table 4 gives an estimate and an %RSE per row and no CV%
    # label anywhere, so the 'ISV on EFF' and 'Additive error' rows are raw
    # NONMEM variances, matching the convention proved for Table 2 in the
    # companion popPK model file. Both are converted to SDs below.
    #
    # These four structural values reproduce the source's own Fig. 4 typical
    # curves to within read-off precision at three independent time points
    # (placebo 11.0% at week 24, dupilumab 26.3% at week 16 and 34.2% at week
    # 24, nemolizumab 34.5% at week 16 and 40.7% at week 24). The vignette
    # runs that check as a regression test.
    # ========================================================================

    e0 <- -4.92
    label("Intercept of the placebo effect on the logit scale at time zero (logit units); expit(-4.92) = 0.72% IGA success at baseline")
    # Table 4 IGA column, row A: -4.92 (%RSE 11.8).

    pmax_placebo <- 2.99
    label("Maximum placebo effect on the logit scale (logit units); the plateau placebo IGA success rate is expit(-4.92 + 2.99) = 12.5%")
    # Table 4 IGA column, row Pmax: 2.99 (%RSE 16.5).

    lkp <- log(0.123)
    label("Log rate constant for the onset of the placebo effect (log 1/week); back-transform Kp = 0.123 1/week, onset half-life 5.6 weeks")
    # Table 4 IGA column, row Kp (/week): 0.123 (%RSE 8.4).

    edrug_nemolizumab <- 1.71
    label("Maximum nemolizumab effect on the logit scale (logit units)")
    # Table 4 IGA column, row Edrug,nemolizumab: 1.71 (%RSE 9.5). Table 4's
    # 'PPM effect on Edrug,nemolizumab' row is a dash for IGA, i.e. the
    # covariate retained for the PP-NRS endpoint was NOT retained here.

    lkdrug_nemolizumab <- log(0.416)
    label("Log rate constant for the onset of the nemolizumab effect (log 1/week); back-transform Kd = 0.416 1/week, onset half-life 1.7 weeks")
    # Table 4 IGA column, row Kd,nemolizumab (/week): 0.416 (%RSE 47.1).

    edrug_dupilumab <- 1.49
    label("Maximum dupilumab effect on the logit scale (logit units)")
    # Table 4 IGA column, row Edrug,dupilumab: 1.49 (%RSE 7.7).

    lkdrug_dupilumab <- log(0.134)
    label("Log rate constant for the onset of the dupilumab effect (log 1/week); back-transform Kd = 0.134 1/week, onset half-life 5.2 weeks")
    # Table 4 IGA column, row Kd,dupilumab (/week): 0.134 (%RSE 52.7). Roughly
    # a third of the nemolizumab rate constant, which is what produces the
    # source's separation between the two drugs in early weeks even though
    # their plateaus are closer together.

    # ---- Between-STUDY variability ------------------------------------------
    # Methods: "The between-study variability is described by eta, which was
    # assumed to follow a normal distribution (mean = 0, variance = omega^2)."
    # ONE eta, added directly to the logit-scale EFF, shared by every arm of a
    # study. Named eta_study_* per the MBMA convention in this package so it
    # cannot be mistaken for a between-SUBJECT eta - there is no
    # between-subject variability in this model at all.
    eta_study_eff ~ 0.226
    # Table 4 IGA column, row 'ISV on EFF': 0.226 (%RSE 62.8). Variance, so
    # the between-study SD on the logit scale is sqrt(0.226) = 0.475.

    # ---- Residual error -----------------------------------------------------
    addSd <- 0.786765
    label("Residual standard deviation multiplier on the arm-level IGA success rate; the per-observation SD is addSd * sqrt(Pr * (1 - Pr) / N_arm)")
    # Table 4 IGA column, row 'Additive error': 0.619 (%RSE 25.4), read as a
    # VARIANCE, so the multiplier is sqrt(0.619) = 0.786765.
    #
    # ARM-SIZE WEIGHTING. The source's displayed residual equation is
    #   success rate (%) = Pr + eps * (Pr * (1 - Pr) / N)^(1/2)
    # i.e. the residual SD is this multiplier times the BINOMIAL standard
    # error of an arm of N subjects. model() encodes the prediction-dependent
    # factor sqrt(Pr * (1 - Pr)) exactly, but N is a property of the study arm
    # being simulated rather than of the model, so the 1 / sqrt(N) factor is
    # left to the user - the same convention as Boucher_2018_naproxen_mbma,
    # whose residual is likewise reweighted per arm downstream. The model as
    # shipped therefore emits the N = 1 residual; DIVIDE THE SIMULATED
    # RESIDUAL BY sqrt(N_arm), or equivalently set addSd to 0.786765 /
    # sqrt(N_arm), to reproduce an arm of N subjects. The derived variable
    # armSeUnit in model() exposes sqrt(Pr * (1 - Pr)) so the rescaling can be
    # done directly from the solve output. Arm sizes for the six included
    # trials are in Table 3 and are repeated in population$notes. The vignette
    # demonstrates the rescaling.
  })

  model({
    # ---- 1. Arm selectors ----------------------------------------------------
    # Constant maximum effect per drug; there is no dose-response.
    onNemolizumab <- (DOSE_NEMOLIZUMAB_MG > 0)
    onDupilumab   <- (DOSE_DUPILUMAB_MG   > 0)

    # ---- 2. Rate constants ---------------------------------------------------
    kp             <- exp(lkp)
    kdNemolizumab  <- exp(lkdrug_nemolizumab)
    kdDupilumab    <- exp(lkdrug_dupilumab)

    # ---- 3. Placebo and drug effects on the logit scale ----------------------
    # f0 = A + Pmax * (1 - exp(-Kp * time)); present in every arm including
    # the active arms, which is why an active arm's plateau is the sum of the
    # placebo plateau and its own drug maximum.
    f0 <- e0 + pmax_placebo * (1 - exp(-kp * time))

    # fdrug = Edrug * (1 - exp(-Kd * time)); at most one term is non-zero
    # because no included trial combined the two drugs.
    fdrug <-
      edrug_nemolizumab * (1 - exp(-kdNemolizumab * time)) * onNemolizumab +
      edrug_dupilumab   * (1 - exp(-kdDupilumab   * time)) * onDupilumab

    # ---- 4. Logit transform (Methods 'Model-Based Meta-analysis') ------------
    # The between-study eta is added on the logit scale, before the transform.
    eff   <- f0 + fdrug + eta_study_eff
    pResp <- expit(eff)

    # ---- 5. Observation and error -------------------------------------------
    # The source's residual SD is addSd * sqrt(Pr * (1 - Pr) / N_arm). rxode2's
    # prop() multiplies its argument by the prediction, so passing the computed
    # weight addSd * sqrt((1 - Pr) / Pr) yields
    #   SD = addSd * sqrt((1 - Pr) / Pr) * Pr = addSd * sqrt(Pr * (1 - Pr)),
    # which is the source's equation with N = 1. See the addSd note in ini()
    # for the per-arm rescaling by 1 / sqrt(N_arm).
    armSeUnit <- sqrt(pResp * (1 - pResp))
    propBinom <- addSd * sqrt((1 - pResp) / pResp)

    # NAMING. The single output is named `Cc` because that is the package's
    # canonical single-output observation variable (R/conventions.R
    # `observationVar`), NOT because it is a concentration - it is an IGA
    # success-rate probability. Same convention as `Goteti_2024_SLE_mbma.R`.
    Cc <- pResp
    Cc ~ prop(propBinom)
  })
}
