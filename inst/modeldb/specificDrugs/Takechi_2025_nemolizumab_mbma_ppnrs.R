Takechi_2025_nemolizumab_mbma_ppnrs <- function() {
  description <- paste0(
    "MBMA. Longitudinal logistic model-based meta-analysis of the Peak ",
    "Pruritus Numerical Rating Scale (PP-NRS) success rate - the proportion ",
    "of patients achieving at least a 4-point improvement - in prurigo ",
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
    "UNLIKE its IGA sibling this endpoint RETAINS one arm-level covariate: ",
    "the proportion of the arm classified as moderate by IGA ",
    "(IGA_MOD_PCT) scales the nemolizumab maximum effect, so an arm enriched ",
    "for severe disease shows a larger itch response. Nemolizumab has the ",
    "faster onset by a wide margin (1.93 vs 0.119 per week), which is the ",
    "source's basis for its rapid-itch-relief claim. Variability is ",
    "BETWEEN-STUDY only (one eta added to the logit-scale effect); the ",
    "simulation scope is the study-arm success rate, NOT an individual ",
    "patient trajectory. Sister models from the same paper: ",
    "modellib('Takechi_2025_nemolizumab') (popPK), ",
    "modellib('Takechi_2025_nemolizumab_ppnrs') (population PD for the weekly ",
    "average PP-NRS score) and ",
    "modellib('Takechi_2025_nemolizumab_mbma_iga') (the IGA-endpoint MBMA)."
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
    "Parameter values: Table 4, 'PP-NRS' column.",
    "Included trials, arm sizes and the per-arm IGA severity split: Table 3.",
    "Covariate functional form: reconstructed from Fig. 4A and 4C - the paper",
    "never prints the covariate equation. See the e_igamod_edrug_nemolizumab",
    "note in ini().",
    sep = " "
  )

  vignette <- "Takechi_2025_nemolizumab"

  # The single random effect is a BETWEEN-STUDY effect added directly to the
  # logit-scale treatment effect EFF, not an IIV on a named fixed-effect
  # parameter, so it cannot follow the popPK eta<transformed-param-name>
  # convention. Declared paper-specific per the Goteti_2024_SLE_mbma
  # precedent, and matching the IGA sibling of this same paper.
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
      "fraction/arm (probability that a patient in the study arm achieves an ",
      "improvement of at least 4 points in the weekly average Peak Pruritus ",
      "Numerical Rating Scale; this output is NOT a drug concentration - the ",
      "slash satisfies checkModelConventions parsing)"
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
        "dupilumab trials also report the Worst Itch NRS rather than the ",
        "PP-NRS; the source treats these as conceptually equivalent while ",
        "cautioning against over-interpretation (Discussion, final ",
        "limitation). The covariate below is NOT applied to dupilumab - ",
        "Table 4 attaches it to Edrug,nemolizumab only."
      ),
      source_name        = "Drug / Dose (Takechi 2025 Table 3, Dupilumab rows)"
    ),
    IGA_MOD_PCT = list(
      description        = paste0(
        "Percentage (0-100) of the study arm whose baseline Investigator's ",
        "Global Assessment was 3 (moderate) rather than 4 (severe). The ",
        "paper's 'PPM'."
      ),
      units              = "%",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "ARM-LEVEL AGGREGATE, not a per-subject characteristic - this is a ",
        "property of the trial arm being simulated. Table 3 prints the ",
        "COMPLEMENT ('IGA 4 at baseline (%)', i.e. the SEVERE share), so ",
        "IGA_MOD_PCT = 100 - that column. Across the 13 arms of Table 3 the ",
        "values are 45.5 47.0 47.4 52.7 53.9 56.3 59.0 60.5 61.0 62.8 64.6 ",
        "70.7 72.0 percent, whose median is exactly 59.0 - the centring ",
        "constant used in model(). The model consumes it as a FRACTION ",
        "(IGA_MOD_PCT / 100), converted inline; the column itself carries a ",
        "percent, matching TUMTP_SQUAM_PCT, PS_ECOG_0_PCT and RACE_ASIAN_PCT. ",
        "A LOWER IGA_MOD_PCT (a more severe arm) gives a LARGER nemolizumab ",
        "effect, because the coefficient is negative and the deviation from ",
        "the median is then negative. Simulating far outside the observed ",
        "45.5-72.0 percent range is extrapolation: the multiplier ",
        "(1 + e * (PPM - 0.590)) reaches zero at PPM = 1.133 and is therefore ",
        "positive over the whole [0, 1] domain, but the source calibrates it ",
        "only over the printed range."
      ),
      source_name        = "IGA4 at baseline (%) (Takechi 2025 Table 3), complemented to 100"
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
      "results'. That severity split is nevertheless the one covariate this ",
      "endpoint retains - see IGA_MOD_PCT."
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
      "arm's reported PP-NRS success rate at one time point, digitised from ",
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
      "The covariates screened for this endpoint were age, sex, race, body ",
      "weight, baseline PP-NRS and IGA scores; only the arm's moderate-IGA ",
      "proportion was retained, and only on the nemolizumab maximum effect."
    )
  )

  ini({
    # ========================================================================
    # All values are Takechi 2025 Table 4, 'PP-NRS' column.
    #
    # STRUCTURE (Methods 'Model-Based Meta-analysis'):
    #   Pr(PP-NRS) = 1 / (1 + exp(-EFF))
    #   EFF        = f0 + fdrug + eta_study
    #   f0         = A + Pmax  * (1 - exp(-Kp * time))
    #   fdrug      = Edrug * (1 - exp(-Kd * time))
    # with Edrug and Kd estimated separately per drug, and the nemolizumab
    # Edrug carrying the arm-level covariate below. Time is in WEEKS: the rate
    # constants are printed as "(/week)" in Table 4 and the source's Figs. 3
    # and 4 use a 0-24 week axis.
    #
    # VARIANCE SCALE. Table 4 gives an estimate and an %RSE per row and no CV%
    # label anywhere, so the 'ISV on EFF' and 'Additive error' rows are raw
    # NONMEM variances, matching the convention proved for Table 2 in the
    # companion popPK model file. Both are converted to SDs below.
    # ========================================================================

    e0 <- -4.38
    label("Intercept of the placebo effect on the logit scale at time zero (logit units); expit(-4.38) = 1.24% PP-NRS success at baseline")
    # Table 4 PP-NRS column, row A: -4.38 (%RSE 13.6).

    pmax_placebo <- 2.64
    label("Maximum placebo effect on the logit scale (logit units); the plateau placebo PP-NRS success rate is expit(-4.38 + 2.64) = 14.9%")
    # Table 4 PP-NRS column, row Pmax: 2.64 (%RSE 21.7).

    lkp <- log(0.277)
    label("Log rate constant for the onset of the placebo effect (log 1/week); back-transform Kp = 0.277 1/week, onset half-life 2.5 weeks")
    # Table 4 PP-NRS column, row Kp (/week): 0.277 (%RSE 31.4). Note this is
    # more than twice the IGA model's 0.123 /week: the itch placebo response
    # sets in considerably faster than the skin-lesion placebo response.

    edrug_nemolizumab <- 1.63
    label("Maximum nemolizumab effect on the logit scale at the reference arm severity (logit units)")
    # Table 4 PP-NRS column, row Edrug,nemolizumab: 1.63 (%RSE 4.5). This is
    # the value AT the covariate reference (IGA_MOD_PCT = 59.0), not an
    # unconditional maximum - see e_igamod_edrug_nemolizumab.

    e_igamod_edrug_nemolizumab <- -1.84
    label("Fractional change in the maximum nemolizumab effect per unit increase in the arm's moderate-IGA fraction (IGA_MOD_PCT/100) away from 0.590 (unitless)")
    # Table 4 PP-NRS column, row 'PPM effect on Edrug,nemolizumab': -1.84
    # (%RSE 17.5). The IGA column prints a dash for this row, i.e. the
    # covariate was tested for that endpoint and NOT retained.
    #
    # FUNCTIONAL FORM IS A RECONSTRUCTION, NOT A QUOTATION. The paper never
    # prints the covariate equation for the MBMA. The form encoded in model()
    # is a fractional deviation centred on the median arm:
    #   Edrug,nemo = 1.63 * (1 + (-1.84) * (IGA_MOD_PCT/100 - 0.590))
    # It was recovered by arithmetic against the source's own Fig. 4, and the
    # two competing readings are falsified at the same three points. Fig. 4's
    # caption states the simulation used "the recommended dose of each drug in
    # Japan", i.e. the Yokozeki 30 mg Q4W arm, whose IGA_MOD_PCT is
    # 100 - 54.5 = 45.5. Substituting gives Edrug = 1.63 * 1.2484 = 2.035 and
    # predicted success rates of 35.9 / 56.5 / 57.2 percent at weeks 4 / 16 /
    # 24, against the ~36 and ~56 percent week-4 and week-16 medians of the
    # Fig. 4C forest plot and the ~57 percent plateau of Fig. 4A. The
    # alternatives miss all three: a
    # normalised POWER function - the form the paper's Methods states for the
    # PopPD analysis - gives 50.4 / 70.2 / 70.8, and a plain additive shift on
    # the logit scale gives 32.4 / 52.7 / 53.4. With NO covariate at all the
    # model gives 27.2 / 46.5 / 47.2, which is the ~10-point miss that makes
    # the covariate load-bearing. The same parameters reproduce the placebo
    # arm (14.9% at week 24 vs ~15% plotted) and the dupilumab arm (46.8% at
    # week 16, 51.9% at week 24 vs ~47% and ~52% plotted) with no covariate,
    # so the nemolizumab gap is specifically this term. The vignette runs all
    # of these as regression tests.

    edrug_dupilumab <- 1.93
    label("Maximum dupilumab effect on the logit scale (logit units)")
    # Table 4 PP-NRS column, row Edrug,dupilumab: 1.93 (%RSE 13.1). No
    # covariate is attached to this parameter in Table 4.

    lkdrug_nemolizumab <- log(1.93)
    label("Log rate constant for the onset of the nemolizumab effect (log 1/week); back-transform Kd = 1.93 1/week, onset half-life 0.36 weeks")
    # Table 4 PP-NRS column, row Kd,nemolizumab (/week): 1.93 (%RSE 34.4).
    # Essentially complete within the first fortnight, which is the source's
    # quantitative basis for its rapid-itch-relief claim.

    lkdrug_dupilumab <- log(0.119)
    label("Log rate constant for the onset of the dupilumab effect (log 1/week); back-transform Kd = 0.119 1/week, onset half-life 5.8 weeks")
    # Table 4 PP-NRS column, row Kd,dupilumab (/week): 0.119 (%RSE 66.7).
    # Sixteen times slower than nemolizumab's, which is a far wider separation
    # than the IGA endpoint shows (0.416 vs 0.134, a factor of three).

    # ---- Between-STUDY variability ------------------------------------------
    # Methods: "The between-study variability is described by eta, which was
    # assumed to follow a normal distribution (mean = 0, variance = omega^2)."
    # ONE eta, added directly to the logit-scale EFF, shared by every arm of a
    # study. Named eta_study_* per the MBMA convention in this package so it
    # cannot be mistaken for a between-SUBJECT eta - there is no
    # between-subject variability in this model at all.
    eta_study_eff ~ 0.105
    # Table 4 PP-NRS column, row 'ISV on EFF': 0.105 (%RSE 73.4). Variance, so
    # the between-study SD on the logit scale is sqrt(0.105) = 0.324 - less
    # than half the IGA model's 0.475, consistent with the covariate having
    # absorbed some of the between-study heterogeneity for this endpoint.

    # ---- Residual error -----------------------------------------------------
    addSd <- 0.859069
    label("Residual standard deviation multiplier on the arm-level PP-NRS success rate; the per-observation SD is addSd * sqrt(Pr * (1 - Pr) / N_arm)")
    # Table 4 PP-NRS column, row 'Additive error': 0.738 (%RSE 8.5), read as a
    # VARIANCE, so the multiplier is sqrt(0.738) = 0.859069.
    #
    # ARM-SIZE WEIGHTING. The source's displayed residual equation is
    #   success rate (%) = Pr + eps * (Pr * (1 - Pr) / N)^(1/2)
    # i.e. the residual SD is this multiplier times the BINOMIAL standard
    # error of an arm of N subjects. model() encodes the prediction-dependent
    # factor sqrt(Pr * (1 - Pr)) exactly, but N is a property of the study arm
    # being simulated rather than of the model, so the 1 / sqrt(N) factor is
    # left to the user - the same convention as Boucher_2018_naproxen_mbma and
    # as this paper's IGA sibling. The model as shipped therefore emits the
    # N = 1 residual; DIVIDE THE SIMULATED RESIDUAL BY sqrt(N_arm), or
    # equivalently set addSd to 0.859069 / sqrt(N_arm), to reproduce an arm of
    # N subjects. The derived variable armSeUnit in model() exposes
    # sqrt(Pr * (1 - Pr)) so the rescaling can be done directly from the solve
    # output. Arm sizes for the six included trials are in Table 3 and are
    # repeated in population$notes.
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

    # ---- 3. Arm-severity covariate on the nemolizumab maximum ----------------
    # The register stores the moderate-IGA share as a PERCENT; the covariate
    # model is written on the FRACTION, so convert inline (same idiom as
    # Franzese_2026_pdl1_nsclc_mbma).
    ppmFrac <- IGA_MOD_PCT / 100

    # Fractional deviation from the median arm of Table 3 (59.0 percent
    # moderate). Reconstructed form - see the e_igamod_edrug_nemolizumab note
    # in ini() for the arithmetic that selects it over a power function and
    # over an additive logit shift.
    edrugNemolizumab <-
      edrug_nemolizumab * (1 + e_igamod_edrug_nemolizumab * (ppmFrac - 0.590))

    # ---- 4. Placebo and drug effects on the logit scale ----------------------
    # f0 = A + Pmax * (1 - exp(-Kp * time)); present in every arm including
    # the active arms, which is why an active arm's plateau is the sum of the
    # placebo plateau and its own drug maximum.
    f0 <- e0 + pmax_placebo * (1 - exp(-kp * time))

    # fdrug = Edrug * (1 - exp(-Kd * time)); at most one term is non-zero
    # because no included trial combined the two drugs.
    fdrug <-
      edrugNemolizumab * (1 - exp(-kdNemolizumab * time)) * onNemolizumab +
      edrug_dupilumab  * (1 - exp(-kdDupilumab   * time)) * onDupilumab

    # ---- 5. Logit transform (Methods 'Model-Based Meta-analysis') ------------
    # The between-study eta is added on the logit scale, before the transform.
    eff   <- f0 + fdrug + eta_study_eff
    pResp <- expit(eff)

    # ---- 6. Observation and error -------------------------------------------
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
    # `observationVar`), NOT because it is a concentration - it is a PP-NRS
    # success-rate probability. Same convention as `Goteti_2024_SLE_mbma.R`
    # and as this paper's IGA sibling.
    Cc <- pResp
    Cc ~ prop(propBinom)
  })
}
