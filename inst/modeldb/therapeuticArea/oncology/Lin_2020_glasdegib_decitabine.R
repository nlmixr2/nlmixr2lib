Lin_2020_glasdegib_decitabine <- function() {
  description <- "Parametric exponential time-to-event model for overall survival in the exploratory pooled BRIGHT AML 1003 Phase 1b + Phase 2 analysis (Lin 2020): older adults with newly diagnosed AML or high-risk MDS ineligible for intensive chemotherapy across three treatment arms (LDAC alone, glasdegib + LDAC, glasdegib + decitabine; n = 162). Hazard h(t) = lambda * (1 - theta_gl_ldac * GL_LDAC_arm) * (1 - theta_gl_dec * GL_DEC_arm), where lambda is the daily death hazard for the LDAC alone reference arm and the two arm indicators encode the proportional hazard reduction associated with each glasdegib-containing combination. Treatment arm was the only covariate retained after SCM backward elimination (alpha < 0.001) from a full model that also included log-transformed baseline percentage of bone marrow blasts and prior hypomethylating-agent use. Original NONMEM run was fit with no estimated IIV and no residual error (likelihood-based parametric TTE)."
  reference <- paste(
    "Lin S, Shaik N, Chan G, Cortes JE, Ruiz-Garcia A.",
    "An evaluation of overall survival in patients with newly diagnosed",
    "acute myeloid leukemia and the relationship with glasdegib treatment",
    "and exposure.",
    "Cancer Chemother Pharmacol. 2020;86(4):451-459.",
    "doi:10.1007/s00280-020-04132-x.",
    "BRIGHT AML 1003: ClinicalTrials.gov NCT01546038.",
    sep = " "
  )
  vignette <- "Lin_2020_glasdegib_AML_overall_survival"
  units <- list(
    time          = "day",
    dosing        = "n/a (no drug-dosing events; treatment arm enters via the CONMED_GLASDEGIB and CONMED_DECITABINE indicators)",
    concentration = "probability (the model output `sur` is a survival probability, not a drug concentration)"
  )

  covariateData <- list(
    CONMED_GLASDEGIB = list(
      description        = "Binary indicator: 1 = subject is in a glasdegib-containing arm (glasdegib + LDAC or glasdegib + decitabine), 0 = subject is in the LDAC alone arm.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (LDAC alone arm; reference category in this exploratory three-arm analysis)",
      notes              = "Time-fixed per subject. Paired with CONMED_DECITABINE to distinguish the two glasdegib-containing arms. Cohort split for this analysis (Lin 2020 exploratory analysis page 6): LDAC alone (CONMED_GLASDEGIB = 0, CONMED_DECITABINE = 0) n = 38 (Phase 2); glasdegib + LDAC (CONMED_GLASDEGIB = 1, CONMED_DECITABINE = 0) n = 117 (Phase 2 + Phase 1b Arm A); glasdegib + decitabine (CONMED_GLASDEGIB = 1, CONMED_DECITABINE = 1) n = 7 (Phase 1b Arm B: 5 AML + 2 MDS). Sum 38 + 117 + 7 = 162.",
      source_name        = "TRT"
    ),
    CONMED_DECITABINE = list(
      description        = "Binary indicator: 1 = subject is in the glasdegib + decitabine arm, 0 = subject is in either the LDAC alone arm or the glasdegib + LDAC arm.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (LDAC backbone, i.e. either LDAC alone or glasdegib + LDAC)",
      notes              = "Time-fixed per subject. Paired with CONMED_GLASDEGIB. Subjects in the glasdegib + decitabine arm (n = 7) have CONMED_GLASDEGIB = 1 and CONMED_DECITABINE = 1. In Lin 2020 the only source of decitabine exposure is BRIGHT AML 1003 Phase 1b Arm B (glasdegib + decitabine combination); there is no decitabine-alone arm in the trial.",
      source_name        = "TRT"
    )
  )

  covariatesDataExcluded <- list(
    LOG_BMBLAST_PCT = list(
      description        = "Log-transformed baseline percentage of bone marrow blasts (%).",
      units              = "log(%)",
      type               = "continuous",
      notes              = "Entered the SCM full model during forward inclusion (alpha < 0.05) but did not survive backward elimination (alpha < 0.001) for the exploratory three-arm analysis (page 6)."
    ),
    PRIOR_HMA = list(
      description        = "Prior hypomethylating-agent treatment indicator (1 = prior HMA, 0 = no prior HMA).",
      units              = "(binary)",
      type               = "binary",
      notes              = "Entered the SCM full model during forward inclusion but did not survive backward elimination for the exploratory three-arm analysis (page 6)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 162L,
    n_studies      = 1L,
    age_range      = "Pooled Phase 1b + Phase 2 BRIGHT AML 1003 cohorts; trial inclusion required age >= 55 years. Per-cohort age distributions per Lin 2020 Table 1 (Phase 2) and Online Resource 2 (full pooled cohort).",
    weight_range   = "Per-cohort weight distributions per Lin 2020 Table 1 (Phase 2 median 78.2 kg, range 47.5-118.0) and Online Resource 2 (full pooled cohort).",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = "newly diagnosed AML (WHO 2008) or high-risk MDS in adults ineligible for intensive chemotherapy. The MDS subset enters via the two Phase 1b Arm B subjects.",
    dose_range     = "LDAC alone (LDAC 20 mg SC BID for 10 days per 28-day cycle, n = 38, Phase 2); Glasdegib 100 mg QD orally + LDAC (n = 117, Phase 2 randomised + Phase 1b Arm A which evaluated 100 or 200 mg QD); Glasdegib 100 or 200 mg QD orally + decitabine (n = 7, Phase 1b Arm B). Glasdegib dose reductions allowed for adverse-event management.",
    regions        = "Multicenter (NCT01546038)",
    notes          = "162 patients pooled across BRIGHT AML 1003 phases: Phase 2 AML cohort (78 glasdegib + LDAC, 38 LDAC alone) plus Phase 1b Arm A (glasdegib + LDAC AML) and Arm B (5 AML + 2 MDS on glasdegib + decitabine). Per the paper, the Phase 1b Arm B sample is small and the CI of the glasdegib + decitabine hazard reduction is wide (-95.0% to -28.6%). Demographic and baseline-characteristic table is published as Online Resource 2 (supplementary; not on disk for this extraction)."
  )

  ini({
    # Exponential time-to-event hazard for overall survival in the
    # exploratory pooled Phase 1b + Phase 2 BRIGHT AML 1003 cohort
    # (n = 162). Final estimates from Lin 2020 Results page 6 and the
    # paper's typeset equation (page 6, equation rendered as
    # "S(t) = exp(-(0.00540 * (1 - 0.480 * glasdegib+LDAC) *
    #                       (1 - 0.618 * glasdegib+decitabine) * t)"):
    #   * Base hazard (LDAC alone reference): 0.00540, RSE 14.1%.
    #   * Glasdegib + LDAC arm: hazard fraction 1 - 0.480 = 0.520 of
    #     LDAC alone.
    #   * Glasdegib + decitabine arm: hazard fraction 1 - 0.618 = 0.382
    #     of LDAC alone, i.e. ~61.8% lower than LDAC alone (page 6).
    # Cross-check on median OS for glasdegib + decitabine:
    #   median OS = ln(2) / (0.00540 * 0.382) = log(2) / 0.002063 = 336
    #   days = 11.0 months, matching the paper's "estimated median OS of
    #   11.1 months for glasdegib + decitabine" (page 6).
    llam_haz       <- log(0.00540); label("Exponential baseline hazard for overall survival in the LDAC alone reference arm (lambda, log(1/day))")  # Lin 2020 Results page 6: lambda = 0.00540, RSE 14.1%
    e_gl_ldac_haz  <- 0.480;        label("Proportional hazard reduction for glasdegib + LDAC vs LDAC alone (theta_gl_ldac, unitless multiplier; hazard fraction = 1 - 0.480 = 0.520)")  # Lin 2020 page 6 equation: 0.480 = 48.0% hazard reduction
    e_gl_dec_haz   <- 0.618;        label("Proportional hazard reduction for glasdegib + decitabine vs LDAC alone (theta_gl_dec, unitless multiplier; hazard fraction = 1 - 0.618 = 0.382)")  # Lin 2020 page 6 equation: 0.618 = 61.8% hazard reduction (95% CI -95.0% to -28.6%)

    # No estimated IIV. Parametric population TTE fit with no subject-
    # level random effects on a small (n = 162) cohort with two binary
    # treatment-arm indicators.

    # No residual error. NONMEM LIKE / F_FLAG branch for parametric TTE.
  })
  model({
    # Map the two canonical drug indicators into the paper's two binary
    # arm flags so the published parameters (0.00540, 0.480, 0.618) are
    # preserved verbatim. The encoding is:
    #   LDAC alone               -> CONMED_GLASDEGIB = 0, CONMED_DECITABINE = 0
    #   Glasdegib + LDAC         -> CONMED_GLASDEGIB = 1, CONMED_DECITABINE = 0
    #   Glasdegib + decitabine   -> CONMED_GLASDEGIB = 1, CONMED_DECITABINE = 1
    gl_ldac_arm <- CONMED_GLASDEGIB * (1 - CONMED_DECITABINE)
    gl_dec_arm  <- CONMED_GLASDEGIB * CONMED_DECITABINE

    # Constant exponential baseline hazard with multiplicative arm effects.
    lam_haz <- exp(llam_haz)
    hazard  <- lam_haz *
      (1 - e_gl_ldac_haz * gl_ldac_arm) *
      (1 - e_gl_dec_haz  * gl_dec_arm)

    # Cumulative hazard and survival probability.
    d/dt(cumhaz) <- hazard
    cumhaz(0)    <- 0
    sur          <- exp(-cumhaz)
  })
}
