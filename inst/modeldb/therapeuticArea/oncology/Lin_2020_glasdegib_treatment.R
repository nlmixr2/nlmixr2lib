Lin_2020_glasdegib_treatment <- function() {
  description <- "Parametric exponential time-to-event model for overall survival in older adults (>= 55 y) with newly diagnosed AML who were ineligible for intensive chemotherapy in BRIGHT AML 1003 Phase 2 (Lin 2020). Treatment-response analysis: hazard h(t) = lambda * (1 + theta_ldac_alone * LDAC_alone), with lambda the daily death hazard for glasdegib + LDAC (reference) and theta_ldac_alone the proportional increase in hazard for LDAC alone. Treatment arm was the only covariate retained after SCM backward elimination (alpha < 0.001); demographics, baseline safety labs, and disease characteristics were not significant. Original NONMEM run was fit with no estimated IIV and no residual error (likelihood-based parametric TTE)."
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
    dosing        = "n/a (no drug-dosing events; the treatment arm enters via the CONMED_GLASDEGIB indicator)",
    concentration = "probability (the model output `sur` is a survival probability, not a drug concentration)"
  )

  covariateData <- list(
    CONMED_GLASDEGIB = list(
      description        = "Binary treatment-arm indicator: 1 = subject is in the glasdegib + LDAC arm, 0 = subject is in the LDAC alone arm.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (LDAC alone arm; LDAC = low-dose cytarabine, the standard-of-care comparator in BRIGHT AML 1003 Phase 2)",
      notes              = "Time-fixed per subject. The paper's published equation uses the complementary indicator LDAC_alone = 1 if LDAC alone, 0 otherwise. Inside model() the canonical column is mapped via LDAC_alone = 1 - CONMED_GLASDEGIB so the published parameters lambda = 0.00253 and theta_ldac_alone = 1.376 are preserved verbatim. BRIGHT AML 1003 Phase 2 cohort split: 78 glasdegib + LDAC (CONMED_GLASDEGIB = 1) and 38 LDAC alone (CONMED_GLASDEGIB = 0). Source dataset column not published; derive from the protocol-defined treatment arm assignment.",
      source_name        = "TRT"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Baseline age (years) at randomisation.",
      units              = "years",
      type               = "continuous",
      notes              = "Tested in SCM forward inclusion (alpha < 0.05); not retained in the final treatment-response model (median 76 years overall, range 58-92, per Table 1)."
    ),
    WT = list(
      description        = "Baseline body weight (kg) at randomisation.",
      units              = "kg",
      type               = "continuous",
      notes              = "Tested in SCM forward inclusion; not retained in the final model (median 78.2 kg overall, range 47.5-118.0, per Table 1)."
    ),
    SEXF = list(
      description        = "Female sex indicator (1 = female, 0 = male).",
      units              = "(binary)",
      type               = "binary",
      notes              = "Tested in SCM forward inclusion; not retained in the final model (29% female overall, per Table 1)."
    ),
    RACE_BLACK = list(
      description        = "Black race indicator (1 = Black, 0 = otherwise).",
      units              = "(binary)",
      type               = "binary",
      notes              = "Tested in SCM forward inclusion; not retained in the final model (1% Black overall, per Table 1)."
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator (1 = Asian, 0 = otherwise).",
      units              = "(binary)",
      type               = "binary",
      notes              = "Tested in SCM forward inclusion; not retained in the final model (2% Asian overall, per Table 1)."
    ),
    ECOG_GE1 = list(
      description        = "Binary baseline ECOG performance status indicator (1 = ECOG >= 1, 0 = ECOG = 0).",
      units              = "(binary)",
      type               = "binary",
      notes              = "Tested in SCM forward inclusion; not retained in the final treatment-response model (89% with ECOG >= 1 overall, per Table 1)."
    ),
    CRCL_COCKCROFT = list(
      description        = "Baseline creatinine clearance (mL/min) by Cockcroft-Gault.",
      units              = "mL/min",
      type               = "continuous",
      notes              = "Tested in SCM forward inclusion; not retained in the final model (median 62.7 mL/min, range 32.5-134.0, per Table 1; most subjects had mild renal impairment per KDOQI)."
    ),
    AST = list(
      description        = "Baseline aspartate transaminase (U/L).",
      units              = "U/L",
      type               = "continuous",
      notes              = "Tested in SCM forward inclusion; not retained in the final model (median 20.5 U/L, range 7.0-111.0, per Table 1)."
    ),
    WBC = list(
      description        = "Baseline white blood cell count (10^9 cells/L).",
      units              = "10^9 cells/L",
      type               = "continuous",
      notes              = "Tested in SCM forward inclusion; not retained in the final model (median 3.1, range 0.4-5850.0, per Table 1)."
    ),
    BMBLAST_PCT = list(
      description        = "Baseline percentage of bone marrow blasts (%).",
      units              = "%",
      type               = "continuous",
      notes              = "Tested in SCM forward inclusion; not retained in the final model (median 44.0%, range 13-99, per Table 1)."
    ),
    PERIPH_BLAST_PCT = list(
      description        = "Baseline percentage of peripheral blasts (%).",
      units              = "%",
      type               = "continuous",
      notes              = "Tested in SCM forward inclusion; not retained in the final model (median 7.0%, range 0-91, per Table 1)."
    ),
    AML_SECONDARY = list(
      description        = "Secondary (vs. de novo) AML indicator (1 = secondary disease, 0 = de novo).",
      units              = "(binary)",
      type               = "binary",
      notes              = "Tested in SCM forward inclusion; not retained in the final model (52% secondary overall, per Table 1)."
    ),
    CYTO_POOR = list(
      description        = "Poor cytogenetic risk indicator (1 = poor, 0 = good/intermediate).",
      units              = "(binary)",
      type               = "binary",
      notes              = "Significant in SCM forward inclusion (alpha < 0.05) but eliminated during backward elimination (alpha < 0.001); not retained in the final model (40% poor risk overall, per Table 1)."
    ),
    PRIOR_HMA = list(
      description        = "Prior hypomethylating-agent treatment indicator (1 = prior HMA, 0 = no prior HMA).",
      units              = "(binary)",
      type               = "binary",
      notes              = "Tested in SCM forward inclusion; not retained in the final model (15% with prior HMA overall, per Table 1)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 116L,
    n_studies      = 1L,
    age_range      = "58-92 years (median 76; BRIGHT AML 1003 Phase 2 inclusion required age >= 55 years and ineligibility for intensive chemotherapy)",
    weight_range   = "47.5-118.0 kg (median 78.2)",
    sex_female_pct = 29,
    race_ethnicity = c(White = 97, Black = 1, Asian = 2),
    disease_state  = "newly diagnosed, previously untreated acute myeloid leukemia (WHO 2008) or high-risk MDS in adults ineligible for intensive chemotherapy",
    dose_range     = "Glasdegib 100 mg orally once daily in 28-day cycles + LDAC (low-dose cytarabine) 20 mg subcutaneously BID for 10 days per 28-day cycle, vs. LDAC alone (same schedule); follow-up up to 4 years",
    regions        = "Multicenter (Phase 1b/2 trial sites; NCT01546038)",
    notes          = "116 patients from BRIGHT AML 1003 Phase 2 (NCT01546038): 78 randomised 2:1 to glasdegib + LDAC and 38 to LDAC alone. Baseline demographics and safety labs per Lin 2020 Table 1; data cut-off 3 January 2017."
  )

  ini({
    # Exponential time-to-event hazard for overall survival (BRIGHT AML 1003
    # Phase 2 treatment-response analysis; n = 116). Final estimates from
    # Lin 2020 Results / page 4: "the base hazard (relative standard error
    # [RSE]) for glasdegib + LDAC was estimated at 0.00253 (13.82%). LDAC
    # alone treatment resulted in an ~138% increase (multiply 1.376 by 1,
    # or by 0 if not LDAC alone treatment) in base hazard (RSE, 37.74%)."
    # The reported equation (page 4 / Fig. 1b annotation):
    #   S(t) = exp(-0.00253 * (1 + 1.376 * LDAC_alone) * t),  t in days.
    # Fig. 1b also annotates S(t) = exp(-0.00601 * t) for the LDAC alone
    # arm, which cross-checks 0.00253 * (1 + 1.376) = 0.00601 /day.
    # The corresponding hazard ratio glasdegib+LDAC vs LDAC alone is
    # 0.00253 / 0.00601 = 0.421, matching the paper's HR 0.42 (95% CI
    # 0.28-0.66).
    llam_haz       <- log(0.00253); label("Exponential baseline hazard for overall survival in the glasdegib + LDAC reference arm (lambda, log(1/day))")  # Lin 2020 Results page 4: lambda = 0.00253, RSE 13.82%
    e_ldac_alone   <- 1.376;        label("Proportional increase in hazard for LDAC alone vs glasdegib + LDAC (theta_ldac_alone, unitless multiplier on lambda)")  # Lin 2020 Results page 4: 1.376, RSE 37.74%; produces HR 0.42 glasdegib+LDAC vs LDAC alone

    # No estimated IIV. The treatment-response TTE model was fit with no
    # subject-level random effects (parametric population TTE on a binary
    # study-arm covariate). No eta* parameters are added.

    # No residual error. The NONMEM $ERROR / $EST framework for parametric
    # TTE uses the LIKE / F_FLAG branch (the likelihood is the survival /
    # event-density itself, not an observation-error model). For forward
    # simulation in nlmixr2lib, `sur` and `hazard` are exposed as derived
    # outputs.
  })
  model({
    # Map the canonical CONMED_GLASDEGIB indicator (1 = glasdegib + LDAC,
    # 0 = LDAC alone) to the paper's "LDAC_alone" indicator so the
    # published parameter values (lambda = 0.00253, theta_ldac_alone =
    # 1.376) remain verbatim.
    ldac_alone <- 1 - CONMED_GLASDEGIB

    # Constant exponential baseline hazard with multiplicative arm effect.
    lam_haz <- exp(llam_haz)
    hazard  <- lam_haz * (1 + e_ldac_alone * ldac_alone)

    # Cumulative hazard and survival probability. Exponential hazard means
    # cumhaz(t) = hazard * t analytically; the ODE form is used for
    # consistency with other nlmixr2lib TTE models (see
    # NA_NA_tte_lognormal, Zecchin_2016_survival).
    d/dt(cumhaz) <- hazard
    cumhaz(0)    <- 0
    sur          <- exp(-cumhaz)
  })
}
