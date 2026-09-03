Patel_2020_sapanisertib_QTcI <- function() {
  description <- paste(
    "Concentration-DeltaQTcI linear mixed-effects PD model for the mTORC1/2",
    "inhibitor sapanisertib (MLN0128 / INK128 / TAK-228) in 44 patients with",
    "advanced solid tumors, from a dedicated QTc (DQT) study of a single",
    "40 mg oral dose. The endpoint is change from time-matched day -1",
    "baseline in the individually rate-corrected QT interval",
    "(DeltaQTcI, msec). The model is:",
    "DeltaQTcI = e0 + slope * CP_SAPANISERTIB_NGML,",
    "with e0 = -0.304 msec and slope = -0.033 msec per ng/mL. The slope is",
    "NEGATIVE and statistically significant (P = 0.002): sapanisertib",
    "concentration is associated with QTcI SHORTENING, not prolongation,",
    "which is the paper's central negative finding.",
    "PD-only model: sapanisertib plasma concentration is supplied as a",
    "time-varying covariate CP_SAPANISERTIB_NGML (ng/mL). The source",
    "publication does not fit a population PK model -- the PK was",
    "characterised by NCA only (geometric mean Cmax 297 ng/mL, AUCinf",
    "2480 ng*h/mL, CL/F 16.1 L/h, t1/2 9.5 h after 40 mg); users wishing to",
    "drive the PD model from a simulated PK source must supply their own",
    "concentration trajectory (no sapanisertib popPK model exists in the",
    "nlmixr2lib registry). Companion model files",
    "Patel_2020_sapanisertib_QTcF.R and Patel_2020_sapanisertib_RR.R report",
    "the same linear structure for the Fridericia-corrected QT interval and",
    "for the RR interval, respectively.",
    sep = " "
  )

  reference <- paste(
    "Patel C, Goel S, Patel MR, Rangachari L, Wilbur JD, Shou Y,",
    "Venkatakrishnan K, Lockhart AC. (2020). Phase 1 Study to Evaluate the",
    "Effect of the Investigational Anticancer Agent Sapanisertib on the QTc",
    "Interval in Patients With Advanced Solid Tumors.",
    "Clinical Pharmacology in Drug Development 9(7):876-888.",
    "doi:10.1002/cpdd.808.",
    sep = " "
  )

  vignette <- "Patel_2020_sapanisertib"

  units <- list(
    time          = "h",
    dosing        = "(none; PD-only model fed by an external sapanisertib plasma-concentration covariate)",
    concentration = "(observation QTcI is the change from time-matched day -1 baseline in the individually rate-corrected QT interval, msec; driving covariate CP_SAPANISERTIB_NGML is in ng/mL)"
  )

  covariateData <- list(
    CP_SAPANISERTIB_NGML = list(
      description        = "Instantaneous sapanisertib plasma concentration at the time of each PD observation, supplied as a time-varying covariate from observed plasma samples or an upstream PK source.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying per event row. Drives the linear concentration-",
        "DeltaQTcI expression DeltaQTcI = e0 + slope *",
        "CP_SAPANISERTIB_NGML.",
        "The slope is reported directly in msec per ng/mL (Patel 2020",
        "Figure 4A in-panel parameter table, x-axis 'Concentration (ng",
        "mL-1)'), so CP_SAPANISERTIB_NGML is supplied in ng/mL with no",
        "in-model unit rescaling.",
        "In Patel 2020 this was the observed sapanisertib plasma",
        "concentration determined by a validated LC-MS/MS method with",
        "deuterated sapanisertib as internal standard (dynamic range",
        "1-1000 ng/mL; accuracy -4.0% to -2.7% bias; precision 3.0-9.5%",
        "CV; Patel 2020 Methods 'Electrocardiogram, PK, and Safety",
        "Assessments'), drawn time-matched to the triplicate ECGs at 0",
        "(<=15 min predose), 0.25, 0.5, 1, 1.5, 2, 2.5, 3, 4, 6, 8, 10,",
        "24 and 48 h after the single 40 mg oral dose.",
        "Reference values observed: geometric mean Cmax 297 ng/mL (%CV",
        "52.9) after 40 mg, median 317 (range 94-976) ng/mL; median Tmax",
        "1.52 h (range 0.50-24.00) (Patel 2020 Table S2). The",
        "concentration-QTc dataset spans roughly 0-976 ng/mL (Patel 2020",
        "Figure 4A x-axis extent). Geometric mean Cmax at the two lower",
        "clinical doses, used by the paper as projection anchors, are",
        "36.9 ng/mL (4 mg once daily) and 235 ng/mL (30 mg once weekly),",
        "both from prior single-agent studies (Patel 2020 Table 4",
        "footnotes a and b).",
        "Set to 0 outside the drug-exposure window (the concentration-",
        "slope term then collapses to the intercept e0)."
      ),
      source_name        = "sapanisertib plasma concentration"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 44L,
    n_studies        = 1L,
    age_range        = "22-79 years",
    age_median       = "59.5 years",
    weight_range     = "49-122 kg",
    weight_median    = "70.6 kg",
    sex_female_pct   = round(100 * 28 / 44, 1),
    race_ethnicity   = c(White = 68, Black = 20, `Not reported` = 11),
    disease_state    = paste(
      "Adults with a radiographically or clinically evaluable advanced",
      "solid tumor and ECOG performance status 0 (30%), 1 (68%) or 2",
      "(2%). 42 of 44 patients (95%) had received prior antineoplastic",
      "therapy. Disease stage was IV or higher in 34 of 44 (77%).",
      "Tumor types: colon and endometrial (5 each, 11% each), kidney (4,",
      "9%), colorectal and sarcoma (2 each, 5% each), and 26 (59%)",
      "other (anal, cervical, gastric, liver, melanoma, non-small cell",
      "lung, ovarian, small cell lung, soft tissue, stomach, thyroid).",
      "Entry required adequate hepatic, renal and hematologic function,",
      "fasting glucose <=120 mg/dL, fasting triglycerides <=300 mg/dL,",
      "and left ventricular ejection fraction within 5 absolute",
      "percentage points of the institutional normal. Patients with",
      "baseline QTcF > 430 msec (men) or > 450 msec (women), congenital",
      "long QT syndrome, diabetes mellitus, gastrointestinal disease, or",
      "significant active cardiovascular or pulmonary disease were",
      "EXCLUDED, so the cohort is pre-selected against QT-vulnerable",
      "patients (Patel 2020 Table 1 and Supplemental Methods",
      "'Patients')."
    ),
    dose_range       = paste(
      "Single 40 mg oral dose of sapanisertib (immediate-release",
      "capsules) on cycle 1 day 1 in the fasted state or >= 2 h after a",
      "light meal -- the single-agent maximum tolerated dose on the",
      "once-weekly schedule and a supratherapeutic exposure relative to",
      "the recommended phase 2 doses. From cycle 1 day 8 patients could",
      "optionally continue sapanisertib 30 mg once weekly in continuous",
      "28-day cycles for up to 12 months; only the single-dose day -1 to",
      "day 3 window contributes to the concentration-QTc analysis",
      "(Patel 2020 Methods 'Study Design and Patients')."
    ),
    regions          = "United States (5 centers: Montefiore Medical Center, Bronx NY; Washington University School of Medicine, St. Louis MO; Stephenson Cancer Center, University of Oklahoma, Oklahoma City OK; Mary Crowley Cancer Research Center, Dallas TX; Florida Cancer Specialists / Sarah Cannon Research Institute, Sarasota FL).",
    notes            = paste(
      "NCT02197572. Open-label, single-arm, dedicated QTc (DQT) study --",
      "a reduced design used in place of a thorough QT (TQT) study",
      "because the toxicity profile of an anticancer agent precludes",
      "dosing healthy volunteers and a placebo / positive-control arm",
      "(Patel 2020 Introduction).",
      "Time-matched baseline design: serial triplicate 12-lead ECGs were",
      "recorded on day -1 by Holter H12+ recorder at the same clock",
      "times as the day 1 postdose schedule, and DeltaQTcI was formed by",
      "subtracting the time-matched day -1 mean QTcI from the day 1 to",
      "day 3 mean QTcI. Each ECG collection was preceded by 5 minutes of",
      "supine rest; triplicates were extracted at 2- to 5-minute",
      "intervals and read centrally (Biotelemetry, Malvern, PA).",
      "Heart-rate correction method: individual (QTcI). Each patient's",
      "own slope b was estimated by linear regression of",
      "log(QT) = log(a) + b * log(RR) over all day -1 QT/RR pairs and",
      "used to form QTcI. QTcI was designated a priori as the primary",
      "correction method, supported by its smaller QTcI-RR residual",
      "slope (0.0068) versus QTcF-RR (0.0139) (Patel 2020 Methods",
      "'Electrocardiogram Analysis and End Points').",
      "ANALYSIS POPULATIONS DIFFER: all 44 patients contributed to the",
      "concentration-QTc analysis that this model file encodes, whereas",
      "the separate repeated-measures mixed-effects LSM analysis (Patel",
      "2020 Table 2) used only the 32 patients with a complete day -1 /",
      "day 1 Holter series. Do not gate this model against Table 2.",
      "Model-fitting software was NONMEM 7.3 (ICON Development",
      "Solutions) with post-processing in R; 95% CIs for the population",
      "parameters were constructed by nonparametric bootstrap (Patel",
      "2020 Methods 'Statistical Analyses'). Residual degrees of freedom",
      "for this DeltaQTcI fit were 409 (Patel 2020 Figure 4A in-panel",
      "table).",
      "The paper reports that the concentration-QTc models are mixed-",
      "effects models (patient as a random effect) but does NOT report",
      "any variance component -- see the ini() block comment and the",
      "vignette Errata."
    )
  )

  ini({
    # ==================================================================
    # Linear concentration-DeltaQTcI model. Both parameters are printed
    # inside the Figure 4A plot panel as a regression summary table
    # (columns: Parameter | Estimate | SE | df | t-value | P-value |
    # 95% Bootstrap CI). That in-panel table is a PRINTED parameter
    # table with full printed-value authority -- it is not a digitised
    # or curve-fitted value. It is invisible to pdftotext / docling
    # because Figure 4 is a raster image; it was read by rendering
    # page 10 of the PDF at 400 dpi.
    #
    #   Parameter | Estimate |    SE | df  | t-value | P-value | 95% boot CI
    #   Intercept |   -0.304 | 1.755 | 409 |  -0.173 |   0.863 | (-3.757, 3.177)
    #   Slope     |   -0.033 | 0.011 | 409 |  -3.137 |   0.002 | (-0.060, -0.013)
    #
    # Source equation (Patel 2020 Results 'Relationship of Sapanisertib
    # Concentration and QTc', Figure 4A):
    #   DeltaQTcI = -0.304 - 0.033 * C_sapanisertib(ng/mL)
    #
    # Transcription was verified by triple redundancy: every P-value
    # reproduces from its t-value and df via 2*pt(-|t|, df)
    # (0.863 -> 0.8627; 0.002 -> 0.00183), and both bootstrap CIs
    # contain their point estimates.
    #
    # INDEPENDENT CORROBORATION. Table 4 tabulates model-predicted
    # DeltaQTcI at three geometric-mean Cmax values (36.9, 235 and
    # 297 ng/mL -> -1.538, -8.159 and -10.231 msec). Three points
    # OVER-DETERMINE a two-parameter line, and they are collinear to
    # R^2 = 0.9999999997, so Table 4 inverts to intercept -0.30475 and
    # slope -0.033422. Both printed sources carry only three decimals,
    # so they are compared by overlap of their rounding intervals, not
    # by exact re-rounding: the slope -0.033422 does round to the
    # printed -0.033, but the intercept -0.30475 rounds to -0.305 and
    # NOT to the printed -0.304. That 0.00075 gap is fully inside the
    # inversion's own uncertainty -- rounding Table 4's three entries
    # to three decimals admits any intercept in [-0.3055, -0.3040]
    # (all 27 corner perturbations of +/-0.0005 enumerated), which
    # overlaps the panel's own [-0.3045, -0.3035]. The two independent
    # printed sources are therefore mutually consistent, and the
    # vignette asserts the interval overlap rather than an equality
    # that does not hold.
    #
    # This file uses the PRINTED 3-decimal panel estimates, so that
    # every value in ini() is a number printed in the paper. The
    # consequence is a <= 0.13 msec offset from Table 4's predictions
    # (e.g. -10.105 vs -10.231 msec at 297 ng/mL), which is entirely
    # slope-rounding: the exact back-solved slope -0.033422 reproduces
    # all three Table 4 entries to 1e-4 msec. Users who need to
    # reproduce Table 4 digit-for-digit should set slope = -0.033422.
    # Both facts are shown in the vignette source-trace section.
    # ==================================================================

    e0 <- -0.304
    label("Linear-mixed-effects intercept e0 on DeltaQTcI (msec)")
    # Patel 2020 Figure 4A in-panel table, 'Intercept' row: estimate
    # -0.304 (SE 1.755, df 409, t -0.173, P = 0.863, 95% bootstrap CI
    # -3.757 to 3.177). NOT statistically different from zero, as
    # expected for a change-from-time-matched-baseline endpoint with a
    # properly constructed day -1 baseline. Deliberately NOT log-
    # transformed: the intercept of a DeltaQTc model is negative here
    # and must be free to take either sign. Canonical bare PD baseline
    # parameter `e0`; same encoding as Mukker_2026_tuvusertib_QTcF.R
    # (e0 = -0.125) and Zhou_2025_fruquintinib_QTcF_M11.R (e0 = 2.37).

    slope <- -0.033
    label("Linear concentration-DeltaQTcI slope (msec per ng/mL)")
    # Patel 2020 Figure 4A in-panel table, 'Slope' row: estimate -0.033
    # (SE 0.011, df 409, t -3.137, P = 0.002, 95% bootstrap CI -0.060
    # to -0.013). Deliberately NOT log-transformed -- the slope is
    # NEGATIVE, so a log parameterisation is impossible. The whole
    # bootstrap CI lies below zero, so the QTcI-shortening association
    # is statistically significant and is a real feature of the fit,
    # not noise around zero. Same bare-`slope` encoding as
    # Mukker_2026_tuvusertib_QTcF.R and Zhou_2025_fruquintinib_*.R.
    # NOTE the printed t-value -3.137 is inconsistent with
    # -0.033 / 0.011 = -3.00 but consistent with the unrounded
    # -0.033422 / 0.010654; both displayed numbers are rounded to
    # 3 decimals, and this is the expected artefact rather than a
    # transcription error.

    # ==================================================================
    # Inter-individual variability: Patel 2020 Methods 'Statistical
    # Analyses' states that "Mixed-effects models were subsequently
    # developed to describe the direct effect of sapanisertib on change
    # from time-matched baseline in RR, QTcI, and QTcF", and the
    # residual df of 409 against 44 patients is consistent with a
    # per-subject random intercept. However the paper does NOT report
    # any variance estimate (no omega, no random-effect SD, no
    # shrinkage) for any of the three concentration-effect models. Per
    # the standing operator policy on unreported IIV, this file omits
    # the eta declarations and ships a typical-value-only model. The
    # vignette Errata documents the gap. The only surviving trace of
    # the variance components is the Table 4 threshold-exceedance
    # probability column and the Figure 4B 90% prediction interval;
    # neither is invertible to a unique variance triple, so no value
    # is guessed here.
    # ==================================================================

    # ==================================================================
    # Residual error: not reported numerically in Patel 2020. Encoded
    # as fixed(0) per the standing operator policy on unreported
    # residual error (the model returns the deterministic typical-
    # value prediction). Same encoding as
    # Darpo_2014_racSotalol_QTcI.R. The vignette Errata documents the
    # gap.
    # ==================================================================
    addSd <- fixed(0)
    label("Additive residual error standard deviation on DeltaQTcI (msec; ZERO - not reported in source)")
  })

  model({
    # ==================================================================
    # Linear concentration-DeltaQTcI prediction (Patel 2020 Figure 4A).
    # CP_SAPANISERTIB_NGML is supplied as a time-varying covariate in
    # ng/mL; the slope is already on the msec-per-ng/mL scale, so no
    # unit rescaling is applied. Output is the change from time-matched
    # day -1 baseline in QTcI, in msec.
    #
    # There are no ODE states: this is a purely algebraic PD model
    # driven by an external concentration column, so the file declares
    # no compartmentData and rxSolve() must NOT be given omega = NA
    # (there are no etas to suppress).
    # ==================================================================
    QTcI <- e0 + slope * CP_SAPANISERTIB_NGML

    QTcI ~ add(addSd)
  })
}
