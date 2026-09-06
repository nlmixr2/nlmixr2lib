Patel_2020_sapanisertib_RR <- function() {
  description <- paste(
    "Concentration-DeltaRR linear mixed-effects PD model for the mTORC1/2",
    "inhibitor sapanisertib (MLN0128 / INK128 / TAK-228) in 44 patients with",
    "advanced solid tumors, from a dedicated QTc (DQT) study of a single",
    "40 mg oral dose. The endpoint is change from time-matched day -1",
    "baseline in the RR interval (the time between two consecutive R waves,",
    "i.e. the reciprocal-scale heart-rate measure; observation variable dRR,",
    "msec). The model is:",
    "dRR = e0 + slope * CP_SAPANISERTIB_NGML,",
    "with e0 = -25.504 msec and slope = 0.147 msec per ng/mL. Unlike the two",
    "companion QTc models the slope here is POSITIVE and strongly",
    "significant (P < 0.001): higher sapanisertib concentration is",
    "associated with RR LENGTHENING, i.e. heart-rate SLOWING. The large",
    "negative intercept means that at zero concentration RR is ~25 msec",
    "below its time-matched baseline, so the net predicted DeltaRR crosses",
    "zero near 174 ng/mL and reaches about +18 msec at the 40 mg geometric",
    "mean Cmax of 297 ng/mL.",
    "PD-only model: sapanisertib plasma concentration is supplied as a",
    "time-varying covariate CP_SAPANISERTIB_NGML (ng/mL). The source",
    "publication does not fit a population PK model -- the PK was",
    "characterised by NCA only (geometric mean Cmax 297 ng/mL, AUCinf",
    "2480 ng*h/mL, CL/F 16.1 L/h, t1/2 9.5 h after 40 mg); users wishing to",
    "drive the PD model from a simulated PK source must supply their own",
    "concentration trajectory (no sapanisertib popPK model exists in the",
    "nlmixr2lib registry). Companion model files",
    "Patel_2020_sapanisertib_QTcI.R and Patel_2020_sapanisertib_QTcF.R",
    "report the same linear structure for the individually rate-corrected",
    "and Fridericia-corrected QT intervals.",
    sep = " "
  )

  reference <- paste(
    "Patel C, Goel S, Patel MR, Rangachari L, Wilbur JD, Shou Y,",
    "Venkatakrishnan K, Lockhart AC. (2020). Phase 1 Study to Evaluate the",
    "Effect of the Investigational Anticancer Agent Sapanisertib on the QTc",
    "Interval in Patients With Advanced Solid Tumors.",
    "Clinical Pharmacology in Drug Development 9(7):876-888.",
    "doi:10.1002/cpdd.808.",
    "Parameter values are from Supplementary Figure S2 of the supplemental",
    "information (file cpdd808-sup-0001-suppmat.docx).",
    sep = " "
  )

  vignette <- "Patel_2020_sapanisertib"

  units <- list(
    time          = "h",
    dosing        = "(none; PD-only model fed by an external sapanisertib plasma-concentration covariate)",
    concentration = "(observation dRR is the change from time-matched day -1 baseline in the RR interval, msec; driving covariate CP_SAPANISERTIB_NGML is in ng/mL)"
  )

  covariateData <- list(
    CP_SAPANISERTIB_NGML = list(
      description        = "Instantaneous sapanisertib plasma concentration at the time of each PD observation, supplied as a time-varying covariate from observed plasma samples or an upstream PK source.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying per event row. Drives the linear concentration-",
        "DeltaRR expression DeltaRR = e0 + slope *",
        "CP_SAPANISERTIB_NGML.",
        "The slope is reported directly in msec per ng/mL (Patel 2020",
        "Supplementary Figure S2 in-panel parameter table, x-axis",
        "'Concentration (ng ml-1)'), so CP_SAPANISERTIB_NGML is supplied",
        "in ng/mL with no in-model unit rescaling.",
        "In Patel 2020 this was the observed sapanisertib plasma",
        "concentration determined by a validated LC-MS/MS method with",
        "deuterated sapanisertib as internal standard (dynamic range",
        "1-1000 ng/mL; Patel 2020 Methods 'Electrocardiogram, PK, and",
        "Safety Assessments'), drawn time-matched to the triplicate ECGs",
        "at 0 (<=15 min predose), 0.25, 0.5, 1, 1.5, 2, 2.5, 3, 4, 6, 8,",
        "10, 24 and 48 h after the single 40 mg oral dose.",
        "Reference values observed: geometric mean Cmax 297 ng/mL (%CV",
        "52.9) after 40 mg, median 317 (range 94-976) ng/mL; median Tmax",
        "1.52 h (range 0.50-24.00) (Patel 2020 Table S2). The",
        "concentration-RR dataset spans roughly 0-976 ng/mL (Patel 2020",
        "Figure S2 x-axis extent). Geometric mean Cmax at the two lower",
        "clinical doses are 36.9 ng/mL (4 mg once daily) and 235 ng/mL",
        "(30 mg once weekly) (Patel 2020 Table 4 footnotes a and b);",
        "note that Table 4 itself tabulates projections for DeltaQTcI",
        "and DeltaQTcF only, NOT for DeltaRR, so this model has no",
        "published point-prediction table to gate against.",
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
      "day 3 window contributes to the concentration-RR analysis",
      "(Patel 2020 Methods 'Study Design and Patients')."
    ),
    regions          = "United States (5 centers: Montefiore Medical Center, Bronx NY; Washington University School of Medicine, St. Louis MO; Stephenson Cancer Center, University of Oklahoma, Oklahoma City OK; Mary Crowley Cancer Research Center, Dallas TX; Florida Cancer Specialists / Sarah Cannon Research Institute, Sarasota FL).",
    notes            = paste(
      "NCT02197572. Open-label, single-arm, dedicated QTc (DQT) study.",
      "Time-matched baseline design: serial triplicate 12-lead ECGs were",
      "recorded on day -1 by Holter H12+ recorder at the same clock",
      "times as the day 1 postdose schedule, and DeltaRR was formed by",
      "subtracting the time-matched day -1 mean RR from the day 1 to",
      "day 3 mean RR.",
      "RR is the reciprocal-scale partner of heart rate; the paper",
      "modelled RR rather than heart rate for the concentration-effect",
      "analysis (Patel 2020 Results 'Relationship of Sapanisertib",
      "Concentration and Heart Rate' and Figure S2) while reporting the",
      "time-course analysis on the heart-rate scale (Table S1).",
      "The TIME course and the CONCENTRATION course of heart rate point",
      "in opposite directions and are not in conflict: heart rate rose",
      "progressively with time, reaching its maximum LSM increase of",
      "+12.4 bpm at 24 h after dosing (Patel 2020 Table S1) -- by which",
      "point sapanisertib concentrations have fallen well below Cmax",
      "given the 9.5 h half-life -- whereas the concentration-effect fit",
      "encoded here finds RR LENGTHENING (heart-rate slowing) with",
      "increasing concentration. The negative intercept (-25.5 msec, a",
      "faster heart rate than baseline at zero concentration) is what",
      "carries the observed net tachycardia. Mean heart-rate values",
      "remained inside the normal range throughout (64.5-88.4 bpm on",
      "treatment versus 68.1-74.0 bpm at baseline).",
      "ANALYSIS POPULATIONS DIFFER: all 44 patients contributed to the",
      "concentration-effect analysis that this model file encodes,",
      "whereas the separate repeated-measures mixed-effects LSM analysis",
      "(Patel 2020 Table S1) used only the 32 patients with a complete",
      "day -1 / day 1 Holter series. Do not gate this model against",
      "Table S1.",
      "Model-fitting software was NONMEM 7.3 (ICON Development",
      "Solutions) with post-processing in R; 95% CIs for the population",
      "parameters were constructed by nonparametric bootstrap (Patel",
      "2020 Methods 'Statistical Analyses'). Residual degrees of freedom",
      "for this DeltaRR fit were 470 (Patel 2020 Figure S2 in-panel",
      "table) -- the largest of the three fits, as RR needs neither a",
      "per-patient QT/RR regression nor a QT measurement.",
      "The paper reports that the concentration-effect models are mixed-",
      "effects models (patient as a random effect) but does NOT report",
      "any variance component -- see the ini() block comment and the",
      "vignette Errata."
    )
  )

  ini({
    # ==================================================================
    # Linear concentration-DeltaRR model. Both parameters are printed
    # inside the Supplementary Figure S2 plot panel as a regression
    # summary table (columns: Parameter | Estimate | SE | df | t-value
    # | P-value | 95% Bootstrap CI). That in-panel table is a PRINTED
    # parameter table with full printed-value authority -- it is not a
    # digitised or curve-fitted value. It is invisible to every
    # text-extraction path (the trimmed markdown of the .docx renders
    # the figure as `<!-- image -->`); it was read by unzipping
    # cpdd808-sup-0001-suppmat.docx and rendering word/media/image4.png
    # at 2x.
    #
    #   Parameter | Estimate |    SE | df  | t-value | P-value | 95% boot CI
    #   Intercept |  -25.504 | 9.453 | 470 |  -2.698 |   0.007 | (-44.154, -7.851)
    #   Slope     |    0.147 | 0.034 | 470 |   4.278 |   0.000 | (0.082, 0.242)
    #
    # Source equation (Patel 2020 Results 'Relationship of Sapanisertib
    # Concentration and Heart Rate', Figure S2):
    #   DeltaRR = -25.504 + 0.147 * C_sapanisertib(ng/mL)
    #
    # Transcription was verified by triple redundancy: the intercept
    # t-value reproduces exactly from estimate / SE
    # (-25.504 / 9.453 = -2.698), every P-value reproduces from its
    # t-value and df via 2*pt(-|t|, df) (0.007 -> 0.0072;
    # 0.000 -> 2.3e-5), and both bootstrap CIs contain their point
    # estimates.
    #
    # Unlike the two QTc models, this fit has NO published
    # point-prediction table (Table 4 covers DeltaQTcI and DeltaQTcF
    # only), so there is no second printed source to cross-check the
    # estimates against. The internal Estimate / SE / t / P / CI
    # redundancy above is the whole of the available verification.
    # ==================================================================

    e0 <- -25.504
    label("Linear-mixed-effects intercept e0 on DeltaRR (msec)")
    # Patel 2020 Figure S2 in-panel table, 'Intercept' row: estimate
    # -25.504 (SE 9.453, df 470, t -2.698, P = 0.007, 95% bootstrap CI
    # -44.154 to -7.851). Significantly BELOW zero and by far the
    # largest-magnitude intercept of the three fits: at zero
    # concentration the RR interval sits ~25 msec under its
    # time-matched baseline, i.e. the heart rate is faster. This is the
    # term that carries the net tachycardia the paper reports on the
    # heart-rate scale (Table S1). Deliberately NOT log-transformed:
    # the intercept is negative and must be free to take either sign.
    # Canonical bare PD baseline parameter `e0`.

    slope <- 0.147
    label("Linear concentration-DeltaRR slope (msec per ng/mL)")
    # Patel 2020 Figure S2 in-panel table, 'Slope' row: estimate 0.147
    # (SE 0.034, df 470, t 4.278, P = 0.000, 95% bootstrap CI 0.082 to
    # 0.242). POSITIVE -- the opposite sign to both QTc slopes -- and
    # the most strongly significant of the six parameters across the
    # three fits (P = 2.3e-5 recomputed). Deliberately NOT
    # log-transformed, for consistency with the companion QTc files
    # where a negative slope makes a log parameterisation impossible;
    # `slope` is kept on the natural scale throughout the three-model
    # set so the three files are directly comparable. Interpretation:
    # RR lengthening / heart-rate slowing with rising concentration,
    # partially offsetting the negative intercept. The net predicted
    # DeltaRR crosses zero at 25.504 / 0.147 = 173.5 ng/mL.

    # ==================================================================
    # Inter-individual variability: Patel 2020 Methods 'Statistical
    # Analyses' states that "Mixed-effects models were subsequently
    # developed to describe the direct effect of sapanisertib on change
    # from time-matched baseline in RR, QTcI, and QTcF", and the
    # residual df of 470 against 44 patients is consistent with a
    # per-subject random intercept. However the paper does NOT report
    # any variance estimate for any of the three concentration-effect
    # models. Per the standing operator policy on unreported IIV, this
    # file omits the eta declarations and ships a typical-value-only
    # model. The vignette Errata documents the gap.
    # ==================================================================

    # ==================================================================
    # Residual error: not reported numerically in Patel 2020. Encoded
    # as fixed(0) per the standing operator policy on unreported
    # residual error (the model returns the deterministic typical-
    # value prediction). The vignette Errata documents the gap.
    # ==================================================================
    addSd <- fixed(0)
    label("Additive residual error standard deviation on DeltaRR (msec; ZERO - not reported in source)")
  })

  model({
    # ==================================================================
    # Linear concentration-DeltaRR prediction (Patel 2020 Figure S2).
    # CP_SAPANISERTIB_NGML is supplied as a time-varying covariate in
    # ng/mL; the slope is already on the msec-per-ng/mL scale, so no
    # unit rescaling is applied. Output is the change from time-matched
    # day -1 baseline in the RR interval, in msec.
    #
    # The observation variable is named `dRR`, the canonical
    # change-from-baseline RR interval registered in
    # inst/references/compartment-names.md alongside this extraction
    # (founding example). The register had no RR entry of either form,
    # and the two nearest families gave opposite guidance: the bare
    # interval names `QTcI` / `QTcF` / `QTcP` are used directly as the
    # observation variable in change-from-baseline concentration-QTc
    # models (Darpo_2014_racSotalol_*.R,
    # Fostvedt_2021_glasdegib_QTcF.R, Mukker_2026_tuvusertib_QTcF.R,
    # Zhou_2025_fruquintinib_*.R), whereas `dHR` exists precisely for a
    # change-from-baseline heart-rate-domain endpoint. Operator ruling
    # (sidecar oare_PMC7586797 request-001 q1, answered 2026-09-02)
    # adopted `dRR`: the bare-name precedent is for ABSOLUTE interval
    # names, and a delta is what `dHR` already encodes, so the
    # heart-rate-domain pair (dHR, dRR) is kept mutually consistent.
    # The scale-incomparability rationale transfers directly -- an
    # absolute RR interval is ~850 msec, this DeltaRR is -25 to
    # +18 msec. The consequence is a deliberate naming asymmetry
    # WITHIN this three-file set: the two QTc siblings keep their bare
    # corrected-interval names (`QTcI`, `QTcF`) per the standing QTc
    # precedent, while the heart-rate endpoint takes the `d` prefix.
    #
    # There are no ODE states: this is a purely algebraic PD model
    # driven by an external concentration column, so the file declares
    # no compartmentData and rxSolve() must NOT be given omega = NA
    # (there are no etas to suppress).
    # ==================================================================
    dRR <- e0 + slope * CP_SAPANISERTIB_NGML

    dRR ~ add(addSd)
  })
}
