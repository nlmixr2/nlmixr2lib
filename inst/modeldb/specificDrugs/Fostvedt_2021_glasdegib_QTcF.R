Fostvedt_2021_glasdegib_QTcF <- function() {
  description <- paste(
    "Population PD model for glasdegib concentration-driven prolongation",
    "of the QT interval corrected for heart rate using Fridericia's",
    "formula (QTcF) in 70 adult patients with advanced cancer pooled",
    "from two phase 1 dose-escalation trials (B1371001 in hematologic",
    "malignancies; B1371002 in solid tumors). The exposure-response",
    "form is a linear mixed-effects model:",
    "QTcF = theta1 + theta2 * (CP_GLASDEGIB_NGML / 1000) + eta1 + W * eps,",
    "with an additive random effect on the intercept (eta1 ~ N(0,",
    "omega2_1)) and a 'thetarized' additive residual error (W is a",
    "fitted scalar; eps ~ N(0, 1)). The covariate analysis (age, sex,",
    "study) retained no covariates; a random effect on the slope was",
    "considered and dropped because shrinkage exceeded 20%. PD-only",
    "model: plasma glasdegib concentration is supplied as a time-varying",
    "covariate CP_GLASDEGIB_NGML (ng/mL). The slope is reported in the",
    "source publication on the microgram-per-mL scale (4.3 msec per",
    "microgram per mL) which is equivalent to 0.0043 msec per ng/mL;",
    "this file keeps the slope in the paper's microgram-per-mL form and",
    "applies the unit conversion `CP_GLASDEGIB_NGML / 1000` inside",
    "model(). The source publication does not fit a population PK model;",
    "users wishing to drive the PD model from a simulated PK source must",
    "supply their own concentration trajectory (no glasdegib popPK model",
    "exists in the nlmixr2lib registry).",
    sep = " "
  )

  reference <- paste(
    "Fostvedt LK, Shaik N, Martinelli G, Wagner AJ, Ruiz-Garcia A.",
    "Exposure-response modeling of the effect of glasdegib on cardiac",
    "repolarization in patients with cancer.",
    "Expert Review of Clinical Pharmacology 2021;14(7):927-935.",
    "doi:10.1080/17512433.2021.1925538.",
    sep = " "
  )

  vignette <- "Fostvedt_2021_glasdegib"

  # Paper-specific etas: the source paper places the random intercept
  # additively on the linear QTcF scale (eta1 ~ N(0, omega2_1) added
  # directly to the typical-value intercept theta1). The encoding here
  # log-transforms theta1 to le0 for canonical positivity, so the linear-
  # scale variance 277.6 msec^2 is translated to a log-normal equivalent
  # omega^2_log = log(1 + (sqrt(277.6) / 423.8)^2) = 0.001544 (CV ~3.93%
  # on the baseline); the typical-value predictions are identical and the
  # IIV envelope is approximated to within ~0.8% for this small CV. See
  # the validation vignette's Assumptions and deviations section for the
  # log-vs-linear-eta accounting. The eta name `etale0` pairs canonically
  # with the `le0` parameter (no paper_specific_etas declaration needed).

  units <- list(
    time          = "h",
    dosing        = "(none; PD-only model fed by an external glasdegib plasma-concentration covariate)",
    concentration = "(observation QTcF is the QT interval corrected for heart rate using Fridericia's formula, msec; driving covariate CP_GLASDEGIB_NGML is in ng/mL)"
  )

  covariateData <- list(
    CP_GLASDEGIB_NGML = list(
      description        = "Instantaneous glasdegib plasma concentration at the time of each PD observation, supplied as a time-varying covariate from observed plasma samples or an upstream PK source.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying per event row. Drives the linear concentration-QTc expression QTcF = e0 + slope * (CP_GLASDEGIB_NGML / 1000).",
        "In Fostvedt 2021 this was the observed glasdegib plasma concentration around the time of the ECG collection (paired within +/-15 minutes of the PK sample; Methods 2.2 and 2.3). The source NM-TRAN column name in the model is CONC.",
        "The slope is reported on the microgram-per-mL scale in the paper (Table 3 footnote: 'The estimates are reported on the microgram scale as the scaling helped the estimation procedure'); the canonical covariate is in ng/mL so the `/ 1000` rescaling lives inside model().",
        "Reference values observed: geometric-mean steady-state Cmax was 1137 ng/mL at 100 mg QD (therapeutic dose) and 2445 ng/mL at 200 mg QD (Fostvedt 2021 Methods 2.3 and Table 4; the Cmax values come from a separate Pfizer phase 2 study NCT01546038 cited as reference [15] in the source paper).",
        "Set to 0 outside the drug-exposure window (the concentration-slope term then collapses to 0)."
      ),
      source_name        = "CONC"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age (years)",
      units       = "years",
      type        = "continuous",
      notes       = "Screened by forward-selection (alpha = 0.05) / backward-elimination (alpha = 0.001) stepwise covariate procedure on the QTcF and QTcS model parameters and not retained (Fostvedt 2021 Methods 2.3 and Results 3.3)."
    ),
    SEXF = list(
      description = "Biological sex (1 = female)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened (paper covariate 'sex') and not retained (Fostvedt 2021 Results 3.3). Of 70 patients in the pooled analysis, 28 (40%) were female (Fostvedt 2021 Table 1)."
    ),
    STUDY = list(
      description = "Source phase 1 study (hematologic malignancies B1371001 vs solid tumors B1371002)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened (paper covariate 'study (hematologic versus solid tumors)') and not retained (Fostvedt 2021 Results 3.3 and Discussion: 'no statistically significant differences were found based on the pre-specified backwards elimination significant criteria of 0.001'). Electrolyte imbalances and CYP3A4 inhibitor / inducer comedications were prohibited by protocol and therefore could not be tested as covariates (Fostvedt 2021 Methods 2.3)."
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 70L,
    n_studies        = 2L,
    n_observations   = 747L,
    age_range        = "25-89 years (Study B1371001 25-89, median 69; Study B1371002 27-76, median 61; pooled median 67)",
    weight_range     = NA_character_,
    sex_female_pct   = 40,
    race_ethnicity   = c(White = 80, Asian = 6, Black = 6, Other = 9),
    disease_state    = "Patients with advanced cancer pooled from two phase 1 dose-escalation trials: hematologic malignancies (Study B1371001, NCT00953758, n = 47) and solid tumor malignancies (Study B1371002, NCT01286467, n = 23). Eligibility excluded patients with screening QTcF > 470 msec; on-study QTcF > 480 msec triggered correction of reversible causes (electrolyte abnormalities, hypoxia, concomitant medications with QTc-prolonging potential).",
    dose_range       = paste(
      "Glasdegib QD oral monotherapy, with PK-ECG pairs collected across a wide dose range (Fostvedt 2021 Section 2.1 and Discussion).",
      "Study B1371001: 5, 10, 20, 40, 80, 120, 180, 270, 400, and 600 mg QD; Cycle 1 was preceded by a single lead-in dose for single-dose PK characterisation, with continuous once-daily dosing in 28-day cycles thereafter. PK-ECG pairs collected pre dose, 1, 4, and 24 h post lead-in dose; Cycle 1 Day 1 pre dose and 1 h post dose; Cycle 1 Day 8 1 h post dose; Cycle 1 Day 15 1 h post dose; Cycle 1 Day 21 pre dose, 1, 2, 4, and 24 h post dose; and Day 1 of every subsequent cycle at 1 h post dose.",
      "Study B1371002: 80, 160, 320, and 640 mg QD. Cycle 1 administered glasdegib for 25 days followed by 3 days off treatment to characterise steady-state elimination; Cycle 2 and beyond continuous QD. Triplicate ECGs scheduled Days 1, 15, and 25 of Cycle 1 and Day 1 of each subsequent cycle, paired with PK pre dose and 2 h post dose; Cycle 1 Day 25 additional PK-ECG pairs at pre dose and 2, 6, and 24 h post dose.",
      "Maximum dose tested (640 mg QD) was more than 6-fold above the selected clinical dose of 100 mg QD."
    ),
    regions          = NA_character_,
    notes            = paste(
      "Concentration-QTc analysis data set: 747 PK-ECG pairs (B1371001 n = 589; B1371002 n = 158) collected within +/-15 minutes of paired PK draws, from 70 patients (Fostvedt 2021 Methods 2.2 and Table 1).",
      "Both studies were open-label, multicentre, phase 1 first-in-patient (B1371001) or solid-tumor phase 1 (B1371002) dose-escalation designs (Fostvedt 2021 Section 2.1).",
      "ECG: triplicate 12-lead measurements with a 10-second rhythm strip, collected approximately 2 to 5 minutes apart at each scheduled timepoint. Baseline QTcF mean (SD) = 421.1 (18.1) msec; on-treatment QTcF mean (SD) = 430.6 (24.0) msec (Fostvedt 2021 Table 2).",
      "Glasdegib has a plasma terminal half-life of approximately 17 hours and is metabolised by CYP3A4 (Fostvedt 2021 Sections 2.4 and 4).",
      "A separate parametric bootstrap was performed by the authors at therapeutic and supratherapeutic Cmax to estimate the 95% CI of the predicted mean QTcF change from baseline (Fostvedt 2021 Section 3.4 and Table 4). A subsequent ICH E14-compliant thorough QT (TQT) study confirmed the conclusions: estimated slope 0.005 msec/(ng/mL) (TQT) vs. 0.00430 msec/(ng/mL) (this E-R analysis) for QTcF; estimated mean prolongation at 2445 ng/mL was 12.09 msec (95% CI 10.03-14.25) in this E-R analysis vs. the TQT 90% CI 11.95-15.49 msec (Fostvedt 2021 Discussion).",
      "The exposure-response analysis was performed using NONMEM (Fostvedt 2021 implicit from the NM-TRAN-style equation and stepwise-covariate procedure)."
    )
  )

  ini({
    # ==================================================================
    # Linear concentration-QTcF model (Fostvedt 2021 Section 2.3 and
    # Table 3 QTcF row). The published equation in the source paper is:
    #   QTcF_ij = theta1 + theta2 * CONC_ij + eta1_i + W * eps_ij
    # with theta2 reported on the microgram-per-mL scale (the slope was
    # estimated with a CONC / 1000 rescaling internal to NONMEM).
    # ==================================================================

    le0 <- log(423.8)
    label("Baseline QTcF intercept theta1 (msec)")
    # Fostvedt 2021 Table 3 QTcF 'Intercept, theta_1' = 423.8 (SE 2.1; RSE 0.5%).
    # Canonical log-transformed PD baseline; the linear-scale typical
    # value of 423.8 msec is recovered as exp(le0) inside model().

    lslope <- log(4.3)
    label("Concentration-QTcF slope theta2 (msec per microgram per mL)")
    # Fostvedt 2021 Table 3 QTcF 'Slope, theta_2' = 4.3 (SE 0.37; RSE 8.6%).
    # Reported on the microgram-per-mL scale per the Table 3 / Figure 2
    # footnote ("The estimates are reported on the microgram scale as the
    # scaling helped the estimation procedure"). Equivalent to 0.0043
    # msec per ng/mL (Fostvedt 2021 Discussion: "the estimated slopes
    # were 0.00430 ... msec/ng/mL"). Paper-specific log-transformed
    # linear concentration-effect slope; the bare value 4.3 is recovered
    # as exp(lslope) inside model(). Established precedent in the
    # registry: Liesenfeld_2006_dabigatran_aPTT.R, deVriesSchultink_2018
    # _anthracycline_troponinT.R.

    # ==================================================================
    # Random effect on the intercept (Fostvedt 2021 Section 2.3:
    # "eta_1 ~ N(0, omega^2_1)" added additively to the linear-scale
    # QTcF intercept; Table 3 QTcF 'Variance of the random effect on
    # the intercept, omega^2_1' = 277.6 msec^2). No random effect on the
    # slope: Section 2.3 reports it was dropped because "the shrinkage
    # was very large (>20%)".
    #
    # Translation from linear-scale variance to log-normal-equivalent
    # omega^2_log:
    #   omega^2_log = log(1 + (sqrt(277.6) / 423.8)^2)
    #               = log(1 + 277.6 / 423.8^2)
    #               = log(1.001546) ~= 0.0015446
    # CV at the linear-scale typical value: sqrt(277.6) / 423.8 ~= 3.93%.
    # The log-vs-linear approximation is within ~0.8% at this CV; the
    # vignette Errata documents the conversion.
    # ==================================================================
    etale0 ~ 0.0015446
    # Fostvedt 2021 Table 3 QTcF omega^2_1 = 277.6 (SE 50.7; RSE 18.3%;
    # shrinkage 4.4%); converted to log-normal equivalent for canonical
    # log-transformed intercept encoding.

    # ==================================================================
    # Residual error (Fostvedt 2021 Section 2.3: W is a fitted scalar
    # and eps ~ N(0, 1), so the residual error is additive on the QTcF
    # scale with SD = W; Table 3 QTcF 'Residual standard deviation, W'
    # = 14.1 msec). This maps directly onto the nlmixr2 canonical
    # additive residual `addSd`.
    # ==================================================================
    addSd <- 14.1
    label("Additive residual error standard deviation on QTcF (msec)")
    # Fostvedt 2021 Table 3 QTcF W = 14.1 (SE 1.1; RSE 8.0%; shrinkage 4.0%).
  })

  model({
    # ==================================================================
    # Linear concentration-QTcF model. CP_GLASDEGIB_NGML is supplied as
    # a time-varying covariate (ng/mL); the slope is on the microgram-
    # per-mL scale (per Fostvedt 2021 Table 3 footnote), so the in-model
    # `/ 1000` performs the ng/mL -> microgram/mL conversion.
    # ==================================================================
    e0    <- exp(le0 + etale0)
    slope <- exp(lslope)

    QTcF <- e0 + slope * (CP_GLASDEGIB_NGML / 1000)
    QTcF ~ add(addSd)
  })
}
