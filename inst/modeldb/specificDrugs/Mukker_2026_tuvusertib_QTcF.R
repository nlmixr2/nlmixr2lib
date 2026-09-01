Mukker_2026_tuvusertib_QTcF <- function() {
  description <- paste(
    "Garnett-type linear mixed-effects concentration-QTc model for the",
    "investigational ATR kinase inhibitor tuvusertib (M1774) in 55",
    "patients with advanced solid tumors enrolled in Part A1 of the",
    "phase I first-in-human study DDRiver Solid Tumors 301",
    "(NCT04170153; dose range 5-270 mg). The endpoint is the change",
    "from baseline in the Fridericia-corrected QT interval",
    "(DeltaQTcF, ms). The published model is:",
    "DeltaQTcF = (theta0 + eta0) + (theta1 + eta1) * C",
    "+ theta2(nominal time) + theta3 * (QTC_BL - 422) + eps,",
    "with theta0 = -0.125 ms, theta1 = 0.00244 ms/(ng/mL), four nominal",
    "post-dose-time intercept shifts (-1.59, 1.11, 2.87, -0.402 ms at 0,",
    "1, 2 and 3 h) attributed by the authors to diurnal variation, and a",
    "centered-baseline-QTcF coefficient theta3 = -0.0980 ms/ms. Random",
    "effects sit additively on both the intercept (SD 5.41 ms) and the",
    "concentration slope (SD 0.0026 ms/(ng/mL)) with an estimated",
    "correlation of -0.00273; residual SD is 7.07 ms. PD-only model:",
    "tuvusertib plasma concentration is supplied as a time-varying",
    "covariate CP_TUVUSERTIB_NGML (ng/mL). No covariate (body weight,",
    "age, race, sex) reached significance. The source publication does",
    "not fit a population PK model -- concentrations enter as observed,",
    "time-matched measurements -- so users wishing to drive the PD model",
    "from a simulated PK source must supply their own concentration",
    "trajectory. Companion model files Mukker_2026_tuvusertib_HR.R (the",
    "concentration-DeltaHR model that establishes the first Garnett",
    "assumption) and Mukker_2026_tuvusertib_hERG.R (the nonclinical in",
    "vitro hERG concentration-response) complete the paper's integrated",
    "nonclinical-and-clinical QTc risk assessment.",
    sep = " "
  )

  reference <- paste(
    "Mukker JK, Yap TA, Tolcher AW, de Bono JS, Plummer R, Grosser G,",
    "van Amsterdam C, Schieferstein H, Witjes H, Diderichsen PM,",
    "Krebs-Brown A, Gao W, Strotmann R, Szucs Z, Gounaris I,",
    "Venkatakrishnan K. An Integrated Nonclinical and Clinical Risk",
    "Assessment of the Effects of Investigational ATRi Tuvusertib on QTc",
    "Interval in Patients With Solid Tumors.",
    "Clinical and Translational Science 2026;19(2):e70496.",
    "doi:10.1111/cts.70496.",
    sep = " "
  )

  vignette <- "Mukker_2026_tuvusertib_QTc"

  units <- list(
    time          = "h",
    dosing        = "(none; PD-only model fed by an external tuvusertib plasma-concentration covariate)",
    concentration = "(observation QTcF is the change from baseline in the Fridericia-corrected QT interval, ms; driving covariate CP_TUVUSERTIB_NGML is in ng/mL)"
  )

  covariateData <- list(
    CP_TUVUSERTIB_NGML = list(
      description        = "Instantaneous tuvusertib plasma concentration at the time of each ECG observation, supplied as a time-varying covariate from observed plasma samples or an upstream PK source.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying per event row. Drives the linear concentration-DeltaQTcF term (theta1 + eta1) * CP_TUVUSERTIB_NGML.",
        "In Mukker 2026 this was the observed tuvusertib plasma concentration in the time-matched PK-ECG dataset: serial blood samples were drawn immediately after each triplicate 12-lead ECG, and concentrations were measured by a validated LC/MS method with a lower limit of quantification of 0.5 ng/mL (Mukker 2026 Methods 2.2.2).",
        "The slope is reported directly in ms per ng/mL (Mukker 2026 Table 2), so CP_TUVUSERTIB_NGML is supplied in ng/mL with no in-model unit rescaling.",
        "Reference values observed: geometric mean steady-state Cmax 1410 ng/mL at 180 mg QD (the recommended dose for expansion), with 2x and 3x that value (2820 and 4230 ng/mL) used for the supratherapeutic projections in Table 3; median 524, P90 1732, P95 2252 and maximum 3290 ng/mL across the C-QTc dataset (Table S2).",
        "Set to 0 outside the drug-exposure window (the concentration-slope term then collapses to 0, leaving the intercept, nominal-time and baseline terms)."
      ),
      source_name        = "tuvusertib plasma concentration (C in Equation 1)"
    ),
    QTC_BL = list(
      description        = "Subject's baseline (time-zero) Fridericia-corrected QT interval, treated as a per-subject time-fixed covariate. Enters the linear-mixed-effects intercept as the centered term e_qtc_bl_e0 * (QTC_BL - 422). Set QTC_BL = 422 ms for the typical subject -- the centered term then collapses to 0.",
      units              = "ms",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Derived from the triplicate 12-lead ECGs recorded at baseline; triplicates were averaged to a single mean value per patient per timepoint and QTcF computed per Fridericia (Mukker 2026 Methods 2.2.2).",
        "Centering reference: unlike Darpo 2014 (where the cohort median was not published and a rounded 390 ms standard had to be assumed), Mukker 2026 PUBLISHES the centering constant. The Equation 1 legend defines QTcF0 as 'the overall mean QTcF_i0, that is, the mean of all the baseline (time zero) QTcF values', and Table 1 reports mean (SD) baseline QTcF = 422 (21.0) ms. The model therefore uses qtc_bl_ref = 422 ms with no assumption.",
        "Eligibility for the parent trial excluded patients with pre-existing QTc prolongation (average QTcF > 450 ms for males, > 470 ms for females; Mukker 2026 Methods 2.2.1), so the covariate range is truncated above.",
        "The QT correction method is not carried in the canonical name (per the QTC_BL register entry); here it is Fridericia."
      ),
      source_name        = "QTcF_i0 / Baseline QTcF"
    ),
    T_LASTDOSE = list(
      description        = "Nominal (protocol-scheduled) time after the most recent tuvusertib dose at which the triplicate ECG was recorded, in hours. Selects one of the four estimated nominal-timepoint intercept shifts.",
      units              = "h",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying per event row. Mukker 2026 Equation 1 carries a term theta2 * Time_j where theta2 is 'a factor estimating time-specific intercepts at nominal times j accounting for diurnal variation'. Table 2 reports four estimates, one per nominal post-dose hour (0, 1, 2, 3 h), each footnoted as 'the estimated difference from the population mean intercept' -- i.e. mean-centered class effects rather than contrasts against an omitted reference level.",
        "Because this is a PD-only model with no dosing records, rxode2's native tad() is unavailable and the post-dose clock must be supplied as a covariate column; see the T_LASTDOSE register entry.",
        "The column carries the NOMINAL scheduled hour, not the actual elapsed time: the class-effect model is defined on the nominal sampling grid. model() bins the supplied value to its nearest scheduled level (< 0.5 h -> 0 h; 0.5-1.5 h -> 1 h; 1.5-2.5 h -> 2 h; >= 2.5 h -> 3 h) so that small protocol-window slippage does not silently fall through to the wrong class.",
        "Every ECG in the analysis dataset falls at one of the four nominal times, so the four indicators are exhaustive and mutually exclusive by construction."
      ),
      source_name        = "nominal time after dose (Time_j in Equation 1)"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened as a candidate covariate on the C-DeltaQTcF model and not retained (p > 0.1; Mukker 2026 Results 3.2.2). Cohort mean (SD) 78.6 (18.0) kg (Table 1). No point estimate is reported anywhere on disk, so the effect cannot be encoded."
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened and not retained (p > 0.1; Mukker 2026 Results 3.2.2). Cohort mean (SD) 61.9 (10.9) years (Table 1)."
    ),
    SEXF = list(
      description = "Biological sex (1 = female)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened (paper covariate 'sex') and not retained (p > 0.1; Mukker 2026 Results 3.2.2). Cohort was 32 female (58.2%) / 23 male (41.8%) (Table 1). Note this is the covariate whose retention drives the headline finding of the Darpo 2014 rac-sotalol companion model in this registry; it did not reach significance here."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as part of the 'race' covariate and not retained (p > 0.1; Mukker 2026 Results 3.2.2). 5 of 55 patients (9.1%) were Asian (Table 1)."
    ),
    RACE_BLACK = list(
      description = "Black or African American race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as part of the 'race' covariate and not retained (p > 0.1; Mukker 2026 Results 3.2.2). 2 of 55 patients (3.6%) were Black or African American (Table 1)."
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 55L,
    n_studies        = 1L,
    age_range        = "mean (SD) 61.9 (10.9) years",
    weight_range     = "mean (SD) 78.6 (18.0) kg",
    sex_female_pct   = 58.2,
    race_ethnicity   = c(White = 76.4, Asian = 9.1, Black_or_African_American = 3.6, Other = 10.9),
    disease_state    = paste(
      "Patients with advanced solid tumors enrolled in Part A1 of DDRiver",
      "Solid Tumors 301 (NCT04170153), an open-label, first-in-human,",
      "multicenter phase I dose-escalation study of tuvusertib monotherapy.",
      "Patients with pre-existing QTc prolongation (average QTcF > 450 ms",
      "for males, > 470 ms for females) were excluded. Mean (SD) baseline",
      "QTcF interval 422 (21.0) ms (Mukker 2026 Table 1)."
    ),
    dose_range       = paste(
      "Tuvusertib 5-270 mg orally. Regimens contributing to the C-QTc",
      "analysis: 5, 10, 20, 40, 80, 130, 180, 220 and 270 mg once daily",
      "(QD) continuous; 180 and 220 mg QD 2 weeks on / 1 week off; and",
      "150 mg twice daily (BID) 4 days on / 3 days off. The maximum",
      "tolerated dose was 180 mg QD continuous and the recommended dose",
      "for expansion 180 mg QD 2 weeks on / 1 week off (Mukker 2026",
      "Abstract, Introduction and Figure S2 legend)."
    ),
    regions          = NA_character_,
    notes            = paste(
      "Inclusion in the C-QTc dataset required 12-lead triplicate ECG readings at baseline plus at least one post-baseline reading, each with a time-matched tuvusertib plasma concentration (Mukker 2026 Methods 2.2.1).",
      "ECG handling: digital triplicate ECGs collected at prespecified timepoints and centrally read by an independent cardiologist (IQVIA). Duplicate and single ECGs were excluded; triplicates were averaged to one mean value per patient per timepoint, and QTcF derived by the Fridericia formula (Mukker 2026 Methods 2.2.2). Nominal post-dose ECG times were 0, 1, 2 and 3 h (Table 2).",
      "The four Garnett assumptions were each checked and met (Mukker 2026 Results 3.2.2 and Supporting Information Results S1): (i) no tuvusertib effect on heart rate -- established by the companion C-DeltaHR model in Mukker_2026_tuvusertib_HR.R; (ii) QTcF independent of heart rate after Fridericia correction (p = 0.11 for the residual QTcF-RR dependence); (iii) no hysteresis between concentration and DeltaQTcF; (iv) no nonlinearity in the C-DeltaQTcF relationship by LOESS.",
      "Key model-based projections used as validation targets: DeltaQTcF_drug (the drug effect corrected for intercept, baseline and time, i.e. slope * C) was 3.44 ms (90% CI 1.39, 5.49) at the 180 mg QD geometric mean steady-state Cmax of 1410 ng/mL, 6.88 (2.78, 10.98) at 2x and 10.32 (4.17, 16.47) at 3x that concentration -- all upper bounds below the 20 ms threshold of concern for oncology drugs (Mukker 2026 Table 3).",
      "Estimation was by linear mixed-effects modeling in R 4.3.1 (lme-class fitting); the C-QTc dataset was assembled in SAS 9.4 (Mukker 2026 Methods 2.2.3).",
      "Data set assembly and the analysis were performed by Certara USA on behalf of the sponsor (EMD Serono / the healthcare business of Merck KGaA)."
    )
  )

  ini({
    # ==================================================================
    # Garnett linear mixed-effects concentration-QTc model.
    # Mukker 2026 Equation 1 (reproduced from the display-equation
    # image CTS-19-e70496-e001.jpg supplied in the article's
    # supplementary-file bundle; the trimmed text renders it only as a
    # 'formula-not-decoded' placeholder):
    #
    #   DeltaQTcF_ij = (theta0 + eta0_i)
    #                + (theta1 + eta1_i) * C_ij
    #                + theta2 * Time_j
    #                + theta3 * (QTcF_i0 - QTcF0)
    #                + eps_ij
    #
    # All fixed-effect values below are from Mukker 2026 Table 2
    # ('Parameter estimates of the final C-DeltaQTcF model').
    # ==================================================================

    e0 <- -0.125
    label("Population mean intercept theta0 on DeltaQTcF (ms)")
    # Mukker 2026 Table 2 'Mean intercept, ms' = -0.125 (95% CI -1.98,
    # 1.73). NOT log-transformed: the intercept of a change-from-
    # baseline QTc model is signed and is in fact negative here.
    # Same encoding decision as Darpo_2014_racSotalol_QTcF.R.

    slope <- 0.00244
    label("Concentration-DeltaQTcF slope theta1 (ms per ng/mL)")
    # Mukker 2026 Table 2 'Tuvusertib plasma concentration slope,
    # ms/(ng/mL)' = 0.00244 (95% CI 0.000708, 0.00417). The CI excludes
    # zero, which the authors cite as evidence of assay sensitivity in
    # the absence of a positive control (Results 3.2.2, Discussion).
    #
    # Deliberately NOT log-transformed, unlike the `lslope` encoding in
    # Darpo_2014_racSotalol_QTcF.R and Fostvedt_2021_glasdegib_QTcF.R.
    # The reason is structural, not cosmetic: this model carries an
    # ADDITIVE random effect on the slope whose SD (0.0026) is LARGER
    # than the typical value itself (0.00244), so individual slopes are
    # legitimately negative for a large share of the population. A
    # log-transform with a log-normal eta would force every individual
    # slope positive and would materially misstate the model. The
    # published equation places eta1 additively on the linear scale and
    # this file reproduces that exactly.

    # ------------------------------------------------------------------
    # Nominal-post-dose-time intercept shifts (theta2 * Time_j).
    # Mukker 2026 Table 2 rows 'Nominal time after dose (0h/1h/2h/3h)',
    # each footnoted (footnote a): 'The fixed-effect estimates for
    # nominal time after the dose represent the estimated difference
    # from the population mean intercept.' They are therefore mean-
    # centered class effects, not contrasts against an omitted
    # reference level -- all four are estimated and all four are
    # carried here. The authors attribute the term to diurnal variation
    # (Methods 2.2.3).
    # ------------------------------------------------------------------

    e_tad0h_e0 <- -1.59
    label("Effect of the 0 h nominal post-dose timepoint on e0 (additive, ms)")
    # Mukker 2026 Table 2 'Nominal time after dose (0h), ms' = -1.59
    # (95% CI -3.34, 0.167).

    e_tad1h_e0 <- 1.11
    label("Effect of the 1 h nominal post-dose timepoint on e0 (additive, ms)")
    # Mukker 2026 Table 2 'Nominal time after dose (1h), ms' = 1.11
    # (95% CI -0.162, 2.38).

    e_tad2h_e0 <- 2.87
    label("Effect of the 2 h nominal post-dose timepoint on e0 (additive, ms)")
    # Mukker 2026 Table 2 'Nominal time after dose (2h), ms' = 2.87
    # (95% CI 1.51, 4.23). The only nominal-time effect whose CI
    # excludes zero.

    e_tad3h_e0 <- -0.402
    label("Effect of the 3 h nominal post-dose timepoint on e0 (additive, ms)")
    # Mukker 2026 Table 2 'Nominal time after dose (3h), ms' = -0.402
    # (95% CI -1.69, 0.884).

    e_qtc_bl_e0 <- -0.0980
    label("Effect of centered baseline QTcF on e0 (additive, ms per ms)")
    # Mukker 2026 Table 2 'Baseline QTcF, ms' = -0.0980 (95% CI -0.176,
    # -0.0196). Negative: patients whose baseline QTcF is above the
    # cohort mean show a smaller change from that baseline (regression
    # to the mean), the same sign and interpretation as the -0.70 ms/ms
    # coefficient in Darpo_2014_racSotalol_QTcF.R.

    # ==================================================================
    # Random effects. Mukker 2026 Methods 2.2.3: 'Patients were
    # included as an additive random effect on both intercept and slope
    # terms', with 'random effects ... normally distributed with mean
    # [0,0] and an unstructured covariance matrix Omega'. Both etas are
    # therefore additive on the linear (ms) scale.
    #
    # Table 2 reports the variance components as STANDARD DEVIATIONS
    # plus a correlation, not as variances:
    #   SD of intercept                              = 5.41   ms
    #   SD of tuvusertib plasma concentration effect = 0.0026 ms/(ng/mL)
    #   Correlation between the two                  = -0.00273
    # nlmixr2 wants a lower-triangular VARIANCE-covariance block, so:
    #   var(eta0)     = 5.41^2               = 29.2681
    #   cov(eta0,eta1)= -0.00273 * 5.41 * 0.0026
    #                                        = -3.840018e-05
    #   var(eta1)     = 0.0026^2             = 6.76e-06
    # Positive-definiteness check: cov^2 = 1.47e-09 is far below
    # var(eta0) * var(eta1) = 1.979e-04, so the block is valid. The
    # correlation is negligible in practice (-0.003); it is carried
    # anyway because the paper reports an unstructured Omega.
    # ==================================================================
    etae0 + etaslope ~ c(
      29.2681,
      -3.840018e-05, 6.76e-06
    )

    # ==================================================================
    # Residual unexplained variability. Mukker 2026 Methods 2.2.3:
    # 'the residual unexplained variability, eps_ij, was normally
    # distributed with mean 0 and variance sigma^2' -- i.e. additive on
    # the DeltaQTcF (ms) scale. Table 2 'Residual SD, ms' = 7.07.
    # ==================================================================
    addSd <- 7.07
    label("Additive residual error standard deviation on DeltaQTcF (ms)")
  })

  model({
    # ==================================================================
    # 1. Per-individual intercept and slope. Both random effects are
    #    additive on the linear scale (Mukker 2026 Equation 1).
    # ==================================================================
    e0_i    <- e0 + etae0
    slope_i <- slope + etaslope

    # ==================================================================
    # 2. Nominal-post-dose-time class effect. T_LASTDOSE carries the
    #    nominal scheduled post-dose hour; it is binned to the nearest
    #    of the four scheduled levels so that protocol-window slippage
    #    in a user-supplied column cannot fall through to the wrong
    #    class. Comparisons evaluate to 0/1 in rxode2, and the four
    #    indicators are mutually exclusive and exhaustive, so the sum
    #    selects exactly one shift.
    # ==================================================================
    i_tad0h <- (T_LASTDOSE < 0.5)
    i_tad1h <- (T_LASTDOSE >= 0.5) * (T_LASTDOSE < 1.5)
    i_tad2h <- (T_LASTDOSE >= 1.5) * (T_LASTDOSE < 2.5)
    i_tad3h <- (T_LASTDOSE >= 2.5)

    tad_effect <-
      e_tad0h_e0 * i_tad0h +
      e_tad1h_e0 * i_tad1h +
      e_tad2h_e0 * i_tad2h +
      e_tad3h_e0 * i_tad3h

    # ==================================================================
    # 3. Centered baseline-QTcF term. The centering constant is the
    #    PUBLISHED cohort mean baseline QTcF, 422 ms (Mukker 2026
    #    Table 1 'Mean (SD) baseline QTcF interval, ms' = 422 (21.0),
    #    which the Equation 1 legend names QTcF0). No assumption is
    #    needed here, in contrast to Darpo_2014_racSotalol_QTcF.R where
    #    the cohort median was unpublished.
    # ==================================================================
    qtc_bl_ref      <- 422
    qtc_bl_centered <- QTC_BL - qtc_bl_ref

    # ==================================================================
    # 4. Predicted change from baseline in QTcF (ms).
    #
    #    Note on reproducing Mukker 2026 Table 3: the quantity tabulated
    #    there is DeltaQTcF_drug, defined as the tuvusertib effect
    #    'corrected for intercept, baseline, and time' -- that is, the
    #    slope term alone. Evaluate this model at the concentration of
    #    interest and subtract its own prediction at zero concentration
    #    to recover it; do not read Table 3 against the full prediction
    #    below. The validation vignette does exactly that.
    # ==================================================================
    QTcF <- e0_i + slope_i * CP_TUVUSERTIB_NGML + tad_effect +
      e_qtc_bl_e0 * qtc_bl_centered

    QTcF ~ add(addSd)
  })
}
