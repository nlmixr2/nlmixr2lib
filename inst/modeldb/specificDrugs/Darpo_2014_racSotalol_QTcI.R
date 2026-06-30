Darpo_2014_racSotalol_QTcI <- function() {
  description <- paste(
    "Concentration-QTcI linear mixed-effects PD model for oral rac-sotalol",
    "in 39 healthy adults (28 men, 11 women) following a single 160 mg",
    "oral dose on day 1, with day-0 ECGs serving as the time-matched",
    "baseline. The endpoint is change from baseline in the individually-",
    "corrected QT interval (DeltaQTcI, ms). The full model is:",
    "DeltaQTcI = e0 + e_sexf_e0 * SEXF",
    "+ (slope + e_sexf_slope * SEXF) * CP_RACSOTALOL_UGML",
    "+ e_qtc_bl_e0 * (QTC_BL - 390),",
    "with male and median-baseline (QTC_BL = 390 ms) reference. Reference",
    "parameter values (male, median-baseline subject): e0 = -3.2 ms, slope",
    "= 23 ms per ug/mL, female increments +11.1 ms on the intercept and",
    "+7 ms per ug/mL on the slope, centred-baseline-QTcI coefficient",
    "-0.70 ms/ms. PD-only model: rac-sotalol plasma concentration is",
    "supplied as a time-varying covariate CP_RACSOTALOL_UGML (ug/mL).",
    "The source publication does not fit a population PK model -- the",
    "pharmacokinetic profile was characterised by NCA only (Cmax 1.4",
    "ug/mL men, 1.8 ug/mL women; AUC_inf 14.6 men, 17.1 women); users",
    "wishing to drive the PD model from a simulated PK source must",
    "supply their own concentration trajectory (no rac-sotalol popPK",
    "model exists in the nlmixr2lib registry). Companion model file",
    "Darpo_2014_racSotalol_QTcF.R reports the same structure for the",
    "Fridericia-corrected QT interval.",
    sep = " "
  )

  reference <- paste(
    "Darpo B, Karnad DR, Badilini F, Florian J, Garnett CE, Kothari S,",
    "Panicker GK, Sarapa N. (2014). Are women more susceptible than men",
    "to drug-induced QT prolongation? Concentration-QTc modelling in a",
    "phase 1 study with oral rac-sotalol.",
    "British Journal of Clinical Pharmacology 77(3):522-531.",
    "doi:10.1111/bcp.12201.",
    sep = " "
  )

  vignette <- "Darpo_2014_racSotalol"

  units <- list(
    time          = "h",
    dosing        = "(none; PD-only model fed by an external rac-sotalol plasma-concentration covariate)",
    concentration = "(observation DeltaQTcI is the change from time-matched day-0 baseline in the individually-corrected QT interval, ms; driving covariate CP_RACSOTALOL_UGML is in ug/mL)"
  )

  covariateData <- list(
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Source paper uses 'gender' as a categorical fixed effect with",
        "the male group as the implicit reference (Darpo 2014 Methods",
        "'Statistical analysis' paragraph and Table 1; the 'Female",
        "gender' row of Table 1 is the additive intercept shift, and",
        "the 'Concentration x Female gender interaction' row is the",
        "additive slope shift).",
        "Cohort composition: 11 women / 28 men out of 39 subjects."
      ),
      source_name        = "gender (paper notation; female = 1)"
    ),
    CP_RACSOTALOL_UGML = list(
      description        = "Instantaneous rac-sotalol plasma concentration at the time of each PD observation, supplied as a time-varying covariate from observed plasma samples or an upstream PK source.",
      units              = "ug/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying per event row. Drives the linear concentration-",
        "DeltaQTcI expression DeltaQTcI = e0 + slope *",
        "CP_RACSOTALOL_UGML (plus sex and centered-baseline covariate",
        "terms).",
        "In Darpo 2014 this was the observed rac-sotalol plasma",
        "concentration determined by validated HPLC (PPD Development,",
        "LLOQ 10 ng/mL, inter-day CV < 8.4%; Darpo 2014 Methods 'Study",
        "outline') at the 15 ECG time points (1, 1.5, 2, 2.5, 3, 3.5,",
        "4, 4.5, 5, 6, 8, 10, 13, 16, 22.5 h post 160 mg oral dose).",
        "Reference values observed: Tmax mean 2.8 h (men) / 2.9 h",
        "(women); Cmax 1.4 ug/mL (men, range 0.9-1.9) / 1.8 ug/mL",
        "(women, range 1.1-2.8); AUC_inf 14.6 (men) / 17.1 (women)",
        "ug/mL*h (Darpo 2014 Results 'Rac-sotalol pharmacokinetic",
        "profile' and Figure 1).",
        "The slope coefficient is reported on the ug/mL scale (Darpo",
        "2014 Table 1 'Plasma concentration of rac-sotalol' row: 'ms",
        "per ug ml-1'), so CP_RACSOTALOL_UGML is supplied directly in",
        "ug/mL with no in-model unit rescaling needed.",
        "Set to 0 outside the drug-exposure window (the concentration-",
        "slope term then collapses to 0)."
      ),
      source_name        = "rac-sotalol plasma concentration"
    ),
    QTC_BL = list(
      description        = "Subject's pre-dose (day-0) individually-corrected QT interval (QTcI) baseline, treated as a per-subject time-fixed covariate. Used to compute the centered baseline-QTcI effect on the linear-mixed-effects intercept: e_qtc_bl_e0 * (QTC_BL - 390). Set QTC_BL = 390 ms for the typical (median-baseline) subject -- the centered term then collapses to 0 and the model returns the pure sex- and concentration-driven DeltaQTcI prediction.",
      units              = "ms",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. The source paper computed each",
        "subject's baseline QTcI as the mean of the day-0 ECG QTcI",
        "values (15 nominal time points, triplicate 10-second strips,",
        "Darpo 2014 Methods 'ECG assessment' and 'Statistical",
        "analysis').",
        "Centering reference: the paper centred the per-subject",
        "baseline QTcI on the population median ('Median-centred",
        "baseline QTcI for the population was calculated using the",
        "median of all baseline QTcI values from all subjects.",
        "Individual median-centred baseline values were obtained by",
        "subtracting each baseline QTcI value from this median value';",
        "Darpo 2014 Methods 'Statistical analysis' paragraph 5). The",
        "specific median value is not quoted in the paper. Reference",
        "values observed: QTcI day 0 ranged 384-394 ms in men (n=28)",
        "and 398-412 ms in women (n=11) (Darpo 2014 Results 'The",
        "effect of rac-sotalol on QTc' paragraph 1). The packaged",
        "model uses the rounded standard 390 ms as the centering",
        "reference (qtc_bl_ref); this matches the male-dominated",
        "cohort median to within the reported per-sex range. The",
        "assumption is documented in the vignette Errata.",
        "For QTcF instead of QTcI, see the companion model file",
        "Darpo_2014_racSotalol_QTcF.R, which uses the same QTC_BL",
        "canonical covariate but interprets it as the QTcF baseline",
        "and applies a 390 ms QTcF centering reference (the paper's",
        "ranges are similar across the two correction methods: 380-396",
        "ms in men, 393-410 ms in women)."
      ),
      source_name        = "baseline QTcI (day-0 individually-corrected QT interval)"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 39L,
    n_studies        = 1L,
    n_observations   = 3456L,
    age_range        = "18-45 years; mean 27 years",
    weight_range     = "47-108 kg; mean 74 kg",
    sex_female_pct   = round(100 * 11 / 39, 1),
    race_ethnicity   = NA_character_,
    disease_state    = paste(
      "Healthy young adult subjects meeting typical inclusion criteria",
      "for phase 1 studies. Studied at the Pharmacia Clinical Research",
      "Unit in Kalamazoo, MI under Good Clinical Practice (Darpo 2014",
      "Methods 'Study outline'). Body mass index range 18-31 kg/m^2",
      "(mean 24)."
    ),
    dose_range       = paste(
      "Single 160 mg oral dose of rac-sotalol (Betapace, Berlex",
      "Laboratories) given at 08:00 in the fasted state on day 1, with",
      "a separate baseline day 0 (no drug) used to construct the",
      "DeltaQTc reference. A 320 mg dose was planned for day 2 in the",
      "parent study but all 11 women were excluded from day-2 dosing",
      "by the protocol discontinuation criterion DeltaQTcF > 60 ms",
      "from the day-1 dose; day-2 data are not part of this PD",
      "extraction (Darpo 2014 Methods 'Study outline' paragraph 2)."
    ),
    regions          = "United States (Pharmacia Clinical Research Unit, Kalamazoo, Michigan).",
    notes            = paste(
      "Pooled day-0 (baseline) + day-1 (160 mg rac-sotalol) ECG data",
      "from 39 subjects: 3456 ECGs (33 missing, 49 with QT not",
      "measurable; 1.9% missing rate). Triplicate 10-second 12-lead",
      "ECG strips at 15 post-dose nominal time points: 1, 1.5, 2,",
      "2.5, 3, 3.5, 4, 4.5, 5, 6, 8, 10, 13, 16 and 22.5 h. ECGs",
      "recorded by Holter at 180 Hz and up-sampled to 1000 Hz",
      "before manual QT/RR measurement on the superimposed median",
      "beat (Darpo 2014 Methods 'ECG assessment').",
      "Heart-rate correction method: individual (QTcI = QT / RR^I),",
      "where I is each subject's own log-linear regression slope of",
      "ln(QT) on ln(RR) using all drug-free baseline-day ECGs",
      "(Darpo 2014 Methods 'Statistical analysis' paragraph 1; QTcI",
      "selected because its regression-line slope on the on-drug",
      "data was 0.016 ms/ms RR -- the lowest among the four tested",
      "correction methods QTcB, QTcF, QTcI, QTcN; Darpo 2014 Results",
      "'Correction for heart rate changes').",
      "Pre-specified base model structure (after AIC selection):",
      "linear concentration-DeltaQTcI relationship with both slope",
      "and intercept as fixed effects with additive between-subject",
      "variability on both (Darpo 2014 Methods 'Statistical",
      "analysis' paragraph 5, choice (i)). Covariates retained in",
      "the final structural model: median-centred baseline QTcI",
      "(on intercept), female gender (on both intercept and slope).",
      "Statistical software was SAS 9.2 (SAS Institute Inc., Cary,",
      "NC, USA)."
    )
  )

  ini({
    # ==================================================================
    # Linear concentration-DeltaQTcI model (Darpo 2014 Methods
    # 'Statistical analysis' paragraph 5 and Figure 6A caption /
    # Table 1 QTcI row). Source equation reported for men:
    #   DeltaQTcI = -3.2 - 0.7 * centered_baseline_QTcI
    #             + 23 * rac-sotalol_concentration_ug_per_mL
    # and for women:
    #   DeltaQTcI =  7.9 - 0.7 * centered_baseline_QTcI
    #             + 30 * rac-sotalol_concentration_ug_per_mL
    # Encoded here in a male-reference + additive-sex-effect form so
    # the male/female contrast is reproduced by setting SEXF = 0 or 1.
    # ==================================================================

    e0 <- -3.2
    label("Linear-mixed-effects intercept e0 on DeltaQTcI, male median-baseline reference (ms)")
    # Darpo 2014 Table 1 QTcI row 'Intercept (ms)' = -3.2 (SE 2.0;
    # P = 0.12). Not log-transformed because the intercept of a
    # DeltaQTc model can take negative values (close-to-zero baseline-
    # day-vs-baseline-day prediction). Canonical bare PD baseline
    # parameter `e0` is registered in parameter-names.md.

    lslope <- log(23)
    label("Linear concentration-DeltaQTcI slope, male reference (ms per ug/mL)")
    # Darpo 2014 Table 1 QTcI row 'Plasma concentration of rac-
    # sotalol (ms per ug ml-1)' = 23 (SE 1.7; P < 0.0001). Log-
    # transformed because the slope is positive (concentration-
    # related QTc prolongation). Bare name inside model() is `slope`.

    # ------------------------------------------------------------------
    # Covariate effects (Darpo 2014 Table 1 QTcI rows). All three are
    # estimated point estimates with reported SE / P-value; not fixed.
    # ------------------------------------------------------------------

    e_sexf_e0 <- 11.1
    label("Effect of female sex on e0 (additive, ms)")
    # Darpo 2014 Table 1 QTcI row 'Female gender (ms)' = 11.1 (SE 3.8;
    # P = 0.004). Reproduces female-reference intercept e0_female =
    # -3.2 + 11.1 = 7.9 ms (Figure 6A caption).

    e_sexf_slope <- 7
    label("Effect of female sex on slope (additive, ms per ug/mL)")
    # Darpo 2014 Table 1 QTcI row 'Concentration x Female gender
    # interaction (ms per ug ml-1)' = 7 (SE 2.8; P = 0.01).
    # Reproduces female-reference slope slope_female = 23 + 7 = 30
    # ms per ug/mL (Figure 6A caption).

    e_qtc_bl_e0 <- -0.70
    label("Effect of centered baseline QTcI on e0 (additive, ms per ms)")
    # Darpo 2014 Table 1 QTcI row 'Centred baseline QTcI (ms)' =
    # -0.70 (SE 0.05; P < 0.0001). The negative coefficient means
    # subjects with a higher-than-median baseline QTcI have a smaller
    # predicted DeltaQTcI (regression-to-the-mean shrinkage of the
    # change-from-baseline at high baseline).

    # ==================================================================
    # Inter-individual variability: Darpo 2014 Methods 'Statistical
    # analysis' paragraph 5 reports that the chosen base model includes
    # 'additive between-subject variability associated with slope and
    # intercept' (option (i) of the three evaluated structures), but
    # the source paper does NOT report the variance estimates for
    # omega_e0 or omega_slope. Per the standing operator policy on
    # unreported IIV, this file omits the eta declarations and ships a
    # typical-value-only model. The vignette Errata documents the gap.
    # ==================================================================

    # ==================================================================
    # Residual error: not reported numerically in Darpo 2014. Encoded
    # as fixed(0) per the standing operator policy on unreported
    # residual error (the model returns the deterministic typical-
    # value prediction). The vignette Errata documents the gap.
    # ==================================================================
    addSd <- fixed(0)
    label("Additive residual error standard deviation on DeltaQTcI (ms; FIXED AT ZERO - not reported in source)")
  })

  model({
    # ==================================================================
    # 1. Per-subject typical-value intercept and slope after additive
    #    sex covariate effects (Darpo 2014 Figure 6A caption).
    # ==================================================================
    e0_i    <- e0 + e_sexf_e0 * SEXF
    slope_i <- exp(lslope) + e_sexf_slope * SEXF

    # ==================================================================
    # 2. Centered baseline-QTcI term. The paper centered each
    #    subject's baseline QTcI on the cohort median (Darpo 2014
    #    Methods 'Statistical analysis' paragraph 5). The specific
    #    median is not quoted in the paper; the packaged model uses
    #    the rounded standard 390 ms, which matches the male-
    #    dominated cohort (28 men, baseline QTcI 384-394 ms; 11
    #    women, 398-412 ms) to within the reported per-sex range.
    #    See vignette Errata.
    # ==================================================================
    qtc_bl_ref <- 390
    qtc_bl_centered <- QTC_BL - qtc_bl_ref

    # ==================================================================
    # 3. Linear concentration-DeltaQTcI prediction. Output is the
    #    change from time-matched day-0 baseline in QTcI, in ms.
    # ==================================================================
    QTcI <- e0_i + slope_i * CP_RACSOTALOL_UGML + e_qtc_bl_e0 * qtc_bl_centered

    QTcI ~ add(addSd)
  })
}
