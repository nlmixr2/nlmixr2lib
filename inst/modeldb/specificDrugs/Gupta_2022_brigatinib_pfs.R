Gupta_2022_brigatinib_pfs <- function() {
  description <- paste0(
    "Cox proportional-hazards exposure-efficacy model for ",
    "progression-free survival by blinded independent review (BIRC) as ",
    "a function of the brigatinib daily AUC time-averaged between the ",
    "last two disease-assessment scans preceding progression or ",
    "censoring, a STATIC (time-fixed) covariate, in adults with ",
    "ALK-inhibitor-naive ALK-positive advanced non-small cell lung ",
    "cancer treated first-line with brigatinib 180 mg orally once daily ",
    "after a 7-day lead-in at 90 mg (Gupta 2022, ALTA-1L brigatinib ",
    "arm, n = 123). The model returns the RELATIVE hazard exp(log(1.03) ",
    "* AUC_BRIG_SCAN), the printed hazard ratio 1.03 (95 percent CI ",
    "1.01-1.05) being per 1 ug*h/mL/day. The relationship is ",
    "statistically significant (p = 0.01) but runs in the ",
    "pharmacologically inconsistent direction -- a higher exposure is ",
    "associated with a HIGHER hazard of progression or death -- and ",
    "Gupta 2022 attributes it to the static metric: patients who stay ",
    "on treatment longer accumulate dose reductions, so a long PFS is ",
    "attached to a lower time-averaged exposure. With the time-varying ",
    "metrics that account for dose changes ",
    "(Gupta_2022_brigatinib_pfs_scan_tv and ",
    "Gupta_2022_brigatinib_pfs_daily_tv) the effect is not significant. ",
    "NO BASELINE HAZARD IS ENCODED: a Cox regression is semiparametric ",
    "and leaves h0(t) unspecified, so the absolute survivor function is ",
    "a quantity the fit never produced rather than an unreported ",
    "parameter. There is no PK layer and no ODE: the exposure metric is ",
    "supplied as a data column, derived in the source analysis from the ",
    "individual Bayesian CL/F of the population PK model packaged as ",
    "Gupta_2021_brigatinib. Twenty companion exposure-response models ",
    "in the Gupta_2022_brigatinib_* family."
  )
  reference <- paste(
    "Gupta N, Reckamp KL, Camidge DR, Kleijn HJ, Ouerdani A, Bellanti F,",
    "Maringwa J, Hanley MJ, Wang S, Zhang P, Venkatakrishnan K. (2022).",
    "Population pharmacokinetic and exposure-response analyses from ALTA-1L:",
    "Model-based analyses supporting the brigatinib dose in ALK-positive NSCLC.",
    "Clin Transl Sci 15:1143-1154.",
    "doi:10.1111/cts.13231.",
    "Individual brigatinib exposures derive from the previously published",
    "population PK model applied without modification by Bayesian",
    "re-estimation; see modellib('Gupta_2021_brigatinib').",
    "Supporting Figures S6 and S7 are Supporting Information files s002 and",
    "s005; Methods S1 is file s007."
  )
  vignette <- "Gupta_2022_brigatinib_exposure_response"
  units <- list(
    time = "n/a (semiparametric Cox relative hazard; the baseline hazard and therefore the time scale are unspecified by the fit)",
    dosing = "n/a (no dose events; exposure enters as the AUC_BRIG_SCAN covariate column)",
    concentration = "hr (relative hazard of progression or death, unitless; also lhr, the log relative hazard)"
  )

  covariateData <- list(
    AUC_BRIG_SCAN = list(
      description = "Individual brigatinib daily AUC averaged over the disease-assessment scan interval -- in the static models, the interval between the last two scans preceding the event, best confirmed response or censoring. Supplied as data: this model has no PK layer, and the source analysis derived it from the individual Bayesian CL/F of the Gupta 2021 population PK model and the actual dosing history (Gupta 2022 Methods S1).",
      units = "ug*h/mL/day",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Daily AUC on day i is the increment of the cumulative AUC,",
        "AUCD(i) = CAUC(i + 1) - CAUC(i), with CAUC(i) = cumulative dose to",
        "day i / individual CL/F (Gupta 2022 Methods S1); the metric is the",
        "mean of AUCD over the scan interval, so dose reductions and",
        "interruptions inside the interval lower it. It is NOT a",
        "steady-state AUC(0-24) at the labelled dose, although for an",
        "uninterrupted 180 mg once-daily course the two coincide",
        "(180 mg / CL/F). Enters LINEARLY and uncentred.",
        "Gupta 2022 Figure 2a gives the PFS analysis-set median (range)",
        "as 15.2 (0, 77.0) and the quartile medians as 9.7, 13.0, 17.1 and",
        "28.3 ug*h/mL/day.",
        "In this model the exposure enters with log hazard ratio log(1.03) per ug*h/mL/day."
      ),
      source_name = "Time-averaged AUC between last two scans (ug.h/mL/day)"
    )
  )

  covariatesDataExcluded <- list(
    ECOG_GE1 = list(
      description = "Baseline ECOG performance status of 1 or more.",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Gupta 2022 Results: in the stepwise covariate analysis only ECOG",
        "performance status was a significant predictor of progression or",
        "death (p < 0.001), patients with ECOG 1 progressing faster than",
        "those with ECOG 0, and after its addition the exposure effect",
        "remained significant. The paper prints neither the ECOG",
        "coefficient nor the ECOG-adjusted exposure hazard ratio, so this",
        "file encodes the reported exposure-only model and the ECOG term",
        "is documented here rather than encoded. ECOG 0 / 1 / 2 in",
        "39.8 / 56.1 / 4.1 percent of the 123 patients (Table 1); how the",
        "five ECOG 2 patients were coded is not stated."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 123L,
    n_studies = 1L,
    n_observations = "123 right-censored BIRC progression-free-survival times, one per patient; 59 events (Figure 2a)",
    age_range = "27-85 years (median 57; Gupta 2022 Table 1, PK-evaluable ALTA-1L brigatinib arm, n = 123)",
    weight_range = "43-111 kg (median 67)",
    sex_female_pct = 48.8,
    race_ethnicity = c(White = 55.3, Asian = 43.1, Other = 1.6),
    disease_state = "ALK-inhibitor-naive ALK-positive advanced non-small cell lung cancer, first-line treatment in the phase 3 ALTA-1L trial (NCT02737501), second interim analysis (data cut-off 28 June 2019). ECOG performance status 0 / 1 / 2 in 39.8 / 56.1 / 4.1 percent; prior chemotherapy 28.5 percent; brain metastases at baseline 29.3 percent; baseline albumin median 41 g/L (range 24-48).",
    dose_range = "Brigatinib 180 mg orally once daily after a 7-day lead-in at 90 mg once daily (28-day cycles); the crizotinib comparator arm (250 mg twice daily) contributes no exposure and is not modelled.",
    regions = "Multiregional (Gupta 2022 Discussion: ALTA-1L was a multiregional clinical trial); per-region counts not reported.",
    notes = paste0(
      "Of 137 intent-to-treat brigatinib patients, 13 had no ",
      "quantifiable brigatinib concentration and one was not dosed, ",
      "leaving 123 for the population PK re-estimation (1069 samples) ",
      "and the exposure-response analyses of PFS, ORR (both by ",
      "blinded independent review) and safety. Median BIRC PFS 24 ",
      "months (95 percent CI 21-NA) in the 123-patient analysis set ",
      "(Figure 2a)."
    )
  )

  ini({
    # ==================================================================
    # Gupta 2022 Results, Exposure-efficacy analysis: static time-averaged
    # AUC between the last two disease assessments, hazard ratio 1.03 (95
    # percent CI 1.01-1.05), p = 0.01, per 1 ug*h/mL/day. The coefficient
    # below is the natural log of the printed hazard ratio, the scale a Cox
    # model estimates on; the paper prints no standard error and no
    # coefficient, so the two-decimal hazard ratio is the only precision
    # available (log(1.03) = 0.02956).
    #
    # NO BASELINE HAZARD IS ENCODED. A Cox regression leaves h0(t)
    # unspecified, so this model returns the relative hazard only; Gupta
    # 2022 Figure 2a gives the Kaplan-Meier PFS curves by exposure quartile
    # as the empirical description of the absolute time course. The
    # regression was fitted with standard survival software in R, with no
    # variance components, so there is no IIV, no residual error and no
    # observation endpoint.
    # ==================================================================
    e_auc_haz <- log(1.03); label("Log hazard ratio for progression or death per 1 ug*h/mL/day increase in AUC_BRIG_SCAN (log scale)") # Gupta 2022 Results, Exposure-efficacy analysis: hazard ratio 1.03 (95 percent CI 1.01-1.05), p = 0.01
  })

  model({
    # Relative hazard against a hypothetical zero-exposure patient. Under
    # proportional hazards the subject hazard is h(t) = h0(t) * hr for any
    # baseline the user supplies; only ratios of hr between two exposures
    # are meaningful.
    lhr <- e_auc_haz * AUC_BRIG_SCAN

    hr <- exp(lhr)
  })
}
