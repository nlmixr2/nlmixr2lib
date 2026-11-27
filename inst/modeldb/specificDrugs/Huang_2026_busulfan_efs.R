Huang_2026_busulfan_efs <- function() {
  description <- paste(
    "Cox proportional-hazards exposure-response model for event-free survival",
    "(EFS; composite of graft failure, disease relapse or nonrelapse mortality)",
    "after allogeneic haematopoietic stem cell transplantation in all 455",
    "Chinese pediatric patients receiving 3- or 4-day IV busulfan-based",
    "myeloablative conditioning (Huang 2026, Table S2a). The exposure driver",
    "AUC_BUSULFAN_FIRST is the first-dose busulfan AUC (dose / individual CL)",
    "from the companion popPK model modellib('Huang_2026_busulfan'); it enters",
    "only through the 4.5 mg*h/L cut point chosen from a restricted cubic",
    "spline, and subjects below it have a 1.63-fold hazard of an EFS event.",
    "No other covariate was retained in the multivariable model. Because the",
    "fit is a semiparametric Cox regression, the baseline hazard is unspecified",
    "and not reported; the model returns the RELATIVE hazard hr against the",
    "AUC >= 4.5 mg*h/L reference group, not an absolute survivor function.",
    "Algebraic and deterministic (no ODE state, no drug input, no IIV, no",
    "residual error)."
  )
  reference <- paste(
    "Huang C, Du X, Xue L, Song L, Liu X, Chen Z, Wang X, Qin M, Sun X, Wang D,",
    "Yang J, Li D, Jiang H, Ding W, Qu Y, Zhu Z, Zhao L, Yuan J, Dong J, Wu D,",
    "Hu S, Miao L. Association of Busulfan Exposure with Event-Free Survival in",
    "a Chinese Pediatric Cohort: A Multicenter, Prospective Cohort Study.",
    "Drug Des Devel Ther. 2026;20:593835. doi:10.2147/DDDT.S593835.",
    "Hazard ratio from Supplementary Table S2a (multivariate column).",
    sep = " "
  )
  vignette <- "Huang_2026_busulfan"
  units <- list(
    time = "n/a (semiparametric Cox regression; the baseline hazard and hence the time scale are left unspecified; follow-up was measured from transplantation)",
    dosing = "n/a (exposure-response model; busulfan exposure enters as the AUC_BUSULFAN_FIRST covariate, not as a dose record)",
    concentration = "hr (hazard of an EFS event relative to AUC_BUSULFAN_FIRST >= 4.5 mg*h/L, unitless); not a drug concentration"
  )

  covariateData <- list(
    AUC_BUSULFAN_FIRST = list(
      description = "First-dose busulfan AUC from time zero to infinity, computed as the first dose divided by the individual empirical-Bayes clearance",
      units = "mg*h/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Bayesian individual CL from modellib('Huang_2026_busulfan') (Methods, 'Pharmacokinetics of Busulfan'). Cohort median 4.5 mg*h/L, range 1.8-10.5 (Results). The model consumes it only through the indicator AUC_BUSULFAN_FIRST < auc_cut.",
      source_name = "AUC first"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 455L,
    n_studies = 1L,
    n_events = 88L,
    age_range = "0.5-17.0 years (median 7.5)",
    weight_range = "5.0-96.5 kg (median 24.5)",
    sex_female_pct = 36.7,
    race_ethnicity = "Chinese",
    disease_state = "Pediatric first allogeneic HSCT with IV busulfan-based myeloablative conditioning (3-day, n = 326; 4-day, n = 129); 73.9% malignant disease. See modellib('Huang_2026_busulfan') for the full cohort description.",
    dose_range = "0.8-1.2 mg/kg per dose, 2 h IV infusion four times daily for 3 or 4 days",
    regions = "China (five pediatric HSCT centres)",
    notes = "88 EFS events (Table S2a); median follow-up 2.2 years. 231 patients had AUC first >= 4.5 mg*h/L (35 events) and 224 had < 4.5 mg*h/L (53 events). Kaplan-Meier 2-year EFS 84.6% vs 76.4% (Figure 4A)."
  )

  ini({
    # Cut point read off the restricted cubic spline (Figure 3A; Methods:
    # 'thresholds ... defined at points at which the HR was greater than 1'),
    # chosen rather than estimated with uncertainty.
    auc_cut <- fixed(4.5); label("First-dose busulfan AUC cut point separating the low-exposure group from the reference (mg*h/L)") # Results 'EFS'; Figure 3A; Table S2a

    # Log-CI midpoint (log(1.06) + log(2.50)) / 2 = 0.4873 matches log(1.63) = 0.4886.
    e_auc_busulfan_first_low_haz <- log(1.63); label("Log hazard ratio for an EFS event, AUC first < 4.5 vs >= 4.5 mg*h/L (log scale; HR 1.63)") # Table S2a multivariate HR 1.63 (95% CI 1.06-2.50), P = 0.02
  })

  model({
    # Reference group: AUC first >= 4.5 mg*h/L (Table S2a).
    aucLow <- (AUC_BUSULFAN_FIRST < auc_cut)

    lhr <- e_auc_busulfan_first_low_haz * aucLow
    hr <- exp(lhr)
  })
}
