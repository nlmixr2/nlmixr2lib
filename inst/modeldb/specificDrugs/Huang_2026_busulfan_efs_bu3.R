Huang_2026_busulfan_efs_bu3 <- function() {
  description <- paste(
    "Cox proportional-hazards exposure-response model for event-free survival",
    "(EFS; composite of graft failure, disease relapse or nonrelapse mortality)",
    "after allogeneic haematopoietic stem cell transplantation in the 326",
    "Chinese pediatric patients receiving 3-day (12-dose) IV busulfan-based",
    "myeloablative conditioning (Huang 2026, Table S2b). The exposure driver",
    "AUC_BUSULFAN_COURSE is the cumulative 3-day busulfan AUC from the",
    "companion popPK model modellib('Huang_2026_busulfan'); it enters only",
    "through the 54 mg*h/L cut point chosen from a restricted cubic spline, and",
    "subjects below it have a 1.73-fold hazard of an EFS event. No other",
    "covariate was retained in the multivariable model. Semiparametric Cox",
    "regression: the baseline hazard is unspecified and not reported, so the",
    "model returns the RELATIVE hazard hr against the cumulative AUC >= 54",
    "mg*h/L reference group. Algebraic and deterministic (no ODE state, no drug",
    "input, no IIV, no residual error)."
  )
  reference <- paste(
    "Huang C, Du X, Xue L, Song L, Liu X, Chen Z, Wang X, Qin M, Sun X, Wang D,",
    "Yang J, Li D, Jiang H, Ding W, Qu Y, Zhu Z, Zhao L, Yuan J, Dong J, Wu D,",
    "Hu S, Miao L. Association of Busulfan Exposure with Event-Free Survival in",
    "a Chinese Pediatric Cohort: A Multicenter, Prospective Cohort Study.",
    "Drug Des Devel Ther. 2026;20:593835. doi:10.2147/DDDT.S593835.",
    "Hazard ratio from Supplementary Table S2b (multivariate column).",
    sep = " "
  )
  vignette <- "Huang_2026_busulfan"
  units <- list(
    time = "n/a (semiparametric Cox regression; the baseline hazard and hence the time scale are left unspecified; follow-up was measured from transplantation)",
    dosing = "n/a (exposure-response model; busulfan exposure enters as the AUC_BUSULFAN_COURSE covariate, not as a dose record)",
    concentration = "hr (hazard of an EFS event relative to AUC_BUSULFAN_COURSE >= 54 mg*h/L, unitless); not a drug concentration"
  )

  covariateData <- list(
    AUC_BUSULFAN_COURSE = list(
      description = "Cumulative busulfan AUC over the whole 3-day conditioning course (12 doses), computed as dose divided by the individual empirical-Bayes clearance summed over the course",
      units = "mg*h/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Bayesian individual CL from modellib('Huang_2026_busulfan') (Methods, 'Pharmacokinetics of Busulfan'). The paper's cAUC3 days; median 54.2 mg*h/L, range 26.5-119.0 in the 3-day cohort (Results). The model consumes it only through the indicator AUC_BUSULFAN_COURSE < cauc_cut.",
      source_name = "cAUC3 days"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 326L,
    n_studies = 1L,
    n_events = 57L,
    age_range = "0.5-17.0 years (median 8.6)",
    weight_range = "5.0-96.5 kg (median 27.0)",
    sex_female_pct = 38.65,
    race_ethnicity = "Chinese",
    disease_state = "Pediatric first allogeneic HSCT with a 3-day IV busulfan-based myeloablative conditioning regimen; 76.7% malignant disease. See modellib('Huang_2026_busulfan') for the full cohort description.",
    dose_range = "0.8-1.2 mg/kg per dose, 2 h IV infusion four times daily for 3 days",
    regions = "China (five pediatric HSCT centres)",
    notes = "57 EFS events (Table S2b). 166 patients had cAUC3 days >= 54 mg*h/L (22 events) and 160 had < 54 mg*h/L (35 events). Kaplan-Meier 2-year EFS 86.4% vs 77.7% (Figure 4B)."
  )

  ini({
    # Cut point read off the restricted cubic spline (Figure 3B), chosen
    # rather than estimated with uncertainty.
    cauc_cut <- fixed(54); label("Cumulative 3-day busulfan AUC cut point separating the low-exposure group from the reference (mg*h/L)") # Results 'EFS'; Figure 3B; Table S2b

    # Log-CI midpoint (log(1.01) + log(2.94)) / 2 = 0.5442 matches log(1.73) = 0.5481.
    e_auc_busulfan_course_low_haz <- log(1.73); label("Log hazard ratio for an EFS event, cAUC3 days < 54 vs >= 54 mg*h/L (log scale; HR 1.73)") # Table S2b multivariate HR 1.73 (95% CI 1.01-2.94), P = 0.04
  })

  model({
    # Reference group: cAUC3 days >= 54 mg*h/L (Table S2b).
    caucLow <- (AUC_BUSULFAN_COURSE < cauc_cut)

    lhr <- e_auc_busulfan_course_low_haz * caucLow
    hr <- exp(lhr)
  })
}
