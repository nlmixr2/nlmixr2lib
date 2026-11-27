Huang_2026_busulfan_efs_bu4 <- function() {
  description <- paste(
    "Cox proportional-hazards exposure-response model for event-free survival",
    "(EFS; composite of graft failure, disease relapse or nonrelapse mortality)",
    "after allogeneic haematopoietic stem cell transplantation in the 129",
    "Chinese pediatric patients receiving 4-day (16-dose) IV busulfan-based",
    "myeloablative conditioning (Huang 2026, Table S2c). The exposure driver",
    "AUC_BUSULFAN_COURSE is the cumulative 4-day busulfan AUC from the",
    "companion popPK model modellib('Huang_2026_busulfan'); it enters through",
    "a 70-90 mg*h/L reference window chosen from a restricted cubic spline.",
    "Subjects below the window have a 2.94-fold and subjects above it a",
    "3.22-fold hazard of an EFS event. No other covariate was retained in the",
    "multivariable model. Semiparametric Cox regression: the baseline hazard is",
    "unspecified and not reported, so the model returns the RELATIVE hazard hr",
    "against the 70-90 mg*h/L reference group. Algebraic and deterministic (no",
    "ODE state, no drug input, no IIV, no residual error)."
  )
  reference <- paste(
    "Huang C, Du X, Xue L, Song L, Liu X, Chen Z, Wang X, Qin M, Sun X, Wang D,",
    "Yang J, Li D, Jiang H, Ding W, Qu Y, Zhu Z, Zhao L, Yuan J, Dong J, Wu D,",
    "Hu S, Miao L. Association of Busulfan Exposure with Event-Free Survival in",
    "a Chinese Pediatric Cohort: A Multicenter, Prospective Cohort Study.",
    "Drug Des Devel Ther. 2026;20:593835. doi:10.2147/DDDT.S593835.",
    "Hazard ratios from Supplementary Table S2c (multivariate column).",
    sep = " "
  )
  vignette <- "Huang_2026_busulfan"
  units <- list(
    time = "n/a (semiparametric Cox regression; the baseline hazard and hence the time scale are left unspecified; follow-up was measured from transplantation)",
    dosing = "n/a (exposure-response model; busulfan exposure enters as the AUC_BUSULFAN_COURSE covariate, not as a dose record)",
    concentration = "hr (hazard of an EFS event relative to AUC_BUSULFAN_COURSE within 70-90 mg*h/L, unitless); not a drug concentration"
  )

  covariateData <- list(
    AUC_BUSULFAN_COURSE = list(
      description = "Cumulative busulfan AUC over the whole 4-day conditioning course (16 doses), computed as dose divided by the individual empirical-Bayes clearance summed over the course",
      units = "mg*h/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Bayesian individual CL from modellib('Huang_2026_busulfan') (Methods, 'Pharmacokinetics of Busulfan'). The paper's cAUC4 days; median 70.5 mg*h/L, range 30.2-109.9 in the 4-day cohort (Results). The model consumes it only through the indicators AUC_BUSULFAN_COURSE < cauc_lo and > cauc_hi.",
      source_name = "cAUC4 days"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 129L,
    n_studies = 1L,
    n_events = 31L,
    age_range = "0.5-16.9 years (median 4.5)",
    weight_range = "6.0-80.0 kg (median 17.0)",
    sex_female_pct = 31.78,
    race_ethnicity = "Chinese",
    disease_state = "Pediatric first allogeneic HSCT with a 4-day IV busulfan-based myeloablative conditioning regimen; 66.7% malignant disease. See modellib('Huang_2026_busulfan') for the full cohort description.",
    dose_range = "0.8-1.2 mg/kg per dose, 2 h IV infusion four times daily for 4 days",
    regions = "China (five pediatric HSCT centres)",
    notes = "31 EFS events (Table S2c). cAUC4 days 70-90 mg*h/L: 45 patients (5 events); < 70: 60 (18 events); > 90: 24 (8 events). Kaplan-Meier 2-year EFS 88.9%, 71.0% and 66.1% (Figure 4C)."
  )

  ini({
    # Window edges read off the restricted cubic spline (Figure 3C), chosen
    # rather than estimated with uncertainty.
    cauc_lo <- fixed(70); label("Lower edge of the cumulative 4-day busulfan AUC reference window (mg*h/L)") # Results 'EFS'; Figure 3C; Table S2c
    cauc_hi <- fixed(90); label("Upper edge of the cumulative 4-day busulfan AUC reference window (mg*h/L)") # Results 'EFS'; Figure 3C; Table S2c

    # Log-CI midpoints: (log(1.09) + log(7.92)) / 2 = 1.0778 vs log(2.94) = 1.0784;
    # (log(1.05) + log(9.85)) / 2 = 1.1681 vs log(3.22) = 1.1694.
    e_auc_busulfan_course_low_haz <- log(2.94); label("Log hazard ratio for an EFS event, cAUC4 days < 70 vs 70-90 mg*h/L (log scale; HR 2.94)") # Table S2c multivariate HR 2.94 (95% CI 1.09-7.92), P = 0.03
    e_auc_busulfan_course_high_haz <- log(3.22); label("Log hazard ratio for an EFS event, cAUC4 days > 90 vs 70-90 mg*h/L (log scale; HR 3.22)") # Table S2c multivariate HR 3.22 (95% CI 1.05-9.85), P = 0.04
  })

  model({
    # Reference group: 70 <= cAUC4 days <= 90 mg*h/L ('lower than 70' and
    # 'higher than 90' are the comparison groups, Results 'EFS').
    caucLow <- (AUC_BUSULFAN_COURSE < cauc_lo)
    caucHigh <- (AUC_BUSULFAN_COURSE > cauc_hi)

    lhr <- e_auc_busulfan_course_low_haz * caucLow +
      e_auc_busulfan_course_high_haz * caucHigh
    hr <- exp(lhr)
  })
}
