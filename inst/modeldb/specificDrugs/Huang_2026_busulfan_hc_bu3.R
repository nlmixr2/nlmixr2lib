Huang_2026_busulfan_hc_bu3 <- function() {
  description <- paste(
    "Cox proportional-hazards exposure-response model for grade 3-4",
    "hemorrhagic cystitis (HC) within 100 days of allogeneic haematopoietic",
    "stem cell transplantation in the 326 Chinese pediatric patients receiving",
    "3-day (12-dose) IV busulfan-based myeloablative conditioning (Huang 2026,",
    "Table S4b). The exposure driver AUC_BUSULFAN_COURSE is the cumulative",
    "3-day busulfan AUC from the companion popPK model",
    "modellib('Huang_2026_busulfan'); subjects below the 54 mg*h/L cut point",
    "have a 0.38-fold hazard of grade 3-4 HC (borderline, P = 0.05). No other",
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
    "Hazard ratio from Supplementary Table S4b (multivariate column).",
    sep = " "
  )
  vignette <- "Huang_2026_busulfan"
  units <- list(
    time = "n/a (semiparametric Cox regression; the baseline hazard and hence the time scale are left unspecified; follow-up was measured from transplantation over 100 days)",
    dosing = "n/a (exposure-response model; busulfan exposure enters as the AUC_BUSULFAN_COURSE covariate, not as a dose record)",
    concentration = "hr (hazard of grade 3-4 hemorrhagic cystitis relative to AUC_BUSULFAN_COURSE >= 54 mg*h/L, unitless); not a drug concentration"
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
    n_events = 22L,
    age_range = "0.5-17.0 years (median 8.6)",
    weight_range = "5.0-96.5 kg (median 27.0)",
    sex_female_pct = 38.65,
    race_ethnicity = "Chinese",
    disease_state = "Pediatric first allogeneic HSCT with a 3-day IV busulfan-based myeloablative conditioning regimen; 76.7% malignant disease. See modellib('Huang_2026_busulfan') for the full cohort description.",
    dose_range = "0.8-1.2 mg/kg per dose, 2 h IV infusion four times daily for 3 days",
    regions = "China (five pediatric HSCT centres)",
    notes = "22 grade 3-4 HC events (Table S4b). cAUC3 days >= 54 mg*h/L: 166 patients (16 events); < 54: 160 (6 events). Grade 3-4 HC occurred only in patients with malignant disease who received ATG (Results 'HC')."
  )

  ini({
    # Cut point shared with the 3-day EFS analysis (Figure 3B), chosen rather
    # than estimated with uncertainty.
    cauc_cut <- fixed(54); label("Cumulative 3-day busulfan AUC cut point separating the low-exposure group from the reference (mg*h/L)") # Results 'HC'; Table S4b

    # Log-CI midpoint (log(0.15) + log(0.98)) / 2 = -0.9587 matches log(0.38) = -0.9676.
    e_auc_busulfan_course_low_haz <- log(0.38); label("Log hazard ratio for grade 3-4 HC, cAUC3 days < 54 vs >= 54 mg*h/L (log scale; HR 0.38)") # Table S4b multivariate HR 0.38 (95% CI 0.15-0.98), P = 0.05
  })

  model({
    # Reference group: cAUC3 days >= 54 mg*h/L (Table S4b).
    caucLow <- (AUC_BUSULFAN_COURSE < cauc_cut)

    lhr <- e_auc_busulfan_course_low_haz * caucLow
    hr <- exp(lhr)
  })
}
