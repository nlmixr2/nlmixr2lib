Huang_2026_busulfan_nrm_bu3 <- function() {
  description <- paste(
    "Cox proportional-hazards exposure-response model for nonrelapse mortality",
    "(NRM) after allogeneic haematopoietic stem cell transplantation in the 326",
    "Chinese pediatric patients receiving 3-day (12-dose) IV busulfan-based",
    "myeloablative conditioning (Huang 2026, Table S6). The exposure driver",
    "AUC_BUSULFAN_COURSE is the cumulative 3-day busulfan AUC from the",
    "companion popPK model modellib('Huang_2026_busulfan'); subjects below the",
    "54 mg*h/L cut point have a 2.48-fold hazard of NRM after adjustment for",
    "age group (3-6, 6-12 and > 12 years against 0-3 years). Semiparametric Cox",
    "regression: the baseline hazard is unspecified and not reported, so the",
    "model returns the RELATIVE hazard hr against a 0-3-year-old with cumulative",
    "AUC >= 54 mg*h/L. Algebraic and deterministic (no ODE state, no drug input,",
    "no IIV, no residual error)."
  )
  reference <- paste(
    "Huang C, Du X, Xue L, Song L, Liu X, Chen Z, Wang X, Qin M, Sun X, Wang D,",
    "Yang J, Li D, Jiang H, Ding W, Qu Y, Zhu Z, Zhao L, Yuan J, Dong J, Wu D,",
    "Hu S, Miao L. Association of Busulfan Exposure with Event-Free Survival in",
    "a Chinese Pediatric Cohort: A Multicenter, Prospective Cohort Study.",
    "Drug Des Devel Ther. 2026;20:593835. doi:10.2147/DDDT.S593835.",
    "Hazard ratios from Supplementary Table S6 (multivariate column).",
    sep = " "
  )
  vignette <- "Huang_2026_busulfan"
  units <- list(
    time = "n/a (semiparametric Cox regression; the baseline hazard and hence the time scale are left unspecified; follow-up was measured from transplantation)",
    dosing = "n/a (exposure-response model; busulfan exposure enters as the AUC_BUSULFAN_COURSE covariate, not as a dose record)",
    concentration = "hr (hazard of NRM relative to a 0-3-year-old with AUC_BUSULFAN_COURSE >= 54 mg*h/L, unitless); not a drug concentration"
  )

  covariateData <- list(
    AUC_BUSULFAN_COURSE = list(
      description = "Cumulative busulfan AUC over the whole 3-day conditioning course (12 doses), computed as dose divided by the individual empirical-Bayes clearance summed over the course",
      units = "mg*h/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Bayesian individual CL from modellib('Huang_2026_busulfan') (Methods, 'Pharmacokinetics of Busulfan'). The paper's cAUC3 days; median 54.2 mg*h/L, range 26.5-119.0 in the 3-day cohort (Results). The model consumes it only through the indicator AUC_BUSULFAN_COURSE < cauc_cut.",
      source_name = "cAUC3 days"
    ),
    AGE = list(
      description = "Postnatal age at transplantation",
      units = "years",
      type = "continuous",
      reference_category = "0-3 years",
      notes = "Enters only through the Table S6 age bands 0-3 (reference), 3-6, 6-12 and > 12 years. Band edges are read as AGE < 3, 3 <= AGE < 6, 6 <= AGE <= 12 and AGE > 12; the paper does not state which band a boundary age falls in.",
      source_name = "Age (years)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 326L,
    n_studies = 1L,
    n_events = 26L,
    age_range = "0.5-17.0 years (median 8.6)",
    weight_range = "5.0-96.5 kg (median 27.0)",
    sex_female_pct = 38.65,
    race_ethnicity = "Chinese",
    disease_state = "Pediatric first allogeneic HSCT with a 3-day IV busulfan-based myeloablative conditioning regimen; 76.7% malignant disease. See modellib('Huang_2026_busulfan') for the full cohort description.",
    dose_range = "0.8-1.2 mg/kg per dose, 2 h IV infusion four times daily for 3 days",
    regions = "China (five pediatric HSCT centres)",
    notes = "26 NRM events (Table S6). cAUC3 days >= 54 mg*h/L: 166 patients (8 events); < 54: 160 (18 events). Age 0-3 years: 33 (7 events); 3-6: 73 (7); 6-12: 148 (6); > 12: 72 (6)."
  )

  ini({
    # Cut point shared with the 3-day EFS analysis (Figure 3B), chosen rather
    # than estimated with uncertainty.
    cauc_cut <- fixed(54); label("Cumulative 3-day busulfan AUC cut point separating the low-exposure group from the reference (mg*h/L)") # Discussion; Table S6

    # Log-CI midpoint (log(1.07) + log(5.78)) / 2 = 0.9110 matches log(2.48) = 0.9083.
    e_auc_busulfan_course_low_haz <- log(2.48); label("Log hazard ratio for NRM, cAUC3 days < 54 vs >= 54 mg*h/L (log scale; HR 2.48)") # Table S6 multivariate HR 2.48 (95% CI 1.07-5.78), P = 0.04

    e_age_3to6_haz <- log(0.39); label("Log hazard ratio for NRM, age 3-6 vs 0-3 years (log scale; HR 0.39)") # Table S6 multivariate HR 0.39 (95% CI 0.14-1.10), P = 0.08
    e_age_6to12_haz <- log(0.16); label("Log hazard ratio for NRM, age 6-12 vs 0-3 years (log scale; HR 0.16)") # Table S6 multivariate HR 0.16 (95% CI 0.05-0.48), P < 0.01
    e_age_over12_haz <- log(0.40); label("Log hazard ratio for NRM, age > 12 vs 0-3 years (log scale; HR 0.40)") # Table S6 multivariate HR 0.40 (95% CI 0.13-1.19), P = 0.10
  })

  model({
    # Reference group: cAUC3 days >= 54 mg*h/L and age 0-3 years (Table S6).
    caucLow <- (AUC_BUSULFAN_COURSE < cauc_cut)
    age3to6 <- (AGE >= 3) * (AGE < 6)
    age6to12 <- (AGE >= 6) * (AGE <= 12)
    ageOver12 <- (AGE > 12)

    lhr <- e_auc_busulfan_course_low_haz * caucLow +
      e_age_3to6_haz * age3to6 +
      e_age_6to12_haz * age6to12 +
      e_age_over12_haz * ageOver12
    hr <- exp(lhr)
  })
}
