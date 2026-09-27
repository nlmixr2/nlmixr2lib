Zhou_2020_paclitaxel_liposomal_neutropenia <- function() {
  description <- paste(
    "Landmark logistic-regression exposure-safety model for neutropenia of",
    "CTCAE v4.0 grade 2 or worse (the source's 'grade > 1') as a linear",
    "function of the individual total-paclitaxel AUC, in adults with",
    "squamous non-small cell lung cancer given paclitaxel liposome",
    "(Lipusu) 175 mg/m^2 as a 3-h infusion followed by platinum",
    "chemotherapy (Zhou 2020, n = 45). There is no PK layer and no ODE:",
    "the exposure metric is supplied as the AUC_PTX data column, which the",
    "source computed from each patient's dose and post hoc clearance from",
    "the companion population PK model Zhou_2020_paclitaxel_liposomal.",
    "The exposure slope was significant at p = 0.0469."
  )
  reference <- paste(
    "Zhou H, Yan J, Chen W, Yang J, Liu M, Zhang Y, Shen X, Ma Y, Hu X,",
    "Wang Y, Du K, Li G. Population Pharmacokinetics and Exposure-Safety",
    "Relationship of Paclitaxel Liposome in Patients With Non-small Cell",
    "Lung Cancer. Front Oncol. 2020;10:1731 (issue dated 5 February 2021).",
    "doi:10.3389/fonc.2020.01731.",
    "Model form in Results 'Exposure-Safety Analysis'; estimates in",
    "Table 3; fitted curve in Figure 5.",
    "Individual exposures derive from the companion population PK model;",
    "see modellib('Zhou_2020_paclitaxel_liposomal').",
    sep = " "
  )
  vignette <- "Zhou_2020_paclitaxel_liposomal"
  units <- list(
    time = "n/a (static landmark exposure-response regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the AUC_PTX covariate column)",
    concentration = "prob_neutropenia_grade2 (probability of CTCAE grade 2 or worse neutropenia, 0-1; also logit_neutropenia_grade2)"
  )

  covariateData <- list(
    AUC_PTX = list(
      description = "Individual area under the total (liposome-encapsulated plus released) plasma paclitaxel concentration-time curve for one treatment cycle, computed from the administered dose and the individual post hoc clearance. Supplied as data: this model has no PK layer.",
      units = "mg*h/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Methods 'Exposure-Safety Analysis': 'The individual post hoc PK",
        "parameters and dosage were used to calculate the AUC.' For the",
        "linear companion PK model this is Dose / CL_i. The units are the",
        "Figure 5 x-axis label (mg*hr/L), and the observed exposures there",
        "span about 7.7 to 17 mg*h/L, consistent with a 210-300 mg dose over",
        "a typical clearance of 21.55 L/h (240 mg / 21.55 L/h = 11.1",
        "mg*h/L). Enters the logit linearly and uncentred.",
        "The source does not say which cycle's AUC was paired with the",
        "neutropenia outcome; see the vignette."
      ),
      source_name = "AUC"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 45L,
    n_studies = 1L,
    age_range = "36-75 years",
    age_median = "59 years",
    weight_range = "45-100 kg",
    weight_median = "71.5 kg",
    sex_female_pct = 8.9,
    race_ethnicity = "Chinese (single centre in Beijing); race not tabulated",
    disease_state = "Squamous non-small cell lung cancer",
    dose_range = "Paclitaxel liposome 175 mg/m^2 (210-300 mg after vial rounding) as a 3-h IV infusion on day 1 of a 3-week cycle, then cisplatin 75 mg/m^2 or carboplatin AUC 4-5 on day 2",
    regions = "China",
    notes = "Same cohort as the companion population PK model (Table 1). Neutropenia of any grade was collected from electronic medical records and graded by CTCAE v4.0; the analysis was not stratified by grade because of sparse data per grade, and dichotomised at grade > 1. Figure 5 plots one event indicator per patient."
  )

  ini({
    # logit[P(NE > 1)] = a + b * AUC      (Results 'Exposure-Safety Analysis')
    # Table 3 prints both coefficients; the fitted curve in Figure 5
    # reproduces them (P = 0.16 at 5 mg*h/L, 0.56 at 10, 0.98 at 20).
    logit_ref <- -3.5008; label("Logit of the grade 2 or worse neutropenia probability at zero paclitaxel AUC (unitless logit)") # Table 3, Neutropenia grade > 1, a = -3.5008 (SE 2.1942), p = 0.1106
    e_auc_logit <- 0.372; label("Change in the log-odds of grade 2 or worse neutropenia per 1 mg*h/L increase in total paclitaxel AUC (1/(mg*h/L))") # Table 3, Neutropenia grade > 1, b = 0.372 (SE 0.1872), p = 0.0469

    # The source likelihood is Bernoulli, so it estimates neither a
    # random effect nor a residual variance. The small additive term
    # below exists only so rxode2 has an error model to attach to the
    # typical-value probability; it is NOT a published quantity.
    addSd_prob_neutropenia_grade2 <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_neutropenia_grade2 <- logit_ref + e_auc_logit * AUC_PTX
    prob_neutropenia_grade2 <- expit(logit_neutropenia_grade2)

    prob_neutropenia_grade2 ~ add(addSd_prob_neutropenia_grade2)
  })
}
