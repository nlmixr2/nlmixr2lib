PonceBobadilla_2023_upadacitinib_maintenance_clinrem <- function() {
  description <- paste0(
    "Landmark logistic exposure-response model for the probability of ",
    "clinical remission per Adapted Mayo score at Week 52 of maintenance ",
    "in adults with ulcerative colitis who responded to 8-week ",
    "upadacitinib 45 mg induction (449 patients from phase 3 Study 3; ",
    "placebo or upadacitinib 15 or 30 mg once daily). Clinical remission ",
    "per Adapted Mayo score is stool-frequency subscore <= 1 and not ",
    "greater than baseline, rectal-bleeding subscore 0 and endoscopic ",
    "subscore <= 1 (OSM Table 2). Drug effect is logarithmic in exposure, ",
    "beta * log(1 + CAV), with CAV the steady-state average plasma ",
    "concentration (ng/mL) from ",
    "modellib('PonceBobadilla_2023_upadacitinib') and 0 for placebo. ",
    "Pre-specified stratification covariates on the logit: Bio-IR status, ",
    "corticosteroid use at maintenance Week 0 and clinical remission at ",
    "maintenance Week 0. No between-subject variability and no residual ",
    "error in the source (Bernoulli likelihood). Model-predicted response ",
    "by arm in Figure 6 of the paper: 13, 40 and 48% for placebo, 15 and ",
    "30 mg. One of nine exposure-response models in the ",
    "PonceBobadilla_2023_upadacitinib_* family."
  )
  reference <- paste(
    "Ponce-Bobadilla AV, Stodtmann S, Eckert D, Zhou W, Liu W, Mohamed M-EF.",
    "Upadacitinib population pharmacokinetics and exposure-response",
    "relationships in ulcerative colitis patients. Clin Pharmacokinet.",
    "2023;62(1):101-112. doi:10.1007/s40262-022-01191-6.",
    "Coefficients from Supplementary Table 7 of the Online Supplementary",
    "Material (OSM), block 'Clinical Remission per Adapted Mayo Score [NRI]'; model form from the OSM Methods",
    "'Exposure-response analysis'. The exposure metric is the steady-state",
    "average concentration from modellib('PonceBobadilla_2023_upadacitinib').",
    sep = " "
  )
  vignette <- "PonceBobadilla_2023_upadacitinib"
  units <- list(
    time = "n/a (static landmark logistic regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "prob_clinrem (probability of clinical remission per Adapted Mayo score at Week 52, 0-1; not a drug concentration)"
  )

  covariateData <- list(
    CAV = list(
      description = "Model-estimated steady-state average upadacitinib plasma concentration over the once-daily dosing interval.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Enters the logit as log(1 + CAV), natural log, uncentered, so the ",
        "intercept is the logit at zero exposure (the placebo arm). The OSM ",
        "Methods compute it per subject as Dose x F_ER / (CL/F x 24 h) from ",
        "the post hoc estimates of ",
        "modellib('PonceBobadilla_2023_upadacitinib'), and set it to 0 for ",
        "every placebo patient. The +1 inside the logarithm makes the ",
        "exposure term vanish at zero, so the unit matters: CAV must be ",
        "supplied in ng/mL. Model-estimated steady-state medians in ",
        "ulcerative colitis (OSM Table 5): 6.69, 14.3, 27.2 and 41.4 ng/mL at ",
        "7.5, 15, 30 and 45 mg once daily."
      ),
      source_name = "CAVG"
    ),
    BIOLOGIC_IR = list(
      description = "Biologic therapy-intolerant or inadequate responder (Bio-IR) indicator (1 = Bio-IR, 0 = not).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Bio-IR, the paper's 'no' category).",
      notes = paste0(
        "Pre-specified stratification variable (OSM Methods). Log-odds shift ",
        "-0.156 for Bio-IR versus non-Bio-IR (OSM Table 7). Non-Bio-IR ",
        "includes biologic-naive patients and patients exposed to a biologic ",
        "without intolerance or inadequate response."
      ),
      source_name = "Bio-IR"
    ),
    CONMED_STEROID = list(
      description = "Systemic corticosteroid use at Week 0 of the maintenance study (the end of induction) (1 = yes, 0 = no).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no corticosteroid use).",
      notes = paste0(
        "Time-fixed pre-specified stratification variable. Log-odds shift ",
        "-0.121 (OSM Table 7, row 'Corticosteroid Use at Week 0 of Study 3 ",
        "(yes versus no)')."
      ),
      source_name = "Corticosteroid Use at Week 0 of Study 3"
    ),
    CLINICAL_REMISSION_BL = list(
      description = "Clinical remission per Adapted Mayo score at Week 0 of the maintenance study (1 = in remission on entering maintenance, 0 = clinical responder not in remission).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not in clinical remission at maintenance Week 0).",
      notes = paste0(
        "Pre-specified stratification variable; remission uses the OSM Table ",
        "2 definition. Log-odds shift 1.02 (OSM Table 7, row 'Clinical ",
        "Remission per Adapted Mayo Score at Week 0 of Study 3 (yes versus ",
        "no)'). Every maintenance patient was an 8-week clinical responder to ",
        "upadacitinib 45 mg once daily, so this indicator separates ",
        "responders that were also remitters."
      ),
      source_name = "Clinical Remission per Adapted Mayo Score at Week 0 of Study 3"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 449L,
    n_studies = 1L,
    age_median = "41 years (OSM Table 8 scaling median)",
    weight_median = "69.85 kg (OSM Table 8 scaling median)",
    sex_female_pct = 38,
    disease_state = "Moderately to severely active ulcerative colitis, clinical responders after 8 weeks of upadacitinib 45 mg once daily induction; 53% biologic-intolerant or inadequate responders in the PK subset (OSM Table 4).",
    dose_range = "Placebo or upadacitinib 15 or 30 mg once daily (extended-release tablet) for 52 weeks.",
    regions = "Multinational; maintenance Study 3 (phase 3, NCT02819635).",
    notes = paste0(
      "Week-52 landmark logistic regression on the first approximately 450 ",
      "45 mg induction responders randomized in Study 3, non-responder ",
      "imputation, glm in R 4.0.2. Placebo patients carry CAV = 0 although ",
      "they received 45 mg during induction. The sex percentage is from the ",
      "maintenance PK subset (OSM Table 4). Other medians (OSM Table 8): CRP ",
      "4.05 mg/L, fecal calprotectin 1658.5 mg/kg, albumin 4.3 g/dL, disease ",
      "duration 6.21 years."
    )
  )

  ini({
    # OSM Table 7, block 'Clinical Remission per Adapted Mayo Score [NRI]' under 'Maintenance Endpoints'. Every
    # coefficient is on the log-odds scale and was estimated; the 95% CIs and
    # p-values are as printed.
    logite0 <- -2.22
    label("Logit of the probability of clinical remission per Adapted Mayo score at Week 52 for a reference placebo patient (unitless)") # OSM Table 7 Intercept = -2.22 (95% CI -2.80, -1.64)
    e_cav_clinrem <- 0.57
    label("Log-odds shift per unit of log(1 + CAV), CAV in ng/mL (unitless)") # OSM Table 7 log(1+CAVG) = 0.57 (95% CI 0.404, 0.737), p < 0.05
    e_biologic_ir_clinrem <- -0.156
    label("Log-odds shift for biologic-intolerant or inadequate responders (unitless)") # OSM Table 7 Bio-IR (yes versus no) = -0.156 (95% CI -0.598, 0.286), p not significant
    e_conmed_steroid_clinrem <- -0.121
    label("Log-odds shift for corticosteroid use at maintenance Week 0 (unitless)") # OSM Table 7 corticosteroid use (yes versus no) = -0.121 (95% CI -0.587, 0.344), p not significant
    e_clinical_remission_bl_clinrem <- 1.02
    label("Log-odds shift for clinical remission at maintenance Week 0 (unitless)") # OSM Table 7 Clinical Remission at Week 0 of Study 3 (yes versus no) = 1.02 (95% CI 0.570, 1.46), p < 0.05

    # The source maximises a Bernoulli likelihood on a 0/1 outcome; this
    # placeholder residual only lets nlmixr2 accept the model for simulation.
    addSd_prob_clinrem <- fixed(0.001)
    label("Placeholder additive residual SD on prob_clinrem; not from the source (unitless)")
  })

  model({
    # OSM Methods: logit P = alpha + f(Cavg) + sum(beta_j * x_j).
    logit_clinrem <- logite0 +
      e_cav_clinrem * log(1 + CAV) +
      e_biologic_ir_clinrem * BIOLOGIC_IR +
      e_conmed_steroid_clinrem * CONMED_STEROID +
      e_clinical_remission_bl_clinrem * CLINICAL_REMISSION_BL
    prob_clinrem <- expit(logit_clinrem)
    prob_clinrem ~ add(addSd_prob_clinrem)
  })
}
