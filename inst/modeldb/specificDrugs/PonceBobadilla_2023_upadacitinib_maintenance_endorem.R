PonceBobadilla_2023_upadacitinib_maintenance_endorem <- function() {
  description <- paste0(
    "Landmark logistic exposure-response model for the probability of ",
    "endoscopic remission at Week 52 of maintenance in adults with ",
    "ulcerative colitis who responded to 8-week upadacitinib 45 mg ",
    "induction (449 patients from phase 3 Study 3; placebo or ",
    "upadacitinib 15 or 30 mg once daily). Endoscopic remission is Mayo ",
    "endoscopic subscore of 0 (OSM Table 2). The selected model is a ",
    "TREATMENT-EFFECT model: both upadacitinib doses carry one shared ",
    "log-odds shift through ON_TREATMENT and no exposure term enters; the ",
    "paper found no exposure-response relationship beyond the treatment ",
    "effect for this endpoint. Pre-specified stratification covariates on ",
    "the logit: Bio-IR status, corticosteroid use at maintenance Week 0 ",
    "and clinical remission at maintenance Week 0. No between-subject ",
    "variability and no residual error in the source (Bernoulli ",
    "likelihood). Not shown in the paper's Figure 6, which plots only the ",
    "endpoints with an exposure-response relationship. One of nine ",
    "exposure-response models in the PonceBobadilla_2023_upadacitinib_* ",
    "family."
  )
  reference <- paste(
    "Ponce-Bobadilla AV, Stodtmann S, Eckert D, Zhou W, Liu W, Mohamed M-EF.",
    "Upadacitinib population pharmacokinetics and exposure-response",
    "relationships in ulcerative colitis patients. Clin Pharmacokinet.",
    "2023;62(1):101-112. doi:10.1007/s40262-022-01191-6.",
    "Coefficients from Supplementary Table 7 of the Online Supplementary",
    "Material (OSM), block 'Endoscopic Remission [NRI]'; model form from the OSM Methods",
    "'Exposure-response analysis'. The exposure metric is the steady-state",
    "average concentration from modellib('PonceBobadilla_2023_upadacitinib').",
    sep = " "
  )
  vignette <- "PonceBobadilla_2023_upadacitinib"
  units <- list(
    time = "n/a (static landmark logistic regression; no time dimension)",
    dosing = "n/a (no dose events; treatment enters as the ON_TREATMENT covariate)",
    concentration = "prob_endorem (probability of endoscopic remission at Week 52, 0-1; not a drug concentration)"
  )

  covariateData <- list(
    BIOLOGIC_IR = list(
      description = "Biologic therapy-intolerant or inadequate responder (Bio-IR) indicator (1 = Bio-IR, 0 = not).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Bio-IR, the paper's 'no' category).",
      notes = paste0(
        "Pre-specified stratification variable (OSM Methods). Log-odds shift ",
        "-0.277 for Bio-IR versus non-Bio-IR (OSM Table 7). Non-Bio-IR ",
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
        "-0.444 (OSM Table 7, row 'Corticosteroid Use at Week 0 of Study 3 ",
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
        "2 definition. Log-odds shift 1.04 (OSM Table 7, row 'Clinical ",
        "Remission per Adapted Mayo Score at Week 0 of Study 3 (yes versus ",
        "no)'). Every maintenance patient was an 8-week clinical responder to ",
        "upadacitinib 45 mg once daily, so this indicator separates ",
        "responders that were also remitters."
      ),
      source_name = "Clinical Remission per Adapted Mayo Score at Week 0 of Study 3"
    ),
    ON_TREATMENT = list(
      description = "Maintenance treatment-arm indicator (1 = upadacitinib 15 or 30 mg once daily, 0 = placebo).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (placebo; the intercept is the placebo log-odds).",
      notes = paste0(
        "Log-odds shift 1.75 for active treatment (OSM Table 7, row 'Active ",
        "treatment'). A treatment-effect model was selected for this endpoint ",
        "over the linear, logarithmic and Emax exposure models (main text ",
        "Section 3.3), so both doses carry the same effect and no exposure ",
        "term enters."
      ),
      source_name = "Active treatment"
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
    # OSM Table 7, block 'Endoscopic Remission [NRI]' under 'Maintenance Endpoints'. Every
    # coefficient is on the log-odds scale and was estimated; the 95% CIs and
    # p-values are as printed.
    logite0 <- -3.1
    label("Logit of the probability of endoscopic remission at Week 52 for a reference placebo patient (unitless)") # OSM Table 7 Intercept = -3.1 (95% CI -3.94, -2.26)
    e_on_treatment_endorem <- 1.75
    label("Log-odds shift for upadacitinib 15 or 30 mg versus placebo (unitless)") # OSM Table 7 Active treatment = 1.75 (95% CI 0.977, 2.52), p < 0.05
    e_biologic_ir_endorem <- -0.277
    label("Log-odds shift for biologic-intolerant or inadequate responders (unitless)") # OSM Table 7 Bio-IR (yes versus no) = -0.277 (95% CI -0.818, 0.263), p not significant
    e_conmed_steroid_endorem <- -0.444
    label("Log-odds shift for corticosteroid use at maintenance Week 0 (unitless)") # OSM Table 7 corticosteroid use (yes versus no) = -0.444 (95% CI -1.04, 0.151), p not significant
    e_clinical_remission_bl_endorem <- 1.04
    label("Log-odds shift for clinical remission at maintenance Week 0 (unitless)") # OSM Table 7 Clinical Remission at Week 0 of Study 3 (yes versus no) = 1.04 (95% CI 0.515, 1.56), p < 0.05

    # The source maximises a Bernoulli likelihood on a 0/1 outcome; this
    # placeholder residual only lets nlmixr2 accept the model for simulation.
    addSd_prob_endorem <- fixed(0.001)
    label("Placeholder additive residual SD on prob_endorem; not from the source (unitless)")
  })

  model({
    # OSM Methods: logit P = alpha + f(Cavg) + sum(beta_j * x_j).
    logit_endorem <- logite0 +
      e_on_treatment_endorem * ON_TREATMENT +
      e_biologic_ir_endorem * BIOLOGIC_IR +
      e_conmed_steroid_endorem * CONMED_STEROID +
      e_clinical_remission_bl_endorem * CLINICAL_REMISSION_BL
    prob_endorem <- expit(logit_endorem)
    prob_endorem ~ add(addSd_prob_endorem)
  })
}
