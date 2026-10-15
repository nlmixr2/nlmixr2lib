PonceBobadilla_2023_upadacitinib_induction_clinresp <- function() {
  description <- paste0(
    "Landmark logistic exposure-response model for the probability of ",
    "clinical response per Adapted Mayo score at the end of 8-week ",
    "induction in adults with moderately to severely active ulcerative ",
    "colitis (1,234 patients from phase 2b Study 1 and phase 3 Studies 2 ",
    "and 4; placebo or upadacitinib 7.5-45 mg once daily). Clinical ",
    "response per Adapted Mayo score is decrease from baseline in the ",
    "Adapted Mayo score of >= 2 points and >= 30%, plus a decrease in the ",
    "rectal-bleeding subscore of >= 1 or an absolute rectal-bleeding ",
    "subscore <= 1 (OSM Table 2). Drug effect is logarithmic in exposure, ",
    "beta * log(1 + CAV), with CAV the steady-state average plasma ",
    "concentration (ng/mL) from ",
    "modellib('PonceBobadilla_2023_upadacitinib') and 0 for placebo. ",
    "Pre-specified stratification covariates on the logit: Bio-IR status, ",
    "baseline corticosteroid use, baseline Adapted Mayo score > 7, and ",
    "study (Study 2 and Study 4 versus Study 1). No between-subject ",
    "variability and no residual error in the source (Bernoulli ",
    "likelihood). Model-predicted response by arm in Figure 5 of the ",
    "paper: 23, 48, 57, 65 and 70% for placebo, 7.5, 15, 30 and 45 mg. ",
    "One of nine exposure-response models in the ",
    "PonceBobadilla_2023_upadacitinib_* family."
  )
  reference <- paste(
    "Ponce-Bobadilla AV, Stodtmann S, Eckert D, Zhou W, Liu W, Mohamed M-EF.",
    "Upadacitinib population pharmacokinetics and exposure-response",
    "relationships in ulcerative colitis patients. Clin Pharmacokinet.",
    "2023;62(1):101-112. doi:10.1007/s40262-022-01191-6.",
    "Coefficients from Supplementary Table 7 of the Online Supplementary",
    "Material (OSM), block 'Clinical Response per Adapted Mayo Score [NRI]'; model form from the OSM Methods",
    "'Exposure-response analysis'. The exposure metric is the steady-state",
    "average concentration from modellib('PonceBobadilla_2023_upadacitinib').",
    sep = " "
  )
  vignette <- "PonceBobadilla_2023_upadacitinib"
  units <- list(
    time = "n/a (static landmark logistic regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "prob_clinresp (probability of clinical response per Adapted Mayo score at Week 8, 0-1; not a drug concentration)"
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
        "-1.04 for Bio-IR versus non-Bio-IR (OSM Table 7). Non-Bio-IR ",
        "includes biologic-naive patients and patients exposed to a biologic ",
        "without intolerance or inadequate response."
      ),
      source_name = "Bio-IR"
    ),
    CONMED_STEROID = list(
      description = "Systemic corticosteroid use at induction baseline (1 = yes, 0 = no).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no corticosteroid use).",
      notes = paste0(
        "Time-fixed pre-specified stratification variable. Log-odds shift ",
        "0.176 (OSM Table 7, row 'Baseline Corticosteroid Use (yes versus ",
        "no)')."
      ),
      source_name = "Baseline Corticosteroid Use"
    ),
    SCORE_AMAYO = list(
      description = "Baseline Adapted Mayo score for ulcerative colitis (stool-frequency + rectal-bleeding + endoscopic subscores, 0-9).",
      units = "(score, 0-9)",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Enters only through the stratification indicator (SCORE_AMAYO > 7), ",
        "computed in model(); reference category is a score <= 7. Log-odds ",
        "shift 0.387 for a baseline score > 7 (OSM Table 7, row 'Baseline ",
        "Adapted Mayo Score (>7 versus <=7)')."
      ),
      source_name = "Baseline Adapted Mayo Score"
    ),
    STUDY_UPA_UC_S2 = list(
      description = "Ulcerative colitis induction Study 2 indicator (1 = Study 2, the phase 3 induction substudy of NCT02819635; 0 = other).",
      units = "(binary)",
      type = "binary",
      reference_category = "0. Study 1 (the phase 2b dose-ranging induction study) is the reference for both study indicators.",
      notes = "Log-odds shift 0.539 (OSM Table 7, row 'Study 2'). Set to 0 to predict a Study 1 patient.",
      source_name = "Study 2"
    ),
    STUDY_UPA_UC_S4 = list(
      description = "Ulcerative colitis induction Study 4 indicator (1 = Study 4, the replicate phase 3 induction study NCT03653026; 0 = other).",
      units = "(binary)",
      type = "binary",
      reference_category = "0. Study 1 (the phase 2b dose-ranging induction study) is the reference for both study indicators.",
      notes = "Log-odds shift 0.581 (OSM Table 7, row 'Study 4'). Set to 0 to predict a Study 1 patient.",
      source_name = "Study 4"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1234L,
    n_studies = 3L,
    age_median = "41 years (OSM Table 8 scaling median)",
    weight_median = "70.35 kg (OSM Table 8 scaling median)",
    sex_female_pct = 38,
    disease_state = "Moderately to severely active ulcerative colitis; 58% biologic-intolerant or inadequate responders in the PK subset (OSM Table 3).",
    dose_range = "Placebo or upadacitinib 7.5, 15, 30 or 45 mg once daily (extended-release tablet) for 8 weeks; 7.5-30 mg only in Study 1.",
    regions = "Multinational; Study 1 (phase 2b, 248 main-group patients), Study 2 (phase 3, 472) and Study 4 (phase 3, 514).",
    notes = paste0(
      "Week-8 landmark logistic regression, one record per patient with ",
      "non-responder imputation, fitted with glm in R 4.0.2. The ",
      "sex percentage is from the ulcerative colitis PK subset (OSM Table 3); ",
      "the exposure-response set is not tabulated separately. Other medians ",
      "(OSM Table 8): CRP 4.78 mg/L, fecal calprotectin 1737 mg/kg, albumin ",
      "4.3 g/dL, disease duration 5.85 years."
    )
  )

  ini({
    # OSM Table 7, block 'Clinical Response per Adapted Mayo Score [NRI]' under 'Induction Endpoints'. Every
    # coefficient is on the log-odds scale and was estimated; the 95% CIs and
    # p-values are as printed.
    logite0 <- -1.41
    label("Logit of the probability of clinical response per Adapted Mayo score at Week 8 for a reference placebo patient (unitless)") # OSM Table 7 Intercept = -1.41 (95% CI -1.83, -1)
    e_cav_clinresp <- 0.598
    label("Log-odds shift per unit of log(1 + CAV), CAV in ng/mL (unitless)") # OSM Table 7 log(1+CAVG) = 0.598 (95% CI 0.517, 0.678), p < 0.05
    e_biologic_ir_clinresp <- -1.04
    label("Log-odds shift for biologic-intolerant or inadequate responders (unitless)") # OSM Table 7 Bio-IR (yes versus no) = -1.04 (95% CI -1.32, -0.764), p < 0.05
    e_conmed_steroid_clinresp <- 0.176
    label("Log-odds shift for baseline corticosteroid use (unitless)") # OSM Table 7 corticosteroid use (yes versus no) = 0.176 (95% CI -0.0918, 0.444), p not significant
    e_score_amayo_clinresp <- 0.387
    label("Log-odds shift for a baseline Adapted Mayo score above 7 (unitless)") # OSM Table 7 Baseline Adapted Mayo Score (>7 versus <=7) = 0.387 (95% CI 0.12, 0.653), p < 0.05
    e_study_upa_uc_s2_clinresp <- 0.539
    label("Log-odds shift for Study 2 versus Study 1 (unitless)") # OSM Table 7 Study 2 = 0.539 (95% CI 0.189, 0.889), p < 0.05
    e_study_upa_uc_s4_clinresp <- 0.581
    label("Log-odds shift for Study 4 versus Study 1 (unitless)") # OSM Table 7 Study 4 = 0.581 (95% CI 0.235, 0.927), p < 0.05

    # The source maximises a Bernoulli likelihood on a 0/1 outcome; this
    # placeholder residual only lets nlmixr2 accept the model for simulation.
    addSd_prob_clinresp <- fixed(0.001)
    label("Placeholder additive residual SD on prob_clinresp; not from the source (unitless)")
  })

  model({
    # OSM Methods: logit P = alpha + f(Cavg) + sum(beta_j * x_j).
    amayo_gt7 <- 0
    if (SCORE_AMAYO > 7) amayo_gt7 <- 1
    logit_clinresp <- logite0 +
      e_cav_clinresp * log(1 + CAV) +
      e_biologic_ir_clinresp * BIOLOGIC_IR +
      e_conmed_steroid_clinresp * CONMED_STEROID +
      e_score_amayo_clinresp * amayo_gt7 +
      e_study_upa_uc_s2_clinresp * STUDY_UPA_UC_S2 +
      e_study_upa_uc_s4_clinresp * STUDY_UPA_UC_S4
    prob_clinresp <- expit(logit_clinresp)
    prob_clinresp ~ add(addSd_prob_clinresp)
  })
}
