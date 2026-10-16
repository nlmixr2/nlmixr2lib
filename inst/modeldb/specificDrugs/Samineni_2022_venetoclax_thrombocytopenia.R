Samineni_2022_venetoclax_thrombocytopenia <- function() {
  description <- paste(
    "Binomial logistic-regression exposure-SAFETY model for grade >= 3",
    "thrombocytopenia as a function of the venetoclax steady-state",
    "nominal-dose AUC (AUCss nominal) plus the retained patient-level",
    "covariates AGE and ECOG_GE2, in patients with non-Hodgkin lymphoma",
    "treated with venetoclax plus R-CHOP in CAVALLI (Samineni 2022 Table",
    "4, final model; n = 216, 46 events/216 patients). The 95% CI of the",
    "exposure slope spans zero (not significant). There is no PK layer",
    "and no ODE: the exposure is supplied as a data column, derived in",
    "the source from the companion population PK model packaged as",
    "Samineni_2022_venetoclax. The curve drawn in Figure 5 differs from",
    "these coefficients; it is packaged separately as",
    "Samineni_2022_venetoclax_thrombocytopenia_base. No random effects",
    "and no residual error (Bernoulli likelihood).",
    sep = " "
  )
  reference <- paste(
    "Samineni D, Huang W, Gibiansky L, Ding H, Zhang R, Li C, Sinha A,",
    "Rajwanshi R, Humphrey K, Bazeos A, Salem AH, Miles D.",
    "Population Pharmacokinetics and Exposure-Response Analyses for",
    "Venetoclax in Combination with R-CHOP in Relapsed/Refractory and",
    "Previously Untreated Patients with Diffuse Large B Cell Lymphoma.",
    "Adv Ther. 2022;39(1):598-618. doi:10.1007/s12325-021-01919-z.",
    "Logistic-regression coefficients from Table 4.",
    sep = " "
  )
  vignette <- "Samineni_2022_venetoclax"
  units <- list(
    time = "n/a (static landmark exposure-response regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the AUC_VENETOCLAX_SS covariate column)",
    concentration = "prob_thrombocytopenia_grade3 (probability of grade >= 3 thrombocytopenia, 0-1; also logit_thrombocytopenia_grade3)"
  )

  covariateData <- list(
    AUC_VENETOCLAX_SS = list(
      description = "Individual venetoclax steady-state nominal-dose AUC over one day, AUCss nominal = Dnom * F1 / (CL/F) (Samineni 2022 Eq. 1), from empirical Bayes post hoc PK parameters, the protocol-assigned dose and day-1 covariates.",
      units = "ug*day/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Supplied as data: this model has no PK layer. The source derived it",
        "from the companion population PK model, packaged as",
        "Samineni_2022_venetoclax. Enters the logit uncentred, per 1",
        "ug*day/mL. The source exposure axis (Figures 4 and 5) runs from",
        "about 0.20 to 1.83 ug*day/mL with exposure quartile bounds near",
        "0.43, 0.55 and 0.68. Computing Eq. 1 from the printed Table 2 values",
        "gives exposures about 1.8-fold higher than that axis; see the",
        "vignette before chaining the PK model into this one.",
        sep = " "
      ),
      source_name = "AUCss nominal"
    ),
    AGE = list(
      description = "Age at baseline.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters the logit uncentred, per year (Table 4 row 'Age'). CAVALLI ages 18-85 years, median 64 (Table 1).",
      source_name = "Age"
    ),
    ECOG_GE2 = list(
      description = "ECOG performance status at baseline >= 2: 1 = yes, 0 = ECOG 0 or 1.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ECOG 0 or 1)",
      notes = "Table 4 row 'ECOG (2+)', footnote b. 35 of 232 CAVALLI patients (15.1%) had ECOG 2.",
      source_name = "ECOG (2+)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 216L,
    n_studies = 1L,
    n_observations = "216 binary records, one per patient",
    age_range = "18-85 years (median 64) in the CAVALLI PK population (Table 1)",
    weight_range = "46.0-182 kg (median 77.0) in the CAVALLI PK population (Table 1)",
    sex_female_pct = 44.4,
    race_ethnicity = c(
      White = 71.1,
      Asian = 2.2,
      Black = 1.7,
      `Native Hawaiian or Pacific Islander` = 1.3,
      Unknown = 23.7
    ),
    disease_state = "Relapsed/refractory and previously untreated non-Hodgkin lymphoma, mostly DLBCL, phase 1b and phase 2 of CAVALLI (GO27878)",
    dose_range = "Venetoclax 400 mg (n = 3), 600 mg (n = 8) or 800 mg (n = 205) once daily for 10 days of each 21-day cycle, with R-CHOP",
    regions = "North America, Western and Eastern Europe and other regions (CAVALLI)",
    notes = "Exposure-safety analysis population. Demographic percentages are those of the 232-patient CAVALLI PK population in Table 1; the paper does not tabulate the 216-patient subset separately. Grade >= 3 adverse events were counted through 30 days after the last venetoclax or CHOP dose or 90 days after the last rituximab dose, whichever was later; serious adverse events up to 13 July 2018."
  )

  ini({
    # Samineni 2022 Table 4, 'grade >= 3 thrombocytopenia' block. Uncentred logit-linear
    # model: logit(p) = intercept + slope * AUC_VENETOCLAX_SS + covariate terms.
    logit_ref <- -8.70; label("Logit of the event probability at zero venetoclax exposure with every covariate at zero (unitless logit)") # Table 4 Intercept -8.70 (SE 1.59, 95% CI -11.82 to -5.57), p = 0.000000049
    e_auc_logit <- 0.59; label("Log-odds of the event per 1 ug*day/mL of venetoclax AUCss nominal (unitless logit)") # Table 4 Slope of AUCss nominal 0.59 (SE 0.75, 95% CI -0.88 to 2.05), p = 0.521
    e_age_logit <- 0.11; label("Log-odds of the event per year of age (unitless logit)") # Table 4 Age 0.11 (SE 0.02, 95% CI 0.06 to 0.15), p = 0.0000013
    e_ecog_ge2_logit <- 1.65; label("Log-odds of the event for ECOG performance status 2 or more versus 0 or 1 (unitless logit)") # Table 4 ECOG (2+) 1.65 (SE 0.45, 95% CI 0.75 to 2.54), p = 0.00029

    # Logistic regression (glm, logit link): Bernoulli likelihood, no random
    # effects and no residual error. The tiny additive residual below exists
    # only so rxode2 has an error model to attach to the probability; it is
    # not a published quantity.
    addSd_prob_thrombocytopenia_grade3 <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_thrombocytopenia_grade3 <- logit_ref + e_auc_logit * AUC_VENETOCLAX_SS +
      e_age_logit * AGE +
      e_ecog_ge2_logit * ECOG_GE2
    prob_thrombocytopenia_grade3 <- expit(logit_thrombocytopenia_grade3)

    prob_thrombocytopenia_grade3 ~ add(addSd_prob_thrombocytopenia_grade3)
  })
}
