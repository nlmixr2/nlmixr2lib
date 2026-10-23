Samineni_2022_venetoclax_febrile_neutropenia <- function() {
  description <- paste(
    "Binomial logistic-regression exposure-SAFETY model for grade >= 3",
    "febrile neutropenia as a function of the venetoclax steady-state",
    "nominal-dose AUC (AUCss nominal) plus the retained patient-level",
    "covariate IPI_GE4, in patients with non-Hodgkin lymphoma treated",
    "with venetoclax plus R-CHOP in CAVALLI (Samineni 2022 Table 4, final",
    "model; n = 216, 68 events/216 patients). The 95% CI of the exposure",
    "slope spans zero (not significant). There is no PK layer and no ODE:",
    "the exposure is supplied as a data column, derived in the source",
    "from the companion population PK model packaged as",
    "Samineni_2022_venetoclax. The curve drawn in Figure 5 differs from",
    "these coefficients; it is packaged separately as",
    "Samineni_2022_venetoclax_febrile_neutropenia_base. No random effects",
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
    concentration = "prob_febrile_neutropenia (probability of grade >= 3 febrile neutropenia, 0-1; also logit_febrile_neutropenia)"
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
    IPI_GE4 = list(
      description = "International Prognostic Index high-risk group at baseline: 1 = IPI 4 or 5, 0 = IPI 0 to 3.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (IPI 0 to 3)",
      notes = "Table 4 row 'IPI (high) [4, 5]'. Methods compare high (4, 5) against non-high (2, 3); the Table 4 footnote says the comparison group was IPI 0, 1 and 2. Table 1 has 26 patients with IPI 0-1, so the reference is read as every non-high patient (IPI 0 to 3). 51 of 232 CAVALLI patients (22.0%) had IPI 4-5.",
      source_name = "IPI (high) [4, 5]"
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
    # Samineni 2022 Table 4, 'grade >= 3 febrile neutropenia' block. Uncentred logit-linear
    # model: logit(p) = intercept + slope * AUC_VENETOCLAX_SS + covariate terms.
    logit_ref <- -0.81; label("Logit of the event probability at zero venetoclax exposure with every covariate at zero (unitless logit)") # Table 4 Intercept -0.81 (SE 0.42, 95% CI -1.63 to 0.01), p = 0.052
    e_auc_logit <- -0.46; label("Log-odds of the event per 1 ug*day/mL of venetoclax AUCss nominal (unitless logit)") # Table 4 Slope of AUCss nominal -0.46 (SE 0.66, 95% CI -1.75 to 0.84), p = 0.0608
    e_ipi_ge4_logit <- 1.20; label("Log-odds of the event for the high-risk IPI group versus the rest (unitless logit)") # Table 4 IPI (high) [4, 5] 1.20 (SE 0.33, 95% CI 0.55 to 1.84), p = 0.00026

    # Logistic regression (glm, logit link): Bernoulli likelihood, no random
    # effects and no residual error. The tiny additive residual below exists
    # only so rxode2 has an error model to attach to the probability; it is
    # not a published quantity.
    addSd_prob_febrile_neutropenia <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_febrile_neutropenia <- logit_ref + e_auc_logit * AUC_VENETOCLAX_SS +
      e_ipi_ge4_logit * IPI_GE4
    prob_febrile_neutropenia <- expit(logit_febrile_neutropenia)

    prob_febrile_neutropenia ~ add(addSd_prob_febrile_neutropenia)
  })
}
