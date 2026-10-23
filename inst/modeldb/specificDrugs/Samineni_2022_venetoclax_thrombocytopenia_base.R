Samineni_2022_venetoclax_thrombocytopenia_base <- function() {
  description <- paste(
    "Binomial logistic-regression exposure-SAFETY curve for grade >= 3",
    "thrombocytopenia as a function of the venetoclax steady-state",
    "nominal-dose AUC alone, as drawn in Samineni 2022 Figure 5d (n =",
    "216, 46 events/216 patients). The paper prints no coefficients for",
    "this curve; the intercept and slope were digitised from the figure",
    "by the maintainers. It does not coincide with the Table 4",
    "coefficients (packaged as",
    "Samineni_2022_venetoclax_thrombocytopenia), and is read as the",
    "exposure-only base fit whose p-value the figure and Table 4 both",
    "print. No PK layer: the exposure is supplied as a data column",
    "derived from Samineni_2022_venetoclax. No random effects and no",
    "residual error (Bernoulli likelihood).",
    sep = " "
  )
  reference <- paste(
    "Samineni D, Huang W, Gibiansky L, Ding H, Zhang R, Li C, Sinha A,",
    "Rajwanshi R, Humphrey K, Bazeos A, Salem AH, Miles D.",
    "Population Pharmacokinetics and Exposure-Response Analyses for",
    "Venetoclax in Combination with R-CHOP in Relapsed/Refractory and",
    "Previously Untreated Patients with Diffuse Large B Cell Lymphoma.",
    "Adv Ther. 2022;39(1):598-618. doi:10.1007/s12325-021-01919-z.",
    "Fitted curve of Figure 5d, digitised (the source prints no coefficients for it).",
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
    # Samineni 2022 Figure 5d fitted line, digitised by the maintainers from the
    # 300-dpi figure image: about 780 pixel columns along the line, converted
    # with the axis-tick calibration and fitted logit-linear in AUC (residual
    # RMSE under 0.008 on the logit scale).
    logit_ref <- -1.375; label("Logit of the event probability at zero venetoclax exposure (unitless logit)") # Figure 5d fitted curve, digitised; p(0) = 0.202
    e_auc_logit <- 0.478; label("Log-odds of the event per 1 ug*day/mL of venetoclax AUCss nominal (unitless logit)") # Figure 5d fitted curve, digitised; figure p = 0.521

    # Logistic regression (glm, logit link): Bernoulli likelihood, no random
    # effects and no residual error. The tiny additive residual below exists
    # only so rxode2 has an error model to attach to the probability; it is
    # not a published quantity.
    addSd_prob_thrombocytopenia_grade3 <- fixed(0.001); label("Placeholder additive residual SD on the event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_thrombocytopenia_grade3 <- logit_ref + e_auc_logit * AUC_VENETOCLAX_SS
    prob_thrombocytopenia_grade3 <- expit(logit_thrombocytopenia_grade3)

    prob_thrombocytopenia_grade3 ~ add(addSd_prob_thrombocytopenia_grade3)
  })
}
