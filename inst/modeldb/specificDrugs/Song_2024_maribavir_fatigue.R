Song_2024_maribavir_fatigue <- function() {
  description <- paste0(
    "Univariate binomial logistic-regression exposure-SAFETY model for the ",
    "probability of fatigue as a treatment-emergent adverse event in solid ",
    "organ and haematopoietic cell transplant recipients with refractory ",
    "cytomegalovirus infection treated with maribavir 400 mg orally twice ",
    "daily in the phase 3 SOLSTICE study (Song 2024, n = 231, 28 patients ",
    "with the event). The probability is expit of a linear predictor in the ",
    "steady-state maribavir Cmax on the last day of exposure (odds ratio ",
    "1.09 per ug/mL, p = 0.000592). The paper publishes this relationship ",
    "only as a fitted curve (Supplementary Figure S1.4); the intercept and ",
    "slope were recovered by digitising that curve. The paper's final fatigue ",
    "model also retained enrolment region, whose coefficients are not ",
    "reported, so this file is the univariate exposure-only fit the figure ",
    "displays. No PK layer and no ODE: the exposure metric is supplied as a ",
    "data column, derived in the source analysis from the population PK ",
    "model packaged as Sun_2023_maribavir. No between-subject random effect ",
    "and no residual error are estimated (Bernoulli likelihood). Five ",
    "companion models in the Song_2024_maribavir_* family."
  )
  reference <- paste(
    "Song IH, Chen G, Hayes S, Farrell C, Jomphe C, Gosselin NH, Sun K.",
    "Population pharmacokinetics and exposure-response relationships of",
    "maribavir in transplant recipients with cytomegalovirus infection.",
    "J Pharmacokinet Pharmacodyn. 2024;51:887-904.",
    "doi:10.1007/s10928-024-09939-2.",
    "Coefficients digitised from the fitted curve of Online Supplementary Material 1, Figure S1.4.",
    sep = " "
  )
  vignette <- "Song_2024_maribavir"
  units <- list(
    time = "n/a (static landmark exposure-response regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CMAX covariate column)",
    concentration = "prob_fatigue (probability of treatment-emergent fatigue, 0-1; also logit_fatigue)"
  )

  covariateData <- list(
    CMAX = list(
      description = "Individual maribavir maximum plasma concentration at steady state on the last day of exposure. Supplied as data: this model has no PK layer, and the source analysis used individual (empirical-Bayes) predictions from the companion maribavir population PK model.",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The digitised slope is per 1 ug/mL. The figure's quartile boundaries",
        "are 7.37, 12.9, 17.6, 22.6 and 47.0 ug/mL; extrapolation outside that",
        "range is not supported by the source. For calibration, the SOLSTICE",
        "refractory geometric mean Cmax.ss is 16.9 ug/mL (Table 3). Enrolment",
        "region, also retained in the paper's final fatigue model, is not",
        "encoded because neither its levels nor its coefficients are",
        "reported. Derive from modellib('Sun_2023_maribavir')."
      ),
      source_name = "Cmax.ss of maribavir"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 231L,
    n_studies = 1L,
    n_observations = "231 binary treatment-emergent-fatigue records, one per patient (28 with the event; Figure S1.4 quartile table 3 + 5 + 5 + 15)",
    age_range = "19-79 years (median 56.0, mean 53.7, SD 13.4); 18 to <45: 55, 45 to <65: 123, >=65: 53 (Table 4)",
    weight_range = "36.1-124 kg (median 74.1, mean 75.7, SD 17.9) (Table 4)",
    sex_female_pct = 37.7,
    race_ethnicity = c(White = 75.8, Black = 12.6, Asian = 3.9, Other = 6.9, Missing = 0.9),
    disease_state = "Transplant recipients with CMV infection refractory to ganciclovir, valganciclovir, foscarnet or cidofovir, with or without genotypic resistance: solid organ transplant 140 (61%), haematopoietic cell transplant 91 (39%). CMV asymptomatic 199, symptomatic 14, tissue-invasive disease 18. Baseline CMV DNA low (<9100 IU/mL) 151, intermediate 66, high (>=91 000 IU/mL) 14. Resistant to prior anti-CMV therapy at baseline 100, not resistant 126, not reported 5. Prior CMV prophylaxis 98 (Tables 4 and S1.2).",
    dose_range = "Maribavir 400 mg orally twice daily for 8 weeks, the SOLSTICE randomized dose",
    regions = "Multicentre phase 3 trial (NCT02931539); enrolment regions are not tabulated for the exposure-response set.",
    notes = paste0(
      "The exposure-response analysis set: 231 of the 235 maribavir-arm ",
      "patients of the phase 3 SOLSTICE study who had PK exposure ",
      "parameters. Adverse events were counted over the 8-week treatment ",
      "phase."
    )
  )

  ini({
    # ================================================================
    # Song 2024 Online Supplementary Material 1, Figure S1.4: the fitted
    # (predicted) probability curve versus Cmax.ss, p = 0.000592.
    #
    #   logit(p) = intercept + slope * CMAX
    #
    # Digitised by the maintainers from the figure raster: axis ticks
    # calibrated, the fitted curve tracked column by column (1338 points
    # between Cmax.ss 7.39 and 47.0 ug/mL), and logit(p) regressed on
    # untransformed Cmax.ss (R^2 = 0.99995, max |residual| 0.11 logit, the
    # largest residuals at the lowest probabilities where one pixel spans
    # the most logit), which confirms the exposure enters linearly and
    # untransformed. Rounded to four significant figures.
    # ================================================================

    logit_ref <- -3.713; label("Logit of the fatigue probability at zero maribavir Cmax (unitless logit)")  # Figure S1.4 fitted curve, digitised; intercept -3.71285 (not printed)

    e_cmax_logit <- 0.08244; label("Log-odds of fatigue per 1 ug/mL increase in maribavir steady-state Cmax (unitless logit per ug/mL)")  # Figure S1.4 fitted curve, digitised; slope 0.0824396 per ug/mL (not printed; p = 0.000592 printed in the panel)

    # Bernoulli likelihood in the source; the fixed additive residual is a
    # placeholder so rxode2 has an error model, not a published quantity.
    addSd_prob_fatigue <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_fatigue <- logit_ref + e_cmax_logit * CMAX

    prob_fatigue <- expit(logit_fatigue)

    prob_fatigue ~ add(addSd_prob_fatigue)
  })
}
