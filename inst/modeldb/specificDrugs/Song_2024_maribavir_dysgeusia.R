Song_2024_maribavir_dysgeusia <- function() {
  description <- paste0(
    "Univariate binomial logistic-regression exposure-SAFETY model for the ",
    "probability of taste disturbance (dysgeusia) as a treatment-emergent ",
    "adverse event in solid organ and haematopoietic cell transplant ",
    "recipients with refractory cytomegalovirus infection treated with ",
    "maribavir 400 mg orally twice daily in the phase 3 SOLSTICE study (Song ",
    "2024, n = 231, 108 patients with the event). The probability is expit of ",
    "a linear predictor in the average maribavir concentration, with a ",
    "NEGATIVE slope (odds ratio 0.923 per ug/mL, p = 0.00202). The paper ",
    "publishes this relationship only as a fitted curve (Supplementary ",
    "Figure S1.3); the intercept and slope were recovered by digitising that ",
    "curve. The paper's final taste-disturbance model also retained ",
    "transplant type and enrolment region, whose coefficients are not ",
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
    "Coefficients digitised from the fitted curve of Online Supplementary Material 1, Figure S1.3.",
    sep = " "
  )
  vignette <- "Song_2024_maribavir"
  units <- list(
    time = "n/a (static landmark exposure-response regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "prob_dysgeusia (probability of treatment-emergent taste disturbance, 0-1; also logit_dysgeusia)"
  )

  covariateData <- list(
    CAV = list(
      description = "Individual maribavir average plasma concentration. Supplied as data: this model has no PK layer, and the source analysis used individual (empirical-Bayes) predictions from the companion maribavir population PK model.",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Song 2024 defines Cavg as the 'average concentration on each study",
        "day' and plots one value per patient; at steady state on the 400 mg",
        "twice-daily regimen the daily average equals the dosing-interval",
        "average (AUC0-tau / 12 h, or AUC_MBV_SS24 / 24 h). The paper does not",
        "state how the per-day values were summarised into one per patient.",
        "The digitised slope is per 1 ug/mL. The figure's quartile boundaries",
        "are 3.50, 7.83, 10.6, 14.6 and 38.3 ug/mL; extrapolation outside that",
        "range is not supported by the source. Derive from",
        "modellib('Sun_2023_maribavir')."
      ),
      source_name = "Cavg of maribavir"
    )
  )

  covariatesDataExcluded <- list(
    TX_HCT = list(
      description = "Haematopoietic cell transplant recipient (1) versus solid organ transplant recipient (0).",
      units = "(binary)",
      type = "binary",
      reference_category = "not reported",
      notes = "Transplant type was retained in Song 2024's final taste-disturbance model together with enrolment region (Results), but neither coefficient nor the region levels are reported, so they cannot be encoded; this file is the univariate fit drawn in Figure S1.3."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 231L,
    n_studies = 1L,
    n_observations = "231 binary treatment-emergent-taste-disturbance records, one per patient (108 with the event; Figure S1.3 quartile table 31 + 33 + 26 + 18)",
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
    # Song 2024 Online Supplementary Material 1, Figure S1.3: the fitted
    # (predicted) probability curve versus Cavg, p = 0.00202.
    #
    #   logit(p) = intercept + slope * CAV
    #
    # Digitised by the maintainers from the figure raster: axis ticks
    # calibrated, the fitted curve tracked column by column (1346 points
    # between Cavg 3.48 and 38.4 ug/mL), and logit(p) regressed on
    # untransformed Cavg (R^2 = 0.99999, max |residual| 0.012 logit), which
    # confirms the exposure enters linearly and untransformed. Rounded to
    # four significant figures.
    # ================================================================

    logit_ref <- 0.8056; label("Logit of the taste-disturbance probability at zero maribavir average concentration (unitless logit)")  # Figure S1.3 fitted curve, digitised; intercept 0.805631 (not printed)

    e_cav_logit <- -0.07977; label("Log-odds of taste disturbance per 1 ug/mL increase in maribavir average concentration (unitless logit per ug/mL)")  # Figure S1.3 fitted curve, digitised; slope -0.0797684 per ug/mL (not printed; p = 0.00202 printed in the panel)

    # Bernoulli likelihood in the source; the fixed additive residual is a
    # placeholder so rxode2 has an error model, not a published quantity.
    addSd_prob_dysgeusia <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_dysgeusia <- logit_ref + e_cav_logit * CAV

    prob_dysgeusia <- expit(logit_dysgeusia)

    prob_dysgeusia ~ add(addSd_prob_dysgeusia)
  })
}
