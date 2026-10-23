Song_2024_maribavir_sae <- function() {
  description <- paste0(
    "Univariate binomial logistic-regression exposure-SAFETY model for the ",
    "probability of any treatment-emergent serious adverse event in solid ",
    "organ and haematopoietic cell transplant recipients with refractory ",
    "cytomegalovirus infection treated with maribavir 400 mg orally twice ",
    "daily in the phase 3 SOLSTICE study (Song 2024, n = 231, 88 patients ",
    "with an event). The probability is expit of a linear predictor in the ",
    "0-24 h steady-state maribavir AUC on the last day of exposure (odds ",
    "ratio 1.19 per 50 ug*h/mL, p = 0.000566). The paper publishes this ",
    "relationship only as a fitted curve (Supplementary Figure S1.5); the ",
    "intercept and slope were recovered by digitising that curve. The ",
    "paper's final safety model also retained baseline CMV resistance and ",
    "baseline CMV DNA, whose coefficients are not reported, so this file is ",
    "the univariate exposure-only fit the figure displays. No PK layer and no ",
    "ODE: the exposure metric is supplied as a data column, derived in the ",
    "source analysis from the population PK model packaged as ",
    "Sun_2023_maribavir. No between-subject random effect and no residual ",
    "error are estimated (Bernoulli likelihood). Five companion models in the ",
    "Song_2024_maribavir_* family."
  )
  reference <- paste(
    "Song IH, Chen G, Hayes S, Farrell C, Jomphe C, Gosselin NH, Sun K.",
    "Population pharmacokinetics and exposure-response relationships of",
    "maribavir in transplant recipients with cytomegalovirus infection.",
    "J Pharmacokinet Pharmacodyn. 2024;51:887-904.",
    "doi:10.1007/s10928-024-09939-2.",
    "Coefficients digitised from the fitted curve of Online Supplementary Material 1, Figure S1.5.",
    sep = " "
  )
  vignette <- "Song_2024_maribavir"
  units <- list(
    time = "n/a (static landmark exposure-response regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the AUC_MBV_SS24 covariate column)",
    concentration = "prob_sae (probability of any treatment-emergent serious adverse event, 0-1; also logit_sae)"
  )

  covariateData <- list(
    AUC_MBV_SS24 = list(
      description = "Individual maribavir area under the plasma concentration-time curve from 0 to 24 h at steady state on the last day of exposure. Supplied as data: this model has no PK layer, and the source analysis used individual (empirical-Bayes) predictions from the companion maribavir population PK model.",
      units = "ug*h/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The digitised slope is per 1 ug*h/mL (no rescaling in model()), equal",
        "to 0.174 per 50 ug*h/mL. This is the 24-h integral, NOT the",
        "dosing-interval AUC0-tau of Song 2024 Table 3 and NOT the AUC_MBV_SS",
        "column of Sun 2025, either of which is half as large at steady state.",
        "Derive from modellib('Sun_2023_maribavir'). The figure's quartile",
        "boundaries are 72.2, 170, 238, 354 and 802 ug*h/mL; extrapolation",
        "outside that range is not supported by the source."
      ),
      source_name = "AUCss of maribavir"
    )
  )

  covariatesDataExcluded <- list(
    DIS_CMV_RES_BL = list(
      description = "CMV genotypically resistant to prior anti-CMV therapy at baseline.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not CMV resistant at baseline)",
      notes = "Retained in Song 2024's final serious-adverse-event model (Results) but its coefficient is not reported anywhere in the paper or supplement, so it cannot be encoded; this file is the univariate fit drawn in Figure S1.5."
    ),
    CMVDNA_HIGH = list(
      description = "Baseline plasma CMV DNA >= 9100 IU/mL.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (low baseline CMV DNA)",
      notes = "Baseline CMV DNA level was retained in Song 2024's final serious-adverse-event model (Results), but neither its coding (pooled or three-level) nor its coefficient is reported, so it cannot be encoded."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 231L,
    n_studies = 1L,
    n_observations = "231 binary treatment-emergent-serious-adverse-event records, one per patient (88 with an event; Figure S1.5 quartile table 17 + 20 + 20 + 31)",
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
    # Song 2024 Online Supplementary Material 1, Figure S1.5: the fitted
    # (predicted) probability curve versus AUCss, p = 0.000566.
    #
    #   logit(p) = intercept + slope * AUC_MBV_SS24
    #
    # Digitised by the maintainers from the figure raster: axis ticks
    # calibrated, the fitted curve tracked column by column (1346 points
    # between AUCss 71.5 and 802.9 ug*h/mL), and logit(p) regressed on
    # untransformed AUCss (R^2 = 0.999993, max |residual| 0.014 logit), which
    # confirms the exposure enters linearly and untransformed. Rounded to
    # four significant figures. The vignette checks the curve against the
    # printed quartile table and the observed event count.
    # ================================================================

    logit_ref <- -1.467; label("Logit of the serious-adverse-event probability at zero maribavir exposure (unitless logit)")  # Figure S1.5 fitted curve, digitised; intercept -1.46667 (not printed)

    e_auc_logit <- 0.003481; label("Log-odds of a serious adverse event per 1 ug*h/mL increase in maribavir 0-24 h steady-state AUC (unitless logit per ug*h/mL)")  # Figure S1.5 fitted curve, digitised; slope 0.00348108 per ug*h/mL (not printed; p = 0.000566 printed in the panel)

    # Bernoulli likelihood in the source; the fixed additive residual is a
    # placeholder so rxode2 has an error model, not a published quantity.
    addSd_prob_sae <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_sae <- logit_ref + e_auc_logit * AUC_MBV_SS24

    prob_sae <- expit(logit_sae)

    prob_sae ~ add(addSd_prob_sae)
  })
}
