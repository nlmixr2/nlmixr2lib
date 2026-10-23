Song_2024_maribavir_cmv_clearance_wk8 <- function() {
  description <- paste0(
    "Binomial logistic-regression exposure-EFFICACY model for the PRIMARY ",
    "endpoint of the phase 3 SOLSTICE study -- confirmed clearance of plasma ",
    "cytomegalovirus (CMV) DNA at week 8 -- in solid organ and haematopoietic ",
    "cell transplant recipients with CMV infection refractory (with or ",
    "without resistance) to ganciclovir, valganciclovir, foscarnet or ",
    "cidofovir, treated with maribavir 400 mg orally twice daily (Song 2024, ",
    "n = 231, the SOLSTICE maribavir arm). The probability of clearance is ",
    "expit of a linear predictor in the 0-24 h steady-state maribavir AUC on ",
    "the last day of exposure, an on-treatment treatment-emergent ",
    "maribavir-resistance mutation, and the baseline CD8+CD69+ pp65 T-cell ",
    "stratum. The exposure slope is NEGATIVE (odds ratio 0.877 per 50 ",
    "ug*h/mL); the source attributes it to resistance mutations clustering in ",
    "the upper exposure quartiles, and it was not significant once those ",
    "patients were excluded. No PK layer and no ODE: the exposure metric is ",
    "supplied as a data column, derived in the source analysis from the ",
    "population PK model packaged as Sun_2023_maribavir. No between-subject ",
    "random effect and no residual error are estimated (Bernoulli ",
    "likelihood). Five companion models in the Song_2024_maribavir_* family."
  )
  reference <- paste(
    "Song IH, Chen G, Hayes S, Farrell C, Jomphe C, Gosselin NH, Sun K.",
    "Population pharmacokinetics and exposure-response relationships of",
    "maribavir in transplant recipients with cytomegalovirus infection.",
    "J Pharmacokinet Pharmacodyn. 2024;51:887-904.",
    "doi:10.1007/s10928-024-09939-2.",
    "Logistic-regression coefficients in Table 5.",
    sep = " "
  )
  vignette <- "Song_2024_maribavir"
  units <- list(
    time = "n/a (static landmark exposure-response regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the AUC_MBV_SS24 covariate column)",
    concentration = "prob_cmv_clearance_wk8 (probability of confirmed CMV clearance of plasma CMV DNA at week 8, 0-1; also logit_cmv_clearance_wk8)"
  )

  covariateData <- list(
    AUC_MBV_SS24 = list(
      description = "Individual maribavir area under the plasma concentration-time curve from 0 to 24 h at steady state on the last day of exposure. Supplied as data: this model has no PK layer, and the source analysis used individual (empirical-Bayes) predictions from the companion maribavir population PK model.",
      units = "ug*h/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Song 2024 estimates the coefficient PER INCREMENT OF 50 ug*h/mL, so",
        "model() divides this column by 50; a user supplying the raw AUC needs no",
        "rescaling of their own. This is the 24-h integral (both doses of the",
        "twice-daily regimen), NOT the dosing-interval AUC0-tau of Song 2024",
        "Table 3 and NOT the AUC_MBV_SS column of Sun 2025, either of which is",
        "half as large at steady state. Derive the column from",
        "modellib('Sun_2023_maribavir'), which encodes the same NONMEM run",
        "(RUN# 171) Song 2024 used; integrate Cc over the 24 h of the last",
        "dosing day. The SOLSTICE analysis-set quartile boundaries are 72.2,",
        "170, 238, 354 and 802 ug*h/mL (Figure S1.5).",
        "Slope -0.131 per 50 ug*h/mL (odds ratio 0.877, 95% CI 0.780-0.982,",
        "p = 0.0233)."
      ),
      source_name = "AUCss of maribavir - increment of 50 ug.h/mL"
    ),
    TE_RESIST_MBV_ONTX = list(
      description = "Treatment-emergent CMV mutation conferring resistance to maribavir, detected during maribavir treatment.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no treatment-emergent maribavir-resistance mutation detected during treatment)",
      notes = paste(
        "47 of 231 patients positive; a further 13 had a mutation first",
        "detected after the end of treatment and are 0 here (Table 4). This is",
        "an ON-STUDY OUTCOME used as a predictor of another on-study outcome,",
        "so the coefficient is a descriptive association, not a quantity",
        "available for prognosis before treatment. All 47 were week-8",
        "non-responders. Coefficient -5.39 (SE 1.42), odds ratio 0.00457: by",
        "far the strongest term. Song 2024 Results: probability of response",
        "1.3% with versus 74.2% without, at CD8+CD69+ < 0.5% and AUCss 200",
        "ug*h/mL."
      ),
      source_name = "Treatment-emergent maribavir resistance on maribavir treatment"
    ),
    CD8_PP65_MID = list(
      description = "Baseline CD8+CD69+ pp65-stimulated T-cell percentage in the >= 0.5% to < 2% stratum.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (baseline CD8+CD69+ < 0.5%, which also holds results reported as 'too few cells for analysis')",
      notes = paste(
        "28 of 231 patients (Table S1.2). One of three indicators decomposing",
        "the four-level baseline CMV-specific T-cell assay; the other two are",
        "CD8_PP65_HI and CD8_PP65_NR, and all three are 0 for the < 0.5%",
        "reference. Coefficient -0.304, not significant (p = 0.535)."
      ),
      source_name = "Baseline CD8+CD69+ cell count, >= 0.5% to < 2%"
    ),
    CD8_PP65_HI = list(
      description = "Baseline CD8+CD69+ pp65-stimulated T-cell percentage in the >= 2% stratum.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (baseline CD8+CD69+ < 0.5%, which also holds results reported as 'too few cells for analysis')",
      notes = paste(
        "29 of 231 patients (Table S1.2). Coefficient +1.19 (p = 0.0371):",
        "probability of response 90.0% versus 74.2% in the reference stratum,",
        "without a resistance mutation and at AUCss 200 ug*h/mL (Results)."
      ),
      source_name = "Baseline CD8+CD69+ cell count, >= 2%"
    ),
    CD8_PP65_NR = list(
      description = "Baseline CD8+CD69+ pp65-stimulated T-cell assay not reported.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (baseline CD8+CD69+ < 0.5%, which also holds results reported as 'too few cells for analysis')",
      notes = paste(
        "46 of 231 patients (20%). A FITTED level with its own coefficient, not",
        "a missing-data code: give an unreported patient CD8_PP65_NR = 1 rather",
        "than imputing them into the reference. A result reported as 'too few",
        "cells for analysis' is NOT 'not reported' -- Table S1.2 footnote a",
        "assigns it to the < 0.5% reference. Coefficient -0.314 (p = 0.447)."
      ),
      source_name = "Baseline CD8+CD69+ cell count, Not reported"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 231L,
    n_studies = 1L,
    n_observations = "231 binary confirmed-CMV-clearance-at-week-8 records, one per patient (131 responders)",
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
      "parameters. Individual exposures were derived from the final ",
      "population PK model (Song 2024 Table 2; the same NONMEM run as ",
      "Sun_2023_maribavir) using the actual dosing history and time-varying ",
      "strong CYP3A inhibitor / inducer use. Table 4 labels its youngest age ",
      "band '18 to < 65', which from the counts and the next band is a ",
      "misprint for 18 to < 45."
    )
  )

  ini({
    # ================================================================
    # Song 2024 Table 5 (primary endpoint).
    #
    #   logit(p) = intercept + slope * (AUC_MBV_SS24 / 50)
    #              + risk-factor terms
    #
    # Every value is the 'Estimate (SE)' column. Covariates are not
    # centred, so the intercept is the logit at zero exposure with every
    # risk factor at its reference level. The Results-text scenario
    # probabilities (74.2%, 1.3%, 90.0% at AUCss 200 ug*h/mL) reproduce from
    # these values and are re-checked in the validation vignette.
    # ================================================================

    logit_ref <- 1.58; label("Logit of the confirmed CMV clearance at week 8 probability at zero maribavir exposure with all risk factors at their reference level (unitless logit)")  # Table 5, Intercept 1.58 (SE 0.413), p < 0.001

    e_auc_logit <- -0.131; label("Log-odds of confirmed CMV clearance at week 8 per 50 ug*h/mL increase in maribavir 0-24 h steady-state AUC (unitless logit)")  # Table 5, AUCss of maribavir - increment of 50 ug.h/mL, -0.131 (SE 0.0588), OR 0.877 (0.780, 0.982), p = 0.0233

    e_te_resist_mbv_ontx_logit <- -5.39; label("Log-odds of confirmed CMV clearance at week 8 for a treatment-emergent maribavir-resistance mutation detected during treatment (unitless logit)")  # Table 5, Treatment-emergent maribavir resistance on maribavir treatment -5.39 (SE 1.42), OR 0.00457, p < 0.001
    e_cd8_pp65_mid_logit <- -0.304; label("Log-odds of confirmed CMV clearance at week 8 for the middle baseline CD8+CD69+ pp65 stimulation stratum versus the lowest stratum (unitless logit)")  # Table 5, Baseline CD8+CD69+ cell count >= 0.5% to < 2%, -0.304 (SE 0.492), p = 0.535; the printed OR 0.739 does not match exp(-0.304) = 0.738, and the estimate column is used
    e_cd8_pp65_hi_logit <- 1.19; label("Log-odds of confirmed CMV clearance at week 8 for the highest baseline CD8+CD69+ pp65 stimulation stratum versus the lowest stratum (unitless logit)")  # Table 5, Baseline CD8+CD69+ cell count >= 2%, 1.19 (SE 0.630), OR 3.29, p = 0.0371
    e_cd8_pp65_nr_logit <- -0.314; label("Log-odds of confirmed CMV clearance at week 8 for an unreported baseline CD8+CD69+ pp65 assay versus the lowest stratum (unitless logit)")  # Table 5, Baseline CD8+CD69+ cell count Not reported, -0.314 (SE 0.415), OR 0.730, p = 0.447

    # The source is a binomial logistic regression (Bernoulli likelihood, no
    # random effects). The tiny fixed additive residual exists only so
    # rxode2 has an error model to attach to the probability; it is NOT a
    # published quantity.
    addSd_prob_cmv_clearance_wk8 <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # Table 5 estimates the exposure coefficient per 50 ug*h/mL increment.
    auc50 <- AUC_MBV_SS24 / 50

    logit_cmv_clearance_wk8 <- logit_ref +
      e_auc_logit * auc50 +
      e_te_resist_mbv_ontx_logit * TE_RESIST_MBV_ONTX +
      e_cd8_pp65_mid_logit * CD8_PP65_MID +
      e_cd8_pp65_hi_logit * CD8_PP65_HI +
      e_cd8_pp65_nr_logit * CD8_PP65_NR

    prob_cmv_clearance_wk8 <- expit(logit_cmv_clearance_wk8)

    prob_cmv_clearance_wk8 ~ add(addSd_prob_cmv_clearance_wk8)
  })
}
