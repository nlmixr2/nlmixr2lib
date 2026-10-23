Song_2024_maribavir_cmv_clearance_wk16 <- function() {
  description <- paste0(
    "Binomial logistic-regression exposure-EFFICACY model for the KEY ",
    "SECONDARY endpoint of the phase 3 SOLSTICE study -- confirmed clearance ",
    "of plasma cytomegalovirus (CMV) DNA AND CMV infection symptom control at ",
    "week 8, maintained through week 16 -- in solid organ and haematopoietic ",
    "cell transplant recipients with refractory CMV infection treated with ",
    "maribavir 400 mg orally twice daily (Song 2024, n = 231, the SOLSTICE ",
    "maribavir arm). The probability of response is expit of a linear ",
    "predictor in the 0-24 h steady-state maribavir AUC on the last day of ",
    "exposure, a treatment-emergent maribavir-resistance mutation detected at ",
    "any time, the baseline CD8+CD69+ pp65 T-cell stratum, an ",
    "intermediate-or-high baseline CMV DNA load and prior CMV prophylaxis. ",
    "The exposure slope is NEGATIVE (odds ratio 0.740 per 50 ug*h/mL). No PK ",
    "layer and no ODE: the exposure metric is supplied as a data column, ",
    "derived in the source analysis from the population PK model packaged as ",
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
    "Logistic-regression coefficients in Online Supplementary Material 1, Table S1.3.",
    sep = " "
  )
  vignette <- "Song_2024_maribavir"
  units <- list(
    time = "n/a (static landmark exposure-response regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the AUC_MBV_SS24 covariate column)",
    concentration = "prob_cmv_clearance_wk16 (probability of confirmed CMV clearance and symptom control at week 8 maintained through week 16, 0-1; also logit_cmv_clearance_wk16)"
  )

  covariateData <- list(
    AUC_MBV_SS24 = list(
      description = "Individual maribavir area under the plasma concentration-time curve from 0 to 24 h at steady state on the last day of exposure. Supplied as data: this model has no PK layer, and the source analysis used individual (empirical-Bayes) predictions from the companion maribavir population PK model.",
      units = "ug*h/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Song 2024 estimates the coefficient PER INCREMENT OF 50 ug*h/mL, so",
        "model() divides this column by 50. This is the 24-h integral (both",
        "doses of the twice-daily regimen), NOT the dosing-interval AUC0-tau of",
        "Song 2024 Table 3 and NOT the AUC_MBV_SS column of Sun 2025, either of",
        "which is half as large at steady state. Derive the column from",
        "modellib('Sun_2023_maribavir'). Slope -0.301 per 50 ug*h/mL (odds",
        "ratio 0.740, 95% CI 0.611-0.870, p < 0.001); observed response rates",
        "fall from 31% to 5% across the exposure quartiles (Results)."
      ),
      source_name = "AUCss of maribavir - increment of 50 ug.h/mL"
    ),
    TE_RESIST_MBV = list(
      description = "Treatment-emergent CMV mutation conferring resistance to maribavir, detected at any time (during treatment or follow-up).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no treatment-emergent maribavir-resistance mutation at any time)",
      notes = paste(
        "60 of 231 patients: 47 detected during treatment and 13 after it",
        "(Table 4). Unlike the primary-endpoint model, which uses the",
        "during-treatment subset (TE_RESIST_MBV_ONTX), this model uses the",
        "ANY-TIME definition, matching an endpoint whose maintenance window",
        "runs into follow-up. An on-study outcome used as a predictor, so the",
        "coefficient is a descriptive association. Coefficient -3.27 (SE",
        "1.39), odds ratio 0.0381."
      ),
      source_name = "Treatment-emergent CMV mutation conferring resistance to maribavir"
    ),
    CD8_PP65_MID = list(
      description = "Baseline CD8+CD69+ pp65-stimulated T-cell percentage in the >= 0.5% to < 2% stratum.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (baseline CD8+CD69+ < 0.5%, which also holds results reported as 'too few cells for analysis')",
      notes = "28 of 231 patients (Table S1.2). Coefficient -0.521, not significant (p = 0.387).",
      source_name = "Baseline CD8+CD69+ cell count, >= 0.5% to < 2%"
    ),
    CD8_PP65_HI = list(
      description = "Baseline CD8+CD69+ pp65-stimulated T-cell percentage in the >= 2% stratum.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (baseline CD8+CD69+ < 0.5%, which also holds results reported as 'too few cells for analysis')",
      notes = "29 of 231 patients (Table S1.2). Coefficient +1.65 (p = 0.002): probability of response 77.9% versus 40.3% in the reference stratum (Results).",
      source_name = "Baseline CD8+CD69+ cell count, >= 2%"
    ),
    CD8_PP65_NR = list(
      description = "Baseline CD8+CD69+ pp65-stimulated T-cell assay not reported.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (baseline CD8+CD69+ < 0.5%, which also holds results reported as 'too few cells for analysis')",
      notes = paste(
        "46 of 231 patients (20%). A FITTED level, not a missing-data code:",
        "give an unreported patient CD8_PP65_NR = 1 rather than imputing them",
        "into the reference. Coefficient -0.107 (p = 0.830)."
      ),
      source_name = "Baseline CD8+CD69+ cell count, Not reported"
    ),
    CMVDNA_HIGH = list(
      description = "Baseline plasma CMV DNA >= 9100 IU/mL (SOLSTICE intermediate or high stratum, pooled).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (low baseline CMV DNA, < 9100 IU/mL)",
      notes = paste(
        "80 of 231 patients (66 intermediate + 14 high, Table 4). Song 2024",
        "does not print the cut-point; the SOLSTICE primary publication (Avery",
        "2022 Clin Infect Dis) defines low as < 9100 IU/mL and pools",
        "intermediate with high as >= 9100 IU/mL, which is the CMVDNA_HIGH",
        "threshold. Coefficient -1.01 (p = 0.030): probability of response",
        "19.7% versus 40.3% (Results)."
      ),
      source_name = "Baseline CMV DNA level intermediate/high"
    ),
    CMV_PROPHY_PRIOR = list(
      description = "Prior use of CMV prophylaxis.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no prior CMV prophylaxis)",
      notes = "98 of 231 patients (Table S1.2). Coefficient -0.818 (p = 0.041): probability of response 23.0% versus 40.3% (Results).",
      source_name = "Prior use of CMV prophylaxis"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 231L,
    n_studies = 1L,
    n_observations = "231 binary key-secondary-endpoint records, one per patient (44 responders, Results)",
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
      "population PK model (the same NONMEM run as Sun_2023_maribavir). ",
      "44 of 231 (19%) met the key secondary endpoint (Results)."
    )
  )

  ini({
    # ================================================================
    # Song 2024 Online Supplementary Material 1, Table S1.3 (key
    # secondary endpoint).
    #
    #   logit(p) = intercept + slope * (AUC_MBV_SS24 / 50)
    #              + risk-factor terms
    #
    # Covariates are not centred. The Results-text scenario probabilities
    # (40.3%, 2.5%, 77.9%, 19.7%, 23.0% at AUCss 200 ug*h/mL) reproduce from
    # these values and are re-checked in the validation vignette.
    # ================================================================

    logit_ref <- 0.814; label("Logit of the key secondary endpoint probability at zero maribavir exposure with all risk factors at their reference level (unitless logit)")  # Table S1.3, Intercept 0.814 (SE 0.543), p = 0.128

    e_auc_logit <- -0.301; label("Log-odds of the key secondary endpoint per 50 ug*h/mL increase in maribavir 0-24 h steady-state AUC (unitless logit)")  # Table S1.3, AUCss of maribavir - increment of 50 ug.h/mL, -0.301 (SE 0.0884), OR 0.740 (0.611, 0.870), p < 0.001

    e_te_resist_mbv_logit <- -3.27; label("Log-odds of the key secondary endpoint for a treatment-emergent maribavir-resistance mutation detected at any time (unitless logit)")  # Table S1.3, Treatment-emergent CMV mutation conferring resistance to maribavir -3.27 (SE 1.39), OR 0.0381, p < 0.001
    e_cd8_pp65_mid_logit <- -0.521; label("Log-odds of the key secondary endpoint for the middle baseline CD8+CD69+ pp65 stimulation stratum versus the lowest stratum (unitless logit)")  # Table S1.3, Baseline CD8+CD69+ cell count >= 0.5% to < 2%, -0.521 (SE 0.624), OR 0.594, p = 0.387
    e_cd8_pp65_hi_logit <- 1.65; label("Log-odds of the key secondary endpoint for the highest baseline CD8+CD69+ pp65 stimulation stratum versus the lowest stratum (unitless logit)")  # Table S1.3, Baseline CD8+CD69+ cell count >= 2%, 1.65 (SE 0.552), OR 5.18, p = 0.002
    e_cd8_pp65_nr_logit <- -0.107; label("Log-odds of the key secondary endpoint for an unreported baseline CD8+CD69+ pp65 assay versus the lowest stratum (unitless logit)")  # Table S1.3, Baseline CD8+CD69+ cell count Not reported, -0.107 (SE 0.499), OR 0.899, p = 0.830
    e_cmvdna_high_logit <- -1.01; label("Log-odds of the key secondary endpoint for intermediate or high baseline CMV DNA (>= 9100 IU/mL) versus low (unitless logit)")  # Table S1.3, Baseline CMV DNA level intermediate/high -1.01 (SE 0.489), OR 0.363, p = 0.030
    e_cmv_prophy_prior_logit <- -0.818; label("Log-odds of the key secondary endpoint for prior CMV prophylaxis (unitless logit)")  # Table S1.3, Prior use of CMV prophylaxis -0.818 (SE 0.409), OR 0.441, p = 0.041

    # Bernoulli likelihood in the source; the fixed additive residual is a
    # placeholder so rxode2 has an error model, not a published quantity.
    addSd_prob_cmv_clearance_wk16 <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # Table S1.3 estimates the exposure coefficient per 50 ug*h/mL increment.
    auc50 <- AUC_MBV_SS24 / 50

    logit_cmv_clearance_wk16 <- logit_ref +
      e_auc_logit * auc50 +
      e_te_resist_mbv_logit * TE_RESIST_MBV +
      e_cd8_pp65_mid_logit * CD8_PP65_MID +
      e_cd8_pp65_hi_logit * CD8_PP65_HI +
      e_cd8_pp65_nr_logit * CD8_PP65_NR +
      e_cmvdna_high_logit * CMVDNA_HIGH +
      e_cmv_prophy_prior_logit * CMV_PROPHY_PRIOR

    prob_cmv_clearance_wk16 <- expit(logit_cmv_clearance_wk16)

    prob_cmv_clearance_wk16 ~ add(addSd_prob_cmv_clearance_wk16)
  })
}
