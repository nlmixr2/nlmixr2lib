Song_2024_maribavir_te_resistance <- function() {
  description <- paste0(
    "Binomial logistic-regression model for the probability of developing a ",
    "treatment-emergent cytomegalovirus (CMV) mutation conferring resistance ",
    "to maribavir (detected during treatment or follow-up) in solid organ and ",
    "haematopoietic cell transplant recipients with refractory CMV infection ",
    "treated with maribavir 400 mg orally twice daily in the phase 3 SOLSTICE ",
    "study (Song 2024, n = 231, 60 events). The probability is expit of a ",
    "linear predictor in the steady-state maribavir trough concentration on ",
    "the last day of exposure, the baseline CMV DNA stratum (low, ",
    "intermediate, high), CMV resistance to prior anti-CMV therapy at ",
    "baseline, the baseline CD4+CD69+ pp65 T-cell stratum, age and transplant ",
    "type. The exposure slope is POSITIVE (odds ratio 1.23 per 2.5 ug/mL): ",
    "higher trough, higher resistance probability. No PK layer and no ODE: ",
    "the exposure metric is supplied as a data column, derived in the source ",
    "analysis from the population PK model packaged as Sun_2023_maribavir. ",
    "No between-subject random effect and no residual error are estimated ",
    "(Bernoulli likelihood). Five companion models in the ",
    "Song_2024_maribavir_* family."
  )
  reference <- paste(
    "Song IH, Chen G, Hayes S, Farrell C, Jomphe C, Gosselin NH, Sun K.",
    "Population pharmacokinetics and exposure-response relationships of",
    "maribavir in transplant recipients with cytomegalovirus infection.",
    "J Pharmacokinet Pharmacodyn. 2024;51:887-904.",
    "doi:10.1007/s10928-024-09939-2.",
    "Logistic-regression coefficients in Table 6.",
    sep = " "
  )
  vignette <- "Song_2024_maribavir"
  units <- list(
    time = "n/a (static landmark exposure-response regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CTROUGH covariate column)",
    concentration = "prob_te_resist_mbv (probability of developing a treatment-emergent maribavir-resistant CMV mutation, 0-1; also logit_te_resist_mbv)"
  )

  covariateData <- list(
    CTROUGH = list(
      description = "Individual maribavir minimum (trough) plasma concentration at steady state on the last day of exposure. Supplied as data: this model has no PK layer, and the source analysis used individual (empirical-Bayes) predictions from the companion maribavir population PK model.",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Song 2024 estimates the coefficient PER INCREMENT OF 2.5 ug/mL, so",
        "model() divides this column by 2.5. Derive the column from",
        "modellib('Sun_2023_maribavir'), which encodes the same NONMEM run",
        "(RUN# 171) Song 2024 used. For calibration, the SOLSTICE refractory",
        "geometric mean Ctrough.ss is 5.13 ug/mL (Table 3). Slope +0.210 per",
        "2.5 ug/mL (odds ratio 1.23, 95% CI 1.06-1.43, p = 0.006). A",
        "complementary analysis of observed predose concentrations found",
        "trough levels higher in patients who later developed a mutation from",
        "week 1 onward, which the source reads as exposure preceding (not",
        "caused by) the mutation."
      ),
      source_name = "Ctrough.ss of maribavir - increment of 2.5 ug/mL"
    ),
    CMVDNA_HIGH = list(
      description = "Baseline plasma CMV DNA >= 9100 IU/mL (SOLSTICE intermediate or high stratum).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (low baseline CMV DNA, < 9100 IU/mL)",
      notes = paste(
        "Used together with the nested CMVDNA_VHIGH to encode Table 6's",
        "three-level categorical (low reference): intermediate is",
        "CMVDNA_HIGH = 1 with CMVDNA_VHIGH = 0, high is both = 1. Cut-points",
        "from the SOLSTICE primary publication (Avery 2022 Clin Infect Dis):",
        "low < 9100, intermediate >= 9100 to < 91 000, high >= 91 000 IU/mL.",
        "151 / 66 / 14 of 231 patients (Table 4). Intermediate coefficient",
        "+0.694 (p = 0.056)."
      ),
      source_name = "Baseline CMV DNA levels, Intermediate"
    ),
    CMVDNA_VHIGH = list(
      description = "Baseline plasma CMV DNA >= 91 000 IU/mL (SOLSTICE high stratum).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (baseline CMV DNA < 91 000 IU/mL)",
      notes = paste(
        "Nested inside CMVDNA_HIGH: a patient with CMVDNA_VHIGH = 1 must also",
        "have CMVDNA_HIGH = 1. 14 of 231 patients. High-versus-low coefficient",
        "+1.52 (p = 0.016): probability 84.5% versus 54.5% for low CMV DNA in",
        "the Results scenario."
      ),
      source_name = "Baseline CMV DNA levels, High"
    ),
    DIS_CMV_RES_BL = list(
      description = "CMV genotypically resistant to prior anti-CMV therapy at baseline.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not CMV resistant at baseline)",
      notes = paste(
        "100 resistant, 126 not resistant and 5 not reported of 231 (Table",
        "S1.2); the source does not state how the 5 unreported patients were",
        "coded. Resistance to the prior agents (ganciclovir / valganciclovir,",
        "foscarnet, cidofovir), not to maribavir. Coefficient -1.08",
        "(p = 0.006): probability 28.9% versus 54.4% in the Results scenario."
      ),
      source_name = "Baseline CMV resistant"
    ),
    CD4_PP65_MID = list(
      description = "Baseline CD4+CD69+ pp65-stimulated T-cell percentage in the >= 0.5% to < 2% stratum.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (baseline CD4+CD69+ < 0.5%, which also holds results reported as 'too few cells for analysis')",
      notes = paste(
        "42 of 231 patients (Table S1.2). Coefficient -1.52 (p = 0.010):",
        "probability 20.7% versus 54.5% in the Results scenario."
      ),
      source_name = "Baseline CD4+CD69+ cell count, >= 0.5% to < 2%"
    ),
    CD4_PP65_HI = list(
      description = "Baseline CD4+CD69+ pp65-stimulated T-cell percentage in the >= 2% stratum.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (baseline CD4+CD69+ < 0.5%, which also holds results reported as 'too few cells for analysis')",
      notes = "12 of 231 patients (Table S1.2). Coefficient -0.0534, not significant (p = 0.945).",
      source_name = "Baseline CD4+CD69+ cell count, >= 2%"
    ),
    CD4_PP65_NR = list(
      description = "Baseline CD4+CD69+ pp65-stimulated T-cell assay not reported.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (baseline CD4+CD69+ < 0.5%, which also holds results reported as 'too few cells for analysis')",
      notes = paste(
        "46 of 231 patients (20%). A FITTED level, not a missing-data code:",
        "give an unreported patient CD4_PP65_NR = 1 rather than imputing them",
        "into the reference. 'Too few cells for analysis' is NOT 'not",
        "reported' -- Table S1.2 footnote a assigns it to the < 0.5%",
        "reference. Coefficient -0.185 (p = 0.666)."
      ),
      source_name = "Baseline CD4+CD69+ cell count, Not reported"
    ),
    AGE = list(
      description = "Age at baseline.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters UNCENTRED, per year: coefficient -0.0328 (odds ratio 0.968,",
        "p = 0.011). The mean predicted probability over the cohort's printed",
        "covariate marginals is 0.262 with the uncentred term, against the",
        "observed 60/231 = 0.260; centring at the median age of 56 would give",
        "about 0.6. Note that the Results-text scenario probabilities (54.5%",
        "etc.) are evaluated with the age term at zero."
      ),
      source_name = "Age"
    ),
    TX_HCT = list(
      description = "Haematopoietic cell transplant recipient (1) versus solid organ transplant recipient (0).",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (HCT recipient) is the source's reference; the published effect is on the complement, solid organ transplant",
      notes = paste(
        "Song 2024 Table 6 estimates +0.970 for 'Solid organ transplant'",
        "against an HCT reference; model() forms the SOT indicator as",
        "(1 - TX_HCT). This identity holds only because the cohort is exactly",
        "140 SOT + 91 HCT = 231 transplant recipients with no non-transplant",
        "subject. Probability 75.9% for SOT versus 54.4% for HCT in the",
        "Results scenario."
      ),
      source_name = "Solid organ transplant"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 231L,
    n_studies = 1L,
    n_observations = "231 binary treatment-emergent-maribavir-resistance records, one per patient (60 events: 47 during treatment, 13 after it)",
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
      "population PK model (the same NONMEM run as Sun_2023_maribavir)."
    )
  )

  ini({
    # ================================================================
    # Song 2024 Table 6.
    #
    #   logit(p) = intercept + slope * (CTROUGH / 2.5)
    #              + b_int * intermediate + b_high * high CMV DNA
    #              + other risk-factor terms
    #
    # Covariates are not centred; AGE enters per year from zero. The
    # Results-text scenario probabilities (54.5%, 70.5%, 84.5%, 28.9%,
    # 20.7%, 75.9% at Ctrough.ss 5 ug/mL with the age term at zero)
    # reproduce from these values and are re-checked in the vignette.
    # ================================================================

    logit_ref <- -0.242; label("Logit of the treatment-emergent maribavir-resistance probability at zero trough concentration and zero age with all risk factors at their reference level (unitless logit)")  # Table 6, Intercept -0.242 (SE 0.711), p = 0.733

    e_ctrough_logit <- 0.210; label("Log-odds of treatment-emergent maribavir resistance per 2.5 ug/mL increase in steady-state trough concentration (unitless logit)")  # Table 6, Ctrough.ss of maribavir - increment of 2.5 ug/mL, 0.210 (SE 0.0762), OR 1.23 (1.06, 1.43), p = 0.006

    e_cmvdna_int_logit <- 0.694; label("Log-odds of treatment-emergent maribavir resistance for intermediate versus low baseline CMV DNA (unitless logit)")  # Table 6, Baseline CMV DNA levels Intermediate 0.694 (SE 0.364), OR 2.00, p = 0.056
    e_cmvdna_vhigh_logit <- 1.52; label("Log-odds of treatment-emergent maribavir resistance for high (>= 91 000 IU/mL) versus low baseline CMV DNA (unitless logit)")  # Table 6, Baseline CMV DNA levels High 1.52 (SE 0.632), OR 4.57, p = 0.016
    e_dis_cmv_res_bl_logit <- -1.08; label("Log-odds of treatment-emergent maribavir resistance for CMV resistant to prior anti-CMV therapy at baseline (unitless logit)")  # Table 6, Baseline CMV resistant -1.08 (SE 0.388), OR 0.341, p = 0.006
    e_cd4_pp65_mid_logit <- -1.52; label("Log-odds of treatment-emergent maribavir resistance for the middle baseline CD4+CD69+ pp65 stimulation stratum versus the lowest stratum (unitless logit)")  # Table 6, Baseline CD4+CD69+ cell count >= 0.5% to < 2%, -1.52 (SE 0.594), OR 0.219, p = 0.010
    e_cd4_pp65_hi_logit <- -0.0534; label("Log-odds of treatment-emergent maribavir resistance for the highest baseline CD4+CD69+ pp65 stimulation stratum versus the lowest stratum (unitless logit)")  # Table 6, Baseline CD4+CD69+ cell count >= 2%, -0.0534 (SE 0.768), OR 0.948, p = 0.945
    e_cd4_pp65_nr_logit <- -0.185; label("Log-odds of treatment-emergent maribavir resistance for an unreported baseline CD4+CD69+ pp65 assay versus the lowest stratum (unitless logit)")  # Table 6, Baseline CD4+CD69+ cell count Not reported, -0.185 (SE 0.428), OR 0.831, p = 0.666
    e_age_logit <- -0.0328; label("Log-odds of treatment-emergent maribavir resistance per year of age (unitless logit per year)")  # Table 6, Age -0.0328 (SE 0.0129), OR 0.968, p = 0.011
    e_sot_logit <- 0.970; label("Log-odds of treatment-emergent maribavir resistance for a solid organ transplant versus haematopoietic cell transplant recipient (unitless logit)")  # Table 6, Solid organ transplant 0.970 (SE 0.413), OR 2.64, p = 0.019

    # Bernoulli likelihood in the source; the fixed additive residual is a
    # placeholder so rxode2 has an error model, not a published quantity.
    addSd_prob_te_resist_mbv <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # Table 6 estimates the exposure coefficient per 2.5 ug/mL increment.
    ctrough25 <- CTROUGH / 2.5

    # Three-level baseline CMV DNA from the two nested indicators:
    # intermediate = CMVDNA_HIGH without CMVDNA_VHIGH.
    cmvdna_int <- CMVDNA_HIGH - CMVDNA_VHIGH

    # Solid organ transplant indicator; exact because the cohort holds only
    # SOT and HCT recipients.
    tx_sot <- 1 - TX_HCT

    logit_te_resist_mbv <- logit_ref +
      e_ctrough_logit * ctrough25 +
      e_cmvdna_int_logit * cmvdna_int +
      e_cmvdna_vhigh_logit * CMVDNA_VHIGH +
      e_dis_cmv_res_bl_logit * DIS_CMV_RES_BL +
      e_cd4_pp65_mid_logit * CD4_PP65_MID +
      e_cd4_pp65_hi_logit * CD4_PP65_HI +
      e_cd4_pp65_nr_logit * CD4_PP65_NR +
      e_age_logit * AGE +
      e_sot_logit * tx_sot

    prob_te_resist_mbv <- expit(logit_te_resist_mbv)

    prob_te_resist_mbv ~ add(addSd_prob_te_resist_mbv)
  })
}
