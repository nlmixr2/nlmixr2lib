Sun_2025_maribavir_anemia <- function() {
  description <- paste0(
    "Binomial logistic-regression exposure-SAFETY model for ",
    "treatment-emergent anemia as a function of the maribavir AUC on the ",
    "DAY OF THE EVENT, in allogeneic hematopoietic cell transplant (HCT) ",
    "recipients with a first asymptomatic cytomegalovirus infection ",
    "treated with maribavir 400 mg orally twice daily (Sun 2025, n = 238, ",
    "the AURORA maribavir arm). Note that this is the DAY-OF-EVENT ",
    "exposure metric (AUCday), not the steady-state metric: Sun 2025's ",
    "abstract reports that STEADY-STATE exposures were not significantly ",
    "associated with any adverse event other than nausea and vomiting, ",
    "whereas all fourteen day-of-event relationships in Figure S3 are ",
    "significant. A day-of-event AUC is partly a consequence of when in ",
    "the treatment course the event happened, so it is the weaker causal ",
    "claim of the two. Enrolling region is retained and, unusually, with ",
    "POSITIVE coefficients -- both non-North-American regions reported ",
    "more treatment-emergent anemia. The two region estimates are the ",
    "subject of a source discrepancy between Table S4 and Figure S3 that ",
    "is resolved in favour of the figure; see the covariate notes and the ",
    "vignette. The exposure slope is +0.0669 per 10 ug*h/mL. The output ",
    "is named prob_anemia_ae rather than prob_anemia because the ",
    "canonical prob_anemia is already registered for a LABORATORY-derived ",
    "CTCAE grade >= 3 anemia endpoint, whereas this one is an ungraded ",
    "investigator-reported adverse-event term and therefore sits at a ",
    "different event rate; the two must not be conflated. There is no PK ",
    "layer and no ODE: the exposure metric is supplied as a data column, ",
    "derived in the source analysis from the companion population PK ",
    "model packaged as Sun_2025_maribavir. No between-subject random ",
    "effect and no residual error are estimated (Bernoulli likelihood). ",
    "Fifteen companion exposure-response models in the ",
    "Sun_2025_maribavir_* family."
  )
  reference <- paste(
    "Sun K, Jomphe C, Gosselin NH, Pheng L, Durairaj C, Hang Y, Bhattacharya I.",
    "Population Pharmacokinetics and Exposure-Response Relationships of Maribavir",
    "in Transplant Recipients With First Episode or Refractory Cytomegalovirus.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(8):1346-1356.",
    "doi:10.1002/psp4.70054.",
    "Logistic-regression coefficients in the Figure S3 parameter tables, Supporting Information file s002; the corresponding odds ratios are also tabulated in Table S4, file s001.",
    sep = " "
  )
  vignette <- "Sun_2025_maribavir"
  units <- list(
    time          = "n/a (static landmark exposure-response regression; no time dimension)",
    dosing        = "n/a (no dose events; exposure enters as the AUC_MBV_DAY covariate column)",
    concentration = "prob_anemia_ae (probability of treatment-emergent anemia, 0-1; also logit_anemia_ae)"
  )

  covariateData <- list(
    AUC_MBV_DAY = list(
      description        = "Individual maribavir area under the plasma concentration-time curve over the calendar day on which the adverse event occurred. Supplied as data: this model has no PK layer, and the source analysis used individual predictions from the companion maribavir population PK model.",
      units              = "ug*h/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Sun 2025 estimates the coefficient PER INCREMENT OF 10 ug*h/mL, so",
        "model() divides this column by 10. This is the DAY-OF-EVENT exposure, not",
        "the steady-state exposure: where a subject had more than one event of the",
        "same type, the source retained the exposure on the day of the MOST SEVERE",
        "event (Methods 2.4), and subjects without the event contribute the",
        "exposure on the corresponding on-treatment day. Do not substitute",
        "AUC_MBV_SS -- that column drives the two EFFICACY models and is a",
        "different quantity. The Figure S3 quartile boundaries run from about",
        "[0, 116) to [313, 951] ug*h/mL, so the observed range spans roughly an",
        "order of magnitude.",
        "In this model the exposure enters with slope 0.0669 per 10 ug*h/mL (odds",
        "ratio 1.07 (1.04-1.10), p < 0.001)."
      ),
      source_name        = "AUCday of maribavir (increment of 10 h.ug/mL)"
    ),
    REGION_ASIAPACIFIC = list(
      description        = "Asia Pacific enrolling-region indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (North America, the reference region, when REGION_EUROPE is also 0)",
      notes              = paste(
        "41 of 238 subjects (17.2%). Paired with REGION_EUROPE against a NORTH",
        "AMERICA reference. Broader than REGION_JAPAN or REGION_EASTASIA: Asia",
        "Pacific as a trial-operations region conventionally spans East Asia,",
        "South-East Asia and Oceania.",
        "Coefficient +1.14, odds ratio 3.12. SOURCE DISCREPANCY: Table S4 assigns",
        "odds ratio 3.12 to EUROPE and 4.62 to Asia Pacific, the reverse of Figure",
        "S3. The figure is preferred because it is internally consistent --",
        "exp(1.14) = 3.13 and exp(1.53) = 4.62 reproduce its own odds ratios --",
        "whereas Table S4 prints odds ratios only and cannot be checked against",
        "itself. Anemia is the only endpoint where the two sources conflict."
      ),
      source_name        = "Enrolling region, Asia Pacific"
    ),
    REGION_EUROPE = list(
      description        = "Europe enrolling-region indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (North America, the reference region, when REGION_ASIAPACIFIC is also 0)",
      notes              = paste(
        "138 of 238 subjects (58.0%). Paired with REGION_ASIAPACIFIC against a",
        "NORTH AMERICA reference -- both indicators 0 -- rather than the",
        "REGION_* family's more usual US reference. Enrolling region here is a",
        "proxy for regional differences in transplant and CMV-management practice,",
        "not for ancestry: Table S2 records 138 European enrollments but only 31",
        "subjects of Asian race against 41 Asia Pacific enrollments, so region and",
        "race disagree within this cohort.",
        "Coefficient +1.53, odds ratio 4.62. See the REGION_ASIAPACIFIC note for",
        "the Table S4 / Figure S3 discrepancy that this pair of values resolves."
      ),
      source_name        = "Enrolling region, Europe"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 238L,
    n_studies      = 1L,
    n_observations = "238 binary anemia records, one per patient",
    age_range      = "12 to <18 years: 1 (0.4%); 18 to <45: 47 (19.7%); 45 to <65: 142 (59.7%); >=65: 48 (20.2%) (Table S2)",
    weight_range   = "Not reported for the AURORA arm; Table 1 reports 'NR' for the AURORA weights specifically. The pooled PK analysis population spans 36.1-141 kg.",
    sex_female_pct = 45.4,
    race_ethnicity = c(Caucasian = 79.8, Asian = 13.0, Black = 3.4, Other = 2.9, Missing = 0.8),
    disease_state  = "Allogeneic HCT recipients (238/238, 100%) with FIRST asymptomatic cytomegalovirus infection after transplant; all were maribavir-susceptible at baseline. Baseline CMV DNA: very low 71, low 127, high 39, missing 1. CMV serostatus D+/R+ 123 (51.7%), D-/R+ 83 (34.9%), D+/R- 18 (7.6%), D-/R- 8 (3.4%). Reason for transplant: acute myeloid leukaemia 88 (37.0%), myelodysplastic syndrome 37 (15.5%), acute lymphocytic leukaemia 21 (8.8%), non-Hodgkin lymphoma 19 (8.0%), other 73 (30.7%). Conditioning: reduced-intensity 116 (48.7%), myeloablative 86 (36.1%), non-myeloablative 32 (13.4%).",
    dose_range     = "Maribavir 400 mg orally twice daily, the AURORA randomized dose",
    regions        = "North America 59 (24.8%), Europe 138 (58.0%), Asia Pacific 41 (17.2%) (Table S2)",
    notes          = paste0(
      "This is the exposure-response analysis population: the maribavir arm of ",
      "the phase 3 AURORA study, a subset of the 930-subject population PK ",
      "analysis population. Individual exposures were derived from the ",
      "companion population PK model, packaged as Sun_2025_maribavir. ",
      "Sun 2025 selected safety endpoints for formal exposure-response analysis by incidence (TEAEs occurring in more than 10% of participants, plus the most severe TEAEs)."
    )
  )

  ini({
    # ================================================================
    # Sun 2025 Figure S3 panel j (Supporting Information file s002).
    #
    #   logit(p) = intercept + slope * (AUC_MBV_DAY / 10)
    #              + risk-factor terms
    #
    # Every value below is the source's 'Estimate (SE)' column. Each is
    # confirmed a second time by its printed odds ratio, which the source
    # defines as the exponentiated estimate; because both columns are
    # rounded, the check is the INTERVAL one (does a real b exist with
    # round(b, 3) == Estimate and round(exp(b), 3) == OR) rather than
    # exp(printed) == printed. The checks are re-run mechanically in the
    # validation vignette. Covariates are not centred except where noted,
    # so the intercept is an anchor at zero exposure, not a reference-
    # patient probability.
    # ================================================================

    # ----- Logit intercept -----
    logit_ref <- -3.84; label("Logit of the anemia probability at zero maribavir exposure with all risk factors at their reference level (unitless logit)")  # Sun 2025 Figure S3 panel j (Supporting Information file s002), Intercept -3.84 (SE 0.590), p < 0.001

    # ----- Exposure effect -----
    e_auc_logit <- 0.0669; label("Log-odds of anemia per 10 ug*h/mL increase in maribavir AUC on the day of the event (unitless logit)")  # Sun 2025 Figure S3 panel j (Supporting Information file s002), Estimate 0.0669 (SE 0.0125), odds ratio 1.07 (1.04-1.10), p < 0.001

    # ----- Risk-factor effects on the logit -----
    e_region_asiapacific_logit <- 1.14; label("Log-odds of anemia for enrolment in Asia Pacific versus North America (unitless logit)")  # Sun 2025 Figure S3 panel j (Supporting Information file s002), Enrolling region, Asia Pacific Estimate 1.14 (SE 0.518), odds ratio 3.12 (1.13-8.61), p = 0.02826
    e_region_europe_logit <- 1.53; label("Log-odds of anemia for enrolment in Europe versus North America (unitless logit)")  # Sun 2025 Figure S3 panel j (Supporting Information file s002), Enrolling region, Europe Estimate 1.53 (SE 0.594), odds ratio 4.62 (1.44-14.8), p = 0.00991

    # ----- No between-subject variability, no residual error -----
    # The source is a binomial logistic regression; a Bernoulli likelihood
    # has no sigma, and the analysis estimates no random effects. The tiny
    # fixed additive residual below exists only so rxode2 has an error
    # model to attach to the typical-value probability; it is NOT a
    # published quantity. See the vignette's Assumptions and deviations.
    addSd_prob_anemia_ae <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # Sun 2025 estimates every coefficient per 10 ug*h/mL increment of the
    # DAY-OF-EVENT AUC, so the exposure column is divided by 10 here.
    auc10 <- AUC_MBV_DAY / 10

    # ----- Linear predictor -----
    logit_anemia_ae <- logit_ref +
      e_auc_logit * auc10 +
      e_region_asiapacific_logit * REGION_ASIAPACIFIC +
      e_region_europe_logit * REGION_EUROPE

    prob_anemia_ae <- expit(logit_anemia_ae)

    # ----- Observation -----
    prob_anemia_ae ~ add(addSd_prob_anemia_ae)
  })
}
