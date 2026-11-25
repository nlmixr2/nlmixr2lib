Zhu_2019_ipatasertib_diarrhea_grade3_nominal <- function() {
  description <- paste0(
    "Binomial logistic-regression exposure-SAFETY model for grade >= 3 ",
    "diarrhea as a function of the ipatasertib steady-state daily AUC ",
    "based on the NOMINAL dose (AUCss; the exploratory analysis of Zhu ",
    "2019 that ignores dose interruptions and reductions) in men with ",
    "metastatic castration-resistant prostate cancer receiving ",
    "ipatasertib 200 or 400 mg once daily or placebo with abiraterone ",
    "and prednisone (A.MARTIN phase II, n = 253). The probability is ",
    "expit(-4.38 + 0.000640 * AUC_IPAT_SS) with the exposure in ",
    "ng*h/mL entering linearly and uncentred, so the intercept is the ",
    "placebo log-odds. There is no PK layer and no ODE: the exposure ",
    "metric is supplied as a data column, derived in the source ",
    "analysis from the individual apparent clearance of an ipatasertib ",
    "population PK model that the paper does not print. No ",
    "between-subject random effect and no residual error are estimated ",
    "(Bernoulli likelihood). Companion fit of the same endpoint on the ",
    "actual-dose exposure metric: ",
    "Zhu_2019_ipatasertib_diarrhea_grade3."
  )
  reference <- paste(
    "Zhu R, Poland B, Wada R, Liu Q, Musib L, Maslyar D, Cho E, Yu W, Ma H,",
    "Jin JY, Budha N.",
    "Exposure-response-based product profile-driven clinical utility index",
    "for ipatasertib dose selection in prostate cancer.",
    "CPT Pharmacometrics Syst Pharmacol. 2019;8(4):240.",
    "doi:10.1002/psp4.12394.",
    "Coefficients are transcribed from Supporting Information Table S2",
    "(Table S1 for the Cox proportional-hazards models); exposure units are",
    "taken from the x-axis of Supporting Figure S3.",
    sep = " "
  )
  vignette <- "Zhu_2019_ipatasertib_exposure_response"
  units <- list(
    time = "n/a (static landmark exposure-response regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the AUC_IPAT_SS covariate column)",
    concentration = "prob_diarrhea_grade3 (probability of grade >= 3 diarrhea, 0-1; also logit_diarrhea_grade3)"
  )

  covariateData <- list(
    AUC_IPAT_SS = list(
      description = "Individual ipatasertib steady-state daily AUC based on the NOMINAL (randomised) dose, AUCss: nominal daily dose divided by the individual apparent clearance from the ipatasertib population PK model, ignoring dose interruptions and reductions. Supplied as data: this model has no PK layer.",
      units = "ng*h/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Zhu 2019 Methods, Additional exploratory E-R analysis. Enters",
        "LINEARLY and uncentred on the linear predictor; placebo patients",
        "carry 0. Because it ignores dose modifications, the fitted",
        "exposure-response is flatter than the AUC_IPAT_EVT fit of the same",
        "endpoint (Zhu 2019 Figure 3, Figure S4) and the model is used for",
        "dose projection WITHOUT a dose-intensity model.",
        "In this model the exposure enters with slope 0.000640 per ng*h/mL."
      ),
      source_name = "AUCss"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 253L,
    n_studies = 1L,
    sex_female_pct = 0,
    disease_state = "Metastatic castration-resistant prostate cancer previously treated with docetaxel chemotherapy; randomised, double-blind phase II part of the phase Ib/II A.MARTIN study (NCT01485861).",
    dose_range = "Ipatasertib 200 mg or 400 mg orally once daily, or placebo, each combined with abiraterone 1000 mg and prednisone 10 mg daily.",
    arms = c(placebo = 82L, ipatasertib_200mg = 87L, ipatasertib_400mg = 84L),
    notes = paste(
      "All 253 randomised patients provided evaluable PK, rPFS and safety",
      "data (Zhu 2019 Results). Placebo patients enter the regressions with",
      "zero ipatasertib exposure. Dose modifications, mainly for adverse",
      "events, occurred in 14.6 percent (12/82) of placebo, 14.9 percent",
      "(13/87) of 200 mg and 39.3 percent (33/84) of 400 mg patients. Age,",
      "weight and race are not reported in Zhu 2019; see the primary",
      "A.MARTIN publication (de Bono 2019, Clin Cancer Res 25:928)."
    )
  )

  ini({
    # ==================================================================
    # Zhu 2019 Table S2 row block 'Diarrhea grade >= 3, AUC ss as exposure metrics':
    #   logit(P(AE = 1)) = alpha + beta * AUC
    # (equation printed beneath Table S2). Estimates with standard error
    # and 95 percent confidence interval:
    #   alpha -4.38  SE 0.661  95% CI -5.87 to -3.24
    #   beta  0.000640  SE 0.000204  95% CI 0.000259 to 0.00107  (P 0.00173)
    # Fitted with R glm(), binomial likelihood, logit link (Zhu 2019
    # Methods, E-R analysis). The exposure unit, ng*h/mL, is the x-axis
    # label of Supporting Figure S3.
    # ==================================================================
    logit_ref <- -4.38; label("Logit of the probability of grade >= 3 diarrhea at zero ipatasertib exposure, i.e. on placebo plus abiraterone (unitless logit)") # Zhu 2019 Table S2, 'Diarrhea grade >= 3, AUC ss as exposure metrics', alpha = -4.38 (SE 0.661)
    e_auc_logit <- 0.000640; label("Change in the log-odds of grade >= 3 diarrhea per 1 ng*h/mL increase in AUC_IPAT_SS (unitless logit)") # Zhu 2019 Table S2, 'Diarrhea grade >= 3, AUC ss as exposure metrics', beta = 0.000640 (SE 0.000204, P 0.00173)

    # The source likelihood is Bernoulli, so it estimates neither a random
    # effect nor a residual variance. The small additive term below exists
    # only so rxode2 has an error model to attach to the typical-value
    # probability; it is NOT a published quantity.
    addSd_prob_diarrhea_grade3 <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_diarrhea_grade3 <- logit_ref + e_auc_logit * AUC_IPAT_SS
    prob_diarrhea_grade3 <- expit(logit_diarrhea_grade3)

    prob_diarrhea_grade3 ~ add(addSd_prob_diarrhea_grade3)
  })
}
