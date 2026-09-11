Chen_2025_hemoporfin_significant_improvement <- function() {
  description <- paste0(
    "Binomial logistic exposure-response model for the probability that a ",
    "blinded expert panel grades the response to hemoporfin photodynamic ",
    "therapy as category I or II (significant improvement, at least 60% ",
    "improvement in port-wine-stain colour blanching) 8 weeks after ",
    "treatment, in 24 Chinese pediatric patients aged 7-14 years (Chen ",
    "2025, NCT03125057). The probability is ",
    "expit(-2.09 + 0.169 * AUC_HEMO) where AUC_HEMO is the individual ",
    "hemoporfin AUC(0-30min) in h*ug/mL derived from the companion ",
    "population PK model (Chen 2025 Table S1). This endpoint did NOT reach ",
    "statistical significance (slope p = 0.505, 95% CI -0.322 to 0.699) and ",
    "the non-significance is the paper's reported result, not a ",
    "transcription gap: of the five endpoints screened, only the ",
    "investigator- and patient-rating endpoints were significant. The ",
    "slope is positive but its confidence interval spans zero. There is no ",
    "PK layer and no ODE -- exposure enters as a covariate column. No ",
    "between-subject random effect and no residual error are estimated ",
    "(Bernoulli likelihood). Companion models: Chen_2025_hemoporfin ",
    "(population PK) and four sibling endpoints in the ",
    "Chen_2025_hemoporfin_* family."
  )
  reference <- paste(
    "Chen R, Zhang B, Tao J, Yao Q, Zhou T, Ma L, Xu Z.",
    "Population Pharmacokinetics and Exposure-Response Relationship of",
    "Hemoporfin in Pediatric Patients With Port-Wine Stain.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(8):1449-1457.",
    "doi:10.1002/psp4.70050. PMCID: PMC12439283.",
    "Coefficients from Table S1 (Data S1); model form from Equation 3.",
    "Exposure driver supplied by modellib('Chen_2025_hemoporfin').",
    sep = " "
  )
  vignette <- "Chen_2025_hemoporfin"
  depends <- c("Chen_2025_hemoporfin")

  units <- list(
    time          = "n/a (static landmark exposure-response regression evaluated 8 weeks after therapy; no time dimension)",
    dosing        = "n/a (no dose events; hemoporfin exposure enters as the covariate AUC_HEMO)",
    concentration = "prob_significant_improvement (probability of expert-panel category I or II, 0-1; also logit_significant_improvement)"
  )

  covariateData <- list(
    AUC_HEMO = list(
      description        = "Individual hemoporfin AUC from 0 to 30 minutes after the start of the intravenous infusion",
      units              = "h*ug/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Chen 2025 Methods 2.8: AUC(0-30min) and Cmax were chosen as the ",
        "efficacy-related exposure metrics because the photodynamic ",
        "treatment window (light applied 10-30 minutes after the start of ",
        "the 20-minute infusion) overlaps the early disposition phase. ",
        "Results 3.5: AUC(0-30min) and Cmax were highly linearly ",
        "correlated, so only AUC(0-30min) was carried into the ",
        "exposure-response analysis. Values are empirical Bayes estimates ",
        "from the companion population PK model, not observed NCA. The ",
        "mean value in this cohort at the empirical 5 mg/kg dose is ",
        "12.4 h*ug/mL (Table S2). NOT centred and NOT scaled -- the ",
        "intercept is an extrapolated anchor at AUC_HEMO = 0, not a ",
        "reference-patient probability."
      ),
      source_name        = "AUC 0-30min"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 24L,
    n_studies      = 1L,
    age_range      = "7-13 years",
    age_median     = "9 years",
    weight_range   = "21-72 kg",
    weight_median  = "32 kg",
    height_range   = "120-164 cm",
    height_median  = "139 cm",
    sex_female_pct = 41.7,
    race_ethnicity = c(Asian = 100),
    disease_state  = "Port-wine stain (congenital capillary malformation) of the head and neck.",
    dose_range     = "Single 5 mg/kg intravenous hemoporfin infusion over 20 minutes, followed by 530 nm LED irradiation at 60 mW/cm^2 (4 patients) or 75 mW/cm^2 (20 patients) from 10 to 30 minutes after the start of infusion.",
    administration = "Intravenous infusion (20 minutes) plus photodynamic therapy",
    regions        = "China",
    notes          = paste0(
      "Chen 2025 Table 1 (pediatric column) and Methods 2.5. Efficacy was ",
      "graded 8 weeks after therapy by a blinded expert panel of one ",
      "dermatologist and two plastic surgeons, from before / after photos ",
      "of the target lesion, on a four-category improvement scale: ",
      "I >= 90%, II 60-90%, III 20-60%, IV < 20%. The modelled event here ",
      "is category I or II. Only the pediatric trial contributed efficacy ",
      "data -- the adult phase I trial was PK-only and involved no ",
      "photodynamic therapy."
    )
  )

  ini({
    # ==================================================================
    # Chen 2025 Table S1 (Data S1), row block "Probability of achieving
    # significant improvement (category I+II)". Model form is Equation 3:
    #   log(P / (1 - P)) = beta0 + beta1 * Exposure metric
    # Fitted by logistic regression in R 3.4.2 (Methods 2.9), separately
    # from the NONMEM population PK run.
    #
    # NEITHER coefficient is significant at the 5% level; the slope's
    # confidence interval spans zero. Retained for completeness because
    # Table S1 prints point estimates. See the vignette's Assumptions and
    # deviations before using this model for prediction.
    # ==================================================================

    # ----- Logit intercept -----
    logit_ref <- -2.09 ; label("Logit of the probability of expert-panel category I or II at AUC_HEMO = 0 h*ug/mL (unitless logit)")  # Table S1, beta0 = -2.09 (95% CI -8.69 to 4.05, p = 0.509; not significant)

    # ----- Exposure effect on the logit -----
    e_auc_hemo_logit <- 0.169 ; label("Log-odds of expert-panel category I or II per 1 h*ug/mL increase in hemoporfin AUC(0-30min) (mL/(ug*h))")  # Table S1, beta1 = 0.169 mL*ug^-1*h^-1 (95% CI -0.322 to 0.699, p = 0.505; not significant)

    # ----- No between-subject variability, no residual error -----
    # The source is a binomial logistic regression; a Bernoulli likelihood
    # has no sigma, and no random effects were estimated. The tiny fixed
    # additive residual below exists only so rxode2 has an error model to
    # attach to the typical-value probability; it is NOT a published
    # quantity. See the vignette's Assumptions and deviations.
    addSd_prob_significant_improvement <- fixed(0.001) ; label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # ----- Linear predictor (Chen 2025 Equation 3, Table S1) -----
    logit_significant_improvement <- logit_ref + e_auc_hemo_logit * AUC_HEMO

    prob_significant_improvement <- expit(logit_significant_improvement)

    # ----- Observation -----
    prob_significant_improvement ~ add(addSd_prob_significant_improvement)
  })
}
