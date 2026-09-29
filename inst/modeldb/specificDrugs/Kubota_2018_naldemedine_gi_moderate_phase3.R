Kubota_2018_naldemedine_gi_moderate_phase3 <- function() {
  description <- paste0(
    "Landmark binomial logistic-regression exposure-SAFETY model relating ",
    "the steady-state daily AUC of naldemedine to the probability of a treatment-emergent gastrointestinal-disorder adverse event of MODERATE OR GREATER severity ",
    "in patients with chronic non-cancer pain and opioid-induced constipation (Kubota 2018, pooled studies 1314V9231 and 1315V9232, n = 993). ",
    "The model is the paper's Eq. (3), a two-parameter logistic in AUCss: ",
    "prob = 1 / (1 + exp(-(a + b * AUCss))) with a = -2.75 and b = 0.0125. ",
    "The 95% CI of the slope includes zero; Kubota 2018 attributes this to an insufficient frequency of moderate and severe events rather than to a genuinely flat relationship, and notes that dose-dependent gastrointestinal adverse reactions had been seen in an earlier study that went up to 3 mg, 15 times the clinical dose. ",
    "There is no PK layer and no ODE: the exposure metric is supplied as the ",
    "AUC_NALD data column, which the source analysis obtained by empirical ",
    "Bayes estimation from the companion population PK model packaged as ",
    "Kubota_2018_naldemedine. No covariates were tested, because no ",
    "prognostic factors of efficacy or safety had been identified in the ",
    "contributing studies. No between-subject random effect and no residual ",
    "error are estimated (Bernoulli likelihood). Six companion ",
    "exposure-response models in the Kubota_2018_naldemedine_* family."
  )
  reference <- paste(
    "Kubota R, Fukumura K, Wajima T.",
    "Population Pharmacokinetics and Exposure-Response Relationships of Naldemedine.",
    "Pharm Res. 2018;35(11):225.",
    "doi:10.1007/s11095-018-2501-7. PMCID: PMC6182381.",
    "Logistic-regression coefficients are Table V(b), rows '1314V9231 and 1315V9232 (Phase 3)' / 'Moderate, Severe';",
    "the logistic form is Eq. (3), spelled out in the Table VI footnotes as",
    "'Probability = 1 / [1 + exp(-a - b x AUCss)]';",
    "model-predicted probabilities against the observed frequencies are",
    "Table VI and Supplemental Table S7.",
    sep = " "
  )
  vignette <- "Kubota_2018_naldemedine"
  units <- list(
    time = "n/a (static landmark exposure-response regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the AUC_NALD covariate column)",
    concentration = "prob_gi_moderate_or_worse (probability of a gastrointestinal disorder of moderate or greater severity, 0-1; also logit_gi_moderate_or_worse)"
  )

  covariateData <- list(
    AUC_NALD = list(
      description = "Individual steady-state area under the naldemedine plasma concentration-time curve over the 24 h once-daily dosing interval. Supplied as data: this model has no PK layer, and the source analysis used empirical Bayes estimates from the companion population PK model.",
      units = "ng*h/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Kubota 2018 estimates the slope per 1 ng*h/mL, so model() uses this",
        "column unscaled. Subjects randomised to PLACEBO were assigned AUCss = 0",
        "and contribute to the fit; the intercept is therefore an anchor at true",
        "zero exposure with a real placebo-arm interpretation, not an",
        "extrapolation.",
        "Observed exposures (Supplemental Table S4, empirical Bayes estimates):",
        "phase 2b 0.1 mg mean 10.50 (median 10.38, range 7.02-13.49), 0.2 mg mean",
        "22.11 (median 20.25, range 12.66-51.36), 0.4 mg mean 43.76 (median 38.18,",
        "range 25.78-74.93); phase 3 0.2 mg mean 27.50 (median 24.47, range",
        "5.669-92.1).",
        "Note that the probabilities Kubota 2018 tabulates in Table VI for the 0.1",
        "and 0.4 mg groups were NOT computed at those observed means. The Table VI",
        "footnotes state that the 0.1 mg AUC was 'assumed to be half' and the 0.4",
        "mg AUC 'double' the 0.2 mg value, giving 11.06 and 44.22 in phase 2b",
        "rather than the observed 10.50 and 43.76. Reproducing Table VI requires",
        "the assumed values; reproducing the fit requires the observed ones.",
        "A reference dose-proportional value for a typical subject can be obtained",
        "from the companion PK model as Dose / (CL/F): 0.2 mg / 9.10 L/h = 21.98",
        "ng*h/mL for the reference subject."
      ),
      source_name = "AUCss"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 993L,
    n_studies = 2L,
    n_observations = "993 binary gastrointestinal disorder of moderate or greater severity records, one per subject",
    age_range = "Not reported separately for the exposure-response subset; the parent population PK analysis population spans 18-90 years",
    weight_range = "Not reported separately for the exposure-response subset; the parent population PK analysis population spans 34.4-188.1 kg",
    sex_female_pct = NA_real_,
    disease_state = "patients with chronic non-cancer pain and opioid-induced constipation",
    dose_range = "Placebo or 0.2 mg naldemedine once daily for 12 weeks",
    regions = "Global",
    notes = paste(
      "The endpoint is the occurrence of a treatment-emergent adverse event in the gastrointestinal-disorders system organ class, dichotomised at the MODERATE-or-greater severity threshold. Observed frequencies (Table IV): placebo 28/548 (5.1%), 0.2 mg 42/445 (9.4%).",
      "The exposure-response population is NOT the population PK analysis",
      "population: it comprises the subjects of the population PK analysis OR of",
      "the placebo group who also had SAFETY data, so placebo subjects who",
      "contributed no naldemedine concentrations are included with AUCss = 0.",
      "Kubota 2018 analysed the phase 2b study and the pooled phase 3 studies",
      "SEPARATELY because the study durations and the responder definitions",
      "differ; this file packages the pooled phase 3 fit only. Its companion is the",
      "matching model for the other study set."
    )
  )

  ini({
    # ==================================================================
    # Kubota 2018 Table V(b), rows '1314V9231 and 1315V9232 (Phase 3)' / 'Moderate, Severe'. The model is the paper's Eq. (3):
    #
    #   Probability(event) = 1 / (1 + exp(-(a + b * AUCss)))
    #
    # written out in the Table VI footnotes as
    # 'Probability = 1 / [1 + exp(-a - b x AUCss)]'.
    #
    # Both coefficients below were confirmed a second time, independently
    # of Table V, by reproducing every model-predicted probability Kubota
    # 2018 tabulates for this endpoint in Table VI and Supplemental Table
    # S7 to the three decimal places printed. That check is re-run
    # mechanically in the validation vignette.
    #
    # AUCss is NOT centred, so the intercept is the logit at zero
    # exposure. Because placebo subjects entered the fit with AUCss = 0,
    # that is a real placebo-arm probability rather than an extrapolation.
    # ==================================================================

    logit_ref <- -2.75; label("Logit of the gastrointestinal disorder of moderate or greater severity probability at zero naldemedine exposure, i.e. the placebo arm (unitless logit)")  # Table V(b), rows '1314V9231 and 1315V9232 (Phase 3)' / 'Moderate, Severe', parameter a = -2.75 (95% CI -3.09 to -2.44)
    e_auc_logit <- 0.0125; label("Log-odds of a gastrointestinal disorder of moderate or greater severity per 1 ng*h/mL increase in naldemedine steady-state daily AUC (unitless logit)")  # Table V(b), rows '1314V9231 and 1315V9232 (Phase 3)' / 'Moderate, Severe', parameter b = 0.0125 (95% CI -0.00169 to 0.0259)

    # ----- No between-subject variability, no residual error -----
    # The source is a binomial logistic regression fitted in SAS 9.2; a
    # Bernoulli likelihood has no sigma, and the analysis estimates no
    # random effects. The tiny fixed additive residual below exists only
    # so rxode2 has an error model to attach to the typical-value
    # probability; it is NOT a published quantity. See the vignette's
    # Assumptions and deviations.
    addSd_prob_gi_moderate_or_worse <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # Linear predictor in the untransformed steady-state daily AUC.
    logit_gi_moderate_or_worse <- logit_ref + e_auc_logit * AUC_NALD

    prob_gi_moderate_or_worse <- expit(logit_gi_moderate_or_worse)

    # ----- Observation -----
    prob_gi_moderate_or_worse ~ add(addSd_prob_gi_moderate_or_worse)
  })
}
