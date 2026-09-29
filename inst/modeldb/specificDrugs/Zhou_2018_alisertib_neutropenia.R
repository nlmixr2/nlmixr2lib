Zhou_2018_alisertib_neutropenia <- function() {
  description <- paste0(
    "Landmark logistic-regression exposure-SAFETY model for CTCAE grade 3 or ",
    "higher neutropenia as a function of the log of the time-averaged daily ",
    "AUC of the Aurora A kinase inhibitor alisertib, in adult cancer patients ",
    "given 5 to 200 mg/day on the 7-days-on, 14-days-off schedule (Zhou 2018, ",
    "n = 591). There is no PK layer and no ODE: the exposure metric is ",
    "supplied as the AUC_ALIS data column, computed in the source analysis ",
    "from each patient's actually administered doses (dose modifications ",
    "included) and their individual CL/F from the companion population PK ",
    "model packaged as Zhou_2018_alisertib, averaged from the start of dosing ",
    "to the onset of the worst grade of the toxicity. Exposure was a ",
    "significant predictor at p < 0.0001. Both coefficients are ",
    "FIGURE-DERIVED: Zhou 2018 tabulates no logistic-regression parameters, ",
    "so they were recovered by digitising the fitted curve of Figure 7 panel ",
    "A, which is exactly logit-linear in log(AUC). The recovered model ",
    "reproduces all three probabilities the paper prints in its Results text ",
    "to within the printed rounding -- see the vignette. Companion models for ",
    "the other two mechanism-related antiproliferative toxicities are ",
    "Zhou_2018_alisertib_stomatitis and Zhou_2018_alisertib_diarrhea."
  )
  reference <- paste(
    "Zhou X, Mould DR, Takubo T, Sheldon-Waniga E, Huebner D, Milton A,",
    "Venkatakrishnan K. Global population pharmacokinetics of the",
    "investigational Aurora A kinase inhibitor alisertib in cancer patients:",
    "rationale for lower dosage in Asia.",
    "Br J Clin Pharmacol. 2018;84(1):35-51. doi:10.1111/bcp.13430.",
    "Exposure-safety methods in Methods 'Exposure-safety analyses'; fitted",
    "relationship in Figure 7 panel A; model-predicted incidences quoted in",
    "Results 'Exposure-safety relationships'.",
    "Individual exposures derive from the companion population PK model; see",
    "modellib('Zhou_2018_alisertib').",
    sep = " "
  )
  vignette <- "Zhou_2018_alisertib"
  units <- list(
    time = "n/a (static landmark exposure-response regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the AUC_ALIS covariate column)",
    concentration = "prob_neutropenia_grade3 (probability of grade 3 or higher neutropenia, 0-1; also logit_neutropenia_grade3)"
  )

  covariateData <- list(
    AUC_ALIS = list(
      description = "Individual time-averaged daily area under the alisertib plasma concentration-time curve, from the start of alisertib dosing to the onset of the worst grade of the toxicity of interest. Supplied as data: this model has no PK layer.",
      units = "umol*h/L/day",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Zhou 2018 computed this per patient from the ACTUALLY ADMINISTERED",
        "doses -- dose modifications included -- and the individual empirical",
        "Bayes estimate of CL/F from the companion population PK model, so it is",
        "a treatment-course average and not a steady-state metric. The model",
        "enters log(AUC_ALIS), as Methods specify a logistic regression on the",
        "log-transformed time-averaged AUC; the Figure 7 curve is exactly",
        "logit-linear in log(AUC), confirming the form.",
        "Figure 7 panel A plots the observed range out to about 160 umol*h/L/day.",
        "Reference exposures quoted in Results: 15.63 (Western patients at 50 mg",
        "twice daily), 23.72 (East Asian patients at 50 mg twice daily) and 14.23",
        "(East Asian patients at the 30 mg twice daily regional recommended",
        "phase 2 dose).",
        "Units are micromolar-hours per day because alisertib concentrations in",
        "this analysis are molar; multiply by the molar mass to obtain mass units."
      ),
      source_name = "Time-averaged AUC (uM x h/day)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 591L,
    n_studies = 10L,
    age_range = "21 to 88 years across the pooled data set",
    age_median = "62 years",
    weight_range = "40.8 to 205.0 kg across the pooled data set",
    weight_median = "74.0 kg",
    sex_female_pct = 45,
    race_ethnicity = c(White = 83, `Asian (all)` = 9, Black = 5, Other = 1, Missing = 1),
    disease_state = "Advanced haematological and nonhaematological malignancies",
    dose_range = "5 to 200 mg/day alisertib on the 7-days-on, 14-days-off schedule",
    regions = "Western countries and the East Asian region (Japan, Singapore, Taiwan, Hong Kong, South Korea)",
    notes = paste0(
      "The exposure-safety analysis population is the subset of the pooled ",
      "alisertib clinical programme that received the 7-day dosing schedule ",
      "and was evaluable for the endpoint; n differs slightly across the three ",
      "toxicities (591 neutropenia, 593 stomatitis, 594 diarrhoea). Baseline ",
      "demographics are those of the combined 671-patient data set, Tables 2 ",
      "and 3 of Zhou 2018. Grade 3 was chosen as the neutropenia cut-off, ",
      "versus grade 2 for stomatitis and diarrhoea, on the differential ",
      "clinical relevance of each event to tolerability and quality of life."
    )
  )

  ini({
    # ================================================================
    #   logit(p) = logit_ref + e_auc_logit * log(AUC_ALIS)
    #
    # FIGURE-DERIVED PARAMETERS. Zhou 2018 publishes no coefficient
    # table for the exposure-safety logistic regressions; Figure 7
    # panel A is the only place the fitted relationship appears. Both
    # values below were recovered by digitising that curve at 400 dpi
    # from the publisher PDF (page 13) and regressing logit(p) on
    # log(AUC): the fit is essentially exact, with a residual SD of
    # 0.0024 on the logit scale over 498 extracted points spanning
    # AUC 2.2 to 160 umol*h/L/day, which is what confirms the assumed
    # logit-linear-in-log(AUC) form rather than assuming it.
    #
    # The digitisation is self-validating. Without using any printed
    # number, the recovered coefficients predict probabilities of
    # 0.387 / 0.456 / 0.372 at the three exposures Zhou 2018 quotes in
    # its Results text, against the printed 0.39 / 0.46 / 0.37 -- every
    # one inside the printed rounding. The vignette re-runs this check.
    # Treat the values as good to roughly the second decimal, not as
    # published point estimates.
    # ================================================================

    logit_ref <- -2.316; label("Logit of the grade 3 or higher neutropenia probability at a time-averaged alisertib AUC of 1 umol*h/L/day (unitless logit)")  # digitised from Zhou 2018 Figure 7 panel A; not tabulated in the source
    e_auc_logit <- 0.675; label("Change in the log-odds of grade 3 or higher neutropenia per unit increase in the natural log of the time-averaged alisertib AUC (unitless logit)")  # digitised from Zhou 2018 Figure 7 panel A; not tabulated in the source

    # The source likelihood is Bernoulli, so it estimates neither a
    # random effect nor a residual variance. The small additive term
    # below exists only so rxode2 has an error model to attach to the
    # typical-value probability; it is NOT a published quantity.
    addSd_prob_neutropenia_grade3 <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_neutropenia_grade3 <- logit_ref + e_auc_logit * log(AUC_ALIS)
    prob_neutropenia_grade3 <- expit(logit_neutropenia_grade3)

    prob_neutropenia_grade3 ~ add(addSd_prob_neutropenia_grade3)
  })
}
