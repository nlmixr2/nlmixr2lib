Shemesh_2020_atezolizumab_aeg34_auc <- function() {
  description <- paste(
    "Binomial logistic-regression exposure-safety model for GRADE 3/4",
    "ADVERSE EVENTS versus Cycle-1 AUC (AUC_ATEZO, ug*day/mL) in patients",
    "with high tissue tumour mutational burden (tTMB >= 16 mutations/Mb)",
    "treated with single-agent intravenous atezolizumab, pooled from",
    "seven trials (Shemesh 2020, n = 167, Figure 3C). The probability is",
    "expit(0.4644 + -0.06491 * (AUC_ATEZO / 1000)). THE EXPOSURE TERM IS",
    "NOT STATISTICALLY SIGNIFICANT (exploratory P = .812); the flat",
    "relationship is the paper's result, not a transcription gap. Shemesh",
    "2020 prints no coefficients: both were recovered from the plotted",
    "fitted curve. There is no PK layer and no ODE; the exposure metric",
    "is supplied as a data column, in the source an empirical-Bayes",
    "prediction from the Stroh 2017 atezolizumab popPK model (packaged as",
    "the fixed intravenous layer of Chan_2025_atezolizumab). One of nine",
    "companion models in the Shemesh_2020_atezolizumab_* family (3",
    "endpoints x 3 Cycle-1 exposure metrics)."
  )
  reference <- paste(
    "Shemesh CS, Chan P, Legrand FA, Shames DS, Das Thakur M, Shi J,",
    "Bailey L, Vadhavkar S, He X, Zhang W, Bruno R. Pan-cancer population",
    "pharmacokinetics and exposure-safety and -efficacy analyses of",
    "atezolizumab in patients with high tumor mutational burden.",
    "Pharmacol Res Perspect. 2020;8(6):e00685. doi:10.1002/prp2.685.",
    "Coefficients recovered from Figure 3C."
  )
  vignette <- "Shemesh_2020_atezolizumab_tmb"
  units <- list(
    time = "n/a (static Cycle-1 landmark regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the AUC_ATEZO data column)",
    concentration = "prob_aeg34 (probability of a grade 3/4 adverse event, 0-1)"
  )

  covariateData <- list(
    AUC_ATEZO = list(
      description = "Model-predicted Cycle-1 atezolizumab AUC over the first 21-day cycle",
      units = "ug*day/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters LINEARLY and UNCENTRED; the model() block divides by 1000",
        "so the slope is per 1000 ug*day/mL, the rescaling used by the",
        "sibling Chan_2025_atezolizumab_* models. Derived in the source",
        "by Bayesian post hoc estimation (NONMEM MAXEVAL = 0) of",
        "individual parameters under the nominal regimen, using the Stroh",
        "2017 phase I popPK model without re-estimation (the same",
        "parameter set is packaged as the fixed intravenous layer of",
        "modellib('Chan_2025_atezolizumab')). tTMB-high distribution",
        "(Shemesh 2020 Table 2): geometric mean 2722 ug*day/mL (CV 67.2",
        "percent); plotted range about 290-5140 ug*day/mL (Figure 3). The",
        "linear-in-AUC form (not log-AUC) was established from the",
        "curvature of the plotted Figure 3D curve; see the vignette."
      ),
      source_name = "Cycle 1 AUC"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 167,
    n_studies = 7,
    age_range = "37-89 years (tTMB-high, Table 1)",
    weight_range = "35.4-149 kg (tTMB-high, Table 1)",
    sex_female_pct = 27.4,
    race_ethnicity = "82.9 percent White (tTMB-high, Table 1)",
    disease_state = "Locally advanced or metastatic solid tumours with high tissue tumour mutational burden (>= 16 mutations/Mb), predominantly NSCLC and urothelial carcinoma",
    dose_range = "Atezolizumab 1200 mg IV every 3 weeks (about 94 percent of patients); some PCD4989g patients received 10, 15 or 20 mg/kg IV every 3 weeks",
    regions = "Multinational",
    notes = paste(
      "tTMB-high subgroup (tTMB >= 16 mutations/Mb, FoundationOne assay)",
      "of the tTMB-evaluable population (986 of 2894 treated patients)",
      "pooled from OAK, POPLAR, BIRCH, FIR, IMvigor210, IMvigor211 and",
      "PCD4989g (Shemesh 2020 Table S1). Baseline, tTMB-high column of",
      "Table 1 (n = 175): median age 65 (37-89) years, median weight 75",
      "(35.4-149) kg, 27.4 percent female, 82.9 percent White, median",
      "albumin 39 g/L, median baseline SLD 58 mm, ADA-positive 29.3",
      "percent. Tumour types (Table S2): NSCLC 83, urothelial 70,",
      "melanoma 12, other 10. Observed endpoint rate: grade 3/4",
      "adverse-event incidence 56.9 percent (Results 3.5, 167 patients).",
      "Cycle-1 exposure metrics were used deliberately 'to avoid",
      "confounding factors on exposure, such as response-dependent",
      "time-varying clearance'."
    )
  )

  ini({
    # ==================================================================
    # Source: Figure 3C of Shemesh 2020, the black 'model-fitted curve'
    # of a logistic regression with the Cycle-1 exposure metric as the
    # only predictor. No coefficient table exists anywhere in the paper
    # or its supplement, so both coefficients are recovered from the
    # drawing.
    #
    # Figure 3C of Shemesh 2020 is embedded in the article PDF as vector
    # graphics, so the fitted curve's end points were read exactly from
    # the drawing coordinates (axis calibration from the tick marks) and
    # the two logistic coefficients solved from them; nothing is
    # eyeballed.
    #
    # End points of the plotted fitted curve (drawn from the lowest to
    # the highest observed exposure): AUC_ATEZO = 336 ug*day/mL ->
    # probability 0.6089; AUC_ATEZO = 5141.3 ug*day/mL -> probability
    # 0.5326.
    #
    # The logit is LINEAR in the untransformed exposure (not in its
    # logarithm), the form that reproduces the drawn curvature of all
    # nine Shemesh 2020 exposure-response panels. The intercept is the
    # logit at AUC_ATEZO = 0, outside the observed range; only
    # differences in the linear predictor across the plotted range are
    # supported by the source.
    # ==================================================================
    logit_ref <- 0.4644; label("Logit of the probability of a grade 3/4 adverse event at AUC_ATEZO = 0 ug*day/mL (unitless logit)")  # digitised from Figure 3C (exact curve end points; see comment block above)
    e_auc_atezo_logit <- -0.06491; label("Log-odds of a grade 3/4 adverse event per 1000 ug*day/mL increase in Cycle-1 AUC (unitless logit)")  # digitised from Figure 3C; printed exploratory P = .812 (not significant)

    # No between-subject variability and no residual error: the source is
    # a Bernoulli-likelihood logistic regression. The placeholder additive
    # term exists only because rxode2 requires an observation declaration;
    # see the vignette Assumptions and deviations.
    addSd_prob_aeg34 <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # Linear predictor, linear and uncentred in the exposure metric. The
    # divisor 1000 converts the ug*day/mL data column into the per-1000
    # unit in which the slope is expressed.
    logit_aeg34 <- logit_ref + e_auc_atezo_logit * (AUC_ATEZO / 1000)

    prob_aeg34 <- expit(logit_aeg34)

    # Deterministic probability of the event. Downstream callers can
    # sample binary outcomes with rbinom(n, 1, prob_aeg34) on the rxSolve
    # output.
    prob_aeg34 ~ add(addSd_prob_aeg34)
  })
}
