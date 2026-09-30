Shemesh_2020_atezolizumab_orr_ctrough <- function() {
  description <- paste(
    "Binomial logistic-regression exposure-efficacy model for OBJECTIVE",
    "RESPONSE (complete or partial response, investigator-assessed RECIST",
    "1.1) versus Cycle-1 Cmin (CTROUGH, ug/mL) in patients with high",
    "tissue tumour mutational burden (tTMB >= 16 mutations/Mb) treated",
    "with single-agent intravenous atezolizumab, pooled from seven trials",
    "(Shemesh 2020, n = 171, Figure 3B). The probability is expit(-0.866",
    "+ 0.0 * (CTROUGH / 10)). THE EXPOSURE TERM IS NOT STATISTICALLY",
    "SIGNIFICANT (exploratory P = .998); the flat relationship is the",
    "paper's result, not a transcription gap. Shemesh 2020 prints no",
    "coefficients: both were recovered from the plotted fitted curve.",
    "There is no PK layer and no ODE; the exposure metric is supplied as",
    "a data column, in the source an empirical-Bayes prediction from the",
    "Stroh 2017 atezolizumab popPK model (packaged as the fixed",
    "intravenous layer of Chan_2025_atezolizumab). One of nine companion",
    "models in the Shemesh_2020_atezolizumab_* family (3 endpoints x 3",
    "Cycle-1 exposure metrics)."
  )
  reference <- paste(
    "Shemesh CS, Chan P, Legrand FA, Shames DS, Das Thakur M, Shi J,",
    "Bailey L, Vadhavkar S, He X, Zhang W, Bruno R. Pan-cancer population",
    "pharmacokinetics and exposure-safety and -efficacy analyses of",
    "atezolizumab in patients with high tumor mutational burden.",
    "Pharmacol Res Perspect. 2020;8(6):e00685. doi:10.1002/prp2.685.",
    "Coefficients recovered from Figure 3B."
  )
  vignette <- "Shemesh_2020_atezolizumab_tmb"
  units <- list(
    time = "n/a (static Cycle-1 landmark regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CTROUGH data column)",
    concentration = "prob_orr_investigator (probability of objective response (CR + PR, investigator-assessed RECIST 1.1), 0-1)"
  )

  covariateData <- list(
    CTROUGH = list(
      description = "Model-predicted Cycle-1 minimum (trough) serum atezolizumab concentration",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "TOTAL serum atezolizumab at the END OF CYCLE 1: Shemesh 2020",
        "Methods, 'cycle 1 Cmin was based on day 1 of cycle 2 pre-dose",
        "samples', i.e. the day-21 trough after a single 1200 mg dose.",
        "Source name 'Cmin', carried by the CTROUGH canonical (source",
        "alias). Enters LINEARLY and UNCENTRED; the model() block divides",
        "by 10 so the slope is per 10 ug/mL, as in",
        "Chan_2025_atezolizumab_pfs. Empirical-Bayes prediction from the",
        "Stroh 2017 phase I popPK model (packaged as the fixed",
        "intravenous layer of modellib('Chan_2025_atezolizumab')).",
        "tTMB-high distribution (Shemesh 2020 Table 2): geometric mean",
        "63.2 ug/mL (CV 200 percent); plotted range about 1-135 ug/mL",
        "(Figure 3B, Figure S6)."
      ),
      source_name = "Cycle 1 Cmin"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 171,
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
      "melanoma 12, other 10. Observed endpoint rate: ORR 29.7 percent",
      "(52 of 175 tTMB-high patients; 171 in the exposure-response set).",
      "Cycle-1 exposure metrics were used deliberately 'to avoid",
      "confounding factors on exposure, such as response-dependent",
      "time-varying clearance'."
    )
  )

  ini({
    # ==================================================================
    # Source: Figure 3B of Shemesh 2020, the black 'model-fitted curve'
    # of a logistic regression with the Cycle-1 exposure metric as the
    # only predictor. No coefficient table exists anywhere in the paper
    # or its supplement, so both coefficients are recovered from the
    # drawing.
    #
    # Figure 3B of Shemesh 2020 is embedded in the article PDF as vector
    # graphics, so the fitted curve's end points were read exactly from
    # the drawing coordinates (axis calibration from the tick marks) and
    # the two logistic coefficients solved from them; nothing is
    # eyeballed.
    #
    # End points of the plotted fitted curve (drawn from the lowest to
    # the highest observed exposure): CTROUGH = 6.1 ug/mL -> probability
    # 0.2961; CTROUGH = 129.4 ug/mL -> probability 0.2961.
    #
    # The logit is LINEAR in the untransformed exposure (not in its
    # logarithm), the form that reproduces the drawn curvature of all
    # nine Shemesh 2020 exposure-response panels. The intercept is the
    # logit at CTROUGH = 0, outside the observed range; only differences
    # in the linear predictor across the plotted range are supported by
    # the source.
    # ==================================================================
    logit_ref <- -0.866; label("Logit of the probability of objective response (CR + PR, investigator-assessed RECIST 1.1) at CTROUGH = 0 ug/mL (unitless logit)")  # digitised from Figure 3B (exact curve end points; see comment block above)
    e_ctrough_logit <- 0.0; label("Log-odds of objective response (CR + PR, investigator-assessed RECIST 1.1) per 10 ug/mL increase in Cycle-1 Cmin (unitless logit)")  # drawn exactly horizontal; slope 0 to drawing resolution, printed exploratory P = .998

    # No between-subject variability and no residual error: the source is
    # a Bernoulli-likelihood logistic regression. The placeholder additive
    # term exists only because rxode2 requires an observation declaration;
    # see the vignette Assumptions and deviations.
    addSd_prob_orr_investigator <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # Linear predictor, linear and uncentred in the exposure metric. The
    # divisor 10 converts the ug/mL data column into the per-10
    # unit in which the slope is expressed.
    logit_orr <- logit_ref + e_ctrough_logit * (CTROUGH / 10)

    prob_orr_investigator <- expit(logit_orr)

    # Deterministic probability of the event. Downstream callers can
    # sample binary outcomes with rbinom(n, 1, prob_orr_investigator) on the rxSolve
    # output.
    prob_orr_investigator ~ add(addSd_prob_orr_investigator)
  })
}
