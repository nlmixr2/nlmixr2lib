Chen_2021_lorlatinib_icorr <- function() {
  description <- paste0(
    "Binomial logistic-regression efficacy model for INTRACRANIAL ",
    "objective response rate (IC-ORR) by independent central review in ",
    "adults with advanced ALK-positive or ROS1-positive non-small cell ",
    "lung cancer, baseline CNS metastasis and at least one prior ALK ",
    "inhibitor, treated with the third-generation ALK/ROS1 tyrosine ",
    "kinase inhibitor lorlatinib (Chen 2021, n = 132, the ",
    "CNS-metastatic subset of the phase I/II study B7461001, ",
    "NCT01970865, phase II expansion cohorts 2-5 at 100 mg orally once ",
    "daily). The probability of intracranial response is ",
    "expit(3.929 - 1.015 * log(ALP) + 0.015 * AMYL) where ALP is ",
    "baseline alkaline phosphatase in U/L and AMYL is baseline serum ",
    "amylase in U/L (Chen 2021 Eq. 5, Table 4). The log is a natural ",
    "log. NO LORLATINIB EXPOSURE TERM IS PRESENT and that absence is ",
    "the paper's headline efficacy result, not a transcription gap: ",
    "the best exposure predictor screened for this endpoint, ",
    "log(Ctrough,P1), was not statistically significant, so lorlatinib ",
    "plasma exposure does not predict intracranial response at the ",
    "doses studied. The two retained covariates act in opposite ",
    "directions -- higher baseline alkaline phosphatase (a marker of ",
    "greater disease burden) lowers the odds of response (odds ratio ",
    "0.363 per e-fold), while higher baseline amylase raises them ",
    "slightly (odds ratio 1.015 per U/L). There is no PK layer and no ",
    "ODE. No between-subject random effect and no residual error are ",
    "estimated (Bernoulli likelihood). Four companion ",
    "exposure-response models in the Chen_2021_lorlatinib_* family."
  )
  reference <- paste(
    "Chen J, Ruiz-Garcia A, James LP, Peltz G, Thurm H, Clancy J, Hibma J.",
    "Lorlatinib exposure-response analyses for safety and efficacy in a",
    "phase I/II trial to support benefit-risk assessment in non-small cell",
    "lung cancer.",
    "Clin Pharmacol Ther. 2021;110(5):1273-1281.",
    "doi:10.1002/cpt.2228.",
    sep = " "
  )
  vignette <- "Chen_2021_lorlatinib_exposure_response"
  units <- list(
    time          = "n/a (static landmark exposure-efficacy regression; no time dimension and no treatment-duration covariate)",
    dosing        = "n/a (no dose events, and no exposure covariate: this endpoint retained no lorlatinib exposure metric)",
    concentration = "prob_icorr (probability of intracranial objective response, 0-1; also logit_icorr)"
  )

  covariateData <- list(
    ALP = list(
      description        = "Baseline serum alkaline phosphatase (BAP).",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Chen 2021 Table 4 writes the unit as IU/L; U/L and IU/L are used",
        "interchangeably for this assay and the canonical spelling is U/L.",
        "Enters on the NATURAL LOG scale and is NOT normalised to any",
        "reference value -- Eq. 5 is a bare -1.015 * Log(BAP), so the",
        "intercept 3.929 is the logit at ALP = 1 U/L and AMYL = 0 U/L.",
        "The log base is not stated in the paper; natural log is",
        "established by reproduction rather than assumed. At the IC-ORR",
        "population medians (ALP 100.50 U/L, AMYL 62.00 U/L) the natural",
        "log gives a predicted probability of 0.545, consistent with the",
        "63% intracranial response rate reported for lorlatinib in",
        "ALK-TKI-pretreated patients, whereas a base-10 log gives an",
        "implausible 0.944. Efficacy-IC-ORR analysis set (N = 132) median",
        "100.50 U/L, range 13.00-1,552.00, mean 142.66 (SD 153.63)",
        "(Chen 2021 Table 2). The strongly right-skewed distribution is",
        "why the paper log-transforms this covariate and not amylase. The",
        "negative coefficient is read by Chen 2021 as reflecting greater",
        "underlying disease burden rather than a pharmacological effect."
      ),
      source_name        = "BAP (baseline alkaline phosphatase)"
    ),
    AMYL = list(
      description        = "Baseline serum amylase (BAMY).",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Chen 2021 Table 4 writes the unit as IU/L; U/L and IU/L are used",
        "interchangeably for this assay and the canonical spelling is U/L.",
        "Enters LINEARLY and untransformed, in contrast to the",
        "log-transformed ALP term in the same equation. NOT centred, so",
        "the odds ratio 1.015 is per 1 U/L. Efficacy-IC-ORR analysis set",
        "(N = 132) median 62.00 U/L, range 13.00-218.00, mean 70.32",
        "(SD 35.86), 4 patients missing (Chen 2021 Table 2). Amylase is",
        "read here as a routine clinical-chemistry marker retained by the",
        "backward-elimination algorithm, not as a mechanistic driver of",
        "intracranial response; Chen 2021 offers no biological",
        "interpretation for it and the effect is small (a full",
        "interquartile shift moves the logit by well under one unit)."
      ),
      source_name        = "BAMY (baseline amylase)"
    )
  )

  covariatesDataExcluded <- list(
    CTROUGH = list(
      description = "Individual lorlatinib trough plasma concentration over cycle 1 (Ctrough,P1).",
      units       = "ng/mL",
      type        = "continuous",
      notes       = paste(
        "The best-fitting lorlatinib exposure metric for this endpoint in",
        "the univariate forward screen -- Chen 2021 Results states 'the",
        "best lorlatinib exposure predictor to be evaluated in the final",
        "model was ... the logarithmic value of Ctrough,P1 for the IC-ORR",
        "analysis' -- but it was NOT statistically significant and does",
        "not appear in Table 4: 'in the E-R efficacy analysis for IC-ORR,",
        "lorlatinib exposure using the metric Ctrough,P1, was not a",
        "statistically significant predictor of achieving IC-ORR.' No",
        "point estimate exists because no exposure term entered the final",
        "model -- this is a published null result, not a reporting gap,",
        "and the correct encoding is the omission of the term rather than",
        "a fixed(0) coefficient. Mean (SD) Ctrough,P1 in this analysis set",
        "was 114.97 (40.28) ng/mL. Chen 2021's Discussion attributes the",
        "null to the homogeneity of the exposure data, nearly all patients",
        "having been dosed at 100 mg q.d."
      )
    ),
    TCHOL = list(
      description = "Baseline total serum cholesterol.",
      units       = "mg/dL",
      type        = "continuous",
      notes       = paste(
        "Screened as a candidate covariate (Chen 2021 Table S4) but not",
        "retained for this endpoint. Efficacy-IC-ORR analysis set median",
        "201.00 mg/dL, range 88.00-321.00, no missing values",
        "(Chen 2021 Table 2)."
      )
    ),
    WT = list(
      description = "Baseline body weight.",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened but not retained. Efficacy-IC-ORR analysis set median",
        "64.25 kg, range 31.80-124.70 (Chen 2021 Table 1)."
      )
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator; 1 = Asian, 0 = other.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened but not retained. Efficacy-IC-ORR analysis set Asian",
        "40/132 (30%), White 72 (55%), Other 4 (3%), Black 0 (0%),",
        "Missing 16 (12%) (Chen 2021 Table 1)."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 132L,
    n_studies      = 1L,
    n_observations = "132 binary intracranial-response records (one per patient; landmark analysis, no repeated measures)",
    age_range      = "median 51.00 years, range 29.00-77.00 (Chen 2021 Table 1, Efficacy-IC-ORR population)",
    weight_range   = "median 64.25 kg, range 31.80-124.70 (Chen 2021 Table 1, Efficacy-IC-ORR population)",
    sex_female_pct = 60.0,
    race_ethnicity = c(White = 55.0, Asian = 30.0, Other = 3.0, Black = 0.0, Missing = 12.0),
    disease_state  = "advanced anaplastic lymphoma kinase (ALK)-positive or c-ROS oncogene 1 (ROS1)-positive non-small cell lung cancer with BASELINE CNS METASTASIS per independent central review; 97% had received at least one prior ALK inhibitor, 89% prior crizotinib, 64% prior chemotherapy and 58% prior CNS radiotherapy; ECOG performance status 0 (46%), 1 (50%), 2 (4%)",
    dose_range     = "lorlatinib 100 mg orally once daily (phase II expansion cohorts 2-5, the labelled dose)",
    regions        = "Multicentre international phase I/II study B7461001 (NCT01970865); regional breakdown not reported in Chen 2021",
    baseline_alkaline_phosphatase = "median 100.50 U/L, range 13.00-1,552.00, mean 142.66 (SD 153.63) (Chen 2021 Table 2, Efficacy-IC-ORR population)",
    baseline_amylase              = "median 62.00 U/L, range 13.00-218.00, mean 70.32 (SD 35.86), 4 missing (Chen 2021 Table 2, Efficacy-IC-ORR population)",
    baseline_intracranial_tumor_size = "median 31.95 mm, range 5.50-129.00, mean 40.52 (SD 27.06), 52 missing (Chen 2021 Table 2)",
    notes          = paste0(
      "This is the ONLY efficacy endpoint from Chen 2021 that can be ",
      "packaged. The paper also analysed systemic objective response ",
      "rate (ORR) in the wider 197-patient ALK-inhibitor-pretreated ",
      "set, but reports that 'none of the tested parameters, including ",
      "the lorlatinib exposure metric Cmax,P1, were significant ",
      "predictors of achieving ORR' and prints no ORR coefficient ",
      "table anywhere, so there is nothing to encode -- see the ",
      "vignette Errata. Model fit statistics for the IC-ORR model ",
      "(Chen 2021 Table 4): deviance difference vs null 15.15341 on 2 ",
      "degrees of freedom, AIC 167.510, log-likelihood -80.75511, ",
      "1 - P = 0.0005122."
    )
  )

  ini({
    # ==================================================================
    # Chen 2021 Table 4 ("Final model for exposure-response analysis:
    # IC-ORR >= 1 prior ALK inhibitors") and the printed regression
    # Eq. 5:
    #
    #   logit(p) = 3.929 - 1.015*Log(BAP) + 0.015*BAMY
    #
    # Table 4 prints BOTH an Estimate block and an OR block for the
    # same three rows, so each value is reported twice by the paper
    # itself and a third time by Eq. 5. BOTH blocks are rounded to 3
    # decimals, so the correct consistency test is an interval one --
    # does a real b exist with round(b, 3) == Estimate AND
    # round(exp(b), 3) == OR -- rather than round(log(OR), 3) ==
    # Estimate, which wrongly treats the OR block as exact. All three
    # rows pass:
    #   Estimate  3.929 admits exp(b) in [50.8291, 50.8800); OR 50.863
    #     is in that band.
    #   Estimate  0.015 admits exp(b) in [1.01506, 1.01606); OR 1.015
    #     admits [1.0145, 1.0155) -- they overlap.
    #   Estimate -1.015 admits exp(b) in [0.362223, 0.362585); OR 0.363
    #     admits [0.3625, 0.3635) -- they overlap on [0.3625, 0.362585).
    # The Log(BAP) row is the tightest of the three and is the one a
    # naive round(log(0.363), 3) = -1.013 check falsely flags; see
    # vignette Errata item 3.
    #
    # Covariates are NOT centred and NOT scaled; the intercept is an
    # extrapolated anchor at ALP = 1 U/L and AMYL = 0 U/L, not a
    # reference-patient probability.
    # ==================================================================

    # ----- Logit intercept -----
    logit_ref <- 3.929      ; label("Logit of the intracranial objective response probability at ALP = 1 U/L and AMYL = 0 U/L (unitless logit)")  # Chen 2021 Table 4, Estimate 3.929 (95% CI 1.0210, 7.0845), z = 2.553, P = 0.0107; OR block 50.863 (2.7759, 1193.3282); Eq. 5

    # ----- Covariate effects on the logit -----
    e_alp_logit  <- -1.015  ; label("Log-odds of intracranial objective response per e-fold increase in baseline alkaline phosphatase (unitless logit)")  # Chen 2021 Table 4, Estimate -1.015, z = -3.015, P = 0.0026; Eq. 5; OR block 0.363 (0.1801, 0.6778), consistent with the Estimate under the interval test above. The printed Estimate CI "(-1.7145 to 0.3889)" drops a minus sign on its upper limit -- it must be (-1.7145, -0.3889), which the OR interval confirms exactly; see vignette Errata item 2.
    e_amyl_logit <-  0.015  ; label("Log-odds of intracranial objective response per 1 U/L increase in baseline serum amylase (unitless logit)")  # Chen 2021 Table 4, Estimate 0.015 (0.0037, 0.0268), z = 2.506, P = 0.0122; OR block 1.015 (1.0037, 1.0272)

    # ----- No between-subject variability, no residual error -----
    # The source is a generalized binomial logistic regression fitted
    # with R's glm(family = "binomial"); a Bernoulli likelihood has no
    # sigma, and the analysis estimates no random effects. The tiny
    # fixed additive residual below exists only so rxode2 has an error
    # model to attach to the typical-value probability; it is NOT a
    # published quantity. See the vignette's Assumptions and deviations.
    addSd_prob_icorr <- fixed(0.001) ; label("Placeholder additive residual SD on the typical-value response probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # ----- Linear predictor (Chen 2021 Eq. 5) -----
    # log() is rxode2's natural log. ALP is log-transformed; AMYL is not.
    # No exposure term: log(Ctrough,P1) was screened and found not
    # significant for this endpoint.
    logit_icorr <- logit_ref +
      e_alp_logit  * log(ALP) +
      e_amyl_logit * AMYL

    prob_icorr <- expit(logit_icorr)

    # ----- Observation -----
    prob_icorr ~ add(addSd_prob_icorr)
  })
}
