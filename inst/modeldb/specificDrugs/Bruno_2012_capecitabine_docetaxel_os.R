Bruno_2012_capecitabine_docetaxel_os <- function() {
  description <- paste(
    "Drug-independent parametric overall-survival (OS) model in second-line",
    "locally advanced / metastatic breast cancer (Bruno 2012, pooled studies",
    "SO14999 and NO16853). Survival time follows a log-normal accelerated",
    "failure-time distribution whose log-scale location is a linear function of",
    "baseline tumor size, the model-predicted relative change from baseline in",
    "tumor size at week 6 (end of cycle 2), ECOG performance status, the number",
    "of metastatic sites and a study effect. The week-6 tumor-size change is the",
    "only treatment-dependent input, which is what makes the model",
    "drug-independent: any tumor-growth-inhibition model can supply it. The",
    "model has no ODE state, no inter-individual random effects and no residual",
    "error; it exposes the survival probability `surv` and the median survival",
    "time `tmed`. Time is in months. Companion models:",
    "modellib('Bruno_2012_capecitabine_docetaxel_tgi') supplies RCFB6_SLD and",
    "modellib('Bruno_2012_capecitabine_docetaxel_pfs') is the matching",
    "progression-free-survival model.",
    sep = " "
  )
  reference <- paste(
    "Bruno R, Lindbom L, Schaedeli Stark F, Chanu P, Gilberg F, Frey N, Claret L.",
    "Simulations to assess phase II noninferiority trials of different doses of",
    "capecitabine in combination with docetaxel for metastatic breast cancer.",
    "CPT Pharmacometrics Syst Pharmacol. 2012;1(12):e19.",
    "doi:10.1038/psp.2012.20.",
    "Parameter estimates are in Table 1; the accelerated-failure-time form is",
    "given in the Supplementary Information.",
    sep = " "
  )
  vignette <- "Bruno_2012_capecitabine_docetaxel_mbc"

  units <- list(
    time          = "month",
    dosing        = "n/a (no dosing events; treatment enters only through the week-6 tumor-size change covariate RCFB6_SLD)",
    concentration = "n/a (the model outputs are a survival probability and a median survival time in months, not a drug concentration)"
  )

  covariateData <- list(
    TUM_SLD = list(
      description        = "Observed baseline tumor size: the sum of the longest diameters of all measurable lesions (WHO criteria, study SO14999) or of the target lesions (RECIST, study NO16853).",
      units              = "mm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Enters the log-scale location linearly and uncentred, so the intercept corresponds to TUM_SLD = 0 mm. Observed range in the pooled 888-patient dataset: 10 to 520 mm. Bruno 2012 Results: survival decreased with increasing baseline tumor size.",
      source_name        = "Baseline tumor size"
    ),
    RCFB6_SLD = list(
      description        = "Relative (fractional) change from baseline in tumor size at week 6, i.e. at the end of treatment cycle 2, computed from the individual predictions of the companion tumor-growth-inhibition model as (SLD(week 6) - SLD(0)) / SLD(0). Negative values are tumor shrinkage.",
      units              = "unitless (fraction of baseline)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject once week 6 has passed. This is the single treatment-dependent predictor in the model and the mechanism by which a tumor-growth-inhibition model is linked to survival. Bruno 2012 derives it from individual predictions of the companion model modellib('Bruno_2012_capecitabine_docetaxel_tgi'); reproducing the published framework therefore requires a two-stage simulation (solve the TGI model, take the week-6 fractional change per subject, then bind it here). Bruno 2012 Results illustrate the effect with expected median survival varying from 11.8 months (95% CI 10.1-14.0) at RCFB6_SLD = +0.30 to 19.3 months (17.3-21.6) at RCFB6_SLD = -0.30.",
      source_name        = "Fractional change in tumor size at week 6"
    ),
    WHO_PS = list(
      description        = "Eastern Cooperative Oncology Group (ECOG) performance-status integer at baseline, identical by construction with the WHO performance status.",
      units              = "(integer score)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Enters linearly with the raw score. Only scores 0 and 1 occur in the Bruno 2012 cohort: 319 of 888 patients (36%) had ECOG 1 rather than 0 (Bruno 2012 Results, model for OS). The reference patient therefore has WHO_PS = 0. Unlike the number-of-metastases and study terms, this covariate is NOT on the 1 / 2 coding: reproducing the published expected-survival values in Bruno 2012 Results requires WHO_PS = 0 for the reference patient.",
      source_name        = "ECOG performance status"
    ),
    MET_GE2 = list(
      description        = "Binary indicator dichotomising the count of baseline metastatic sites at 2: 1 = more than one metastatic site, 0 = one or fewer metastatic sites.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (one or fewer metastatic sites)",
      notes              = "Time-fixed per subject. 687 of 888 patients (77%) had more than one metastatic site (Bruno 2012 Results, model for OS). Bruno 2012 enters the covariate on a 1 / 2 code rather than 0 / 1: the Table 1 footnote states the tabulated -0.200 corresponds to <=1 metastasis and that for >1 metastases the estimate is -0.2 x 2 = -0.4. The model body therefore multiplies the coefficient by (1 + MET_GE2).",
      source_name        = "Number of metastases"
    ),
    STUDY_NO16853 = list(
      description        = "Study indicator for the randomized phase II noninferiority study NO16853: 1 = NO16853, 0 = the pivotal phase III study SO14999.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (study SO14999)",
      notes              = "Time-fixed per subject. Survival was longer in NO16853 than in SO14999; Bruno 2012 Results state the prognostic factors in the model could not explain the difference and attribute it to a change in standard of care between the two studies, and report no interaction between the study effect and the tumor-shrinkage effect. Entered on a 1 / 2 code as for the metastases term: the Table 1 footnote states the tabulated 0.131 corresponds to SO14999 and that for NO16853 the estimate is 0.131 x 2 = 0.262, so the model body multiplies the coefficient by (1 + STUDY_NO16853).",
      source_name        = "Study effect"
    )
  )

  covariatesDataExcluded <- list(
    ER_POS = list(
      description = "Estrogen-receptor status of the tumor (positive versus negative).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but NOT retained in the final model. Bruno 2012 Results: estrogen- and progesterone-receptor status were missing in 220 of 888 patients (25%) and so could only be tested on the subset with the information available. Patients whose tumors expressed at least one receptor survived significantly longer (P < 0.001), and there was no interaction between that effect and the effect of the predicted week-6 fractional change in tumor size, but the effect was dropped from the final model because of the missingness. No point estimate is published, so no coefficient can be transcribed."
    ),
    PR_POS = list(
      description = "Progesterone-receptor status of the tumor (positive versus negative).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but NOT retained in the final model; tested jointly with ER_POS as the composite 'both negative versus at least one positive' covariate. See the ER_POS entry for the reason it was dropped. No point estimate is published."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 888L,
    n_studies      = 2L,
    n_events       = 556L,
    disease_state  = "locally advanced / metastatic breast cancer, second line after anthracycline pretreatment",
    dose_range     = "capecitabine 825 or 1,250 mg/m^2 twice daily on days 1-14 of each 3-week cycle plus docetaxel 75 mg/m^2 on day 1 of each cycle; study SO14999 additionally contributed a single-agent docetaxel every-3-weeks arm (Bruno 2012 does not state that arm's docetaxel dose)",
    notes          = "Pooled database of 888 patients (463 from SO14999, 425 from NO16853), of whom 556 (63%) died during the observation period. Median (95% CI) survival was 14.8 months (13.6-16.0). ECOG performance status 1 rather than 0 in 319 patients (36%); more than one metastatic site in 687 patients (77%). Estrogen- and progesterone-receptor status were missing in 220 patients (25%). The log-normal distribution best described survival time, as previously observed in colorectal cancer (Claret 2009). Parameters were estimated with the censorReg function in S-PLUS 8.0, not NONMEM. Bruno 2012 does not tabulate age, weight, sex or race for this cohort."
  )

  ini({
    # Table 1 "Survival model parameter estimates (RSE, %)". The Table 1
    # footnote states "Estimates correspond to values in month", i.e. the
    # linear predictor below is log(median survival time in months). Bruno 2012
    # Methods: "The probability density function that best described the
    # observed survival time was selected among normal, log normal, Weibull,
    # logistic, loglogistic, exponential, and extreme"; Results: "The log-normal
    # distribution best described the survival time".
    ltmed0 <- 3.02
    label("Intercept beta_0: log median OS (log-months) at the reference covariate values")   # Bruno 2012 Table 1: beta_0 = 3.02 (RSE 5%, 95% CI 2.72 to 3.31, P < 0.001)

    e_tumsld_tmed <- -0.00231
    label("Effect of baseline tumor size on log median OS (per mm)")                          # Bruno 2012 Table 1: beta_1 = -0.00231 (RSE 19%, 95% CI -0.00316 to -0.00146, P < 0.001, delta-2LL +27.97)

    e_rcfb6sld_tmed <- -0.801
    label("Effect of the week-6 relative change from baseline in tumor size on log median OS (per unit fraction)") # Bruno 2012 Table 1: beta_2 = -0.801 (RSE 17%, 95% CI -1.060 to -0.541, P < 0.001, delta-2LL +36.16)

    e_whops_tmed <- -0.352
    label("Effect of ECOG performance status on log median OS (per unit score)")              # Bruno 2012 Table 1: beta_3 = -0.352 (RSE 17%, 95% CI -0.468 to -0.236, P < 0.001, delta-2LL +34.92)

    e_metge2_tmed <- -0.200
    label("Effect of the 1/2-coded number-of-metastases indicator on log median OS")          # Bruno 2012 Table 1: beta_4 = -0.200 (RSE 36%, 95% CI -0.339 to -0.060, P = 0.00502, delta-2LL +7.89); footnote c: the estimate corresponds to <=1 metastasis and for >1 metastases is -0.2 x 2 = -0.4

    e_studyno16853_tmed <- 0.131
    label("Effect of the 1/2-coded study indicator on log median OS")                         # Bruno 2012 Table 1: beta_5 = 0.131 (RSE 44%, 95% CI 0.017 to 0.245, P = 0.0239, delta-2LL +5.16); footnote c: the estimate corresponds to study SO14999 and for study NO16853 is 0.131 x 2 = 0.262

    sdlogt <- 0.773
    label("Scale sigma of the log-normal survival distribution (SD of log survival time)")    # Bruno 2012 Table 1: sigma = 0.773 (random variability; no RSE or CI reported)
  })
  model({
    # Log-normal accelerated-failure-time model. Supplementary Information:
    # log(T_i) = beta_0 + beta_1*x_1i + ... + beta_p*x_pi + sigma*epsilon_i with
    # epsilon_i standard normal, so exp(mu_logt) is the MEDIAN survival time.
    #
    # The number-of-metastases and study terms use Bruno 2012's 1 / 2 coding
    # (Table 1 footnote c), reproduced here as (1 + indicator). ECOG performance
    # status enters on its raw 0 / 1 score, NOT on the 1 / 2 coding: only that
    # combination reproduces the expected median survival of 11.8 months at
    # RCFB6_SLD = +0.30 and 19.3 months at -0.30 quoted in Bruno 2012 Results
    # (reference patient: baseline tumor size 75 mm, ECOG 0, >1 metastatic site,
    # study NO16853).
    mu_logt <- ltmed0 +
      e_tumsld_tmed * TUM_SLD +
      e_rcfb6sld_tmed * RCFB6_SLD +
      e_whops_tmed * WHO_PS +
      e_metge2_tmed * (1 + MET_GE2) +
      e_studyno16853_tmed * (1 + STUDY_NO16853)

    # Median survival time, months.
    tmed <- exp(mu_logt)

    # Survival function of the log-normal distribution. The 1e-6 offset keeps
    # log(t) finite at t = 0, where the expression already evaluates to 1.
    surv <- 1 - pnorm((log(t + 1e-6) - mu_logt) / sdlogt)

    # Cumulative hazard, exposed for convenience when sampling event times.
    cumhaz <- -log(surv)
  })
}
