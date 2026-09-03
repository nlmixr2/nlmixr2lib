Bruno_2012_capecitabine_docetaxel_pfs <- function() {
  description <- paste(
    "Drug-independent parametric progression-free-survival (PFS) model in",
    "second-line locally advanced / metastatic breast cancer (Bruno 2012,",
    "pooled studies SO14999 and NO16853). Time to progression or death follows",
    "a log-normal accelerated failure-time distribution whose log-scale location",
    "is a linear function of baseline tumor size, the model-predicted relative",
    "change from baseline in tumor size at week 6 (end of cycle 2), ECOG",
    "performance status and a study effect. Unlike the companion overall-survival",
    "model this one does not retain a number-of-metastases term. The week-6",
    "tumor-size change is the only treatment-dependent input, which is what makes",
    "the model drug-independent. The model has no ODE state, no",
    "inter-individual random effects and no residual error; it exposes the",
    "progression-free probability `surv` and the median PFS time `tmed`. Time is",
    "in months. Companion models:",
    "modellib('Bruno_2012_capecitabine_docetaxel_tgi') supplies RCFB6_SLD and",
    "modellib('Bruno_2012_capecitabine_docetaxel_os') is the matching",
    "overall-survival model.",
    sep = " "
  )
  reference <- paste(
    "Bruno R, Lindbom L, Schaedeli Stark F, Chanu P, Gilberg F, Frey N, Claret L.",
    "Simulations to assess phase II noninferiority trials of different doses of",
    "capecitabine in combination with docetaxel for metastatic breast cancer.",
    "CPT Pharmacometrics Syst Pharmacol. 2012;1(12):e19.",
    "doi:10.1038/psp.2012.20.",
    "Parameter estimates are in Table 2; the accelerated-failure-time form is",
    "given in the Supplementary Information.",
    sep = " "
  )
  vignette <- "Bruno_2012_capecitabine_docetaxel_mbc"

  units <- list(
    time          = "month",
    dosing        = "n/a (no dosing events; treatment enters only through the week-6 tumor-size change covariate RCFB6_SLD)",
    concentration = "n/a (the model outputs are a progression-free probability and a median PFS time in months, not a drug concentration)"
  )

  covariateData <- list(
    TUM_SLD = list(
      description        = "Observed baseline tumor size: the sum of the longest diameters of all measurable lesions (WHO criteria, study SO14999) or of the target lesions (RECIST, study NO16853).",
      units              = "mm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Enters the log-scale location linearly and uncentred, so the intercept corresponds to TUM_SLD = 0 mm. Observed range in the pooled 888-patient dataset: 10 to 520 mm. Bruno 2012 Results: baseline tumor size influenced PFS the same way it influenced survival.",
      source_name        = "Baseline tumor size"
    ),
    RCFB6_SLD = list(
      description        = "Relative (fractional) change from baseline in tumor size at week 6, i.e. at the end of treatment cycle 2, computed from the individual predictions of the companion tumor-growth-inhibition model as (SLD(week 6) - SLD(0)) / SLD(0). Negative values are tumor shrinkage.",
      units              = "unitless (fraction of baseline)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject once week 6 has passed. This is the single treatment-dependent predictor in the model and the mechanism by which a tumor-growth-inhibition model is linked to PFS. Bruno 2012 derives it from individual predictions of the companion model modellib('Bruno_2012_capecitabine_docetaxel_tgi'); reproducing the published framework therefore requires a two-stage simulation (solve the TGI model, take the week-6 fractional change per subject, then bind it here). Bruno 2012 Results illustrate the effect with expected PFS varying from 3.8 months (95% CI 3.2-4.5) at RCFB6_SLD = +0.30 to 7.7 months (7.1-8.5) at RCFB6_SLD = -0.30.",
      source_name        = "Fractional change in tumor size at week 6"
    ),
    WHO_PS = list(
      description        = "Eastern Cooperative Oncology Group (ECOG) performance-status integer at baseline, identical by construction with the WHO performance status.",
      units              = "(integer score)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Enters linearly with the raw score. Only scores 0 and 1 occur in the Bruno 2012 cohort: 319 of 888 patients (36%) had ECOG 1 rather than 0. The reference patient therefore has WHO_PS = 0. Unlike the study term, this covariate is NOT on the 1 / 2 coding: reproducing the published expected-PFS values in Bruno 2012 Results requires WHO_PS = 0 for the reference patient.",
      source_name        = "ECOG performance status"
    ),
    STUDY_NO16853 = list(
      description        = "Study indicator for the randomized phase II noninferiority study NO16853: 1 = NO16853, 0 = the pivotal phase III study SO14999.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (study SO14999)",
      notes              = "Time-fixed per subject. PFS was longer in NO16853 than in SO14999; Bruno 2012 Results note that some of this may reflect the different progression-assessment criteria in the two studies (WHO in SO14999, RECIST in NO16853) and state that the clinical-trial simulations were conditioned on NO16853. Entered on a 1 / 2 code: the Table 2 footnote states the tabulated 0.216 corresponds to SO14999 and that for NO16853 the estimate is 0.216 x 2 = 0.432, so the model body multiplies the coefficient by (1 + STUDY_NO16853).",
      source_name        = "Study effect"
    )
  )

  covariatesDataExcluded <- list(
    MET_GE2 = list(
      description = "Binary indicator dichotomising the count of baseline metastatic sites at 2: 1 = more than one metastatic site, 0 = one or fewer.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but NOT retained in the final PFS model. Bruno 2012 tested the same covariate set on OS and PFS; the number-of-metastases term survived backward elimination for OS (Table 1, beta_4 = -0.200) but does not appear in the final PFS model (Table 2), so no coefficient can be transcribed here."
    ),
    ER_POS = list(
      description = "Estrogen-receptor status of the tumor (positive versus negative).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but NOT retained. Estrogen- and progesterone-receptor status were missing in 220 of 888 patients (25%) and were tested only on the subset with the information available; no point estimate is published."
    ),
    PR_POS = list(
      description = "Progesterone-receptor status of the tumor (positive versus negative).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but NOT retained; tested jointly with ER_POS as the composite 'both negative versus at least one positive' covariate. No point estimate is published."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 888L,
    n_studies      = 2L,
    n_events       = 790L,
    disease_state  = "locally advanced / metastatic breast cancer, second line after anthracycline pretreatment",
    dose_range     = "capecitabine 825 or 1,250 mg/m^2 twice daily on days 1-14 of each 3-week cycle plus docetaxel 75 mg/m^2 on day 1 of each cycle; study SO14999 additionally contributed a single-agent docetaxel every-3-weeks arm (Bruno 2012 does not state that arm's docetaxel dose)",
    notes          = "Pooled database of 888 patients (463 from SO14999, 425 from NO16853), of whom 790 (89%) experienced progression or died during the observation time. Median (95% CI) PFS was 5.8 months (5.5-6.3). ECOG performance status 1 rather than 0 in 319 patients (36%). Tumor response was assessed with WHO criteria in SO14999 and RECIST in NO16853, which Bruno 2012 notes may contribute to the study effect on PFS. The log-normal distribution best described the PFS time. Parameters were estimated with the censorReg function in S-PLUS 8.0, not NONMEM. Bruno 2012 does not tabulate age, weight, sex or race for this cohort."
  )

  ini({
    # Table 2 "Progression-free survival model parameter estimates (RSE, %)".
    # The Table 2 footnote states "Estimates correspond to values in month", so
    # the linear predictor below is log(median PFS time in months). Bruno 2012
    # Results: "The log-normal distribution best described the PFS time."
    ltmed0 <- 1.38
    label("Intercept beta_0: log median PFS (log-months) at the reference covariate values")  # Bruno 2012 Table 2: beta_0 = 1.38 (RSE 7%, 95% CI 1.20 to 1.56, P < 0.001)

    e_tumsld_tmed <- -0.00165
    label("Effect of baseline tumor size on log median PFS (per mm)")                         # Bruno 2012 Table 2: beta_1 = -0.00165 (RSE 26%, 95% CI -0.00248 to -0.00082, P < 0.001, delta-2LL +14.89)

    e_rcfb6sld_tmed <- -1.18
    label("Effect of the week-6 relative change from baseline in tumor size on log median PFS (per unit fraction)") # Bruno 2012 Table 2: beta_2 = -1.18 (RSE 11%, 95% CI -1.43 to -0.92, P < 0.001, delta-2LL +77.66)

    e_whops_tmed <- -0.195
    label("Effect of ECOG performance status on log median PFS (per unit score)")             # Bruno 2012 Table 2: beta_3 = -0.195 (RSE 29%, P < 0.001, delta-2LL +11.61). The printed 95% CI "(-0.306; 0.083)" does not bracket the estimate; the RSE implies (-0.306; -0.084), i.e. the upper limit lost its minus sign in typesetting. See the vignette Errata.

    e_studyno16853_tmed <- 0.216
    label("Effect of the 1/2-coded study indicator on log median PFS")                        # Bruno 2012 Table 2: beta_4 = 0.216 (RSE 26%, 95% CI 0.107 to 0.325, P < 0.001, delta-2LL +15.09); footnote c: the estimate corresponds to study SO14999 and for study NO16853 is 0.216 x 2 = 0.432

    sdlogt <- 0.799
    label("Scale sigma of the log-normal PFS distribution (SD of log PFS time)")              # Bruno 2012 Table 2: sigma = 0.799 (random variability; no RSE or CI reported)
  })
  model({
    # Log-normal accelerated-failure-time model. Supplementary Information:
    # log(T_i) = beta_0 + beta_1*x_1i + ... + beta_p*x_pi + sigma*epsilon_i with
    # epsilon_i standard normal, so exp(mu_logt) is the MEDIAN PFS time.
    #
    # The study term uses Bruno 2012's 1 / 2 coding (Table 2 footnote c),
    # reproduced here as (1 + indicator). ECOG performance status enters on its
    # raw 0 / 1 score: only that combination reproduces the expected PFS of
    # 3.8 months at RCFB6_SLD = +0.30 and 7.7 months at -0.30 quoted in Bruno
    # 2012 Results (reference patient: baseline tumor size 75 mm, ECOG 0, study
    # NO16853).
    mu_logt <- ltmed0 +
      e_tumsld_tmed * TUM_SLD +
      e_rcfb6sld_tmed * RCFB6_SLD +
      e_whops_tmed * WHO_PS +
      e_studyno16853_tmed * (1 + STUDY_NO16853)

    # Median PFS time, months.
    tmed <- exp(mu_logt)

    # Progression-free (survival) function of the log-normal distribution. The
    # 1e-6 offset keeps log(t) finite at t = 0, where the expression already
    # evaluates to 1.
    surv <- 1 - pnorm((log(t + 1e-6) - mu_logt) / sdlogt)

    # Cumulative hazard, exposed for convenience when sampling event times.
    cumhaz <- -log(surv)
  })
}
