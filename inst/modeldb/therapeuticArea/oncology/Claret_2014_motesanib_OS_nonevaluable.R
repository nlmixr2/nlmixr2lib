Claret_2014_motesanib_OS_nonevaluable <- function() {
  description <- paste(
    "Parametric overall-survival (log-normal accelerated failure time)",
    "model for the TUMOR-SIZE-NONEVALUABLE patients of the phase III",
    "MONET1 study of carboplatin/paclitaxel (CP) plus motesanib 125 mg",
    "q.d. vs. CP plus placebo in first-line advanced nonsquamous",
    "non-small cell lung cancer (NSCLC). These are the 156 of 1,090",
    "enrolled patients (14 percent) for whom no tumor-size response",
    "metric could be predicted -- 115 had no postbaseline tumor",
    "assessment and 41 had only nontarget lesions -- so the companion",
    "TTG-driven survival model cannot be applied to them. Claret et al.",
    "fitted them a separate model in which baseline serum albumin is the",
    "ONLY independent prognostic factor; there is no drug-effect term, no",
    "tumor-size term and no ODE. The model exists so that trial",
    "simulations built on modellib('Claret_2014_motesanib_tumorsize_OS')",
    "can be corrected for the nonevaluable fraction, exactly as the",
    "source paper did. Outputs are the survivor function sur(t), the",
    "hazard, the cumulative hazard and the individual median survival",
    "time. No IIV and no residual error are reported: the AFT scale",
    "parameter is the model's only source of between-subject variability.",
    sep = " "
  )
  reference <- paste(
    "Claret L, Bruno R, Lu J-F, Sun Y-N, Hsu C-P.",
    "Exploratory modeling and simulation to support development of",
    "motesanib in Asian patients with non-small cell lung cancer based on",
    "MONET1 study results.",
    "Clin Pharmacol Ther. 2014;95(4):446-451.",
    "doi:10.1038/clpt.2014.11.",
    "Parameter estimates: Supplementary Table S2",
    "('Parameter Estimates from the Final OS Model for Non-evaluable",
    "Patients'); parameter covariance matrix: Supplementary Table S4.",
    "The companion model for TS-evaluable patients is",
    "modellib('Claret_2014_motesanib_tumorsize_OS').",
    sep = " "
  )
  vignette <- "Claret_2014_motesanib_nsclc"
  units <- list(
    time          = "day",
    dosing        = "n/a (no drug-dosing events; this model has no drug-effect term -- survival depends only on baseline serum albumin)",
    concentration = "probability (the model output sur is a survival probability, not a drug concentration)"
  )

  covariateData <- list(
    ALB = list(
      description        = "Baseline serum albumin concentration. The only independent prognostic factor retained in the final OS model for TS-nonevaluable patients, and the strongest factor in the corresponding univariate Cox screen (P < 0.0001).",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Required input; NOT centred -- the coefficient applies to the raw g/L value, so mu_OS = 0.681 + 0.108 * ALB.",
        "Supplementary Table S2 labels the row 'Baseline albumin (g/L)', i.e. SI units, which is also the register's canonical unit for ALB, so no conversion is applied.",
        "The units are self-checking: at a typical 38 g/L the model gives a median OS of exp(0.681 + 0.108 * 38) = 120 days, consistent with the poor outcome expected of patients with no evaluable postbaseline tumor assessment;",
        "reading the same coefficient as per g/dL would give exp(0.681 + 0.108 * 3.8) = 3 days, which is impossible.",
        "Claret 2014 Results ('Nonevaluable patients') reports that an imbalance in baseline albumin could account for the better outcome in the CP-plus-placebo arm than in the CP-plus-motesanib arm among these patients.",
        sep = " "
      ),
      source_name        = "Baseline albumin"
    )
  )

  # Covariates screened in the Claret 2014 univariate Cox analysis of the
  # nonevaluable patients but NOT retained in the final model. Documentation
  # only; checkModelConventions() does not require these to appear in model().
  covariatesDataExcluded <- list(
    ALT = list(
      description = "Baseline alanine aminotransferase.",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Claret 2014 Results ('Nonevaluable patients'): univariate Cox P = 0.0041, but baseline albumin was the only factor with independent predictive value in the multivariate model."
    ),
    AST = list(
      description = "Baseline aspartate aminotransferase.",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Claret 2014 Results ('Nonevaluable patients'): univariate Cox P < 0.02; not retained."
    ),
    LDH = list(
      description = "Baseline lactate dehydrogenase.",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Claret 2014 Results ('Nonevaluable patients'): univariate Cox P < 0.02; not retained."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 146L,
    n_studies      = 1L,
    age_range      = "adults with advanced nonsquamous NSCLC; the age distribution is reported in the MONET1 primary publication (Scagliotti 2012, J Clin Oncol 30:2829-2836) rather than in Claret 2014",
    weight_range   = "not reported in Claret 2014",
    sex_female_pct = NA_real_,
    race_ethnicity = "21 of the 219 Asian MONET1 patients resampled for the virtual phase III simulations were TS-nonevaluable and were simulated with this model.",
    disease_state  = "Advanced (stage IIIB with pleural effusion or stage IV) or recurrent nonsquamous NSCLC, previously untreated, whose tumor-size response could not be modelled: 115 patients had no postbaseline tumor assessment and 41 had only nontarget lesions.",
    dose_range     = "Carboplatin/paclitaxel (CP) plus motesanib 125 mg orally once daily, vs. CP plus placebo. Treatment arm is NOT a term in this model.",
    regions        = "MONET1, an international, randomized, placebo-controlled, double-blind phase III study (so described in the title of Claret 2014's reference 13, Scagliotti 2012, J Clin Oncol 30:2829-2836). Claret 2014 does not print a trial-registry identifier, so none is recorded here.",
    biomarkers     = "Overall survival (OS, days); baseline serum albumin (g/L) as the sole prognostic covariate.",
    notes          = paste(
      "n_subjects is taken as 146 from the Supplementary Table S2 title, 'Parameter Estimates from the Final OS Model for Non-evaluable Patients (n = 146)'.",
      "The main text instead reports 156 nonevaluable patients (88 in the motesanib arm and 68 in the placebo arm; 115 with no postbaseline tumor assessment plus 41 with only nontarget lesions, which also sums to 156).",
      "The paper does not reconcile the two figures; the most likely explanation is that 10 of the 156 lacked a baseline albumin measurement and so could not enter the albumin regression.",
      "The 156 figure is the one used for the 14 percent nonevaluable correction applied in the simulations. See the vignette Errata.",
      sep = "\n"
    )
  )

  ini({
    # ==================================================================
    # Parametric accelerated failure time (AFT) regression with a
    # log-normal survival-time distribution, fitted with survreg() in R.
    # Claret 2014 Supplementary Table S2. Survival time is analysed in
    # DAYS (Table S2 footnote a). The parameter covariance matrix is
    # Supplementary Table S4.
    #
    # No IIV and no residual error: the AFT scale parameter is the only
    # source of between-subject survival variability, and Table S2
    # reports point estimates with standard errors only.
    # ==================================================================
    mu_os_int_logdays <- 0.681 ; label("Intercept of mu_OS on the log(days) AFT scale")                        # Table S2: Intercept = 0.681 (SE 0.881, Z 0.773, P 0.44)
    e_alb_mu_os       <- 0.108 ; label("Slope of baseline serum albumin on mu_OS (per g/L, uncentred)")        # Table S2: Baseline albumin (g/L) = 0.108 (SE 0.0237, Z 4.581, P < 0.00001)
    sigma_os_log      <- 0.502 ; label("log(scale) of the log-normal AFT survival distribution (unitless)")    # Table S2: Log(scale) = 0.502 (SE 0.0702, Z 7.153, P < 0.00001); "Scale, standard deviation of log(OS)"
  })

  model({
    # ---- Linear predictor on the log-survival-time scale ----
    # Albumin enters UNCENTRED, in g/L, exactly as reported.
    mu_os_logdays <- mu_os_int_logdays + e_alb_mu_os * ALB
    sigma_os      <- exp(sigma_os_log)

    # ---- Log-normal AFT survivor function ----
    # S(t) = 1 - Phi((log t - mu) / sigma); del_t keeps log(t) finite at
    # t = 0, and sur is floored inside the quotients so the hazard and
    # cumulative hazard stay finite in the far tail.
    del_t     <- 1e-6
    t_pos     <- t + del_t
    z_os      <- (log(t_pos) - mu_os_logdays) / sigma_os
    sur       <- 1 - pnorm(z_os)
    pdf_os    <- exp(-0.5 * z_os * z_os) / (t_pos * sigma_os * sqrt(2 * pi))
    hazard    <- pdf_os / (sur + 1e-30)
    cumhazard <- -log(sur + 1e-30)

    # Median survival time in days: exp(mu_OS), since the median of a
    # log-normal is the exponential of the mean of its logarithm.
    median_os <- exp(mu_os_logdays)
  })
}
attr(Claret_2014_motesanib_OS_nonevaluable, "message") <-
  "Log-normal AFT overall-survival model for the tumor-size-NONEVALUABLE MONET1 patients (Claret 2014, n=146). Baseline serum albumin (ALB, g/L, uncentred) is the sole prognostic factor; there is no drug-effect term, no tumor-size term, no ODE, no IIV and no residual error. Outputs: mu_os_logdays, sigma_os, median_os, sur, hazard, cumhazard. Use it to correct trial simulations built on modellib('Claret_2014_motesanib_tumorsize_OS') for the 14 percent nonevaluable fraction, as the source paper did."
Claret_2014_motesanib_OS_nonevaluable
