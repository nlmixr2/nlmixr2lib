Claret_2014_DCR_OS_nsclc <- function() {
  description <- paste(
    "Parametric overall-survival (OS) model for first-line advanced non-small cell",
    "lung cancer (NSCLC), linking the week-8 RECIST disease control rate (DCR) to",
    "survival. Accelerated failure time (AFT) log-normal regression fit with the R",
    "survreg function to pooled individual-level data from 774 evaluable patients in",
    "two bevacizumab-plus-chemotherapy studies: E4599 (Phase III, 878 Western",
    "patients, bevacizumab + carboplatin/paclitaxel versus carboplatin/paclitaxel)",
    "and SAiL (Phase IV single arm, 198 Chinese patients). The location parameter of",
    "log(OS) is a linear function of two retained binary covariates: week-8 disease",
    "control (CR + PR + SD versus PD) and baseline hypoalbuminaemia (< 3.5 g/dL).",
    "Bevacizumab treatment, Chinese ethnicity and ECOG performance status were all",
    "significant in univariate Cox screening but were eliminated in backward",
    "selection -- disease control fully captured the bevacizumab effect on OS, so the",
    "authors describe the final OS model as drug independent. The model is algebraic",
    "and deterministic (no ODE state, no PK input, no IIV, no residual error): it",
    "returns the survivor function directly. The companion progression-free-survival",
    "model from the same paper is modellib('Claret_2014_DCR_PFS_nsclc')."
  )
  reference <- paste(
    "Claret L, Gupta M, Han K, Joshi A, Sarapa N, He J, Powell B, Bruno R.",
    "Prediction of overall survival or progression free survival by disease control",
    "rate at week 8 is independent of ethnicity: Western versus Chinese patients with",
    "first-line non-small cell lung cancer treated with chemotherapy with or without",
    "bevacizumab.",
    "J Clin Pharmacol. 2014;54(3):253-257. doi:10.1002/jcph.191.",
    "Parameter estimates are from Table 2 (Parameter Estimates of the OS Model).",
    "The companion PFS model from the same paper (Table 3) is available as",
    "modellib('Claret_2014_DCR_PFS_nsclc').",
    sep = " "
  )
  vignette <- "Claret_2014_DCR_OS_PFS_nsclc"

  units <- list(
    time = "month",
    dosing = "n/a (no PK input; the treatment effect enters only through the week-8 disease-control covariate RESP_DCR)",
    concentration = "probability (the model output sur_os is the probability of surviving beyond time t, not a drug concentration)"
  )

  covariateData <- list(
    RESP_DCR = list(
      description        = "Binary week-8 disease-control indicator: 1 = best RECIST response between weeks 3 and 14 (the assessment closest to week 8) was complete response, partial response or stable disease; 0 = progressive disease.",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "The paper's primary drug-effect metric. Claret 2014 Methods, 'Response Metrics':",
        "patients needed at least one response assessment between 3 and 14 weeks of treatment",
        "and the response closest to week 8 was used; 774 of 1,076 patients (629 of 878 Western,",
        "145 of 198 Chinese) were evaluable. Observed week-8 DCR by arm, computed from the",
        "Claret 2014 Table 1 response-category counts: Western chemotherapy alone",
        "(2 CR + 67 PR + 152 SD) / 320 = 0.691; Western bevacizumab + chemotherapy",
        "(4 + 133 + 128) / 309 = 0.858; Chinese bevacizumab + chemotherapy",
        "(2 + 67 + 69) / 145 = 0.952. The paper notes the Chinese DCR (95%) exceeded the",
        "Western bevacizumab DCR (85%). Additive on the log-time scale:",
        "mu_os = mu_os_int + e_dcr_mu_os * RESP_DCR + e_alb_lt35_mu_os * ALB_LT35.",
        "The landmark week is model-specific; record it here rather than assuming it,",
        "because other papers in the RESP_ family use an end-of-induction re-baseline instead."
      ),
      source_name        = "DCR"
    ),
    ALB_LT35 = list(
      description        = "Binary baseline hypoalbuminaemia indicator: 1 = baseline serum albumin below 35 g/L (3.5 g/dL); 0 = 35 g/L or above.",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "Claret 2014 Methods, 'Patient Characteristics', reports the threshold in",
        "US-convention units as '< 3.5 g/dL vs >= 3.5 g/dL, yes/no'; 3.5 g/dL = 35 g/L in the",
        "register's canonical SI unit, hence the canonical name ALB_LT35. The paper describes",
        "low albumin as the poor-prognostic category and as 'a marker of inflammation, a",
        "well-established prognostic factor of OS in cancer' (Discussion).",
        "Observed prevalence of low albumin: 29% in the Western study, 6% in the Chinese study",
        "(Claret 2014 Methods). The much lower Chinese prevalence is one of the two model",
        "covariates the authors credit with explaining the longer survival observed in the",
        "Chinese study, the other being the higher Chinese DCR.",
        "Additive on the log-time scale with a negative coefficient (shorter survival).",
        "Derive from a continuous albumin column with ALB_LT35 = as.integer(ALB < 35) when ALB",
        "is in g/L, or as.integer(ALB_gdL < 3.5) when the source column is in g/dL.",
        "The paper reports only the dichotomised covariate, so the continuous canonical ALB",
        "cannot reconstruct this effect."
      ),
      source_name        = "Albumin < 3.5 (g/dL)"
    )
  )

  # Covariates the source paper screened but did NOT retain in the final OS
  # model. Documented here rather than in covariateData so they carry their
  # provenance without triggering a "declared but not referenced" convention
  # warning: none of them appears in model().
  covariatesDataExcluded <- list(
    DRUG_BEV = list(
      description = "Bevacizumab add-on treatment-arm indicator (1 = bevacizumab + chemotherapy, 0 = chemotherapy alone).",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Significant in univariate Cox screening (deviance 13, Claret 2014 Results,",
        "'Overall Survival Model') but eliminated in backward selection: 'DCR captured the",
        "effect of bevacizumab as bevacizumab treatment was no longer significant in the",
        "multivariate model.' The Discussion draws the consequence explicitly: 'This model can",
        "be assumed to be drug independent (as there is no effect of bevacizumab treatment in",
        "the final model).' Retained in the companion PFS model, where DCR did NOT fully",
        "capture the bevacizumab effect -- see modellib('Claret_2014_DCR_PFS_nsclc')."
      )
    ),
    ECOG_GE1 = list(
      description = "Baseline ECOG performance status dichotomised at 1 (1 = ECOG >= 1, 0 = ECOG 0).",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Claret 2014 Methods dichotomises ECOG performance status as '0 vs >= 1'.",
        "Significant in univariate Cox screening (deviance 8, Claret 2014 Results) but",
        "eliminated in backward selection; only DCR and albumin survived. Observed prevalence",
        "of ECOG >= 1: 60% in the Western study, 78% in the Chinese study (Claret 2014",
        "Methods). No point estimate is published for its effect, so it cannot be reinstated",
        "from this paper."
      )
    ),
    RACE_CHINESE = list(
      description = "Chinese-ethnicity indicator, tested as the paper's ethnicity covariate (1 = Chinese SAiL cohort, 0 = Western E4599 cohort).",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Significant in univariate Cox screening (deviance 19, Claret 2014 Results) but",
        "eliminated in backward selection, and the paper reports no DCR-by-ethnicity",
        "interaction. This is the paper's central negative finding and the reason for its",
        "title: 'The covariates in the model explained the ethnicity differences as Chinese",
        "ethnicity was no longer significant in the final model and there was no interaction",
        "between DCR and Chinese ethnicity, that is, ethnicity did not impact the link between",
        "DCR and survival.' Absence of this covariate from ini() is therefore a modelling",
        "result, not an omission. No point estimate is published, so it cannot be reinstated.",
        "The canonical name RACE_CHINESE is used because the source contrast breaks out the",
        "Chinese SAiL cohort against the Western E4599 cohort as its own indicator; the paper",
        "does not resolve finer ethnic strata within either study."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 774L,
    n_studies      = 2L,
    age_range      = "adults with previously untreated advanced NSCLC; per-cohort age distributions are reported in the underlying trial publications (Sandler 2006 E4599; Crino 2010 SAiL), not in Claret 2014",
    weight_range   = "not reported in Claret 2014",
    sex_female_pct = NA_real_,
    race_ethnicity = "Western (E4599, United States) and Chinese (SAiL China subset); the paper treats ethnicity as a two-level Western-versus-Chinese contrast and does not resolve finer strata",
    disease_state  = "First-line locally advanced, metastatic or recurrent non-small cell lung cancer",
    dose_range     = "n/a (no PK input). E4599: bevacizumab 15 mg/kg IV Q3W plus carboplatin/paclitaxel, versus carboplatin/paclitaxel alone. SAiL: bevacizumab 15 mg/kg IV Q3W plus investigator-choice chemotherapy, single arm.",
    regions        = "Western (E4599, Phase III) and China (SAiL, Phase IV single arm)",
    biomarkers     = "Survival endpoint: overall survival (OS), months. Time-fixed covariates: week-8 RECIST disease control (RESP_DCR) and baseline hypoalbuminaemia (ALB_LT35).",
    notes          = paste(
      "n_subjects = 774 is the EVALUABLE population actually used to fit the survival model",
      "(Claret 2014 Methods, 'Response Metrics': 629 of 878 Western [72%] and 145 of 198",
      "Chinese [73%] patients had at least one response assessment between weeks 3 and 14).",
      "The enrolled populations were 878 + 198 = 1,076.",
      "",
      "Evaluable-population arm sizes and week-8 response categories (Claret 2014 Table 1):",
      "  Western, chemotherapy alone      : CR 2 (0.6%), PR 67 (21%), SD 152 (48%), PD 99 (31%); n = 320",
      "  Western, bevacizumab + chemo     : CR 4 (1.3%), PR 133 (43%), SD 128 (41%), PD 44 (14%); n = 309",
      "  Chinese, bevacizumab + chemo     : CR 2 (1%),   PR 67 (46%),  SD 69 (48%),  PD 7 (5%);   n = 145",
      "",
      "Baseline prognostic-factor prevalences (Claret 2014 Methods, 'Patient Characteristics'):",
      "  low albumin (< 3.5 g/dL): 29% Western, 6% Chinese",
      "  ECOG performance status >= 1: 60% Western, 78% Chinese",
      "",
      "Published outcome anchors for validating simulations:",
      "  E4599 median OS 12.3 months (bevacizumab + chemo) vs 10.3 months (chemo), HR 0.79",
      "  E4599 median PFS 6.2 months vs 4.5 months, HR 0.66",
      "  SAiL median OS 18.9 months; median time to progressive disease 8.3 months",
      "  Predictive check, Western evaluable patients (n = 629): model-predicted bevacizumab",
      "  OS hazard ratio 0.84 (95% prediction interval 0.71-0.98) versus observed 0.77",
      "",
      "Censoring for the paper's posterior predictive checks was simulated by drawing each",
      "patient's study duration from a uniform distribution between 8 and 30 months",
      "(Claret 2014 Methods, 'Survival Model Development'), consistent with the minimum and",
      "maximum time each patient stayed in the study. Reproduce that censoring scheme when",
      "comparing a simulated hazard ratio against the published one.",
      sep = "\n"
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Accelerated failure time (AFT) log-normal survival regression,
    # Claret 2014 J Clin Pharmacol 54(3):253-257 Table 2.
    #
    # Claret 2014 Methods, 'Survival Model Development', names both the
    # functional form and the software: "A parametric survival regression
    # model (using the survreg function in R version 2.13.1) was developed
    # that describes the clinical endpoint (OS or PFS) distribution as a
    # function of covariates", with the distribution selected among
    # normal, lognormal, Weibull, logistic, log-logistic and exponential;
    # Results, 'Overall Survival Model' reports "A lognormal distribution
    # had the best likelihood to describe the OS distribution".
    #
    # survreg's log-normal parameterisation is the AFT form
    #   log(T) = mu + sigma * eps,  eps ~ N(0, 1)
    # so T is log-normal with median exp(mu), and Table 2's "Log (scale)"
    # row is log(sigma). Survivor function:
    #   S(t) = 1 - Phi((log(t) - mu) / sigma).
    #
    # SIGN RECOVERY. The minus signs on the Albumin and Log(scale) rows of
    # Table 2 are typeset as U+2212 and are dropped by both the trimmed
    # markdown conversion and `pdftotext -layout`. They are recovered
    # unambiguously from the table's own Wald statistics, since z equals
    # estimate divided by SE for every row:
    #     1.855 / 0.0613 = +30.3  vs published z = 30.2
    #     0.939 / 0.0646 = +14.5  vs published z = 14.5
    #    -0.270 / 0.0601 =  -4.49 vs published z = -4.5   (sign required)
    #    -0.370 / 0.0308 = -12.0  vs published z = -12.0  (sign required)
    # The negative albumin coefficient is independently corroborated by the
    # paper's prose: "the probability to survive decreases in patients with
    # low baseline albumin (< 3.5 g/dL; disease severity)".
    #
    # No IIV and no residual error: this is a population-level parametric
    # time-to-event regression, and all subject-level variability is
    # carried by sigma_os, the AFT scale. Table 2 reports point estimates
    # and standard errors only.
    # ------------------------------------------------------------------

    mu_os_int <- 1.855; label("Intercept of the log-normal AFT location parameter for OS (log-months) in a reference patient: progressive disease at week 8 and baseline albumin >= 3.5 g/dL")  # Table 2 Intercept = 1.855 (SE 0.0613, z 30.2, P < .00001)

    e_dcr_mu_os <- 0.939; label("Additive effect of week-8 disease control (RESP_DCR) on the OS log-normal AFT location parameter (log-months); median OS is exp(0.939) = 2.56-fold longer with disease control")  # Table 2 DCR = 0.939 (SE 0.0646, z 14.5, P < .00001)

    e_alb_lt35_mu_os <- -0.270; label("Additive effect of baseline albumin < 3.5 g/dL (ALB_LT35) on the OS log-normal AFT location parameter (log-months); median OS is exp(-0.270) = 0.763-fold, i.e. 23.7% shorter")  # Table 2 "Albumin < 3.5 (g/dL)" = -0.270 (SE 0.0601, z -4.5, P < .00001); minus sign recovered from z = estimate/SE and corroborated by the Results prose

    lsigma_os <- -0.370; label("log(scale) of the OS log-normal AFT distribution; sigma_os = exp(-0.370) = 0.691 is the standard deviation of log(OS) in log-months")  # Table 2 "Log (scale)" = -0.370 (SE 0.0308, z -12.0, P < .00001); minus sign recovered from z = estimate/SE
  })

  model({
    # ------------------------------------------------------------------
    # 1. Location parameter of log(OS), in log-months. Linear in the two
    #    retained binary covariates (Claret 2014 Table 2). Bevacizumab
    #    treatment, Chinese ethnicity and ECOG performance status are
    #    deliberately absent -- they were eliminated in backward selection;
    #    see covariatesDataExcluded.
    # ------------------------------------------------------------------
    mu_os <- mu_os_int + e_dcr_mu_os * RESP_DCR + e_alb_lt35_mu_os * ALB_LT35

    # Median survival time in months, exp(mu_os), because the log-normal
    # median is exp of the AFT location parameter. This is the quantity the
    # paper's covariate-effect statements are expressed against.
    median_os <- exp(mu_os)

    # ------------------------------------------------------------------
    # 2. AFT scale on the log-time axis.
    # ------------------------------------------------------------------
    sigma_os <- exp(lsigma_os)

    # ------------------------------------------------------------------
    # 3. Log-normal survivor function S(t) = 1 - Phi((log t - mu)/sigma),
    #    plus the density and hazard for users driving a time-to-event
    #    simulation from this model. del_t keeps log(t) finite at t = 0,
    #    where S is 1 by definition; the epsilon in the hazard denominator
    #    keeps h(t) finite in the far tail where S underflows to 0. Same
    #    numerical device as Struemper_2025_tumorsize_OS_nsclc.
    # ------------------------------------------------------------------
    del_t <- 1e-6
    t_pos <- t + del_t
    z_os  <- (log(t_pos) - mu_os) / sigma_os

    sur_os    <- 1 - pnorm(z_os)
    pdf_os    <- exp(-0.5 * z_os * z_os) / (t_pos * sigma_os * sqrt(2 * pi))
    hazard_os <- pdf_os / (sur_os + 1e-30)

    # Cumulative hazard, for users assembling a joint or competing-risk
    # simulation on top of this endpoint model.
    cumhazard_os <- -log(sur_os + 1e-30)
  })
}
attr(Claret_2014_DCR_OS_nsclc, "message") <- paste(
  "Log-normal AFT overall-survival model for first-line advanced NSCLC (Claret 2014,",
  "n = 774 evaluable across E4599 and SAiL). log(OS in months) = 1.855 + 0.939 * RESP_DCR",
  "- 0.270 * ALB_LT35 + sigma * N(0, 1) with sigma = exp(-0.370) = 0.691.",
  "Algebraic and deterministic: no ODE state, no PK input, no IIV, no residual error.",
  "Outputs: mu_os, median_os (months), sigma_os, sur_os (S(t)), hazard_os, cumhazard_os.",
  "Bevacizumab treatment, ECOG performance status and Chinese ethnicity were screened and",
  "eliminated (see covariatesDataExcluded); the authors describe the final OS model as drug",
  "independent because week-8 disease control fully captured the bevacizumab effect.",
  "Companion PFS model: modellib('Claret_2014_DCR_PFS_nsclc')."
)
Claret_2014_DCR_OS_nsclc
