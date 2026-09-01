Claret_2014_DCR_PFS_nsclc <- function() {
  description <- paste(
    "Parametric progression-free-survival (PFS) model for first-line advanced",
    "non-small cell lung cancer (NSCLC), linking the week-8 RECIST disease control",
    "rate (DCR) to PFS. Accelerated failure time (AFT) log-normal regression fit with",
    "the R survreg function to pooled individual-level data from 774 evaluable",
    "patients in two bevacizumab-plus-chemotherapy studies: E4599 (Phase III, 878",
    "Western patients, bevacizumab + carboplatin/paclitaxel versus",
    "carboplatin/paclitaxel) and SAiL (Phase IV single arm, 198 Chinese patients).",
    "The location parameter of log(PFS) is a linear function of two retained binary",
    "covariates: week-8 disease control (CR + PR + SD versus PD) and bevacizumab",
    "treatment. Unlike the companion overall-survival model, disease control did NOT",
    "fully capture the bevacizumab effect here, so a separate bevacizumab term is",
    "retained and the authors state the PFS model cannot be assumed drug independent.",
    "Baseline albumin, ECOG performance status and Chinese ethnicity were all",
    "significant in univariate Cox screening but were eliminated in backward",
    "selection. The model is algebraic and deterministic (no ODE state, no PK input,",
    "no IIV, no residual error): it returns the survivor function directly. The",
    "companion overall-survival model from the same paper is",
    "modellib('Claret_2014_DCR_OS_nsclc')."
  )
  reference <- paste(
    "Claret L, Gupta M, Han K, Joshi A, Sarapa N, He J, Powell B, Bruno R.",
    "Prediction of overall survival or progression free survival by disease control",
    "rate at week 8 is independent of ethnicity: Western versus Chinese patients with",
    "first-line non-small cell lung cancer treated with chemotherapy with or without",
    "bevacizumab.",
    "J Clin Pharmacol. 2014;54(3):253-257. doi:10.1002/jcph.191.",
    "Parameter estimates are from Table 3 (Parameter Estimates of the PFS Model).",
    "The companion OS model from the same paper (Table 2) is available as",
    "modellib('Claret_2014_DCR_OS_nsclc').",
    sep = " "
  )
  vignette <- "Claret_2014_DCR_OS_PFS_nsclc"

  units <- list(
    time = "month",
    dosing = "n/a (no PK input; bevacizumab enters only as the binary treatment-arm covariate DRUG_BEV)",
    concentration = "probability (the model output sur_pfs is the probability of being alive and progression-free beyond time t, not a drug concentration)"
  )

  covariateData <- list(
    RESP_DCR = list(
      description        = "Binary week-8 disease-control indicator: 1 = best RECIST response between weeks 3 and 14 (the assessment closest to week 8) was complete response, partial response or stable disease; 0 = progressive disease.",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "The paper's primary drug-effect metric, identical in definition to the one used by",
        "the companion OS model but carrying a different coefficient here. Claret 2014",
        "Methods, 'Response Metrics': patients needed at least one response assessment between",
        "3 and 14 weeks of treatment and the response closest to week 8 was used; 774 of 1,076",
        "patients were evaluable. Observed week-8 DCR by arm, computed from the Claret 2014",
        "Table 1 response-category counts: Western chemotherapy alone (2 CR + 67 PR + 152 SD)",
        "/ 320 = 0.691; Western bevacizumab + chemotherapy (4 + 133 + 128) / 309 = 0.858;",
        "Chinese bevacizumab + chemotherapy (2 + 67 + 69) / 145 = 0.952.",
        "Additive on the log-time scale:",
        "mu_pfs = mu_pfs_int + e_dcr_mu_pfs * RESP_DCR + e_bev_mu_pfs * DRUG_BEV.",
        "Note the DCR effect on PFS (1.378) is markedly larger than on OS (0.939), consistent",
        "with the week-8 response assessment being a much more direct read-out of progression",
        "than of death.",
        "There is an intentional partial confounding between this covariate and DRUG_BEV:",
        "bevacizumab raises the probability that RESP_DCR equals 1 AND additionally lengthens",
        "PFS at a given RESP_DCR. Both paths must be simulated to reproduce the paper's",
        "published bevacizumab hazard ratio."
      ),
      source_name        = "DCR"
    ),
    DRUG_BEV = list(
      description        = "Binary bevacizumab add-on treatment-arm indicator: 1 = bevacizumab plus chemotherapy, 0 = chemotherapy backbone alone.",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "Retained in the final PFS model with a positive coefficient (longer PFS on",
        "bevacizumab). Claret 2014 Results, 'PFS Model': 'DCR did not fully capture the effect",
        "of bevacizumab on PFS', and the Discussion draws the consequence: 'the PFS model",
        "cannot be assumed to be drug independent and it has less value than the OS model.'",
        "This is the structural difference between the two models in the paper -- the",
        "companion OS model has no bevacizumab term at all.",
        "Reference arm is ACTIVE chemotherapy, not placebo: E4599 randomised bevacizumab +",
        "carboplatin/paclitaxel against carboplatin/paclitaxel alone. Set DRUG_BEV = 1 for all",
        "SAiL (Chinese) subjects, which was a single-arm bevacizumab + chemotherapy study,",
        "and for the E4599 bevacizumab arm; set 0 for the E4599 chemotherapy-alone arm.",
        "Because the Chinese cohort is entirely bevacizumab-treated, this covariate is",
        "perfectly confounded with study in that cohort; the bevacizumab contrast is",
        "identified from the randomised Western arms only."
      ),
      source_name        = "Bevacizumab"
    )
  )

  # Covariates the source paper screened but did NOT retain in the final PFS
  # model. Documented here rather than in covariateData so they carry their
  # provenance without triggering a "declared but not referenced" convention
  # warning: none of them appears in model().
  covariatesDataExcluded <- list(
    ALB_LT35 = list(
      description = "Baseline hypoalbuminaemia indicator (1 = serum albumin below 35 g/L / 3.5 g/dL, 0 = at or above).",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Significant in univariate Cox screening for PFS (deviance 11, P < .0001; Claret 2014",
        "Results, 'PFS Model') but eliminated in backward selection -- only DCR and",
        "bevacizumab treatment survived. This is the mirror image of the OS model, where",
        "albumin was retained and bevacizumab was dropped; see",
        "modellib('Claret_2014_DCR_OS_nsclc'). No PFS point estimate is published for the",
        "albumin effect, so it cannot be reinstated from this paper."
      )
    ),
    ECOG_GE1 = list(
      description = "Baseline ECOG performance status dichotomised at 1 (1 = ECOG >= 1, 0 = ECOG 0).",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Claret 2014 Methods dichotomises ECOG performance status as '0 vs >= 1'. The weakest",
        "of the screened covariates for PFS (univariate Cox deviance 2.1, P = .0385; Claret",
        "2014 Results) and eliminated in backward selection, which used a P < .01 cut-off.",
        "Observed prevalence of ECOG >= 1: 60% Western, 78% Chinese (Claret 2014 Methods).",
        "No point estimate is published, so it cannot be reinstated."
      )
    ),
    RACE_CHINESE = list(
      description = "Chinese-ethnicity indicator, tested as the paper's ethnicity covariate (1 = Chinese SAiL cohort, 0 = Western E4599 cohort).",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Significant in univariate Cox screening (deviance 12; Claret 2014 Results, 'PFS",
        "Model') but eliminated in backward selection, and the paper reports no",
        "DCR-by-ethnicity interaction: 'There was no interaction between DCR and Chinese",
        "ethnicity, that is, ethnicity did not impact the link between DCR and PFS.'",
        "This is the paper's central negative finding and the reason for its title, so the",
        "absence of this covariate from ini() is a modelling result, not an omission.",
        "No point estimate is published, so it cannot be reinstated."
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
    biomarkers     = "Survival endpoint: progression-free survival (PFS), months. Time-fixed covariates: week-8 RECIST disease control (RESP_DCR) and bevacizumab treatment arm (DRUG_BEV).",
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
      "  E4599 median PFS 6.2 months (bevacizumab + chemo) vs 4.5 months (chemo), HR 0.66",
      "  E4599 median OS 12.3 months vs 10.3 months, HR 0.79",
      "  SAiL median time to progressive disease 8.3 months; median OS 18.9 months",
      "  Predictive check, Western evaluable patients (n = 629): model-predicted bevacizumab",
      "  PFS hazard ratio 0.59 (95% prediction interval 0.49-0.72) versus observed 0.58",
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
    # Claret 2014 J Clin Pharmacol 54(3):253-257 Table 3.
    #
    # Claret 2014 Methods, 'Survival Model Development', names both the
    # functional form and the software: "A parametric survival regression
    # model (using the survreg function in R version 2.13.1) was developed
    # that describes the clinical endpoint (OS or PFS) distribution as a
    # function of covariates"; Results, 'PFS Model' reports "A lognormal
    # distribution had the best likelihood to describe the PFS
    # distribution".
    #
    # survreg's log-normal parameterisation is the AFT form
    #   log(T) = mu + sigma * eps,  eps ~ N(0, 1)
    # so T is log-normal with median exp(mu), and Table 3's "Log (scale)"
    # row is log(sigma). Survivor function:
    #   S(t) = 1 - Phi((log(t) - mu) / sigma).
    #
    # SIGN RECOVERY. The minus sign on the Log(scale) row of Table 3 is
    # typeset as U+2212 and is dropped by both the trimmed markdown
    # conversion and `pdftotext -layout`. It is recovered unambiguously
    # from the table's own Wald statistics, since z equals estimate
    # divided by SE for every row:
    #     0.466 / 0.0477 =  +9.77 vs published z = 9.8
    #     1.378 / 0.0525 = +26.2  vs published z = 26.2
    #     0.232 / 0.0429 =  +5.41 vs published z = 5.4    (POSITIVE)
    #    -0.641 / 0.0282 = -22.7  vs published z = -22.7  (sign required)
    # The Bevacizumab row is positive: longer PFS on bevacizumab, matching
    # the paper's prose that "the probability of progression or death
    # decreases in patients with disease control and with bevacizumab
    # treatment" and the published protective hazard ratio of 0.59.
    #
    # No IIV and no residual error: this is a population-level parametric
    # time-to-event regression, and all subject-level variability is
    # carried by sigma_pfs, the AFT scale. Table 3 reports point estimates
    # and standard errors only.
    # ------------------------------------------------------------------

    mu_pfs_int <- 0.466; label("Intercept of the log-normal AFT location parameter for PFS (log-months) in a reference patient: progressive disease at week 8 and chemotherapy alone")  # Table 3 Intercept = 0.466 (SE 0.0477, z 9.8, P < .00001)

    e_dcr_mu_pfs <- 1.378; label("Additive effect of week-8 disease control (RESP_DCR) on the PFS log-normal AFT location parameter (log-months); median PFS is exp(1.378) = 3.97-fold longer with disease control")  # Table 3 DCR = 1.378 (SE 0.0525, z 26.2, P < .00001)

    e_bev_mu_pfs <- 0.232; label("Additive effect of bevacizumab treatment (DRUG_BEV) on the PFS log-normal AFT location parameter (log-months), over and above the effect captured by disease control; median PFS is exp(0.232) = 1.26-fold longer")  # Table 3 Bevacizumab = 0.232 (SE 0.0429, z 5.4, P < .00001)

    lsigma_pfs <- -0.641; label("log(scale) of the PFS log-normal AFT distribution; sigma_pfs = exp(-0.641) = 0.527 is the standard deviation of log(PFS) in log-months")  # Table 3 "Log (scale)" = -0.641 (SE 0.0282, z -22.7, P < .00001); minus sign recovered from z = estimate/SE
  })

  model({
    # ------------------------------------------------------------------
    # 1. Location parameter of log(PFS), in log-months. Linear in the two
    #    retained binary covariates (Claret 2014 Table 3). Baseline
    #    albumin, Chinese ethnicity and ECOG performance status are
    #    deliberately absent -- they were eliminated in backward selection;
    #    see covariatesDataExcluded.
    # ------------------------------------------------------------------
    mu_pfs <- mu_pfs_int + e_dcr_mu_pfs * RESP_DCR + e_bev_mu_pfs * DRUG_BEV

    # Median PFS in months, exp(mu_pfs), because the log-normal median is
    # exp of the AFT location parameter. This is the quantity the paper's
    # covariate-effect statements are expressed against.
    median_pfs <- exp(mu_pfs)

    # ------------------------------------------------------------------
    # 2. AFT scale on the log-time axis.
    # ------------------------------------------------------------------
    sigma_pfs <- exp(lsigma_pfs)

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
    z_pfs <- (log(t_pos) - mu_pfs) / sigma_pfs

    sur_pfs    <- 1 - pnorm(z_pfs)
    pdf_pfs    <- exp(-0.5 * z_pfs * z_pfs) / (t_pos * sigma_pfs * sqrt(2 * pi))
    hazard_pfs <- pdf_pfs / (sur_pfs + 1e-30)

    # Cumulative hazard, for users assembling a joint or competing-risk
    # simulation on top of this endpoint model.
    cumhazard_pfs <- -log(sur_pfs + 1e-30)
  })
}
attr(Claret_2014_DCR_PFS_nsclc, "message") <- paste(
  "Log-normal AFT progression-free-survival model for first-line advanced NSCLC (Claret 2014,",
  "n = 774 evaluable across E4599 and SAiL). log(PFS in months) = 0.466 + 1.378 * RESP_DCR",
  "+ 0.232 * DRUG_BEV + sigma * N(0, 1) with sigma = exp(-0.641) = 0.527.",
  "Algebraic and deterministic: no ODE state, no PK input, no IIV, no residual error.",
  "Outputs: mu_pfs, median_pfs (months), sigma_pfs, sur_pfs (S(t)), hazard_pfs, cumhazard_pfs.",
  "Unlike the companion OS model, week-8 disease control did NOT fully capture the bevacizumab",
  "effect, so DRUG_BEV is retained and the model is NOT drug independent. Baseline albumin,",
  "ECOG performance status and Chinese ethnicity were screened and eliminated (see",
  "covariatesDataExcluded). Companion OS model: modellib('Claret_2014_DCR_OS_nsclc')."
)
Claret_2014_DCR_PFS_nsclc
