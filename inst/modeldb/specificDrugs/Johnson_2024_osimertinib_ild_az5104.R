Johnson_2024_osimertinib_ild_az5104 <- function() {
  description <- paste0(
    "Penalized (Firth) logistic-regression exposure-safety model for ",
    "interstitial lung disease (ILD) and ILD-like events in patients ",
    "with EGFR-mutated non-small cell lung cancer treated with ",
    "osimertinib in the first-line or second-or-later-line setting, ",
    "driven by exposure to the active metabolite AZ5104 (Johnson 2024, ",
    "final Model 3 of Table S7). The probability of an ILD-like event ",
    "is expit(a + b * log(AUC_AZ5104) + c * RACE_JAPANESE + ",
    "d * RACE_ASIAN), where AUC_AZ5104 is the individual steady-state ",
    "AZ5104 AUC in nmol*h/L derived from the companion population PK ",
    "model Johnson_2025_osimertinib, and the race indicators are ",
    "Japanese and Asian non-Japanese against a non-Asian reference. ",
    "There is no PK layer and no ODE: exposure is supplied as data. ",
    "Osimertinib and AZ5104 exposures are highly correlated and were ",
    "never entered together; the parent-driven alternative is ",
    "Johnson_2024_osimertinib_ild_osimertinib."
  )
  reference <- paste(
    "Johnson M, Lin YW, Schmidt H, Sunnaker M, Van Maanen E, Huang X,",
    "Rukazenkov Y, Tomkinson H, Vishwanathan K. Exposure-response",
    "modelling of osimertinib in patients with non-small cell lung",
    "cancer. Br J Clin Pharmacol. 2024;90(12):3263-3276.",
    "doi:10.1111/bcp.16199. PMCID PMC11602960.",
    sep = " "
  )
  vignette <- "Johnson_2024_osimertinib_exposure_response"

  units <- list(
    time = "n/a (static landmark exposure-response model; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the AUC_AZ5104 covariate column)",
    concentration = "prob_ild (probability of an ILD or ILD-like event, 0-1)"
  )

  covariateData <- list(
    AUC_AZ5104 = list(
      description = "Individual steady-state AUC of the osimertinib active metabolite AZ5104 over the once-daily dosing interval, derived from the individual apparent AZ5104 clearance of the population PK model and supplied as data (this model has no PK layer).",
      units = "nmol*h/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Johnson 2024 Methods 'Estimation of plasma exposure': AUCss of osimertinib and/or AZ5104 derived from the individual apparent clearances of the population PK analysis (reference 23 = Johnson 2025, modellib('Johnson_2025_osimertinib')) and the first dose received. At steady state the molar AZ5104 AUC is fm * dose / (CLmetabolite/F)_i with the fixed conversion fraction fm = 0.25, i.e. AUC_AZ5104 = 0.25 * dose_mg / cl_az5104 * 1e6 / 499.62 using the cl_az5104 (L/h) of Johnson_2025_osimertinib. Unit nM*h, matching the osimertinib AUCss axes of Johnson 2024 Figures 2 and 3. Enters the logit on the natural-log scale, uncentred.",
      source_name = "AUCss (AZ5104)"
    ),
    RACE_JAPANESE = list(
      description = "Japanese race indicator (1 = Japanese, 0 = otherwise).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the reference group of the race covariates is non-Asian: White, Black, Hispanic/Latino and other non-Asian patients).",
      notes = "Additive on the logit (Johnson 2024 Supplementary Methods, categorical covariate term beta_k,l * U_k,l,i). Mutually exclusive with RACE_ASIAN in this model: a Japanese patient has RACE_JAPANESE = 1 and RACE_ASIAN = 0.",
      source_name = "Japanese"
    ),
    RACE_ASIAN = list(
      description = "Asian non-Japanese race indicator (1 = Asian but not Japanese, 0 = otherwise).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the reference group of the race covariates is non-Asian).",
      notes = "In this model RACE_ASIAN EXCLUDES Japanese patients, who carry RACE_JAPANESE instead (the Wade_2015_certolizumab convention of pairing the two indicators as mutually exclusive). It includes Chinese patients and the other Asian groups. Note the difference from the companion PK model Johnson_2025_osimertinib, whose RACE_ASIAN_OTH excludes Chinese patients as well. Retained in the final model although its likelihood-ratio gain was small (decrease in OBJ of 1.5 vs Model 2; Johnson 2024 Results).",
      source_name = "Asian (non-Japanese)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline.",
      units = "years",
      type = "continuous",
      notes = "Screened as a single covariate added to the final model; not retained (decrease in OBJ 0.362, p = 0.547; Johnson 2024 Table S9). Weight (p = 0.829), body mass index (p = 0.516), sex (p = 0.458), WHO performance status 0 vs >= 1 (p = 0.313), history of hypertension (p = 0.948), line of therapy first vs later (p = 0.673), nicotine use (p = 0.254) and statin use (p = 0.278) were screened the same way and also dropped; no point estimates are reported."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1364L,
    n_studies = 5L,
    age_range = "25-91 years (median 62; Johnson 2024 Table S3)",
    weight_range = "29-122 kg (median 60-62; Johnson 2024 Table S3)",
    sex_female_pct = 64.7,
    race_ethnicity = c(
      Japanese = 18.6,
      Chinese = 20.3,
      `Asian non-Chinese non-Japanese` = 23.7,
      White = 27.5,
      Other = 9.9
    ),
    disease_state = "EGFR mutation-positive advanced non-small cell lung cancer (first-line, or T790M-positive after prior EGFR-TKI therapy)",
    dose_range = "osimertinib 20-240 mg once daily (AURA) and 80 mg once daily (AURA extension, AURA2, AURA3, FLAURA)",
    regions = "Multinational; 62% Asian",
    notes = "Same analysis set as Johnson_2024_osimertinib_ild_osimertinib: studies AURA, AURA extension, AURA2, AURA3 and FLAURA (Johnson 2024 Table S1). n_subjects is the osimertinib-treated first-line (338) plus second-or-later-line (1026) population; the fit used the subset with a computable AUCss, whose size Johnson 2024 does not state. Percentages are the first-line plus second-or-later-line columns of Table S3. The final model was used to predict ILD in the adjuvant ADAURA study (N = 325) as an external check (Figure 2E)."
  )

  ini({
    # Johnson 2024 Table S7, Model 3 (final AZ5104 ILD model); same
    # functional form as the osimertinib model of Table 1, with C_i the
    # natural log of the AZ5104 AUCss in nM*h.
    logit_ref           <- -7.46;  label("Logit of the ILD probability for a non-Asian patient at AUC_AZ5104 = 1 nmol*h/L (unitless logit)")    # Johnson 2024 Table S7, Model 3: Intercept -7.46 (95% CI -10.6, -4.33)
    e_auc_ild           <- 0.559;  label("Log-odds of ILD per e-fold increase in AZ5104 AUCss (unitless logit)")                                 # Johnson 2024 Table S7, Model 3: Exposure 0.559 (95% CI 0.135, 0.969)
    e_race_japanese_ild <- 1.17;   label("Log-odds of ILD for Japanese vs non-Asian patients (unitless logit)")                                  # Johnson 2024 Table S7, Model 3: Japanese 1.17 (95% CI 0.541, 1.83)
    e_race_asian_ild    <- -0.551; label("Log-odds of ILD for Asian non-Japanese vs non-Asian patients (unitless logit)")                        # Johnson 2024 Table S7, Model 3: Asian (non-Japanese) -0.551 (95% CI -1.37, 0.223)

    # Binary landmark endpoint with a penalized Bernoulli likelihood and no
    # random effect. The placeholder additive residual exists only so that
    # rxode2 accepts an observation declaration.
    addSd_prob_ild <- fixed(0.001); label("Placeholder additive residual SD on the ILD probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # Natural-log AUCss, uncentred: logit_ref is the logit at 1 nM*h, far
    # below the observed range. Johnson 2024 quotes expit(logit_ref) =
    # 0.06% as its 'AUCss = 0' placebo extrapolation.
    logit_ild <- logit_ref + e_auc_ild * log(AUC_AZ5104) +
      e_race_japanese_ild * RACE_JAPANESE + e_race_asian_ild * RACE_ASIAN
    prob_ild <- expit(logit_ild)

    prob_ild ~ add(addSd_prob_ild)
  })
}
