Johnson_2024_osimertinib_ild_osimertinib <- function() {
  description <- paste0(
    "Penalized (Firth) logistic-regression exposure-safety model for ",
    "interstitial lung disease (ILD) and ILD-like events in patients ",
    "with EGFR-mutated non-small cell lung cancer treated with ",
    "osimertinib in the first-line or second-or-later-line setting ",
    "(Johnson 2024, final Model 3 of Table 1). The probability of an ",
    "ILD-like event is expit(a + b * log(AUC_OSIMERTINIB) + ",
    "c * RACE_JAPANESE + d * RACE_ASIAN), where AUC_OSIMERTINIB is the ",
    "individual steady-state osimertinib AUC in nmol*h/L (dose divided ",
    "by the individual apparent clearance of the companion population ",
    "PK model Johnson_2025_osimertinib) and the race indicators are ",
    "Japanese and Asian non-Japanese against a non-Asian reference. ",
    "There is no PK layer and no ODE: exposure is supplied as data. ",
    "A sibling model driven by the AZ5104 metabolite AUC is ",
    "Johnson_2024_osimertinib_ild_az5104."
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
    dosing = "n/a (no dose events; exposure enters as the AUC_OSIMERTINIB covariate column)",
    concentration = "prob_ild (probability of an ILD or ILD-like event, 0-1)"
  )

  covariateData <- list(
    AUC_OSIMERTINIB = list(
      description = "Individual steady-state osimertinib AUC over the once-daily dosing interval, AUCss = dose / (CL/F)_i, supplied as data (this model has no PK layer).",
      units = "nmol*h/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Johnson 2024 Methods 'Estimation of plasma exposure': AUCss = dose/clearance from the individual apparent clearances of the population PK analysis (reference 23 = Johnson 2025, modellib('Johnson_2025_osimertinib')); for the safety analysis the FIRST dose a patient received was used. The unit is nM*h, the axis label of Johnson 2024 Figure 2A and 2C; Johnson_2025_osimertinib returns cl in L/h, so AUC_OSIMERTINIB = dose_mg / cl * 1e6 / 499.62. Enters the logit on the natural-log scale, uncentred, so a zero AUC (placebo) gives log(0) = -Inf and a probability of exactly 0. Johnson 2024 Figure 2A quartile boundaries of the analysis set: 6584.57, 9689.56, 12910.9 and 21683.3 nM*h.",
      source_name = "AUCss (osimertinib)"
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
      notes = "In this model RACE_ASIAN EXCLUDES Japanese patients, who carry RACE_JAPANESE instead (the Wade_2015_certolizumab convention of pairing the two indicators as mutually exclusive). It includes Chinese patients and the other Asian groups (Johnson 2024 Table S2 grouped race: Chinese 21.7% and Asian non-Chinese non-Japanese 23.7% of the osimertinib-treated patients). Note the difference from the companion PK model Johnson_2025_osimertinib, whose RACE_ASIAN_OTH excludes Chinese patients as well. Retained in the final model although its likelihood-ratio gain was small (decrease in OBJ of 1.6 vs Model 2; Johnson 2024 Results).",
      source_name = "Asian (non-Japanese)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline.",
      units = "years",
      type = "continuous",
      notes = "Screened as a single covariate added to the final model; not retained (decrease in OBJ 0.397, p = 0.529; Johnson 2024 Table S8). Weight (p = 0.660), body mass index (p = 0.389), sex (p = 0.495), WHO performance status 0 vs >= 1 (p = 0.273), history of hypertension (p = 0.989), line of therapy first vs later (p = 0.701), nicotine use (p = 0.286) and statin use (p = 0.258) were screened the same way and also dropped; no point estimates are reported."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1364L,
    n_studies = 5L,
    n_events = "51 ILD-like events across the four osimertinib AUCss quartiles of Johnson 2024 Figure 2A (6, 12, 10 and 23)",
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
    notes = "Studies AURA, AURA extension, AURA2, AURA3 and FLAURA (Johnson 2024 Table S1). n_subjects is the osimertinib-treated first-line (338) plus second-or-later-line (1026) population; the fit used the subset with a computable AUCss, whose size Johnson 2024 does not state. Percentages are the first-line plus second-or-later-line columns of Table S3. The final model was then used to predict ILD in the adjuvant ADAURA study (N = 325) as an external check (Figure 2E). Comparator EGFR-TKI and chemotherapy patients were not included in the fit."
  )

  ini({
    # Johnson 2024 Table 1, Model 3 (final osimertinib ILD model).
    # Supplementary Methods: log(p / (1 - p)) = alpha + beta_0 * C_i plus
    # additive categorical-covariate terms. C_i is the natural log of the
    # AUCss in nM*h: the Model 3 estimates reproduce the predicted-mean
    # curves of Figure 2C on that scale (see the vignette), and the
    # LVEF analysis states the same 'log-transformed' AUCss.
    logit_ref           <- -9.35;  label("Logit of the ILD probability for a non-Asian patient at AUC_OSIMERTINIB = 1 nmol*h/L (unitless logit)")    # Johnson 2024 Table 1, Model 3: Intercept -9.35 (95% CI -13.8, -4.76)
    e_auc_ild           <- 0.628;  label("Log-odds of ILD per e-fold increase in osimertinib AUCss (unitless logit)")                                 # Johnson 2024 Table 1, Model 3: Exposure 0.628 (95% CI 0.149, 1.09)
    e_race_japanese_ild <- 1.14;   label("Log-odds of ILD for Japanese vs non-Asian patients (unitless logit)")                                       # Johnson 2024 Table 1, Model 3: Japanese 1.14 (95% CI 0.510, 1.80)
    e_race_asian_ild    <- -0.567; label("Log-odds of ILD for Asian non-Japanese vs non-Asian patients (unitless logit)")                             # Johnson 2024 Table 1, Model 3: Asian (non-Japanese) -0.567 (95% CI -1.38, 0.206)

    # Binary landmark endpoint with a penalized Bernoulli likelihood and no
    # random effect. The placeholder additive residual exists only so that
    # rxode2 accepts an observation declaration.
    addSd_prob_ild <- fixed(0.001); label("Placeholder additive residual SD on the ILD probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # Natural-log AUCss, uncentred: logit_ref is the logit at 1 nM*h, far
    # below the observed range. Johnson 2024 quotes expit(logit_ref) =
    # 0.009% as its 'AUCss = 0' placebo extrapolation.
    logit_ild <- logit_ref + e_auc_ild * log(AUC_OSIMERTINIB) +
      e_race_japanese_ild * RACE_JAPANESE + e_race_asian_ild * RACE_ASIAN
    prob_ild <- expit(logit_ild)

    prob_ild ~ add(addSd_prob_ild)
  })
}
