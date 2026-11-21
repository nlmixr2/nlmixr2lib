Sangro_2023_nivolumab_ipilimumab_orr <- function() {
  description <- paste0(
    "Logistic-regression exposure-efficacy model for objective tumor ",
    "response (complete or partial response per RECIST 1.1, by blinded ",
    "independent central review) in sorafenib-experienced adults with ",
    "advanced hepatocellular carcinoma receiving nivolumab monotherapy ",
    "(0.1-10 mg/kg every 2 weeks) or nivolumab plus ipilimumab ",
    "(Sangro 2023, n = 328, cohorts 1, 2 and 4 of CheckMate 040). The ",
    "probability of response is expit(-1.419 - 0.00841 * CAV_NIVOLUMAB ",
    "+ 0.3712 * log(CAV_IPILIMUMAB)), where both columns are the ",
    "model-predicted average serum concentrations over the FIRST dosing ",
    "interval (Cavg1, ug/mL) and the ipilimumab log term is 0 for ",
    "patients who received no ipilimumab (CAV_IPILIMUMAB = 0). Nivolumab ",
    "exposure enters linearly and is not significant (odds ratio 0.9916 ",
    "per ug/mL); ipilimumab exposure is the response driver (odds ratio ",
    "1.45 per unit natural-log Cavg1). No other covariate was retained. ",
    "There is no PK layer and no ODE: the exposures are supplied as data ",
    "and were derived in the source from updated versions of the ",
    "nivolumab (Bajaj 2017) and ipilimumab (Sanghavi 2020) population PK ",
    "models. No random effect and no residual error are estimated."
  )
  reference <- paste(
    "Sangro B, Yau T, El-Khoueiry AB, Kudo M, Shen Y, Tschaika M, Roy A,",
    "Feng Y, Gao L, Aras U. Exposure-response analysis for nivolumab plus",
    "ipilimumab combination therapy in patients with advanced",
    "hepatocellular carcinoma (CheckMate 040). Clin Transl Sci.",
    "2023;16(8):1445-1457. doi:10.1111/cts.13544. PMCID PMC10432868.",
    "Coefficients are transcribed from Table 2. The first-dose exposures",
    "came from updated versions of the population PK models packaged as",
    "modellib('Bajaj_2017_nivolumab') and",
    "modellib('Sanghavi_2020_ipilimumab').",
    sep = " "
  )
  vignette <- "Sangro_2023_nivolumab_ipilimumab_hcc"

  units <- list(
    time = "n/a (static landmark exposure-response model; no time dimension)",
    dosing = "n/a (no dose events; exposures enter as the CAV_NIVOLUMAB and CAV_IPILIMUMAB columns)",
    concentration = "prob_orr_central (probability of objective tumor response by blinded independent central review, 0-1)"
  )

  covariateData <- list(
    CAV_NIVOLUMAB = list(
      description = "Model-predicted nivolumab average serum concentration over the first dosing interval (Cavg1)",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Averaging window: the FIRST dosing interval, i.e. AUC over",
        "days 0-14 / 14 for every-2-weeks nivolumab and days 0-21 / 21",
        "for the every-3-weeks induction of arms A and B. Sangro 2023",
        "Methods: the updated population PK models 'were applied to",
        "derive nivolumab and ipilimumab time-averaged concentrations",
        "following the first dose (Cavg1) ... based upon the assigned",
        "treatment regimen', using empirical Bayes estimates. Enters",
        "LINEARLY and UNCENTRED (Table 2 row 'Nivolumab Cavg1, ug/mL';",
        "the Results state that linear-scale nivolumab exposure was",
        "retained after log and linear forms were compared). Table S2",
        "geometric means: 3.8 (nivolumab < 3 mg/kg), 24.8 (3 mg/kg",
        "every 2 weeks), 86.8 (10 mg/kg), 6.5 (arm A), 21.5 (arm B),",
        "23.0 (arm C) ug/mL."
      ),
      source_name = "Nivolumab Cavg1"
    ),
    CAV_IPILIMUMAB = list(
      description = "Model-predicted ipilimumab average serum concentration over the first dosing interval (Cavg1); 0 for nivolumab monotherapy",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Averaging window: the FIRST dosing interval, i.e. AUC over",
        "days 0-21 / 21 for every-3-weeks ipilimumab (arms A and B) and",
        "days 0-42 / 42 for every-6-weeks ipilimumab (arm C). Set to 0",
        "for patients on nivolumab monotherapy (Table S2 lists 0.0 for",
        "every monotherapy group). Enters on the NATURAL-LOG scale",
        "(Table 2 row 'Log-transformed ipilimumab Cavg1, ug/mL'), and",
        "the log term is set to 0 when CAV_IPILIMUMAB = 0. Neither the",
        "zero handling nor the log base is printed; both are pinned by",
        "Figure 2 (monotherapy odds ratio about 1; natural log puts the",
        "arm A geometric-mean odds ratio at 3.38, inside the plotted",
        "2.92-3.70, where log10 would give 1.85); see the vignette.",
        "Geometric means: 16.7 / 5.1 / 1.7 ug/mL in arms A / B / C",
        "(Results text) or 17.5 / 5.8 / 2.3 ug/mL (Table S2)."
      ),
      source_name = "Ipilimumab Cavg1"
    )
  )

  covariatesDataExcluded <- list(
    ECOG_GE1 = list(
      description = "Baseline ECOG performance status indicator (1 = ECOG 1, 0 = ECOG 0)",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened in the objective-response logistic model but not",
        "retained: 'No other covariates were found to have a significant",
        "effect on OTR; therefore, they were not included in the final",
        "model' (Sangro 2023 Results). No point estimate is reported. The",
        "same applies to the other screened covariates: HCC etiology,",
        "baseline tumor-cell PD-L1 status (< 1% vs >= 1%), baseline",
        "tumor size, extrahepatic spread or vascular invasion, and",
        "baseline nivolumab clearance."
      )
    ),
    PDL1_TUM_POS = list(
      description = "Baseline tumor-cell PD-L1 positivity indicator (1 = >= 1%, 0 = < 1%)",
      units = "(binary)",
      type = "binary",
      notes = "Screened in the objective-response model but not retained (Sangro 2023 Methods and Results); no point estimate reported."
    ),
    TUMSZ = list(
      description = "Baseline tumor size",
      units = "mm",
      type = "continuous",
      notes = "Screened in the objective-response model but not retained (Sangro 2023 Methods and Results); no point estimate reported."
    ),
    CLI = list(
      description = "Baseline posthoc nivolumab clearance",
      units = "L/h",
      type = "continuous",
      notes = "Screened in the objective-response model but not retained (Sangro 2023 Results). It IS retained in the companion overall-survival model, Sangro_2023_nivolumab_ipilimumab_os."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 328L,
    n_studies = 1L,
    n_observations = "328 binary best-overall-response records (one per patient)",
    age_range = "18-83 years (median 62)",
    weight_range = "38-130 kg (median 67)",
    sex_female_pct = 21,
    race_ethnicity = c(White = 41, Asian = 46, Black = 3, Other = 11),
    disease_state = "advanced hepatocellular carcinoma, sorafenib-experienced (progressed on or intolerant to sorafenib); Child-Pugh A (98%) or B; etiology HBV 41%, HCV 20%, HBV/HCV coinfected 2%, uninfected 37%; extrahepatic spread or vascular invasion 86%",
    dose_range = paste(
      "Nivolumab monotherapy 0.1, 0.3, 1, 3 or 10 mg/kg every 2 weeks",
      "(n = 182; cohorts 1 and 2); arm A nivolumab 1 mg/kg + ipilimumab",
      "3 mg/kg every 3 weeks x 4 then nivolumab 240 mg every 2 weeks",
      "(n = 49); arm B nivolumab 3 mg/kg + ipilimumab 1 mg/kg every",
      "3 weeks x 4 then nivolumab 240 mg every 2 weeks (n = 49); arm C",
      "nivolumab 3 mg/kg every 2 weeks + ipilimumab 1 mg/kg every",
      "6 weeks (n = 48)"
    ),
    regions = "CheckMate 040 (NCT01658878) international sites (Asia, Europe, North America)",
    notes = paste(
      "Baseline demographics from Sangro 2023 Table 1; treatment",
      "assignment from Table S1; observed objective response by group",
      "from Table S2 (0.22, 0.14 and 0 for nivolumab < 3, 3 and",
      "10 mg/kg monotherapy; 0.33, 0.31 and 0.31 for arms A, B and C).",
      "ECOG performance status 0 in 60% and 1 in 40%. Data cutoff March",
      "2019. A post hoc exploratory analysis; the monotherapy and",
      "combination cohorts were not randomised against each other."
    )
  )

  ini({
    # ==================================================================
    # Sangro 2023 Table 2, 'Parameter estimates of exposure-response for
    # BICR-assessed objective tumor response'. Every row is printed as
    # an estimate on the logit scale with its SE and odds ratio, so
    # nothing here is back-solved or digitised.
    #
    #   Intercept                          -1.419    SE 0.3522  OR 0.242
    #   Nivolumab Cavg1, ug/mL             -0.00841  SE 0.01229 OR 0.9916
    #   Log-transformed ipilimumab Cavg1    0.3712   SE 0.1263  OR 1.45
    #
    # The logistic regression has an exact Bernoulli likelihood with no
    # random effect and no residual error. rxode2 requires an observation
    # declaration, so the deterministic probability is emitted with a
    # tiny placeholder additive residual, mirroring
    # Babel_2026_telisotuzumab_orr.R.
    # ==================================================================
    logit_ref <- -1.419; label("Logit of objective response with nivolumab Cavg1 = 0 and no ipilimumab log term (unitless logit)") # Sangro 2023 Table 2, Intercept: -1.419 (SE 0.3522; OR 0.242)
    e_cav_nivolumab_logit <- -0.00841; label("Log-odds of objective response per ug/mL of nivolumab first-dose Cavg1 (unitless logit per ug/mL)") # Sangro 2023 Table 2, Nivolumab Cavg1 row: -0.00841 (SE 0.01229; OR 0.9916)
    e_cav_ipilimumab_logit <- 0.3712; label("Log-odds of objective response per unit natural-log ipilimumab first-dose Cavg1 (unitless logit)") # Sangro 2023 Table 2, Log-transformed ipilimumab Cavg1 row: 0.3712 (SE 0.1263; OR 1.45)

    addSd_prob_orr_central <- fixed(0.001); label("Placeholder additive residual SD on the response probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # The ipilimumab exposure term is the natural log of Cavg1 for
    # patients who received ipilimumab and 0 for nivolumab monotherapy
    # (CAV_IPILIMUMAB = 0), where the log is undefined. The source
    # prints neither choice; Figure 2 plots each arm's odds ratio against a
    # nivolumab-only reference, which this handling reproduces (see the
    # vignette for the arm-by-arm check).
    log_cav_ipilimumab <- 0
    if (CAV_IPILIMUMAB > 0) log_cav_ipilimumab <- log(CAV_IPILIMUMAB)

    logit_orr_central <- logit_ref +
      e_cav_nivolumab_logit * CAV_NIVOLUMAB +
      e_cav_ipilimumab_logit * log_cav_ipilimumab
    prob_orr_central <- expit(logit_orr_central)

    # Deterministic probability of objective response by blinded
    # independent central review. Sample binary outcomes with
    # rbinom(n, 1, prob_orr_central) on the rxSolve output.
    prob_orr_central ~ add(addSd_prob_orr_central)
  })
}
