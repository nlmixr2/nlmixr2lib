Sangro_2023_nivolumab_ipilimumab_os <- function() {
  description <- paste0(
    "Cox proportional-hazards exposure-efficacy model for OVERALL ",
    "SURVIVAL in sorafenib-experienced adults with advanced ",
    "hepatocellular carcinoma receiving nivolumab monotherapy ",
    "(0.1-10 mg/kg every 2 weeks) or nivolumab plus ipilimumab ",
    "(Sangro 2023, n = 328, cohorts 1, 2 and 4 of CheckMate 040). The ",
    "model returns the RELATIVE hazard exp(0.07747 * log(CAV_NIVOLUMAB) ",
    "- 0.1539 * log(CAV_IPILIMUMAB) + 0.4444 * ECOG_GE1 + 0.05068 * CL ",
    "+ 0.7962 * log(LDHR)), where CAV_NIVOLUMAB and CAV_IPILIMUMAB are ",
    "the model-predicted average serum concentrations over the FIRST ",
    "dosing interval (Cavg1, ug/mL), the ipilimumab log term is 0 for ",
    "patients who received no ipilimumab, and CL is the baseline posthoc ",
    "nivolumab clearance in mL/h (supplied as CLI in L/h). Higher ",
    "ipilimumab exposure lowers the hazard (HR 0.857 per unit natural-log ",
    "Cavg1); the nivolumab exposure term is not significant (HR 1.081). ",
    "NO BASELINE HAZARD IS ENCODED: a Cox regression leaves h0(t) ",
    "unspecified, so only ratios of hr between two covariate settings ",
    "are meaningful. There is no PK layer and no ODE: exposures and ",
    "clearance are supplied as data, derived in the source from updated ",
    "versions of the nivolumab (Bajaj 2017) and ipilimumab (Sanghavi ",
    "2020) population PK models. No random effect and no residual error ",
    "are estimated."
  )
  reference <- paste(
    "Sangro B, Yau T, El-Khoueiry AB, Kudo M, Shen Y, Tschaika M, Roy A,",
    "Feng Y, Gao L, Aras U. Exposure-response analysis for nivolumab plus",
    "ipilimumab combination therapy in patients with advanced",
    "hepatocellular carcinoma (CheckMate 040). Clin Transl Sci.",
    "2023;16(8):1445-1457. doi:10.1111/cts.13544. PMCID PMC10432868.",
    "Coefficients are transcribed from Table 3. The first-dose exposures",
    "and baseline nivolumab clearance came from updated versions of the",
    "population PK models packaged as modellib('Bajaj_2017_nivolumab') and",
    "modellib('Sanghavi_2020_ipilimumab').",
    sep = " "
  )
  vignette <- "Sangro_2023_nivolumab_ipilimumab_hcc"

  units <- list(
    time = "n/a (semiparametric Cox regression; the baseline hazard and hence the time scale are left unspecified by the method)",
    dosing = "n/a (no dose events; exposures enter as the CAV_NIVOLUMAB and CAV_IPILIMUMAB columns)",
    concentration = "hr (relative hazard of death, unitless; also lhr, the linear predictor)"
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
        "for the every-3-weeks induction of arms A and B (Sangro 2023",
        "Methods; empirical Bayes predictions from the updated nivolumab",
        "population PK model under the assigned regimen). Enters on the",
        "NATURAL-LOG scale, uncentred (Table 3 row 'Log-transformed",
        "nivolumab Cavg1, ug/mL'). Every patient received nivolumab, so",
        "the column is always positive. Table S2 geometric means: 3.8",
        "(nivolumab < 3 mg/kg), 24.8 (3 mg/kg every 2 weeks), 86.8",
        "(10 mg/kg), 6.5 (arm A), 21.5 (arm B), 23.0 (arm C) ug/mL."
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
        "for nivolumab monotherapy. Enters on the NATURAL-LOG scale",
        "(Table 3 row 'Log-transformed ipilimumab Cavg1, ug/mL'), and",
        "the log term is set to 0 when CAV_IPILIMUMAB = 0. Neither the",
        "zero handling nor the log base is printed; both are pinned by",
        "Figure 3 (monotherapy hazard ratio about 1; natural log puts",
        "the arm A geometric-mean hazard ratio at 0.58, inside the",
        "plotted 0.56-0.60); see the vignette."
      ),
      source_name = "Ipilimumab Cavg1"
    ),
    ECOG_GE1 = list(
      description = "Baseline ECOG performance status indicator (1 = ECOG 1, 0 = ECOG 0)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ECOG 0)",
      notes = paste(
        "Table 3 row 'ECOG performance status (>0)': 0.4444 (HR 1.56).",
        "CheckMate 040 enrolled ECOG 0 or 1 only, so ECOG > 0 is",
        "ECOG = 1 (40% of the analysis population, Table 1)."
      ),
      source_name = "ECOG performance status"
    ),
    CLI = list(
      description = "Baseline (first-dose) posthoc nivolumab clearance from the updated nivolumab population PK model",
      units = "L/h",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Supply in the canonical L/h; the model multiplies by 1000",
        "because the Table 3 coefficient is per mL/h ('Baseline",
        "nivolumab CL, mL/h': 0.05068, HR 1.052 per mL/h). Linear and",
        "uncentred. Derived in the source from the empirical Bayes",
        "estimates of the updated nivolumab population PK model (Bajaj",
        "2017 lineage); the Bajaj 2017 reference-patient baseline CL is",
        "9.4 mL/h (0.0094 L/h). Populate it by simulating",
        "modellib('Bajaj_2017_nivolumab') and taking the time-zero",
        "clearance, or set every subject to the typical value for a",
        "typical-patient calculation."
      ),
      source_name = "Baseline nivolumab CL"
    ),
    LDHR = list(
      description = "Baseline serum lactate dehydrogenase divided by the laboratory upper limit of normal",
      units = "(unitless ratio)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters as its natural log (Table 3 row 'Log ratio baseline",
        "LDH/ULN': 0.7962, HR 2.217 per unit log ratio), so a patient at",
        "the upper limit of normal (LDHR = 1) contributes nothing. The",
        "log base is not printed; natural log is assumed, matching the",
        "paper's other 'log-transformed' terms whose base the figures",
        "pin."
      ),
      source_name = "Log ratio baseline LDH/ULN"
    )
  )

  covariatesDataExcluded <- list(
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      notes = paste(
        "Screened in the overall-survival Cox model but not retained",
        "(Sangro 2023 Methods lists it among the candidate covariates;",
        "the final model in Table 3 does not carry it). No point",
        "estimate is reported. The same applies to HCC etiology,",
        "baseline alpha-fetoprotein and baseline tumor size."
      )
    ),
    TUMSZ = list(
      description = "Baseline tumor size",
      units = "mm",
      type = "continuous",
      notes = "Screened in the overall-survival model but not retained (Sangro 2023 Methods and Table 3); no point estimate reported."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 328L,
    n_studies = 1L,
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
      "assignment from Table S1. Minimum follow-up 28 months (cohort 4)",
      "and 31 / 27 months (cohorts 1 / 2); data cutoff March 2019.",
      "ECOG performance status 0 in 60% and 1 in 40%. The model-predicted",
      "survival agreed with observed Kaplan-Meier curves (Figure S1)",
      "except for a deviation after week 12 that the authors attribute",
      "to the first-dose ipilimumab exposure under-representing",
      "continuous ipilimumab dosing in arm C."
    )
  )

  ini({
    # ==================================================================
    # Sangro 2023 Table 3, 'Parameter estimates of exposure-response for
    # overall survival'. Every row is printed as a coefficient on the
    # log-hazard scale with its SE and hazard ratio, so nothing here is
    # back-solved or digitised.
    #
    #   Log-transformed nivolumab Cavg1    0.07747  SE 0.1004   HR 1.081
    #   Log-transformed ipilimumab Cavg1  -0.1539   SE 0.06909  HR 0.8574
    #   ECOG performance status (>0)       0.4444   SE 0.1400   HR 1.56
    #   Baseline nivolumab CL, mL/h        0.05068  SE 0.01682  HR 1.052
    #   Log ratio baseline LDH/ULN         0.7962   SE 0.2047   HR 2.217
    #
    # NO BASELINE HAZARD IS ENCODED. A Cox regression is semiparametric:
    # h0(t) is not an unreported parameter but a quantity the fit never
    # produced. The model returns the relative hazard hr only, matching
    # Chan_2025_atezolizumab_os.R. There is no IIV, no residual error
    # and no observation endpoint.
    # ==================================================================
    e_cav_nivolumab_haz <- 0.07747; label("Log hazard ratio for death per unit natural-log nivolumab first-dose Cavg1 (log scale)") # Sangro 2023 Table 3, Log-transformed nivolumab Cavg1 row: 0.07747 (SE 0.1004; HR 1.081)
    e_cav_ipilimumab_haz <- -0.1539; label("Log hazard ratio for death per unit natural-log ipilimumab first-dose Cavg1 (log scale)") # Sangro 2023 Table 3, Log-transformed ipilimumab Cavg1 row: -0.1539 (SE 0.06909; HR 0.8574)
    e_ecog_ge1_haz <- 0.4444; label("Log hazard ratio for death for ECOG performance status 1 versus 0 (log scale)") # Sangro 2023 Table 3, ECOG performance status (>0) row: 0.4444 (SE 0.1400; HR 1.56)
    e_cli_haz <- 0.05068; label("Log hazard ratio for death per mL/h of baseline nivolumab clearance (log scale)") # Sangro 2023 Table 3, Baseline nivolumab CL (mL/h) row: 0.05068 (SE 0.01682; HR 1.052)
    e_ldhr_haz <- 0.7962; label("Log hazard ratio for death per unit natural-log of baseline LDH / ULN (log scale)") # Sangro 2023 Table 3, Log ratio baseline LDH/ULN row: 0.7962 (SE 0.2047; HR 2.217)
  })

  model({
    # The ipilimumab exposure term is the natural log of Cavg1 for
    # patients who received ipilimumab and 0 for nivolumab monotherapy
    # (CAV_IPILIMUMAB = 0), where the log is undefined. The source
    # prints neither choice; Figure 3 plots each arm's hazard ratio against a
    # nivolumab-only reference, which this handling reproduces (see the
    # vignette for the arm-by-arm check).
    log_cav_ipilimumab <- 0
    if (CAV_IPILIMUMAB > 0) log_cav_ipilimumab <- log(CAV_IPILIMUMAB)

    # CLI is supplied in L/h; Table 3 reports the coefficient per mL/h.
    cli_mlh <- CLI * 1000

    # Linear predictor and relative hazard against a hypothetical patient
    # with nivolumab Cavg1 = 1 ug/mL, no ipilimumab, ECOG 0, clearance 0
    # and LDH at the ULN. Only ratios of hr between two covariate
    # settings are meaningful; under proportional hazards the subject
    # hazard is h0(t) * hr for any baseline h0(t) the user supplies.
    lhr <- e_cav_nivolumab_haz * log(CAV_NIVOLUMAB) +
      e_cav_ipilimumab_haz * log_cav_ipilimumab +
      e_ecog_ge1_haz * ECOG_GE1 +
      e_cli_haz * cli_mlh +
      e_ldhr_haz * log(LDHR)

    hr <- exp(lhr)
  })
}
