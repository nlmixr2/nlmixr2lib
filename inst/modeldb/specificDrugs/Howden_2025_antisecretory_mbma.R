Howden_2025_antisecretory_mbma <- function() {
  description <- paste(
    "MBMA. Study-arm-level sigmoid Emax model on the logit scale linking the",
    "steady-state percentage of 24 h with intragastric pH > 4 (holding-time",
    "ratio, HTR) of an antisecretory regimen (histamine H2-receptor",
    "antagonist, proton pump inhibitor or potassium-competitive acid blocker)",
    "to the erosive esophagitis healing rate at week 4 or week 8, pooled from",
    "82 pH-metry studies and 104 randomized healing trials (Howden 2025).",
    "Between-drug and between-treatment random effects are scaled by each",
    "arm's standard error of the logit healing rate; the model simulates",
    "arm-level healing rates, not individual patients."
  )
  reference <- paste(
    "Howden CW, Scarpignato C, Vaezi MF, Jacob R, Wolle M, Hunt RH.",
    "Mathematical model of the relationship between pH holding time and",
    "erosive esophagitis healing rates. CPT Pharmacometrics Syst Pharmacol.",
    "2025;14(1):28-35. doi:10.1002/psp4.13235"
  )
  vignette <- "Howden_2025_antisecretory_mbma"
  units <- list(
    time = "week",
    dosing = "none",
    concentration = paste(
      "percent (arm erosive esophagitis healing rate, exposed as ee_healing;",
      "the model observation logit_ee_healing is that rate on the",
      "logit-of-fraction scale and is NOT a drug concentration)"
    )
  )

  covariateData <- list(
    GASTRIC_PH4_TIME_PCT = list(
      description = paste(
        "Steady-state (days 5-8 of daily dosing) mean percentage of a 24 h",
        "period with intragastric pH > 4 (pH > 4 holding-time ratio) for the",
        "arm's drug and dose regimen, the model's exposure driver"
      ),
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Regimen-level, not measured in the healing-trial patients: the",
        "authors pooled 24 h pH-metry studies in healthy volunteers by drug",
        "and dose (weighted mean, weights = per-protocol size) and merged the",
        "pooled value onto every healing-trial arm of the same regimen",
        "(Methods 'Data analysis'). Range 10.2% (placebo) to 96.6%",
        "(vonoprazan 40 mg OD) in Table 1. Holding time in hours =",
        "GASTRIC_PH4_TIME_PCT / 100 * 24."
      ),
      source_name = "PH4"
    ),
    SE_ARM = list(
      description = paste(
        "Standard error of the arm's pooled erosive esophagitis healing rate",
        "on the logit scale"
      ),
      units = "(unitless, logit scale)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Supplied per observation. Multiplies BOTH the between-drug random",
        "effect and the between-treatment (residual) random effect",
        "(Supplementary 'Final Model Code': IPRED = MU + EEHR_SE*BDV_ +",
        "EEHR_SE*BTV_). The paper's simulations used 0.25 (high precision,",
        "N about 300), 0.5 (average, about 100 patients) and 1.0 (low",
        "precision)."
      ),
      source_name = "EEHR_SE"
    ),
    TRT_APPROVED = list(
      description = paste(
        "Regimen-approval indicator: 1 = the drug and dose regimen is",
        "approved for erosive esophagitis in the US or the European Union,",
        "0 = not approved"
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "1 (approved)",
      notes = paste(
        "Raises EC50 (by -46.5%, i.e. lowers it) when 0 and multiplies the",
        "between-treatment SD by 1 + 2.05 when 1 (Supplementary 'Final Model",
        "Code' APPR). Non-approved regimens are footnoted in Table 1."
      ),
      source_name = "APPR"
    ),
    TRT_H2RA = list(
      description = paste(
        "Drug-class indicator: 1 = histamine H2-receptor antagonist arm",
        "(famotidine, nizatidine, ranitidine), 0 otherwise"
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (proton pump inhibitor or placebo)",
      notes = paste(
        "NONMEM CLASS = 3. Multiplies EC50 by 1 + 0.413 and the",
        "between-treatment SD by 1 - 0.71. Mutually exclusive with TRT_PCAB."
      ),
      source_name = "CLASS"
    ),
    TRT_PCAB = list(
      description = paste(
        "Drug-class indicator: 1 = potassium-competitive acid blocker arm",
        "(vonoprazan, tegoprazan), 0 otherwise"
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (proton pump inhibitor or placebo)",
      notes = paste(
        "NONMEM CLASS = 5. Multiplies EC50 by 1 - 0.308. Mutually exclusive",
        "with TRT_H2RA."
      ),
      source_name = "CLASS"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 43417,
    n_studies = 186,
    age_range = ">= 18 years (literature-search restriction)",
    disease_state = paste(
      "Erosive esophagitis (healing RCTs, 43,417 patients in 196 arms);",
      "intragastric pH data from healthy volunteers (82 phase I studies,",
      "4297 participants in 201 arms)"
    ),
    dose_range = paste(
      "Placebo; famotidine 20 mg BID, nizatidine 150 mg BID, ranitidine",
      "150 mg QID / 300 mg OD; esomeprazole, lansoprazole, omeprazole,",
      "pantoprazole, rabeprazole, azeloprazole, dexlansoprazole and",
      "ilaprazole at 10-90 mg OD or BID; tegoprazan 50-100 mg OD;",
      "vonoprazan 10-40 mg OD (Table 1)"
    ),
    regions = "Not reported (published literature, MEDLINE and Embase to 31 March 2023)",
    notes = paste(
      "Model-based meta-analysis of published study-arm summaries; the",
      "observation is an arm's pooled healing rate at week 4 or week 8",
      "(intention-to-treat). Approved and non-approved regimens were both",
      "used in model development (Table 1). n_subjects counts the",
      "healing-trial patients; n_studies counts both datasets (82 pH",
      "studies + 104 healing RCTs)."
    )
  )

  ini({
    # Supplementary 'Final Model Code' (NONMEM $PRED) and Table S2. Values
    # are the Table S2 final estimates; the $THETA block of the printed code
    # carries a different set of numbers (e.g. EEHR0 10.32, EPH50 24.74),
    # which reproduce neither Table S3 nor the Results probabilities.
    logite0 <- log(0.123 / (1 - 0.123))
    label("Logit of the placebo (HTR = 0) erosive esophagitis healing rate (unitless)") # Table S2 E0 = 12.3% (RSE 68.9%)
    logitemax <- log(0.959 / (1 - 0.959))
    label("Logit of the maximum erosive esophagitis healing rate (unitless)") # Table S2 Emax = 95.9% (RSE 2.3%)
    lec50_wk8 <- log(26.9)
    label("Log of the pH > 4 holding-time ratio giving half-maximal effect at week 8 (log %)") # Table S2 E_HTR,50 Week 8 = 26.9% (RSE 21.6%)
    lec50_wk4 <- log(45.7)
    label("Log of the pH > 4 holding-time ratio giving half-maximal effect at week 4 (log %)") # Table S2 E_HTR,50 Week 4 = 45.7% (RSE 22.0%)
    e_trt_approved_ec50 <- -0.465
    label("Fractional change in EC50 for a NON-approved regimen, applied when TRT_APPROVED = 0 (fraction)") # Table S2 E_HTR,50 'Non-approval effect' = -46.5% (RSE 4.5%)
    e_trt_pcab_ec50 <- -0.308
    label("Fractional change in EC50 for a potassium-competitive acid blocker (fraction)") # Table S2 E_HTR,50 'P-CAB effect' = -30.8% (RSE 45.1%)
    e_trt_h2ra_ec50 <- 0.413
    label("Fractional change in EC50 for an H2-receptor antagonist (fraction)") # Table S2 E_HTR,50 'H2RA effect' = 41.3% (RSE 21.8%)
    lhill <- log(1.31)
    label("Log of the Emax-model shape (Hill) coefficient (unitless)") # Table S2 Gamma = 1.31 (RSE 33.5%)
    bdv_pct <- 8.1
    label("Between-drug variability: healing-rate increase above 50% for a +1 SD drug effect at SE_ARM = 1 (percentage points)") # Table S2 BDV '% at HR=50%, SE=1' = 8.1 (RSE 27.8%)
    btv_pct <- 6.86
    label("Between-treatment variability: healing-rate increase above 50% for a +1 SD treatment effect at SE_ARM = 1 (percentage points)") # Table S2 BTV '% at HR=50%, SE=1' = 6.86 (RSE 25.8%)
    e_trt_approved_btv <- 2.05
    label("Fractional change in the between-treatment SD for an approved regimen, applied when TRT_APPROVED = 1 (fraction)") # Table S2 BTV 'Non-approval effect' = 205% (RSE 48.7%); applied on APPR = 1 per the Final Model Code
    e_trt_h2ra_btv <- -0.71
    label("Fractional change in the between-treatment SD for an H2-receptor antagonist (fraction)") # Table S2 BTV 'H2RA effect' = -71% (RSE 13.0%)

    # $OMEGA (1 FIX): the drug-level random effect is a standard normal
    # scaled in model() by the BDV logit SD and the arm's SE.
    eta_drug_mu ~ fixed(1)

    # $SIGMA (1 FIX): the treatment-level (residual) random effect is a
    # standard normal scaled in model() by the BTV logit SD and the arm's SE.
    addSd <- fixed(1)
    label("Unit residual SD; scaled per record by the between-treatment logit SD and SE_ARM (unitless)") # $SIGMA (1 FIX)
  })

  model({
    # Week-4 vs week-8 assessment (NONMEM WEEK.EQ.1 = week 4). Observations
    # are defined only at time = 4 or time = 8 weeks.
    wk4 <- 0
    if (time < 6) wk4 <- 1
    ec50 <- exp(lec50_wk8 * (1 - wk4) + lec50_wk4 * wk4) *
      (1 + e_trt_approved_ec50 * (1 - TRT_APPROVED)) *
      (1 + e_trt_pcab_ec50 * TRT_PCAB) *
      (1 + e_trt_h2ra_ec50 * TRT_H2RA)
    hill <- exp(lhill)

    # Sigmoid Emax on the logit scale (MU in the Final Model Code)
    mu <- logite0 +
      GASTRIC_PH4_TIME_PCT^hill * (logitemax - logite0) / (ec50^hill + GASTRIC_PH4_TIME_PCT^hill)

    # Random-effect SDs on the logit scale: BDV/BTV are the healing-rate
    # percentage points above 50% for a +1 SD effect at SE = 1, so
    # logit SD = logit((value + 50) / 100).
    sdDrug <- log((bdv_pct + 50) / (50 - bdv_pct))
    sdTrt <- log((btv_pct + 50) / (50 - btv_pct)) *
      (1 + e_trt_approved_btv * TRT_APPROVED) *
      (1 + e_trt_h2ra_btv * TRT_H2RA)

    logit_ee_healing <- mu + SE_ARM * sdDrug * eta_drug_mu
    ee_healing <- 100 * expit(logit_ee_healing)

    sdRes <- addSd * sdTrt * SE_ARM
    logit_ee_healing ~ add(sdRes)
  })
}
