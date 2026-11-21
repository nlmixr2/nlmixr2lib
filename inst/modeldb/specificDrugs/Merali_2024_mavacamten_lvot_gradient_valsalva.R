Merali_2024_mavacamten_lvot_gradient_valsalva <- function() {
  description <- paste(
    "Exposure-response (efficacy) model for the Valsalva left ventricular",
    "outflow tract (LVOT) peak gradient in adults with obstructive",
    "hypertrophic cardiomyopathy given mavacamten or placebo (Merali 2024).",
    "The gradient is the individual baseline value VLVOTg0, times an",
    "exponential placebo time course that declines to a fixed 87.67% of",
    "baseline, times exp(-kMAVA * Cavg168), where Cavg168 (CAV) is the average",
    "mavacamten plasma concentration over the 168 h before the echocardiogram,",
    "supplied from a population PK model such as Chang_2024_mavacamten. The",
    "observed baseline gradient enters VLVOTg0 and kMAVA, and baseline",
    "NT-proBNP enters VLVOTg0, as power terms; correlated log-normal IIV on",
    "VLVOTg0 and kMAVA and an exponential residual error. The paired safety",
    "model is Merali_2024_mavacamten_lvef."
  )
  reference <- paste(
    "Merali S, Salinger DH, Palmisano M, Sehnert AJ, Thanneer N, Back H,",
    "Seroogy JD, Gretler DD, Roy A, Perera V. Recommendation of mavacamten",
    "posology by model-based analyses in adults with obstructive",
    "hypertrophic cardiomyopathy. CPT Pharmacometrics Syst Pharmacol.",
    "2024;13(9):1448-1461. doi:10.1002/psp4.13138. PMID: 38695527.",
    "PMCID: PMC11533099."
  )
  vignette <- "Merali_2024_mavacamten"
  units <- list(
    time = "week",
    dosing = "mg (mavacamten, given to the companion PK model that supplies CAV)",
    concentration = "ng/mL (CAV, the Cavg168 exposure covariate)",
    lvot_gradient_valsalva = "mmHg"
  )

  covariateData <- list(
    CAV = list(
      description = "Average mavacamten plasma concentration over the 168 h (1 week) before the echocardiographic assessment (Cavg168).",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Merali 2024 Methods 'Model development': a backward-rolling 1-week",
        "average, updated weekly, from the individual predictions of the",
        "sponsor population PK model (the same final model as",
        "modellib('Chang_2024_mavacamten')). Enters uncentred inside",
        "exp(-kMAVA * CAV); 0 for placebo. Time-varying: it follows dose",
        "titrations and interruptions."
      ),
      source_name = "Cavg168"
    ),
    LVOT_GRADIENT_VALSALVA_BL = list(
      description = "Observed Valsalva LVOT peak gradient at baseline (before the first study dose).",
      units = "mmHg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power covariate on VLVOTg0 and on kMAVA, centred at 69 mmHg",
        "(Merali 2024 Table 2 'log(VLVOTg0/69)'; Methods 'baseline VLVOTg on",
        "both VLVOTg0 and kMAVA'). Patients rolling from a parent study into",
        "an extension had this value reset to the extension baseline",
        "(Methods 'Data'). Analysis-set mean (SD) 75.1 (34.0) mmHg",
        "(Table S2)."
      ),
      source_name = "VLVOTg0 (baseline VLVOTg; VLVOT_BL in Data S2)"
    ),
    NTPROBNP = list(
      description = "Baseline serum NT-proBNP.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power covariate on VLVOTg0 centred at 0.736 ng/mL (Table 2",
        "'log(NT-proBNP/736)'; the Table 2 footnote prints the reference as",
        "'736 ng/mL' but Table S2 reports NT-proBNP in ng/L (mean 1280), so",
        "the reference is 736 ng/L = 0.736 ng/mL). Supply values in ng/mL."
      ),
      source_name = "NT-proBNP"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 272L,
    n_studies = 4L,
    n_observations = 3597L,
    age_range = "Adults (not reported for the analysis set)",
    weight_range = "Not reported for the analysis set",
    sex_female_pct = 40.8,
    disease_state = "Symptomatic obstructive hypertrophic cardiomyopathy",
    dose_range = "Placebo or mavacamten 2.5, 5, 10 or 15 mg orally once daily, titrated on drug concentration and/or echocardiography",
    regions = "Multinational",
    baseline_vlvotg = "mean 75.1 mmHg (SD 34.0)",
    baseline_lvef = "mean 74% (SD 5.9)",
    notes = paste(
      "Patients with obstructive HCM in PIONEER-HCM (NCT02842242),",
      "EXPLORER-HCM (NCT03470545), PIONEER-OLE (NCT03496168) and the",
      "EXPLORER-LTE cohort of MAVA-LTE (NCT03723655) (Merali 2024 Methods",
      "'Data'; Results). Beta-blocker use 75.4%, non-dihydropyridine calcium",
      "channel blocker use 25.7%, NYHA class II 71.7% / III 28.3%, pacemaker",
      "9.9% (Table S2). Fitted with nlme 3.1-137 in R 4.0."
    )
  )

  ini({
    # Reference patient: baseline VLVOTg 69 mmHg and baseline NT-proBNP
    # 736 ng/L (Merali 2024 Table 2 footnote). The paper reports the
    # log-scale estimates of VLVOTg0 and kMAVA directly.
    lrbase <- 4.22; label("Log typical baseline Valsalva LVOT gradient VLVOTg0 for the reference patient (log mmHg)") # Table 2 'Log(VLVOTg0.REF), log(mm Hg)' 4.22 (SE 0.0248)
    e_lvot_gradient_valsalva_bl_rbase <- 0.521; label("Power exponent of baseline VLVOTg/69 on VLVOTg0 (unitless)") # Table 2 'Log(VLVOTg0) ~ log(VLVOTg0/69)' 0.521 (SE 0.0319)
    e_ntprobnp_rbase <- 0.108; label("Power exponent of baseline NT-proBNP/736 ng/L on VLVOTg0 (unitless)") # Table 2 'Log(VLVOTg0) ~ log(NT-proBNP/736)' 0.108 (SE 0.0152)
    lslope <- -5.66; label("Log typical exposure-response rate constant kMAVA on Cavg168 (log 1/(ng/mL))") # Table 2 'Log(kMAVA.REF), log(1/(ng/mL))' -5.66 (SE 0.0343)
    e_lvot_gradient_valsalva_bl_slope <- 0.438; label("Power exponent of baseline VLVOTg/69 on kMAVA (unitless)") # Table 2 'Log(kMAVA) ~ log(VLVOTg0/69)' 0.438 (SE 0.0506)
    lkpbo <- log(0.0373); label("Placebo-effect decay rate constant kPBO (1/week)") # Table 2 'kPBO, 1/wk' 0.0373 (SE 0.0205)
    pbomax <- fixed(0.8767); label("Fraction of baseline VLVOTg remaining at the maximal placebo effect VLVOTgPBOMAX (unitless)") # Table 2 'VLVOTgPBOMAX' 0.8767; footnote: estimated on placebo data only (SE 0.0316) and fixed in the final model

    # Table 2 reports %CV for the log-normal IIVs; omega^2 = log(1 + CV^2).
    etalrbase + etalslope ~ c(
      log(1 + 0.245^2),
      0.010 * sqrt(log(1 + 0.245^2) * log(1 + 0.433^2)), log(1 + 0.433^2)
    ) # Table 2 '%CV (VLVOTg0.REF)' 24.5, '%CV (kMAVA.REF)' 43.3, 'Corr (VLVOTg0, kMAVA)' 0.010

    expSd_lvot_gradient_valsalva <- 0.468; label("Exponential (log-scale) residual SD on the Valsalva LVOT gradient (unitless)") # Table 2 'Residual SD (exponential error)' 0.468
  })

  model({
    # Covariates are the OBSERVED baseline values of the patient
    rbase <- exp(
      lrbase + etalrbase +
        e_lvot_gradient_valsalva_bl_rbase * log(LVOT_GRADIENT_VALSALVA_BL / 69) +
        e_ntprobnp_rbase * log(NTPROBNP / 0.736)
    )
    slope <- exp(
      lslope + etalslope +
        e_lvot_gradient_valsalva_bl_slope * log(LVOT_GRADIENT_VALSALVA_BL / 69)
    )
    kpbo <- exp(lkpbo)

    # Placebo (background-treatment) time course; t in weeks from baseline
    pbo <- pbomax + (1 - pbomax) * exp(-kpbo * t)

    lvot_gradient_valsalva <- rbase * pbo * exp(-slope * CAV)

    lvot_gradient_valsalva ~ lnorm(expSd_lvot_gradient_valsalva)
  })
}
