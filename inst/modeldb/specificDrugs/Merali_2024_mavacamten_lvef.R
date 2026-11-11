Merali_2024_mavacamten_lvef <- function() {
  description <- paste(
    "Exposure-response (safety) model for left ventricular ejection fraction",
    "(LVEF) in adults with obstructive or nonobstructive hypertrophic",
    "cardiomyopathy given mavacamten or placebo (Merali 2024). LVEF is the",
    "individual baseline value LVEF0 times (1 - kMAVA * Cavg168^qq), where",
    "Cavg168 (CAV) is the average mavacamten plasma concentration over the",
    "168 h before the echocardiogram, supplied from a population PK model such",
    "as Chang_2024_mavacamten; there is no placebo effect. The observed",
    "baseline LVEF, nonobstructive HCM and female sex enter LVEF0; the",
    "observed baseline LVEF and studies other than EXPLORER-HCM enter kMAVA,",
    "and studies other than EXPLORER-HCM shift the exposure exponent qq.",
    "Correlated log-normal IIV on LVEF0 and kMAVA and a log-normal residual",
    "error. The paired efficacy model is",
    "Merali_2024_mavacamten_lvot_gradient_valsalva."
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
    lvef = "%"
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
        "modellib('Chang_2024_mavacamten')). Enters uncentred as",
        "kMAVA * CAV^qq; 0 for placebo. Time-varying: it follows dose",
        "titrations and interruptions."
      ),
      source_name = "Cavg168"
    ),
    LVEF_BL = list(
      description = "Observed left ventricular ejection fraction at baseline (before the first study dose).",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power covariate on LVEF0 and on kMAVA, centred at 74% (Merali 2024",
        "Table 2 'log(LVEF0/74)'). Patients rolling from a parent study into",
        "an extension had this value reset to the extension baseline",
        "(Methods 'Data'). Analysis-set mean (SD) 73% (6.4) (Table S2)."
      ),
      source_name = "LVEF0 (baseline LVEF; BL.LVEF in Data S2)"
    ),
    DIS_HCM_NONOBSTRUCTIVE = list(
      description = "Nonobstructive hypertrophic cardiomyopathy indicator, 1 = nonobstructive HCM, 0 = obstructive HCM.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (obstructive HCM)",
      notes = "Merali 2024 Table 2 'LVEF0 ~ nonobstructive HCM'; 59 of the 331 patients (17.8%; Table S2); MAVERICK-HCM is the nonobstructive-HCM study of the pool.",
      source_name = "Population (obstructive HCM vs nonobstructive HCM)"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Merali 2024 Table 2 'LVEF0 ~ female'; 43.5% female (Table S2).",
      source_name = "Sex"
    ),
    STUDY_EXPLORERHCM = list(
      description = "EXPLORER-HCM study indicator, 1 = record from the phase 3 EXPLORER-HCM study (NCT03470545), 0 = any other study of the analysis.",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (EXPLORER-HCM); the paper's reference patient is an EXPLORER-HCM patient",
      notes = paste(
        "Merali 2024 Table 2 'all studies except EXPLORER-HCM' shifts log",
        "kMAVA by +4.29 and qq by -0.64, so the effects are applied to",
        "(1 - STUDY_EXPLORERHCM). Records from PIONEER-HCM, MAVERICK-HCM,",
        "PIONEER-OLE and MAVA-LTE carry 0, including EXPLORER-HCM patients",
        "once they roll over into the EXPLORER-LTE cohort of MAVA-LTE."
      ),
      source_name = "all studies except EXPLORER-HCM"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 331L,
    n_studies = 5L,
    n_observations = 4175L,
    age_range = "Adults (not reported for the analysis set)",
    weight_range = "Not reported for the analysis set",
    sex_female_pct = 43.5,
    disease_state = "Symptomatic obstructive hypertrophic cardiomyopathy (272 patients) or nonobstructive hypertrophic cardiomyopathy (59 patients)",
    dose_range = "Placebo or mavacamten 2.5, 5, 10 or 15 mg orally once daily, titrated on drug concentration and/or echocardiography",
    regions = "Multinational",
    baseline_lvef = "mean 73% (SD 6.4)",
    notes = paste(
      "All five studies: PIONEER-HCM (NCT02842242), MAVERICK-HCM",
      "(NCT03442764), EXPLORER-HCM (NCT03470545), PIONEER-OLE (NCT03496168)",
      "and MAVA-LTE (NCT03723655) (Merali 2024 Methods 'Data'; Results).",
      "Sex not recorded for 1 patient. Beta-blocker use 74.0%, NYHA class II",
      "72.8% / III 27.2% (Table S2). Fitted with nlme 3.1-137 in R 4.0."
    )
  )

  ini({
    # Reference patient: male EXPLORER-HCM patient with obstructive HCM and
    # baseline LVEF 74% (Merali 2024 Table 2 footnote). The model is fitted
    # on the log scale (Results: log-normal random effects on LVEF0 and
    # kMAVA, log-normal RUV); the Table 2 rows labelled '%' are log-scale
    # quantities times 100, so each is divided by 100 here (see the vignette
    # Assumptions and deviations).
    lrbase <- log(74.0); label("Typical baseline LVEF0 for the reference patient (%)") # Table 2 'LVEF0.REF, %' 74.0 (SE 0.243)
    e_lvef_bl_rbase <- 0.336; label("Power exponent of baseline LVEF/74 on LVEF0 (unitless)") # Table 2 'LVEF0 ~ log(LVEF0/74), %' 33.6 / 100
    e_dis_hcm_nonobstructive_rbase <- -0.0324; label("Additive shift in log LVEF0 for nonobstructive vs obstructive HCM (log units)") # Table 2 'LVEF0 ~ nonobstructive HCM, %' -3.24 / 100
    e_sexf_rbase <- 0.0140; label("Additive shift in log LVEF0 for female vs male (log units)") # Table 2 'LVEF0 ~ female, %' 1.40 / 100
    lslope <- -12.8; label("Log typical exposure-response slope kMAVA on Cavg168^qq for EXPLORER-HCM (log units)") # Table 2 'Log(kMAVA.REF), 1/(ng/mL)' -12.8 (SE 0.901)
    e_study_explorerhcm_slope <- 4.29; label("Additive shift in log kMAVA for studies other than EXPLORER-HCM (applied when STUDY_EXPLORERHCM = 0; log units)") # Table 2 'Log(kMAVA.REF) ~ all studies except EXPLORER-HCM' 4.29 (SE 0.869)
    e_lvef_bl_slope <- 1.73; label("Power exponent of baseline LVEF/74 on kMAVA (unitless)") # Table 2 'Log(kMAVA.REF) ~ log(LVEF0/74)' 1.73 (SE 0.418)
    pow <- 1.69; label("Exponent qq on Cavg168 for EXPLORER-HCM (unitless)") # Table 2 'qqREF' 1.69 (SE 0.146)
    e_study_explorerhcm_pow <- -0.64; label("Additive shift in qq for studies other than EXPLORER-HCM (applied when STUDY_EXPLORERHCM = 0; unitless)") # Table 2 'qqREF ~ all studies except EXPLORER-HCM' -0.64 (SE 0.142)

    # LVEF0 IIV is reported as an SD ('SD (LVEF0.REF), %' 2.38 -> 0.0238 on
    # the log scale); kMAVA IIV as %CV, omega^2 = log(1 + CV^2).
    etalrbase + etalslope ~ c(
      0.0238^2,
      -0.062 * 0.0238 * sqrt(log(1 + 0.616^2)), log(1 + 0.616^2)
    ) # Table 2 'SD (LVEF0.REF), %' 2.38, '%CV (kMAVA.REF)' 61.6, 'Corr (LVEF0, kMAVA)' -0.062

    expSd_lvef <- 0.048; label("Log-normal residual SD on LVEF (unitless)") # Table 2 'Residual SD' 4.80 / 100
  })

  model({
    # Covariates are the OBSERVED baseline values of the patient; the study
    # effects act on records outside EXPLORER-HCM.
    non_explorer <- 1 - STUDY_EXPLORERHCM
    rbase <- exp(
      lrbase + etalrbase +
        e_lvef_bl_rbase * log(LVEF_BL / 74) +
        e_dis_hcm_nonobstructive_rbase * DIS_HCM_NONOBSTRUCTIVE +
        e_sexf_rbase * SEXF
    )
    slope <- exp(
      lslope + etalslope +
        e_study_explorerhcm_slope * non_explorer +
        e_lvef_bl_slope * log(LVEF_BL / 74)
    )
    qq <- pow + e_study_explorerhcm_pow * non_explorer

    lvef <- rbase * (1 - slope * CAV^qq)

    lvef ~ lnorm(expSd_lvef)
  })
}
