Wu_2025_mavacamten_lvef <- function() {
  description <- paste(
    "Exposure-response (safety) model for left ventricular ejection fraction",
    "(LVEF) in adults with hypertrophic cardiomyopathy given mavacamten or",
    "placebo, refitted with Chinese EXPLORER-CN patients added to the Merali",
    "2024 analysis set (Wu 2025). LVEF is the individual baseline value LVEF0",
    "times (1 - kMAVA * Cavg168^qq), where Cavg168 (CAV) is the average",
    "mavacamten plasma concentration over the 168 h before the",
    "echocardiogram, supplied from a population PK model such as",
    "Chang_2024_mavacamten; there is no placebo effect. The observed baseline",
    "LVEF, nonobstructive HCM, the septal-reduction-therapy patient type,",
    "female sex and enrolment in China enter LVEF0; the septal-reduction-",
    "therapy patient type and enrolment in China lower kMAVA. Log-normal IIV",
    "on LVEF0 and kMAVA and a log-normal residual error. The paired efficacy",
    "model is Wu_2025_mavacamten_lvot_gradient_valsalva."
  )
  reference <- paste(
    "Wu X, Puli S, Chen N, Tian Z, Hsu P, Sun J, Lyu C, Merali S, Zhang J.",
    "Model-informed recommendation of mavacamten posology for Chinese adults",
    "with obstructive hypertrophic cardiomyopathy. CPT Pharmacometrics Syst",
    "Pharmacol. 2025;14(4):751. doi:10.1002/psp4.13312. PMID: 39871653.",
    "PMCID: PMC12001263."
  )
  vignette <- "Wu_2025_mavacamten"
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
        "Wu 2025 Methods 'Observed data' and Table 2 equation; derived from",
        "the population PK model of Merali 2024 / Chang 2024, which Wu 2025",
        "kept unchanged (modellib('Chang_2024_mavacamten')). Enters uncentred",
        "as kMAVA * CAV^qq; 0 for placebo. Time-varying: it follows dose",
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
        "Power covariate on LVEF0 centred at 72.82% (Wu 2025 Table 2",
        "'LVEF0 ~ log(LVEF0/72.82)'; the footnote's reference LVEF of 74% is",
        "the Merali 2024 value and is not used)."
      ),
      source_name = "LVEF0 (baseline LVEF)"
    ),
    DIS_HCM_NONOBSTRUCTIVE = list(
      description = "Nonobstructive hypertrophic cardiomyopathy indicator, 1 = nonobstructive HCM, 0 = obstructive HCM.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (obstructive HCM)",
      notes = "Wu 2025 Table 2 'LVEF0 ~ nonobstructive HCM'; one level of the patient-type covariate (nHCM, oHCM, SRT) of Table 1.",
      source_name = "nonobstructive HCM (patient type nHCM)"
    ),
    DIS_HCM_SRT = list(
      description = "Septal-reduction-therapy (SRT) patient-type indicator, 1 = patient of the 'SRT' type, 0 = obstructive HCM (oHCM).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (obstructive HCM)",
      notes = paste(
        "Wu 2025 Table 2 'LVEF0 ~ SRT' and 'Log(kMAVA.REF) ~ SRT'; one level",
        "of the patient-type covariate (nHCM, oHCM, SRT) of Table 1. The",
        "paper does not define the SRT patient type; it most plausibly denotes",
        "patients with obstructive HCM referred for septal reduction therapy",
        "(the VALOR-HCM population). Mutually exclusive with",
        "DIS_HCM_NONOBSTRUCTIVE."
      ),
      source_name = "SRT (patient type)"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Wu 2025 Table 2 'LVEF0 ~ female'.",
      source_name = "Sex"
    ),
    REGION_CHINA = list(
      description = "China enrolment indicator, 1 = patient enrolled in a Chinese study (EXPLORER-CN), 0 = enrolled in the global studies.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (global studies of Merali 2024)",
      notes = paste(
        "Wu 2025 Table 2 'LVEF0 ~ Chinese' and 'Log(kMAVA.REF) ~ Chinese'.",
        "Methods call the covariate 'Chinese study population'; all 81",
        "Chinese patients came from EXPLORER-CN, enrolled in China."
      ),
      source_name = "Chinese"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 412L,
    n_studies = 6L,
    age_range = "Adults (not reported for the analysis set)",
    weight_range = "Not reported for the analysis set",
    disease_state = "Symptomatic obstructive or nonobstructive hypertrophic cardiomyopathy",
    dose_range = "Placebo or mavacamten 2.5, 5, 10 or 15 mg orally once daily, titrated on drug concentration and/or echocardiography",
    regions = "Multinational plus China",
    notes = paste(
      "The Merali 2024 LVEF analysis set (331 patients with obstructive or",
      "nonobstructive HCM from PIONEER-HCM, MAVERICK-HCM, EXPLORER-HCM,",
      "PIONEER-OLE and MAVA-LTE) plus 81 Chinese patients with obstructive HCM",
      "from EXPLORER-CN (NCT05174416), who added 519 LVEF observations (Wu",
      "2025 Results). n_subjects is that sum; Wu 2025 does not tabulate the",
      "pooled analysis set and does not name the study that contributed the",
      "SRT patient type."
    )
  )

  ini({
    # Reference patient: male patient with obstructive HCM, baseline LVEF
    # 72.82%, outside China (Wu 2025 Table 2 row labels and footnote). The
    # model is fitted on the log scale (Table 2 note: log-normal IIV and
    # RUV, log-transformed model equations); the rows labelled '%' other
    # than LVEF0.REF are log-scale quantities times 100, so each is divided
    # by 100 here. 'Log(kMAVA.REF)' prints the back-transformed value.
    lrbase <- log(73.1); label("Typical baseline LVEF0 for the reference patient (%)") # Table 2 'LVEF0.REF, %' 73.1 (RSE 0.389%)
    e_lvef_bl_rbase <- 0.315; label("Power exponent of baseline LVEF/72.82 on LVEF0 (unitless)") # Table 2 'LVEF0 ~ log(LVEF0/72.82), %' 0.315 (RSE 10.6%)
    e_dis_hcm_nonobstructive_rbase <- -0.0567; label("Additive shift in log LVEF0 for nonobstructive vs obstructive HCM (log units)") # Table 2 'LVEF0 ~ nonobstructive HCM, %' -5.67 / 100
    e_sexf_rbase <- 0.0212; label("Additive shift in log LVEF0 for female vs male (log units)") # Table 2 'LVEF0 ~ female, %' 2.12 / 100
    e_dis_hcm_srt_rbase <- -0.0617; label("Additive shift in log LVEF0 for the SRT patient type vs obstructive HCM (log units)") # Table 2 'LVEF0 ~ SRT, %' -6.17 / 100
    e_region_china_rbase <- 0.0754; label("Additive shift in log LVEF0 for patients enrolled in China (log units)") # Table 2 'LVEF0 ~ Chinese, %' 7.54 / 100
    lslope <- log(0.0000615); label("Log typical exposure-response slope kMAVA on Cavg168^qq (log units)") # Table 2 'Log(kMAVA.REF), 1/(ng/mL)' 0.0000615 (RSE 63.8%)
    e_dis_hcm_srt_slope <- -0.469; label("Additive shift in log kMAVA for the SRT patient type vs obstructive HCM (log units)") # Table 2 'Log(kMAVA.REF) ~ SRT (%)' -46.9 / 100
    e_region_china_slope <- -0.871; label("Additive shift in log kMAVA for patients enrolled in China (log units)") # Table 2 'Log(kMAVA.REF) ~ Chinese (%)' -87.1 / 100
    pow <- 1.20; label("Exponent qq on Cavg168 (unitless)") # Table 2 'qqREF' 1.20 (RSE 8.37%)

    # LVEF0 IIV is reported as an SD ('SD (LVEF0.REF), %' 3.86 -> 0.0386 on
    # the log scale); kMAVA IIV as %CV, omega^2 = log(1 + CV^2). No
    # correlation is reported for the LVEF model.
    etalrbase ~ 0.0386^2 # Table 2 'SD (LVEF0.REF), %' 3.86 / 100
    etalslope ~ log(1 + 0.762^2) # Table 2 '%CV (kMAVA.REF), 1/(ng/mL)' 76.2

    expSd_lvef <- 0.0461; label("Log-normal residual SD on LVEF (unitless)") # Table 2 'Residual SD' 4.61 / 100
  })

  model({
    # Covariates are the OBSERVED baseline values of the patient
    rbase <- exp(
      lrbase + etalrbase +
        e_lvef_bl_rbase * log(LVEF_BL / 72.82) +
        e_dis_hcm_nonobstructive_rbase * DIS_HCM_NONOBSTRUCTIVE +
        e_sexf_rbase * SEXF +
        e_dis_hcm_srt_rbase * DIS_HCM_SRT +
        e_region_china_rbase * REGION_CHINA
    )
    slope <- exp(
      lslope + etalslope +
        e_dis_hcm_srt_slope * DIS_HCM_SRT +
        e_region_china_slope * REGION_CHINA
    )
    qq <- pow

    lvef <- rbase * (1 - slope * CAV^qq)

    lvef ~ lnorm(expSd_lvef)
  })
}
