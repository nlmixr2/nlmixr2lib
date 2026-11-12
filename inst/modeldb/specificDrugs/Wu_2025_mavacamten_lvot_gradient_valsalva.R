Wu_2025_mavacamten_lvot_gradient_valsalva <- function() {
  description <- paste(
    "Exposure-response (efficacy) model for the Valsalva left ventricular",
    "outflow tract (LVOT) peak gradient in adults with obstructive",
    "hypertrophic cardiomyopathy, refitted with Chinese EXPLORER-CN patients",
    "added to the Merali 2024 analysis set (Wu 2025). The gradient is the",
    "individual baseline value VLVOTg0, times an exponential placebo time",
    "course that declines to a fixed 85.1% of baseline, times",
    "exp(-kMAVA * Cavg168), where Cavg168 (CAV) is the average mavacamten",
    "plasma concentration over the 168 h before the echocardiogram, supplied",
    "from a population PK model such as Chang_2024_mavacamten. The observed",
    "baseline gradient enters VLVOTg0 and kMAVA, and baseline NT-proBNP enters",
    "VLVOTg0, as power terms; enrolment in China raises VLVOTg0 and lowers",
    "kMAVA. Correlated log-normal IIV on VLVOTg0 and kMAVA and an exponential",
    "residual error. The paired safety model is Wu_2025_mavacamten_lvef."
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
    lvot_gradient_valsalva = "mmHg"
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
        "inside exp(-kMAVA * CAV); 0 for placebo. Time-varying: it follows",
        "dose titrations and interruptions."
      ),
      source_name = "Cavg168"
    ),
    LVOT_GRADIENT_VALSALVA_BL = list(
      description = "Observed Valsalva LVOT peak gradient at baseline (before the first study dose).",
      units = "mmHg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power covariate on VLVOTg0 and on kMAVA, centred at 74.1 mmHg (Wu",
        "2025 Table 2 'log(VLVOTg0/74.1)' and footnote 'baseline VLVOTg of",
        "74.1 mmHg')."
      ),
      source_name = "VLVOTg0 (baseline VLVOTg)"
    ),
    NTPROBNP = list(
      description = "Baseline serum NT-proBNP.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power covariate on VLVOTg0 centred at 0.674 ng/mL (Table 2",
        "'log(NT-proBNP/674)'; footnote 'baseline NT-proBNP of 674 pg/mL').",
        "Supply values in ng/mL."
      ),
      source_name = "NT-proBNP"
    ),
    REGION_CHINA = list(
      description = "China enrolment indicator, 1 = patient enrolled in a Chinese study (EXPLORER-CN), 0 = enrolled in the global studies.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (global studies of Merali 2024)",
      notes = paste(
        "Wu 2025 Table 2 'Log(VLVOTg0) ~ log(Chinese)' 0.297 and 'Chinese on",
        "kMAVA' -0.351, both additive on the log-transformed parameter",
        "(footnote). Methods call the covariate 'Chinese study population';",
        "all 81 Chinese patients came from EXPLORER-CN, enrolled in China."
      ),
      source_name = "Chinese"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 353L,
    n_studies = 5L,
    age_range = "Adults (not reported for the analysis set)",
    weight_range = "Not reported for the analysis set",
    disease_state = "Symptomatic obstructive hypertrophic cardiomyopathy",
    dose_range = "Placebo or mavacamten 2.5, 5, 10 or 15 mg orally once daily, titrated on drug concentration and/or echocardiography",
    regions = "Multinational plus China",
    notes = paste(
      "The Merali 2024 VLVOTg analysis set (272 patients with obstructive",
      "HCM from PIONEER-HCM, EXPLORER-HCM, PIONEER-OLE and MAVA-LTE) plus 81",
      "Chinese patients with obstructive HCM from EXPLORER-CN (NCT05174416),",
      "who added 514 VLVOTg observations (Wu 2025 Results). n_subjects is",
      "that sum; Wu 2025 does not tabulate the pooled analysis set. 13% of",
      "the EXPLORER-CN mavacamten group were CYP2C19 poor metabolizers",
      "(Methods). Fitted with nonlinear mixed-effects methods as in",
      "Merali 2024."
    )
  )

  ini({
    # Reference patient: baseline VLVOTg 74.1 mmHg and baseline NT-proBNP
    # 674 pg/mL, outside China (Wu 2025 Table 2 footnote). The two
    # 'Log(...REF)' rows print the back-transformed value (77.8 mmHg cannot
    # be a log), so their logs are taken here.
    lrbase <- log(77.8); label("Log typical baseline Valsalva LVOT gradient VLVOTg0 for the reference patient (log mmHg)") # Table 2 'Log(VLVOTg0.REF), log(mmHg)' 77.8 (RSE 2.3%)
    e_lvot_gradient_valsalva_bl_rbase <- 0.365; label("Power exponent of baseline VLVOTg/74.1 on VLVOTg0 (unitless)") # Table 2 'Log(VLVOTg0) ~ log(VLVOTg0/74.1)' 0.365 (RSE 15.9%)
    e_ntprobnp_rbase <- 0.12; label("Power exponent of baseline NT-proBNP/674 pg/mL on VLVOTg0 (unitless)") # Table 2 'Log(VLVOTg0) ~ log(NT-proBNP/674)' 0.12 (RSE 15.8%)
    e_region_china_rbase <- 0.297; label("Additive shift in log VLVOTg0 for patients enrolled in China (log units)") # Table 2 'Log(VLVOTg0) ~ log(Chinese)' 0.297 (RSE 27.9%)
    lslope <- log(0.00291); label("Log typical exposure-response rate constant kMAVA on Cavg168 (log 1/(ng/mL))") # Table 2 'Log(kMAVA.REF), log(1/(ng/mL))' 0.00291 (RSE 4.6%)
    e_lvot_gradient_valsalva_bl_slope <- 0.31; label("Power exponent of baseline VLVOTg/74.1 on kMAVA (unitless)") # Table 2 'Log(kMAVA) ~ log(VLVOTg0/74.1)' 0.31 (RSE 34.4%)
    e_region_china_slope <- -0.351; label("Additive shift in log kMAVA for patients enrolled in China (log units)") # Table 2 'Chinese on kMAVA' -0.351 (RSE 26.4%)
    lkpbo <- log(0.178); label("Placebo-effect decay rate constant kPBO (1/week)") # Table 2 'kPBO, 1/wk' 0.178 (RSE NA)
    pbomax <- fixed(0.851); label("Fraction of baseline VLVOTg remaining at the maximal placebo effect VLVOTgPBOMAX (unitless)") # Table 2 'VLVOTgPBOMAX' 0.851; footnote: estimated with placebo data only and fixed in the final model

    # Table 2 reports %CV for the log-normal IIVs; omega^2 = log(1 + CV^2).
    etalrbase + etalslope ~ c(
      log(1 + 0.37^2),
      0.449 * sqrt(log(1 + 0.37^2) * log(1 + 0.658^2)), log(1 + 0.658^2)
    ) # Table 2 '%CV (VLVOTg0.REF)' 37, '%CV (kMAVA.REF)' 65.8, 'Corr (VLVOTg0, kMAVA)' 0.449

    expSd_lvot_gradient_valsalva <- 0.434; label("Exponential (log-scale) residual SD on the Valsalva LVOT gradient (unitless)") # Table 2 'Residual SD (exponential error)' 0.434
  })

  model({
    # Covariates are the OBSERVED baseline values of the patient
    rbase <- exp(
      lrbase + etalrbase +
        e_lvot_gradient_valsalva_bl_rbase * log(LVOT_GRADIENT_VALSALVA_BL / 74.1) +
        e_ntprobnp_rbase * log(NTPROBNP / 0.674) +
        e_region_china_rbase * REGION_CHINA
    )
    slope <- exp(
      lslope + etalslope +
        e_lvot_gradient_valsalva_bl_slope * log(LVOT_GRADIENT_VALSALVA_BL / 74.1) +
        e_region_china_slope * REGION_CHINA
    )
    kpbo <- exp(lkpbo)

    # Placebo (background-treatment) time course; t in weeks from baseline
    pbo <- pbomax + (1 - pbomax) * exp(-kpbo * t)

    lvot_gradient_valsalva <- rbase * pbo * exp(-slope * CAV)

    lvot_gradient_valsalva ~ lnorm(expSd_lvot_gradient_valsalva)
  })
}
