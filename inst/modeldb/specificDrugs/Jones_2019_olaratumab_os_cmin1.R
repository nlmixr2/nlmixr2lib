Jones_2019_olaratumab_os_cmin1 <- function() {
  description <- "Parametric time-to-event exposure-response model for OVERALL SURVIVAL in adults with locally advanced or metastatic soft tissue sarcoma treated with doxorubicin with or without olaratumab 15 mg/kg IV on days 1 and 8 of 21-day cycles (Jones 2019; randomised phase 2 trial JGDG, NCT01185964, n = 133). Constant (exponential) baseline hazard of 0.00205 deaths per day, multiplied by (1 + 0.862 * ECOG_GE1) * (1 - 0.583 * LINE_1L) and by an inhibitory sigmoid Imax function of the olaratumab trough serum concentration at the end of cycle 1, CTROUGH (Imax 0.765, EC50 66.1 ug/mL, Hill 8 fixed). No PK layer: CTROUGH is a time-fixed data column (0 for the doxorubicin-alone arm), the individual prediction from the population PK model Mo_2018_olaratumab. Outputs `hazard` (deaths per day), `cumhaz`, `sur` (survival probability), `hr_drug` (hazard ratio for olaratumab vs none) and `medtte` (closed-form median survival, days). Sister models: Jones_2019_olaratumab_os_cavg, Jones_2019_olaratumab_pfs_cmin1, Jones_2019_olaratumab_pfs_cavg."
  reference <- paste(
    "Jones RL, Mo G, Baldwin JR, Peterson PM, Ilaria RL, Conti I,",
    "Cronier DM, Tap WD.",
    "Exposure-response relationship of olaratumab for survival outcomes and",
    "safety when combined with doxorubicin in patients with soft tissue",
    "sarcoma. Cancer Chemother Pharmacol. 2019;83(1):191.",
    "doi:10.1007/s00280-018-3723-4. PMCID PMC6373189.",
    "Estimates from Table 2 (OS, Cmin1 model).",
    sep = " "
  )
  vignette <- "Jones_2019_olaratumab_exposure_response"
  units <- list(
    time = "day",
    dosing = "n/a (no drug-dosing events; the exposure driver is supplied as the CTROUGH data covariate, in ug/mL)",
    concentration = "probability (the model output `sur` is the overall survival probability, not a drug concentration)"
  )

  compartmentData <- list(
    cumhaz = list(
      analyte = "cumulative hazard of death",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    CTROUGH = list(
      description = "Individual predicted olaratumab trough serum concentration at the end of the first 21-day treatment cycle (Cmin1)",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Required input, time-fixed per subject; set to 0 for patients not receiving olaratumab (doxorubicin-alone arm), which gives hr_drug = 1. Total serum olaratumab at the end of cycle 1 (day 21, before the cycle-2 dose) of 15 mg/kg IV on days 1 and 8, i.e. the intended exposure before any dose reduction. Individual estimate from the population PK model of Mo 2018 (modellib('Mo_2018_olaratumab')). Enters through the inhibitory sigmoid Imax function 1 - imax * CTROUGH^hill / (ec50^hill + CTROUGH^hill), uncentred. Jones 2019 Table 4: range 12.3-188.1 ug/mL in 62 patients, quartile boundaries 62.8, 86.9 and 105.6 ug/mL.",
      source_name = "Cmin1"
    ),
    ECOG_GE1 = list(
      description = "Baseline ECOG performance status >= 1 indicator (1 = ECOG >= 1, 0 = ECOG 0)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Multiplies the baseline hazard by (1 + e_ecog_ge1_haz * ECOG_GE1); Jones 2019 Results: 'Patients with an ECOG PS >= 1 have a predicted 86.2% increase in the baseline hazard'.",
      source_name = "EGRP (ECOG group)"
    ),
    LINE_1L = list(
      description = "No prior line of systemic treatment for advanced disease (1 = no prior line, 0 = one or more prior lines)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Multiplies the baseline hazard by (1 + e_line1l_haz * LINE_1L); Jones 2019 Results: 'patients who had received no prior lines of treatment had a 58.3% decrease in the baseline hazard'. Randomisation was stratified by 0 versus >= 1 prior lines.",
      source_name = "PRVTRT (number of prior treatments)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 133L,
    n_studies = 1L,
    disease_state = "Histologically confirmed locally advanced or metastatic soft tissue sarcoma (Jones 2019 Methods; randomised open-label phase 2 trial JGDG, NCT01185964).",
    dose_range = "Olaratumab 15 mg/kg IV on days 1 and 8 plus doxorubicin 75 mg/m2 on day 1 of each 21-day cycle (n = 66), or doxorubicin 75 mg/m2 alone (n = 67); doxorubicin for up to 8 cycles, olaratumab continued until progression or toxicity.",
    age_range = "not reported in Jones 2019",
    weight_range = "not reported in Jones 2019",
    sex_female_pct = NA_real_,
    regions = "not reported in Jones 2019 (multicenter)",
    notes = "Survival data for all 133 randomised patients; individual olaratumab exposures available for 62 patients of the experimental arm (Table 4). Observed medians in this trial, quoted in the Jones 2019 Introduction from Tap 2016: OS 26.5 vs 14.7 months (HR 0.46) and PFS 6.6 vs 4.1 months (HR 0.67) for olaratumab plus doxorubicin vs doxorubicin alone. Covariates screened (Table 1): ECOG group, tumour size, histology, age group, sex, race, body weight, prior treatment, haemoglobin, albumin; only ECOG group and prior treatment were retained, both on the baseline hazard."
  )

  ini({
    # Jones 2019 Table 2, Overall survival, 'Cmin1 model' column
    # (population estimate, %SEE). Time unit is days: the VPC x-axes in
    # Figures 1a and 2a are 'Time (Days)', and ln(2) / 0.00604 = 115 days
    # = 3.8 months reproduces the 4.1-month doxorubicin-arm median PFS.
    lhaz_base <- log(0.00205); label("Baseline hazard of death for ECOG 0 with >= 1 prior line and no olaratumab (1/day)") # Table 2 OS BaseHazard, Cmin1 model = 0.00205 (2.18 %SEE)
    limax <- log(0.765); label("Maximum fractional decrease of the hazard by olaratumab (unitless)") # Table 2 OS EMAX, Cmin1 model = 0.765 (8.63 %SEE)
    lec50 <- log(66.1); label("Cmin1 giving half the maximal hazard decrease, ECmin150 (ug/mL)") # Table 2 OS ECmin150 = 66.1 ug/mL (12.1 %SEE)
    lhill <- fixed(log(8)); label("Hill coefficient of the olaratumab effect (unitless)") # Table 2 OS Hill = 8 (fixed)
    e_ecog_ge1_haz <- 0.862; label("Fractional change in baseline hazard for ECOG >= 1 (unitless)") # Table 2 OS EGRPBase, Cmin1 model = 0.862 (42.1 %SEE)
    e_line1l_haz <- -0.583; label("Fractional change in baseline hazard for no prior line of treatment (unitless)") # Table 2 OS PRVTRTBase, Cmin1 model = -0.583 (15.9 %SEE)
  })
  model({
    haz_base <- exp(lhaz_base)
    imax <- exp(limax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # Inhibitory sigmoid Imax effect of olaratumab as a fractional
    # decrease of the hazard (Jones 2019 Methods and Figure 1b).
    hr_drug <- 1 - imax * CTROUGH^hill / (ec50^hill + CTROUGH^hill)

    hazard <- haz_base * (1 + e_ecog_ge1_haz * ECOG_GE1) * (1 + e_line1l_haz * LINE_1L) * hr_drug

    d/dt(cumhaz) <- hazard
    cumhaz(0) <- 0
    sur <- exp(-cumhaz)

    # Constant hazard: the median survival time solves exp(-hazard * t) = 0.5.
    medtte <- log(2) / hazard
  })
}
