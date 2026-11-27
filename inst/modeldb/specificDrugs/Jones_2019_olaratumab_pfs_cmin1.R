Jones_2019_olaratumab_pfs_cmin1 <- function() {
  description <- "Parametric time-to-event exposure-response model for PROGRESSION-FREE SURVIVAL in adults with locally advanced or metastatic soft tissue sarcoma treated with doxorubicin with or without olaratumab 15 mg/kg IV on days 1 and 8 of 21-day cycles (Jones 2019; randomised phase 2 trial JGDG, NCT01185964, n = 133). Constant (exponential) baseline hazard of 0.00604 progression-or-death events per day, multiplied by an inhibitory sigmoid Imax function of the olaratumab trough serum concentration at the end of cycle 1, CTROUGH (Imax 0.571, EC50 82.0 ug/mL, Hill 8 fixed); no patient covariate was retained. No PK layer: CTROUGH is a time-fixed data column (0 for the doxorubicin-alone arm), the individual prediction from the population PK model Mo_2018_olaratumab. Outputs `hazard` (events per day), `cumhaz`, `sur` (progression-free survival probability), `hr_drug` (hazard ratio for olaratumab vs none) and `medtte` (closed-form median PFS, days). Sister models: Jones_2019_olaratumab_pfs_cavg, Jones_2019_olaratumab_os_cmin1, Jones_2019_olaratumab_os_cavg."
  reference <- paste(
    "Jones RL, Mo G, Baldwin JR, Peterson PM, Ilaria RL, Conti I,",
    "Cronier DM, Tap WD.",
    "Exposure-response relationship of olaratumab for survival outcomes and",
    "safety when combined with doxorubicin in patients with soft tissue",
    "sarcoma. Cancer Chemother Pharmacol. 2019;83(1):191.",
    "doi:10.1007/s00280-018-3723-4. PMCID PMC6373189.",
    "Estimates from Table 3 (PFS, Cmin1 model).",
    sep = " "
  )
  vignette <- "Jones_2019_olaratumab_exposure_response"
  units <- list(
    time = "day",
    dosing = "n/a (no drug-dosing events; the exposure driver is supplied as the CTROUGH data covariate, in ug/mL)",
    concentration = "probability (the model output `sur` is the progression-free survival probability, not a drug concentration)"
  )

  compartmentData <- list(
    cumhaz = list(
      analyte = "cumulative hazard of progression or death",
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
    notes = "Survival data for all 133 randomised patients; individual olaratumab exposures available for 62 patients of the experimental arm (Table 4). Tumour response was assessed every 6 weeks by RECIST 1.1. Observed medians in this trial, quoted in the Jones 2019 Introduction from Tap 2016: PFS 6.6 vs 4.1 months (HR 0.67) for olaratumab plus doxorubicin vs doxorubicin alone. None of the Table 1 covariates (ECOG group, tumour size, histology, age group, sex, race, body weight, prior treatment, haemoglobin, albumin) was retained for PFS."
  )

  ini({
    # Jones 2019 Table 3, Progression-free survival, 'Cmin1 model' column
    # (population estimate, %SEE). The typeset %SEE values '2-52',
    # '15-4' and '6-15' are read as 2.52, 15.4 and 6.15 (a dash printed
    # in place of the decimal point). Time unit is days (Figure 2a x-axis
    # 'Time (Days)').
    lhaz_base <- log(0.00604); label("Baseline hazard of progression or death without olaratumab (1/day)") # Table 3 PFS BaseHazard, Cmin1 model = 0.00604 (2.52 %SEE)
    limax <- log(0.571); label("Maximum fractional decrease of the hazard by olaratumab (unitless)") # Table 3 PFS EMAX, Cmin1 model = 0.571 (15.4 %SEE)
    lec50 <- log(82.0); label("Cmin1 giving half the maximal hazard decrease, ECmin150 (ug/mL)") # Table 3 PFS ECmin150 = 82.0 ug/mL (6.15 %SEE)
    lhill <- fixed(log(8)); label("Hill coefficient of the olaratumab effect (unitless)") # Table 3 PFS Hill = 8 (fixed)
  })
  model({
    haz_base <- exp(lhaz_base)
    imax <- exp(limax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # Inhibitory sigmoid Imax effect of olaratumab as a fractional
    # decrease of the hazard (Jones 2019 Methods and Figure 2b).
    hr_drug <- 1 - imax * CTROUGH^hill / (ec50^hill + CTROUGH^hill)

    hazard <- haz_base * hr_drug

    d/dt(cumhaz) <- hazard
    cumhaz(0) <- 0
    sur <- exp(-cumhaz)

    # Constant hazard: the median survival time solves exp(-hazard * t) = 0.5.
    medtte <- log(2) / hazard
  })
}
