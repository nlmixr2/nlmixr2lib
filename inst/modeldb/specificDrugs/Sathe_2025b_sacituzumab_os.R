Sathe_2025b_sacituzumab_os <- function() {
  description <- paste0(
    "Cox proportional-hazards exposure-efficacy model for OVERALL ",
    "SURVIVAL in adults with metastatic triple-negative breast cancer ",
    "treated with sacituzumab govitecan 8 or 10 mg/kg on days 1 and 8 of ",
    "21-day cycles (Sathe 2025b, n = 277, pooled IMMU-132-01 and ",
    "ASCENT). The log relative hazard is -0.0123 * CAV + 0.000622 * LDH, ",
    "where CAV is the individual average serum concentration of TOTAL ",
    "ANTIBODY (tAB) from the first dose to death or censoring, in ug/mL, ",
    "and LDH is baseline serum lactate dehydrogenase in IU/L. Outputs: ",
    "lhr, hr and prob_os_12mo, the probability of being alive at 12 ",
    "months, computed with a 12-month baseline cumulative hazard ",
    "BACK-SOLVED from the two median-LDH predictions of Table S12 at the ",
    "Table S3 median LDH of 277 IU/L; a Cox fit leaves the baseline ",
    "hazard nonparametric and the paper does not print it. There is no ",
    "PK layer and no ODE; exposure derives from the companion population ",
    "PK model Sathe_2024_sacituzumab. One of six exposure-response models ",
    "in the Sathe_2025b_sacituzumab_* family."
  )
  reference <- paste(
    "Sathe AG, Diderichsen PM, Fauchet F, Phan SC, Girish S, Othman AA.",
    "Exposure-Response Analyses of Sacituzumab Govitecan Efficacy and",
    "Safety in Patients With Metastatic Triple-Negative Breast Cancer.",
    "Clin Pharmacol Ther. 2025;117(2):570-578. doi:10.1002/cpt.3495.",
    "PMCID PMC11739744. Coefficients from Supplementary Table S10; the",
    "12-month baseline cumulative hazard is back-solved from Table S12.",
    sep = " "
  )
  vignette <- "Sathe_2025b_sacituzumab_exposure_response"

  units <- list(
    time = "months (prob_os_12mo is a single landmark at 12 months; lhr and hr are time-invariant under proportional hazards)",
    dosing = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "prob_os_12mo (probability of being alive at 12 months, 0-1); also hr (relative hazard, unitless) and lhr (its log)"
  )

  covariateData <- list(
    CAV = list(
      description = "Individual average serum concentration of sacituzumab govitecan TOTAL ANTIBODY (tAB, hRS7 conjugated or unconjugated to SN-38) from the first dose to death or censoring. Supplied as data: this model has no PK layer.",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Sathe 2025b Eq. 1 with t_event from Table S1 (time of event, or 21 days after the last SG dose if no event was reported). CAVG_tAB was the most significant of nine exposure metrics for OS (Table S9, -2LL drop 67.9). Reproduce the column from the total-antibody output (Cc_tab) of modellib('Sathe_2024_sacituzumab'). Enters LINEARLY and UNCENTRED. Table 2 OS exposure quartiles: medians 88.7, 125, 153 and 204 ug/mL (overall range 50.6-451); Table S12 footnote: mean 125 ug/mL at 8 mg/kg and 148 ug/mL at 10 mg/kg.",
      source_name = "CAVG_tAB"
    ),
    LDH = list(
      description = "Baseline serum lactate dehydrogenase",
      units = "IU/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters LINEARLY per IU/L and UNCENTRED (Sathe 2025b Table S10, 'Slope for LDH (per IU/L)'). Table S3: median 277 IU/L (range 107-5930) overall, 241 in IMMU-132-01 and 283 in ASCENT; one ASCENT patient was missing. The abstract quotes 'median lactate dehydrogenase (227 IU/L)'; that is a digit transposition, since the median of the pooled set must lie between the two study medians (241 and 283).",
      source_name = "LDH"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 277L,
    n_studies = 2L,
    age_range = "27-82 years (median 54)",
    weight_range = "37.2-132 kg (median 66.7)",
    sex_female_pct = 99.3,
    race_ethnicity = c(White = 80.5, Black = 11.2, Asian = 4.0, Other = 4.3),
    disease_state = "metastatic triple-negative breast cancer, previously treated",
    dose_range = "sacituzumab govitecan 8, 10 or 12 mg/kg starting dose IV on days 1 and 8 of 21-day cycles (10 mg/kg in about 92%)",
    regions = "North America 67.9%, Europe 32.1%",
    notes = "191 deaths (69.0%; Table S4). The proportional-hazards assumption was formally violated for CAVG_tAB, but Schoenfeld and Martingale residuals showed limited bias and the model was accepted (Supplementary Section 1). Prior treatment, body weight, age, race, region, number of prior lines, ECOG performance status, Trop-2 H score and UGT1A1 genotype were screened and not retained."
  )

  ini({
    # Sathe 2025b Supplementary Table S10, OS rows (log-hazard scale;
    # positive = higher risk of death).
    e_cav_os_haz <- -0.0123; label("Log hazard ratio for death per 1 ug/mL increase in CAVG_tAB (per ug/mL)") # Sathe 2025b Table S10, OS Slope for CAVG_tAB (per ug/mL): -0.0123 [-0.0161 to -0.00854], p < 0.001
    e_ldh_os_haz <- 0.000622; label("Log hazard ratio for death per 1 IU/L increase in baseline LDH (per IU/L)") # Sathe 2025b Table S10, OS Slope for LDH (per IU/L): 0.000622 [0.000411-0.000832], p < 0.001

    # Baseline (all covariates zero) cumulative hazard of death at 12
    # months, NOT printed. Back-solved from the Table S12 median-LDH
    # rows as -log(S12) / exp(lhr) at LDH = 277 IU/L and the footnoted
    # mean CAVG_tAB: 3.315 (125 ug/mL, S = 0.429) and 3.300 (148 ug/mL,
    # S = 0.530); mean 3.31. With this value the lower- and upper-quartile
    # rows imply LDH quartiles of about 200 and 446 IU/L.
    ch0_os_12mo <- 3.31; label("Baseline cumulative hazard of death at 12 months for the all-zero covariate reference (unitless)") # back-solved from Sathe 2025b Table S12 median-LDH rows at LDH 277 IU/L; no printed value exists
  })

  model({
    lhr <- e_cav_os_haz * CAV + e_ldh_os_haz * LDH
    hr <- exp(lhr)

    # Probability of being alive at the 12-month landmark. Valid at 12
    # months only: the baseline hazard at other times is not recoverable
    # from the paper.
    prob_os_12mo <- exp(-ch0_os_12mo * hr)
  })
}
