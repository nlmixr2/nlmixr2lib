Shahraz_2026_savolitinib_orr_lowmet <- function() {
  description <- paste(
    "Exposure-efficacy logistic regression for objective response to",
    "savolitinib plus osimertinib in the SAVANNAH subgroup 'Low biomarker, no prior chemo or IO'",
    "(Shahraz 2026, Figure 3c): the probability of an objective response is",
    "expit(intercept + slope * log(CTROUGH)), with CTROUGH the individual",
    "steady-state savolitinib trough concentration (Cmin,ss) supplied as",
    "data. The coefficients are digitised from the plotted curve because",
    "the paper prints none. No PK layer, no ODE and no observation",
    "endpoint; prob_orr is an algebraic output."
  )
  reference <- paste(
    "Shahraz A, Derippe T, Kode KK, Allmond J, Nelson S, Igwegbe I,",
    "Sawant-Basak A, Vishwanathan K, Miah K. Population Pharmacokinetics and",
    "Exposure-Response Analyses of Savolitinib Alone or in Combination With",
    "Osimertinib in Patients With Locally Advanced or Metastatic Solid",
    "Tumors. Clin Transl Sci. 2026;19:e70700. doi:10.1111/cts.70700.",
    "PMCID PMC13504015.",
    sep = " "
  )
  vignette <- "Shahraz_2026_savolitinib"
  units <- list(
    time = "n/a (static landmark exposure-response model; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CTROUGH covariate column)",
    concentration = "prob_orr (probability of an objective response, 0-1)"
  )

  covariateData <- list(
    CTROUGH = list(
      description = "Individual steady-state savolitinib trough plasma concentration (Cmin,ss), supplied as data (this model has no PK layer).",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "TOTAL savolitinib, steady state, from the individual empirical Bayes estimates of the population PK model (modellib('Shahraz_2026_savolitinib')) at each patient's most prevalent dose and regimen (Shahraz 2026 Methods); the minimum over the dosing interval of an oral drug, so Cmin,ss and the pre-dose trough coincide. Enters the logit on the natural-log scale, uncentred: values must be > 0 (log(0) = -Inf gives a probability of 0 for a positive slope). Plotted range of Figure 3: about 0.7-2500 ng/mL.",
      source_name = "Cmin,ss"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 87L,
    n_studies = 1L,
    n_responders = 12L,
    disease_state = "EGFR-mutated, MET-overexpressed and/or amplified locally advanced or metastatic non-small cell lung cancer after progression on osimertinib; SAVANNAH subgroup: low MET biomarker status (MET IHC3+/>= 50% and/or FISH5+ but not IHC3+/>= 90% and/or FISH10+), all lines of therapy, excluding patients with prior chemotherapy or immunotherapy",
    dose_range = "Savolitinib 300 mg once daily, 600 mg once daily or 300 mg twice daily (most prevalent dose) with osimertinib 80 mg once daily",
    regions = "Multinational (SAVANNAH, NCT03778229)",
    notes = "Shahraz 2026 Figure 3c: 87 patients (12 responders); Cmin,ss quartiles Q1-Q4 had 3/22, 2/22, 2/21 and 5/22 responders (13.6%, 9.1%, 9.5%, 22.7%). The 24 patients on savolitinib plus placebo were excluded from the exposure-efficacy analysis."
  )

  ini({
    logit_ref_orr <- -2.161
    label("Logit of the objective-response probability at Cmin,ss = 1 ng/mL (unitless logit)") # Digitised by the maintainers from the Figure 3c curve (logistic in log Cmin,ss; RMS about 0.002 on the probability scale)
    e_logctrough_orr <- 0.1150
    label("Change in the objective-response logit per unit natural log of Cmin,ss in ng/mL (unitless logit)") # Digitised by the maintainers from the Figure 3c curve
  })

  model({
    prob_orr <- expit(logit_ref_orr + e_logctrough_orr * log(CTROUGH))
  })
}
