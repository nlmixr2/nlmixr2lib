Ojara_2023_paclitaxel_os <- function() {
  description <- paste(
    "Weibull parametric time-to-event model of overall survival (OS) in",
    "paclitaxel/platinum-treated advanced non-small cell lung cancer",
    "(CEPAC-TDM phase III study, 221 patients with all seven baseline blood",
    "biomarkers). Hazard h(t) = lambda * alpha * (lambda * t)^(alpha - 1) *",
    "exp(E_RS8 * (RS8 - 80.0)) * exp(E_CRP * (CRP - 23.5)), where RS8 is",
    "the tumor size at week 8 relative to baseline (percent) and CRP the",
    "baseline serum C-reactive protein (mg/L). RS8 enters through the",
    "RCFB8_SLD covariate (fractional change, RS8 = 100 * (1 + RCFB8_SLD));",
    "in the paper it is derived from the companion paclitaxel-exposure",
    "driven tumor growth inhibition model (Ojara_2023_paclitaxel_tgi). No",
    "interindividual variability. Time 0 is treatment start."
  )
  reference <- paste(
    "Ojara FW, Henrich A, Frances N, Nassar YM, Huisinga W, Hartung N,",
    "Geiger K, Holdenrieder S, Joerger M, Kloft C. A prognostic baseline",
    "blood biomarker and tumor growth kinetics integrated model in",
    "paclitaxel/platinum treated advanced non-small cell lung cancer",
    "patients. CPT Pharmacometrics Syst Pharmacol. 2023;12(11):1714-1725.",
    "doi:10.1002/psp4.12937"
  )
  vignette <- "Ojara_2023_paclitaxel_nsclc"
  units <- list(
    time = "h",
    dosing = "none (no drug input; prognostic survival model)",
    concentration = "not applicable"
  )

  covariateData <- list(
    CRP = list(
      description = "Baseline (cycle 1 day 1) serum C-reactive protein, standard assay (latex-enhanced turbidimetry on a Roche Cobas C501 analyzer).",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Per-subject, baseline only. Centered at the median 23.5 mg/L (Table S2; NONMEM code in Appendix S1) inside exp(e_crp_haz * (CRP - 23.5)). 5th/95th percentiles 1.18 / 151 mg/L (Table S5, Figure 5). The paper notes that CRP values are analyzer-dependent and the model should be re-evaluated if a different method is used.",
      source_name = "BASECRP"
    ),
    RCFB8_SLD = list(
      description = "Tumor size (RECIST 1.1 sum of target-lesion diameters) at week 8 relative to baseline, expressed as the fractional change from baseline (TS(week 8) - TS(0)) / TS(0). The paper's RS8 is the same quantity as a percent of baseline: RS8 = 100 * (1 + RCFB8_SLD).",
      units = "unitless (fraction of baseline)",
      type = "continuous",
      reference_category = NULL,
      notes = "Per-subject, time-fixed. Model-derived in the paper: the individual week-8 tumor size from the paclitaxel AUCcycle-driven TGI model (Ojara_2023_paclitaxel_tgi), also for patients who died before week 8. Median RS8 80.0% (RCFB8_SLD = -0.200), the centering value of the hazard term; 5th/95th percentiles 42.4% / 91.4% (Table S5).",
      source_name = "RS8"
    )
  )

  compartmentData <- list(
    cumhaz = list(
      analyte = "cumulative hazard of death",
      units = "(unitless)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 221L,
    n_studies = 1L,
    age_range = "41-78 years (median 64.0)",
    weight_range = "BSA 1.34-2.49 m^2 (median 1.86)",
    sex_female_pct = 33.5,
    race_ethnicity = "Not reported.",
    disease_state = "Newly diagnosed advanced non-small cell lung cancer (stage IIIB 14.9%, IV 85.1%; non-squamous adenocarcinoma 76.0%, squamous 24.0%; ECOG 0/1/2 52.5/40.3/7.2%), first-line paclitaxel plus carboplatin (82.8%) or cisplatin (17.2%) every 3 weeks for up to six cycles.",
    dose_range = "Paclitaxel 3-h IV infusion q3w, 200 mg/m^2 (BSA-guided arm, 108 patients) or PK- and neutropenia-guided doses (PK-guided arm, 113 patients).",
    regions = "Switzerland and Germany (CEPAC-TDM, NCT01326767).",
    notes = "Subset of the 365 CEPAC-TDM patients with baseline measurements of all seven blood biomarkers (CRP, IL-6, LDH, CEA, CYFRA 21-1, TSH, NLR); 155 deaths (70%) and 66 right-censored (Results; Table S1). Covariates screened but not retained in the final model are listed in covariatesDataExcluded."
  )

  covariatesDataExcluded <- list(
    IL6 = list(
      description = "Baseline serum interleukin-6 concentration.",
      units = "pg/mL",
      type = "continuous",
      notes = "Significant only in univariate analysis (p < 0.05); not retained after stepwise covariate modeling. Median 2.57 pg/mL (Table S2)."
    ),
    LDH = list(
      description = "Baseline lactate dehydrogenase activity.",
      units = "U/L",
      type = "continuous",
      notes = "Significant only in univariate analysis (p < 0.05); not retained. Median 178 U/L (Table S2)."
    ),
    NLR = list(
      description = "Baseline neutrophil-to-lymphocyte ratio.",
      units = "(ratio)",
      type = "continuous",
      notes = "Significant only in univariate analysis (p < 0.05); not retained. Median 4.23 (Table S2)."
    ),
    TUMSZ = list(
      description = "Baseline tumor size (RECIST 1.1 sum of target-lesion diameters).",
      units = "mm",
      type = "continuous",
      notes = "Significant only in univariate analysis (p < 0.05); not retained. Median 8.40 cm = 84.0 mm (Table S1)."
    ),
    CEA = list(
      description = "Baseline carcinoembryonic antigen concentration.",
      units = "ng/mL",
      type = "continuous",
      notes = "Not significantly associated with OS. Median 7.10 ng/mL (Table S2)."
    ),
    CYFRA21_1 = list(
      description = "Baseline cytokeratin-19 fragment (CYFRA 21-1) concentration.",
      units = "ng/mL",
      type = "continuous",
      notes = "Not significantly associated with OS. Median 6.50 ng/mL (Table S2)."
    ),
    TSH = list(
      description = "Baseline thyroid stimulating hormone.",
      units = "mU/L",
      type = "continuous",
      notes = "Not significantly associated with OS. Median 1.33 mU/L (Table S2)."
    ),
    AGE = list(
      description = "Age at baseline.",
      units = "years",
      type = "continuous",
      notes = "Not significantly associated with OS."
    ),
    SEXF = list(
      description = "Female sex indicator.",
      units = "(binary)",
      type = "binary",
      notes = "Not significantly associated with OS."
    )
  )

  ini({
    llam_haz <- log(0.00010); label("Weibull hazard scale factor lambda_OS (1/h)") # Table 2, Weibull covariate model: lambda_OS = 0.00010 1/h (RSE 7.30%)
    lalfa_haz <- log(1.29); label("Weibull hazard shape factor alpha (unitless)") # Table 2, Weibull covariate model: alpha = 1.29 (RSE 6.90%)
    e_crp_haz <- 0.00744; label("Effect of baseline CRP on the log hazard (per mg/L, centered at 23.5 mg/L)") # Table 2, Weibull covariate model: CRP on lambda = 0.00744 (RSE 16.9%)
    e_rcfb8_haz <- 0.0233; label("Effect of week-8 relative tumor size RS8 on the log hazard (per percent of baseline, centered at 80.0%)") # Table 2, Weibull covariate model: RS8 on lambda = 0.0233 (RSE 26.6%)
  })

  model({
    lam_haz <- exp(llam_haz)
    alfa_haz <- exp(lalfa_haz)

    # RS8 in percent of baseline, as in the paper's NONMEM data column.
    rs8 <- 100 * (1 + RCFB8_SLD)

    # Equation 5 / Appendix S1 $DES. The 1e-6 offset is the NONMEM DEL term,
    # which keeps (lam_haz * t)^(alfa_haz - 1) finite at t = 0.
    del <- 1e-6
    hazard <- lam_haz * alfa_haz * (lam_haz * t + del)^(alfa_haz - 1) *
      exp(e_crp_haz * (CRP - 23.5)) * exp(e_rcfb8_haz * (rs8 - 80.0))
    d/dt(cumhaz) <- hazard
    sur <- exp(-cumhaz)
  })
}
