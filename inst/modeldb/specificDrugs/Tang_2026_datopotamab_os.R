Tang_2026_datopotamab_os <- function() {
  description <- paste0(
    "Cox proportional-hazards exposure-efficacy model for OVERALL ",
    "SURVIVAL in adults with unresectable or metastatic hormone ",
    "receptor-positive, HER2-negative breast cancer receiving ",
    "datopotamab deruxtecan (Dato-DXd) 6 mg/kg every 3 weeks in ",
    "TROPION-Breast01 (Tang 2026, n = 352). The model returns the ",
    "RELATIVE hazard of death hr = exp(-0.0031 * AUC_DATODXD), where ",
    "AUC_DATODXD is the individual Cycle-1 Dato-DXd AUC in ug*day/mL; ",
    "higher exposure is associated with longer survival. Dato-DXd AUC1 ",
    "was the only covariate retained by stepwise selection. NO BASELINE ",
    "HAZARD IS ENCODED: a Cox regression is semiparametric and leaves ",
    "h0(t) unspecified, so only RATIOS of hr between two exposures are ",
    "meaningful. There is no PK layer and no ODE: exposure is supplied ",
    "as data and was derived in the source from empirical Bayes ",
    "estimates of the companion population PK model packaged as ",
    "Hong_2025_datopotamab. No between-subject random effect and no ",
    "residual error are estimated. One of ten exposure-response models ",
    "in the Tang_2026_datopotamab_* family."
  )
  reference <- paste(
    "Tang Z, Lim K, Jiang Y, Denduluri N, Rokutanda N, Pan Y, Hong Y, Ren S,",
    "Vajjah P, Zhou D. Exposure-Response Analyses of Datopotamab Deruxtecan",
    "(Dato-DXd) in Patients with Hormone Receptor (HR) Positive, Human",
    "Epidermal Growth Factor Receptor 2 (HER2)-Negative Breast Cancer.",
    "J Clin Pharmacol. 2026;66(4):245-258. doi:10.1002/jcph.70188.",
    "The coefficient is printed in the Results equation for overall survival.",
    "Individual Cycle-1 AUC values derive from the population PK model of",
    "Hong 2025 (doi:10.1002/psp4.70118); see modellib('Hong_2025_datopotamab').",
    sep = " "
  )
  vignette <- "Tang_2026_datopotamab_exposure_response"
  units <- list(
    time = "n/a (semiparametric Cox regression; the baseline hazard and hence the time scale are left unspecified by the method)",
    dosing = "n/a (no dose events; exposure enters as the AUC_DATODXD data column)",
    concentration = "hr (relative hazard of death, unitless; also lhr, the linear predictor)"
  )

  covariateData <- list(
    AUC_DATODXD = list(
      description = "Model-predicted Cycle-1 (day 0 to day 21) area under the plasma concentration-time curve of intact datopotamab deruxtecan",
      units = "ug*day/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters LINEARLY and UNCENTRED: Tang 2026 Methods give the",
        "exposure term as theta_k * PK with no normalisation (only the",
        "non-exposure continuous covariates are divided by their median).",
        "Derived in the source from empirical Bayes estimates of the",
        "Hong 2025 population PK model, NOT observed. TROPION-Breast01",
        "distribution: median 650, 5th percentile 467 and 95th percentile",
        "864 ug*day/mL (Tang 2026 Results); quartile medians 519, 607, 688",
        "and 779 ug*day/mL (Tang 2026 Table S5)."
      ),
      source_name = "Dato-DXd AUC1"
    )
  )

  covariatesDataExcluded <- list(
    TUMSZ = list(
      description = "Baseline tumor size (sum of target-lesion diameters)",
      units = "mm",
      type = "continuous",
      notes = "Tang 2026 Figure 1B shows an apparent positive relationship between baseline tumor size and overall survival, and Table S3 lists tumor size among the screened continuous covariates, but the Results state 'No other covariates were identified as significant' for overall survival; no coefficient is reported. Retained instead in the companion progression-free-survival model Tang_2026_datopotamab_pfs."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 352L,
    n_studies = 1L,
    age_range = "29-86 years (median 55)",
    weight_range = "35.6-141 kg (median 62)",
    sex_female_pct = 98.6,
    race_ethnicity = "50% White",
    disease_state = "unresectable or metastatic hormone receptor-positive, HER2-negative breast cancer after prior endocrine-based therapy and chemotherapy",
    dose_range = "datopotamab deruxtecan 6 mg/kg IV every 3 weeks",
    regions = "38.1% Europe; also Japan, USA, mainland China and rest of world",
    notes = paste(
      "Efficacy- and PK-evaluable patients of TROPION-Breast01 (data",
      "cut-off 29 April 2024); 8 patients without PK records were",
      "excluded. 83% had received a prior CDK4/6 inhibitor. All patients",
      "received a single dose level, so the exposure range reflects",
      "between-patient PK variability only and Tang 2026 cautions that",
      "the exposure-survival association may be confounded with",
      "disease-related covariates (Table S5)."
    )
  )

  ini({
    # ==================================================================
    # Tang 2026 Results, 'Exposure-Efficacy Analyses':
    #   Overall survival: h(t) = h0(t) * e^(-0.0031 * Dato-DXd AUC1)
    # The coefficient is on the log-hazard scale, the scale a Cox model
    # estimates on, so nothing here is back-solved or digitised. The
    # printed forest-plot hazard ratios are an internal check:
    # exp(-0.0031 * (864 - 650)) = 0.515 and exp(-0.0031 * (467 - 650))
    # = 1.763, against the printed 0.515 and 1.76.
    #
    # NO BASELINE HAZARD IS ENCODED. A Cox regression is semiparametric:
    # h0(t) is left unspecified by the method, so it is not an unreported
    # parameter but a quantity the fit never produced. This model
    # returns the RELATIVE hazard hr only, matching the pattern of
    # Chan_2025_atezolizumab_os.R. The regression was fitted in R 4.3
    # with no variance components, so there is no IIV and no residual
    # error and no observation endpoint is declared.
    # ==================================================================
    e_auc_datodxd_haz <- -0.0031; label("Log hazard ratio for death per 1 ug*day/mL increase in Cycle-1 Dato-DXd AUC (log scale)") # Tang 2026 Results, overall-survival CPH equation: e^(-0.0031 Dato-DXd AUC1)
  })

  model({
    # Linear predictor against a hypothetical patient with AUC1 = 0.
    # Under proportional hazards hr is constant in time, so the subject
    # hazard is h(t) = h0(t) * hr for any baseline h0(t) the user
    # supplies; only ratios of hr between two exposures are meaningful.
    lhr <- e_auc_datodxd_haz * AUC_DATODXD
    hr <- exp(lhr)
  })
}
