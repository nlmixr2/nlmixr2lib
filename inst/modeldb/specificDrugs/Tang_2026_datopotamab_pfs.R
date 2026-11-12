Tang_2026_datopotamab_pfs <- function() {
  description <- paste0(
    "Cox proportional-hazards model for PROGRESSION-FREE SURVIVAL (per ",
    "independent central review) in adults with unresectable or ",
    "metastatic hormone receptor-positive, HER2-negative breast cancer ",
    "receiving datopotamab deruxtecan (Dato-DXd) 6 mg/kg every 3 weeks ",
    "in TROPION-Breast01 (Tang 2026, n = 352). The model returns the ",
    "RELATIVE hazard of progression or death hr = exp(0.28 * TUMSZ / ",
    "55), where TUMSZ is baseline tumor size in mm normalised by the ",
    "cohort median of 55 mm. NO EXPOSURE TERM: Dato-DXd Cycle-1 AUC ",
    "entered in forward addition but was removed in backward ",
    "elimination (p = 0.0011 against the p < 0.001 retention ",
    "criterion), leaving baseline tumor size as the only covariate; ",
    "this is a baseline-prognostic model, packaged because it is the ",
    "final model the paper reports for this endpoint. NO BASELINE ",
    "HAZARD IS ENCODED: a Cox regression is semiparametric and leaves ",
    "h0(t) unspecified, so only RATIOS of hr between two covariate ",
    "settings are meaningful. No between-subject random effect and no ",
    "residual error are estimated. One of ten exposure-response models ",
    "in the Tang_2026_datopotamab_* family."
  )
  reference <- paste(
    "Tang Z, Lim K, Jiang Y, Denduluri N, Rokutanda N, Pan Y, Hong Y, Ren S,",
    "Vajjah P, Zhou D. Exposure-Response Analyses of Datopotamab Deruxtecan",
    "(Dato-DXd) in Patients with Hormone Receptor (HR) Positive, Human",
    "Epidermal Growth Factor Receptor 2 (HER2)-Negative Breast Cancer.",
    "J Clin Pharmacol. 2026;66(4):245-258. doi:10.1002/jcph.70188.",
    "The coefficient is printed in the Results equation for progression-free",
    "survival; the median normalisation follows the Methods CPH equation.",
    sep = " "
  )
  vignette <- "Tang_2026_datopotamab_exposure_response"
  units <- list(
    time = "n/a (semiparametric Cox regression; the baseline hazard and hence the time scale are left unspecified by the method)",
    dosing = "n/a (no dose events; the model has no exposure term)",
    concentration = "hr (relative hazard of progression or death, unitless; also lhr, the linear predictor)"
  )

  covariateData <- list(
    TUMSZ = list(
      description = "Baseline tumor size (sum of target-lesion diameters)",
      units = "mm",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters LINEARLY after division by the TROPION-Breast01 median of",
        "55 mm, per the Tang 2026 Methods CPH equation (theta_j * x_j /",
        "x_median). The Results print the coefficient as",
        "'e^(0.28 Baseline tumor size)', and the published hazard ratios",
        "fix the scaling: exp(0.28 * (138/55 - 1)) = 1.53 and",
        "exp(0.28 * (17/55 - 1)) = 0.824 against the printed 1.52 and",
        "0.825 for the 95th and 5th percentiles; read per mm (no",
        "normalisation) the same coefficient would give hazard ratios of",
        "order 10^10. TROPION-Breast01 range 10-338 mm (Table S2)."
      ),
      source_name = "Baseline tumor size"
    )
  )

  covariatesDataExcluded <- list(
    AUC_DATODXD = list(
      description = "Model-predicted Cycle-1 area under the plasma concentration-time curve of intact datopotamab deruxtecan",
      units = "ug*day/mL",
      type = "continuous",
      notes = "Selected in forward addition but removed in backward elimination (p = 0.0011 > 0.001); Tang 2026 reports no coefficient for it in the progression-free-survival model. Figure 1D shows an apparent positive Kaplan-Meier relationship by AUC1 quartile, and the Discussion attributes part of it to confounding with baseline tumor size and albumin (Table S5)."
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
    notes = "Efficacy- and PK-evaluable patients of TROPION-Breast01 (data cut-off 29 April 2024). Baseline tumor size median 55 mm, 5th percentile 17 mm, 95th percentile 138 mm (Tang 2026 Results)."
  )

  ini({
    # ==================================================================
    # Tang 2026 Results, 'Exposure-Efficacy Analyses':
    #   progression-free survival: h(t) = h0(t) * e^(0.28 Baseline tumor size)
    # with the Methods equation dividing continuous covariates by their
    # median (55 mm in TROPION-Breast01). The coefficient is printed to
    # two significant figures only; the two published hazard ratios
    # back-solve to 0.2775 and 0.2784, both of which round to 0.28, so
    # the printed value is used unchanged.
    #
    # NO BASELINE HAZARD IS ENCODED; this model returns the RELATIVE
    # hazard hr only. No IIV, no residual error and no observation
    # endpoint (Cox regression in R 4.3).
    # ==================================================================
    tumsz_ref <- fixed(55); label("Median baseline tumor size used to normalise TUMSZ (mm)") # Tang 2026 Results: 'compared to the median tumor size of 55 mm'
    e_tumsz_haz <- 0.28; label("Log hazard ratio for progression or death per unit of median-normalised baseline tumor size (log scale)") # Tang 2026 Results, PFS CPH equation: e^(0.28 Baseline tumor size)
  })

  model({
    # The 'reference' patient for ratios is the median patient
    # (TUMSZ = 55 mm, lhr = 0.28); hr itself is relative to a
    # hypothetical TUMSZ = 0.
    lhr <- e_tumsz_haz * TUMSZ / tumsz_ref
    hr <- exp(lhr)
  })
}
