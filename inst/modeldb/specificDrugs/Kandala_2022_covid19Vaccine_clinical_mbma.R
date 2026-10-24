Kandala_2022_covid19Vaccine_clinical_mbma <- function() {
  description <- paste0(
    "MBMA (clinical only). COVID-19 vaccine immune correlate of protection: ",
    "the clinical MBMA, the third of three efficacy-scale models in Kandala ",
    "2022. A standalone sigmoid (Hill) relating the convalescent-normalized ",
    "serum neutralizing titre (SN/HCS ratio) to clinical vaccine efficacy ",
    "(VE %), fitted to clinical data alone. The Hill coefficient was fixed to ",
    "1 (its log-scale 95% CI contained 0). The 50% protective SN titre in ",
    "humans is 21% of the geometric-mean human convalescent titre (Figure 5, ",
    "ET50 = 0.21). Weighted (trial-size) generalized least squares with an ",
    "additive residual; no random effects, so it simulates typical vaccine ",
    "efficacy and is NOT an individual-level model. The two sibling efficacy ",
    "models are modellib('Kandala_2022_covid19Vaccine_calibrated_mbma') and ",
    "modellib('Kandala_2022_covid19Vaccine_joint_mbma')."
  )

  reference <- paste(
    "Kandala B, Plock N, Chawla A, Largajolli A, Robey S, Watson K, et al.",
    "Accelerating model-informed decisions for COVID-19 vaccine candidates",
    "using a model-based meta-analysis approach.",
    "EBioMedicine. 2022;84:104264. doi:10.1016/j.ebiom.2022.104264.",
    sep = " "
  )

  vignette <- "Kandala_2022_covid19Vaccine"

  units <- list(
    time = paste0(
      "not applicable -- static study-level meta-analysis with no time course."
    ),
    dosing = paste0(
      "not applicable -- vaccine exposure enters through the SN_HCS_RATIO ",
      "covariate, not through dosing records."
    ),
    concentration = paste0(
      "ve is the predicted clinical vaccine efficacy in percent. The output ",
      "Cc is unused."
    )
  )

  covariateData <- list(
    SN_HCS_RATIO = list(
      description = paste0(
        "Study-arm ratio of the vaccine SN titre to the same-assay geometric-",
        "mean human convalescent serum (HCS) titre. Dimensionless; the sole ",
        "driver of the clinical sigmoid."
      ),
      units = "(fold; dimensionless ratio)",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Kandala 2022 Equation 12 and Figure 5 x-axis, computed from Table 3 ",
        "columns 'SN Titre (GM)' and 'HCS (GM)'. Observed clinical values ",
        "span roughly 0.045 to 3.9; predictions outside that range ",
        "extrapolate beyond the calibration data."
      ),
      source_name = "SN/HCS"
    )
  )

  population <- list(
    species = "human",
    n_studies = 6L,
    age_range = paste0(
      "Adults; Ad26.COV2.S and BNT162b2 each contribute separate younger-",
      "adult and older-adult (65-85 years) efficacy strata."
    ),
    disease_state = paste0(
      "PCR-confirmed symptomatic COVID-19 of any severity; efficacy is the ",
      "relative incidence-rate reduction versus placebo."
    ),
    dose_range = paste0(
      "Six vaccine candidates (BNT162b2, mRNA-1273, Ad26.COV2.S, 2019nCoV, ",
      "Gam-COVID-Vac, ChAdOx1 nCoV-19); Phase 1/2 immunogenicity matched to ",
      "Phase 3 efficacy from identical regimens (Kandala 2022 Table 3)."
    ),
    regions = "International (published Phase 1/2 and Phase 3 trials)",
    notes = paste0(
      "Model-based meta-analysis: the unit of observation is a trial-level ",
      "efficacy with its matched SN/HCS ratio. Efficacy observations were ",
      "weighted by the number of participants in each trial. External ",
      "validation (BBV152, CoronaVac) and delta-variant predictions are in ",
      "the vignette; this model reproduces the Table 4 external predictions ",
      "exactly from the typical-value curve."
    )
  )

  # ==========================================================================
  # ini(): Kandala 2022 Figure 5 parameter table (clinical MBMA). The Hill
  # coefficient was fixed to 1 because the 95% CI of its log contained 0
  # (log-Hill estimate 0.29, 95% CI [-0.88, 1.48]; exp(0) = 1). VEmax is a 0-1
  # fraction (estimate 1.00). Reported 95% CIs are in the trailing comments.
  # The additive residual scale is not tabulated, so it is a placeholder.
  # ==========================================================================
  ini({
    vemax <- 1.00
    label("Maximum vaccine efficacy (fraction)")  # Kandala 2022 Figure 5 table, row 'VEmax': 1.00, 95% CI [0.94, 1.07]

    et50 <- 0.21
    label("SN/HCS titre ratio for 50% clinical vaccine efficacy in humans (dimensionless)")  # Kandala 2022 Figure 5 table, row 'ET50,human': 0.21, 95% CI [0.15, 0.30]

    hill <- fixed(1)
    label("Hill coefficient of the clinical-efficacy sigmoid (unitless)")  # Kandala 2022 Results: log-Hill 0.29, 95% CI [-0.88, 1.48] contained 0, so Hill fixed to exp(0) = 1

    addSd <- fixed(0.001)
    label("Placeholder additive residual SD on VE (percent); NOT a published value")  # not from source -- clinical-model residual scale is unreported; see vignette Assumptions and deviations
  })

  model({
    # ======================================================================
    # Clinical sigmoid (Kandala 2022 Equation 12):
    #   VE(%) = 100 * VEmax * (SN/HCS)^g / ((SN/HCS)^g + ET50^g)
    # VEmax is a 0-1 fraction, so the factor of 100 puts VE on the percent
    # scale the paper plots (Figure 5).
    # ======================================================================
    frac <- SN_HCS_RATIO^hill / (SN_HCS_RATIO^hill + et50^hill)
    ve <- 100 * vemax * frac

    ve ~ add(addSd)
  })
}
