Tang_2026_datopotamab_dose_interruption <- function() {
  description <- paste0(
    "Logistic-regression exposure-safety model for a treatment-emergent ",
    "adverse event leading to dose interruption in adults with breast or ",
    "lung cancer treated with datopotamab deruxtecan (Dato-DXd) every 3 ",
    "weeks (Tang 2026; TROPION-PanTumor01 and TROPION-Breast01 only, the ",
    "two studies that recorded this endpoint). The probability is ",
    "expit(-2.733 + 0.002195 * AUC_DATODXD), where AUC_DATODXD is the ",
    "individual Cycle-1 intact Dato-DXd AUC in ug*day/mL. No covariate ",
    "other than exposure was retained. There is no PK layer and no ODE: ",
    "exposure is supplied as data and was derived in the source from ",
    "the companion population PK model packaged as ",
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
    "Coefficients are transcribed from Supporting Information Table S4.",
    sep = " "
  )
  vignette <- "Tang_2026_datopotamab_exposure_response"
  units <- list(
    time = "n/a (static landmark exposure-response model; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the AUC_DATODXD covariate column)",
    concentration = "prob_dose_interruption (probability of a TEAE leading to dose interruption, 0-1)"
  )

  covariateData <- list(
    AUC_DATODXD = list(
      description = "Model-predicted Cycle-1 (day 0 to day 21) area under the plasma concentration-time curve of intact datopotamab deruxtecan",
      units = "ug*day/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters LINEARLY and UNCENTRED. Derived in the source from empirical Bayes estimates of the Hong 2025 population PK model (modellib('Hong_2025_datopotamab')). Tang 2026 Figure 2C shows medians near 440, 650 and 910 ug*day/mL at 4, 6 and 8 mg/kg.",
      source_name = "Dato-DXd AUC1"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 647L,
    n_studies = 2L,
    disease_state = "advanced solid tumours: non-small cell lung cancer, HR-positive/HER2-negative breast cancer and triple-negative breast cancer",
    dose_range = "datopotamab deruxtecan 0.27-10 mg/kg IV every 3 weeks (TROPION-PanTumor01) and 6 mg/kg (TROPION-Breast01)",
    notes = "Subset of the 1081-patient pooled safety set restricted to TROPION-PanTumor01 (295) and TROPION-Breast01 (352), the two studies in which a TEAE leading to dose interruption was analysed (Tang 2026 Table 2 footnote); 647 is the sum of the two study counts in Table S2, the paper does not print the subset size. Observed incidence 8% (4/50), 24% (117/485) and 35% (29/82) at 4, 6 and 8 mg/kg (Table 2). Pooled-set demographics: see Tang_2026_datopotamab_teae_grade3."
  )

  ini({
    # Tang 2026 Supporting Information Table S4, 'Dose Interruption'
    # rows, printed on the log-odds scale.
    logit_ref <- -2.73329; label("Logit of the event probability at zero Cycle-1 Dato-DXd AUC (unitless logit)") # Tang 2026 Table S4, Dose Interruption (Intercept): -2.73329 [-3.36; -2.11]
    e_auc_datodxd_logit <- 0.002195; label("Log-odds of a TEAE leading to dose interruption per 1 ug*day/mL increase in Cycle-1 Dato-DXd AUC (unitless logit)") # Tang 2026 Table S4, Dose Interruption Dato-DXd AUC1: 0.002195 [0.00138; 0.00300]

    addSd_prob_dose_interruption <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_dose_interruption <- logit_ref + e_auc_datodxd_logit * AUC_DATODXD
    prob_dose_interruption <- expit(logit_dose_interruption)

    prob_dose_interruption ~ add(addSd_prob_dose_interruption)
  })
}
