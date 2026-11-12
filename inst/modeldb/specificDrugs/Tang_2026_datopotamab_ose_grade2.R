Tang_2026_datopotamab_ose_grade2 <- function() {
  description <- paste0(
    "Logistic-regression exposure-safety model for grade 2 or worse ",
    "ocular surface event (OSE; MedDRA v26.0 preferred-term list) in ",
    "adults with breast or non-small cell lung cancer treated with ",
    "datopotamab deruxtecan (Dato-DXd) 0.27-10 mg/kg every 3 weeks ",
    "(Tang 2026, n = 1081 pooled from TROPION-PanTumor01, -Lung01, ",
    "-Lung05 and -Breast01). The probability is expit(-4.748 + ",
    "0.07011 * CAV), where CAV is the individual average plasma ",
    "concentration of INTACT Dato-DXd (ug/mL) up to the end of the ",
    "cycle of the event. No covariate other than exposure was ",
    "retained. There is no PK layer and no ODE: exposure is supplied ",
    "as data and was derived in the source from the companion ",
    "population PK model packaged as Hong_2025_datopotamab. No ",
    "between-subject random effect and no residual error are ",
    "estimated. One of ten exposure-response models in the ",
    "Tang_2026_datopotamab_* family."
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
    dosing = "n/a (no dose events; exposure enters as the CAV covariate column)",
    concentration = "prob_ocular_surface_event_grade2 (probability of grade 2 or worse ocular surface event, 0-1)"
  )

  covariateData <- list(
    CAV = list(
      description = "Individual average plasma concentration of INTACT datopotamab deruxtecan, cumulative AUC divided by the time from the first dose to the end of the cycle of the event (or, without an event, to the cycle of the last treatment or the data cut-off, whichever came first)",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Intact conjugate (Dato-DXd), in ug/mL; the companion Tang_2026_datopotamab_teae_grade3 and _sae models carry the PAYLOAD DXd Cavg in ng/mL under the same column name. Enters LINEARLY and UNCENTRED. Derived in the source from empirical Bayes estimates of the Hong 2025 population PK model (modellib('Hong_2025_datopotamab'), Cc in ug/mL), accounting for dose interruptions and reductions; endpoint-specific because the averaging window ends at the event. Tang 2026 Figure 2D shows the analysis range spanning roughly 15-60 ug/mL with medians near 22, 31 and 42 ug/mL at 4, 6 and 8 mg/kg. Tang 2026 cautions that Cavg derived this way may be biased.",
      source_name = "Dato-DXd Cavg"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1081L,
    n_studies = 4L,
    age_range = "26-86 years (median 60)",
    weight_range = "35.6-156 kg (median 64.2)",
    sex_female_pct = 67.2,
    race_ethnicity = "47.2% White, 40% Asian, 1.8% Black, 7.3% multiple or other, 2.9% unknown or not reported (Table S2)",
    disease_state = "advanced or metastatic non-small cell lung cancer (59.5%), HR-positive/HER2-negative breast cancer (36.4%) or triple-negative breast cancer (4.1%)",
    dose_range = "datopotamab deruxtecan 0.27-10 mg/kg IV every 3 weeks; 919 of 1081 (85%) at 6 mg/kg, 50 at 4 mg/kg and 82 at 8 mg/kg",
    regions = "USA 29.7%, Europe 27.8%, Japan 17.8%, mainland China 4.5%, rest of world 20.2%",
    notes = "Safety analysis set (data cut-off 29 April 2024). Observed incidence 8% (4/50), 10% (88/919) and 13% (11/82) at 4, 6 and 8 mg/kg (Table 2)."
  )

  ini({
    # Tang 2026 Supporting Information Table S4, 'OSE Grade 2+' rows,
    # printed on the log-odds scale; nothing is back-solved.
    logit_ref <- -4.74822; label("Logit of the event probability at zero Dato-DXd Cavg (unitless logit)") # Tang 2026 Table S4, OSE Grade 2+ (Intercept): -4.74822 [-5.45; -4.04]
    e_cav_logit <- 0.070107; label("Log-odds of grade 2 or worse ocular surface event per 1 ug/mL increase in Dato-DXd Cavg (unitless logit)") # Tang 2026 Table S4, OSE Grade 2+ Dato-DXd Cavg: 0.070107 [0.0528; 0.0874]

    # Bernoulli-likelihood logistic regression (R 4.3): no variance
    # components. The placeholder additive term exists only because
    # rxode2 requires an observation declaration.
    addSd_prob_ocular_surface_event_grade2 <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)") # not from source; see vignette Assumptions and deviations
  })

  model({
    logit_ocular_surface_event_grade2 <- logit_ref + e_cav_logit * CAV
    prob_ocular_surface_event_grade2 <- expit(logit_ocular_surface_event_grade2)

    prob_ocular_surface_event_grade2 ~ add(addSd_prob_ocular_surface_event_grade2)
  })
}
