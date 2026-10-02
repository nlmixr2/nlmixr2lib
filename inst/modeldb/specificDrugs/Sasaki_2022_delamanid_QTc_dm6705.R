Sasaki_2022_delamanid_QTc_dm6705 <- function() {
  description <- "Linear mixed-effects concentration-QTc model relating the time-matched change from baseline in the Bazett-corrected QT interval (DeltaQTcB, ms) to the plasma concentration of DM-6705, the major metabolite of delamanid, in children and adolescents (0.67-17 years) with multidrug-resistant tuberculosis: DeltaQTcB = (theta0 + eta0) + (theta1 + eta1) * C + theta2 * (QTcB0 - mean QTcB0), with additive random effects on the intercept and slope. The DM-6705 slope (0.0613 ms per ng/mL, 90% CI 0.016-0.107) was the significant concentration-QTc relationship of the analysis. PD-only; the concentration is supplied as a covariate, for example from Sasaki_2022_delamanid."
  reference <- paste(
    "Sasaki T, Svensson EM, Wang X, Wang Y, Hafkin J, Karlsson MO, Mallikaarjun S.",
    "Population Pharmacokinetic and Concentration-QTc Analysis of Delamanid in",
    "Pediatric Participants with Multidrug-Resistant Tuberculosis.",
    "Antimicrob Agents Chemother. 2022;66(2):e01608-21.",
    "doi:10.1128/aac.01608-21"
  )
  vignette <- "Sasaki_2022_delamanid"
  units <- list(
    time = "h",
    dosing = "(none; PD-only model fed by an external DM-6705 plasma-concentration covariate)",
    concentration = "(observation QTc is the change from baseline in the Bazett-corrected QT interval, ms; driving covariate CP_DM6705_NGML is in ng/mL)"
  )

  covariateData <- list(
    CP_DM6705_NGML = list(
      description = "Plasma concentration of DM-6705 at the time of each ECG observation",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying. In Sasaki 2022 this was the observed, time-matched DM-6705 concentration (UPLC-MS/MS, LLOQ 1.00 ng/mL); for simulation it can be taken from the Cc_dm6705 output of Sasaki_2022_delamanid. The slope is reported in ms per ng/mL (Table S5), so no unit rescaling is applied. Set to 0 to recover the drug-free prediction.",
      source_name = "C (DM-6705 concentration)"
    ),
    QTC_BL = list(
      description = "Subject's baseline Bazett-corrected QT interval (QTcB), time-fixed",
      units = "ms",
      type = "continuous",
      reference_category = NULL,
      notes = "Bazett correction. Baseline is the triplicate-averaged QTcB at the time-matched day -1 timepoint (Sasaki 2022 Methods and Table S1). Enters as e_qtc_bl_e0 * (QTC_BL - 421). The paper defines the centering constant QTc0 as the mean of all baseline QTc values but does not print it; 421 ms is derived from the baseline QTcB-versus-RR regression of Figure S5 (QTcB = 0.03 * RR + 403.167, a least-squares line that passes through the mean point) at the mean RR of about 600 ms implied by the Discussion's baseline heart rate of 'approximately 100 bpm'. With theta2 = 0.0309, a 10-ms error in this constant shifts the prediction by 0.3 ms.",
      source_name = "QTc_ij0"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 37L,
    n_studies = 2L,
    age_range = "0.67-17 years",
    age_median = "mean (SD) 6.36 (5.19) years; median not reported",
    weight_range = "not reported; mean (SD) by age group 39 (4.59), 24.9 (6.79), 14.2 (3.2) and 9.76 (1.83) kg",
    weight_median = "mean (SD) 19.2 (11.6) kg; median not reported",
    sex_female_pct = 51.4,
    race_ethnicity = c(Asian = 67.6, Black = 5.41, Other = 27),
    disease_state = "Multidrug-resistant tuberculosis, on an optimized background regimen that could include QT-prolonging drugs (fluoroquinolones, clofazimine, macrolides)",
    dose_range = "Oral delamanid 5 mg QD to 100 mg BID by age and weight group (Table S1)",
    regions = "Philippines (67.6%), South Africa (32.4%)",
    notes = "354 QT measurements with time-matched delamanid / DM-6705 concentrations from 37 participants of trials 232 and 233 (Sasaki 2022 Results). No placebo arm, so the endpoint is DeltaQTcB, not a placebo-corrected DeltaDeltaQTcB. Bazett rather than Fridericia correction was chosen because QTcB was less dependent on RR at baseline (Figure S5)."
  )

  ini({
    # Linear mixed-effects model (Sasaki 2022 Methods):
    #   DeltaQTc_ijk = (theta0 + eta0_i) + (theta1 + eta1_i) * C_ijk
    #                  + theta2 * (QTc_ij0 - QTc0) + eps
    # Fixed effects from Table S5 ('Linear Mixed Effects Model for
    # DM-6705/DeltaQTcB').
    e0 <- 0.923
    label("Population mean intercept theta0 on DeltaQTcB (ms)") # Table S5 'Intercept (ms)' = 0.923 (SE 1.61; 90% CI -2.22, 4.07)
    slope <- 0.0613
    label("DM-6705 concentration-DeltaQTcB slope theta1 (ms per ng/mL)") # Table S5 'Slope (ms/[ng/mL])' = 0.0613 (SE 0.0231; 90% CI 0.016, 0.107)
    e_qtc_bl_e0 <- 0.0309
    label("Effect of centered baseline QTcB on DeltaQTcB (ms per ms)") # Table S5 'Baseline QTcB effect' = 0.0309 (SE 0.0728)

    # Random effects are normally distributed and additive (Methods);
    # Table S5 reports variances and no covariance.
    etae0 ~ 59.9 # Table S5 'Random effect for the intercept' = 59.9 (variance)
    etaslope ~ 0.00446 # Table S5 'Random effect for the slope' = 0.00446 (variance)

    addSd <- 13.190906
    label("Additive residual error on DeltaQTcB (ms)") # Table S5 'Residual Variability (shown as variance)' = 174; SD = sqrt(174)
  })

  model({
    e0_i <- e0 + etae0
    slope_i <- slope + etaslope

    # Centering constant QTc0 (mean baseline QTcB), derived from Figure S5;
    # see covariateData$QTC_BL$notes.
    qtc_bl_ref <- 421

    QTc <- e0_i + slope_i * CP_DM6705_NGML + e_qtc_bl_e0 * (QTC_BL - qtc_bl_ref)

    QTc ~ add(addSd)
  })
}
