Sasaki_2022_delamanid_QTc_parent <- function() {
  description <- "Linear mixed-effects concentration-QTc model relating the time-matched change from baseline in the Bazett-corrected QT interval (DeltaQTcB, ms) to the plasma concentration of parent delamanid in children and adolescents (0.67-17 years) with multidrug-resistant tuberculosis: DeltaQTcB = (theta0 + eta0) + (theta1 + eta1) * C + theta2 * (QTcB0 - mean QTcB0), with additive random effects on the intercept and slope. The delamanid slope (0.00792 ms per ng/mL, 90% CI -0.00132 to 0.0172) was not statistically significant; the companion DM-6705 model (Sasaki_2022_delamanid_QTc_dm6705) is the one the paper used for its QTc projections. PD-only; the concentration is supplied as a covariate, for example from Sasaki_2022_delamanid."
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
    dosing = "(none; PD-only model fed by an external delamanid plasma-concentration covariate)",
    concentration = "(observation QTc is the change from baseline in the Bazett-corrected QT interval, ms; driving covariate CP_DELAMANID_NGML is in ng/mL)"
  )

  covariateData <- list(
    CP_DELAMANID_NGML = list(
      description = "Plasma concentration of delamanid at the time of each ECG observation",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying. In Sasaki 2022 this was the observed, time-matched delamanid concentration (UPLC-MS/MS, LLOQ 1.00 ng/mL); for simulation it can be taken from the Cc output of Sasaki_2022_delamanid. The slope is reported in ms per ng/mL (Table S4), so no unit rescaling is applied. Set to 0 to recover the drug-free prediction.",
      source_name = "C (delamanid concentration)"
    ),
    QTC_BL = list(
      description = "Subject's baseline Bazett-corrected QT interval (QTcB), time-fixed",
      units = "ms",
      type = "continuous",
      reference_category = NULL,
      notes = "Bazett correction. Baseline is the triplicate-averaged QTcB at the time-matched day -1 timepoint (Sasaki 2022 Methods and Table S1). Enters as e_qtc_bl_e0 * (QTC_BL - 421). The paper defines the centering constant QTc0 as the mean of all baseline QTc values but does not print it; 421 ms is derived from the baseline QTcB-versus-RR regression of Figure S5 (QTcB = 0.03 * RR + 403.167, a least-squares line that passes through the mean point) at the mean RR of about 600 ms implied by the Discussion's baseline heart rate of 'approximately 100 bpm'. With theta2 = 0.0318, a 10-ms error in this constant shifts the prediction by 0.3 ms.",
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
    notes = "354 QT measurements with time-matched delamanid / DM-6705 concentrations from 37 participants of trials 232 and 233 (Sasaki 2022 Results). No placebo arm, so the endpoint is DeltaQTcB, not a placebo-corrected DeltaDeltaQTcB."
  )

  ini({
    # Linear mixed-effects model (Sasaki 2022 Methods):
    #   DeltaQTc_ijk = (theta0 + eta0_i) + (theta1 + eta1_i) * C_ijk
    #                  + theta2 * (QTc_ij0 - QTc0) + eps
    # Fixed effects from Table S4 ('Linear Mixed Effects Model for
    # Delamanid/DeltaQTcB').
    e0 <- 1.47
    label("Population mean intercept theta0 on DeltaQTcB (ms)") # Table S4 'Intercept (ms)' = 1.47 (SE 1.67; 90% CI -1.80, 4.75)
    slope <- 0.00792
    label("Delamanid concentration-DeltaQTcB slope theta1 (ms per ng/mL)") # Table S4 'Slope (ms/[ng/mL])' = 0.00792 (SE 0.00471; 90% CI -0.00132, 0.0172)
    e_qtc_bl_e0 <- 0.0318
    label("Effect of centered baseline QTcB on DeltaQTcB (ms per ms)") # Table S4 'Baseline QTcB effect' = 0.0318 (SE 0.0739)

    # Random effects are normally distributed and additive (Methods);
    # Table S4 reports variances and no covariance.
    etae0 ~ 63.5 # Table S4 'Random effect for the intercept' = 63.5 (variance)
    etaslope ~ 0.0000937 # Table S4 'Random effect for the slope' = 0.0000937 (variance)

    addSd <- 13.416408
    label("Additive residual error on DeltaQTcB (ms)") # Table S4 'Residual Variability (shown as variance)' = 180; SD = sqrt(180)
  })

  model({
    e0_i <- e0 + etae0
    slope_i <- slope + etaslope

    # Centering constant QTc0 (mean baseline QTcB), derived from Figure S5;
    # see covariateData$QTC_BL$notes.
    qtc_bl_ref <- 421

    QTc <- e0_i + slope_i * CP_DELAMANID_NGML + e_qtc_bl_e0 * (QTC_BL - qtc_bl_ref)

    QTc ~ add(addSd)
  })
}
