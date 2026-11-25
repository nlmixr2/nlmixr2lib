Schlachter_2026b_atogepant <- function() {
  description <- paste(
    "Exposure-response model for monthly migraine days (MMD) under oral",
    "atogepant, a calcitonin gene-related peptide (CGRP) receptor antagonist,",
    "in adults with episodic or chronic migraine pooled across one phase 2b/3",
    "and three phase 3 studies (Schlachter 2026). The predicted MMD is the",
    "individual baseline MMD multiplied by one minus a combined placebo plus",
    "atogepant effect that approaches its plateau exponentially in time with",
    "half-time T50. The combined effect is built on a shifted-logit scale",
    "(back-transformed as 2 * expit(eff) - 1, so it lies in (-1, 1)) as an",
    "additive placebo effect plus an Emax function of the individual",
    "steady-state atogepant trough concentration (nM) from the published",
    "population PK model. Baseline MMD lowers the atogepant Emax and lengthens",
    "T50; failure of two to four classes of conventional oral migraine",
    "preventives lowers the combined effect and shortens T50; Asian race",
    "lowers the combined effect. Residual error is additive on MMD with",
    "separate magnitudes for episodic and chronic migraine. The source",
    "compressed simulated MMD into the 0-28 day range after simulation",
    "(compressed normal distribution); that clamp is not part of the",
    "likelihood and is applied in the vignette, not in the model."
  )
  reference <- paste(
    "Schlachter L, Beck D, Boinpally RR, Stodtmann S. Exposure-Response",
    "Modeling of Monthly Migraine Days for Efficacy of Atogepant in Patients",
    "With Episodic or Chronic Migraine. CPT Pharmacometrics Syst Pharmacol.",
    "2026;15(1):e70154. doi:10.1002/psp4.70154.",
    "Individual atogepant trough concentrations come from the population PK",
    "model of Schlachter et al. (2026) Clin Pharmacokinet 65(1):149-164,",
    "doi:10.1007/s40262-025-01566-5; see",
    "modellib('Schlachter_2026_atogepant').",
    sep = " "
  )
  vignette <- "Schlachter_2026b_atogepant"

  # Non-canonical residual-SD names carried by this paper: the additive
  # residual error has separate magnitudes for episodic and chronic migraine
  # (Table 4 'Additive error episodic' / 'Additive error chronic'; Data S1
  # EPS(1) / EPS(2) selected by IND). Same shape as Pohl_2022_linzagolix.R.
  paper_specific_residual_sds <- c("addSdEpisodic", "addSdChronic")

  units <- list(
    time = "month",
    dosing = "n/a (exposure-response model; no dose events)",
    concentration = "nM",
    response = "migraine days/month"
  )

  covariateData <- list(
    CTROUGH = list(
      description = "Individual predicted steady-state atogepant minimum (trough) plasma concentration",
      units = "nM",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "TOTAL atogepant plasma concentration at the end of the steady-state",
        "dosing interval, an empirical-Bayes prediction from the population PK",
        "model of Schlachter 2026 (Clin Pharmacokinet; modellib",
        "'Schlachter_2026_atogepant', which predicts in ng/mL -- divide by",
        "about 0.6035 ng/mL per nM, atogepant molar mass 603.5 g/mol, to",
        "obtain nM). Time-fixed per subject. 0 for placebo patients.",
        "Enters through a hyperbolic Emax on the shifted-logit scale.",
        "Table 2 medians (5th, 95th percentiles): 10 mg QD 4.34 (1.94,",
        "10.2), 30 mg QD 14.7 (7.32, 34.3), 60 mg QD 28.9 (10.2, 67.0),",
        "30 mg BID 72.9 (33.1, 236), 60 mg BID 167 (69.4, 480) nM."
      ),
      source_name = "CMIN"
    ),
    PLACEBO = list(
      description = "Placebo-arm membership indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (atogepant arm)",
      notes = paste(
        "PLACEBO = 1 when the Data S1 column TRT is 0. Selects the",
        "placebo-only branch of the combined effect, which carries the",
        "placebo-effect random effect (ETA(1)) and no drug term; the",
        "atogepant branch carries the Emax random effect (ETA(2)) and no",
        "placebo random effect (Data S1 'IF (TRT .EQ. 0) EFF = LGPEFF +",
        "ETA(1) + EFFCOV')."
      ),
      source_name = "TRT"
    ),
    MIGRAINE_DAYS_BL = list(
      description = "Baseline number of monthly migraine days",
      units = "days/month",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed per subject. Plays two roles: it is the individual",
        "baseline that the combined effect scales (the prediction at time 0",
        "equals MIGRAINE_DAYS_BL exactly), and, centred at the pooled",
        "median 9.00 days/month (Table 1; Data S1 'BSL - 9.00'), it is a",
        "linear covariate on the atogepant Emax and an exponential",
        "covariate on T50."
      ),
      source_name = "BSL"
    ),
    MIGRAINE_ORAL_PREVENTIVE_IR = list(
      description = "Two to four classes of conventional oral migraine preventive treatments failed",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not previously failed)",
      notes = paste(
        "MIGRAINE_ORAL_PREVENTIVE_IR = 1 when the Data S1 column TRTF is 2",
        "('previously failed treatments' = yes). In the source cohort it",
        "coincides exactly with the ELEVATE study (all 280 ELEVATE patients",
        "and no patient of the other three studies; Table 1). Lowers the",
        "combined placebo + atogepant effect additively on the shifted-logit",
        "scale and multiplies T50 by 0.628."
      ),
      source_name = "TRTF"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian)",
      notes = "RACE_ASIAN = 1 when the Data S1 column ASIAN is 2. Lowers the combined placebo + atogepant effect additively on the shifted-logit scale.",
      source_name = "ASIAN"
    ),
    DIS_MIGRAINE_CHRONIC = list(
      description = "Chronic-migraine patient indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (episodic migraine)",
      notes = paste(
        "DIS_MIGRAINE_CHRONIC = 1 when the Data S1 column IND is 2 (chronic",
        "migraine; the PROGRESS study). Selects the chronic-migraine",
        "additive residual error. Indication was also screened as a",
        "covariate on the combined effect, Emax and T50 and not retained",
        "(Table S1B); in the structural model chronic and episodic migraine",
        "differ only through MIGRAINE_DAYS_BL."
      ),
      source_name = "IND"
    )
  )

  covariatesDataExcluded <- list(
    PRIOR_MIGRAINE_PREVENTIVE = list(
      description = "Prior preventive migraine medication (yes/no)",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened in the exploratory stepwise search; added on the combined",
        "effect in forward step 2 (delta OFV 9.35) but removed in backward",
        "elimination (p = 0.00223 > 0.001; Table S1C, Discussion). Not",
        "retained. Source column PRIMED."
      )
    ),
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened on the combined effect, Emax and T50 (Table S1B); not significant. Source column WTKG."
    ),
    SEXF = list(
      description = "Biological sex",
      units = "(binary)",
      type = "binary",
      notes = "Screened on the combined effect, Emax and T50 (Table S1B); not significant. Source column SEX."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on the combined effect, Emax and T50 (Table S1B); not significant. Source column AGE."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 2544L,
    n_studies = 4L,
    age_range = "18-74 years",
    age_median = "41 years",
    weight_range = "35.8-226 kg",
    weight_median = "76.2 kg",
    sex_female_pct = 88.1,
    race_ethnicity = c(
      White = 78.4,
      Black = 11.7,
      Asian = 8.1,
      AmIndAlaskaNative = 0.5,
      NativeHawaiianPacIsl = 0.1,
      Multiple = 1.2
    ),
    disease_state = paste(
      "Adults with episodic migraine (1948 patients; 4-14 headache days per",
      "month) or chronic migraine (596 patients; at least 15 headache days",
      "per month, more than 8 of them migraine days). Baseline MMD median 9.00",
      "(range 3.23-28.0) days/month."
    ),
    dose_range = paste(
      "Placebo (792), atogepant 10 mg QD (306), 30 mg QD (405), 60 mg QD (703),",
      "30 mg BID (251) or 60 mg BID (87), oral, for 12 weeks."
    ),
    regions = "Multinational (phase 2b/3 EM study, ADVANCE, PROGRESS, ELEVATE).",
    notes = paste(
      "Demographics from Table 1 of Schlachter 2026 (CPT PSP). 280 patients",
      "(11%, all from ELEVATE) had failed two to four classes of conventional",
      "oral preventive treatments; 1604 (63%) had received prior preventive",
      "migraine medication. Observation unit is one month of diary-recorded",
      "migraine days (normalised to 28 days) over a 12-week treatment period."
    )
  )

  ini({
    # Final exposure-response estimates: Schlachter 2026 Table 4 (compressed
    # normal distribution, first-order estimation); structure from Data S1.
    # pbo, emax and the effect covariates act on the shifted-logit scale of
    # the combined effect, e1 = 2 * expit(eff) - 1.
    ec50 <- 0.412; label("Atogepant trough concentration at half-maximal effect EC50 (nM)") # Table 4 'Atogepant Cmin EC50 (nM)' 0.412 (RSE 60.0%); Data S1 THETA(1)
    pbo <- 0.852; label("Placebo effect on the shifted-logit scale (unitless)") # Table 4 'Placebo effect' 0.852 (RSE 4.62%); Data S1 THETA(2) LGPEFF
    lt50 <- log(0.407); label("Time to half of the combined placebo and atogepant effect T50 (month)") # Table 4 'TEF50 (month)' 0.407 (RSE 3.24%); Data S1 THETA(3)
    emax <- 0.665; label("Maximal atogepant effect on the shifted-logit scale at 9 baseline migraine days (unitless)") # Table 4 'Atogepant Cmin Emax' 0.665 (RSE 8.43%); Data S1 THETA(4)
    e_migraine_days_bl_emax <- -0.0166; label("Additive change in atogepant Emax per baseline migraine day above 9 (per day/month)") # Table 4 'Baseline MMD on atogepant Emax' -0.0166 (RSE 43.8%); Data S1 THETA(5)
    e_migraine_oral_preventive_ir_eff <- -0.410; label("Additive shift of the combined effect for failed oral preventive treatments (shifted-logit units)") # Table 4 'Previously failed treatments on combined placebo & atogepant effect' -0.410 (RSE 18.5%); Data S1 THETA(6)
    e_migraine_days_bl_t50 <- log(1.02); label("Log-scale change in T50 per baseline migraine day above 9 (per day/month)") # Table 4 'Baseline MMD on TEF50' 1.02 = exp(THETA(7)), 95% CI 1.02-1.03; Data S1 'THETA(7) * (BSL - 9.00)' inside EXP()
    e_migraine_oral_preventive_ir_t50 <- log(0.628); label("Log-scale T50 multiplier for failed oral preventive treatments (unitless)") # Table 4 'Previously failed treatments on TEF50' 0.628 = exp(THETA(8)), 95% CI 0.513-0.768 symmetric about log(0.628); Data S1 THETA(8) inside EXP()
    e_race_asian_eff <- -0.326; label("Additive shift of the combined effect for Asian race (shifted-logit units)") # Table 4 'Asian race on combined placebo & atogepant effect' -0.326 (RSE 28.1%); Data S1 THETA(9)

    # Additive IIV on the shifted-logit scale (Data S1 '; ETA(1): add',
    # '; ETA(2): add'). Table 4 prints omega^2: %CV column = sqrt(exp(omega^2) - 1).
    etapbo ~ 0.793 # Table 4 'IIV on placebo effect' 0.793 (110 %CV, shrinkage 12.3%); placebo patients only
    etaemax ~ 1.67 # Table 4 'IIV on atogepant Emax' 1.67 (208 %CV, shrinkage 17.6%); atogepant patients only

    # Table 4 prints the residual VARIANCES (days^2); SD = sqrt(variance).
    addSdEpisodic <- sqrt(3.34); label("Additive residual SD for episodic migraine (days/month)") # Table 4 'Additive error episodic (days2)' 3.34 (RSE 0.835%); Data S1 SIGMA1
    addSdChronic <- sqrt(7.98); label("Additive residual SD for chronic migraine (days/month)") # Table 4 'Additive error chronic (days2)' 7.98 (RSE 1.99%); Data S1 SIGMA2
  })

  model({
    # Covariate terms, centred at the pooled median baseline of 9.00 MMD
    # (Data S1 'BSL - 9.00').
    bl_centred <- MIGRAINE_DAYS_BL - 9
    eff_cov <- e_migraine_oral_preventive_ir_eff *
      MIGRAINE_ORAL_PREVENTIVE_IR +
      e_race_asian_eff * RACE_ASIAN

    # Individual atogepant Emax (additive IIV) and hyperbolic drug effect of
    # the steady-state trough concentration.
    emax_i <- emax + e_migraine_days_bl_emax * bl_centred + etaemax
    drug_eff <- emax_i * CTROUGH / (CTROUGH + ec50)

    # Combined effect on the shifted-logit scale. Placebo patients carry the
    # placebo random effect and no drug term; atogepant patients carry the
    # drug term (with its Emax random effect) and no placebo random effect.
    pbo_i <- pbo + etapbo
    eff <- PLACEBO * pbo_i + (1 - PLACEBO) * (pbo + drug_eff) + eff_cov

    # Back-transform to a fractional shift in (-1, 1).
    e1 <- 2 * expit(eff) - 1

    # Time to half of the combined effect.
    t50 <- exp(
      lt50 +
        e_migraine_days_bl_t50 * bl_centred +
        e_migraine_oral_preventive_ir_t50 * MIGRAINE_ORAL_PREVENTIVE_IR
    )

    # Predicted monthly migraine days (Data S1 IPRED).
    migraine_days <- MIGRAINE_DAYS_BL * (1 - e1 * (1 - exp(-time * log(2) / t50)))

    addSd_ind <- addSdEpisodic *
      (1 - DIS_MIGRAINE_CHRONIC) +
      addSdChronic * DIS_MIGRAINE_CHRONIC
    migraine_days ~ add(addSd_ind)
  })
}
