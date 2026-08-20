Shigetome_2025_paroxetine_madrs <- function() {
  description <- "Population PK-PD model for the time course of depression-severity improvement under paroxetine in 50 Japanese adults with major depressive disorder, fitted to Montgomery-Asberg Depression Rating Scale (MADRS) scores at baseline and weeks 1, 2, 4 and 6. The endpoint madrsenh is the paper's enhancement rate, the percentage reduction in MADRS from baseline, described by an Emax model in treatment duration: madrsenh = Emax * Time / (ET50 + Time). Cumulative paroxetine exposure over the first week of treatment (AUC_PAROX, predicted from the companion popPK model Shigetome_2025_paroxetine) enters Emax through the power term -(AUC_PAROX / 2610.43)^-8.56, so higher first-week exposure raises the attainable improvement; the MADRS score at week 1 (SCORE_MADRS) enters ET50 through +(SCORE_MADRS / 28.64)^2.38, so less severe week-1 depression brings the improvement forward. There is no PK compartment and no dosing event in this model: drug exposure is carried entirely by the AUC_PAROX covariate, and the independent variable is treatment duration in weeks. Proportional inter-individual variability on both Emax and ET50, additive residual error."

  reference <- paste(
    "Shigetome K, Egashira T, Tomita T, Higa N, Iwashita K, Morita K, Nishimura M,",
    "Kaneko T, Maeda H, Yamada KD, Kajiwara-Morita A, Oniki K, Yasui-Furukori N,",
    "Saruwatari J.",
    "Effect of Cumulative Exposure on the Efficacy of Paroxetine:",
    "A Population Pharmacokinetic-Pharmacodynamic and Machine Learning Analyses.",
    "CPT Pharmacometrics Syst Pharmacol. 2025 Jun;14(6):1119-1127.",
    "doi:10.1002/psp4.70032.",
    "First-week exposure AUC_PAROX is predicted by the companion popPK model;",
    "see modellib('Shigetome_2025_paroxetine').",
    sep = " "
  )
  vignette <- "Shigetome_2025_paroxetine"
  units <- list(
    time          = "week (treatment duration since paroxetine initiation; the Emax model is a function of time, not of a concurrent concentration)",
    dosing        = "not applicable (no dose events; cumulative paroxetine exposure enters through the AUC_PAROX covariate)",
    concentration = "percent / patient (output madrsenh is the enhancement rate, the percentage reduction in MADRS score from baseline, NOT a drug concentration; the slash is only to satisfy checkModelConventions unit parsing)"
  )

  covariateData <- list(
    AUC_PAROX = list(
      description        = "Cumulative area under the plasma paroxetine concentration-time curve from treatment initiation to the end of the first week of treatment (AUC 0-1 week)",
      units              = "ng*h/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per patient. Not measured directly: predicted for each patient from the companion popPK model (Shigetome_2025_paroxetine) using that patient's dosing history over week 1, then carried into this model as a covariate. Enters Emax as the subtracted power term (AUC_PAROX / 2610.43)^e_auc_parox_emax with the normalizing constant 2610.43 ng*h/mL, which is the population value in this 50-patient cohort (observed means 2764.9 in the remission group and 2472.2 in the non-remission group, Table S5). Because the exponent is negative, the subtracted term shrinks as exposure rises, so Emax increases with first-week exposure. The term is steep: it reaches Emax = 0 near AUC_PAROX = 1541 ng*h/mL and is essentially saturated above the normalizing value, so simulations outside the observed 1642-5814 ng*h/mL range extrapolate badly. Among the eight paroxetine exposure indices screened (measured trough, popPK-predicted trough, Cmax, single-dose AUC, and cumulative AUC to weeks 1, 2, 4 and 6), only AUC 0-1 week was retained (Table S3).",
      source_name        = "AUC_W1"
    ),
    SCORE_MADRS = list(
      description        = "Montgomery-Asberg Depression Rating Scale total score, ten items each rated 0-6 by a physician at interview",
      units              = "(SCORE_MADRS units, 0-60 score)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "This model uses the week-1 value (the MADRS score one week after starting paroxetine, mean 28.6, SD 11.3, range 4-56; Shigetome 2025 Table 2). Time-fixed once week 1 has passed. Enters ET50 as the added power term (SCORE_MADRS / 28.64)^e_score_madrs_et50 with the normalizing constant 28.64 (the cohort value), so a lower week-1 score shortens ET50 and brings the improvement forward. The baseline MADRS score is NOT this covariate: it defines the endpoint instead, since the enhancement rate is the percentage reduction from baseline. Baseline MADRS was screened on both Emax and ET50 and was not retained (Table S3).",
      source_name        = "MADRS_W1"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 50L,
    n_studies      = 1L,
    age_range      = "20-70 years",
    age_mean       = "46.6 years (SD 13.7)",
    weight_range   = "40-85 kg",
    weight_mean    = "56.3 kg (SD 10.4)",
    sex_female_pct = 62.0,
    race_ethnicity = c(Asian = 100),
    disease_state  = "major depressive disorder (DSM-IV), antidepressant-naive at paroxetine initiation; baseline MADRS 39.6 (SD 9.7, range 12-53)",
    dose_range     = "10-40 mg/day orally (mean 32.2 mg/day, SD 11.0); initial dose 10-20 mg/day then weekly increases of 10 mg/day to the maximum tolerated dose",
    regions        = "Japan (Hirosaki University Hospital and Dokkyo Medical University School of Medicine)",
    timepoints     = "MADRS assessed at baseline and 1, 2, 4 and 6 weeks after starting paroxetine; observed enhancement rates 25.9% (SD 30.3), 45.7% (SD 27.4), 57.5% (SD 28.3) and 64.2% (SD 30.6) at weeks 1, 2, 4 and 6",
    notes          = "The 50 patients are the subset of the 179-patient PK cohort (Shigetome_2025_paroxetine) who had serial MADRS assessments and a pre-treatment Temperament and Character Inventory. Baseline demographics are Shigetome 2025 Table 2. Remission (MADRS < 10 at week 6) was reached by 26 of the 50. All seven TCI personality scales and daylight hours were screened on Emax and ET50 and none was retained (Table S3)."
  )

  covariatesDataExcluded <- list(
    DAYLIGHT_H = list(
      description = "Cumulative daylight hours over the treatment period (mean 120.6 h, SD 43.4, range 54.3-173.2; Shigetome 2025 Table 2)",
      units       = "h",
      type        = "continuous",
      notes       = "Screened on both Emax and ET50 and not retained (Shigetome 2025 Table S3: univariate dOFV 0.303 on Emax and -0.193 on ET50, neither reaching the p < 0.1 forward-selection threshold). No point estimate is reported, so no covariate effect can be encoded."
    ),
    SCORE_TCI = list(
      description = "Temperament and Character Inventory scores: novelty seeking, harm avoidance, reward dependence, persistence, self-directedness, cooperativeness and self-transcendence (Shigetome 2025 Table 2)",
      units       = "(TCI subscale points)",
      type        = "continuous",
      notes       = "All seven subscales were screened on both Emax and ET50 and none was retained (Shigetome 2025 Table S3). The paper notes that week-1 MADRS is itself positively associated with harm avoidance and negatively with self-directedness, which it offers as the reason the personality effects reported in its earlier work do not appear here. Reward dependence was, separately, one of four features selected by the machine-learning arm of the paper; that analysis is a classifier, not part of this popPK/PD model."
    )
  )

  ini({
    # ================================================================
    # Final popPK/PD model (Shigetome 2025 Results section 3.3 printed
    # equations, estimates in Table S4 of the Supporting Information;
    # NONMEM control stream in Text S3):
    #
    #   EFF  = Emax * Time / (ET50 + Time)
    #   Emax = 90.9  - (AUC_0-1week / 2610.43)^-8.56
    #   ET50 = 0.984 + (MADRS_W1    / 28.64)^2.38
    #
    # NONMEM $PRED (Text S3), which also fixes the variability structure:
    #   EMAX = (THETA(1) - (AUC_W1  / 2610.43) ** THETA(4)) * (1 + ETA(1))
    #   EC50 = (THETA(2) + (MADRS_W1/   28.64) ** THETA(3)) * (1 + ETA(2))
    #   Y    = EMAX * WEEK / (EC50 + WEEK) + EPS(1)
    # (the stream's local name EC50 is the paper's ET50, a treatment
    # duration in weeks, not a concentration.)
    #
    # Emax and ET50 are NOT log-transformed here, because the covariate
    # terms are subtracted from / added to them on the natural scale;
    # a log parameterisation cannot express that structure. Emax can
    # legitimately go negative at low first-week exposure.
    # ================================================================

    emax <- 90.9
    label("Maximum enhancement rate in depression severity at the reference first-week exposure (percentage reduction in MADRS from baseline)")  # Shigetome 2025 Table S4, Emax = 90.9 (RSE 9.3%, bootstrap median 91.5, 95% CI 74.7-111)

    e_auc_parox_emax <- -8.56
    label("Power exponent on (AUC_PAROX / 2610.43) in the term subtracted from Emax (unitless)")  # Shigetome 2025 Table S4, effect of AUC 0-1 week on Emax = -8.56 (RSE 3.36%, bootstrap 95% CI -12.7 to -2.20)

    et50 <- 0.984
    label("Treatment duration reaching half of Emax at zero week-1 MADRS, the intercept of ET50 (week)")  # Shigetome 2025 Table S4, ET50 = 0.984 (RSE 45.5%, bootstrap median 0.931, 95% CI 0.101-1.98)

    e_score_madrs_et50 <- 2.38
    label("Power exponent on (SCORE_MADRS / 28.64) in the term added to ET50 (unitless)")  # Shigetome 2025 Table S4, effect of MADRS score at week 1 on ET50 = 2.38 (RSE 43.3%, bootstrap 95% CI 0.666-6.66)

    # ================================================================
    # Inter-individual variability: PROPORTIONAL on both Emax and ET50
    # (Results 3.3, "the proportional error model best accounted for
    # the inter-individual errors of Emax and ET50"; NONMEM Text S3
    # writes "* (1 + ETA(n))"). The ini() values are the NONMEM OMEGA
    # variances of the etas, reported in a column headed "%CV" in
    # Table S4. They are variances, not CVs: the same table reports
    # the additive residual as 190 for an endpoint measured in percent,
    # which is only interpretable as a variance (SD 13.8 percentage
    # points), and only the variance reading of the popPK OMEGA
    # reproduces the observed concentration spread in the companion
    # PK model. See the vignette Errata section for the full argument.
    #
    # Note that proportional IIV on ET50 with a variance of 0.486
    # (SD 0.697) puts about 8% of simulated subjects at a negative
    # ET50; this is the authors' structure reproduced faithfully, and
    # simulations should screen for it.
    # ================================================================
    etaemax ~ 0.150  # Shigetome 2025 Table S4, inter-individual variability on Emax (RSE 34.1%, shrinkage 16.0%, bootstrap 95% CI 0.042-0.258)
    etaet50 ~ 0.486  # Shigetome 2025 Table S4, inter-individual variability on ET50 (RSE 55.1%, shrinkage 32.0%, bootstrap 95% CI 0.0315-1.46)

    # ================================================================
    # Residual error: additive on the enhancement rate (Results 3.3,
    # "the residual error was best with an additional error model";
    # NONMEM Text S3 "Y = ... + EPS(1)"). SIGMA = 190 is a variance in
    # percent-squared, so the additive SD is sqrt(190) = 13.78
    # percentage points.
    # ================================================================
    addSd <- 13.78
    label("Additive residual error SD on the enhancement rate (percentage points)")  # Shigetome 2025 Table S4, additive error 190 (variance; RSE 23.5%, shrinkage 16.3%, bootstrap 95% CI 112-316); sqrt(190) = 13.784

    # Note: an inter-individual-variability shrinkage of 32% on ET50 is
    # flagged by the authors as slightly high but acceptable
    # (Results 3.3).
  })

  model({
    # Individual maximum enhancement rate: the reference Emax reduced by
    # a steep power term in cumulative first-week paroxetine exposure,
    # then scaled by the proportional eta.
    emax_ind <- (emax - (AUC_PAROX / 2610.43)^e_auc_parox_emax) * (1 + etaemax)

    # Individual time to half-maximal enhancement: the intercept plus a
    # power term in the week-1 MADRS score, scaled by the proportional
    # eta. Units are weeks, matching the model time unit.
    et50_ind <- (et50 + (SCORE_MADRS / 28.64)^e_score_madrs_et50) * (1 + etaet50)

    # Enhancement rate in depression severity, i.e. the percentage
    # reduction in MADRS from baseline. At time 0 the fraction is 0 so
    # madrsenh = 0; as time grows large madrsenh approaches emax_ind.
    madrsenh <- emax_ind * time / (et50_ind + time)
    madrsenh ~ add(addSd)
  })
}
