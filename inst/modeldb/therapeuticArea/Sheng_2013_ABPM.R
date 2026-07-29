Sheng_2013_ABPM <- function() {
  description <- "Cyclic-fluctuation (circadian rhythm) model for 24-h ambulatory blood pressure monitoring (24-h ABPM) in Chinese patients with mild-to-moderate essential hypertension during the placebo run-in period of four antihypertensive drug clinical trials (Sheng 2013). Predicts systolic blood pressure (SBP) and diastolic blood pressure (DBP) simultaneously as the sum of (a) a rhythm-adjusted 24-h mean and (b) two cosine harmonics with periods 0.5 day (12-h harmonic) and 1 day (24-h harmonic). The two phase-shift parameters PHS1 (12-h harmonic) and PHS2 (24-h harmonic) are shared across SBP and DBP because the authors found the per-output phase-shift estimates were close enough to combine; the four amplitudes (one per output x harmonic) and the two baselines (one per output) remain output-specific. No drug input and no compartments -- a baseline-BP circadian model intended for combination with a drug PK/PD model when simulating antihypertensive trials."

  reference <- paste(
    "Sheng YC, Wang K, Xu L, Yang J, He YC, Zheng QS.",
    "A cyclic fluctuation model for 24-h ambulatory blood pressure monitoring",
    "in Chinese patients with mild to moderate hypertension.",
    "Acta Pharmacologica Sinica. 2013 Aug;34(8):1043-1051.",
    "doi:10.1038/aps.2013.45.",
    sep = " "
  )
  vignette <- "Sheng_2013_ABPM"

  # Etas pair additively with paper-symbol amplitudes (a11/a21/a12/a22) and
  # phase shifts (phs1/phs2) rather than with log-transformed parents.
  # The published variance magnitudes (Table 2: ~5-7 mmHg^2 for amplitudes,
  # ~1-11 day^2 for phase shifts) make sense only as additive variances on
  # the linear scale -- a log-normal interpretation would give CV > 1000%,
  # which is incompatible with the bootstrap CIs the paper reports. See the
  # vignette Assumptions and deviations section for the full reconciliation.
  paper_specific_etas <- c(
    "etaa11", "etaa21", "etaa12", "etaa22",
    "etaphs1", "etaphs2"
  )

  units <- list(
    time          = "day",
    dosing        = "n/a (baseline BP rhythm model with no drug input)",
    concentration = "mmHg (blood pressure; outputs SBP and DBP)"
  )

  # No covariates retained in the final model (Sheng 2013 Results, page 1046:
  # 'The covariates, including gender, age, weight and BMI, did not
  # significantly increase the goodness of fit.'). The four covariates the
  # authors screened are documented in covariatesDataExcluded so the
  # provenance is preserved without triggering convention warnings.
  covariateData <- list()

  covariatesDataExcluded <- list(
    SEXF = list(
      description        = "Female-sex indicator (1 = female, 0 = male). Screened by Sheng 2013 against all eight structural parameters (Base1, Base2, A11, A21, A12, A22, PHS1, PHS2) and not retained in the final model.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Screened in the covariate-evaluation step of Sheng 2013 (Methods 'Covariate model selection', Results 'The covariates ... did not significantly increase the goodness of fit'); not retained in the final structural model. Recorded here only to preserve the provenance of the covariate screen.",
      source_name        = "Sex (male/female)"
    ),
    AGE = list(
      description        = "Subject age. Screened by Sheng 2013 against all eight structural parameters; not retained in the final model.",
      units              = "year",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened (Methods 'Covariate model selection'); not retained. Cohort medians 49-52 years across the four trials (Table 1).",
      source_name        = "Age (year)"
    ),
    WT = list(
      description        = "Body weight. Screened by Sheng 2013 against all eight structural parameters; not retained in the final model.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened (Methods 'Covariate model selection'); not retained. Cohort medians 66-74 kg across the four trials (Table 1).",
      source_name        = "Weight (kg)"
    ),
    BMI = list(
      description        = "Body mass index. Screened by Sheng 2013 against all eight structural parameters; not retained in the final model.",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened (Methods 'Covariate model selection'); not retained. Cohort medians 24-26 kg/m^2 across the four trials (Table 1).",
      source_name        = "BMI (kg/m^2)"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 38L,
    n_studies       = 1L,
    age_range       = "35-69 years",
    age_median      = "51.6 years (mean +/- 7.83 SD; Sheng 2013 Table 1, Study 1)",
    weight_range    = "45-85 kg",
    weight_median   = "65.5 kg (mean +/- 9.27 SD)",
    sex_female_pct  = 55.3,
    race_ethnicity  = "Han Chinese (Sheng 2013 Methods; all subjects enrolled in Chinese centres)",
    disease_state   = "Mild-to-moderate essential hypertension (mean sitting SBP/DBP 140-179/90-109 mmHg). Major exclusions: significant cardiovascular, hepatic, or renal disease; type 1 diabetes or uncontrolled type 2 diabetes; pregnancy potential; current anticonvulsants or antidepressants; drug-abuse history.",
    dose_range      = "Not applicable -- the parameter estimates come from the 2-week placebo run-in period (no active antihypertensive administered).",
    regions         = "China",
    notes           = "Final-model estimates come from Study 1 (38 patients, 2110 SBP and 2110 DBP measurements at the end of the 2-week placebo run-in period; Sheng 2013 Table 2). The structural model was subsequently re-estimated on three additional Chinese antihypertensive-trial cohorts (Study 2: n=42; Study 3: n=25; Study 4: n=29) and on the pooled four-study dataset (n=134), all yielding very similar parameter estimates (Sheng 2013 Table 3). Each subject contributed 24-h ABPM with measurements every 15 min from 8 AM-10 PM and every 30 min from 10 PM-8 AM, recorded with a SunTech Medical Instruments ABP monitor. The four-study population was 56-69 patients per cohort with 50-55 mean age and similar weight / BMI ranges; the full pooled cohort is summarised in Sheng 2013 Table 1.",
    average_sbp_baseline_mmHg = "141 +/- 15.6 (Study 1 Table 1)",
    average_dbp_baseline_mmHg = "90.4 +/- 10.2 (Study 1 Table 1)",
    n_sbp_observations       = 2110L,
    n_dbp_observations       = 2110L,
    nonmem_method            = "FOCE-I (NONMEM 7.2 via Wings for NONMEM WFN720; Perl-speaks-NONMEM 3.5.3 for bootstrap)"
  )

  ini({
    # =====================================================================
    # Final cyclic-fluctuation model (Sheng 2013 Results, page 1045):
    #
    #   BP_k(t) = Base_k
    #           + A1k * cos[ 2*pi * (t - PHS1) / 0.5 ]   <-- 12-h harmonic
    #           + A2k * cos[   pi * (t - PHS2) / 0.5 ]   <-- 24-h harmonic
    #
    # with k = 1 -> SBP and k = 2 -> DBP. The denominator 0.5 in both cosine
    # arguments gives the 12-h harmonic a period of 2*pi / (2*pi/0.5) = 0.5
    # day and the 24-h harmonic a period of 2*pi / (pi/0.5) = 1 day. The
    # phase shifts PHS1 and PHS2 are shared across SBP and DBP because Sheng
    # 2013 found the per-output estimates were sufficiently close to
    # combine ('Based on a graphical inspection of the raw data, equal
    # values of phase-shift parameters (i.e., PHS_i1 = PHS_i2) for both
    # diastolic and systolic blood pressure measurements were also used',
    # Methods Structural model). All point estimates below are from Sheng
    # 2013 Table 2 (Study 1, the developmental dataset; 38 subjects, 2110
    # SBP and 2110 DBP observations).
    # =====================================================================

    # ---- Rhythm-adjusted 24-h mean (mmHg). Log-transformed so the
    # log-normal IIV expression in model() (rbase = exp(lrbase + eta))
    # matches the paper's P_i = P_tv * exp(eta_i) ('Statistical models',
    # page 1045).
    lrbase_SBP <- log(140)    ; label("Rhythm-adjusted 24-h mean SBP (mmHg)")  # Sheng 2013 Table 2, row Base1 (estimate 140, SEM 1.8%; bootstrap median 140, 95% CI 136-145)
    lrbase_DBP <- log(89.5)   ; label("Rhythm-adjusted 24-h mean DBP (mmHg)")  # Sheng 2013 Table 2, row Base2 (estimate 89.5, SEM 1.8%; bootstrap median 89.6, 95% CI 86.7-92.6)

    # ---- Cosine amplitudes (mmHg). Kept on the linear scale (no log
    # transform) because Sheng 2013 reports the IIV variances on the linear
    # additive scale -- see paper_specific_etas comment at the top of the
    # function body.
    a11 <- 7.52              ; label("Amplitude of 12-h-harmonic cosine, SBP (mmHg)")  # Sheng 2013 Table 2, row A11 (estimate 7.52, SEM 8.7%; bootstrap median 7.52, 95% CI 6.37-8.79)
    a21 <- 8.64              ; label("Amplitude of 24-h-harmonic cosine, SBP (mmHg)")  # Sheng 2013 Table 2, row A21 (estimate 8.64, SEM 9.1%; bootstrap median 8.58, 95% CI 6.89-10.19)
    a12 <- 5.61              ; label("Amplitude of 12-h-harmonic cosine, DBP (mmHg)")  # Sheng 2013 Table 2, row A12 (estimate 5.61, SEM 7.6%; bootstrap median 5.60, 95% CI 4.70-6.50)
    a22 <- 6.27              ; label("Amplitude of 24-h-harmonic cosine, DBP (mmHg)")  # Sheng 2013 Table 2, row A22 (estimate 6.27, SEM 9.3%; bootstrap median 6.32, 95% CI 5.21-7.50)

    # ---- Shared phase shifts (days). Linear (un-transformed) because
    # phase shifts can be of either sign; PHS1 (-0.652) is the 12-h-harmonic
    # phase, PHS2 (3.61) the 24-h-harmonic. Both are time offsets in the
    # paper's time-unit-is-day convention; cosines are 2*pi-periodic so
    # equivalent phase shifts differ by multiples of the corresponding
    # period (0.5 day for PHS1, 1 day for PHS2).
    phs1 <- -0.652           ; label("Phase shift of 12-h-harmonic cosine, shared SBP/DBP (day)")  # Sheng 2013 Table 2, row PHS1 (estimate -0.652, SEM 1.3%; bootstrap median -0.651, 95% CI -0.668 to 0.862)
    phs2 <-  3.61            ; label("Phase shift of 24-h-harmonic cosine, shared SBP/DBP (day)")  # Sheng 2013 Table 2, row PHS2 (estimate 3.61, SEM 0.8%; bootstrap median 3.61, 95% CI 3.50-3.69)

    # ---- IIV. Sheng 2013 states 'log-normally distributed' IIV with
    # P_i = P_tv * exp(eta_i), eta_i ~ N(0, omega^2) ('Statistical models',
    # page 1045). The published Table 2 omega^2 values are interpretable
    # as log-scale variances only for the two baselines (0.079, 0.12 ->
    # CV approx 28-35%, plausible). For amplitudes (5.63 mmHg^2 etc.) and
    # phase shifts (10.9 day^2, 1.24 day^2) a log-normal interpretation
    # gives CV > 100% which is incompatible with the narrow bootstrap CIs
    # the same table reports. The only interpretation consistent with both
    # the published point estimates and the bootstrap CIs is that the
    # amplitude and phase-shift IIV is additive on the linear scale, i.e.
    # P_i = P_tv + eta_i with omega^2 in the units of the parameter itself
    # (mmHg^2 or day^2). This is the encoding used in model() below.
    etalrbase_SBP ~ 0.079    # Sheng 2013 Table 2, row omega_Base1 (estimate 0.079, SEM 20.6%; bootstrap median 0.079, 95% CI 0.061-0.094)
    etalrbase_DBP ~ 0.12     # Sheng 2013 Table 2, row omega_Base2 (estimate 0.12, SEM 21.3%; bootstrap median 0.12, 95% CI 0.094-0.146)
    etaa11        ~ 5.63     # Sheng 2013 Table 2, row omega_A11 (estimate 5.63 mmHg^2, SEM 31.8%; bootstrap median 5.47, 95% CI 3.42-7.14); additive IIV on linear scale -- see narrative above
    etaa21        ~ 4.73     # Sheng 2013 Table 2, row omega_A21 (estimate 4.73 mmHg^2, SEM 33.0%; bootstrap median 4.74, 95% CI 3.12-6.91); additive IIV on linear scale
    etaa12        ~ 6.58     # Sheng 2013 Table 2, row omega_A12 (estimate 6.58 mmHg^2, SEM 37.6%; bootstrap median 6.39, 95% CI 3.20-8.70); additive IIV on linear scale
    etaa22        ~ 6.71     # Sheng 2013 Table 2, row omega_A22 (estimate 6.71 mmHg^2, SEM 35.6%; bootstrap median 6.36, 95% CI 4.05-8.83); additive IIV on linear scale
    etaphs1       ~ 10.9     # Sheng 2013 Table 2, row omega_PHS1 (estimate 10.9 day^2, SEM 27.5%; bootstrap median 10.05, 95% CI 7.0-13.3); additive IIV on linear scale (cosine is 2*pi-periodic so any draw is valid)
    etaphs2       ~ 1.24     # Sheng 2013 Table 2, row omega_PHS2 (estimate 1.24 day^2, SEM 37.6%; bootstrap median 1.17, 95% CI 0.59-1.6); additive IIV on linear scale

    # ---- Residual error. Sheng 2013 ('Statistical models') describes an
    # additive normal residual Y_ij = IPRE_ij + epsilon, epsilon ~ N(0,
    # sigma^2). Table 2 reports sigma values with the column unit '(mmHg)';
    # the reported values 12.85 and 9.11 are the additive residual SDs in
    # mmHg (the row unit label is the SD's natural unit; the bootstrap CI
    # 12.1-13.7 mmHg for SBP brackets the same scale). Encoded here as the
    # per-output additive SD via `~ add(addSd_<output>)` in model().
    addSd_SBP <- 12.85       ; label("Intraindividual residual SD for SBP (mmHg)")  # Sheng 2013 Table 2, row sigma_SBP (estimate 12.85, SEM 6.5%; bootstrap median 12.9, 95% CI 12.1-13.7)
    addSd_DBP <-  9.11       ; label("Intraindividual residual SD for DBP (mmHg)")  # Sheng 2013 Table 2, row sigma_DBP (estimate 9.11, SEM 7.9%; bootstrap median 9.10, 95% CI 8.46-9.79)
  })

  model({
    # 1. Individual structural parameters.
    # Baselines: log-normal IIV (multiplicative; paper's stated form).
    rbase_SBP <- exp(lrbase_SBP + etalrbase_SBP)
    rbase_DBP <- exp(lrbase_DBP + etalrbase_DBP)
    # Amplitudes: additive IIV on linear scale (see ini() comment).
    a11i <- a11 + etaa11
    a21i <- a21 + etaa21
    a12i <- a12 + etaa12
    a22i <- a22 + etaa22
    # Phase shifts: additive IIV on linear scale (cosine is periodic).
    phs1i <- phs1 + etaphs1
    phs2i <- phs2 + etaphs2

    # 2. Cosine arguments. The denominator 0.5 in the paper's published
    # equation makes the first cosine the 12-h harmonic (period 0.5 day)
    # and the second the 24-h harmonic (period 1 day). t carries the
    # paper's time-unit-is-day convention; t = 0 is treated as midnight
    # (the 24-h-harmonic peak at t = PHS2 mod 1 ~ 0.61 day ~ 2:30 PM matches
    # the typical daytime BP peak, confirming the midnight origin).
    arg_12h <- 2 * pi * (t - phs1i) / 0.5
    arg_24h <-     pi * (t - phs2i) / 0.5

    # 3. Predictions (algebraic; no ODE / no compartments because the model
    # is a baseline-BP circadian rhythm without dosing).
    SBP <- rbase_SBP + a11i * cos(arg_12h) + a21i * cos(arg_24h)
    DBP <- rbase_DBP + a12i * cos(arg_12h) + a22i * cos(arg_24h)

    # 4. Observation and error.
    SBP ~ add(addSd_SBP)
    DBP ~ add(addSd_DBP)
  })
}
