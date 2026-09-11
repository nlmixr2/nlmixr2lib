Zou_2026_pembrolizumab_qol_mbma <- function() {
  description <- "MBMA. Longitudinal model-based meta-analysis of patient-reported quality of life (EORTC QLQ-C30 Global Health Status / Quality of Life, GHS/QoL) under pembrolizumab versus control, fit to study-arm-level aggregate data digitised from 20 published oncology clinical trials (19 unique trials; 36 arms; 410 arm-timepoint observations) across melanoma, NSCLC, TNBC, HNSCC, urothelial, colorectal and endometrial cancer. The QoL trajectory is a semi-mechanistic sum of an asymptotic toxicity decline and a linear long-term improvement, written on the 0-1 fraction scale as Q(t) = E0 - Emax * (1 - exp(-Kp * t)) + SLP * t, where E0 (baseline) and Emax (maximal toxicity) are logit-transformed to keep the predicted score inside 0-100 and SLP is log-transformed. Treatment is the only retained covariate and acts on BOTH Emax (logit scale, -0.758) and SLP (log scale, -1.4) with pembrolizumab as the reference: the control arm has both less early toxicity and a markedly slower long-term QoL improvement. Variability is two-level and MBMA-specific: between-study variability (BSV) on logit-E0 (SD 0.667) and logit-Emax (SD 0.783), and between-treatment-arm variability (BTAV) on logit-E0 (SD 0.867) and log-SLP (SD 0.871). The BTAV term on E0 and the residual error are BOTH divided by sqrt(N_ARM), the number of participants contributing to the arm-level mean, so arm size is a required model input rather than a downstream weighting; supply it in the N_ARM column. Suitable simulation scope is study-arm-mean QoL trajectories over roughly 0-100 weeks, NOT individual-patient QoL scores. Parameter values are Zou 2026 Table S1 (final model estimates) and the model structure is the Monolix 2024R1 control stream reproduced verbatim in Supplementary Codes section A."

  reference <- paste(
    "Zou Y, Sun Y, Ravva S, Wagner LI, Zhou J.",
    "A Model-Based Meta-Analysis of Pembrolizumab Effects on Patient-Reported",
    "Quality of Life: Advancing Patient-Centered Oncology Drug Development.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15:e70106.",
    "doi:10.1002/psp4.70106.",
    sep = " "
  )
  vignette <- "Zou_2026_pembrolizumab_qol"

  units <- list(
    time          = "week",
    dosing        = "none",
    concentration = "score (arm-mean EORTC QLQ-C30 GHS/QoL on the 0-100 scale, exposed as qol_ghs; the model observation logit_qol_ghs is that score on the logit-of-fraction scale and is NOT a drug concentration)"
  )

  covariateData <- list(
    TRT = list(
      description        = "Study-arm treatment indicator: 0 = pembrolizumab arm (pembrolizumab monotherapy or pembrolizumab in combination with another agent), 1 = control arm (the comparator regimen of the original trial -- chemotherapy, targeted therapy, placebo plus standard of care, or placebo alone, depending on trial design).",
      units              = "(categorical)",
      type               = "categorical",
      reference_category = "0 (pembrolizumab)",
      notes              = "MBMA study-arm-level indicator (a property of the trial arm, not of an individual patient). Zou 2026 Section 2.2 defines exactly two arm groups and pools every comparator regimen into the single 'control' level, so this column cannot distinguish chemotherapy from placebo control. Reference is pembrolizumab, not placebo -- this is the reverse of the usual popPK convention and the sign of both retained covariate effects follows from it (Table S1 row labels read 'Effect of control arm on ..., pembrolizumab as reference'). The paper's Equation 5 writes the effect as P_ik = theta_P * exp(Cov_trt) with Cov_trt = 0 for pembrolizumab and theta_P_control for control; the Monolix [INDIVIDUAL] block applies that shift on the LOGIT scale for the logit-normal Emax and on the LOG scale for the log-normal SLP, which is how it is encoded here.",
      source_name        = "TRT (Monolix [COVARIATE] block, categories 'Pembrolizumab' and 'Placebo'; Zou 2026 Supplementary Codes section B)"
    ),
    N_ARM = list(
      description        = "Number of participants contributing to the study-arm-level QoL mean at that observation; the arm sample size used as the meta-analytic weight.",
      units              = "participants",
      type               = "count",
      reference_category = NULL,
      notes              = "MBMA weighting regressor, supplied per observation row rather than estimated. Zou 2026 uses it in TWO places, which is why it is a covariate here rather than a downstream scaling applied after the solve (the pattern used by Mercier_2014_tramadol_tapentadol_mbma and Chen_2025_methotrexate_*_mbma, where only the residual is weighted). (1) The between-treatment-arm random effect on baseline QoL is divided by sqrt(N_ARM) -- Supplementary Codes section A, 'tE0RE = tE0 + etaBSVE0 + etaBTAVE0/sqrt(NOC)' -- so a large arm's mean baseline sits closer to its study's baseline, as an arm mean's standard error should. (2) The residual SD is divided by sqrt(N_ARM), equivalent to the paper's Equation 3 variance sigma^2 / N_ijk. Baseline arm sizes in the fitted database span 26 to 1098 participants (Zou 2026 Table 1). Must be strictly positive; the model divides by sqrt(N_ARM) twice.",
      source_name        = "NOC (Monolix regressor; Zou 2026 Supplementary Codes sections A and B) / N_ijk (Zou 2026 Equations 2 and 3)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Study-arm mean patient age.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened by the full fixed-effects covariate model of Zou 2026 Section 2.5 on E0, Emax and SLP using the centred exponential form P_i = theta_P * exp(theta_Age * AGE / Age_Ref) (Equation 6), but NOT retained: Section 3.2 states 'No additional significant covariates were identified during model development' beyond treatment. No point estimate and no Age_Ref value is published, so the effect cannot be encoded."
    ),
    SEXF_PCT = list(
      description        = "Study-arm percentage of enrolled participants who are female.",
      units              = "%",
      type               = "continuous",
      reference_category = NULL,
      notes              = "The paper collects and screens the arm-level proportion of MALE patients (Monolix column MaleP); SEXF_PCT = 100 - MaleP. Screened on E0, Emax and SLP per Zou 2026 Section 2.5 and not retained. Recorded on the female-percentage orientation so it matches the individual-level SEXF canonical (1 = female) and the RACE_ASIAN_PCT / PS_ECOG_0_PCT arm-level percentage family; the transformation from the paper's column is stated here so the provenance is not lost."
    ),
    PS_ECOG_0_PCT = list(
      description        = "Study-arm percentage of enrolled participants with an ECOG performance status of 0 at baseline.",
      units              = "%",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Monolix column ECOG0P. Screened on E0, Emax and SLP per Zou 2026 Section 2.5 and not retained. Uses the registered arm-level percentage canonical PS_ECOG_0_PCT rather than a per-subject ECOG indicator, because this MBMA carries the arm's ECOG-0 fraction."
    ),
    DIS_STAGE4_PCT = list(
      description        = "Study-arm percentage of enrolled participants with stage IV disease at baseline.",
      units              = "%",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Monolix column DSIVP. Screened on E0, Emax and SLP per Zou 2026 Section 2.5 and not retained. The Monolix [FILEINFO] header also carries DSIP, DSIIP and DSIIIP (stage I, II and III percentages) but only DSIVP is declared as a covariate in [CONTENT], so only the stage IV fraction entered the screen. Follows the DIS_CHD_PERCENT / TUMTP_SQUAM_PCT arm-level prevalence-percentage family; documentation only, so not registered in inst/references/covariate-columns.md."
    ),
    TUMTP = list(
      description        = "Study-arm tumour type, a seven-level categorical: colorectal, endometrial, HNSCC, melanoma, NSCLC, TNBC, urothelial.",
      units              = "(categorical)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Monolix column DIS, declared with exactly those seven categories in the [COVARIATE] block. Screened on E0, Emax and SLP per Zou 2026 Section 2.5 and not retained. The paper's Limitations paragraph is explicit that this is a power problem rather than evidence of no effect: 'the dataset lacked sufficient power to compare pembrolizumab's QoL benefits across different tumor types.' Recorded as the single source categorical rather than as seven TUMTP_<type> indicator columns because no per-type effect was estimated; a future extraction that does estimate per-type effects should use the registered TUMTP_MEL / TUMTP_NSCLC / TUMTP_CRC / TUMTP_BLADDER indicators."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 11326L,
    n_studies      = 20L,
    n_arms         = 36L,
    n_observations = 410L,
    age_range      = "collected per arm but not reported in the publication (screened as the AGE covariate and not retained)",
    weight_range   = "not reported",
    sex_female_pct = NA_real_,
    race_ethnicity = "not reported",
    disease_state  = "adults with advanced / metastatic solid tumours enrolled in trials of pembrolizumab: melanoma (4 study entries), NSCLC (8), TNBC (3), urothelial carcinoma (2), HNSCC (1), colorectal cancer (1) and endometrial cancer (1)",
    dose_range     = "not modelled -- this MBMA has no exposure term; pembrolizumab arms are pooled across monotherapy and combination regimens and no dose or concentration enters the model",
    regimens       = "pembrolizumab monotherapy or pembrolizumab plus another agent (pooled into one 'pembrolizumab' level) versus the trial's comparator arm (chemotherapy, targeted therapy, placebo plus standard of care, or placebo alone; pooled into one 'control' level)",
    timepoints     = "arm-mean EORTC QLQ-C30 GHS/QoL reported at 3 to 24 timepoints per arm (Zou 2026 Table 1); at least three measurements per arm were required for inclusion",
    regions        = "not reported",
    notes          = "MBMA at the study-arm level: each modelled data point is the arm-mean QoL score in one trial arm at one timepoint, weighted by the arm sample size N_ARM. 228 of the 410 observations are from pembrolizumab arms and 182 from control arms (Zou 2026 Section 3.1). n_studies = 20 counts Table 1 rows; there are 19 unique trials because KEYNOTE-054 (NCT02362594) contributes two rows -- Bottomley 2021 and Buhrer 2024 report the same trial (Table 1 footnote a), with identical baseline arm sizes of 514 pembrolizumab and 505 control. n_arms = 36 counts the 19 rows with non-zero pembrolizumab data plus the 17 rows with non-zero control data. n_subjects = 11326 is a DERIVED lower bound, not a published figure -- the paper reports no participant total. It is the sum of the Table 1 baseline sample sizes (7248 pembrolizumab + 5097 control = 12345) minus the 1019 participants of the second KEYNOTE-054 row. Table 1 footnote a states that the two KEYNOTE-054 entries 'were published with different population inclusion criteria and follow-up durations; therefore, both datasets were included', so the two rows overlap heavily but are not the identical cohort and the exact unique-participant count is indeterminate from what is published. Of the 16 study entries reporting a between-group comparison, 8 reported improved QoL with pembrolizumab and 8 reported no significant or clinically meaningful difference; the model recovers a treatment effect in the second subset too (Zou 2026 Figure 4C). Sources of data were published figures and tables digitised with WebPlotDigitizer. The model is intended for simulating arm-mean QoL trajectories and is NOT suitable for individual-subject simulation. Zou 2026 also reports Wilcoxon comparisons of the empirical-Bayes SLP and Emax between arms (Figure 4B/4C); those are post-hoc statistics on the fitted parameters, not additional model components."
  )

  ini({
    # ================================================================
    # Structural model (Zou 2026 Equation 1, and Supplementary Codes
    # section A which is the fitted Monolix 2024R1 control stream):
    #
    #   Q(t) = E0 - Emax * (1 - exp(-Kp * t)) + SLP * t
    #
    # Q is on the 0-1 FRACTION scale; the reported QoL score is
    # 100 * Q. E0 and Emax are logit-normal and SLP is log-normal
    # (Monolix [INDIVIDUAL] DEFINITION block), so the ini() entries
    # below hold logit(E0), logit(Emax) and log(SLP).
    #
    # All values are Zou 2026 Table S1, "Value" column (final model
    # parameter estimates; the Bootstrap median and 10th/90th
    # percentile columns are quoted in the comments for context).
    # ================================================================

    logite0 <- log(0.656 / (1 - 0.656))
    label("Baseline arm-mean QoL score on the logit-of-fraction scale; back-transforms to a GHS/QoL score of 65.6 on the 0-100 scale (unitless)")  # Zou 2026 Table S1 'Baseline QoL score (E0)' = 0.656 (RSE 5.16%; bootstrap median 0.656, 10th/90th 0.641-0.67). Value is the 0-1 fraction E0FE_pop; logit applied here because Monolix declares E0FE as distribution=logitNormal.

    logitemax <- log(0.0268 / (1 - 0.0268))
    label("Maximal treatment-related toxicity decrement in QoL, on the logit-of-fraction scale; back-transforms to 0.0268, i.e. 2.68 points on the 0-100 GHS/QoL scale (unitless)")  # Zou 2026 Table S1 'Maximal toxicity reducing QoL (Emax)' = 0.0268 (RSE 36.2%; bootstrap median 0.023, 10th/90th 0.011-0.04). Value is the 0-1 fraction EmaxFE_pop; logit applied here because Monolix declares EmaxFE as distribution=logitNormal.

    lkel <- log(0.0705)
    label("First-order rate constant for the onset of the toxicity decline; back-transforms to 0.0705 /week, so the toxicity reaches half its maximum at ln(2)/Kp = 9.8 weeks (log(1/week))")  # Zou 2026 Table S1 'Toxicity offset rate (Kp, 1/week)' = 0.0705 (RSE 30.4%; bootstrap median 0.069, 10th/90th 0.036-0.526). Named lkel rather than lkp to follow the Mercier_2014_tramadol_tapentadol_mbma precedent for the rate constant of an exponential approach-to-plateau in an MBMA time course, and to avoid colliding with the lkp_<tissue> partition-coefficient family.

    lslp <- log(9.98e-4)
    label("Long-term linear QoL improvement rate on the 0-1 fraction scale; back-transforms to 9.98e-4 /week, i.e. 0.0998 GHS/QoL points per week (log(1/week))")  # Zou 2026 Table S1 'QoL improvement rate (SLP, 1/week)' = 9.98e-4 (RSE 24.1%; bootstrap 10th/90th 5.44e-4-1.498e-3). Table S1's Value and percentile columns are on the 1/week scale reached AFTER the control stream's 'SLP2 = SLP*0.001' rescaling, while its Bootstrap-median cell (0.94) is the unscaled Monolix SLP_pop; 0.94 * 0.001 = 9.4e-4 falls inside the quoted 5.44e-4-1.498e-3 interval, which confirms the reading. This file works directly on the 1/week scale and drops the 0.001 bookkeeping factor -- a constant factor shifts log(SLP) without changing its log-scale SD or the log-scale covariate effect, so eta_arm_slp and e_trt_slp are unaffected.

    # ================================================================
    # Covariate model (Zou 2026 Equation 5 and Section 3.2).
    # Treatment is the ONLY retained covariate and acts on Emax and
    # SLP. Pembrolizumab is the reference (TRT = 0); the estimates
    # below are the shift applied when TRT = 1 (control arm).
    # ================================================================

    e_trt_emax <- -0.758
    label("Effect of the control arm on maximal toxicity Emax relative to pembrolizumab, additive on the logit scale (unitless)")  # Zou 2026 Table S1 'Effect of control arm on Emax, pembrolizumab as reference (theta_Emax_PBO)' = -0.758 (RSE 32.6%; bootstrap median -0.642, 10th/90th -1.399 to -0.006). Applied on the LOGIT scale per the Monolix [INDIVIDUAL] block, where EmaxFE is distribution=logitNormal with coefficient={0, beta_EmaxFE_TRT_Placebo}. Negative: control arms carry less early toxicity, consistent with Section 3.3 (Emax significantly greater under pembrolizumab, p = 0.005).

    e_trt_slp <- -1.4
    label("Effect of the control arm on the QoL improvement rate SLP relative to pembrolizumab, additive on the log scale (unitless)")  # Zou 2026 Table S1 'Effect of control arm on SLP, pembrolizumab as reference (theta_SLP_PBO)' = -1.4 (RSE 32.2%; bootstrap median -1.371, 10th/90th -2.076 to -0.768). Applied on the LOG scale per the Monolix [INDIVIDUAL] block, where SLP is distribution=logNormal with coefficient={0, beta_SLP_TRT_Placebo}. exp(-1.4) = 0.247, so the control arm improves at about a quarter of the pembrolizumab rate -- the paper's headline long-term QoL benefit (Section 3.3, p < 0.0001).

    # ================================================================
    # Two-level MBMA random effects (Zou 2026 Section 2.4 and
    # Section 3.2). BSV is between STUDY; BTAV is between TREATMENT
    # ARM within study. Table S1 reports these as STANDARD
    # DEVIATIONS under the heading "Inter-individual variability
    # (standard deviation)", so the ini() values below are their
    # squares. Table S1 additionally lists BSVEmax, BSVE0 and BTAE0
    # as "0 Fixed"; those rows are the Monolix etaBSVEmax_pop /
    # etaBSVE0_pop / etaBTAVE0_pop TYPICAL VALUES held at 0, i.e.
    # the ordinary zero-mean constraint on a random effect, not
    # separate estimable parameters.
    # ================================================================
    eta_study_e0    ~ 0.444889  # Zou 2026 Table S1 'Standard deviation of BSVE0 (Omega_BSVE0)' = 0.667 (RSE 16.1%, shrinkage 29.9%); variance = 0.667^2 = 0.444889. Between-STUDY, on logit(E0).
    eta_study_emax  ~ 0.613089  # Zou 2026 Table S1 'Standard deviation of BSVEmax (Omega_BSVEmax)' = 0.783 (RSE 31.7%, shrinkage 17.3%); variance = 0.783^2 = 0.613089. Between-STUDY, on logit(Emax).
    eta_arm_e0      ~ 0.751689  # Zou 2026 Table S1 'Standard deviation of BTAE0 (Omega_BTAE0)' = 0.867 (RSE 18.4%, shrinkage 20%); variance = 0.867^2 = 0.751689. Between-TREATMENT-ARM, on logit(E0), and divided by sqrt(N_ARM) in model() per Supplementary Codes section A.
    eta_arm_slp     ~ 0.758641  # Zou 2026 Table S1 'Standard deviation of SLP (Omega_SLP)' = 0.871 (RSE 16.1%, shrinkage 25.8%); variance = 0.871^2 = 0.758641. Between-TREATMENT-ARM (Monolix varlevel=id*occ with sd=gamma_SLP), on log(SLP). Section 3.2: 'BTAV was added to the QoL improvement rate (SLP) without logit transformation.'

    # ================================================================
    # Residual error. Zou 2026 Equation 3 gives
    #   trans(y_ijk) = logit(E_ijk) + eps_ijk,
    #   eps_ijk ~ N(0, sigma^2 / N_ijk)
    # i.e. a constant residual on the logit scale whose variance is
    # inversely proportional to the arm sample size. Monolix fits
    # this by pre-multiplying both the data and the prediction by
    # sqrt(NOC) (Equation 2 and 'pred = logit(Q)*sqrt(NOC)') and
    # using a constant error model with SD a. The two formulations
    # are algebraically identical; this file uses the Equation 3
    # form, applying addSd / sqrt(N_ARM) directly in model() so the
    # observation stays on the interpretable logit scale.
    # ================================================================
    addSd <- 1.1
    label("Constant residual SD on the logit-of-fraction QoL scale for an arm of one participant; the per-observation SD is addSd / sqrt(N_ARM) (unitless)")  # Zou 2026 Table S1 'Additive residual error (a)' = 1.1 (RSE 3.95%; bootstrap median 1.068, 10th/90th 0.891-1.276)
  })

  model({
    # Study-arm inputs supplied per row:
    #   TRT    -- 0 = pembrolizumab arm, 1 = control arm
    #   N_ARM  -- number of participants contributing to the arm mean

    # 1. Baseline QoL for this study arm. Monolix Supplementary Codes
    #    section A:
    #      tE0   = logit(E0FE)
    #      tE0RE = tE0 + etaBSVE0 + etaBTAVE0/sqrt(NOC)
    #      E0    = expit(tE0RE)
    #    The between-arm term is divided by sqrt(N_ARM) because an arm
    #    mean's standard error shrinks with arm size; the between-study
    #    term is not.
    btavE0  <- eta_arm_e0 / sqrt(N_ARM)
    e0Logit <- logite0 + eta_study_e0 + btavE0
    e0Arm   <- expit(e0Logit)

    # 2. Maximal toxicity for this study arm. Treatment acts on the
    #    logit scale (Monolix logit-normal EmaxFE with a TRT
    #    coefficient); only between-STUDY variability was retained.
    emaxLogit <- logitemax + e_trt_emax * TRT + eta_study_emax
    emaxArm   <- expit(emaxLogit)

    # 3. Toxicity onset rate. No covariate and no random effect were
    #    retained on Kp (Monolix: no-variability, no covariate).
    kel <- exp(lkel)

    # 4. QoL improvement rate for this study arm. Treatment acts on
    #    the log scale (Monolix log-normal SLP with a TRT
    #    coefficient) and between-ARM variability was retained.
    slp <- exp(lslp + e_trt_slp * TRT + eta_arm_slp)

    # 5. QoL trajectory on the 0-1 fraction scale (Zou 2026 Eq. 1).
    qFrac <- e0Arm - emaxArm * (1 - exp(-kel * time)) + slp * time

    # 6. Saturation guard, verbatim from Supplementary Codes section A
    #    ("adding a saturation to avoid taking logit(0) (undefined)").
    #    qSat is what the logit below is taken of -- see the vignette
    #    Assumptions and deviations section: the supplement's own
    #    `pred` line applies logit to the UNSATURATED Q, which is what
    #    the guard exists to prevent, and which returns NaN once a
    #    high-baseline, high-slope arm carries Q past 1 at long times.
    qSat <- min(max(qFrac, 0.01), 0.99)

    # 7. Reported QoL score on the published 0-100 GHS/QoL scale. This
    #    is the linear-scale counterpart of the modelled observation
    #    and is what Zou 2026 Figures 1, 3 and 4A plot.
    qol_ghs <- 100 * qSat

    # 8. Observation on the logit scale with the arm-size-weighted
    #    residual of Zou 2026 Eq. 3 (SD = addSd / sqrt(N_ARM)). The
    #    residual is additive on the LOGIT scale, not on the 0-100
    #    score scale -- that is the whole point of the transform, and
    #    is what keeps a simulated arm mean inside 0-100.
    logit_qol_ghs  <- logit(qSat)
    addSdArm       <- addSd / sqrt(N_ARM)
    logit_qol_ghs ~ add(addSdArm)
  })
}
