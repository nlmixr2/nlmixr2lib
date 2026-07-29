Vargo_2014_statins_ezetimibe_mbma <- function() {
  description <- "MBMA. Literature-based meta-analysis dose-response model for percent change in low-density lipoprotein cholesterol (LDL-C) from baseline for six statins (atorvastatin, fluvastatin, lovastatin, pravastatin, rosuvastatin, simvastatin), ezetimibe, and statin-plus-ezetimibe combination therapy in adult dyslipidemia. Operates at the study-arm level over 245 trials (1,267 study-arm data points, 106,808 patients). Algebraic Emax-Hill (sigmoid) dose-response with statin-specific ED50 and shared sigmoidicity n=0.417 across statins; ezetimibe sigmoidicity is fixed to 1. Statin Emax depends on study-arm baseline LDL-C, baseline triglycerides, percentage with coronary heart disease (CHD), and binary cohort indicators for acute coronary syndrome (DIS_ACS) and heterozygous familial hypercholesterolemia (HeFH). Combination therapy is modelled via a sub-additive interaction coefficient gamma=0.523 (at maximal monotherapy effect the combined LDL-C reduction is about 7 percent smaller than the sum of the two monotherapies). Fluvastatin and lovastatin twice-daily and extended-release formulations multiply the statin ED50 by a fixed ratio (0.645 for fluvastatin; 0.59 for lovastatin). Between-study variances for Emax and ED50 were fixed to zero in the source paper, so the model has no eta IIV; the residual SD describes study-arm-mean variability and the suitable simulation scope is study-arm-mean percent change in LDL-C, not individual-subject concentrations."

  reference <- paste(
    "Vargo R, Adewale A, Behm MO, Mandema J, Kerbusch T.",
    "Prediction of clinical irrelevance of PK differences in atorvastatin",
    "using PK/PD models derived from literature-based meta-analyses.",
    "Clin Pharmacol Ther. 2014 Jul;96(1):101-109.",
    "doi:10.1038/clpt.2014.66.",
    sep = " "
  )
  vignette <- "Vargo_2014_statins_ezetimibe_mbma"
  units <- list(
    time          = "day (placeholder; the dose-response is steady-state and time-independent)",
    dosing        = "mg/day (per-arm daily dose of each drug supplied as DOSE_* covariate columns; the model is an MBMA dose-response and does not consume rxode2 dose events)",
    concentration = "fraction/fraction (signed fractional change in LDL-C from baseline; e.g. -0.65 = 65 percent LDL-C reduction). Output Cc is NOT a drug concentration; the slash in the unit string is to satisfy checkModelConventions parsing."
  )

  covariateData <- list(
    CONMED_ATORVASTATIN_DOSE = list(
      description        = "Per-arm atorvastatin daily dose (mg/day; 0 if atorvastatin is not in the regimen).",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "MBMA study-arm-level covariate (study-arm mean daily dose). The canonical register in inst/references/covariate-columns.md is for individual-level pop-PK covariates and does not directly fit MBMA study-arm-level dosing columns; this column mirrors the multi-drug Sadouki_2025 precedent (drug-specific dose / presence covariates documented inline rather than registered). Vargo 2014 Table 2 atorvastatin dose range across the 99 atorvastatin-monotherapy trials was 2.5-80 mg/day (median 20 mg/day).",
      source_name        = "Atorvastatin Dose (Vargo 2014 Table 2)"
    ),
    CONMED_FLV_DOSE = list(
      description        = "Per-arm fluvastatin daily dose (mg/day; 0 if fluvastatin is not in the regimen).",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "MBMA study-arm-level covariate. Vargo 2014 Table 2 fluvastatin dose range across the 33 fluvastatin-monotherapy trials was 2.5-60 mg/day (median 40 mg/day).",
      source_name        = "Fluvastatin Dose (Vargo 2014 Table 2)"
    ),
    CONMED_LOV_DOSE = list(
      description        = "Per-arm lovastatin daily dose (mg/day; 0 if lovastatin is not in the regimen).",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "MBMA study-arm-level covariate. Vargo 2014 Table 2 lovastatin dose range across the 26 lovastatin-monotherapy trials was 10-80 mg/day (median 40 mg/day).",
      source_name        = "Lovastatin Dose (Vargo 2014 Table 2)"
    ),
    CONMED_PRV_DOSE = list(
      description        = "Per-arm pravastatin daily dose (mg/day; 0 if pravastatin is not in the regimen).",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "MBMA study-arm-level covariate. Vargo 2014 Table 2 pravastatin dose range across the 67 pravastatin-monotherapy trials was 5-160 mg/day (median 20 mg/day).",
      source_name        = "Pravastatin Dose (Vargo 2014 Table 2)"
    ),
    CONMED_RSV_DOSE = list(
      description        = "Per-arm rosuvastatin daily dose (mg/day; 0 if rosuvastatin is not in the regimen).",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "MBMA study-arm-level covariate. Vargo 2014 Table 2 rosuvastatin dose range across the 43 rosuvastatin-monotherapy trials was 1-80 mg/day (median 10 mg/day).",
      source_name        = "Rosuvastatin Dose (Vargo 2014 Table 2)"
    ),
    CONMED_SMV_DOSE = list(
      description        = "Per-arm simvastatin daily dose (mg/day; 0 if simvastatin is not in the regimen).",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "MBMA study-arm-level covariate. Vargo 2014 Table 2 simvastatin dose range across the 91 simvastatin-monotherapy trials was 2.5-160 mg/day (median 20 mg/day).",
      source_name        = "Simvastatin Dose (Vargo 2014 Table 2)"
    ),
    CONMED_EZT_DOSE = list(
      description        = "Per-arm ezetimibe daily dose (mg/day; 0 if ezetimibe is not in the regimen).",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "MBMA study-arm-level covariate. Vargo 2014 Table 2 ezetimibe dose range across the 10 ezetimibe-monotherapy trials was 0.25-10 mg/day (median 10 mg/day). The covariates were not tested on the ezetimibe ED50 because the 10 mg dose (for which most data are available) is at the ezetimibe Emax.",
      source_name        = "Ezetimibe Dose (Vargo 2014 Table 2)"
    ),
    FORM_FLV_BID_XR = list(
      description        = "Indicator that the fluvastatin regimen is twice-daily or extended-release.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (once-daily immediate-release fluvastatin, the reference regimen)",
      notes              = "1 if the fluvastatin arm used b.i.d. dosing or an extended-release formulation; 0 otherwise. The model multiplies the fluvastatin ED50 by 0.645 when this flag is 1 (Vargo 2014 Table 3 row 'ED50,fluvastatin (b.i.d.|XR)/ED50,fluvastatin'; same ratio used for either regimen because the paper found b.i.d. and XR ED50 estimates were similar). The flag is meaningful only when CONMED_FLV_DOSE > 0; set to 0 for non-fluvastatin arms.",
      source_name        = "fluvastatin formulation/regimen (Vargo 2014 Figure 1)"
    ),
    FORM_LOV_BID_XR = list(
      description        = "Indicator that the lovastatin regimen is twice-daily or extended-release.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (once-daily immediate-release lovastatin, the reference regimen)",
      notes              = "1 if the lovastatin arm used b.i.d. dosing or an extended-release formulation; 0 otherwise. The model multiplies the lovastatin ED50 by 0.59 when this flag is 1 (Vargo 2014 Table 3 row 'ED50,lovastatin (b.i.d.|XR)/ED50,lovastatin'). Meaningful only when CONMED_LOV_DOSE > 0.",
      source_name        = "lovastatin formulation/regimen (Vargo 2014 Figure 1)"
    ),
    LDLC = list(
      description        = "Study-arm mean baseline (pre-treatment) LDL-C concentration.",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-arm baseline LDL-C used in the statin Emax covariate equation (Vargo 2014 Eq 4): the log-ratio (LDLC/180) modifies Emax_statin. Vargo 2014 Table 2 LDL-C ranged 106-349 mg/dL across treatment groups; overall median 181 mg/dL. The reference 180 mg/dL is the centring constant used by the paper and corresponds to the typical-patient definition in the Results section. Source paper uses 'LDL.base'.",
      source_name        = "LDL.base (Vargo 2014 Table 2 / Eq 4)"
    ),
    TRIG = list(
      description        = "Study-arm mean baseline (pre-treatment) triglyceride concentration.",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-arm baseline triglyceride used in the statin Emax covariate equation (Vargo 2014 Eq 4): the log-ratio (TRIG/180) modifies Emax_statin. Vargo 2014 Table 2 TG ranged 59-660 mg/dL across treatment groups; overall median 168 mg/dL. The reference 180 mg/dL is the centring constant used by the paper. Source paper uses 'TG.base'.",
      source_name        = "TG.base (Vargo 2014 Table 2 / Eq 4)"
    ),
    DIS_CHD_PERCENT = list(
      description        = "Percentage of patients in the study arm with coronary heart disease (CHD); range 0-100.",
      units              = "%",
      type               = "continuous",
      reference_category = NULL,
      notes              = "MBMA study-arm-level covariate; the canonical register in inst/references/covariate-columns.md is for individual-level pop-PK covariates and does not directly fit MBMA study-arm aggregate-percentage columns. DIS_CHD_PERCENT enters Vargo 2014 Eq 4 with linear coefficient Emax,4 = -0.000649 per percentage point (so a 24-percent CHD arm reduces statin Emax by 0.000649 x 24 = 0.016). Default 0 in simulation = healthy adult cohort; the paper's typical-patient definition (Figure 2 caption) uses 24 percent CHD.",
      source_name        = "CHD% (Vargo 2014 Eq 4 / Table 3)"
    ),
    DIS_ACS = list(
      description        = "Indicator that the study arm enrolled patients with acute coronary syndrome.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-DIS_ACS arm)",
      notes              = "MBMA study-arm-level covariate; 1 if the arm enrolled an DIS_ACS cohort, 0 otherwise. Vargo 2014 Eq 4 gives an additive shift of -0.117 on Emax_statin in DIS_ACS arms (greater statin LDL-C lowering in DIS_ACS patients). Encoded as DIS_ACS in the Vargo source equation (note the Table 3 row label 'Emax,3 (DIS_ACS)' is a transcription typo; the equation in Eq 4 distinguishes this coefficient from the LDL.base and TG.base coefficients).",
      source_name        = "DIS_ACS=yes (Vargo 2014 Eq 4 / Table 3)"
    ),
    DIS_HEFH = list(
      description        = "Indicator that the study arm enrolled patients with heterozygous familial hypercholesterolemia.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-HeFH arm)",
      notes              = "MBMA study-arm-level covariate; 1 if the arm enrolled a HeFH cohort, 0 otherwise. Vargo 2014 Eq 4 gives an additive shift of +0.127 on Emax_statin in HeFH arms (smaller statin LDL-C lowering in HeFH patients). Note Table 3 row label 'Emax,3 (HeFH)' is a transcription typo; the equation in Eq 4 distinguishes this coefficient from the DIS_ACS coefficient.",
      source_name        = "HeFH=yes (Vargo 2014 Eq 4 / Table 3)"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 106808L,
    n_studies       = 245L,
    n_data_points   = 1267L,
    age_range       = "26-77 years (overall range across treatment groups per Vargo 2014 Table 2)",
    age_median      = "57 years (overall median across treatment groups per Vargo 2014 Table 2)",
    disease_state   = "adults with dyslipidemia (including arms with coronary heart disease, acute coronary syndrome, heterozygous familial hypercholesterolemia, and a placebo reference) -- mixed cohorts pooled at the study-arm level",
    dose_range      = "atorvastatin 2.5-80 mg/day; fluvastatin 2.5-60 mg/day; lovastatin 10-80 mg/day; pravastatin 5-160 mg/day; rosuvastatin 1-80 mg/day; simvastatin 2.5-160 mg/day; ezetimibe 0.25-10 mg/day (per Vargo 2014 Table 2)",
    regimens        = "Primarily once-daily immediate-release. Fluvastatin twice-daily in 6 trials (295 patients) and extended-release in 10 trials (1,558 patients); lovastatin twice-daily in 5 trials (3,562 patients) and extended-release in 2 trials (314 patients); pravastatin twice-daily in 6 trials (493 patients). Same fixed ED50 ratio applied for b.i.d. and XR formulations of fluvastatin (0.645) and lovastatin (0.59); no significant b.i.d./XR effect detected for pravastatin.",
    baseline_ldlc   = "median 181 mg/dL (range across treatment groups 106-349 mg/dL per Vargo 2014 Table 2)",
    baseline_hdlc   = "median 48 mg/dL (range 24-77 mg/dL per Vargo 2014 Table 2)",
    baseline_trig   = "median 168 mg/dL (range 59-660 mg/dL per Vargo 2014 Table 2)",
    sex_female_pct  = NA_real_,
    race_ethnicity  = "Asian vs non-Asian was tested as a covariate but not retained in the final model.",
    regions         = "International; trials from the public literature plus FDA and EMA-registered submissions.",
    notes           = "MBMA at the study-arm level: each data point is the mean LDL-C response in a group of patients at a particular time point in a single trial arm, with the associated normally distributed variance weighted by arm sample size. The model is intended for simulating study-arm-mean percent change in LDL-C and is NOT suitable for individual-subject simulation. Between-trial variances on Emax and ED50 were estimated and found non-significant; the paper fixed them to zero. Steady-state effect was confirmed after at least 4 weeks of treatment so the model has no time component (time did not significantly modify Emax). See Vargo 2014 Methods and Table 2 for the cohort breakdown by treatment."
  )

  ini({
    # ============================================================
    # Statin Emax intercept and covariate coefficients (Vargo 2014 Eq 4,
    # Table 3). The paper writes Emax_statin = Emax,1 + Emax,2*log(LDL/180)
    # + Emax,3*log(TG/180) + Emax,4*CHD% + Emax,5*DIS_ACS + Emax,6*HeFH
    # (the Table 3 row labels collapse coefficients 3-6 to 'Emax,3'; the
    # equation in Eq 4 uses distinct symbols; the implementation below uses
    # six distinct coefficients matching the values in Table 3).
    #
    # Sign convention: Emax,1 = -0.758 in Table 3 (signed, negative = LDL-C
    # reduction). In Eq 4 all covariate terms are added on this signed
    # scale, so Emax_statin is the signed fractional reduction in LDL-C
    # the statin can achieve. In Eq 3, however, f_statin and f_ezetimibe
    # appear as MAGNITUDES so that 1 - gamma*f_statin gives sub-additive
    # combination (paper Results: "reduced by ~7% of the sum of the two
    # maximal monotherapy effects" -- with signed values the same formula
    # gives MORE reduction than the sum, contradicting the prose).
    # The model() block evaluates Eq 4 on the signed scale and takes
    # abs() to convert to magnitudes before evaluating Eq 3.
    # ============================================================

    emax_statin_int <- -0.758
    label("Statin Emax intercept (Emax,1; fractional change in LDL-C, signed)")  # Vargo 2014 Table 3

    e_ldlc_emax_statin <- -0.14
    label("LDL.base coefficient on statin Emax (per natural-log unit of LDL/180)")  # Vargo 2014 Table 3 Emax,2 (LDL.base)

    e_trig_emax_statin <- 0.0506
    label("TG.base coefficient on statin Emax (per natural-log unit of TG/180)")  # Vargo 2014 Table 3 Emax,3 (TG.base)

    e_chd_emax_statin <- -0.000649
    label("CHD% coefficient on statin Emax (per percentage point of CHD)")  # Vargo 2014 Table 3 Emax,4 (CHD%)

    e_acs_emax_statin <- -0.117
    label("DIS_ACS=yes additive shift on statin Emax (signed fractional)")  # Vargo 2014 Table 3 Emax,5 (DIS_ACS)

    e_hefh_emax_statin <- 0.127
    label("HeFH=yes additive shift on statin Emax (signed fractional)")  # Vargo 2014 Table 3 Emax,6 (HeFH)

    # ============================================================
    # Statin ED50 (mg/day) -- one per statin (Vargo 2014 Table 3).
    # Encoded on log scale so the population-typical predictions are
    # back-transformed via exp(led50_*) in model().
    # ============================================================
    led50_atv <- log(15.2)
    label("Atorvastatin ED50 (mg/day)")  # Vargo 2014 Table 3 ED50,atorvastatin

    led50_flv <- log(347)
    label("Fluvastatin ED50 (mg/day; once-daily immediate-release reference)")  # Vargo 2014 Table 3 ED50,fluvastatin

    led50_lov <- log(114)
    label("Lovastatin ED50 (mg/day; once-daily immediate-release reference)")  # Vargo 2014 Table 3 ED50,lovastatin

    led50_prv <- log(145)
    label("Pravastatin ED50 (mg/day)")  # Vargo 2014 Table 3 ED50,pravastatin

    led50_rsv <- log(4.96)
    label("Rosuvastatin ED50 (mg/day)")  # Vargo 2014 Table 3 ED50,rosuvastatin

    led50_smv <- log(36.7)
    label("Simvastatin ED50 (mg/day)")  # Vargo 2014 Table 3 ED50,simvastatin

    # ============================================================
    # Twice-daily / extended-release ED50 ratios (Vargo 2014 Table 3).
    # Fixed to the literature ratio; no IIV.
    # ============================================================
    ratio_red_flv <- fixed(0.645)
    label("Fluvastatin b.i.d./XR ED50 ratio (relative to q.d. IR; FIXED)")  # Vargo 2014 Table 3 ED50,fluvastatin (b.i.d.|XR)/ED50,fluvastatin

    ratio_red_lov <- fixed(0.59)
    label("Lovastatin b.i.d./XR ED50 ratio (relative to q.d. IR; FIXED)")  # Vargo 2014 Table 3 ED50,lovastatin (b.i.d.|XR)/ED50,lovastatin

    # ============================================================
    # Sigmoidicity (shared across statins; ezetimibe sigmoidicity is
    # FIXED to 1 by the paper Methods).
    # ============================================================
    ln_statin <- log(0.417)
    label("Statin sigmoidicity (Hill exponent; shared across statins)")  # Vargo 2014 Table 3 N

    # ============================================================
    # Ezetimibe Emax and ED50 (Vargo 2014 Table 3).
    # ============================================================
    emax_ezt <- -0.184
    label("Ezetimibe Emax (signed fractional change in LDL-C)")  # Vargo 2014 Table 3 Emax (ezetimibe)

    led50_ezt <- log(0.228)
    label("Ezetimibe ED50 (mg/day)")  # Vargo 2014 Table 3 ED50,ezetimibe

    ln_ezt <- fixed(log(1))
    label("Ezetimibe sigmoidicity (Hill exponent; FIXED to 1 per paper Methods)")  # Vargo 2014 Methods, ezetimibe dose-response paragraph: "with the sigmoidicity factor (n) fixed to 1"

    # ============================================================
    # Statin + ezetimibe interaction coefficient gamma (Vargo 2014 Eq 3
    # / Table 3). The combination model is f_combo = f_statin +
    # f_ezetimibe * (1 - gamma * f_statin) where f_* are magnitudes.
    # gamma > 0 implies sub-additive interaction at non-zero statin dose.
    # ============================================================
    gamma_se <- 0.523
    label("Statin x ezetimibe interaction coefficient (Eq 3)")  # Vargo 2014 Table 3 gamma

    # ============================================================
    # Between-study variances. The paper estimated between-trial omega
    # for Emax and ED50 of both statins and ezetimibe; none were
    # statistically significant and all were fixed to zero in the final
    # model. The MBMA therefore has no eta IIV. No eta parameters in
    # ini() per the paper's final-model structure (Vargo 2014 Results,
    # 'Between-trial variances were estimated ... fixed to zero').
    # ============================================================

    # ============================================================
    # Residual error. Vargo 2014 Eq 1: epsilon_ijt ~ N(0, sigma^2 / N_ij)
    # with sigma the SD of the percentage change from baseline. The paper
    # reports sigma = 0.154 (Table 3) on the fractional scale (Methods:
    # 'fractional change (%/100)'). At study-arm level the sample-size
    # weighting is applied by treating sigma/sqrt(N) as the effective
    # observation noise; this model exposes the unweighted sigma and
    # leaves the per-arm reweighting to downstream simulation code.
    # Residual rho (within-arm correlation between time points) = 0.667
    # describes correlation across repeated observation time points in
    # the same arm; the model has no time component so rho does not
    # apply here. The model emits a single steady-state observation per
    # arm so rho is informational only and not encoded.
    # ============================================================
    addSd <- 0.154
    label("Additive residual SD on signed fractional change in LDL-C (study-arm scale; unweighted sigma per Eq 1)")  # Vargo 2014 Table 3 sigma
  })

  model({
    # ----- Reference constants used inside derived quantities -----
    # Reference baseline LDL-C and TG (mg/dL) for the log-ratio covariate
    # terms in Eq 4: the paper centres on 180 mg/dL ('a patient population
    # with a mean baseline LDL-C of 180 mg/dL, a mean baseline triglyceride
    # level of 180 mg/dL').
    ref_ldlc_mgdl <- 180
    ref_trig_mgdl <- 180

    # ----- Statin-specific structural rates -----
    n_statin <- exp(ln_statin)
    n_ezt    <- exp(ln_ezt)

    # ED50 per statin, with the fluvastatin and lovastatin b.i.d./XR
    # ratio applied multiplicatively when the regimen flag is 1.
    # FORM_FLV_BID_XR / FORM_LOV_BID_XR are binary so the algebraic form
    # (1 + (ratio - 1) * flag) yields ratio when flag=1 and 1 when flag=0.
    ed50_atv <- exp(led50_atv)
    ed50_flv <- exp(led50_flv) * (1 + (ratio_red_flv - 1) * FORM_FLV_BID_XR)
    ed50_lov <- exp(led50_lov) * (1 + (ratio_red_lov - 1) * FORM_LOV_BID_XR)
    ed50_prv <- exp(led50_prv)
    ed50_rsv <- exp(led50_rsv)
    ed50_smv <- exp(led50_smv)
    ed50_ezt <- exp(led50_ezt)

    # ----- Statin Emax (Eq 4) on the signed fractional scale -----
    # Negative values mean LDL-C reduction.
    emax_statin_signed <- emax_statin_int +
      e_ldlc_emax_statin * log(LDLC / ref_ldlc_mgdl) +
      e_trig_emax_statin * log(TRIG / ref_trig_mgdl) +
      e_chd_emax_statin  * DIS_CHD_PERCENT +
      e_acs_emax_statin  * DIS_ACS +
      e_hefh_emax_statin * DIS_HEFH

    # Convert to MAGNITUDE (positive) for use in Eq 3. The paper's
    # Methods says ezetimibe Emax is also reported as a signed negative
    # value (-0.184); both magnitudes are positive here.
    emax_statin_mag <- -emax_statin_signed
    emax_ezt_mag    <- -emax_ezt

    # ----- Per-statin dose-response (Eq 2) on the magnitude scale -----
    # Each f_<statin> is non-zero only when its DOSE_<statin> > 0.
    # Standard MBMA usage: only one statin column is non-zero per study
    # arm (the trials in Vargo 2014 Table 2 are single-statin), so the
    # sum collapses to the active statin's term. If a future user codes
    # an arm with two statins simultaneously the model returns the
    # additive sum of their individual Emax contributions, which is
    # outside the paper's calibration range.
    f_atv <- emax_statin_mag * (CONMED_ATORVASTATIN_DOSE ^ n_statin) /
             ((CONMED_ATORVASTATIN_DOSE ^ n_statin) + (ed50_atv ^ n_statin))
    f_flv <- emax_statin_mag * (CONMED_FLV_DOSE ^ n_statin) /
             ((CONMED_FLV_DOSE ^ n_statin) + (ed50_flv ^ n_statin))
    f_lov <- emax_statin_mag * (CONMED_LOV_DOSE ^ n_statin) /
             ((CONMED_LOV_DOSE ^ n_statin) + (ed50_lov ^ n_statin))
    f_prv <- emax_statin_mag * (CONMED_PRV_DOSE ^ n_statin) /
             ((CONMED_PRV_DOSE ^ n_statin) + (ed50_prv ^ n_statin))
    f_rsv <- emax_statin_mag * (CONMED_RSV_DOSE ^ n_statin) /
             ((CONMED_RSV_DOSE ^ n_statin) + (ed50_rsv ^ n_statin))
    f_smv <- emax_statin_mag * (CONMED_SMV_DOSE ^ n_statin) /
             ((CONMED_SMV_DOSE ^ n_statin) + (ed50_smv ^ n_statin))

    f_statin <- f_atv + f_flv + f_lov + f_prv + f_rsv + f_smv

    # ----- Ezetimibe dose-response (Eq 2 with n=1) -----
    f_ezetimibe <- emax_ezt_mag * (CONMED_EZT_DOSE ^ n_ezt) /
                   ((CONMED_EZT_DOSE ^ n_ezt) + (ed50_ezt ^ n_ezt))

    # ----- Combination (Eq 3) on the magnitude scale -----
    f_combination_mag <- f_statin + f_ezetimibe * (1 - gamma_se * f_statin)

    # ----- Output: signed fractional change in LDL-C from baseline ----
    # Cc is overloaded here as the single-output observation per nlmixr2lib
    # convention; it is NOT a drug concentration but the signed fractional
    # change in LDL-C from pretreatment baseline (e.g. Cc = -0.65 means a
    # 65 percent LDL-C reduction). Multiply by 100 downstream for percent
    # display.
    Cc <- -f_combination_mag

    Cc ~ add(addSd)
  })
}
