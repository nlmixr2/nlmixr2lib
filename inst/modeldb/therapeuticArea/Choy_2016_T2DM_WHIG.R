Choy_2016_T2DM_WHIG <- function() {
  description <- paste0(
    "Semi-mechanistic disease-progression model for type 2 diabetes ",
    "(WHIG: Weight, HbA1c, Insulin, Glucose). Body-weight turnover under ",
    "diet+exercise and placebo drives insulin sensitivity, which together ",
    "with beta-cell function drives a fasting serum insulin (FSI) and ",
    "fasting plasma glucose (FPG) homeostatic feedback (steady-state ",
    "quadratic). FPG and an additional postprandial-glucose contribution ",
    "feed a three-compartment transit chain producing total HbA1c. Built ",
    "from the placebo arm of NCT00236600 (181 obese newly-diagnosed adults ",
    "with T2DM on diet+exercise for 67 weeks)."
  )
  reference <- paste(
    "Choy S, Kjellsson MC, Karlsson MO, de Winter W.",
    "Weight-HbA1c-Insulin-Glucose Model for Describing Disease Progression",
    "of Type 2 Diabetes.",
    "CPT Pharmacometrics Syst Pharmacol. 2016;5(1):11-19.",
    "doi:10.1002/psp4.12051"
  )
  vignette <- "Choy_2016_T2DM_WHIG"
  units <- list(
    time          = "day",
    dosing        = "n/a (placebo / lifestyle intervention only)",
    concentration = "weight (kg); FSI (uIU/mL); FPG (mmol/L); HbA1c (%)"
  )

  # No subject-level covariates modify a parameter. The model is fully driven
  # by time and a fixed start-of-active-treatment marker (tTRT). The run-in
  # boundary at tTRT = 42 days is encoded structurally in model(), so no
  # subject indicator column is needed in the event table.
  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 181L,
    n_studies      = 1L,
    age_range      = "18-75 years",
    age_median     = NA_character_,
    weight_range   = "72-159 kg (baseline)",
    weight_median  = "104 kg (baseline)",
    sex_female_pct = 100 * 114 / 181,
    race_ethnicity = NA_character_,
    disease_state  = paste0(
      "Obese (BMI 27-50 kg/m^2), newly diagnosed treatment-naive type 2 ",
      "diabetes mellitus; baseline median FSI 17.8 uIU/mL (range 3.3-79.5), ",
      "FPG 7.6 mmol/L (range 5.0-14.2), HbA1c 6.7% (range 5.3-9.1)."
    ),
    dose_range     = paste0(
      "Placebo arm only: 6-week placebo run-in followed by 60-week ",
      "placebo + lifestyle (individualised energy-deficient diet, ",
      "behavioural modification, physical activity) maintenance phase."
    ),
    regions        = "Sweden",
    notes          = paste0(
      "Placebo arm of NCT00236600 (topiramate T2DM weight-loss study). ",
      "Mean weight change at end of study -4.1 kg; baseline beta-cell ",
      "function 61% of normal; baseline insulin sensitivity 25% of normal."
    )
  )

  # ----- ini(): all parameters from Table 1 of Choy et al. 2016 -----
  # Variability interpretation per Table 1 footnote a/b (paper, p.15):
  #   - Log-normally distributed parameters report CV% in the table; the
  #     packaged ini() variance is log(1 + CV^2).
  #   - Normally distributed parameters listed under footnote b (is0, b0, rb,
  #     efde, efp, efup) have their CV column "reported as absolute values"
  #     rather than as CV%. The cleanest reading consistent with reasonable
  #     posterior magnitudes is: for parameters whose mean is on the logit
  #     scale (is0, b0, rb) the column value IS the SD on that logit scale;
  #     for parameters whose mean is on a percent-effect scale (efde, efp,
  #     efup) the column value IS the absolute SD expressed as a percentage,
  #     so SD = column/100 on the same percent scale as the mean. Both
  #     readings give biologically sensible IIV magnitudes; see the vignette
  #     Assumptions and deviations section for details.
  #   - Supplementary Appendix S2 (variance-covariance correlation matrix,
  #     paper p.15 footnote a) is not on disk; the packaged model uses a
  #     diagonal OMEGA structure. Simulation reproduces typical-value
  #     predictions exactly but cannot exactly reproduce the published 95%
  #     prediction-interval widths in Figures 3 and 4 (which depend on the
  #     missing correlations).
  ini({
    # ----- Weight model (Choy 2016 Eq. 1-4; Table 1) -----
    t_half_wgt <- 96.9     ; label("Half-life of the weight compartment (day)")  # Choy 2016 Table 1
    lblwt      <- log(104) ; label("Typical baseline weight (kg)")               # Choy 2016 Table 1

    # ----- Insulin sensitivity (Choy 2016 Eq. 5-6) -----
    is0         <- 1.1         ; label("Baseline insulin sensitivity, logit scale (unitless); IS_baseline = 1/(1+exp(is0))")  # Choy 2016 Table 1; 1.1 -> 25% of normal
    lscaleefs   <- log(0.0514) ; label("Log scaling factor of change in weight on insulin sensitivity (1/kg)")                # Choy 2016 Table 1

    # ----- beta-cell function and natural disease progression (Choy 2016 Eq. 7) -----
    b0 <- -0.446 ; label("Baseline beta-cell function, logit scale (unitless); B_baseline = 1/(1+exp(b0))")  # Choy 2016 Table 1; -0.446 -> 61% of normal
    rb <-  0.209 ; label("Rate of baseline beta-cell function decrease (logits per year)")                   # Choy 2016 Table 1

    # ----- Treatment effect on beta-cell function (Choy 2016 Eq. 8-9) -----
    lefbmax <- log(0.171) ; label("Log maximal relative increase of beta-cell function during the EFB peak (unitless)")  # Choy 2016 Table 1
    sefbi   <- -3.69      ; label("Shape parameter for the logistic increase of the beta-cell treatment effect (unitless)")  # Choy 2016 Table 1
    sefbd   <-  8.05      ; label("Shape parameter for the logistic decrease of the beta-cell treatment effect (unitless)")  # Choy 2016 Table 1
    lefb50  <- log(190)   ; label("Log time at half of the EFB logistic decline (day)")                                       # Choy 2016 Table 1

    # ----- Weight-input treatment effects (Choy 2016 text under Weight change) -----
    efde <- 4.08 ; label("Effect of diet+exercise on weight input, % reduction (active from t = 0)")     # Choy 2016 Table 1
    efp  <- 2.28 ; label("Effect of placebo on weight input, % reduction (active from t = tTRT)")        # Choy 2016 Table 1
    efup <- 2.99 ; label("Counter-effect on weight input, % increase per year (grows linearly with t)")  # Choy 2016 Table 1

    # ----- HbA1c model (Choy 2016 Eq. 13-17) -----
    lkin_hba1c <- log(0.0129) ; label("Log rate constant for HbA1c-compartment production from FPG (%/d per mmol/L)")  # Choy 2016 Table 1
    lppg       <- log(0.0709) ; label("Log residual HbA1c production rate independent of FPG (%/d)")                   # Choy 2016 Table 1
    scaleppg   <- 0.963       ; label("Scaling factor on PPG for t > 0 (unitless)")                                    # Choy 2016 Table 1
    lmtt       <- log(38.9)   ; label("Log mean transit time of the HbA1c transit chain (day)")                        # Choy 2016 Table 1

    # ----- IIV (diagonal OMEGA; S2 correlation matrix not on disk) -----
    etalblwt     ~ 0.02109   # Choy 2016 Table 1; log-normal, CV 14.6% -> log(1 + 0.146^2)
    etais0       ~ 0.09303   # Choy 2016 Table 1; normal SD 0.305 on the logit scale (footnote b)
    etalscaleefs ~ 0.37570   # Choy 2016 Table 1; log-normal, CV 67% -> log(1 + 0.67^2)
    etab0        ~ 1.96000   # Choy 2016 Table 1; normal SD 1.4 on the logit scale (footnote b)
    etalefbmax   ~ 0.22158   # Choy 2016 Table 1; log-normal, CV 49.9% -> log(1 + 0.499^2)
    etalefb50    ~ 0.11453   # Choy 2016 Table 1; log-normal, CV 34.9% -> log(1 + 0.349^2)
    etarb        ~ 0.04410   # Choy 2016 Table 1; normal SD 0.21 on the logit-rate scale (footnote b)
    etalppg      ~ 0.02358   # Choy 2016 Table 1; log-normal, CV 15.4% -> log(1 + 0.154^2)
    etaefde      ~ 0.12674   # Choy 2016 Table 1; normal SD 0.356 percent-points (footnote b; column/100)
    etaefp       ~ 0.16160   # Choy 2016 Table 1; normal SD 0.402 percent-points (footnote b)
    etaefup      ~ 0.55354   # Choy 2016 Table 1; normal SD 0.744 percent-points (footnote b)

    # ----- Proportional residual errors (Choy 2016 Table 1) -----
    # Table 1 reports the NONMEM $SIGMA variance; the prop() SD is sqrt(variance).
    propSd_WGT   <- sqrt(0.00919) ; label("Proportional residual SD on weight (fraction)")   # Choy 2016 Table 1; sqrt(0.00919) = 0.0959
    propSd_FSI   <- sqrt(0.262)   ; label("Proportional residual SD on FSI (fraction)")      # Choy 2016 Table 1; sqrt(0.262) = 0.512
    propSd_FPG   <- sqrt(0.0688)  ; label("Proportional residual SD on FPG (fraction)")      # Choy 2016 Table 1; sqrt(0.0688) = 0.262
    propSd_HbA1c <- sqrt(0.0241)  ; label("Proportional residual SD on HbA1c (fraction)")    # Choy 2016 Table 1; sqrt(0.0241) = 0.155
  })

  model({
    # ===== Fixed structural / HOMA constants =====
    # tTRT: start of active treatment phase = end of 6-week run-in.
    # Methods (Study design): the placebo step (efp) sets in at week 6, the
    # diet+exercise step (efde) sets in at week 0.
    tTRT <- 42  # day

    # KinFSI/KoutFSI = 7.8 uIU/mL: fixed by the healthy-subject HOMA2
    # convention (FSI_SS = 7.8 uIU/mL at FPG_SS = 4.5 mmol/L when B = 1,
    # EFB = 1, IS = 1). Choy 2016 Methods (FSI-FPG homeostatic feedback
    # model); Wallace 2004 HOMA2.
    KinKoutFSI <- 7.8

    # KinFPG/KoutFPG = 35.1 (mmol/L * uIU/mL): derived as
    # 4.5 mmol/L (healthy FPG_SS) * 7.8 uIU/mL (healthy FSI_SS) = 35.1.
    # Choy 2016 Methods (FSI-FPG homeostatic feedback model).
    KinKoutFPG <- 35.1

    # FPG floor below which fasting hepatic glucose output no longer
    # stimulates insulin secretion (HOMA convention; Choy 2016 cites Levy
    # 1998 and Wallace 2004).
    FPG_floor <- 3.5  # mmol/L

    # ===== Individual parameters =====
    blwt_i       <- exp(lblwt + etalblwt)
    is0_i        <- is0 + etais0
    scaleefs_i   <- exp(lscaleefs + etalscaleefs)
    b0_i         <- b0 + etab0
    efbmax_i     <- exp(lefbmax + etalefbmax)
    efb50_i      <- exp(lefb50 + etalefb50)
    rb_i         <- rb + etarb
    kin_hba1c_i  <- exp(lkin_hba1c)
    ppg_i        <- exp(lppg + etalppg)
    mtt_i        <- exp(lmtt)
    efde_i       <- efde + etaefde
    efp_i        <- efp  + etaefp
    efup_i       <- efup + etaefup

    # ===== Weight model (Eq. 1-4) =====
    # EFW_t: net effect on weight input, normalised to 1 at t = 0.
    #   efde_i (% reduction) acts from t = 0 onwards (D&E step).
    #   efp_i  (% reduction) acts from t = tTRT onwards (placebo step).
    #   efup_i (%/year increase) grows linearly with t (time in days).
    efdeP_step <- 1 - efde_i / 100 - (efp_i / 100) * (t >= tTRT)
    efup_t     <- 1 + (efup_i / 100) * (t / 365)
    EFW_t      <- efup_t * efdeP_step

    # Weight turnover. At steady state with EFW_t = 1, weight = blwt_i.
    kout_wgt <- log(2) / t_half_wgt
    d/dt(weight) <- kout_wgt * (EFW_t * blwt_i - weight)
    weight(0)    <- blwt_i

    WGT  <- weight
    dwgt <- weight - blwt_i  # absolute change from individual baseline

    # ===== Insulin sensitivity (Eq. 5-6) =====
    # EFS = 1 - scaleefs * dwgt: weight loss (dwgt < 0) raises IS.
    EFS <- 1 - scaleefs_i * dwgt
    # IS_baseline uses the "inverse logit" convention of Choy 2016 (Methods,
    # FSI-FPG homeostatic feedback model): IS = 1/(1+exp(is0)).
    # At is0 = 1.1 this gives 0.25 = 25% of normal, matching paper Results.
    IS_baseline <- 1 / (1 + exp(is0_i))
    IS          <- IS_baseline * EFS

    # ===== beta-cell function and disease progression (Eq. 7) =====
    # Natural logistic decline of beta-cell function: with rb > 0 the logit
    # grows over time and B decreases.
    B_logit <- b0_i + rb_i * (t / 365)
    B       <- 1 / (1 + exp(B_logit))

    # ===== Treatment effect on beta-cell function (Eq. 8-9) =====
    # EFB = 1 + efbmax * EFBI(t) * EFBD(t).
    # EFBI: logistic increase centred at t = tTRT (rises from ~0 to ~1 as
    # t crosses tTRT; sefbi is negative so the sigmoid is increasing in t).
    # EFBD: logistic decrease centred at t = efb50 (falls from ~1 to ~0 as
    # t crosses efb50; sefbd is positive so the sigmoid is decreasing in t).
    EFBI <- 1 / (1 + exp(sefbi * (t - tTRT)))
    EFBD <- 1 / (1 + exp(sefbd * (t - efb50_i)))
    EFB  <- 1 + efbmax_i * EFBI * EFBD

    # Combined beta-cell scaling for FSI production.
    Beff <- B * EFB

    # ===== FSI-FPG homeostatic feedback (Eq. 10-12; QSS quadratic) =====
    # Under the paper's short-term steady-state assumption,
    #   FSI_SS = 7.8 * Beff * (FPG - FPG_floor)
    #   FPG_SS = 35.1 / (IS * FSI_SS)
    # which (eliminating FSI) gives the quadratic
    #   FPG^2 - FPG_floor * FPG - KinKoutFPG / (KinKoutFSI * IS * Beff) = 0.
    K_eff <- IS * Beff
    qC    <- KinKoutFPG / (KinKoutFSI * K_eff)
    FPG   <- (FPG_floor + sqrt(FPG_floor * FPG_floor + 4 * qC)) / 2
    FSI   <- KinKoutFSI * Beff * (FPG - FPG_floor)

    # ===== HbA1c transit chain (Eq. 13-17) =====
    # Three transit compartments with shared kout = 3 / MTT. The total
    # HbA1c is the sum across compartments (Eq. 13). PPG contribution is
    # scaled by scaleppg for t > 0 (diet+exercise reduces the postprandial
    # fraction).
    kout_hba1c <- 3 / mtt_i
    ppg_eff    <- ppg_i * (1 - (1 - scaleppg) * (t > 0))

    d/dt(transit1) <- (kin_hba1c_i * FPG + ppg_eff) - kout_hba1c * transit1
    d/dt(transit2) <- kout_hba1c * (transit1 - transit2)
    d/dt(transit3) <- kout_hba1c * (transit2 - transit3)

    # Steady-state initial conditions evaluated at t = 0 with the t = 0
    # baseline FPG and the unscaled PPG (the run-in is short relative to
    # MTT ~ 39 d, so the chain starts at its pre-study steady state).
    B_baseline_init  <- 1 / (1 + exp(b0_i))
    BLFPG_K          <- IS_baseline * B_baseline_init
    BLFPG_qC         <- KinKoutFPG / (KinKoutFSI * BLFPG_K)
    BLFPG            <- (FPG_floor + sqrt(FPG_floor * FPG_floor + 4 * BLFPG_qC)) / 2
    BLA1             <- (kin_hba1c_i * BLFPG + ppg_i) / kout_hba1c
    transit1(0) <- BLA1
    transit2(0) <- BLA1
    transit3(0) <- BLA1

    HbA1c <- transit1 + transit2 + transit3

    # ===== Observations and residual error =====
    WGT   ~ prop(propSd_WGT)
    FSI   ~ prop(propSd_FSI)
    FPG   ~ prop(propSd_FPG)
    HbA1c ~ prop(propSd_HbA1c)
  })
}
