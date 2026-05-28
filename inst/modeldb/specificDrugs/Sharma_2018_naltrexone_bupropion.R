Sharma_2018_naltrexone_bupropion <- function() {
  description <- "Dose- and time-dependent population pharmacodynamic (DTPD) body-weight model for the naltrexone/bupropion fixed-dose combination (Contrave) in obese and overweight adults under lifestyle intervention, based on 4591 subjects pooled from six Contrave clinical trials (placebo and active-treatment arms). Indirect-response body-weight model with linear NHANES-derived disease progression, inverse-Bateman lifestyle-intervention stimulation of body-weight loss (kout), and a combined Emax dose- and time-dependent inhibitory drug effect; diabetes (T2DM) and race covariates on key parameters. Does not include the linked Markov dropout layer (Tr10, Tr01, Tr12, Tr02) of Table 3, and does not include the PPPD concentration-driven variant (whose underlying naltrexone/bupropion PopPK model is an unpublished internal Takeda report)."
  reference <- "Sharma VD, Combes FP, Vakilynejad M, Lahu G, Lesko LJ, Trame MN. Model-Based Approach to Predict Adherence to Protocol During Antiobesity Trials. The Journal of Clinical Pharmacology. 2018;58(2):240-253. doi:10.1002/jcph.994"
  vignette <- "Sharma_2018_naltrexone_bupropion"
  units <- list(
    time = "week",
    dosing = "mg",
    concentration = "kg"  # body weight (the single-output "Cc" canonical carries kg, not a drug concentration; see model() comment)
  )

  covariateData <- list(
    T2DM = list(
      description = "Type-2 diabetes mellitus comorbidity indicator at study entry",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "1 = obese subject with type-2 diabetes mellitus (from Sharma 2018 study 6, NB-304, which was the only included trial enrolling T2DM subjects -- 10.9% of the pooled 4591-subject analysis population); 0 = obese subject without diabetes (studies 1-5). Time-fixed per subject. Diabetes was implemented as a fixed covariate on Emax, kout, kpro, and baseline BW per Sharma 2018 Results (Covariate Analysis).",
      source_name = "T2DM"
    ),
    RACE_BLACK = list(
      description = "Black / African American race indicator at study entry",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "1 = Black / African American (n = 848, 18.5% of the pooled study 1-6 population per Sharma 2018 Table 1); 0 = not Black. Time-fixed per subject. Modulates the maximal fractional kout increase from lifestyle intervention (DSTIM): Black DSTIM = 10.9% vs the White/Asian reference 18.5% (Sharma 2018 Table 2). With RACE_BLACK = 0 and RACE_OTHER = 0 the model reproduces the White/Asian-pooled DSTIM = 18.5% (the paper merged White and Asian into a single reference group because of the small Asian subgroup, n = 49, 1.1%).",
      source_name = "Race == 'Black'"
    ),
    RACE_OTHER = list(
      description = "Race category 'Other' indicator (paper-pooled Pacific Islander / Native Hawaiian / Native American / Alaska Native / Other) at study entry",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "1 = subject in the paper-pooled 'others' race category (n = 152, 3.3% of the pooled study 1-6 population per Sharma 2018 Table 1; per Sharma 2018 Covariate Analysis the small Pacific Islander, Native Hawaiian, Native American / Alaska Native, and 'others' subgroups were lumped together into a single 'other' indicator because each was too small individually for meaningful effect estimation); 0 = not in the paper-pooled 'others' category. Time-fixed per subject. Modulates DSTIM: Other DSTIM = 9.73% vs the White/Asian reference 18.5% (Sharma 2018 Table 2). With RACE_BLACK = 0 and RACE_OTHER = 0 the model reproduces the White/Asian-pooled reference DSTIM.",
      source_name = "Race == 'Other' (paper-pooled)"
    ),
    DOSE_NAL_MGD = list(
      description = "Naltrexone daily dose at the current observation time (mg/day)",
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying covariate input carrying the daily naltrexone dose driving the dose-Hill term NAL/(ED50_NAL + NAL) in the combined Emax drug effect (Sharma 2018 Eq. 4). The full Contrave maintenance dose is 32 mg/day naltrexone (4 tablets/day of the 8 mg / 90 mg formulation); the placebo arm is 0 mg/day. Titration weeks 1-4 ramp from 8 -> 16 -> 24 -> 32 mg/day (one tablet added per week); supply the current week's daily dose. Set to 0 to disable the drug effect (e.g. for the placebo arm, or for use of this model purely for placebo+LSI BW-loss simulation).",
      source_name = "NAL (Sharma 2018 Eq. 4 -- daily naltrexone dose in mg)"
    ),
    DOSE_BUP_MGD = list(
      description = "Bupropion daily dose at the current observation time (mg/day)",
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying covariate input carrying the daily bupropion dose driving the dose-Hill term BUP/(ED50_BUP + BUP) in the combined Emax drug effect (Sharma 2018 Eq. 4). The full Contrave maintenance dose is 360 mg/day bupropion (4 tablets/day of the 8 mg / 90 mg formulation); the placebo arm is 0 mg/day. Titration weeks 1-4 ramp from 90 -> 180 -> 270 -> 360 mg/day (one tablet added per week); supply the current week's daily dose. Set to 0 to disable the drug effect.",
      source_name = "BUP (Sharma 2018 Eq. 4 -- daily bupropion dose in mg)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 4591,
    n_studies = 6,
    age_range = "adults; pooled median 46 y [SD 11.3]; per-study medians 44-55 y",
    weight_range = "pooled median 99 kg [SD 15.8]; per-study medians 94-104 kg",
    sex_female_pct = 82.6,
    race_ethnicity = c(White = 77.2, Asian = 1.1, Black = 18.5, Other = 3.3),
    disease_state = "Obese or overweight adults (BMI 27-45 kg/m^2; overweight 2.4%, class I obesity 38.0%, class II 36.2%, class III 23.3%); study 6 (NB-304, 10.9% of the analysis population) enrolled obese subjects with type-2 diabetes mellitus, the remaining 5 studies enrolled obese non-diabetic subjects with or without other comorbidities (dyslipidemia, controlled hypertension)",
    dose_range = "Contrave (naltrexone + bupropion) fixed-dose combination tablets 8 mg / 90 mg; titrated to a maintenance dose of 32 mg naltrexone + 360 mg bupropion per day (4 tablets/day); placebo and active arms",
    regions = "United States (six Contrave registration trials OT-101, NB-201, NB-301, NB-302, NB-303, NB-304)",
    notes = "Baseline characteristics per Sharma 2018 Table 1 (pooled across all six trials and stratified per trial). For model development only 5 trials (1, 2, 4, 5, 6) were used; trial NB-301 (study 3) was held out as the external evaluation dataset and is reported only as a validation cohort (Sharma 2018 Methods -- 'For model development only 5 studies were used, and BW data from study NB-301 were kept aside as external model evaluation data set'). 21,488 observed BW measurements over 20-65 weeks. All subjects received lifestyle intervention (hypocaloric diet -500 to -1500 kcal/day plus increased physical activity >= 30-min walk x 3/week). The original analysis was run in Monolix 4.3.3 (SAEM)."
  )

  ini({
    # ----- Structural parameters (Sharma 2018 Table 2, "DTPD Model" column, nondiabetic / White-Asian reference) -----

    # Indirect-response BW-loss rate constant (per week).
    lkout <- log(0.0543)
    label("First-order BW loss rate constant kout for nondiabetic reference (1/week)")  # Sharma 2018 Table 2 DTPD column, nondiabetic subjects: 0.0543 [RSE 3%]

    # LSI inverse-Bateman onset rate constant (per week).
    lkrel <- log(0.0344)
    label("LSI onset rate constant krel (1/week)")  # Sharma 2018 Table 2 DTPD column: 0.0344 [RSE 6%]

    # LSI inverse-Bateman loss rate constant (per week).
    lkde <- log(0.0792)
    label("LSI loss rate constant kde (1/week)")  # Sharma 2018 Table 2 DTPD column: 0.0792 [RSE 2%]

    # Maximal fractional increase in kout from lifestyle intervention, White / Asian reference.
    ldstim <- log(0.185)
    label("Maximal fractional kout increase from LSI, DSTIM, for White/Asian reference (fraction)")  # Sharma 2018 Table 2 DTPD column, White and Asian: 18.5% [RSE 3%]

    # Maximal drug-induced BW loss, nondiabetic reference.
    lemax <- log(4.69)
    label("Maximal drug-induced BW loss Emax for nondiabetic reference (kg)")  # Sharma 2018 Table 2 DTPD column, nondiabetic subjects: 4.69 [RSE 3%]

    # Median effective time for the time-driven half-maximal drug effect.
    let50 <- log(9.14)
    label("Time at half-maximal drug effect ET50 (week)")  # Sharma 2018 Table 2 DTPD column "T50 (week-1)": 9.14 [RSE 1%]
    # Units note: Sharma 2018 Table 2 column heading prints "T50 (week-1)" but ET50 enters Eq. (4) as the time-Hill denominator t + ET50 with t in weeks -- so ET50 has units of weeks, not 1/week. The "-1" in the table header is a typesetting artefact (the same "(week-1)" appears in the kout / krel / kde rows which are genuinely 1/week).

    # Median effective doses of naltrexone and bupropion (combined Emax dose-Hill).
    led50_nal <- log(54.6)
    label("Naltrexone half-maximal dose ED50_NAL (mg)")  # Sharma 2018 Table 2 DTPD column "ED50, NAL": 54.6 [RSE 12%]
    led50_bup <- log(645)
    label("Bupropion half-maximal dose ED50_BUP (mg)")   # Sharma 2018 Table 2 DTPD column "ED50, BUP": 645 [RSE 33%]

    # Baseline body weight at t = 0 for nondiabetic reference.
    lbw0 <- log(98.4)
    label("Baseline body weight BW0 for nondiabetic reference (kg)")  # Sharma 2018 Table 2 DTPD column, nondiabetic baseline BW: 98.4 [RSE 0%]

    # Disease-progression BW gain rate, nondiabetic. FIXED at the NHANES literature value (Sharma 2018 Methods).
    lkpro_kgy <- fixed(log(0.7))
    label("BW disease-progression rate kpro for nondiabetic (kg/year, FIXED at NHANES literature value)")
    # Sharma 2018 Methods "Population Dose- and Time-Dependent Pharmacodynamic Model": "kpro was fixed to 0.7 kg per year for subjects without T2DM but was estimated for T2DM subjects to evaluate the effect of diabetes mellitus on the disease progression."
    # Table 2 DTPD column lists kpro = 0.7 [FIXED] in kg/y for nondiabetic. Discussion text inconsistently writes "0.7 kg/week" but that is contradicted by Methods, Table 2, and biological plausibility (0.7 kg/week would be 36 kg/year, not the NHANES adult linear BW-gain rate of ~0.7 kg/year).
    # Stored in kg/year and converted to kg/week in model() for consistency with the per-week time axis.

    # ----- Covariate effects -----
    # Race effects on DSTIM (multiplicative fractional shifts; reference = White/Asian).
    e_race_black_dstim <- 0.109 / 0.185 - 1
    label("Multiplicative shift on DSTIM for Black race vs White/Asian reference (fraction)")
    # Sharma 2018 Table 2 DTPD column: Black DSTIM = 10.9% vs reference 18.5% -> 0.109/0.185 - 1 = -0.411 (-41.1% relative)
    e_race_other_dstim <- 0.0973 / 0.185 - 1
    label("Multiplicative shift on DSTIM for 'Other' race vs White/Asian reference (fraction)")
    # Sharma 2018 Table 2 DTPD column: Other DSTIM = 9.73% vs reference 18.5% -> 0.0973/0.185 - 1 = -0.474 (-47.4% relative)

    # Diabetes effects (multiplicative fractional shifts; reference = nondiabetic).
    e_t2dm_emax <- 3.74 / 4.69 - 1
    label("Multiplicative shift on Emax for T2DM (fraction)")
    # Sharma 2018 Table 2 DTPD column: diabetic Emax = 3.74 vs nondiabetic 4.69 -> -0.2026 (-20.3% relative)
    e_t2dm_kout <- 0.0381 / 0.0543 - 1
    label("Multiplicative shift on kout for T2DM (fraction)")
    # Sharma 2018 Table 2 DTPD column: diabetic kout = 0.0381 vs nondiabetic 0.0543 -> -0.2984 (-29.8% relative)
    e_t2dm_kpro <- 2.7 / 0.7 - 1
    label("Multiplicative shift on kpro for T2DM (fraction)")
    # Sharma 2018 Table 2 DTPD column: diabetic kpro = 2.7 kg/y [RSE 11%] vs FIXED nondiabetic 0.7 -> 2.857 (+286% relative). Estimated for the diabetic cohort only.
    e_t2dm_bw0 <- 103 / 98.4 - 1
    label("Multiplicative shift on baseline BW for T2DM (fraction)")
    # Sharma 2018 Table 2 DTPD column: diabetic baseline BW = 103 vs nondiabetic 98.4 -> +0.04675 (+4.7% relative)

    # ----- IIV (Sharma 2018 Table 2 DTPD IIV column; reported as CV%) -----
    # Variance on the log scale: omega^2 = log(CV^2 + 1).
    # Per Sharma 2018 Results 'Population Dose- and Time-Dependent Pharmacodynamic Model': "IIV was added exponentially on all final model parameter estimates except on kde and kpro because no IIV was tested statistically significant on these parameters." kde and kpro therefore have no eta.
    etalkout      ~ log(0.038^2 + 1)   # 3.8% CV  -> 0.001443  # Sharma 2018 Table 2 DTPD IIV kout: 3.8% [RSE 4%]
    etalkrel      ~ log(1.77^2 + 1)    # 177% CV  -> 1.41908   # Sharma 2018 Table 2 DTPD IIV krel: 177% [RSE 3%]
    etaldstim     ~ log(0.177^2 + 1)   # 17.7% CV -> 0.030840  # Sharma 2018 Table 2 DTPD IIV DSTIM: 17.7% [RSE 3%]
    etalemax      ~ log(0.526^2 + 1)   # 52.6% CV -> 0.244338  # Sharma 2018 Table 2 DTPD IIV Emax: 52.6% [RSE 13%]
    etalet50      ~ log(0.125^2 + 1)   # 12.5% CV -> 0.015511  # Sharma 2018 Table 2 DTPD IIV T50: 12.5% [RSE 6%]
    etaled50_nal  ~ log(2.07^2 + 1)    # 207% CV  -> 1.66495   # Sharma 2018 Table 2 DTPD IIV ED50_NAL: 207% [RSE 6%]
    etaled50_bup  ~ log(6.07^2 + 1)    # 607% CV  -> 3.63284   # Sharma 2018 Table 2 DTPD IIV ED50_BUP: 607% [RSE 6%]
    etalbw0       ~ log(0.157^2 + 1)   # 15.7% CV -> 0.024399  # Sharma 2018 Table 2 DTPD IIV baseline BW: 15.7% [RSE 1%]

    # ----- Residual error (additive on body weight in kg; Sharma 2018 Results "Residual unexplained variability in the model was best described using an additive error model") -----
    addSd <- 1.35
    label("Additive residual SD on body weight (kg)")  # Sharma 2018 Table 2 DTPD column "Additive error (kg)": 1.35 [RSE 1%]
  })

  model({
    # Sharma 2018 DTPD model: indirect-response body-weight model with disease progression (kpro),
    # inverse-Bateman LSI stimulation of kout, and a combined Emax dose- and time-dependent drug effect.
    # Time axis: weeks. State variable `central` holds BWprog,1 from the paper (eq. 3), and the
    # single-output observation `Cc` carries BWprog,2 = BWprog,1 - E_DTPD from eq. (4).
    # `Cc` here is a body-weight prediction (kg), not a drug concentration -- single-output naming
    # per nlmixr2lib convention; the canonical observation name is uniform across drug-PK and
    # endogenous / biomarker single-output models (cf. Bisaso 2014 plasma albumin).

    # ----- Unit-conversion constants -----
    weeks_per_year <- 365.25 / 7   # 52.1786 weeks/year

    # ----- Typical-value parameters with covariate modulation (linear-fraction multiplicative) -----
    kout_tv   <- exp(lkout)  * (1 + e_t2dm_kout * T2DM)
    krel_tv   <- exp(lkrel)
    kde_tv    <- exp(lkde)
    dstim_tv  <- exp(ldstim) * (1 + e_race_black_dstim * RACE_BLACK) *
                               (1 + e_race_other_dstim * RACE_OTHER)
    emax_tv   <- exp(lemax)  * (1 + e_t2dm_emax * T2DM)
    et50_tv   <- exp(let50)
    ed50nal_tv <- exp(led50_nal)
    ed50bup_tv <- exp(led50_bup)
    bw0_tv    <- exp(lbw0)   * (1 + e_t2dm_bw0 * T2DM)
    # kpro is stored in kg/year and converted to kg/week for use with the per-week time axis.
    kpro_kgwk_tv <- exp(lkpro_kgy) * (1 + e_t2dm_kpro * T2DM) / weeks_per_year

    # ----- Individual parameters (log-normal IIV; kde and kpro carry no eta per the source) -----
    kout    <- kout_tv    * exp(etalkout)
    krel    <- krel_tv    * exp(etalkrel)
    kde     <- kde_tv
    dstim   <- dstim_tv   * exp(etaldstim)
    emax    <- emax_tv    * exp(etalemax)
    et50    <- et50_tv    * exp(etalet50)
    ed50nal <- ed50nal_tv * exp(etaled50_nal)
    ed50bup <- ed50bup_tv * exp(etaled50_bup)
    bw0     <- bw0_tv     * exp(etalbw0)
    kpro_kgwk <- kpro_kgwk_tv

    # ----- Baseline body weight with disease progression (Sharma 2018 Eq. 1) -----
    # BW_DP(t) = BW_baseline + kpro * t        [kg]
    bwdp <- bw0 + kpro_kgwk * t

    # ----- Steady-state input to the indirect-response model (Sharma 2018 Eq. 2) -----
    # kin = kout * BW_DP(t)                    [kg/week]
    # Note: kin is time-varying because BW_DP is time-varying through kpro.
    kin <- kout * bwdp

    # ----- LSI inverse-Bateman stimulation term (Sharma 2018 Eq. 3, inner bracket) -----
    # LSI(t) = (kde / (kde - krel)) * (exp(-krel * t) - exp(-kde * t))
    # Equals 0 at t = 0 (no LSI effect at study entry), rises to a peak around t = ln(kde/krel) /
    # (kde - krel) (= ~18.6 weeks at typical values 0.0792 / 0.0344), then decays as the
    # diet-and-exercise effect dissipates.
    lsi <- (kde / (kde - krel)) * (exp(-krel * t) - exp(-kde * t))

    # ----- Drug effect (Sharma 2018 Eq. 4 inner) -----
    # E_DTPD = (Emax * t / (t + ET50)) * (NAL/(ED50_NAL + NAL) + BUP/(ED50_BUP + BUP))     [kg]
    # The time-Hill factor t/(t + ET50) = 0 at t = 0 (no instantaneous drug effect at study entry)
    # and saturates to 1 for t >> ET50, giving a gradual ramp-up of the drug effect over weeks.
    # DOSE_NAL_MGD and DOSE_BUP_MGD are time-varying covariate inputs carrying the current daily
    # naltrexone and bupropion doses in mg/day; set both to 0 to simulate the placebo arm.
    time_factor <- t / (t + et50)
    drug_term   <- DOSE_NAL_MGD / (ed50nal + DOSE_NAL_MGD) +
                   DOSE_BUP_MGD / (ed50bup + DOSE_BUP_MGD)
    e_dtpd      <- emax * time_factor * drug_term

    # ----- ODE for the LSI-modulated BW state (Sharma 2018 Eq. 3) -----
    # dBWprog,1/dt = kin - kout * BWprog,1 * (1 + DSTIM * LSI(t))    [kg/week]
    # The state `central` holds BWprog,1; LSI multiplies kout above its baseline value so the
    # net effect of LSI is to enhance BW loss transiently.
    d/dt(central) <- kin - kout * central * (1 + dstim * lsi)

    # Initial condition: BWprog,1(0) = baseline body weight (kpro contribution is 0 at t = 0)
    central(0) <- bw0

    # ----- Observation: BWprog,2 = BWprog,1 - E_DTPD (Sharma 2018 Eq. 4 outer) -----
    Cc <- central - e_dtpd
    Cc ~ add(addSd)
  })
}
