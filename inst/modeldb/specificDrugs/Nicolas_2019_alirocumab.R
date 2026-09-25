Nicolas_2019_alirocumab <- function() {
  description <- "Sequential population PK/PD model for alirocumab and serum LDL cholesterol in healthy volunteers and adults with hypercholesterolemia (Nicolas 2019, Part II). A type IV (stimulation-of-loss) indirect response model with a Hill coefficient links total alirocumab concentration to LDL-C elimination; the PK layer is the Michaelis-Menten target-mediated approximation of Martinez 2019 (Part I), carried fixed as the concentration driver. Ten covariates act on the four PD parameters (Emax, EC50, Hill coefficient and Kout)."
  reference <- "Nicolas X, Djebli N, Rauch C, Brunet A, Hurbin F, Martinez JM, Fabre D. Population Pharmacokinetic/Pharmacodynamic Analysis of Alirocumab in Healthy Volunteers or Hypercholesterolemic Subjects Using an Indirect Response Model to Predict Low-Density Lipoprotein Cholesterol Lowering: Support for a Biologics License Application Submission: Part II. Clin Pharmacokinet. 2019;58(1):115-130. doi:10.1007/s40262-018-0670-5"
  vignette <- "Nicolas_2019_alirocumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. The three PK states are inherited from Martinez 2019
  # (Part I); `ldl` is the Part II turnover pool and is carried as a
  # CONCENTRATION rather than an amount, because the paper parameterises the
  # indirect response model directly on the measured serum LDL-C concentration
  # (mg/dL) with no volume term anywhere in the model.
  compartmentData <- list(
    depot = list(analyte = "alirocumab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "alirocumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "alirocumab", units = "mg", specimen = "serum", verified = TRUE),
    ldl = list(analyte = "low-density lipoprotein cholesterol", units = "mg/dL", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    LDLC = list(
      description = "Individual pre-treatment (baseline) serum LDL cholesterol concentration",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Dual role, as anticipated by the LDLC register entry. (1) It is the initial condition of the `ldl` turnover state. (2) It anchors the zero-order production rate through the drug-free steady state of the type IV indirect response model, kin = kout * LDLC. Nicolas 2019 tabulates Kout, EC50, Emax and the Hill coefficient but no kin and no typical baseline (Table 2), and 'baseline LDL-C levels' appear in the Sect. 2.5 covariate-screening list without being retained, so the subject's own observed baseline is the only self-consistent source of kin. Because the model is linear in the state, dividing through by the baseline removes it: the PERCENT change from baseline predicted by this model is exactly independent of LDLC, which is why Nicolas 2019 reports its derived endpoints (DLDL-Cmax, DLDL-Ctrough) only as percentages. The vignette asserts that invariance.",
      source_name = "BSLDLC"
    ),
    PCSK9 = list(
      description = "Free (unbound) serum PCSK9 concentration at baseline",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Nicolas 2019 calls this FBSPCSK9 and computes it as the mean of all measurements taken before the first alirocumab dose, including screening (Sect. 2.3). Time-fixed per subject. Two effects, both additive linear deviations from the 265 ng/mL median of the pooled data set (Table 2 footnotes c and d): on Emax, + 0.00156 per ng/mL, applied OUTSIDE the multiplicative sex / age / weight chain; on the Hill coefficient, + 0.00340 per ng/mL. Observed 5th-50th-95th percentiles 126 / 265 / 501 ng/mL (Sect. 3.5). Distinct from TPCSK9_BASE (total, i.e. free plus drug-bound) - the pooled baseline means differ roughly two-fold, 282 vs 675 ng/mL (ESM Table 2).",
      source_name = "FBSPCSK9"
    ),
    TPCSK9 = list(
      description = "Total (free plus alirocumab-bound) serum PCSK9 concentration, time-varying",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Nicolas 2019 calls this TPCSK9. Additive linear-deviation effect on Emax, + 0.000331 per ng/mL relative to the 3340 ng/mL median of the pooled data set, applied INSIDE the multiplicative sex / age / weight chain (Table 2 footnote c). This is a genuinely time-varying, treatment-elevated quantity: anti-PCSK9 antibody binding slows clearance of the PCSK9 pool, so the pooled on-treatment median (3340 ng/mL) is roughly five-fold the pooled pre-treatment median of TPCSK9_BASE (644 ng/mL). Observed 5th-95th percentiles 491-6340 ng/mL (Sect. 3.5). The model does NOT generate this trajectory - it must be supplied as a data column, from measurement or from a mechanistic total-PCSK9 model such as Djebli_2017_alirocumab.R.",
      source_name = "TPCSK9"
    ),
    TPCSK9_BASE = list(
      description = "Total (free plus drug-bound) serum PCSK9 concentration at baseline",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Nicolas 2019 calls this TBSPCSK9. Time-fixed per subject, computed as the mean of all pre-first-dose measurements including screening (Sect. 2.3). Additive linear effect on EC50, + 0.00219 mg/L per ng/mL, entered UNCENTRED - Table 2 footnote b writes EC50 = (1.44 + 0.00219 * TBSPCSK9) * 1.21^HDSTATIN, so the 1.44 mg/L typical value is the EC50 extrapolated to zero total baseline PCSK9 and is never realised in the observed range. Observed 5th-50th-95th percentiles 355 / 644 / 1130 ng/mL (Sect. 3.5 and Sect. 4), giving EC50 = 2.22-3.91 mg/L without a high-dose statin.",
      source_name = "TBSPCSK9"
    ),
    SEXF = list(
      description = "Biological sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Nicolas 2019 codes SEX = 0 for male and 1 for female (Table 2 footnote c), which is the canonical SEXF orientation - no transformation is needed. Multiplicative effect on Emax of 0.703^SEXF, i.e. Emax is 29.7% lower in women than in men at otherwise identical covariates (Sect. 3.5). The factor multiplies only the TPCSK9-adjusted Emax core, not the free-baseline-PCSK9 or statin additive terms; ESM Table 3 pins that placement (see the model-file comment above emax_tv). Male 62.3% of the pooled data set (ESM Table 2).",
      source_name = "SEX"
    ),
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on the Emax core, (AGE/60)^0.415 (Table 2 footnote c); 60 years is the median of the pooled data set. Observed 5th-95th percentiles 37-75 years, spanning a 26.4-34.1% increase in Emax (Sect. 3.5). Pooled mean (SD) 58.2 (11.7) years (ESM Table 2).",
      source_name = "AGE"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on the Emax core, (WEIGHT/82.5)^0.313 (Table 2 footnote c); 82.5 kg is the median of the pooled Nicolas 2019 data set. Note this differs from the 82.9 kg reference of the Part I PK model (Martinez_2019_alirocumab.R) - the two analyses pooled the same 13 studies but different record sets (13,717 alirocumab concentrations vs 14,346 LDL-C values), so each reports its own median. WT is NOT a covariate on the inherited PK layer's clearance here, because those PK parameters enter this model fixed at their Part I typical values; a user simulating exposure-matched cohorts should drive the Part I model directly. Observed 5th-95th percentiles 58.1-119 kg (Sect. 3.5); pooled mean (SD) 85.0 (18.4) kg (ESM Table 2).",
      source_name = "WEIGHT"
    ),
    CONMED_STATIN = list(
      description = "Concomitant statin (HMG-CoA reductase inhibitor) coadministration, any agent at any dose",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (alirocumab given alone)",
      notes = "Nicolas 2019 codes STATIN = 0 if alirocumab was given alone and 1 if coadministered with a statin (Table 2 footnote c), with no dose or molecule restriction - a BROADER definition than the Martinez 2019 Part I CONMED_STATIN, which excluded rosuvastatin >= 20 mg/day and atorvastatin >= 40 mg/day. Additive effect on Emax of +0.408, applied outside the multiplicative chain; +16.8% at median covariates in men (Sect. 3.5). 2588/2799 (92.5%) of the pooled data set were on a statin (ESM Table 2). Compatible with, and not redundant against, CONMED_STATIN_HI: a high-intensity subject carries both = 1.",
      source_name = "STATIN"
    ),
    CONMED_STATIN_HI = list(
      description = "Concomitant high-intensity statin coadministration",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no statin, or a low- or moderate-intensity statin regimen)",
      notes = "Nicolas 2019 calls this HDSTATIN and defines it as rosuvastatin >= 20 mg/day or atorvastatin >= 40 mg/day (Table 2 footnote b), which is exactly the 2018 ACC/AHA high-intensity stratum used by the CONMED_STATIN_LI / _MI / _HI register entry, so no per-model departure from that table has to be recorded. Multiplicative effect on EC50 of 1.21^CONMED_STATIN_HI, i.e. +20.6% (Sect. 3.5). Reference category 0 pools no-statin with the low- and moderate-intensity strata, because the paper fits only the high-dose contrast; CONMED_STATIN_LI and CONMED_STATIN_MI are therefore not referenced by this model. 1305/2799 (46.6%) of the pooled data set (ESM Table 2).",
      source_name = "HDSTATIN"
    ),
    DIS_HEALTHY = list(
      description = "Healthy-volunteer cohort indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (patient with familial or non-familial hypercholesterolemia)",
      notes = "INVERTED alias. Nicolas 2019 codes DISST = 0 for a healthy volunteer and 1 for a patient (Table 2 footnote a); the canonical DIS_HEALTHY runs the other way, so DIS_HEALTHY = 1 - DISST. The published two-level expression Kout = 0.00395*(1-DISST) + 0.00997*DISST is therefore written here as Kout = 0.00395*DIS_HEALTHY + 0.00997*(1-DIS_HEALTHY): the 0.00395 /h 'typical value' row of Table 2 is the HEALTHY-VOLUNTEER Kout, and the 0.00997 /h 'effect of DISST' row is the PATIENT Kout in full, not an increment on it. Patients turn LDL-C over 2.52-fold faster (Sect. 3.5). Only 150/2799 (5.36%) of the pooled data set were healthy volunteers, and none of them received a statin while 2588/2649 (97.7%) of the patients did, so disease state and statin use are strongly collinear here (Sect. 4).",
      source_name = "DISST"
    )
  )

  # Covariates screened in Nicolas 2019 Sect. 2.5 (16 tested in total) but NOT
  # retained in the final model, and therefore not referenced in model().
  covariatesDataExcluded <- list(
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened as a demographic covariate (Sect. 2.5) and used to stratify the derived-endpoint summaries (Fig. 4), but not retained on any PD parameter. Pooled mean (SD) 29.5 (5.42) kg/m^2 (ESM Table 2)."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened (Sect. 2.5), not retained. Pooled mean (SD) 41.7 (3.34) g/L (ESM Table 2)."
    ),
    CONMED_EZE = list(
      description = "Concomitant ezetimibe coadministration",
      units = "(binary)",
      type = "binary",
      notes = "Screened as one of the 'relevant background therapies' (Sect. 2.5), not retained. 457/2799 (16.3%) of the pooled data set (ESM Table 2)."
    ),
    CONMED_FIBRATE = list(
      description = "Concomitant fibrate coadministration",
      units = "(binary)",
      type = "binary",
      notes = "Screened as one of the 'relevant background therapies' (Sect. 2.5), not retained. 130/2799 (4.64%) of the pooled data set (ESM Table 2)."
    ),
    FPCSK9 = list(
      description = "Free (unbound) serum PCSK9 concentration, time-varying",
      units = "ng/mL",
      type = "continuous",
      notes = "Nicolas 2019 screened free PCSK9 'both at baseline and time varying' (Sect. 2.5). Only the BASELINE free value was retained, and it is carried in covariateData as PCSK9; the time-varying free value was not retained on any PD parameter. It IS a covariate of the Part I PK model (Martinez_2019_alirocumab.R, additive effect on Km), but those PK parameters enter this model fixed at their Part I typical values, so nothing here references it."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 2799L,
    n_observations = 14346L,
    n_studies = 13L,
    phases = "Phase I (5 studies), II (4 studies), and III (4 studies: ODYSSEY MONO, COMBO II, FH I and LONG TERM)",
    age_mean_sd = "58.2 (11.7) years; median 60, 5th-95th percentiles 37-75",
    weight_mean_sd = "85.0 (18.4) kg; median 82.5, 5th-95th percentiles 58.1-119",
    bmi_mean_sd = "29.5 (5.42) kg/m^2",
    sex_female_pct = 37.7,
    disease_state = "2649 (94.6%) patients with heterozygous familial or non-familial hypercholesterolemia, several with established coronary heart disease or a risk equivalent, not adequately controlled on a maximally tolerated statin dose or statin-intolerant; 150 (5.36%) healthy volunteers with elevated LDL-C.",
    dose_range = "IV 0.3-12 mg/kg single dose (one phase I study, n = 30). SC 50-300 mg, single or repeated Q2W or Q4W for up to 2 years. The phase III regimens that produced the derived endpoints are 75 mg Q2W and 150 mg Q2W.",
    regions = "Multi-regional pool of 13 Sanofi / Regeneron trials, including two Japanese cohorts (NCT01448317 phase I, NCT01812707 phase II).",
    co_medication = "Any statin 2588 (92.5%); high-intensity statin 1305 (46.6%); low-intensity statin 1283 (45.8%); ezetimibe 457 (16.3%); any fibrate 130 (4.64%); alirocumab monotherapy 161 (5.75%).",
    pcsk9_baseline = "Free baseline PCSK9 mean (SD) 282 (119) ng/mL, median 265, 5th-95th percentiles 126-501. Total baseline PCSK9 mean (SD) 675 (247) ng/mL, median 644, 5th-95th percentiles 355-1130. Time-varying (on-treatment) total PCSK9 median 3340 ng/mL, 5th-95th percentiles 491-6340.",
    notes = "Baseline characteristics are from ESM Table 2; the covariate percentiles used as model reference values are from Table 2 footnotes a-d and Sect. 3.5. LDL-C was computed by the Friedewald formula throughout and is invalid where fasting triglycerides exceeded 400 mg/dL (ESM Supplementary Methods). Records were excluded for missing LDL-C, missing alirocumab concentration, or a concentration below the limit of quantification; missing covariates were imputed by last observation carried forward (Sect. 2.3)."
  )

  ini({
    # ---------------------------------------------------------------------
    # PK LAYER - INHERITED, NOT ESTIMATED HERE.
    #
    # Nicolas 2019 is the second step of a sequential (two-step) analysis: the
    # individual PK parameters of the Part I Michaelis-Menten target-mediated
    # population PK model were used to predict total alirocumab concentrations,
    # which then entered this PD model as a fixed driver (Sect. 2.4). Every
    # value below is therefore a Part I estimate carried unchanged and is
    # wrapped in fixed(); the source is Martinez 2019 Table 2, transcribed in
    # Martinez_2019_alirocumab.R, and NOT Nicolas 2019.
    #
    # Reference covariates for the PK layer: WT = 82.9 kg, AGE = 60 years,
    # CONMED_STATIN = 0, FPCSK9 = 72.9 ng/mL. Because the covariate effects are
    # NOT re-estimated here, the PK layer is carried at those Part I reference
    # values; see the model file's covariateData[[WT]] notes.
    #
    # Rates reported in the paper per hour are converted to per day (x24) for
    # the units$time = "day" convention shared with the Part I model file.
    lka <- fixed(log(7.68e-3 * 24))
    label("Part I: first-order SC absorption rate Ka (1/day)")
    lcl <- fixed(log(0.0124 * 24))
    label("Part I: linear clearance CLL at reference covariates (L/day)")
    lvc <- fixed(log(3.19))
    label("Part I: central volume of distribution V2 (L)")
    lvp <- fixed(log(2.79))
    label("Part I: peripheral volume of distribution V3 at reference age (L)")
    lq <- fixed(log(0.0185 * 24))
    label("Part I: intercompartmental clearance Q (L/day)")
    lvmax <- fixed(log(0.183 * 24))
    label("Part I: maximum Michaelis-Menten elimination rate Vmax (mg/day)")
    lkm <- fixed(log(7.73))
    label("Part I: Michaelis-Menten constant Km at reference free PCSK9 (mg/L)")
    ltlag <- fixed(log(0.641 / 24))
    label("Part I: SC absorption lag time (day)")
    logitfdepot <- fixed(log(0.862 / (1 - 0.862)))
    label("Part I: logit of SC bioavailability F (unitless; F_pop = 0.862)")

    # ---------------------------------------------------------------------
    # PD LAYER - Nicolas 2019 Table 2, 'Final model with covariates' column.
    #
    # Type IV indirect response model (stimulation of the loss of response),
    # Sect. 2.4 and Fig. 1. Kout is reported per hour and converted to per day
    # (x24). EC50 is in mg/L, matching the mg/L serum concentration produced by
    # the PK layer. Emax and the Hill coefficient are dimensionless.
    lkout <- log(0.00395 * 24)
    label("Typical Kout in healthy volunteers, DIS_HEALTHY = 1 (1/day)")
    # Nicolas 2019 Table 2 'Typical value of Kout' 0.00395 /h; footnote a shows
    # this is the DISST = 0 (healthy volunteer) level, not a pooled typical
    # value, because Kout is coded as a two-level switch rather than as a
    # reference value with a multiplicative shift.
    lec50 <- log(1.44)
    label("Typical EC50 extrapolated to TPCSK9_BASE = 0 and no high-intensity statin (mg/L)")
    # Nicolas 2019 Table 2 'Typical value of EC50' 1.44 mg/L; footnote b enters
    # TBSPCSK9 UNCENTRED, so 1.44 is an intercept outside the observed range.
    lemax <- log(2.43)
    label("Typical Emax at median TPCSK9, male, median age and weight (unitless)")
    # Nicolas 2019 Table 2 'Typical value of Emax' 2.43. The tabulated 95% CI
    # for this row, '1.63-1.92', is a duplicate of the Hill-coefficient row's
    # CI and is not used here; the bootstrap 95% CI for Emax is 2.24-2.91.
    lhill <- log(1.78)
    label("Typical Hill coefficient at median free baseline PCSK9 (unitless)")
    # Nicolas 2019 Table 2 'Typical value of c' 1.78, footnote d.

    # Covariate effects on the PD parameters - Nicolas 2019 Table 2 and the
    # four displayed equations of Sect. 3.2 (reproduced verbatim in the Table 2
    # footnotes a-d).
    e_dis_healthy_kout <- 0.00997 * 24
    label("Kout in patients, DIS_HEALTHY = 0 (1/day)")
    # Nicolas 2019 Table 2 'Effect of DISST on Kout' 0.00997 /h. Footnote a:
    # Kout = 0.00395*(1-DISST) + 0.00997*DISST, so this is the patient Kout in
    # full, NOT an increment added to the healthy-volunteer value.
    e_tpcsk9_base_ec50 <- 0.00219
    label("Additive slope of total baseline PCSK9 on EC50 (mg/L per ng/mL, uncentred)")
    # Nicolas 2019 Table 2 'Effect of TBSPCSK9 on EC50'.
    e_conmed_statin_hi_ec50 <- 1.21
    label("Multiplicative factor on EC50 for a high-intensity statin (unitless)")
    # Nicolas 2019 Table 2 'Effect of HDSTATIN on EC50'; footnote b applies it
    # as 1.21^HDSTATIN.
    e_tpcsk9_emax <- 0.000331
    label("Additive slope of (TPCSK9 - 3340 ng/mL) on the Emax core (per ng/mL)")
    # Nicolas 2019 Table 2 'Effect of TPCSK9 on Emax' 0.000331. The Sect. 3.2
    # displayed equation prints '0.0003331'; that extra digit is a typesetting
    # error - Table 2, the Table 2 footnote c, the bootstrap column
    # (0.000328, 95% CI 0.000266-0.000418) and every Emax value in ESM Table 3
    # all reproduce only with 0.000331.
    e_sexf_emax <- 0.703
    label("Multiplicative factor on the Emax core for female sex (unitless)")
    # Nicolas 2019 Table 2 'Effect of SEX on Emax'; footnote c applies it as
    # 0.703^SEX with SEX = 1 for female.
    e_age_emax <- 0.415
    label("Power exponent of AGE/60 years on the Emax core (unitless)")
    # Nicolas 2019 Table 2 'Effect of AGE on Emax'.
    e_wt_emax <- 0.313
    label("Power exponent of WT/82.5 kg on the Emax core (unitless)")
    # Nicolas 2019 Table 2 'Effect of WEIGHT on Emax'.
    e_pcsk9_emax <- 0.00156
    label("Additive slope of (free baseline PCSK9 - 265 ng/mL) on Emax (per ng/mL)")
    # Nicolas 2019 Table 2 'Effect of FBSPCSK9 on Emax'.
    e_conmed_statin_emax <- 0.408
    label("Additive effect of concomitant statin on Emax (unitless)")
    # Nicolas 2019 Table 2 'Effect of STATIN on Emax'.
    e_pcsk9_hill <- 0.00340
    label("Additive slope of (free baseline PCSK9 - 265 ng/mL) on the Hill coefficient (per ng/mL)")
    # Nicolas 2019 Table 2 'Effect of FBSPCSK9 on c'; footnote d.

    # ---------------------------------------------------------------------
    # Inter-individual variability.
    #
    # PK layer: inherited from Martinez 2019 Table 2 and held fixed, for the
    # same reason as the PK thetas. The V3 / Km block was estimated with a
    # correlation of -0.793, so cov = -0.793 * sqrt(0.0735 * 0.298) = -0.11738.
    etalcl ~ fixed(0.232)
    etalvc ~ fixed(0.589)
    etalvp + etalkm ~ fixed(c(0.0735, -0.11738, 0.298))
    etalogitfdepot ~ fixed(1.060)

    # PD layer: exponential inter-individual variability on all four PD
    # parameters (Sect. 3.1). The Table 2 'Estimate (CV %)' column carries the
    # NONMEM variance with the percentage in parentheses equal to
    # 100 * sqrt(variance), not 100 * sqrt(exp(variance) - 1): 0.113 -> 33.6,
    # 0.123 -> 35.1, 0.420 -> 64.8, 0.296 -> 54.4, which reproduces the
    # tabulated 33.7 / 35.1 / 54.4 and the 64.8% quoted for Emax in Sect. 3.4.
    # (The Table 2 Emax IIV cell reads '65.8'; Sect. 3.4 reads 64.8, and 64.8
    # is what sqrt(0.420) gives.) The variances are transcribed directly, so
    # the discrepancy does not propagate.
    etalkout ~ 0.113
    etalec50 ~ 0.123
    etalemax ~ 0.420
    etalhill ~ 0.296

    # Residual variability on LDL-C: combined additive plus proportional
    # (Sect. 3.1, Table 2 'Final model with covariates').
    addSd_ldl <- 5.21
    label("LDL-C additive residual error (mg/dL)")
    propSd_ldl <- 0.224
    label("LDL-C proportional residual error (fraction)")
  })
  model({
    # -------------------------------------------------------------------
    # 1. PK layer (Part I, Martinez 2019), carried fixed. The covariate terms
    #    of the Part I model are held at their reference values, so the typical
    #    values below are the Part I typical values; only the Part I IIV is
    #    retained, to give the PD layer a realistic exposure spread.
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q <- exp(lq)
    ka <- exp(lka)
    vmax <- exp(lvmax)
    km <- exp(lkm + etalkm)
    fdepot <- expit(logitfdepot + etalogitfdepot)
    lag <- exp(ltlag)

    Cc <- central / vc

    # -------------------------------------------------------------------
    # 2. Individual PD parameters (Nicolas 2019 Sect. 3.2 displayed equations,
    #    identical to the Table 2 footnotes a-d), with exponential IIV.

    # Kout: two-level switch on disease state. Written with DIS_HEALTHY, the
    # canonical inverse of the paper's DISST (see covariateData).
    kout_tv <- exp(lkout) * DIS_HEALTHY + e_dis_healthy_kout * (1 - DIS_HEALTHY)
    kout <- kout_tv * exp(etalkout)

    # EC50: TPCSK9_BASE enters UNCENTRED and additively, then a multiplicative
    # high-intensity-statin factor.
    ec50_tv <- (exp(lec50) + e_tpcsk9_base_ec50 * TPCSK9_BASE) *
      e_conmed_statin_hi_ec50^CONMED_STATIN_HI
    ec50 <- ec50_tv * exp(etalec50)

    # Emax. The parenthesisation below is the paper's as printed, and is NOT
    # the only reading its narrative supports - it is pinned by ESM Table 3,
    # which tabulates Emax for all four sex x statin cells at three TPCSK9
    # levels. All twelve of those cells reproduce to three significant figures
    # only when the free-baseline-PCSK9 and statin terms sit OUTSIDE the sex
    # factor: e.g. the sex = 1, statin = 1, median-TPCSK9 cell is 2.12, which
    # is 2.43*0.703 + 0.408 = 2.116 and not (2.43 + 0.408)*0.703 = 1.995.
    # (Sect. 3.5's free-baseline-PCSK9 Emax range, '1.56-2.62' and '1.97-3.21',
    # is the one place in the paper that does NOT reconcile with this form; it
    # matches a sex-outermost reading instead. ESM Table 3 and the two printed
    # equations outvote it - see the vignette's Errata section.)
    emax_tv <- (exp(lemax) + e_tpcsk9_emax * (TPCSK9 - 3340)) *
      e_sexf_emax^SEXF *
      (AGE / 60)^e_age_emax *
      (WT / 82.5)^e_wt_emax +
      e_pcsk9_emax * (PCSK9 - 265) +
      e_conmed_statin_emax * CONMED_STATIN
    emax <- emax_tv * exp(etalemax)

    # Hill coefficient (the paper's gamma).
    hill_tv <- exp(lhill) + e_pcsk9_hill * (PCSK9 - 265)
    hill <- hill_tv * exp(etalhill)

    # -------------------------------------------------------------------
    # 3. Zero-order LDL-C production. Nicolas 2019 tabulates no kin and no
    #    typical baseline; the drug-free steady state of the type IV indirect
    #    response model, LDLC = kin / kout, supplies it from the subject's own
    #    observed pre-treatment LDL-C.
    kin <- kout * LDLC

    # -------------------------------------------------------------------
    # 4. ODEs. The three PK states are the Part I structure: two-compartment
    #    disposition with first-order lagged SC absorption and parallel linear
    #    plus Michaelis-Menten elimination from the central compartment.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot -
      (cl / vc) * central -
      vmax * Cc / (km + Cc) -
      (q / vc) * central +
      (q / vp) * peripheral1
    d/dt(peripheral1) <- (q / vc) * central - (q / vp) * peripheral1

    # LDL-C turnover: type IV indirect response, drug stimulating the loss.
    d/dt(ldl) <- kin - kout * (1 + emax * Cc^hill / (ec50^hill + Cc^hill)) * ldl
    ldl(0) <- LDLC

    # SC bioavailability and absorption lag apply to the depot; an IV dose
    # bypasses the depot via cmt = central on the dose record.
    f(depot) <- fdepot
    alag(depot) <- lag

    # Only LDL-C carries a residual error. Alirocumab concentration enters the
    # Part II analysis as an individual PREDICTION from Part I (Sect. 2.4), so
    # it is noise-free by construction here and Cc is an algebraic observable
    # rather than an endpoint.
    ldl ~ add(addSd_ldl) + prop(propSd_ldl)
  })
}
