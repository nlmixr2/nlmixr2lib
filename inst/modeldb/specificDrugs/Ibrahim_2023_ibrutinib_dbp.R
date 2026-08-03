Ibrahim_2023_ibrutinib_dbp <- function() {
  description <- paste(
    "Ibrutinib-induced diastolic blood pressure (dBP) turnover model in patients with chronic lymphocytic leukemia,",
    "from the phase Ib/II PCYC-1102 study. A zero-order production / first-order loss turnover model preceded by a",
    "single transit compartment reproduces the slow onset of ibrutinib-associated hypertension; ibrutinib stimulates",
    "the zero-order dBP production rate through an Emax function of the daily AUC(0-24). The transit rate constant",
    "and the turnover rate constant are equal, both set to (n + 1) / MTT with n = 1 transit compartment. Baseline age",
    "lowers the dBP baseline. Structurally identical to the sister systolic model but fitted independently, with a",
    "smaller Emax (0.0694 vs 0.113), a lower AUC50 (63.1 vs 91.7 h*ng/mL), a roughly 2-fold slower onset",
    "(MTT 161 vs 79.9 days), no random effect on MTT, and the age effect placed on the baseline rather than on MTT.",
    "Ibrutinib enters only through the time-varying covariate AUC_IBRU; there are no drug-dosing events in this",
    "model. Sister model files from the same paper: modellib('Ibrahim_2023_ibrutinib_sbp'),",
    "modellib('Ibrahim_2023_ibrutinib_leukocyte_spd'), modellib('Ibrahim_2023_ibrutinib_competing_risk')."
  )
  reference <- paste(
    "Ibrahim EIK, Karlsson MO, Friberg LE.",
    "Assessment of ibrutinib scheduling on leukocyte, lymph node size and blood pressure dynamics in chronic",
    "lymphocytic leukemia through pharmacokinetic-pharmacodynamic models.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(9):1305-1318.",
    "doi:10.1002/psp4.13010.",
    "Open Access under CC BY-NC.",
    "Structural equations from main-text equations 2 and 3; parameter values from the authors' own nlmixr control",
    "stream (Supporting Information S4, run300) cross-checked against Table 2.",
    sep = " "
  )
  vignette <- "Ibrahim_2023_ibrutinib"
  units <- list(
    time          = "day",
    dosing        = "n/a (no drug-dosing events; ibrutinib exposure enters as the time-varying covariate AUC_IBRU)",
    concentration = "diastolic blood pressure in mmHg (not a drug concentration)"
  )

  covariateData <- list(
    AUC_IBRU = list(
      description        = "Daily 0-24 h area under the ibrutinib plasma concentration-time curve, AUC(0-24).",
      units              = "h*ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying: the value tracks the patient's current daily ibrutinib dose level and drops to 0 during",
        "treatment interruptions. Ibrutinib enters this model ONLY through this column -- the model contains no drug",
        "compartment and no dosing events. Ibrahim 2023 derived per-subject AUC(0-24) by integrating the individual",
        "post hoc profiles of the two-compartment ibrutinib population PK model of Marostica et al. (Cancer Chemother",
        "Pharmacol. 2015;75(1):111-121), which is NOT part of nlmixr2lib. Enters the dBP-production stimulation Emax",
        "function as Emax * AUC_IBRU / (AUC50 + AUC_IBRU) with AUC50 = 63.1 h*ng/mL (Ibrahim 2023 Table 2) -- about",
        "twice the corresponding IAUC50 of the pBtk efficacy model (34.1 h*ng/mL)."
      ),
      source_name        = "DAILYAUC"
    ),
    AGE = list(
      description        = "Baseline age.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject (baseline value). Enters the dBP baseline in power form",
        "dbpbase = exp(lrbase_dbp + etalrbase_dbp + e_age_base_dbp * log(AGE / 63)), equivalently",
        "dBPbaseline = 69.7 * (AGE / 63)^-0.204, with reference age 63 years (Ibrahim 2023 Table 2 footnote b:",
        "dBPbaseline = e^(LN(69.7) - 0.204 * LN(Age/63))). The population mean baseline age was 62.4 (SD 9.9) years",
        "(Supplementary Table S1). Note the contrast with the sister systolic model, where age acts on the mean",
        "transit time rather than on the baseline: 'The baseline dBP was inversely correlated with age, as previously",
        "described, and older patients had lower MTTsBP values' (Ibrahim 2023 Discussion)."
      ),
      source_name        = "LNAGE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 120L,
    n_studies      = 1L,
    age_range      = "mean 62.4 (SD 9.9) years",
    weight_range   = "mean 82.3 (SD 17) kg",
    sex_female_pct = 24.2,
    race_ethnicity = NULL,
    disease_state  = "Chronic lymphocytic leukemia (CLL); 20.8% treatment-naive, 79.2% relapsed/refractory",
    dose_range     = "ibrutinib 420 mg once daily (n = 94) or 840 mg once daily (n = 38) in PCYC-1102",
    regions        = "United States (PCYC-1102, phase Ib/II)",
    notes          = paste(
      "Baseline demographics from Ibrahim 2023 Supplementary Table S1. The blood-pressure dataset contained 2413",
      "paired sBP and dBP measurements (Ibrahim 2023 Appendix S1 section 1). At baseline 0.25% of patients had",
      "dBP >= 90 mmHg; after 2 years of the approved 420 mg/day schedule the model predicted 7.83%",
      "(Ibrahim 2023 Results 'Evaluation of the de-escalation dosing schedules'). The effect of antihypertensive",
      "drugs was tested as a time-varying binary covariate but was NOT retained in the final model",
      "(Ibrahim 2023 Patients and Methods 'PK-blood pressure model'). Estimation used FOCEI in nlmixr 2.0.6."
    )
  )

  covariatesDataExcluded <- list(
    CONMED_ANTIHYPERTENSIVE = list(
      description = "Antihypertensive co-medication indicator.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened but NOT retained in the final model. Ibrahim 2023 Patients and Methods 'PK-blood pressure model'",
        "tested the effect of antihypertensive drugs as a time-varying binary covariate, 'either as a step function",
        "or as a function in which the antihypertensive effect gradually evolved'; neither form survived the",
        "backward-elimination criterion and no coefficient is reported, so no value can be encoded."
      )
    ),
    SEXF = list(
      description = "Subject sex indicator: 1 = female, 0 = male.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened but NOT retained in the final dBP model. Ibrahim 2023 Patients and Methods 'Covariate analysis'",
        "evaluated baseline age and gender in the blood pressure and competing risk models; only the age effect on",
        "the dBP baseline was retained (Table 2). No coefficient is reported for sex."
      )
    )
  )

  ini({
    # ------------------------------------------------------------------
    # All log-scale THETA values are the final estimates taken verbatim
    # from the authors' own nlmixr control stream (Ibrahim 2023 Supporting
    # Information S4, run300), which carries more significant figures than
    # the rounded back-transformed values printed in Table 2. Every value
    # below back-transforms to its Table 2 row.
    # ------------------------------------------------------------------
    lrbase_dbp <- 4.24388453
    label("Log baseline diastolic blood pressure (mmHg)")  # S4 run300 tdbpbas = 4.24388453 -> 69.7; Table 2 dBPbaseline = 69.7 (95% CI 68.5-70.9)
    lmtt_dbp   <- 5.08188777
    label("Log mean transit time of the dBP transit/turnover chain (day)")  # S4 run300 tmtt = 5.08188777 -> 161; Table 2 MTTdBP = 161 (95% CI 106-244)
    lemax_dbp  <- -2.66779698
    label("Log maximum fractional ibrutinib stimulation of dBP production (unitless)")  # S4 run300 temax = -2.66779698 -> 0.0694; Table 2 Emax,dBP = 0.0694 (95% CI 0.0482-0.1)
    lauc50_dbp <- 4.14476079
    label("Log daily AUC(0-24) giving 50% of the maximum dBP stimulation (h*ng/mL)")  # S4 run300 tec50 = 4.14476079 -> 63.1; Table 2 AUC50,dBP = 63.1 (95% CI 9.75-408)

    e_age_base_dbp <- -0.20429000
    label("Power coefficient of baseline age (referenced to 63 years) on the dBP baseline")  # S4 run300 eff_age_bas = -0.20429000; Table 2 'Coefficient of baseline age effect on dBPbaseline' = -0.204 (95% CI -0.298 to -0.111)

    ntransit_dbp <- fixed(2)
    label("(n + 1) with n = 1 transit compartment; ktr = kout = ntransit_dbp / MTT (unitless)")  # S4 run300 nn = 2; Ibrahim 2023 main text 'ktr,iBP = (n + 1) / MTTiBP' with a single transit compartment

    # --- IIV ($OMEGA variances, S4 run300) ------------------------------
    # Table 2 reports these as CV% via CV = sqrt(exp(omega^2) - 1); the
    # control stream reports the variances directly, which is what ini()
    # needs, so the variances are used here. Unlike the sister systolic
    # model, run300 has NO random effect on MTT -- Table 2 correspondingly
    # reports no CV% for MTTdBP.
    etalrbase_dbp ~ 0.008102388  # S4 run300 eta.dbpbas; Table 2 dBPbaseline CV% = 9.02
    etalemax_dbp  ~ 0.7615876    # S4 run300 eta.emax;   Table 2 Emax,dBP CV% = 107

    # --- residual error --------------------------------------------------
    # Table 2 footnote d: 'Additive residual unexplained variability (RUV)
    # model was implemented on log transformed data'. The authors' code
    # writes `dbp_ipred = log(dbp + 1E-6); dbp_ipred ~ add(add.sd1)`.
    # Additive-on-log is exactly nlmixr2's `lnorm()` residual, so the output
    # stays in mmHg and carries a log-normal residual instead of being
    # emitted on the log scale.
    expSd_dbp <- 0.08883375
    label("Log-scale (log-normal) residual error SD for diastolic blood pressure (unitless on the log scale)")  # S4 run300 add.sd1 = 0.08883375; Table 2 'RUV dBP (%) = 8.9'
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual parameters. Baseline age enters the dBP BASELINE in
    #    power form (Ibrahim 2023 Table 2 footnote b; reference age 63
    #    years) -- contrast the sister systolic model, where the age effect
    #    sits on the mean transit time instead.
    # ------------------------------------------------------------------
    dbpbase   <- exp(lrbase_dbp + etalrbase_dbp + e_age_base_dbp * log(AGE / 63))
    mtt_dbp   <- exp(lmtt_dbp)
    emax_dbp  <- exp(lemax_dbp + etalemax_dbp)
    auc50_dbp <- exp(lauc50_dbp)

    # ------------------------------------------------------------------
    # 2. Derived rate constants. The turnover rate constant was assumed
    #    equal to the transit rate constant (Ibrahim 2023 main text,
    #    'PK-blood pressure model' paragraph).
    # ------------------------------------------------------------------
    ktr_dbp  <- ntransit_dbp / mtt_dbp
    kout_dbp <- ktr_dbp
    kin_dbp  <- dbpbase * kout_dbp

    # Emax stimulation of the zero-order dBP production rate.
    effauc <- emax_dbp * AUC_IBRU / (auc50_dbp + AUC_IBRU)

    # ------------------------------------------------------------------
    # 3. ODE system (Ibrahim 2023 equations 2 and 3).
    # ------------------------------------------------------------------
    d/dt(transit1) <- kin_dbp * (1 + effauc) - ktr_dbp * transit1
    transit1(0)    <- dbpbase

    d/dt(dbp) <- ktr_dbp * transit1 - kout_dbp * dbp
    dbp(0)    <- dbpbase

    # ------------------------------------------------------------------
    # 4. Observation. `dbp` is the registered canonical diastolic-blood-
    #    pressure PD state (Hansson 2013 sunitinib) and is itself the
    #    observation variable, mirroring the systolic sibling `sbp`.
    # ------------------------------------------------------------------
    dbp ~ lnorm(expSd_dbp)
  })
}
