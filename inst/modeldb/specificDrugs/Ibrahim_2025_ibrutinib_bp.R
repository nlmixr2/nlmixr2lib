Ibrahim_2025_ibrutinib_bp <- function() {
  description <- paste(
    "Ibrutinib-induced blood pressure model in patients with chronic lymphocytic leukemia, describing systolic (sBP)",
    "and diastolic (dBP) blood pressure SIMULTANEOUSLY, pooled across the phase Ib/II PCYC-1102 and phase III",
    "PCYC-1115 studies. Each endpoint is a zero-order production / first-order loss turnover compartment preceded by a",
    "single transit compartment, reproducing the slow onset of ibrutinib-associated hypertension; ibrutinib stimulates",
    "the zero-order production rate through an Emax function of the daily AUC(0-24). The transit rate constant and the",
    "turnover rate constant are equal, both set to (n + 1) / MTT with n = 1 transit compartment. Fitting sBP and dBP",
    "together let the authors estimate a SHARED maximum stimulatory effect (Emax,BP) and a SHARED potency (AUC50,BP)",
    "across the two endpoints, with a single shared random effect on Emax, while the mean transit times stay",
    "endpoint-specific and reproduce the observed delay of the diastolic response relative to the systolic one",
    "(126 vs 53.6 days). Baseline age shortens the systolic mean transit time. Ibrutinib enters only through the",
    "time-varying covariate AUC_IBRU; there are no drug-dosing events in this model. Sister model file from the same",
    "paper: modellib('Ibrahim_2025_ibrutinib_cll'). Predecessor framework, which fitted the two endpoints in separate",
    "runs: modellib('Ibrahim_2023_ibrutinib_sbp'), modellib('Ibrahim_2023_ibrutinib_dbp')."
  )
  reference <- paste(
    "Ibrahim EIK, Friberg LE.",
    "Optimizing ibrutinib posology in chronic lymphocytic leukemia using a semi-mechanistic pharmacometric framework.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(12):2186-2198.",
    "doi:10.1002/psp4.70124.",
    "Open Access under CC BY-NC.",
    "Structural equations transcribed from the authors' own RxODE control stream",
    "(Data S1, PSP-2025-0220-s02.docx, `run6023_bp`); parameter values from Table 2.",
    "Model structure inherited from Ibrahim EIK, Karlsson MO, Friberg LE.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(9):1305-1318. doi:10.1002/psp4.13010;",
    "see modellib('Ibrahim_2023_ibrutinib_sbp').",
    sep = " "
  )
  vignette <- "Ibrahim_2025_ibrutinib"
  # `sbp` and `dbp` are canonical PD states. The two upstream transit
  # compartments are declared here because this is the first library model to
  # carry TWO parallel transit chains, one per blood-pressure endpoint: the
  # canonical bare `transit1` is a single chain and cannot name both. The
  # names follow the paired-output suffix idiom already used for
  # `circ_<celltype>` (compartment-names.md), applied to the registered
  # `transit<n>` chain. Declared as paper-specific rather than minted as a new
  # canonical so the operator can decide in review whether a paired-output
  # transit suffix deserves registration; see the vignette Errata.
  paper_specific_compartments <- c("sbp_transit1", "dbp_transit1")
  units <- list(
    time          = "day",
    dosing        = "n/a (no drug-dosing events; ibrutinib exposure enters as the time-varying covariate AUC_IBRU)",
    concentration = "systolic and diastolic blood pressure in mmHg (neither output is a drug concentration)"
  )

  compartmentData <- list(
    sbp_transit1 = list(analyte = "systolic blood pressure", units = "mmHg", specimen = "not applicable", verified = TRUE),
    sbp          = list(analyte = "systolic blood pressure", units = "mmHg", specimen = "not applicable", verified = TRUE),
    dbp_transit1 = list(analyte = "diastolic blood pressure", units = "mmHg", specimen = "not applicable", verified = TRUE),
    dbp          = list(analyte = "diastolic blood pressure", units = "mmHg", specimen = "not applicable", verified = TRUE)
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
        "compartment and no dosing events. Ibrahim 2025 computed daily AUC(0-24) from the individual ibrutinib plasma",
        "concentrations using the two-compartment population PK model of Marostica et al. (Cancer Chemother Pharmacol.",
        "2015;75(1):111-121), which is NOT part of nlmixr2lib. Enters the blood-pressure production stimulation as",
        "Emax,BP * AUC_IBRU / (AUC50,BP + AUC_IBRU) with AUC50,BP = 62.3 h*ng/mL (Ibrahim 2025 Table 2) -- about twice",
        "the IAUC50 of the companion efficacy model (28.4 h*ng/mL; modellib('Ibrahim_2025_ibrutinib_cll')), which is",
        "why the paper's de-escalation schedules shed hypertension risk faster than they shed efficacy."
      ),
      source_name        = "auc"
    ),
    AGE = list(
      description        = "Baseline age.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject (baseline value). Enters the SYSTOLIC mean transit time in power form",
        "mtt_sbp = exp(lmtt_sbp + etalmtt_sbp + e_age_mtt_sbp * log(AGE / 70)), equivalently",
        "MTT_sBP = 53.6 * (AGE / 70)^-3.87. Older patients reach the ibrutinib-elevated sBP steady state faster.",
        "There is no age effect on the diastolic mean transit time (Ibrahim 2025 Table 2).",
        "REFERENCE AGE. The authors' RxODE code (Data S1, PSP-2025-0220-s02.docx) writes the term UNCENTERED, as",
        "`effsbp_age_mtt*LNAGE` with LNAGE = log(AGE). Taken literally that would make the tabulated MTT_sBP of",
        "53.6 days the value at AGE = 1 year, which is not interpretable; the tabulated value must instead be the",
        "typical value at the population's reference age, so the term is centred here. 70 years is used because it is",
        "both the mean baseline age of the n = 246 analysis population (Table S1, 70 +/- 8.9 years) and the age at",
        "which the paper ran every dose-optimization simulation (Methods 2.4: 'in 1000 virtual patients ... at an age",
        "of 70 years'). The immediate predecessor paper did exactly this and stated it explicitly -- Ibrahim 2023",
        "Table 2 footnote c gives MTTsBP = e^(LN(79.9) - 5.04*LN(Age/63)) with 63 years the mean baseline age of that",
        "population. See the vignette Errata."
      ),
      source_name        = "LNAGE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 246L,
    n_studies      = 2L,
    age_range      = "mean 70 (SD 8.9) years",
    disease_state  = "Chronic lymphocytic leukemia (CLL); 151 (61%) treatment-naive, 95 (39%) relapsed/refractory",
    dose_range     = "ibrutinib 420 mg once daily (n = 94) or 840 mg once daily (n = 38) in PCYC-1102; 420 mg once daily in PCYC-1115",
    regions        = "United States and international (PCYC-1102 phase Ib/II; PCYC-1115 phase III)",
    notes          = paste(
      "Baseline demographics from Ibrahim 2025 Table S1 (Data S1, PSP-2025-0220-s01.docx), which reports only age and",
      "CLL group for the pooled n = 246 analysis population. Data were obtained through the Yale University Open Data",
      "Access (YODA) Project 2020-4386. Grade 2 hypertension is defined in the paper as sBP >= 140 mmHg or",
      "dBP >= 90 mmHg, and Grade 3 as sBP >= 160 mmHg or dBP >= 100 mmHg; these thresholds drive the",
      "toxicity-adjusted dosing simulations. Estimation used FOCEI in nlmixr 2.0.6; simulation used RxODE 1.1.1."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Values are the final estimates from Ibrahim 2025 Table 2. The
    # authors' Data S1 supplies the RxODE structure but NOT the fitted
    # THETA vector, so every value below is the back-transformed Table 2
    # point estimate. Each comment names the Table 2 row.
    #
    # SHARED PARAMETERS. Data S1 run6023_bp builds both endpoints' drug
    # effect from the SAME theta and the SAME random effect
    # (`dbpemax <- exp(tbpemax + eta.bpemax)` and
    #  `sbpemax <- exp(tbpemax + eta.bpemax)`; likewise `tbpec50` for both
    # EC50s). This is the coupling the paper highlights -- Results 3.2:
    # "enabling the estimation of a shared maximum ibrutinib stimulatory
    # effect (Emax = 0.0654) and the AUC0-24 at which 50% of the maximum
    # stimulation is achieved (EC50 = 62.3 h.ng.mL-1)". It is why sBP and
    # dBP belong in ONE model file rather than the two separate files the
    # 2023 predecessor used.
    # ------------------------------------------------------------------

    # --- baselines -----------------------------------------------------
    lrbase_sbp <- log(128)
    label("Log baseline systolic blood pressure (mmHg)")  # Table 2 sBPbaseline = 128 (95% CI 125-131)
    lrbase_dbp <- log(69.7)
    label("Log baseline diastolic blood pressure (mmHg)")  # Table 2 dBPbaseline = 69.7 (95% CI 67.8-71.6)

    # --- mean transit times --------------------------------------------
    lmtt_sbp <- log(53.6)
    label("Log mean transit time of the systolic transit/turnover chain (day)")  # Table 2 MTTsBP = 53.6 (95% CI 36.4-78.9)
    lmtt_dbp <- log(126)
    label("Log mean transit time of the diastolic transit/turnover chain (day)")  # Table 2 MTTdBP = 126 (95% CI 97.8-162)

    # --- shared ibrutinib effect ---------------------------------------
    lemax_bp  <- log(0.0654)
    label("Log maximum fractional ibrutinib stimulation of blood-pressure production, shared by sBP and dBP (unitless)")  # Table 2 Emax,BP = 0.0654 (95% CI 0.0476-0.0898)
    lauc50_bp <- log(62.3)
    label("Log daily AUC(0-24) giving 50% of the maximum blood-pressure stimulation, shared by sBP and dBP (h*ng/mL)")  # Table 2 AUC50,BP = 62.3 (95% CI 20.7-188)

    # --- covariate effect ----------------------------------------------
    e_age_mtt_sbp <- -3.87
    label("Power coefficient of baseline age (referenced to 70 years) on the systolic mean transit time")  # Table 2 'Coefficient of baseline age effect on MTTsBP' = -3.87 (95% CI -7.07 to -0.664); centred at 70 y, see covariateData$AGE

    ntransit_bp <- fixed(2)
    label("(n + 1) with n = 1 transit compartment; ktr = kout = ntransit_bp / MTT (unitless)")  # Data S1 run6023_bp `koutdbp = 2/dbpmtt`, `koutsbp = 2/sbpmtt`; Methods 2.3 'a turnover model combined with a single transit compartment'

    # --- IIV ------------------------------------------------------------
    # Table 2 reports IIV as CV%; ini() needs variances, so each CV% is
    # inverted with omega^2 = ln(1 + CV^2) and the arithmetic is shown.
    # Table 2 gives no CV% for MTTdBP or AUC50,BP, so neither carries a
    # random effect -- matching Data S1 run6023_bp, where `dbpmtt` and
    # `dbpec50` / `sbpec50` are built without an eta.
    etalrbase_sbp ~ 0.0077142  # Table 2 sBPbaseline CV% = 8.8 -> ln(1 + 0.088^2) = 0.0077142
    etalrbase_dbp ~ 0.0093649  # Table 2 dBPbaseline CV% = 9.7 -> ln(1 + 0.097^2) = 0.0093649
    etalmtt_sbp   ~ 0.455913   # Table 2 MTTsBP CV% = 76 -> ln(1 + 0.76^2) = 0.455913
    etalemax_bp   ~ 0.673139   # Table 2 Emax,BP CV% = 98 -> ln(1 + 0.98^2) = 0.673139; SHARED by sBP and dBP (Data S1 run6023_bp uses eta.bpemax in both)

    # --- residual error --------------------------------------------------
    # Table 2 footnote b: "An additive residual unexplained variability (RUV)
    # model was implemented on log-transformed data." The authors' code
    # writes `sbp_ipred = log(sbp + 1E-6) + prop.sd2`. Additive-on-log is
    # exactly nlmixr2's `lnorm()` residual, so both outputs stay in mmHg and
    # carry a log-normal residual instead of being emitted on the log scale.
    expSd_sbp <- 0.086
    label("Log-scale (log-normal) residual error SD for systolic blood pressure (unitless on the log scale)")  # Table 2 'RUV sBP (%) = 8.6'
    expSd_dbp <- 0.093
    label("Log-scale (log-normal) residual error SD for diastolic blood pressure (unitless on the log scale)")  # Table 2 'RUV dBP (%) = 9.3'
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual parameters. Baseline age enters the SYSTOLIC mean
    #    transit time in power form, centred at 70 years (see
    #    covariateData$AGE for why the printed uncentred form is centred).
    # ------------------------------------------------------------------
    sbpbase <- exp(lrbase_sbp + etalrbase_sbp)
    dbpbase <- exp(lrbase_dbp + etalrbase_dbp)

    mtt_sbp <- exp(lmtt_sbp + etalmtt_sbp + e_age_mtt_sbp * log(AGE / 70))
    mtt_dbp <- exp(lmtt_dbp)

    # Shared maximum effect and shared potency, with a single shared random
    # effect on Emax (Data S1 run6023_bp).
    emax_bp  <- exp(lemax_bp + etalemax_bp)
    auc50_bp <- exp(lauc50_bp)

    # ------------------------------------------------------------------
    # 2. Derived rate constants. The turnover rate constant was assumed
    #    equal to the transit rate constant, both (n + 1) / MTT with a
    #    single transit compartment.
    # ------------------------------------------------------------------
    ktr_sbp  <- ntransit_bp / mtt_sbp
    kout_sbp <- ktr_sbp
    kin_sbp  <- sbpbase * kout_sbp

    ktr_dbp  <- ntransit_bp / mtt_dbp
    kout_dbp <- ktr_dbp
    kin_dbp  <- dbpbase * kout_dbp

    # Emax stimulation of the zero-order production rates. Both endpoints
    # use the same emax_bp and auc50_bp.
    eff_sbp <- emax_bp * AUC_IBRU / (auc50_bp + AUC_IBRU)
    eff_dbp <- emax_bp * AUC_IBRU / (auc50_bp + AUC_IBRU)

    # ------------------------------------------------------------------
    # 3. ODE system (Data S1 run6023_bp).
    # ------------------------------------------------------------------
    d/dt(sbp_transit1) <- kin_sbp * (1 + eff_sbp) - ktr_sbp * sbp_transit1
    sbp_transit1(0)    <- sbpbase

    d/dt(sbp) <- ktr_sbp * sbp_transit1 - kout_sbp * sbp
    sbp(0)    <- sbpbase

    d/dt(dbp_transit1) <- kin_dbp * (1 + eff_dbp) - ktr_dbp * dbp_transit1
    dbp_transit1(0)    <- dbpbase

    d/dt(dbp) <- ktr_dbp * dbp_transit1 - kout_dbp * dbp
    dbp(0)    <- dbpbase

    # ------------------------------------------------------------------
    # 4. Observations. `sbp` and `dbp` are the registered canonical
    #    blood-pressure PD states and are themselves the observation
    #    variables.
    # ------------------------------------------------------------------
    sbp ~ lnorm(expSd_sbp)
    dbp ~ lnorm(expSd_dbp)
  })
}
