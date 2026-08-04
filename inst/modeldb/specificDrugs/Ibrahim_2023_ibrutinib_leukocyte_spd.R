Ibrahim_2023_ibrutinib_leukocyte_spd <- function() {
  description <- paste(
    "QSP / semi-mechanistic PK-PD. Integrated ibrutinib exposure-response model for leukocyte count and lymph-node",
    "burden (sum of the products of perpendicular diameters, SPD) in chronic lymphocytic leukemia (CLL), from the",
    "phase Ib/II PCYC-1102 study. Five ODEs: one relative-quantity turnover compartment for phosphorylated Bruton",
    "tyrosine kinase (pBtk) whose zero-order production is inhibited by an Imax function of the daily ibrutinib",
    "AUC(0-24), plus four CLL cell subpopulations that reproduce the Calissano CLL life cycle -- two proliferating",
    "stroma-attached colonies in lymphoid tissue with different detachment sensitivities (cll_subpop1, cll_subpop2),",
    "a released/activated tissue pool that exits to blood (cll_subpop3), and a resting peripheral-blood pool",
    "(cll_bld). Free (unphosphorylated) Btk drives four simultaneous drug effects: inhibition of CLL proliferation,",
    "ibrutinib-induced apoptosis of cll_subpop3, enhanced detachment of cll_subpop1/cll_subpop2 from the stroma, and",
    "blockade of homing from blood back to tissue -- the last of which produces the characteristic treatment-related",
    "lymphocytosis. Resistance is an exponential decay of the proliferation and apoptosis effects with time",
    "(t1/2 = 761 days). Two observed outputs: leukocyte count (10^9 cells/L) and total SPD (cm^2). Ibrutinib enters",
    "only through the time-varying covariate AUC_IBRU; there are no drug-dosing events in this model. Sister model",
    "files from the same paper: modellib('Ibrahim_2023_ibrutinib_sbp'), modellib('Ibrahim_2023_ibrutinib_dbp'),",
    "modellib('Ibrahim_2023_ibrutinib_competing_risk')."
  )
  reference <- paste(
    "Ibrahim EIK, Karlsson MO, Friberg LE.",
    "Assessment of ibrutinib scheduling on leukocyte, lymph node size and blood pressure dynamics in chronic",
    "lymphocytic leukemia through pharmacokinetic-pharmacodynamic models.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(9):1305-1318.",
    "doi:10.1002/psp4.13010.",
    "Open Access under CC BY-NC.",
    "Structural equations from Appendix S2 equations S4-S7 and main-text equation 1;",
    "parameter values from the authors' own nlmixr control stream (Supporting Information S4, run100)",
    "cross-checked against Table 1.",
    sep = " "
  )
  vignette <- "Ibrahim_2023_ibrutinib"
  paper_specific_compartments <- c(
    "pbtk", "cll_subpop1", "cll_subpop2", "cll_subpop3", "cll_bld"
  )
  units <- list(
    time          = "day",
    dosing        = "n/a (no drug-dosing events; ibrutinib exposure enters as the time-varying covariate AUC_IBRU)",
    concentration = "leukocyte count in 10^9 cells/L; SPD in cm^2 (neither output is a drug concentration)"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    pbtk        = list(analyte = "phosphorylated Bruton tyrosine kinase (pBtk)", units = NA_character_, specimen = "administration site", verified = FALSE),
    cll_subpop1 = list(analyte = "CLL cell subpopulation 1", units = NA_character_, specimen = "lymph", verified = FALSE),
    cll_subpop2 = list(analyte = "CLL cell subpopulation 2", units = NA_character_, specimen = "lymph", verified = FALSE),
    cll_subpop3 = list(analyte = "CLL cell subpopulation 3", units = NA_character_, specimen = "blood cell", verified = FALSE),
    cll_bld     = list(analyte = "CLL cells in blood", units = NA_character_, specimen = "blood cell", verified = FALSE)
  )

  covariateData <- list(
    AUC_IBRU = list(
      description        = "Daily 0-24 h area under the ibrutinib plasma concentration-time curve, AUC(0-24).",
      units              = "h*ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying: the value tracks the patient's current daily ibrutinib dose level (420, 280, or 140 mg/day in",
        "the schedules simulated by the paper; 420 and 840 mg/day in the PCYC-1102 study itself) and drops to 0",
        "during treatment interruptions. Ibrutinib enters this model ONLY through this column -- the model contains",
        "no drug compartment and no dosing events. Ibrahim 2023 derived per-subject AUC(0-24) by integrating the",
        "individual post hoc profiles of the two-compartment ibrutinib population PK model of Marostica et al.",
        "(Cancer Chemother Pharmacol. 2015;75(1):111-121), which is NOT part of nlmixr2lib; downstream users must",
        "supply AUC(0-24) from that model, from another on-disk ibrutinib popPK model, or from observed exposure.",
        "Enters the pBtk production-inhibition Imax function as AUC_IBRU / (IAUC50 + AUC_IBRU) with",
        "IAUC50 = 34.1 h*ng/mL (Ibrahim 2023 Table 1)."
      ),
      source_name        = "DAILYAUC"
    ),
    TUM_IGHV_MUT = list(
      description        = "Immunoglobulin heavy-chain variable region (IGHV) mutational status of the CLL clone: 1 = IGHV-mutated, 0 = IGHV-unmutated.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (IGHV-unmutated)",
      notes              = paste(
        "Time-fixed per subject. Selects the baseline SPD typical value:",
        "spdbase = spdbase_unm * (1 - TUM_IGHV_MUT) + spdbase_m * TUM_IGHV_MUT.",
        "IGHV-unmutated patients were estimated to have a 2.5-fold higher baseline SPD than IGHV-mutated patients",
        "(48.9 vs 19.5 cm^2; Ibrahim 2023 Results 'PK-SPD-leukocyte model' and Table 1). IGHV-mutated is the",
        "prognostically favourable group in CLL, so the indicator polarity matters: 1 = mutated = smaller baseline",
        "lymph-node burden. Missing categorical covariates (n <= 6) were imputed using the mode",
        "(Ibrahim 2023 Appendix S1 section 2), which is why the source column is named IMPIGVHMS ('imputed IGHV",
        "mutational status')."
      ),
      source_name        = "IMPIGVHMS"
    ),
    LINE_1L = list(
      description        = "First-line-therapy indicator: 1 = treatment-naive (TN) at baseline, 0 = relapsed/refractory (R/R).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (relapsed/refractory)",
      notes              = paste(
        "Time-fixed per subject. Selects the normal (non-leukaemic) peripheral-blood leukocyte number:",
        "nrmbld = nrmbld_tn * LINE_1L + nrmbld_rr * (1 - LINE_1L). Treatment-naive patients were estimated to have",
        "a 1.9-fold higher normal leukocyte number than relapsed/refractory patients (32.9 vs 17.3 x10^9 cells;",
        "Ibrahim 2023 Results 'PK-SPD-leukocyte model' and Table 1). POLARITY WARNING: the source column TRTARM is",
        "the INVERSE of this canonical -- the authors' nlmixr code (Supporting Information S4, run100) writes",
        "bldnrm = bldnrm_tn*(1-iarm) + bldnrm_rr*iarm with iarm = TRTARM, i.e. TRTARM = 0 is treatment-naive.",
        "Convert on ingestion with LINE_1L = 1 - TRTARM."
      ),
      source_name        = "TRTARM"
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
    dose_range     = "ibrutinib 420 mg once daily (n = 94) or 840 mg once daily (n = 38) in PCYC-1102; only the 120 patients with both leukocyte count and SPD measurements were analysed",
    regions        = "United States (PCYC-1102, phase Ib/II)",
    notes          = paste(
      "Baseline demographics from Ibrahim 2023 Supplementary Table S1. Follow-up to a maximum of 2.4 years",
      "(median 1.7 years), with periods of treatment interruption and dose reduction for adverse events.",
      "The analysis dataset contained 2374 ibrutinib plasma concentrations, 2434 leukocyte counts, 507 SPD",
      "measurements and 2413 paired sBP/dBP measurements; 11 patients died and 22 dropped out",
      "(Ibrahim 2023 Appendix S1 section 1). Estimation used FOCEI in nlmixr 2.0.6."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # All log-scale THETA values below are the final estimates taken
    # verbatim from the authors' own nlmixr control stream (Ibrahim 2023
    # Supporting Information S4, run100 `ini({...})`), which carries more
    # significant figures than the rounded back-transformed values printed
    # in Table 1. Each comment gives the control-stream symbol and the
    # matching Table 1 row so both sources can be audited.
    #
    # Note on the two Box-Cox-transformed baselines (CLLbld,0 and
    # SPDbaseline): exp(theta) is 64.35 and 48.72 whereas Table 1 reports
    # 64.1 and 48.9. The small difference is expected -- Table 1 reports the
    # typical value implied by the Box-Cox random-effect distribution while
    # the control stream reports the raw structural THETA. Every other
    # parameter back-transforms to its Table 1 value exactly.
    # ------------------------------------------------------------------

    # --- pBtk turnover -------------------------------------------------
    pbtk_base   <- fixed(1)
    label("Baseline phosphorylated-Btk level (relative quantity; 1 = 100% active)")  # Table 1 'pBTKbaseline (%) = 100 FIX'; S4 run100 base_btk = 1
    lkout_pbtk  <- 0.2986911
    label("Log turnover (elimination) rate constant of pBtk (1/day)")  # S4 run100 tkout = 0.2986911 -> 1.35; Table 1 kout,pBTK = 1.35 (95% CI 0.965-1.88)
    imax_pbtk   <- fixed(1)
    label("Maximum inhibition of pBtk production by ibrutinib (unitless)")  # Table 1 'IMAX = 1 FIX'
    liauc50_pbtk <- 3.5299932
    label("Log daily AUC(0-24) giving 50% of maximum pBtk-production inhibition (h*ng/mL)")  # S4 run100 tic50 = 3.5299932 -> 34.1; Table 1 IAUC50,pBTK = 34.1 (95% CI 16.5-70.7)

    # --- CLL cell baselines --------------------------------------------
    lrbase_cllbld  <- 4.1644135
    label("Log baseline CLL cell number in peripheral blood (10^9 cells)")  # S4 run100 tcbldbas = 4.1644135; Table 1 CLLbld,0 = 64.1 (95% CI 37-112.2)
    lrbase_clltiss <- 7.6449086
    label("Log baseline total proliferating CLL cell number in lymphoid tissues (10^9 cells)")  # S4 run100 tclnbas = 7.6449086 -> 2090; Table 1 CLLtiss,baseline = 2090 (95% CI 1320-3320)
    lrbase_spd_unm <- 3.8860955
    label("Log baseline SPD in IGHV-unmutated patients (cm^2)")  # S4 run100 tspdbas_unm = 3.8860955; Table 1 SPDunm,baseline = 48.9 (95% CI 31.5-75.1)
    lrbase_spd_m   <- 2.9725082
    label("Log baseline SPD in IGHV-mutated patients (cm^2)")  # S4 run100 tspdbas_m = 2.9725082; Table 1 SPDm,baseline = 19.5 (95% CI 6.68-57.4)
    lnrmbld_tn     <- 3.4941826
    label("Log normal (non-leukaemic) leukocyte number in peripheral blood, treatment-naive patients (10^9 cells)")  # S4 run100 tbldnrm_tn = 3.4941826 -> 32.9; Table 1 NRMbld,TN = 32.9 (95% CI 13.2-81.9)
    lnrmbld_rr     <- 2.8503700
    label("Log normal (non-leukaemic) leukocyte number in peripheral blood, relapsed/refractory patients (10^9 cells)")  # S4 run100 tbldnrm_rr = 2.8503700 -> 17.3; Table 1 NRMbld,RR = 17.3 (95% CI 12.1-24.8)
    lnrmln         <- 0.7843116
    label("Log normal (non-leukaemic) lymph-node size (cm^2)")  # S4 run100 tlnss = 0.7843116 -> 2.19; Table 1 NRMLN = 2.19 (95% CI 1.25-3.82)
    vbld           <- fixed(5)
    label("Blood volume used to convert peripheral-blood CLL cell number to leukocyte concentration (L)")  # Table 1 'VBld (L) = 5 FIX'; S4 run100 divides by 5

    # --- CLL cell kinetics ---------------------------------------------
    lkp      <- -5.4811627
    label("Log proliferation rate constant of CLL cells (1/day)")  # S4 run100 tkp1 = -5.4811627 -> 0.00416; Table 1 kp = 0.00416 (95% CI 0.00301-0.00577)
    lkh      <- -0.7571987
    label("Log homing rate constant of CLL cells from peripheral blood back to lymphoid tissue (1/day)")  # S4 run100 tkh = -0.7571987 -> 0.469; Table 1 kh = 0.469 (95% CI 0.273-0.805)
    lkd_bld  <- -4.8761701
    label("Log natural death rate constant of CLL cells in peripheral blood (1/day)")  # S4 run100 tkdb1 = -4.8761701 -> 0.00763 (t1/2 = 90.9 d); Table 1 kd,bld = 0.00763 (95% CI 0.00319-0.0182)
    lkd_tiss <- -1.7320378
    label("Log ibrutinib-induced death rate constant of CLL cells in lymphoid tissue (1/day)")  # S4 run100 tkdt1 = -1.7320378 -> 0.177 (t1/2 = 3.92 d); Table 1 kd,tiss = 0.177 (95% CI 0.0972-0.322)

    # --- ibrutinib effect slopes ---------------------------------------
    slp1  <- fixed(1)
    label("Slope of the ibrutinib inhibitory effect on proliferation (kp) and homing (kh) (unitless)")  # Table 1 'slp1 = 1 FIX'; S4 run100 slp1 = 1
    lslp2 <- 2.4260423
    label("Log slope of the ibrutinib stimulatory effect on detachment of cll_subpop1 from stroma (unitless)")  # S4 run100 tslp2 = 2.4260423 -> 11.3; Table 1 slp2 = 11.3 (95% CI 7.32-17.5); 1 + 11.3 = 12.3-fold enhancement
    lslp3 <- -2.3192759
    label("Log slope of the ibrutinib stimulatory effect on detachment of cll_subpop2 from stroma (unitless)")  # S4 run100 tslp3 = -2.3192759 -> 0.0983; Table 1 slp3 = 0.0983 (95% CI 0.0423-0.229); 1 + 0.0983 = 1.1-fold enhancement

    # --- subpopulation fractions (logit scale) --------------------------
    logit_frc1 <- -0.1152389
    label("Logit of frc1, the fraction of stroma-attached CLL cells that belong to the slow-detaching second subpopulation (unitless)")  # S4 run100 poplogit_fr1 = -0.1152389 -> frc1 = 0.471; Table 1 frc1 = 0.471 (95% CI 0.326-0.621)
    logit_frc2 <- -1.3357462
    label("Logit of frc2, the fraction of CLLtiss,baseline that belongs to the released third subpopulation (unitless)")  # S4 run100 poplogit_fr2 = -1.3357462 -> frc2 = 0.208; Table 1 frc2 = 0.208 (95% CI 0.123-0.331)

    # --- resistance -----------------------------------------------------
    llambda_dec <- -7.0012733
    label("Log exponential decay constant of the ibrutinib proliferation/apoptosis effect over time (1/day)")  # S4 run100 tlmbd = -7.0012733 -> 0.000911 (t1/2 = 761 d); Table 1 lambda_dec = 0.000911 (95% CI 0.000489-0.0017)

    # --- Box-Cox shape for the two baseline random effects --------------
    boxcox_base <- -0.1931613
    label("Box-Cox shape parameter shared by the CLLbld,0 and SPDbaseline random effects (unitless)")  # S4 run100 lambda1 = -0.1931613; Table 1 'Box-Cox for CLLbld,baseline and SPDbaseline random effects' = -0.193 (95% CI -0.398-0.0114)

    # --- IIV ($OMEGA variances, S4 run100) ------------------------------
    # Table 1 reports these as CV% via CV = sqrt(exp(omega^2) - 1); the
    # control stream reports the variances directly, which is what ini()
    # needs, so the variances are used here.
    etalrbase_cllbld  ~ 3.986072    # S4 run100 eta.cbldbas; Table 1 CLLbld,0 CV% = 727
    etalrbase_clltiss ~ 0.2703597   # S4 run100 eta.clnbas;  Table 1 CLLtiss,baseline CV% = 55.7
    etalnrmbld        ~ 0.3730409   # S4 run100 eta.bldnrm;  Table 1 NRMbld CV% = 67.2 (shared by the TN and R/R typical values)
    etalkd_bld        ~ 0.6146217   # S4 run100 eta.kdb1;    Table 1 kd,bld CV% = 92.1
    etalkh            ~ 1.67375     # S4 run100 eta.kh;      Table 1 kh CV% = 208
    etalkd_tiss       ~ 1.572316    # S4 run100 eta.kdt1;    Table 1 kd,tiss CV% = 195
    etalkp            ~ 0.5444774   # S4 run100 eta.kp1;     Table 1 kp CV% = 85.1
    etalogit_frc1     ~ 1.286532    # S4 run100 eta.frc1;    Table 1 frc1 CV% = 162
    etalogit_frc2     ~ 1.41655     # S4 run100 eta.frc2;    Table 1 frc2 CV% = 177
    etaliauc50_pbtk   ~ 2.281496    # S4 run100 eta.ic50;    Table 1 IAUC50,pBTK CV% = 297
    etalkout_pbtk     ~ 2.190465    # S4 run100 eta.kout;    Table 1 kout,pBTK CV% = 282
    etallambda_dec    ~ 0.8146673   # S4 run100 eta.lmbd;    Table 1 lambda_dec CV% = 112

    # Correlated block: baseline SPD and normal lymph-node size.
    # Table 1 footnote b: 'The 63% estimated correlation between omega^2 for
    # SPDbaseline and NRMLN'; check: 0.5053674 / sqrt(1.0127872 * 0.6420107) = 0.627.
    etalrbase_spd + etalnrmln ~ c(1.0127872,
                                  0.5053674, 0.6420107)  # S4 run100 eta.spdbas + eta.lnss; Table 1 SPDbaseline CV% = 132, NRMLN CV% = 94.9

    # --- residual error --------------------------------------------------
    # Table 1 footnote d: 'Additive residual unexplained variability (RUV)
    # model was implemented on log transformed data'. The authors' code
    # writes `ln_spd_ipred = log(tot_spd + 1E-6); ln_spd_ipred ~ add(add.sd1)`.
    # Additive-on-log is exactly nlmixr2's `lnorm()` residual, so the outputs
    # below stay in their natural units (cm^2 and 10^9 cells/L) and carry a
    # log-normal residual instead of being emitted on the log scale.
    expSd_tumorSpd  <- 0.2043994
    label("Log-scale (log-normal) residual error SD for total SPD (unitless on the log scale)")  # S4 run100 add.sd1 = 0.2043994; Table 1 'RUV SPD (%) = 20'
    expSd_leukocyte <- 0.2208791
    label("Log-scale (log-normal) residual error SD for leukocyte count (unitless on the log scale)")  # S4 run100 add.sd2 = 0.2208791; Table 1 'RUV leukocyte count (%) = 22'
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual parameters.
    #    The two baseline parameters carry Box-Cox-transformed random
    #    effects (Ibrahim 2023 S4 run100: etacbldbasbox / etaspdbasbox).
    # ------------------------------------------------------------------
    etabox_cllbld <- ((exp(etalrbase_cllbld)^boxcox_base) - 1) / boxcox_base
    etabox_spd    <- ((exp(etalrbase_spd)^boxcox_base) - 1) / boxcox_base

    cllbld0     <- exp(lrbase_cllbld + etabox_cllbld)
    clltiss0    <- exp(lrbase_clltiss + etalrbase_clltiss)
    nrmbld_tn   <- exp(lnrmbld_tn + etalnrmbld)
    nrmbld_rr   <- exp(lnrmbld_rr + etalnrmbld)
    spdbase_unm <- exp(lrbase_spd_unm + etabox_spd)
    spdbase_m   <- exp(lrbase_spd_m + etabox_spd)
    nrmln       <- exp(lnrmln + etalnrmln)

    kp       <- exp(lkp + etalkp)
    kh       <- exp(lkh + etalkh)
    kd_bld   <- exp(lkd_bld + etalkd_bld)
    kd_tiss  <- exp(lkd_tiss + etalkd_tiss)
    slp2     <- exp(lslp2)
    slp3     <- exp(lslp3)

    kout_pbtk   <- exp(lkout_pbtk + etalkout_pbtk)
    iauc50_pbtk <- exp(liauc50_pbtk + etaliauc50_pbtk)
    lambda_dec  <- exp(llambda_dec + etallambda_dec)

    odds_frc1 <- exp(logit_frc1 + etalogit_frc1)
    frc1      <- odds_frc1 / (1 + odds_frc1)
    odds_frc2 <- exp(logit_frc2 + etalogit_frc2)
    frc2      <- odds_frc2 / (1 + odds_frc2)

    # ------------------------------------------------------------------
    # 2. Covariate models (Ibrahim 2023 Results 'PK-SPD-leukocyte model').
    #    LINE_1L = 1 - TRTARM (see covariateData notes for the polarity).
    # ------------------------------------------------------------------
    nrmbld  <- nrmbld_tn * LINE_1L + nrmbld_rr * (1 - LINE_1L)
    spdbase <- spdbase_unm * (1 - TUM_IGHV_MUT) + spdbase_m * TUM_IGHV_MUT

    # ------------------------------------------------------------------
    # 3. Derived quantities.
    # ------------------------------------------------------------------
    # pBtk zero-order production and its Imax inhibition by ibrutinib
    # (Ibrahim 2023 equation 1).
    kin_pbtk <- kout_pbtk * pbtk_base
    effauc   <- imax_pbtk * AUC_IBRU / (iauc50_pbtk + AUC_IBRU)

    # Free (dephosphorylated) Btk is the pharmacological driver of every
    # downstream CLL effect (Appendix S2: pBTKeff = pBTKbaseline - pBTK).
    pbtk_eff <- pbtk_base - pbtk

    # Resistance: the proliferation and tissue-apoptosis effects decay
    # exponentially with time (Appendix S2: resist = exp(-lambda_dec * t)).
    resist <- exp(-lambda_dec * t)

    # Detachment rate constant was assumed equal to the proliferation rate
    # constant (Ibrahim 2023 Table 1 footnote c: 'kdtch was assumed to be
    # equal to kp'; S4 run100 krdt = kp1).
    kdtch <- kp

    eff_prol  <- slp1 * pbtk_eff * resist    # inhibition of CLL proliferation
    eff_death <- kd_tiss * pbtk_eff * resist # ibrutinib-induced death of cll_subpop3
    eff_home  <- slp1 * pbtk_eff             # blockade of homing blood -> tissue
    eff_dtch1 <- slp2 * pbtk_eff             # enhanced detachment of cll_subpop1
    eff_dtch2 <- slp3 * pbtk_eff             # enhanced detachment of cll_subpop2

    # Scaling factors between the SPD (cm^2) and cell-count (10^9 cells)
    # measurement scales (Ibrahim 2023 Results 'PK-SPD-leukocyte model';
    # S4 run100 sc1 / sc2).
    sc_cells_spd <- spdbase / clltiss0   # SC_cells-SPD: cell count -> SPD
    sc_spd_cells <- clltiss0 / spdbase   # SC_SPD-cells: SPD -> cell count

    # Pseudo-steady-state redistribution rate constant in the absence of
    # ibrutinib (Ibrahim 2023 main text, the kdist definition).
    cll_subpop3_base <- frc2 * clltiss0
    kdist            <- (kh + kd_bld) * (cllbld0 / cll_subpop3_base)

    # ------------------------------------------------------------------
    # 4. ODE system.
    #    pBtk: Ibrahim 2023 equation 1.
    #    CLL subpopulations: Appendix S2 equations S4-S7.
    # ------------------------------------------------------------------
    d/dt(pbtk) <- kin_pbtk * (1 - effauc) - kout_pbtk * pbtk
    pbtk(0)    <- pbtk_base

    # eq. S4 -- stroma-attached colony 1 (fast detachment)
    d/dt(cll_subpop1) <- kp * (1 - eff_prol) * cll_subpop1 -
      kdtch * (1 + eff_dtch1) * cll_subpop1
    cll_subpop1(0)    <- spdbase * (1 - frc2) * (1 - frc1)

    # eq. S5 -- stroma-attached colony 2 (slow detachment; slp2/slp3 = 115)
    d/dt(cll_subpop2) <- kp * (1 - eff_prol) * cll_subpop2 -
      kdtch * (1 + eff_dtch2) * cll_subpop2
    cll_subpop2(0)    <- spdbase * (1 - frc2) * frc1

    # eq. S6 -- released / activated tissue pool that exits to blood
    d/dt(cll_subpop3) <- kp * (1 - eff_prol) * cll_subpop3 +
      kdtch * (1 + eff_dtch1) * cll_subpop1 +
      kdtch * (1 + eff_dtch2) * cll_subpop2 +
      kh * (1 - eff_home) * cll_bld * sc_cells_spd -
      (kdist + eff_death) * cll_subpop3
    cll_subpop3(0)    <- spdbase * frc2

    # eq. S7 -- resting CLL pool in peripheral blood
    d/dt(cll_bld) <- kdist * cll_subpop3 * sc_spd_cells -
      kh * (1 - eff_home) * cll_bld -
      kd_bld * cll_bld
    cll_bld(0)    <- cllbld0

    # ------------------------------------------------------------------
    # 5. Observations (Ibrahim 2023 Figure 1 caption).
    #    Leukocyte count = (peripheral-blood CLL cells + normal leukocytes)
    #    divided by blood volume. SPD = the three tissue subpopulations
    #    plus the normal lymph-node size.
    # ------------------------------------------------------------------
    leukocyte <- (cll_bld + nrmbld) / vbld
    tumorSpd  <- cll_subpop1 + cll_subpop2 + cll_subpop3 + nrmln

    # Total CLL burden on the cell-count scale: the three lymphoid-tissue
    # subpopulations are carried on the SPD scale, so they are converted with
    # SC_SPD-cells before being added to the peripheral-blood pool. This is
    # the quantity summarised in Ibrahim 2023 Figure 4 ('relative change of
    # the individual predictions of total CLL count from their baselines').
    totalCll <- cll_bld +
      (cll_subpop1 + cll_subpop2 + cll_subpop3) * sc_spd_cells

    # Btk occupancy (%), the paper's reported target-engagement metric
    # (Ibrahim 2023 Results 'Evaluation of the de-escalation dosing
    # schedules': occupancy = (1 - pBtk) * 100).
    btkOccupancy <- (pbtk_base - pbtk) * 100

    leukocyte ~ lnorm(expSd_leukocyte)
    tumorSpd  ~ lnorm(expSd_tumorSpd)
  })
}
