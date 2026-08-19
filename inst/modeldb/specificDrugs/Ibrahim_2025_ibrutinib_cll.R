Ibrahim_2025_ibrutinib_cll <- function() {
  description <- paste(
    "QSP / semi-mechanistic PK-PD. Extended ibrutinib exposure-response model for chronic lymphocytic leukemia (CLL)",
    "cell dynamics, describing four biomarkers simultaneously: leukocyte count, lymphocyte count, lymph-node burden",
    "(sum of the products of perpendicular diameters, SPD) and spleen volume, pooled across the phase Ib/II PCYC-1102",
    "and phase III PCYC-1115 studies. Five ODEs: one relative-quantity turnover compartment for phosphorylated Bruton",
    "tyrosine kinase (pBtk) whose zero-order production is inhibited by an Imax function of the daily ibrutinib",
    "AUC(0-24), plus four CLL cell subpopulations -- two proliferating stroma-attached colonies in lymphoid tissue with",
    "different detachment sensitivities (cll_subpop1, cll_subpop2), a released tissue pool that exits to blood",
    "(cll_subpop3), and a resting peripheral-blood pool (cll_bld). Free (unphosphorylated) Btk drives four",
    "simultaneous drug effects: inhibition of CLL proliferation, ibrutinib-induced apoptosis of cll_subpop3, enhanced",
    "detachment of cll_subpop1/cll_subpop2 from the stroma, and blockade of homing from blood back to tissue -- the",
    "last of which produces the characteristic treatment-related lymphocytosis. This 2025 re-estimation carries all",
    "four CLL subpopulations on the cell-count scale and converts to SPD and spleen volume only at observation, and",
    "adds treatment-naive (TN) versus relapsed/refractory (R/R) status as a binary covariate on the pBtk turnover rate,",
    "baseline SPD, peripheral-blood CLL death rate, baseline blood CLL count, normal leukocyte count, and the presence",
    "of acquired resistance. Ibrutinib enters only through the time-varying covariate AUC_IBRU; there are no drug-dosing",
    "events in this model. Sister model file from the same paper:",
    "modellib('Ibrahim_2025_ibrutinib_bp'). Venetoclax combination extension:",
    "modellib('Ibrahim_2025_ibrutinib_venetoclax'). Predecessor framework:",
    "modellib('Ibrahim_2023_ibrutinib_leukocyte_spd')."
  )
  reference <- paste(
    "Ibrahim EIK, Friberg LE.",
    "Optimizing ibrutinib posology in chronic lymphocytic leukemia using a semi-mechanistic pharmacometric framework.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(12):2186-2198.",
    "doi:10.1002/psp4.70124.",
    "Open Access under CC BY-NC.",
    "Structural equations transcribed from the authors' own RxODE control stream",
    "(Data S1, PSP-2025-0220-s02.docx, `run8634_eff`); observation equations from main-text equations 2-5;",
    "parameter values from Table 1.",
    "Model structure inherited from Ibrahim EIK, Karlsson MO, Friberg LE.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(9):1305-1318. doi:10.1002/psp4.13010;",
    "see modellib('Ibrahim_2023_ibrutinib_leukocyte_spd').",
    sep = " "
  )
  vignette <- "Ibrahim_2025_ibrutinib"
  paper_specific_compartments <- c(
    "pbtk", "cll_subpop1", "cll_subpop2", "cll_subpop3", "cll_bld"
  )
  units <- list(
    time          = "day",
    dosing        = "n/a (no drug-dosing events; ibrutinib exposure enters as the time-varying covariate AUC_IBRU)",
    concentration = "leukocyte and lymphocyte counts in 10^9 cells/L; SPD in cm^2; spleen volume in cc (no output is a drug concentration)"
  )

  compartmentData <- list(
    pbtk        = list(analyte = "phosphorylated Bruton tyrosine kinase (pBtk)", units = "relative quantity (1 = 100% of baseline)", specimen = "not applicable", verified = TRUE),
    cll_subpop1 = list(analyte = "CLL cells, stroma-attached proliferating colony 1 (fast detachment)", units = "10^9 cells", specimen = "lymph", verified = TRUE),
    cll_subpop2 = list(analyte = "CLL cells, stroma-attached proliferating colony 2 (slow detachment)", units = "10^9 cells", specimen = "lymph", verified = TRUE),
    cll_subpop3 = list(analyte = "CLL cells, released lymphoid-tissue pool exiting to blood", units = "10^9 cells", specimen = "lymph", verified = TRUE),
    cll_bld     = list(analyte = "CLL cells, resting peripheral-blood pool", units = "10^9 cells", specimen = "blood cell", verified = TRUE)
  )

  covariateData <- list(
    AUC_IBRU = list(
      description        = "Daily 0-24 h area under the ibrutinib plasma concentration-time curve, AUC(0-24).",
      units              = "h*ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying: the value tracks the patient's current daily ibrutinib dose level (420, 280 or 140 mg/day in",
        "the schedules simulated by the paper; 420 and 840 mg/day in PCYC-1102 and 420 mg/day in PCYC-1115) and drops",
        "to 0 during treatment interruptions. Ibrutinib enters this model ONLY through this column -- the model",
        "contains no drug compartment and no dosing events. Ibrahim 2025 computed daily AUC(0-24) from the individual",
        "ibrutinib plasma concentrations using the two-compartment population PK model of Marostica et al.",
        "(Cancer Chemother Pharmacol. 2015;75(1):111-121), which has linear elimination and a sequential",
        "zero-/first-order absorption and is NOT part of nlmixr2lib; downstream users must supply AUC(0-24) from that",
        "model, from another ibrutinib popPK model, or from observed exposure. Enters the pBtk production-inhibition",
        "Imax function as AUC_IBRU / (IAUC50 + AUC_IBRU) with IAUC50 = 28.4 h*ng/mL (Ibrahim 2025 Table 1)."
      ),
      source_name        = "auc"
    ),
    LINE_1L = list(
      description        = "First-line-therapy indicator: 1 = treatment-naive (TN) at baseline, 0 = relapsed/refractory (R/R).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (relapsed/refractory)",
      notes              = paste(
        "Time-fixed per subject. This is the single covariate of the 2025 re-estimation and it acts on six places in",
        "the model (Ibrahim 2025 Results section 3.1 and Table 1):",
        "(1) pBtk turnover rate -- R/R patients have a 1.76-fold higher kout,pBtk, equivalently a 76% longer pBtk",
        "half-life in TN (check: ln2/0.524 = 1.32 d in TN vs ln2/(0.524*1.76) = 0.752 d in R/R, ratio 1.76);",
        "(2) baseline SPD -- 1.58-fold higher in R/R;",
        "(3) peripheral-blood CLL death rate -- R/R multiplier 0.571, i.e. 43% lower in R/R, equivalently a 43%",
        "shorter peripheral CLL cell half-life in TN (check: ln2/0.0153 = 45.3 d in TN vs ln2/(0.0153*0.571) = 79.3 d",
        "in R/R, so TN is 42.9% shorter);",
        "(4) baseline peripheral-blood CLL count -- separate typical values, 208 (TN) vs 53.2 (R/R) x10^9 cells,",
        "a 3.91-fold ratio;",
        "(5) normal leukocyte count -- separate typical values, 37.3 (TN) vs 19.9 (R/R) x10^9 cells, a 1.87-fold ratio;",
        "(6) acquired ibrutinib resistance -- the exponential decay of the proliferation and tissue-apoptosis effects",
        "is active in R/R patients ONLY and is switched off in TN patients (see the lambda_dec ini() comment).",
        "POLARITY WARNING: the source column `arm` is the INVERSE of this canonical -- the authors' RxODE code",
        "(Data S1, PSP-2025-0220-s02.docx) writes `iarm <- arm` and then",
        "`cbldbas <- cbldbas_tn*(1-iarm) + cbldbas_rr*(iarm)`, i.e. arm = 0 is treatment-naive.",
        "Convert on ingestion with LINE_1L = 1 - arm."
      ),
      source_name        = "arm"
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
      "CLL group for the pooled n = 246 analysis population; no weight, sex or race breakdown is given in the 2025",
      "paper, so those fields are omitted rather than carried over from the 2023 predecessor (whose population is a",
      "subset). The analysis included the 120 PCYC-1102 patients with ibrutinib concentrations plus leukocyte count",
      "and SPD, and the 126 PCYC-1115 patients with ibrutinib concentrations plus leukocyte or lymphocyte count and",
      "SPD or spleen volume. Data were obtained through the Yale University Open Data Access (YODA) Project",
      "2020-4386. Estimation used FOCEI in nlmixr 2.0.6; simulation used RxODE 1.1.1."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Values are the final estimates from Ibrahim 2025 Table 1. Unlike the
    # 2023 predecessor, the authors' Data S1 supplies the RxODE structure
    # but NOT the fitted THETA vector, so every value below is the
    # back-transformed Table 1 point estimate rather than a control-stream
    # log-scale THETA. Each comment names the Table 1 row.
    #
    # TIME UNIT. The model integrates in DAYS (kout,pBtk, kh, kd,bld and
    # kd,tiss are all reported in day^-1 and are used directly against `t`
    # in the authors' code). Table 1 reports kp and lambda_dec in month^-1
    # for readability; both are converted here with 30 days per month. That
    # conversion is not stated in the paper and is an assumption -- see the
    # vignette Errata. It is supported by the 2023 predecessor, which
    # reported the same parameters in day^-1: kp = 0.00416/day there vs
    # 0.124/30 = 0.004133/day here (ratio 0.99), comfortably inside the
    # "0.7- to 1.6-fold" difference from the previous publication that the
    # 2025 Discussion reports for all re-estimated PD parameters.
    # ------------------------------------------------------------------

    # --- pBtk turnover -------------------------------------------------
    pbtk_base    <- fixed(1)
    label("Baseline phosphorylated-Btk level (relative quantity; 1 = 100% active)")  # Table 1 'pBtkbaseline (%) = 100 FIX'; Data S1 run8634_eff base_btk = 1
    lkout_pbtk   <- log(0.524)
    label("Log turnover (elimination) rate constant of pBtk in treatment-naive patients (1/day)")  # Table 1 kout,pBtk,TN = 0.524 (95% CI 0.496-0.554)
    e_rr_kout_pbtk <- log(1.76)
    label("Log fold-change in the pBtk turnover rate constant in relapsed/refractory patients (unitless)")  # Table 1 'R/R effect on kout,pBtk = 1.76' (95% CI 0.631-4.90); enters as exp(lkout_pbtk + e_rr_kout_pbtk) in R/R
    imax_pbtk    <- fixed(1)
    label("Maximum inhibition of pBtk production by ibrutinib (unitless)")  # Table 1 'IMAX = 1 FIX'
    liauc50_pbtk <- log(28.4)
    label("Log daily AUC(0-24) giving 50% of maximum pBtk-production inhibition (h*ng/mL)")  # Table 1 IAUC50,pBtk = 28.4 (95% CI 20.4-39.5)

    # --- CLL cell and biomarker baselines -------------------------------
    lrbase_cllbld_tn  <- log(208)
    label("Log baseline CLL cell number in peripheral blood, treatment-naive patients (10^9 cells)")  # Table 1 CLLbld,baseline,TN = 208 (95% CI 145-299)
    lrbase_cllbld_rr  <- log(53.2)
    label("Log baseline CLL cell number in peripheral blood, relapsed/refractory patients (10^9 cells)")  # Table 1 CLLbld,baseline,R/R = 53.2 (95% CI 32.5-87.3)
    lrbase_clltiss    <- log(3279)
    label("Log baseline total proliferating CLL cell number in lymphoid tissues (10^9 cells)")  # Table 1 CLLtiss,baseline = 3279 (95% CI 2420-4430)
    lrbase_spd        <- log(24.3)
    label("Log baseline SPD in treatment-naive patients (cm^2)")  # Table 1 SPDbaseline,TN = 24.3 (95% CI 17.6-33.6)
    e_rr_rbase_spd    <- log(1.58)
    label("Log fold-change in baseline SPD in relapsed/refractory patients (unitless)")  # Table 1 'R/R effect on SPDbaseline = 1.58' (95% CI 1.00-2.49)
    lrbase_spleen     <- log(314)
    label("Log baseline spleen volume (cc)")  # Table 1 SPLEENbaseline = 314 (95% CI 214-463)
    lleuknrm_tn       <- log(37.3)
    label("Log normal (non-leukaemic) leukocyte number in peripheral blood, treatment-naive patients (10^9 cells)")  # Table 1 LEUKnrm,TN = 37.3 (95% CI 32.6-42.6)
    lleuknrm_rr       <- log(19.9)
    label("Log normal (non-leukaemic) leukocyte number in peripheral blood, relapsed/refractory patients (10^9 cells)")  # Table 1 LEUKnrm,R/R = 19.9 (95% CI 14-28.3)
    lspdnrm           <- log(2.70)
    label("Log normal (non-leukaemic) lymph-node size (cm^2)")  # Table 1 SPDnrm = 2.70 (95% CI 1.98-3.67)
    lspleennrm        <- log(252)
    label("Log normal (non-leukaemic) spleen volume (cc)")  # Table 1 SPLEENnrm = 252 (95% CI 200-317)
    vbld              <- fixed(5)
    label("Blood volume used to convert peripheral-blood cell numbers to counts per litre (L)")  # Table 1 'Vbld (L) = 5 FIX'; Data S1 run8634_eff divides by 5
    logit_flymphocyte <- log(0.310 / (1 - 0.310))
    label("Logit of the fraction of the normal lymphocyte count relative to the normal leukocyte count (unitless)")  # Table 1 f_lymphocyte = 0.310 (95% CI 0.253-0.374); Data S1 codes it on the logit scale

    # --- CLL cell kinetics ---------------------------------------------
    lkp        <- log(0.124 / 30)
    label("Log proliferation rate constant of CLL cells (1/day)")  # Table 1 kp = 0.124 month^-1 (95% CI 0.102-0.153); /30 days = 0.004133/day
    lkh        <- log(0.573)
    label("Log homing rate constant of CLL cells from peripheral blood back to lymphoid tissue (1/day)")  # Table 1 kh = 0.573 (95% CI 0.485-0.678)
    lkd_bld    <- log(0.0153)
    label("Log natural death rate constant of CLL cells in peripheral blood, treatment-naive patients (1/day)")  # Table 1 kd,bld,TN = 0.0153 (95% CI 0.0119-0.0196)
    e_rr_kd_bld <- log(0.571)
    label("Log fold-change in the peripheral-blood CLL death rate constant in relapsed/refractory patients (unitless)")  # Table 1 'R/R effect on kd,bld = 0.571' (95% CI 0.440-0.742); 43% lower in R/R
    lkd_tiss   <- log(0.162)
    label("Log ibrutinib-induced death rate constant of CLL cells in lymphoid tissue (1/day)")  # Table 1 kd,tiss = 0.162 (95% CI 0.119-0.22)

    # --- ibrutinib effect slopes ---------------------------------------
    slp1  <- fixed(1)
    label("Slope of the ibrutinib inhibitory effect on proliferation (kp) and homing (kh) (unitless)")  # Table 1 'slp1 = 1 FIX'; Data S1 run8634_eff slp1 = 1
    lslp2 <- log(12.0)
    label("Log slope of the ibrutinib stimulatory effect on detachment of cll_subpop1 from stroma (unitless)")  # Table 1 slp2 = 12.0 (95% CI 9.79-14.7); 1 + 12.0 = 13-fold enhancement at full Btk inhibition
    lslp3 <- log(0.140)
    label("Log slope of the ibrutinib stimulatory effect on detachment of cll_subpop2 from stroma (unitless)")  # Table 1 slp3 = 0.140 (95% CI 0.106-0.183); 1 + 0.140 = 1.14-fold enhancement at full Btk inhibition

    # --- subpopulation fractions (logit scale) --------------------------
    logit_f1 <- log(0.483 / (1 - 0.483))
    label("Logit of f1, the fraction of stroma-attached CLL cells that belong to the slow-detaching second subpopulation (unitless)")  # Table 1 f1 = 0.483 (95% CI 0.444-0.522); Data S1 codes it on the logit scale
    logit_f2 <- log(0.246 / (1 - 0.246))
    label("Logit of f2, the fraction of CLLtiss,baseline that belongs to the released third subpopulation (unitless)")  # Table 1 f2 = 0.246 (95% CI 0.202-0.297); Data S1 codes it on the logit scale

    # --- resistance -----------------------------------------------------
    # TABLE-LABEL TRANSPOSITION. Table 1 prints 'lambda_dec,TN = 0.0230
    # (month^-1)' and 'lambda_dec,R/R = 0 FIX'. Four independent statements
    # in the same paper say the opposite, i.e. that resistance is absent in
    # TN and present in R/R:
    #   (a) Data S1 run8634_eff: `resist = exp(-lmbd*iarm*t)`, and the same
    #       file fixes iarm = 1 for R/R (`cbldbas <- cbldbas_tn*(1-iarm) +
    #       cbldbas_rr*(iarm)`, `kout <- exp(tkout_tn + tkouteff_rr*(iarm))`).
    #       With iarm = 0 in TN the decay term collapses to exp(0) = 1, so
    #       TN patients carry NO resistance and the single estimated lmbd
    #       belongs to R/R.
    #   (b) Abstract: "with no evidence of ibrutinib resistance in TN patients".
    #   (c) Results 3.1: "Resistance to ibrutinib was not apparent in the TN
    #       patients."
    #   (d) Discussion: "the absence of resistance development to ibrutinib
    #       within the analyzed timeframe in TN patients, compared to R/R
    #       patients".
    # The executable code and three prose statements are taken over the two
    # table row labels, so the estimate is assigned to R/R below and TN is
    # fixed to no resistance. See the vignette Errata.
    llambda_dec <- log(0.0230 / 30)
    label("Log exponential decay constant of the ibrutinib proliferation/apoptosis effect in relapsed/refractory patients (1/day)")  # Table 1 lambda_dec = 0.0230 month^-1 (95% CI 0.0109-0.0484), labelled 'TN' but belonging to R/R (see comment above); /30 days = 0.000767/day

    # --- Box-Cox shape for the baseline blood-CLL random effects --------
    boxcox_cllbld <- -0.199
    label("Box-Cox shape parameter for the baseline peripheral-blood CLL count random effects (unitless)")  # Table 1 'Box-Cox for CLLbld,baseline random effects' = -0.199 (95% CI -0.329 to -0.0684)

    # --- shared-eta scaling for kd,bld ----------------------------------
    # Table 1 footnote b: "100% correlation between omega^2 for
    # CLLtiss,baseline and kd,bld,TN". Data S1 run8634_eff implements that
    # perfect correlation by REUSING the CLLtiss,baseline random effect,
    # scaled: `kdb1 <- exp(tkdb1_tn + tkdb1eff_rr*(iarm) + eta.clnbas*sc_kdb)`.
    # sc_kdb is not tabulated, but it is pinned by the two reported CV%:
    #   omega(CLLtiss,baseline) = sqrt(ln(1 + 1.52^2)) = 1.0942  (CV 152%)
    #   omega(kd,bld,TN)        = sqrt(ln(1 + 1.12^2)) = 0.9016  (CV 112%)
    #   sc_kdb = 0.9016 / 1.0942 = 0.824
    sc_kdb <- fixed(0.824)
    label("Scaling of the shared CLLtiss,baseline random effect onto kd,bld, reproducing their 100% correlation (unitless)")  # derived from Table 1 CV% for CLLtiss,baseline (152%) and kd,bld,TN (112%); Table 1 footnote b

    # --- IIV ------------------------------------------------------------
    # Table 1 reports IIV as CV%. Footnote a: CV = sqrt(exp(omega^2) - 1)
    # for every parameter EXCEPT CLLbld,baseline, where CV = sqrt(omega^2)
    # (the Box-Cox-transformed baselines). ini() needs variances, so each
    # CV% is inverted below and the arithmetic is shown in the comment.
    etalkout_pbtk       ~ 3.6997     # Table 1 kout,pBtk CV% = 628 -> ln(1 + 6.28^2) = 3.6997
    etalrbase_cllbld_tn ~ 1.7689     # Table 1 CLLbld,baseline,TN CV% = 133 -> footnote a: omega^2 = 1.33^2 = 1.7689
    etalrbase_cllbld_rr ~ 4.1616     # Table 1 CLLbld,baseline,R/R CV% = 204 -> footnote a: omega^2 = 2.04^2 = 4.1616
    etalogit_flymphocyte ~ 0.37085   # Table 1 f_lymphocyte CV% = 67 -> ln(1 + 0.67^2) = 0.37085
    etalrbase_clltiss   ~ 1.19726    # Table 1 CLLtiss,baseline CV% = 152 -> ln(1 + 1.52^2) = 1.19726; also drives kd,bld via sc_kdb
    etalrbase_spleen    ~ 0.78300    # Table 1 SPLEENbaseline CV% = 109 -> ln(1 + 1.09^2) = 0.78300
    etalleuknrm_tn      ~ 0.115556   # Table 1 LEUKnrm,TN CV% = 35 -> ln(1 + 0.35^2) = 0.115556
    etalleuknrm_rr      ~ 0.417684   # Table 1 LEUKnrm,R/R CV% = 72 -> ln(1 + 0.72^2) = 0.417684
    etalspleennrm       ~ 0.316354   # Table 1 SPLEENnrm CV% = 61 -> ln(1 + 0.61^2) = 0.316354
    etalogit_f1         ~ 2.74188    # Table 1 f1 CV% = 381 -> ln(1 + 3.81^2) = 2.74188
    etalogit_f2         ~ 1.09463    # Table 1 f2 CV% = 141 -> ln(1 + 1.41^2) = 1.09463
    etalkp              ~ 0.75308    # Table 1 kp CV% = 106 -> ln(1 + 1.06^2) = 0.75308
    etalkh              ~ 3.72440    # Table 1 kh CV% = 636 -> ln(1 + 6.36^2) = 3.72440
    etaliauc50_pbtk     ~ 2.99621    # Table 1 IAUC50,pBtk CV% = 436 -> ln(1 + 4.36^2) = 2.99621
    etalkd_tiss         ~ 1.64131    # Table 1 kd,tiss CV% = 204 -> ln(1 + 2.04^2) = 1.64131
    etallambda_dec      ~ 1.20639    # Table 1 lambda_dec CV% = 153 -> ln(1 + 1.53^2) = 1.20639

    # Inter-radiologist variability on the normal lymph-node size. Table 1
    # footnote d: "94% estimated inter-radiologist variability in SPDnrm".
    # Data S1 run8634_eff carries this as three reader-specific random
    # effects selected by a RAD reader-ID column
    # (`irv.lnss <- irv_lnss*(irv_rad1*eta.lnssrad1 + irv_rad2*eta.lnssrad2 +
    # irv_rad3*eta.lnssrad3)`). The three readers draw from one common
    # distribution and exactly one applies to any given reading, so the
    # marginal distribution is reproduced by the single random effect below
    # and no reader-ID covariate column is required. See the vignette Errata.
    etalspdnrm_rad      ~ 0.63313    # Table 1 footnote d inter-radiologist variability in SPDnrm = 94% -> ln(1 + 0.94^2) = 0.63313

    # Correlated block: baseline SPD and normal lymph-node size.
    # Table 1 footnote b: "65% estimated correlation between omega^2 for
    # SPDbaseline and SPDnrm".
    #   omega^2(SPDbaseline) = ln(1 + 1.30^2) = 0.98954   (CV 130%)
    #   omega^2(SPDnrm)      = ln(1 + 0.68^2) = 0.380093  (CV 68%)
    #   covariance = 0.65 * sqrt(0.98954 * 0.380093) = 0.398635
    etalrbase_spd + etalspdnrm ~ c(0.98954,
                                   0.398635, 0.380093)  # Table 1 SPDbaseline,TN CV% = 130, SPDnrm CV% = 68, footnote b correlation 65%

    # --- residual error --------------------------------------------------
    # Table 1 footnote e: "An additive residual unexplained variability (RUV)
    # model was implemented on log-transformed data." The authors' code
    # writes e.g. `wbc_bld_ipred = log(wbc_bld + 1E-6) + prop.sd1`.
    # Additive-on-log is exactly nlmixr2's `lnorm()` residual, so the outputs
    # below stay in their natural units and carry a log-normal residual
    # instead of being emitted on the log scale.
    expSd_leukocyte     <- 0.20
    label("Log-scale (log-normal) residual error SD for leukocyte count (unitless on the log scale)")  # Table 1 'RUV leukocyte count (%) = 20'
    expSd_lymphocyte    <- 0.21
    label("Log-scale (log-normal) residual error SD for lymphocyte count (unitless on the log scale)")  # Table 1 'RUV lymphocyte count (%) = 21'
    expSd_tumorSpd      <- 0.27
    label("Log-scale (log-normal) residual error SD for total SPD (unitless on the log scale)")  # Table 1 'RUV SPD (%) = 27'
    expSd_spleenVolume  <- 0.13
    label("Log-scale (log-normal) residual error SD for spleen volume (unitless on the log scale)")  # Table 1 'RUV spleen volume (%) = 13'
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual parameters.
    #    The two baseline peripheral-blood CLL counts carry Box-Cox
    #    transformed random effects (Data S1 run8634_eff:
    #    etacbldbasbox_tn / etacbldbasbox_rr).
    # ------------------------------------------------------------------
    etabox_cllbld_tn <- ((exp(etalrbase_cllbld_tn)^boxcox_cllbld) - 1) / boxcox_cllbld
    etabox_cllbld_rr <- ((exp(etalrbase_cllbld_rr)^boxcox_cllbld) - 1) / boxcox_cllbld

    cllbld0_tn <- exp(lrbase_cllbld_tn + etabox_cllbld_tn)
    cllbld0_rr <- exp(lrbase_cllbld_rr + etabox_cllbld_rr)

    clltiss0   <- exp(lrbase_clltiss + etalrbase_clltiss)
    spleenbase <- exp(lrbase_spleen + etalrbase_spleen)
    spleennrm  <- exp(lspleennrm + etalspleennrm)
    leuknrm_tn <- exp(lleuknrm_tn + etalleuknrm_tn)
    leuknrm_rr <- exp(lleuknrm_rr + etalleuknrm_rr)

    # Normal lymph-node size carries both between-subject and
    # inter-radiologist variability (Table 1 footnote d). The two random
    # effects are applied as a product of exponentials rather than as
    # `exp(theta + eta1 + eta2)`, which rxode2 cannot mu-reference
    # ("currently do not theta + eta1 + eta2"); the two forms are
    # algebraically identical.
    spdnrm <- exp(lspdnrm + etalspdnrm) * exp(etalspdnrm_rad)

    kp      <- exp(lkp + etalkp)
    kh      <- exp(lkh + etalkh)
    kd_tiss <- exp(lkd_tiss + etalkd_tiss)
    slp2    <- exp(lslp2)
    slp3    <- exp(lslp3)

    # The two fixed effects are parenthesised together in every
    # covariate-adjusted parameter below: rxode2's mu-referencing rejects
    # two bare population parameters in one expression ("2+ single
    # population parameters in a single mu-referenced expression"), and the
    # parentheses resolve it without changing the arithmetic.
    kout_pbtk_tn   <- exp(lkout_pbtk + etalkout_pbtk)
    kout_pbtk_rr   <- exp((lkout_pbtk + e_rr_kout_pbtk) + etalkout_pbtk)
    iauc50_pbtk    <- exp(liauc50_pbtk + etaliauc50_pbtk)
    lambda_dec     <- exp(llambda_dec + etallambda_dec)

    # kd,bld reuses the CLLtiss,baseline random effect, scaled, which is how
    # the authors encoded the 100% correlation of Table 1 footnote b.
    kd_bld_tn <- exp(lkd_bld + etalrbase_clltiss * sc_kdb)
    kd_bld_rr <- exp((lkd_bld + e_rr_kd_bld) + etalrbase_clltiss * sc_kdb)

    odds_f1 <- exp(logit_f1 + etalogit_f1)
    f1      <- odds_f1 / (1 + odds_f1)
    odds_f2 <- exp(logit_f2 + etalogit_f2)
    f2      <- odds_f2 / (1 + odds_f2)

    odds_flymphocyte <- exp(logit_flymphocyte + etalogit_flymphocyte)
    f_lymphocyte     <- odds_flymphocyte / (1 + odds_flymphocyte)

    # ------------------------------------------------------------------
    # 2. Covariate model (Ibrahim 2025 Table 1; Data S1 run8634_eff).
    #    LINE_1L = 1 - arm, so LINE_1L = 1 is treatment-naive.
    # ------------------------------------------------------------------
    cllbld0  <- cllbld0_tn * LINE_1L + cllbld0_rr * (1 - LINE_1L)
    leuknrm  <- leuknrm_tn * LINE_1L + leuknrm_rr * (1 - LINE_1L)
    kout_pbtk <- kout_pbtk_tn * LINE_1L + kout_pbtk_rr * (1 - LINE_1L)
    kd_bld    <- kd_bld_tn * LINE_1L + kd_bld_rr * (1 - LINE_1L)
    spdbase   <- exp((lrbase_spd + e_rr_rbase_spd * (1 - LINE_1L)) + etalrbase_spd)

    # ------------------------------------------------------------------
    # 3. Derived quantities.
    # ------------------------------------------------------------------
    # pBtk zero-order production and its Imax inhibition by ibrutinib.
    kin_pbtk <- kout_pbtk * pbtk_base
    effauc   <- imax_pbtk * AUC_IBRU / (iauc50_pbtk + AUC_IBRU)

    # Free (dephosphorylated) Btk is the pharmacological driver of every
    # downstream CLL effect (Data S1 run8634_eff writes this as `1 - btk`).
    pbtk_eff <- pbtk_base - pbtk

    # Resistance: the proliferation and tissue-apoptosis effects decay
    # exponentially with time in relapsed/refractory patients only
    # (Data S1 run8634_eff `resist = exp(-lmbd*iarm*t)` with iarm = 1 in
    # R/R; see the llambda_dec ini() comment for the Table 1 label
    # transposition).
    resist <- exp(-lambda_dec * (1 - LINE_1L) * t)

    # Detachment rate constant was assumed equal to the proliferation rate
    # constant (Table 1 footnote c: "kdtch was assumed to be equal to kp";
    # Data S1 run8634_eff `krdt <- kp1`).
    kdtch <- kp

    eff_prol  <- slp1 * resist * pbtk_eff     # inhibition of CLL proliferation
    eff_death <- kd_tiss * resist * pbtk_eff  # ibrutinib-induced death of cll_subpop3
    eff_home  <- slp1 * pbtk_eff              # blockade of homing blood -> tissue
    eff_dtch1 <- slp2 * pbtk_eff              # enhanced detachment of cll_subpop1
    eff_dtch2 <- slp3 * pbtk_eff              # enhanced detachment of cll_subpop2

    # Pseudo-steady-state redistribution rate constant from the released
    # tissue pool into blood, set so the system is at steady state at
    # baseline in the absence of ibrutinib (Data S1 run8634_eff:
    # `crdlnbas = frc2*clnbas; krd = (kh + kdb1)*(cbldbas/crdlnbas)`).
    cll_subpop3_base <- f2 * clltiss0
    kdist            <- (kh + kd_bld) * (cllbld0 / cll_subpop3_base)

    # ------------------------------------------------------------------
    # 4. ODE system (Data S1 run8634_eff). All four CLL states are carried
    #    on the cell-count scale; the SPD and spleen-volume observations
    #    are obtained by scaling at the observation step (equations 2-3),
    #    so no cell<->SPD conversion appears inside the ODEs.
    # ------------------------------------------------------------------
    d/dt(pbtk) <- kin_pbtk * (1 - effauc) - kout_pbtk * pbtk
    pbtk(0)    <- pbtk_base

    # stroma-attached colony 1 (fast detachment; slp2)
    d/dt(cll_subpop1) <- kp * (1 - eff_prol) * cll_subpop1 -
      kdtch * (1 + eff_dtch1) * cll_subpop1
    cll_subpop1(0)    <- clltiss0 * (1 - f2) * (1 - f1)

    # stroma-attached colony 2 (slow detachment; slp3)
    d/dt(cll_subpop2) <- kp * (1 - eff_prol) * cll_subpop2 -
      kdtch * (1 + eff_dtch2) * cll_subpop2
    cll_subpop2(0)    <- clltiss0 * (1 - f2) * f1

    # released / activated tissue pool that exits to blood
    d/dt(cll_subpop3) <- kp * (1 - eff_prol) * cll_subpop3 +
      kdtch * (1 + eff_dtch1) * cll_subpop1 +
      kdtch * (1 + eff_dtch2) * cll_subpop2 +
      kh * (1 - eff_home) * cll_bld -
      (kdist + eff_death) * cll_subpop3
    cll_subpop3(0)    <- clltiss0 * f2

    # resting CLL pool in peripheral blood
    d/dt(cll_bld) <- kdist * cll_subpop3 -
      kh * (1 - eff_home) * cll_bld -
      kd_bld * cll_bld
    cll_bld(0)    <- cllbld0

    # ------------------------------------------------------------------
    # 5. Observations (Ibrahim 2025 equations 2-5).
    # ------------------------------------------------------------------
    cll_tiss_total <- cll_subpop1 + cll_subpop2 + cll_subpop3

    # eq. 2 -- SPD (cm^2)
    tumorSpd <- cll_tiss_total * (spdbase / clltiss0) + spdnrm

    # eq. 3 -- spleen volume (cc)
    spleenVolume <- cll_tiss_total * (spleenbase / clltiss0) + spleennrm

    # eq. 4 -- leukocyte count (10^9 cells/L)
    leukocyte <- (cll_bld + leuknrm) / vbld

    # eq. 5 -- lymphocyte count (10^9 cells/L)
    lymphocyte <- (cll_bld + leuknrm * f_lymphocyte) / vbld

    # Total CLL burden on the cell-count scale, the quantity the paper's
    # response-depth criteria are evaluated against.
    totalCll <- cll_bld + cll_tiss_total

    # Btk occupancy (%), the target-engagement metric of the framework.
    btkOccupancy <- (pbtk_base - pbtk) * 100

    leukocyte    ~ lnorm(expSd_leukocyte)
    lymphocyte   ~ lnorm(expSd_lymphocyte)
    tumorSpd     ~ lnorm(expSd_tumorSpd)
    spleenVolume ~ lnorm(expSd_spleenVolume)
  })
}
