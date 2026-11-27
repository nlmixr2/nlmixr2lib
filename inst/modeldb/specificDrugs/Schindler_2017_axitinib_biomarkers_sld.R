Schindler_2017_axitinib_biomarkers_sld <- function() {
  description <- "Joint biomarker and tumour-size model for oral axitinib in Japanese adults with cytokine-refractory metastatic renal cell carcinoma (Schindler 2017). Daily axitinib AUC (AUC_AXITINIB, ug*h/L) drives four indirect-response biomarker models: Imax inhibition of VEGF degradation (with a linear drug-independent increase in VEGF production over time) and Imax / sigmoid-Imax inhibition of sVEGFR-1, sVEGFR-2 and sVEGFR-3 production (Imax fixed to 1; sVEGFR-2 and sVEGFR-3 share one AUC50; sVEGFR-1, -2 and -3 share one AUC50 random effect). Tumour size (sum of longest diameters, SLD) grows first-order and shrinks in proportion to the model-predicted relative change in sVEGFR-3 from baseline, with an exponentially waning (resistance) drug effect; the observed baseline SLD (TUM_SLD) is the initial condition with the residual-error-scaled B2 baseline random effect. A logistic dropout model for tumour assessments (per-day probability driven by progressive disease, current SLD, AUC and time) is carried as the algebraic output p_dropout_day for simulation. Biomarkers have log-normal residual error; SLD has proportional residual error with IIV on its magnitude."
  reference <- paste(
    "Schindler E, Amantea MA, Karlsson MO, Friberg LE.",
    "A Pharmacometric Framework for Axitinib Exposure, Efficacy, and Safety",
    "in Metastatic Renal Cell Carcinoma Patients.",
    "CPT Pharmacometrics Syst Pharmacol. 2017;6(6):373-382.",
    "doi:10.1002/psp4.12193.",
    "AUC_AXITINIB is Dose_daily / (CL/F) with the individual CL/F from the",
    "axitinib popPK model of Rini et al. J Clin Pharmacol 2013;53(5):491-504",
    "(doi:10.1002/jcph.73; modellib('Rini_2013_axitinib')).",
    "Sister model files from the same paper:",
    "modellib('Schindler_2017_axitinib_dbp'),",
    "modellib('Schindler_2017_axitinib_os'),",
    "modellib('Schindler_2017_axitinib_skit').",
    sep = " "
  )
  vignette <- "Schindler_2017_axitinib"
  paper_specific_compartments <- c("svegfr1", "tumor_nadir")
  units <- list(
    time = "day",
    dosing = "mg",
    concentration = "pg/mL (VEGF, sVEGFR-1, sVEGFR-2, sVEGFR-3); mm (tumour sum of longest diameters)"
  )

  compartmentData <- list(
    vegf = list(
      analyte = "vascular endothelial growth factor (VEGF)",
      units = "pg/mL",
      specimen = "plasma",
      verified = TRUE
    ),
    svegfr1 = list(
      analyte = "soluble VEGF receptor 1 (sVEGFR-1)",
      units = "pg/mL",
      specimen = "plasma",
      verified = TRUE
    ),
    svegfr2 = list(
      analyte = "soluble VEGF receptor 2 (sVEGFR-2)",
      units = "pg/mL",
      specimen = "plasma",
      verified = TRUE
    ),
    svegfr3 = list(
      analyte = "soluble VEGF receptor 3 (sVEGFR-3)",
      units = "pg/mL",
      specimen = "plasma",
      verified = TRUE
    ),
    tumor_size = list(
      analyte = "tumour sum of longest diameters (RECIST 1.0)",
      units = "mm",
      specimen = "not applicable",
      verified = TRUE
    ),
    tumor_nadir = list(
      analyte = "running minimum of the model-predicted tumour sum of longest diameters",
      units = "mm",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    AUC_AXITINIB = list(
      description = "Daily axitinib area under the plasma concentration-time curve, AUCdaily = Dose_daily / (CL/F), with Dose_daily the total daily dose actually taken (dose increases and reductions included) and CL/F the individual empirical Bayes apparent clearance from the published axitinib popPK model. Set to 0 on days off therapy.",
      units = "ug*h/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying (step function that changes with the daily dose). Schindler 2017 Methods 'Pharmacokinetics' sets AUCdaily to 0 on days off therapy because axitinib's typical half-life is short (3 h in Japanese patients). Observed range 31.95-1861 ug*h/L (Results 'Biomarkers models'). CL/F came from the Rini 2013 popPK model (modellib('Rini_2013_axitinib')). Supplied in the source datasets as the column AUC.",
      source_name = "AUC"
    ),
    TUM_SLD = list(
      description = "Observed baseline tumour sum of longest diameters (RECIST 1.0), per subject, time-fixed. Initial condition of the tumour-size state through the B2 baseline method.",
      units = "mm",
      type = "continuous",
      reference_category = NULL,
      notes = "Schindler 2017 Results 'Tumor size model': 'The observed baseline tumor size (SLD(0)) was included as a covariate ... acknowledging the same RUV as for postbaseline observations (B2 method)'. Supplement 2 control stream: IBASE = OBASE + ETA(4) * W1 with W1 = THETA(4) * EXP(ETA(5)) * OBASE and OMEGA(4) = 1 FIX. Supplied as the column BSLD.",
      source_name = "BSLD"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 64L,
    n_studies = 1L,
    age_range = "adults (age not reported in Schindler 2017)",
    weight_range = "not reported in Schindler 2017",
    sex_female_pct = NA_real_,
    race_ethnicity = c(Japanese = 100),
    disease_state = "Cytokine-refractory metastatic renal cell carcinoma (mRCC)",
    dose_range = "Axitinib 5 mg PO twice daily in 4-week cycles, titrated up by 2-3 mg twice daily to at most 10 mg twice daily (n = 5) or reduced to 2 mg twice daily (n = 41) on tolerability.",
    regions = "Japan",
    biomarkers = "VEGF, sVEGFR-1, sVEGFR-2, sVEGFR-3 (n = 436 observations each; median follow-up 168 days, range 32-624); tumour SLD (n = 476; median follow-up 337 days, range 36-731).",
    notes = "Single-arm, open-label, multicentre phase II study NCT00569946 (Tomita et al. Eur J Cancer 2011;47:2592-2602). Median treatment duration 51 weeks (range 1.7-104). Four sVEGFR-1 and four sVEGFR-3 values below the limit of quantification were omitted. Schindler 2017 Table 1 summarises assessments; baseline demographics are not tabulated in the paper."
  )

  ini({
    # ------------------------------------------------------------------
    # Joint biomarker model: Schindler 2017 Table 2. The same values appear
    # as FIX in the Supplement 2 tumour-size control stream, which fitted
    # SLD with the biomarker population parameters held at these estimates
    # (PPP&D approach, Methods 'Tumor size model').
    # ------------------------------------------------------------------
    lrbase_vegf <- log(65.0); label("Baseline VEGF (pg/mL)") # Table 2 'Base' VEGF = 65.0 (RSE 7.8%)
    lrbase_svegfr1 <- log(83.5); label("Baseline sVEGFR-1 (pg/mL)") # Table 2 'Base' sVEGFR-1 = 83.5 (RSE 2.9%)
    lrbase_svegfr2 <- log(8850); label("Baseline sVEGFR-2 (pg/mL)") # Table 2 'Base' sVEGFR-2 = 8,850 (RSE 2.8%)
    lrbase_svegfr3 <- log(19500); label("Baseline sVEGFR-3 (pg/mL)") # Table 2 'Base' sVEGFR-3 = 19,500 (RSE 6.5%)

    lmrt_vegf <- log(0.722); label("Mean residence time of VEGF, 1/kout (day)") # Table 2 'MRT' VEGF = 0.722 (RSE 25%)
    lmrt_svegfr1 <- log(0.624); label("Mean residence time of sVEGFR-1, 1/kout (day)") # Table 2 'MRT' sVEGFR-1 = 0.624 (RSE 69%; SIR 95% CI 0.0444-1.58)
    lmrt_svegfr2 <- log(19.7); label("Mean residence time of sVEGFR-2, 1/kout (day)") # Table 2 'MRT' sVEGFR-2 = 19.7 (RSE 17%)
    lmrt_svegfr3 <- log(5.76); label("Mean residence time of sVEGFR-3, 1/kout (day)") # Table 2 'MRT' sVEGFR-3 = 5.76 (RSE 12%)

    limax <- fixed(log(1)); label("Maximal fractional drug inhibition Imax, shared by all four biomarkers (unitless)") # Table 2 'Imax' = 1 FIX for all four biomarkers

    lec50_vegf <- log(354); label("AUCdaily giving half-maximal inhibition of VEGF degradation, AUC50 (ug*h/L)") # Table 2 'AUC50' VEGF = 354 (RSE 13%)
    lec50_svegfr1 <- log(1380); label("AUCdaily giving half-maximal inhibition of sVEGFR-1 production, AUC50 (ug*h/L)") # Table 2 'AUC50' sVEGFR-1 = 1,380 (RSE 13%)
    lec50_svegfr23 <- log(717); label("AUCdaily giving half-maximal inhibition of sVEGFR-2 and sVEGFR-3 production, shared AUC50 (ug*h/L)") # Table 2 'AUC50' sVEGFR-2 and sVEGFR-3 = 717 (RSE 8.6%), footnote c 'Common AUC50 parameter for sVEGFR-2 and 3'

    lhill_svegfr2 <- log(0.733); label("Hill coefficient of the sVEGFR-2 sigmoid Imax model (unitless)") # Table 2 'gamma' sVEGFR-2 = 0.733 (RSE 16%); 1 FIX for the other three biomarkers

    ldp_slope_vegf <- log(0.00178); label("Slope of the linear drug-independent increase in VEGF production (1/day)") # Table 2 'alpha' VEGF = 0.650 /year (RSE 28%); Supplement 2 SLOPE = THETA(18)/1000 with THETA(18) = 1.78, i.e. 0.00178 /day = 0.650 /year

    # IIV (biomarkers). Table 2 reports IIV as %CV = 100 * sqrt(omega^2);
    # the Supplement 2 $OMEGA variances reproduce every Table 2 value
    # under that convention (for example sqrt(0.184) = 0.429 -> 43%).
    etalrbase_vegf ~ 0.184 # Table 2 IIV 'Base' VEGF 43% CV; Supplement 2 $OMEGA 0.184
    etalrbase_svegfr1 ~ 0.0296 # Table 2 IIV 'Base' sVEGFR-1 17% CV; Supplement 2 $OMEGA 0.0296
    etalrbase_svegfr2 ~ 0.0232 # Table 2 IIV 'Base' sVEGFR-2 15% CV; Supplement 2 $OMEGA 0.0232
    etalrbase_svegfr3 ~ 0.237 # Table 2 IIV 'Base' sVEGFR-3 49% CV; Supplement 2 $OMEGA 0.237
    etalmrt_svegfr2 ~ 0.565 # Table 2 IIV 'MRT' sVEGFR-2 75% CV; Supplement 2 $OMEGA 0.565
    etalec50_vegf ~ 0.156 # Table 2 IIV 'AUC50' VEGF 39% CV; Supplement 2 $OMEGA 0.156
    etalec50_svegfr ~ 0.206 # Table 2 IIV 'AUC50' sVEGFR-1, -2, -3 45% CV, footnote b 'common variability term'; Supplement 2 $OMEGA 0.206 (ETA(9) on EC501, EC502 and EC503)
    etaldp_slope_vegf ~ 0.751 # Table 2 IIV 'alpha' VEGF 87% CV; Supplement 2 $OMEGA 0.751

    # Biomarker residual error, additive on the log scale (Table 2
    # footnote d 'Expressed as standard deviation on log-scale').
    expSd_vegf <- 0.376; label("Biomarker-specific log-scale residual SD for VEGF") # Table 2 'RUV' VEGF = 0.376 (RSE 5.9%); Supplement 2 $SIGMA 0.142 = 0.376^2
    expSd_svegfr1 <- 0.193; label("Biomarker-specific log-scale residual SD for sVEGFR-1") # Table 2 'RUV' sVEGFR-1 = 0.193 (RSE 5.3%); Supplement 2 $SIGMA 0.0371
    expSd_svegfr2 <- 0.162; label("Biomarker-specific log-scale residual SD for sVEGFR-2") # Table 2 'RUV' sVEGFR-2 = 0.162 (RSE 14%); Supplement 2 $SIGMA 0.0262
    expSd_svegfr3 <- 0.263; label("Biomarker-specific log-scale residual SD for sVEGFR-3") # Table 2 'RUV' sVEGFR-3 = 0.263 (RSE 6.5%); Supplement 2 $SIGMA 0.0689
    expSd_common <- 0.0593; label("Log-scale residual SD shared by the four biomarkers sampled at the same time") # Table 2 'Common RUV' = 0.0593 (RSE 26%), footnote e; Supplement 2 $SIGMA 0.00351

    # ------------------------------------------------------------------
    # Tumour-size model: Schindler 2017 Table 3 'Tumor size model'. Rate
    # constants are printed per week; model() converts to per day as the
    # Supplement 2 control stream does (KG = THETA(1)/1000/7, etc.).
    # ------------------------------------------------------------------
    lkg <- log(0.00361); label("First-order tumour growth rate constant KG (1/week)") # Table 3 'KG' = 0.00361 /week (RSE 1.8%); estimated with an informative NWPRI prior 3.74e-3 /week (6.64% RSE) from Claret 2015
    k_svegfr3 <- -0.174; label("Tumour-size reduction rate constant per unit relative change in sVEGFR-3, ksVEGFR-3 (1/week)") # Table 3 'ksVEGFR-3' = -0.174 /week (RSE 15%); negative because sVEGFR-3rel(t) is negative
    llambda <- log(0.101); label("Rate constant of the exponential decline of the drug effect (resistance appearance), lambda (1/week)") # Table 3 'lambda' = 0.101 /week (RSE 18%)
    propSd_tumor_size <- 0.105; label("Proportional residual error on tumour SLD (fraction)") # Table 3 'RUV (%)' = 10.5 (RSE 8.2%)

    # IIV (tumour). Same %CV = 100 * sqrt(omega^2) convention as Table 2.
    etalkg ~ 2.56 # Table 3 IIV 'KG' 160% CV -> 1.60^2
    etallambda ~ 0.5184 # Table 3 IIV 'lambda' 72% CV -> 0.72^2
    etapropSd_tumor_size ~ 0.1225 # Table 3 IIV 'RUV' 35% CV -> 0.35^2; Supplement 2 W = IPRED * THETA(4) * EXP(ETA(5))
    etarbase_tumor ~ fixed(1) # B2 baseline method: Supplement 2 IBASE = OBASE + ETA(4) * W1, $OMEGA 1 FIX

    # ------------------------------------------------------------------
    # Tumour-assessment dropout model: Schindler 2017 Table 3 'Dropout
    # model' and Supplementary Material 1 Equations S2-S4.
    # ------------------------------------------------------------------
    dropout_int <- -6.11; label("Intercept of the logistic per-day dropout model, theta0 (logit)") # Table 3 'theta0' = -6.11 (RSE 7.4%)
    e_pd_dropout <- 1.22; label("Effect of progressive disease (20% increase in SLD from nadir) on the logit of the per-day dropout probability") # Table 3 'thetaPD' = 1.22 (RSE 22%)
    e_sld_dropout <- 0.00282; label("Effect of SLD at the time of evaluation on the logit of the per-day dropout probability (1/mm)") # Table 3 'thetaSLD' = 0.00282 /mm (RSE 31%)
    e_auc_dropout <- -0.00529; label("Effect of AUC_AXITINIB on the logit of the per-day dropout probability (L/(ug*h))") # Table 3 'thetaAUC' = -0.00529 L*h^-1*ug^-1 (RSE 18%)
    e_time_dropout <- 0.00371; label("Effect of time since start of study on the logit of the per-day dropout probability (1/day)") # Table 3 'thetaTime' = 0.00371 /day (RSE 45%)
  })

  model({
    # Individual biomarker parameters (Supplement 2 $PK).
    imax <- exp(limax)
    rbase_vegf <- exp(lrbase_vegf + etalrbase_vegf)
    rbase_svegfr1 <- exp(lrbase_svegfr1 + etalrbase_svegfr1)
    rbase_svegfr2 <- exp(lrbase_svegfr2 + etalrbase_svegfr2)
    rbase_svegfr3 <- exp(lrbase_svegfr3 + etalrbase_svegfr3)
    mrt_vegf <- exp(lmrt_vegf)
    mrt_svegfr1 <- exp(lmrt_svegfr1)
    mrt_svegfr2 <- exp(lmrt_svegfr2 + etalmrt_svegfr2)
    mrt_svegfr3 <- exp(lmrt_svegfr3)
    ec50_vegf <- exp(lec50_vegf + etalec50_vegf)
    ec50_svegfr1 <- exp(lec50_svegfr1 + etalec50_svegfr)
    ec50_svegfr23 <- exp(lec50_svegfr23 + etalec50_svegfr)
    hill_svegfr2 <- exp(lhill_svegfr2)
    dp_slope_vegf <- exp(ldp_slope_vegf + etaldp_slope_vegf)

    kout_vegf <- 1 / mrt_vegf
    kout_svegfr1 <- 1 / mrt_svegfr1
    kout_svegfr2 <- 1 / mrt_svegfr2
    kout_svegfr3 <- 1 / mrt_svegfr3
    # VEGF production rises linearly with time (Eq. 1, Rin * (1 + alpha * t)).
    kin_vegf <- rbase_vegf * (1 + dp_slope_vegf * t) * kout_vegf
    kin_svegfr1 <- rbase_svegfr1 * kout_svegfr1
    kin_svegfr2 <- rbase_svegfr2 * kout_svegfr2
    kin_svegfr3 <- rbase_svegfr3 * kout_svegfr3

    # Drug effects (Eqs. 1-2; Supplement 2 $DES EFF, EFF1, EFF2, EFF3).
    eff_vegf <- imax * AUC_AXITINIB / (ec50_vegf + AUC_AXITINIB)
    eff_svegfr1 <- imax * AUC_AXITINIB / (ec50_svegfr1 + AUC_AXITINIB)
    eff_svegfr2 <- imax * AUC_AXITINIB^hill_svegfr2 /
      (ec50_svegfr23^hill_svegfr2 + AUC_AXITINIB^hill_svegfr2)
    eff_svegfr3 <- imax * AUC_AXITINIB / (ec50_svegfr23 + AUC_AXITINIB)

    vegf(0) <- rbase_vegf
    svegfr1(0) <- rbase_svegfr1
    svegfr2(0) <- rbase_svegfr2
    svegfr3(0) <- rbase_svegfr3
    # Eq. 1: axitinib inhibits VEGF loss; Eq. 2: axitinib inhibits sVEGFR production.
    d/dt(vegf) <- kin_vegf - kout_vegf * (1 - eff_vegf) * vegf
    d/dt(svegfr1) <- kin_svegfr1 * (1 - eff_svegfr1) - kout_svegfr1 * svegfr1
    d/dt(svegfr2) <- kin_svegfr2 * (1 - eff_svegfr2) - kout_svegfr2 * svegfr2
    d/dt(svegfr3) <- kin_svegfr3 * (1 - eff_svegfr3) - kout_svegfr3 * svegfr3

    # Tumour size (Eq. 3). Rate constants printed per week -> per day.
    kg <- exp(lkg + etalkg) / 7
    lambda <- exp(llambda + etallambda) / 7
    kd_svegfr3 <- k_svegfr3 / 7
    rsd_tumor_size <- propSd_tumor_size * exp(etapropSd_tumor_size)
    svegfr3_rel <- (svegfr3 - rbase_svegfr3) / rbase_svegfr3
    dsld <- kg * tumor_size - kd_svegfr3 * svegfr3_rel * exp(-lambda * t) * tumor_size
    # B2 baseline: the observed baseline carries the same residual error as
    # post-baseline observations (Supplement 2 IBASE = OBASE + ETA(4) * W1).
    tumor_size(0) <- TUM_SLD * (1 + etarbase_tumor * rsd_tumor_size)
    d/dt(tumor_size) <- dsld

    # Running nadir of the predicted SLD, for the progressive-disease flag
    # of the dropout model (PD = SLD >= 1.2 * nadir, Methods 'Tumor size model').
    # The nadir relaxes onto SLD at a fast rate (1000 /day) whenever SLD is
    # below it and is otherwise constant; a continuous right-hand side keeps
    # the solver from stalling on a switched derivative.
    tumor_nadir(0) <- TUM_SLD * (1 + etarbase_tumor * rsd_tumor_size)
    d/dt(tumor_nadir) <- -1000 * (tumor_nadir - tumor_size) * (tumor_size < tumor_nadir)

    # Logistic dropout (Supplementary Material 1 Eqs. S2-S3): per-day dropout
    # probability mu. Over an interval of dT days, P = 1 - (1 - mu)^dT (Eq. S4).
    pd_flag <- (tumor_size >= 1.2 * tumor_nadir)
    dropout_lp <- dropout_int + e_pd_dropout * pd_flag + e_sld_dropout * tumor_size +
      e_time_dropout * t + e_auc_dropout * AUC_AXITINIB
    p_dropout_day <- expit(dropout_lp)

    # Observations. Biomarkers were modelled as log(concentration) with a
    # biomarker-specific plus a common additive log-scale error; the common
    # term is shared by the four biomarkers sampled together (NONMEM L2), so
    # each marginal log-scale SD is sqrt(specific^2 + common^2).
    lsd_vegf <- sqrt(expSd_vegf^2 + expSd_common^2)
    lsd_svegfr1 <- sqrt(expSd_svegfr1^2 + expSd_common^2)
    lsd_svegfr2 <- sqrt(expSd_svegfr2^2 + expSd_common^2)
    lsd_svegfr3 <- sqrt(expSd_svegfr3^2 + expSd_common^2)
    vegf ~ lnorm(lsd_vegf)
    svegfr1 ~ lnorm(lsd_svegfr1)
    svegfr2 ~ lnorm(lsd_svegfr2)
    svegfr3 ~ lnorm(lsd_svegfr3)
    tumor_size ~ prop(rsd_tumor_size)
  })
}
