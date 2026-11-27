Schindler_2017_axitinib_dbp <- function() {
  description <- "Indirect-response model of axitinib-induced increase in diastolic blood pressure (dBP) during the first treatment cycle in Japanese adults with cytokine-refractory metastatic renal cell carcinoma (Schindler 2017). Daily axitinib AUC (AUC_AXITINIB, ug*h/L) stimulates dBP production through an Emax model parameterised by Emax and the initial slope S0 = Emax / AUC50 (AUC50 = 155 ug*h/L); the loss rate is 1 / MRT. Baseline dBP carries a Box-Cox-transformed random effect; residual error is additive. Fitted to first-cycle data only, because new antihypertensive treatment was prescribed after dBP elevation."
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
    "modellib('Schindler_2017_axitinib_biomarkers_sld'),",
    "modellib('Schindler_2017_axitinib_os'),",
    "modellib('Schindler_2017_axitinib_skit').",
    sep = " "
  )
  vignette <- "Schindler_2017_axitinib"
  units <- list(time = "day", dosing = "mg", concentration = "mmHg (diastolic blood pressure)")

  compartmentData <- list(
    dbp = list(analyte = "diastolic blood pressure", units = "mmHg", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    AUC_AXITINIB = list(
      description = "Daily axitinib area under the plasma concentration-time curve, AUCdaily = Dose_daily / (CL/F), with Dose_daily the total daily dose actually taken and CL/F the individual empirical Bayes apparent clearance from the published axitinib popPK model. Set to 0 on days off therapy.",
      units = "ug*h/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying (step function that changes with the daily dose). AUCdaily predicted dBP better than the daily dose (dOFV = 7.86, Schindler 2017 Results 'Diastolic BP model'). CL/F came from the Rini 2013 popPK model (modellib('Rini_2013_axitinib')). Supplied in the source dataset as the column AUC.",
      source_name = "AUC"
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
    dose_range = "Axitinib 5 mg PO twice daily starting dose (first 4-week cycle modelled).",
    regions = "Japan",
    biomarkers = "Diastolic blood pressure, n = 308 observations from cycle 1 (pre-dose, days 8, 15, 22) and day 1 of cycle 2; median follow-up 29 days (range 21-29).",
    notes = "Phase II study NCT00569946 (Tomita et al. Eur J Cancer 2011;47:2592-2602). Only first-cycle data were modelled because protocol-mandated antihypertensive treatment followed dBP elevation (Schindler 2017 Methods 'Diastolic blood pressure model'; Table 1 footnote c)."
  )

  ini({
    # Schindler 2017 Table 3 'Diastolic blood pressure model'. The
    # Supplement 2 dBP control stream carries the same values in $THETA.
    lrbase_dbp <- log(78.9); label("Typical baseline diastolic blood pressure dBP0 (mmHg)") # Table 3 'dBP0' = 78.9 mmHg (RSE 1.4%)
    boxcox_dbp <- -5.42; label("Box-Cox shape parameter of the baseline dBP random effect (unitless)") # Table 3 'ShapedBP0' = -5.42 (RSE 42%)
    lmrt_dbp <- log(4.92); label("Mean residence time of the dBP response, 1/kout (day)") # Table 3 'MRTdBP' = 4.92 days (RSE 19%)
    lemax_dbp <- log(0.197); label("Maximal fractional stimulation of dBP production, Emax (unitless)") # Table 3 'Emax,dBP' = 0.197 (RSE 14%)
    lslope_dbp <- log(0.00127); label("Initial slope of the Emax model, S0 = Emax / AUC50 (L/(ug*h))") # Table 3 'S0,dBP' = 0.00127 L*h^-1*ug^-1 (RSE 50%; SIR 95% CI 0.609-3.14 e-3); Supplement 2 S0BP = THETA(3)/1000 with THETA(3) = 1.27

    # IIV on baseline dBP: Table 3 6.7% CV; Supplement 2 $OMEGA 0.00453
    # (sqrt = 0.0673). No IIV on the drug-effect parameters.
    etalrbase_dbp ~ 0.00453 # Table 3 IIV 'dBP0' 6.7% CV; Supplement 2 $OMEGA 0.00453

    addSd_dbp <- 6.13; label("Additive residual error on dBP (mmHg)") # Table 3 'RUV (mmHg)' = 6.13 (RSE 7.0%); Supplement 2 $SIGMA 37.6 = 6.13^2
  })

  model({
    # Box-Cox transformed baseline random effect (Supplement 2 $PK):
    # PHI = EXP(ETA(1)); ETATR = (PHI^shape - 1) / shape; BASBP = THETA(1) * EXP(ETATR).
    phi_dbp <- exp(etalrbase_dbp)
    etatr_dbp <- (phi_dbp^boxcox_dbp - 1) / boxcox_dbp
    rbase_dbp <- exp(lrbase_dbp + etatr_dbp)

    mrt_dbp <- exp(lmrt_dbp)
    emax_dbp <- exp(lemax_dbp)
    slope_dbp <- exp(lslope_dbp)
    kout_dbp <- 1 / mrt_dbp
    kin_dbp <- rbase_dbp * kout_dbp

    # Eq. 4: Emax stimulation of production, parameterised by Emax and S0.
    eff_dbp <- emax_dbp * slope_dbp * AUC_AXITINIB / (emax_dbp + slope_dbp * AUC_AXITINIB)

    dbp(0) <- rbase_dbp
    d/dt(dbp) <- kin_dbp * (1 + eff_dbp) - kout_dbp * dbp

    dbp ~ add(addSd_dbp)
  })
}
