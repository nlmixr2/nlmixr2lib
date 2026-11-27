Schindler_2017_axitinib_os <- function() {
  description <- "Parametric overall-survival (OS) model for Japanese adults with cytokine-refractory metastatic renal cell carcinoma on oral axitinib (Schindler 2017). The hazard of death has a log-logistic baseline (scale beta0, shape gamma) multiplied by exp(betaSLD * SLD(t)), where SLD(t) is the model-predicted tumour sum of longest diameters. SLD(t) is simulated in-model exactly as in the authors' OS control stream: daily axitinib AUC (AUC_AXITINIB, ug*h/L) inhibits sVEGFR-3 production (Imax model), and the relative change in sVEGFR-3 from baseline shrinks a first-order-growing tumour with an exponentially waning effect. The sVEGFR-3 and tumour-size parameters and their IIV are fixed to the final estimates of the joint biomarker / tumour-size model (modellib('Schindler_2017_axitinib_biomarkers_sld')); the OS fit took them as individual empirical Bayes estimates. A competing log-logistic hazard describes censoring. Outputs are the survival probabilities sur_os and sur_cens."
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
    "modellib('Schindler_2017_axitinib_dbp'),",
    "modellib('Schindler_2017_axitinib_skit').",
    sep = " "
  )
  vignette <- "Schindler_2017_axitinib"
  paper_specific_compartments <- c("cumhaz_os", "cumhaz_cens")
  units <- list(
    time = "day",
    dosing = "mg",
    concentration = "probability (sur_os and sur_cens are survival probabilities); pg/mL (sVEGFR-3); mm (tumour SLD)"
  )

  compartmentData <- list(
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
    cumhaz_os = list(
      analyte = "cumulative hazard of death",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    ),
    cumhaz_cens = list(
      analyte = "cumulative hazard of censoring",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    AUC_AXITINIB = list(
      description = "Daily axitinib area under the plasma concentration-time curve, AUCdaily = Dose_daily / (CL/F), with Dose_daily the total daily dose and CL/F the individual empirical Bayes apparent clearance from the published axitinib popPK model. Set to 0 on days off therapy.",
      units = "ug*h/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying. For the OS analysis the dose was extrapolated from the last recorded dose, assuming treatment continued until death or censoring (Schindler 2017 Methods 'Overall survival model'). CL/F came from the Rini 2013 popPK model (modellib('Rini_2013_axitinib')). Supplied in the source dataset as the column AUC.",
      source_name = "AUC"
    ),
    TUM_SLD = list(
      description = "Baseline tumour sum of longest diameters (RECIST 1.0), per subject, time-fixed; the initial condition of the tumour-size state.",
      units = "mm",
      type = "continuous",
      reference_category = NULL,
      notes = "The Supplement 2 OS control stream initialises the tumour state at the individual empirical Bayes baseline SLD (IBASS) from the joint biomarker / tumour-size fit. That fit used the observed baseline with a residual-error-sized random effect (B2 method), so the observed baseline is used here; supply an individual model-predicted baseline instead when one is available.",
      source_name = "IBASS"
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
    dose_range = "Axitinib 5 mg PO twice daily in 4-week cycles, titrated 2-10 mg twice daily on tolerability.",
    regions = "Japan",
    biomarkers = "Overall survival: 16 deaths and 48 censored among 64 patients; median follow-up 457 days (range 85-781).",
    notes = "Phase II study NCT00569946 (Tomita et al. Eur J Cancer 2011;47:2592-2602). OS follow-up was limited to the axitinib treatment period (maximum 112 weeks) to avoid confounding by subsequent therapies (Schindler 2017 Discussion)."
  )

  ini({
    # ------------------------------------------------------------------
    # sVEGFR-3 and tumour-size sub-models: final estimates of the joint
    # biomarker / tumour-size model (Schindler 2017 Tables 2-3), held fixed
    # in the OS fit (the Supplement 2 OS stream reads them as individual
    # EBEs IBAS3, IEC53, IBASS, IKG, ILAM, with MRT3 = 5.76 and
    # KD3 = 174/7/1000 hard-coded).
    # ------------------------------------------------------------------
    lrbase_svegfr3 <- fixed(log(19500)); label("Baseline sVEGFR-3 (pg/mL)") # Table 2 'Base' sVEGFR-3 = 19,500
    lmrt_svegfr3 <- fixed(log(5.76)); label("Mean residence time of sVEGFR-3, 1/kout (day)") # Table 2 'MRT' sVEGFR-3 = 5.76; Supplement 2 OS stream MRT3 = 5.76
    limax <- fixed(log(1)); label("Maximal fractional inhibition of sVEGFR-3 production, Imax (unitless)") # Table 2 'Imax' = 1 FIX; Supplement 2 OS stream EMAX = 1
    lec50_svegfr23 <- fixed(log(717)); label("AUCdaily giving half-maximal inhibition of sVEGFR-3 production, AUC50 (ug*h/L)") # Table 2 'AUC50' sVEGFR-2 and -3 = 717
    lkg <- fixed(log(0.00361)); label("First-order tumour growth rate constant KG (1/week)") # Table 3 'KG' = 0.00361 /week
    k_svegfr3 <- fixed(-0.174); label("Tumour-size reduction rate constant per unit relative change in sVEGFR-3, ksVEGFR-3 (1/week)") # Table 3 'ksVEGFR-3' = -0.174 /week; Supplement 2 OS stream KD3 = 174/7/1000 /day with the sign carried by -(-CFBV3*KD3)
    llambda <- fixed(log(0.101)); label("Rate constant of the exponential decline of the drug effect (resistance appearance), lambda (1/week)") # Table 3 'lambda' = 0.101 /week

    etalrbase_svegfr3 ~ fixed(0.237) # Table 2 IIV 'Base' sVEGFR-3 49% CV; Supplement 2 $OMEGA 0.237
    etalec50_svegfr ~ fixed(0.206) # Table 2 IIV 'AUC50' sVEGFR-1, -2, -3 45% CV; Supplement 2 $OMEGA 0.206
    etalkg ~ fixed(2.56) # Table 3 IIV 'KG' 160% CV -> 1.60^2
    etallambda ~ fixed(0.5184) # Table 3 IIV 'lambda' 72% CV -> 0.72^2

    # ------------------------------------------------------------------
    # Overall survival: Schindler 2017 Table 3 'Overall survival model',
    # Eqs. 5-6 (time in days, as in the Supplement 2 OS stream).
    # ------------------------------------------------------------------
    beta0_os <- 7.09; label("Scale parameter of the log-logistic baseline hazard of death, beta0 (psi = exp(-beta0), 1/day)") # Table 3 'beta0' = 7.09 (RSE 3.2%)
    lgamma_os <- log(0.298); label("Shape parameter of the log-logistic baseline hazard of death, gamma (unitless)") # Table 3 'gamma' = 0.298 (RSE 22%)
    e_sld_os <- 0.0115; label("Coefficient of the model-predicted SLD(t) on the log hazard of death (1/mm)") # Table 3 'betaSLD' = 0.0115 /mm (RSE 17%)

    # Censoring hazard: values printed only in the Supplement 2 OS control
    # stream $THETA (Results 'Overall survival model' reports only that a
    # competing log-logistic function described censoring).
    beta0_cens <- 6.21; label("Scale parameter of the log-logistic censoring hazard, beta0 (psi = exp(-beta0), 1/day)") # Supplement 2 OS stream $THETA (0,6.21) 'Beta0_censoring'
    lgamma_cens <- log(0.162); label("Shape parameter of the log-logistic censoring hazard, gamma (unitless)") # Supplement 2 OS stream $THETA (0,0.162) 'Gamma_censoring'
  })

  model({
    # sVEGFR-3 (Supplement 2 OS stream $PK / $DES).
    imax <- exp(limax)
    rbase_svegfr3 <- exp(lrbase_svegfr3 + etalrbase_svegfr3)
    mrt_svegfr3 <- exp(lmrt_svegfr3)
    ec50_svegfr23 <- exp(lec50_svegfr23 + etalec50_svegfr)
    kout_svegfr3 <- 1 / mrt_svegfr3
    kin_svegfr3 <- rbase_svegfr3 * kout_svegfr3
    eff_svegfr3 <- imax * AUC_AXITINIB / (ec50_svegfr23 + AUC_AXITINIB)
    svegfr3(0) <- rbase_svegfr3
    d/dt(svegfr3) <- kin_svegfr3 * (1 - eff_svegfr3) - kout_svegfr3 * svegfr3

    # Tumour size (Eq. 3); rate constants printed per week -> per day.
    kg <- exp(lkg + etalkg) / 7
    lambda <- exp(llambda + etallambda) / 7
    kd_svegfr3 <- k_svegfr3 / 7
    svegfr3_rel <- (svegfr3 - rbase_svegfr3) / rbase_svegfr3
    tumor_size(0) <- TUM_SLD
    d/dt(tumor_size) <- kg * tumor_size - kd_svegfr3 * svegfr3_rel * exp(-lambda * t) * tumor_size

    # Log-logistic hazards (Eqs. 5-6; Supplement 2 OS stream DADT(3), DADT(4)).
    psi_os <- exp(-beta0_os)
    gamma_os <- exp(lgamma_os)
    psi_cens <- exp(-beta0_cens)
    gamma_cens <- exp(lgamma_cens)
    hazard_os <- psi_os^(1 / gamma_os) * t^(1 / gamma_os - 1) /
      (gamma_os * (1 + (psi_os * t)^(1 / gamma_os))) * exp(e_sld_os * tumor_size)
    hazard_cens <- psi_cens^(1 / gamma_cens) * t^(1 / gamma_cens - 1) /
      (gamma_cens * (1 + (psi_cens * t)^(1 / gamma_cens)))

    cumhaz_os(0) <- 0
    cumhaz_cens(0) <- 0
    # Numerical guard, not part of the published model: in the high-KG tail of
    # the IIV the tumour grows without bound once the drug effect wanes, and
    # exp(betaSLD * SLD) overflows the solver. The integrated hazard is capped
    # at 1e6 /day; at that rate survival falls by a factor exp(-1000), i.e. to
    # exactly 0 in double precision, within 0.001 day, so sur_os is unchanged.
    d/dt(cumhaz_os) <- min(hazard_os, 1e6)
    d/dt(cumhaz_cens) <- hazard_cens
    sur_os <- exp(-cumhaz_os)
    sur_cens <- exp(-cumhaz_cens)
  })
}
