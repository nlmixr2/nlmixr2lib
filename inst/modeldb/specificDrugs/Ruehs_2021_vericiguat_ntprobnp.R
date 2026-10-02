Ruehs_2021_vericiguat_ntprobnp <- function() {
  description <- paste(
    "PK/PD turnover model for the effect of oral vericiguat on plasma",
    "NT-proBNP in adults with worsening chronic heart failure and left",
    "ventricular ejection fraction < 45% (phase II SOCRATES-REDUCED). The PK",
    "layer is the paper's PK base model (one compartment, first-order",
    "absorption, correlated IIV on CL/F and V/F, no covariates), which the",
    "authors used to generate the exposure for the PK/PD analysis. NT-proBNP",
    "follows a turnover model whose zero-order production rises with the",
    "log of the individual baseline and is inhibited linearly by the 24-h",
    "vericiguat AUC, and whose first-order elimination is itself inhibited by",
    "NT-proBNP (Emax fixed at 0.95). Without drug the model is not at",
    "equilibrium at baseline, so NT-proBNP drifts downward on standard of care",
    "alone (the placebo arm). Residual error is additive on the log scale."
  )

  reference <- paste(
    "Ruehs H, Klein D, Frei M, Grevel J, Austin R, Becker C, Roessig L,",
    "Pieske B, Garmann D, Meyer M.",
    "Population Pharmacokinetics and Pharmacodynamics of Vericiguat in",
    "Patients with Heart Failure and Reduced Ejection Fraction.",
    "Clin Pharmacokinet. 2021;60(11):1407-1421.",
    "doi:10.1007/s40262-021-01024-y"
  )

  vignette <- "Ruehs_2021_vericiguat"

  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ug/L (vericiguat Cc); pg/mL (NT-proBNP)"
  )

  # NT-proBNP turnover state: first sighting of this PD state in the
  # library, so it is paper-specific until a second model of the same class
  # arrives.
  paper_specific_compartments <- c("ntprobnp")

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "vericiguat", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "vericiguat", units = "mg", specimen = "plasma", verified = TRUE),
    ntprobnp = list(analyte = "NT-proBNP", units = "pg/mL", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 432,
    n_studies = 1,
    age_median = "68 years (PK covariate reference value)",
    disease_state = paste(
      "Worsening chronic heart failure with left ventricular ejection",
      "fraction < 45%, on guideline-directed standard of care"
    ),
    dose_range = paste(
      "Oral vericiguat once daily for 12 weeks, target doses 1.25, 2.5, 5",
      "and 10 mg (the 5 and 10 mg arms started at 2.5 mg and were",
      "up-titrated at weeks 2 and 4), or placebo"
    ),
    regions = "Multinational (SOCRATES-REDUCED, NCT01951625)",
    baseline_ntprobnp = paste(
      "94.1-69,720 pg/mL; quartile boundaries 1559, 3000 and 6246 pg/mL",
      "(Ruehs 2021 Section 2.6)"
    ),
    notes = paste(
      "Ruehs 2021 Section 3.1 and Fig. 1: 2347 eligible NT-proBNP samples",
      "from 432 patients (all arms including placebo). The PK layer was",
      "developed on 454 patients (363 vericiguat-treated, 3376 samples)."
    )
  )

  ini({
    # PK base model, Ruehs 2021 Table 1 (no covariates, F = 1).
    lka <- log(1.5); label("Absorption rate constant ka (1/h)") # Table 1: ka = 1.5 1/h (RSE 9.53%)
    lcl <- log(1.3); label("Apparent clearance CL/F (L/h)") # Table 1: CL/F = 1.3 L/h (RSE 2.09%)
    lvc <- log(38.9); label("Apparent volume of distribution V/F (L)") # Table 1: V/F = 38.9 L (RSE 2.18%)

    # NT-proBNP turnover model, Ruehs 2021 Table 3 and Eqs. 4-6. Rate
    # constants are per day as printed and converted to per hour in model().
    lkin <- log(77.4); label("Typical NT-proBNP production rate TVkin at the median baseline (pg/mL/day)") # Table 3: TVkin = 77.4 pg/mL/day (bootstrap 95% CI 46.4-126.9)
    lrbase <- log(3140); label("Typical baseline NT-proBNP (pg/mL)") # Table 3: [NT-proBNP]baseline = 3140 pg/mL (95% CI 2857-3482)
    e_rbase_kin <- 0.347; label("Effect of log baseline NT-proBNP on kin (per log unit)") # Table 3: theta kin,NT-proBNP = 0.347 (95% CI 0.304-0.405); Eq. 6
    lkout_max <- log(0.157); label("Maximum first-order NT-proBNP elimination rate constant kout_max (1/day)") # Table 3: kout_max = 0.157 1/day (95% CI 0.086-0.330)
    emax <- fixed(0.95); label("Maximum fractional inhibition of kout by NT-proBNP (unitless)") # Table 3: Emax = 0.95 FIX
    lec50 <- log(439); label("NT-proBNP concentration giving half-maximal inhibition of kout (pg/mL)") # Table 3: EC50 = 439 pg/mL (95% CI 233-792)
    e_auc_kin <- 0.0176; label("Linear inhibition of kin by the 24-h vericiguat AUC, ATRT (L/h/mg)") # Table 3: ATRT = 0.0176 L/h/mg (printed CI 0.024-0.360 excludes the estimate)

    # PK IIV, Table 1: omega^2 = log(CV^2 + 1) from the printed CVs (footnote
    # a), CL/F 38.47% -> 0.13802, V/F 28.23% -> 0.076677, ka 102.99% ->
    # 0.72304; covariance = 0.70 x sqrt(0.13802 x 0.076677) = 0.07201.
    etalcl + etalvc ~ c(0.13802, 0.07201, 0.076677) # Table 1: IIV CL/F 38.47%, V/F 28.23%, correlation 0.70
    etalka ~ 0.72304 # Table 1: IIV ka CV 102.99%

    # PD IIV, Table 3 (variances as printed).
    etalkin ~ 0.163 # Table 3: omega^2 (IIV kin) = 0.163
    etalrbase ~ 0.953 # Table 3: omega^2 (IIV NT-proBNP baseline) = 0.953

    propSd <- 0.2602; label("Proportional residual error on vericiguat Cc (fraction)") # Table 1: proportional error 26.02% = SQRT(SIGMA^2) x 100
    addSd <- 7.81; label("Additive residual error on vericiguat Cc (ug/L)") # Table 1: additive error 7.81 ug/L = SQRT(SIGMA^2)
    expSd_ntprobnp <- 0.38079; label("Additive residual error on log NT-proBNP (SD, log scale)") # Table 3: sigma^2 = 0.145 on log-transformed NT-proBNP; SD = sqrt(0.145)
  })

  model({
    # PK base model
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Vericiguat exposure driving the PD: the 24-h AUC in mg*h/L (ATRT is
    # per mg*h/L). The paper updated a precomputed daily AUC at 24-h
    # intervals; here it is the continuous 24 x C(t), whose average over
    # a dosing interval equals the dosing-interval AUC.
    auc24 <- 24 * central / vc

    # NT-proBNP turnover (Eqs. 4-6). rbase is the individual baseline and the
    # initial condition; 8.03 is ln of the median baseline (~3070 pg/mL).
    rbase <- exp(lrbase + etalrbase)
    kin <- exp(lkin + etalkin) * (1 + e_rbase_kin * (log(rbase) - 8.03)) *
      (1 - e_auc_kin * auc24) / 24
    ec50 <- exp(lec50)
    kout <- exp(lkout_max) * (1 - emax * ntprobnp / (ec50 + ntprobnp)) / 24

    ntprobnp(0) <- rbase
    d/dt(ntprobnp) <- kin - kout * ntprobnp

    # Dose in mg, volume in L -> mg/L; x 1000 -> ug/L.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
    ntprobnp ~ lnorm(expSd_ntprobnp)
  })
}
