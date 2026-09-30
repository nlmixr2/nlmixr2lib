Kloosterboer_2020_pipamperone <- function() {
  description <- paste(
    "One-compartment population pharmacokinetic model for oral pipamperone in",
    "children and adolescents (5.6-17.7 years) with behavioural problems,",
    "pooled from a Dutch prospective observational trial and a German",
    "therapeutic-drug-monitoring service (Kloosterboer 2020). First-order",
    "absorption with ka fixed at 2 /h; apparent clearance and apparent volume",
    "scale allometrically with body weight at fixed exponents 0.75 and 1",
    "referenced to 70 kg. Inter-patient variability on CL/F only. Combined",
    "additive + proportional residual error with an extra additive error for",
    "samples quantified by HPLC-UV, and a linear conversion of the plasma",
    "prediction onto the dried-blood-spot scale for DBS samples."
  )
  reference <- paste(
    "Kloosterboer SM, Egberts KM, de Winter BCM, van Gelder T, Gerlach M,",
    "Hillegers MHJ, Dieleman GC, Bahmany S, Reichart CG, van Daalen E,",
    "Kouijzer MEJ, Dierckx B, Koch BCP (2020). Pipamperone population",
    "pharmacokinetics related to effectiveness and side effects in children",
    "and adolescents. Clin Pharmacokinet 59(11):1393-1405.",
    "doi:10.1007/s40262-020-00894-y"
  )
  vignette <- "Kloosterboer_2020_pipamperone"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  compartmentData <- list(
    depot = list(analyte = "pipamperone", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "pipamperone", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight; allometric scaling of CL/F (exponent 0.75) and V/F (exponent 1), referenced to 70 kg",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Kloosterboer 2020 Section 2.4.1 and Table 2 footnote b: allometric",
        "scaling with fixed exponents 1 for V and 0.75 for CL; parameter",
        "estimates are reported per 70 kg. Model-building cohort median",
        "50.4 kg (range 24.8-100.4; Table 1). Weight was recorded at the time",
        "of blood sampling in the Dutch trial (time-varying allowed)."
      ),
      source_name = "bodyweight"
    ),
    ASSAY_LCMSMS = list(
      description = "Bioanalytical method indicator: 1 = liquid chromatography-mass spectrometry (Dutch samples), 0 = HPLC-UV (German TDM samples)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Per-sample. Direction: 1 = LC-MS (the Dutch UHPLC-MS plasma and DBS",
        "methods, LLOQ 1.5 ug/L; Section 2.2), 0 = HPLC-UV (German serum/plasma",
        "method, LLOQ 8 ug/L, linear 2-1050 ug/L). When 0 the extra HPLC-UV",
        "additive residual error (26.6 ug/L; Table 2) is added in variance to",
        "the common additive error. Affects the residual error only."
      ),
      source_name = "analytical method (LCMS vs HPLC-UV)"
    ),
    SAMPLE_DBS = list(
      description = "Dried-blood-spot sample indicator: 1 = DBS, 0 = venous plasma/serum",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Per-sample. Only Dutch samples were collected by DBS (Section 2.2).",
        "When 1 the observed DBS concentration is predicted as",
        "cal_slope_dbs * Cplasma + cal_int_dbs (Table 2, 'DBS correction:",
        "y = ax + b'). Whole-blood venous samples were not collected, so",
        "SAMPLE_WHOLEBLOOD is not used."
      ),
      source_name = "sampling method (DBS vs venepuncture)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 30,
    n_studies = 2,
    n_observations = 70,
    age_range = "5.6-17.7 years (median 13.0)",
    weight_range = "24.8-100.4 kg (median 50.4)",
    sex_female_pct = 30,
    disease_state = paste(
      "Children and adolescents treated with pipamperone for behavioural",
      "problems: autism spectrum disorder 73.3%, ADHD 40%, mental retardation",
      "30%, schizophrenia-spectrum 6.7%, conduct disorder 3.3% (Table 1;",
      "more than one diagnosis possible)."
    ),
    dose_range = "Oral tablet or oral solution, 12-400 mg/day (median 45 mg/day), once to five times daily",
    regions = "Netherlands (SPACe trial, n = 8) and Germany (Wuerzburg TDM service, n = 22)",
    notes = paste(
      "Model-building group of Kloosterboer 2020 Table 1. An external",
      "validation group of 21 German TDM patients (33 concentrations) was not",
      "used for estimation. Psychiatric comedication was common (60% other",
      "antipsychotics)."
    )
  )

  ini({
    lka <- fixed(log(2)); label("Absorption rate constant (1/h)") # Table 2, 'Ka' = 2, footnote a 'Fixed'; Section 2.4.1 'fixed at 2/h, based on the previous literature [21]' (Table 2 prints the unit as 'L/h', a typo for 1/h)
    lcl <- log(22.1); label("Apparent clearance CL/F for a 70 kg patient (L/h)") # Table 2, 'CL/F (L/h/70 kg)' = 22.1 (RSE 12%)
    lvc <- log(416); label("Apparent volume of distribution V/F for a 70 kg patient (L)") # Table 2, 'V/F (L/70 kg)' = 416 (RSE 32%)
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL/F (unitless)") # Table 2 footnote b, 'exponent ... 0.75 for CL'; Section 3.1.2 'fixed exponents'
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on V/F (unitless)") # Table 2 footnote b, 'exponent 1 for V'

    etalcl ~ 0.041153 # Table 2, 'IPV CL' = 20.5%; log(1 + 0.205^2) = 0.041153 (shrinkage 34%)

    cal_slope_dbs <- 0.33; label("Slope a of the plasma-to-DBS conversion y = a*x + b (unitless)") # Table 2, 'DBS correction: y = ax + b', 'a' = 0.33 (RSE 8%)
    cal_int_dbs <- 3.90; label("Intercept b of the plasma-to-DBS conversion y = a*x + b (ug/L)") # Table 2, 'DBS correction: y = ax + b', 'b (ug/L)' = 3.90 (RSE 15%)

    addSd <- 0.21; label("Additive residual error (ug/L)") # Table 2, 'Additional error (ug/L)' = 0.21 (RSE 1%)
    propSd <- 0.39; label("Proportional residual error (fraction)") # Table 2, 'Proportional error' = 0.39 (RSE 19%)
    addSd_hplc <- 26.6; label("Extra additive residual error for HPLC-UV samples (ug/L)") # Table 2, 'Additional error HPLC-UV (ug/L)' = 26.6 (RSE 40%); Section 3.1.3
  })
  model({
    # Individual PK parameters; allometric scaling on body weight (70 kg reference)
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc) * (WT / 70)^e_wt_vc
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Plasma concentration: dose in mg, volume in L -> mg/L; x1000 -> ug/L
    Cplasma <- 1000 * central / vc
    # DBS samples are predicted on the DBS scale via the linear conversion of Table 2
    Cc <- (1 - SAMPLE_DBS) * Cplasma + SAMPLE_DBS * (cal_slope_dbs * Cplasma + cal_int_dbs)

    # Extra additive error for HPLC-UV samples, added in variance to the common additive error
    addSdCc <- sqrt(addSd^2 + (1 - ASSAY_LCMSMS) * addSd_hplc^2)
    Cc ~ add(addSdCc) + prop(propSd)
  })
}
