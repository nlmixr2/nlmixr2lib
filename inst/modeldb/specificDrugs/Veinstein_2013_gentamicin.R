Veinstein_2013_gentamicin <- function() {
  description <- "One-compartment population PK model for intravenous gentamicin in critically ill adult ICU patients with acute kidney injury undergoing 4-hour intermittent hemodialysis (n=10, all male; 6 mg/kg infused over 30 min, with hemodialysis starting 30 min after the end of the infusion; Veinstein 2013). Disposition is parameterised in terms of non-hemodialysis (interdialytic body) clearance, an additive hemodialysis-arm clearance, and volume of distribution. The dialysis arm is gated on/off by the time-varying RRT_HEMODIAL_ACTIVE covariate. Body weight enters the model as a linear (exponent = 1) structural scaler on all three parameters because the published Table 4 estimates are reported per kg; weight was tested as an explicit covariate on V and not retained."
  reference <- "Veinstein A, Venisse N, Badin J, Pinsard M, Robert R, Dupuis A. Gentamicin in hemodialyzed critical care patients: early dialysis after administration of a high dose should be considered. Antimicrob Agents Chemother. 2013;57(2):977-982. doi:10.1128/AAC.01762-12"
  vignette <- "Veinstein_2013_gentamicin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Veinstein 2013 Table 3: actual body weight 54.5-102 kg (mean 75.6 +/- 15.8 across the 10 subjects). Used as a linear (exponent = 1) structural scaler on CL_NHD, CL_HD, and V because the paper reports the population typical values in per-kg units (Table 4 footnotes a and b). The exponent is fixed at 1 because the paper's parameterisation literally normalises by body weight; the typical-value form lcl <- log(per-kg-value) plus cl <- exp(lcl + etalcl) * WT reproduces Table 4 exactly. Weight as an estimable covariate effect on V was tested in the model-building step and not retained (Results, Population PK/PD analysis: 'including weight or ideal body weight in the model combined as factors influencing V did not improve the model fit').",
      source_name        = "WT"
    ),
    RRT_HEMODIAL_ACTIVE = list(
      description        = "Hemodialysis-active indicator (1 during a dialysis session, 0 otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (interdialytic / no dialysis running)",
      notes              = "Time-varying within subject. Gates the additive hemodialysis-arm clearance cl_hemodialysis: the dialyzer contribution is added to the interdialytic body clearance only while a session is running. Veinstein 2013 protocol (Methods, Experimental design): a 4-h intermittent-hemodialysis session was started 30 min after the end of the 30-min 6 mg/kg gentamicin infusion, so RRT_HEMODIAL_ACTIVE = 1 from t = 1 h to t = 5 h after the start of the infusion. The hemodialysis apparatus was a Gambro AK 200 Ultra S with a Toray B3 polymethylmethacrylate dialyzer, blood flow 200-300 mL/min, mean session length 236 +/- 13 min. The estimated typical-value cl_hemodialysis lumps the dialyzer-mediated clearance into a single THETA rather than parameterising it as a Michaels-equation function of blood and dialysate flow rates (cf. Liesenfeld 2013 dabigatran).",
      source_name        = "DIAL"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 10L,
    n_studies        = 1L,
    age_range        = "approx 54-75 years (mean 64.5 +/- 10.1; Results, Patients' characteristics)",
    weight_range     = "54.5-102 kg (mean 72.7 +/- 16.4 on admission; mean 75.6 +/- 15.8 for the dialysed cohort in Table 3)",
    sex_female_pct   = 0,
    race_ethnicity   = "Not reported (single-centre medical ICU at the University Hospital of Poitiers, France)",
    disease_state    = "Critically ill adult ICU patients with acute kidney injury requiring intermittent hemodialysis and suffering from a community-acquired or nosocomial infection requiring treatment with gentamicin. Severity scores: SOFA 3-15 (median approx 11), SAPS II 34-65 (median approx 49). 8/10 received mechanical ventilation; 8/10 received pressor amines (Table 2). Infections included pyelonephritis (n=2), VAP / mediastinitis / fasciitis / lower-limb ischemia / septic thrombophlebitis (nosocomial, n=5), and angiocholitis / peritonitis (community-acquired, n=2). Pathogens identified included Escherichia coli, Serratia marcescens, MSSA, Pseudomonas aeruginosa, and Enterococcus faecium; MICs for gentamicin where available ranged 0.5-2 mg/L (Table 2).",
    dose_range       = "Gentamicin 6 mg/kg actual body weight, IV infusion over 30 min via syringe pump (Schering-Plough SAS), administered 30 min before the start of a 4-h intermittent hemodialysis session. Per-subject administered doses 300-600 mg with infused doses 360 +/- 99 mg on average (Table 3). Concurrent antibacterial therapy varied per indication (piperacillin-tazobactam, ceftriaxone, vancomycin, levofloxacin, amoxicillin, oxacillin, teicoplanin, metronidazole; Table 2).",
    regions          = "France (single-centre medical ICU, University Hospital of Poitiers)",
    renal_function   = "Acute kidney injury requiring intermittent hemodialysis; baseline renal function not separately quantified. Estimated post-dialysis non-hemodialysis clearance CL_NHD in the published cohort was 4.7-22.7 mL/min (Table 3), substantially below normal native CrCL.",
    notes            = "Baseline demographics from Veinstein 2013 Results (Patients' characteristics) and Tables 2-3. Hemodialysis apparatus: Gambro AK 200 Ultra S with a Toray B3 polymethylmethacrylate dialyzer (Toray Medical Co.); blood flow 200-300 mL/min (mean 283 +/- 20); mean session length 236 +/- 13 min. Bioanalytical method: cloned enzyme donor immunoassay (CEDIA; Microgenics / Thermo Scientific) on a modular Roche analyzer; LLOQ 0.24 ug/mL; between-run imprecision 2.1-4.0%."
  )

  ini({
    # Structural PK parameters (Veinstein 2013 Table 4, final-model column).
    # Typical values are reported per kg of body weight (Table 4 footnotes a
    # and b); weight enters the model as a linear (exponent = 1) structural
    # scaler in model(). Unit conversion mL/min/kg -> L/h/kg = x60/1000:
    #   CL_NHD = 0.1205 mL/min/kg = 0.0072300 L/h/kg
    #   CL_HD  = 0.955  mL/min/kg = 0.0573000 L/h/kg
    #   V      = 0.201  L/kg     (already in L/kg, no conversion)
    lcl              <- log(0.007230); label("Non-hemodialysis (interdialytic) body clearance CL_NHD (L/h/kg)") # Veinstein 2013 Table 4: CL_NHD = 0.1205 mL/min/kg (RSE 16%)
    lcl_hemodialysis <- log(0.057300); label("Hemodialysis-arm clearance CL_HD (L/h/kg)")                       # Veinstein 2013 Table 4: CL_HD  = 0.955  mL/min/kg (RSE 14%)
    lvc              <- log(0.201);    label("Volume of distribution V (L/kg)")                                 # Veinstein 2013 Table 4: V      = 0.201  L/kg      (RSE 12%)

    # Inter-individual variability (Veinstein 2013 Table 4 footnote c).
    # Exponential random-effect model, %CV reported on the natural scale.
    # Variance: omega^2 = log(CV^2 + 1).
    #   CL_NHD : 48% CV -> omega^2 = log(0.48^2 + 1) = 0.20734
    #   CL_HD  : 41% CV -> omega^2 = log(0.41^2 + 1) = 0.15538
    #   V      : 35% CV -> omega^2 = log(0.35^2 + 1) = 0.11556
    etalcl              ~ 0.20734  # Veinstein 2013 Table 4 (IIV CL_NHD, 48% CV; RSE 38%)
    etalcl_hemodialysis ~ 0.15538  # Veinstein 2013 Table 4 (IIV CL_HD,  41% CV; RSE 56%)
    etalvc              ~ 0.11556  # Veinstein 2013 Table 4 (IIV V,      35% CV; RSE 59%)

    # Residual error: proportional only (Veinstein 2013 Results, Population
    # PK/PD analysis: "A zero-order input one-compartment model and a
    # proportional-error model provided the best fit to the data and were
    # chosen as the best models").
    propSd <- 0.11; label("Proportional residual error (fraction)")  # Veinstein 2013 Table 4 (11% CV residual; RSE 40%)
  })

  model({
    # Individual PK parameters. The published Table 4 typical values are
    # reported per kg; multiplying by actual body weight returns absolute
    # L/h (clearances) and L (volume) for a given subject. The (WT)^1 scaler
    # is structural, not estimated as a separate exponent.
    cl              <- exp(lcl              + etalcl)              * WT  # L/h
    cl_hemodialysis <- exp(lcl_hemodialysis + etalcl_hemodialysis) * WT  # L/h
    vc              <- exp(lvc              + etalvc)              * WT  # L

    # Total apparent clearance (Veinstein 2013 Methods, Population PK
    # modeling equation: dC/dt = R0/V - [(CL_HD + CL_NHD)/V] * C during the
    # dialysis session, reducing to dC/dt = R0/V - [CL_NHD/V] * C
    # interdialytically). The dialysis arm cl_hemodialysis is added to the
    # body baseline cl only while a hemodialysis session is running,
    # encoded via the RRT_HEMODIAL_ACTIVE time-varying covariate.
    cl_total <- cl + RRT_HEMODIAL_ACTIVE * cl_hemodialysis

    kel <- cl_total / vc

    # One-compartment IV (zero-order input via the dosing infusion rate);
    # no depot. Dose is delivered into 'central' with a duration / rate set
    # by the input event.
    d/dt(central) <- -kel * central

    # Dose in mg, vc in L -> central / vc has units mg/L (matches paper
    # Table 3 concentrations reported in mg/liter).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
