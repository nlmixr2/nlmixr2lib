Dolder_2017_lysergide <- function() {
  description <- paste(
    "One-compartment pharmacokinetic model of oral lysergide (LSD,",
    "d-lysergic acid diethylamide) in healthy adults (Dolder 2017), with",
    "first-order absorption without lag time and first-order elimination.",
    "Clearance and volume are apparent (CL/F, Vd/F); only the oral route was",
    "studied. Dolder 2017 fitted every subject individually in Phoenix",
    "WinNonlin (no population model) and reported the geometric means of the",
    "individual estimates separately for the 100 ug (n = 24) and 200 ug",
    "(n = 16) studies; DOSE_HIGH selects between the two parameter sets. The",
    "model carries no between-subject variability and the residual error,",
    "which was not reported, is fixed at zero. Doses are ug of LSD and",
    "concentrations ng/mL. The effect-compartment sigmoid-Emax PK-PD models",
    "are Dolder_2017_lysergide_anydrugeffect, _gooddrugeffect,",
    "_baddrugeffect, _hr, _bodytemp, _dbp and _sbp."
  )
  reference <- paste(
    "Dolder PC, Schmid Y, Steuer AE, Kraemer T, Rentsch KM, Hammann F, Liechti ME",
    "(2017). Pharmacokinetics and Pharmacodynamics of Lysergic Acid",
    "Diethylamide in Healthy Subjects. Clin Pharmacokinet 56(10):1219-1230.",
    "doi:10.1007/s40262-017-0513-9"
  )
  vignette <- "Dolder_2017_lysergide"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "lysergide", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "lysergide", units = "ug", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    DOSE_HIGH = list(
      description = "Dose group: 1 = 200 ug LSD (Study 2, NCT01878942), 0 = 100 ug LSD (Study 1, NCT02308969).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (100 ug single oral dose)",
      notes = paste(
        "Time-fixed per subject. Dolder 2017 fitted every subject individually",
        "and summarised the estimates separately for the two studies, which",
        "gave 100 ug (n = 24) and 200 ug (n = 16) to different participants",
        "and measured plasma LSD in different laboratories. The indicator",
        "selects the 100 ug or the 200 ug row of Table 1 (PK) and Table 2",
        "(PD); dose group and study are therefore confounded. Only the two",
        "studied doses are represented: the model does not interpolate",
        "between or extrapolate beyond them."
      ),
      source_name = "Dose (Table 1 and Table 2 row label: 100 ug / 200 ug)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 40L,
    n_studies = 2L,
    age_range = "25-60 years",
    age_mean = "33 +/- 11 years (100 ug study); 29 +/- 6 years (200 ug study)",
    weight_range = "52-98 kg",
    weight_mean = "68 +/- 8 kg (100 ug study); 72 +/- 12 kg (200 ug study)",
    sex_female_pct = 50,
    disease_state = "healthy volunteers",
    dose_range = "Single oral dose of 100 ug (n = 24) or 200 ug (n = 16) LSD (d-lysergic acid diethylamide hydrate capsules).",
    regions = "Switzerland (University Hospital Basel)",
    notes = paste(
      "Two placebo-controlled, double-blind, two-period crossover studies at",
      "the University Hospital Basel (washout of at least 7 days); 11 of 40",
      "had used a hallucinogen before. Plasma sampled to 24 h. Pharmacodynamic",
      "effects are differences from placebo in the same subject at the same",
      "time point."
    )
  )

  ini({
    # Table 1: geometric means (95% CI) of the individual one-compartment
    # Phoenix WinNonlin fits, reported separately for the 100 ug (n = 24) and
    # 200 ug (n = 16) studies. The elimination rate constant follows as
    # cl / vc = 12.3 / 46 = 0.267 1/h and 9.9 / 37 = 0.268 1/h (Table 1
    # lambda 0.27 1/h for both doses).
    lka <- fixed(log(1.4)); label("First-order absorption rate constant k01, 100 ug (1/h)") # Table 1 100 ug: k01 = 1.4 (1.2-4.1) 1/h
    lka_highdose <- fixed(log(1.2)); label("First-order absorption rate constant k01, 200 ug (1/h)") # Table 1 200 ug: k01 = 1.2 (0.68-4.6) 1/h
    lcl <- fixed(log(12.3)); label("Apparent clearance CL/F, 100 ug (L/h)") # Table 1 100 ug: CL/F = 12.3 (7.8-24) L/h
    lcl_highdose <- fixed(log(9.9)); label("Apparent clearance CL/F, 200 ug (L/h)") # Table 1 200 ug: CL/F = 9.9 (8.3-12.8) L/h
    lvc <- fixed(log(46)); label("Apparent volume of distribution Vd/F, 100 ug (L)") # Table 1 100 ug: Vd = 46 (35-76) L
    lvc_highdose <- fixed(log(37)); label("Apparent volume of distribution Vd/F, 200 ug (L)") # Table 1 200 ug: Vd = 37 (32-46) L

    propSd <- fixed(0); label("Proportional residual error on plasma LSD (fraction; not reported)") # Not reported; fixed to 0
  })

  model({
    # DOSE_HIGH selects the 200 ug row of Table 1 (and Table 2), 0 the 100 ug row.
    hi <- DOSE_HIGH
    lo <- 1 - DOSE_HIGH
    ka <- exp(lo * lka + hi * lka_highdose)
    cl <- exp(lo * lcl + hi * lcl_highdose)
    vc <- exp(lo * lvc + hi * lvc_highdose)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in ug and volume in L give ug/L = ng/mL, the reported unit.
    Cc <- central / vc

    Cc ~ prop(propSd)
  })
}
