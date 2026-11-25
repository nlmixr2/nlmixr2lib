Dolder_2017_lysergide_dbp <- function() {
  description <- paste(
    "Effect-compartment sigmoid-Emax PK-PD model of",
    "the increase in diastolic blood pressure (mmHg) over placebo after oral lysergide (LSD) in healthy",
    "adults (Dolder 2017). The PK is the one-compartment model with",
    "first-order absorption and elimination of Dolder_2017_lysergide; the",
    "effect-site concentration follows plasma with the first-order",
    "equilibration rate constant ke0, and the effect is",
    "Emax * Ce^gamma / (EC50^gamma + Ce^gamma) with zero baseline, the",
    "effect being the difference from placebo at the same time point.",
    "Dolder 2017 fitted every subject individually in Phoenix WinNonlin",
    "with the individual PK frozen (no population model) and reported the",
    "PD parameters as arithmetic means of the individual estimates and the",
    "PK parameters as geometric means, separately for the 100 ug and",
    "200 ug studies; DOSE_HIGH selects between the two parameter sets.",
    "The model carries no between-subject variability and the residual",
    "errors, which were not reported, are fixed at zero."
  )
  reference <- paste(
    "Dolder PC, Schmid Y, Steuer AE, Kraemer T, Rentsch KM, Hammann F, Liechti ME",
    "(2017). Pharmacokinetics and Pharmacodynamics of Lysergic Acid",
    "Diethylamide in Healthy Subjects. Clin Pharmacokinet 56(10):1219-1230.",
    "doi:10.1007/s40262-017-0513-9"
  )
  vignette <- "Dolder_2017_lysergide"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL", dbpcfb = "mmHg, difference from placebo")

  compartmentData <- list(
    depot = list(analyte = "lysergide", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "lysergide", units = "ug", specimen = "plasma", verified = TRUE),
    effect = list(analyte = "lysergide", units = "ng/mL", specimen = "not applicable", verified = TRUE)
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
      "time point.",
      "Diastolic blood pressure did not change enough to estimate the PD",
      "parameters in one to three subjects (Results)."
    )
  )

  ini({
    # PK layer: the Table 1 geometric means, as in Dolder_2017_lysergide. The
    # individual PK estimates were held fixed for the PK-PD fits (Methods).
    lka <- fixed(log(1.4)); label("First-order absorption rate constant k01, 100 ug (1/h)") # Table 1 100 ug: k01 = 1.4 (1.2-4.1) 1/h
    lka_highdose <- fixed(log(1.2)); label("First-order absorption rate constant k01, 200 ug (1/h)") # Table 1 200 ug: k01 = 1.2 (0.68-4.6) 1/h
    lcl <- fixed(log(12.3)); label("Apparent clearance CL/F, 100 ug (L/h)") # Table 1 100 ug: CL/F = 12.3 (7.8-24) L/h
    lcl_highdose <- fixed(log(9.9)); label("Apparent clearance CL/F, 200 ug (L/h)") # Table 1 200 ug: CL/F = 9.9 (8.3-12.8) L/h
    lvc <- fixed(log(46)); label("Apparent volume of distribution Vd/F, 100 ug (L)") # Table 1 100 ug: Vd = 46 (35-76) L
    lvc_highdose <- fixed(log(37)); label("Apparent volume of distribution Vd/F, 200 ug (L)") # Table 1 200 ug: Vd = 37 (32-46) L

    propSd <- fixed(0); label("Proportional residual error on plasma LSD (fraction; not reported)") # Not reported; fixed to 0

    # PD layer: Table 2 'Diastolic blood pressure increase' rows, arithmetic means +/- SD of the individual
    # estimates (PK-PD link model). Emax upper limit 50 mmHg in the fits (Methods).
    emax <- fixed(23); label("Maximal effect Emax, 100 ug (mmHg)") # Table 2 'Diastolic blood pressure increase' 100 ug: Emax = 23 +/- 14
    emax_highdose <- fixed(18); label("Maximal effect Emax, 200 ug (mmHg)") # Table 2 'Diastolic blood pressure increase' 200 ug: Emax = 18 +/- 11
    lec50 <- fixed(log(0.9)); label("Effect-site concentration at half-maximal effect EC50, 100 ug (ng/mL)") # Table 2 'Diastolic blood pressure increase' 100 ug: EC50 = 0.9 +/- 0.6 ng/mL
    lec50_highdose <- fixed(log(1.6)); label("Effect-site concentration at half-maximal effect EC50, 200 ug (ng/mL)") # Table 2 'Diastolic blood pressure increase' 200 ug: EC50 = 1.6 +/- 0.9 ng/mL
    lhill <- fixed(log(2)); label("Sigmoidicity (Hill) coefficient gamma, 100 ug (unitless)") # Table 2 'Diastolic blood pressure increase' 100 ug: gamma = 2.0 +/- 1.6
    lhill_highdose <- fixed(log(3.5)); label("Sigmoidicity (Hill) coefficient gamma, 200 ug (unitless)") # Table 2 'Diastolic blood pressure increase' 200 ug: gamma = 3.5 +/- 1.6
    lke0 <- fixed(log(2.6)); label("Effect-compartment equilibration rate constant ke0, 100 ug (1/h)") # Table 2 'Diastolic blood pressure increase' 100 ug: ke0 = 2.6 +/- 1.9 1/h
    lke0_highdose <- fixed(log(3.4)); label("Effect-compartment equilibration rate constant ke0, 200 ug (1/h)") # Table 2 'Diastolic blood pressure increase' 200 ug: ke0 = 3.4 +/- 1.9 1/h

    addSd_dbpcfb <- fixed(0); label("Additive residual error on dbpcfb (mmHg; not reported)") # Not reported; fixed to 0
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

    emax_i <- lo * emax + hi * emax_highdose
    ec50 <- exp(lo * lec50 + hi * lec50_highdose)
    hill <- exp(lo * lhill + hi * lhill_highdose)
    ke0 <- exp(lo * lke0 + hi * lke0_highdose)

    # Effect site (ng/mL): dCe/dt = ke0 * (Cp - Ce) (Supplementary Figure S1).
    d/dt(effect) <- ke0 * (Cc - effect)

    dbpcfb <- emax_i * effect^hill / (ec50^hill + effect^hill)

    Cc ~ prop(propSd)
    dbpcfb ~ add(addSd_dbpcfb)
  })
}
