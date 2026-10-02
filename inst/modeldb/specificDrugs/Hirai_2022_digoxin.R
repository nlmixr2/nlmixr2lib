Hirai_2022_digoxin <- function() {
  description <- paste0(
    "One-compartment population PK model with first-order absorption for ",
    "oral digoxin in 391 Japanese adults with atrial fibrillation and heart ",
    "failure, fitted to routine steady-state trough serum concentrations ",
    "(Hirai 2022). The absorption rate constant (1.0 1/h) and the apparent ",
    "volume of distribution (6.0 L/kg, scaled linearly by body weight) were ",
    "fixed from the literature; only the apparent oral clearance was ",
    "estimated. CL/F scales as a power of Cockcroft-Gault creatinine ",
    "clearance normalised to 60 mL/min (capped at 120 mL/min) and falls by a ",
    "fractional 23.8% with concurrent amiodarone. Exponential ",
    "between-subject variability on CL/F and a multiplicative (proportional) ",
    "residual error. Companion Japanese digoxin trough model: ",
    "Komatsu_2015_digoxin."
  )
  reference <- paste(
    "Hirai T, Kasai H, Naganuma M, Hagiwara N, Shiga T.",
    "Population pharmacokinetic analysis and dosage recommendations for",
    "digoxin in Japanese patients with atrial fibrillation and heart failure",
    "using real-world data.",
    "BMC Pharmacol Toxicol. 2022;23:14.",
    "doi:10.1186/s40360-022-00552-y."
  )
  vignette <- "Hirai_2022_digoxin"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL"
  )

  covariateData <- list(
    CRCL = list(
      description = paste0(
        "Creatinine clearance estimated with the Cockcroft-Gault equation ",
        "(Hirai 2022 Methods, Data collection), in raw mL/min and NOT ",
        "BSA-normalised."
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Normalised to 60 mL/min in the power model (CRCL/60)^0.409 ",
        "(Hirai 2022 Results final-model equation). Values above 120 mL/min ",
        "were replaced with 120 mL/min to avoid overestimating renal ",
        "clearance (Methods, Population pharmacokinetic model development); ",
        "the cap is applied inside model(), so the column should carry the ",
        "uncapped Cockcroft-Gault value. Cohort median 56.5 [IQR 40.7-75.6] ",
        "mL/min (Table 1)."
      ),
      source_name = "CLCR"
    ),
    CONMED_AMIO = list(
      description = paste0(
        "1 = concurrent amiodarone, 0 = no amiodarone. 63 of 391 patients ",
        "(16%) were on amiodarone (Hirai 2022 Table 1)."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concurrent amiodarone)",
      notes = paste0(
        "Same polarity as the source ('if amiodarone' = 1). Enters CL/F as ",
        "the fractional factor (1 + e_amio_cl * CONMED_AMIO) with ",
        "e_amio_cl = -0.238 (Table 3), i.e. amiodarone lowers CL/F by 23.8%. ",
        "Amiodarone dose was not modelled (Discussion)."
      ),
      source_name = "amiodarone"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Used only to scale the literature-fixed apparent volume of ",
        "distribution linearly, Vd/F = 6.0 L/kg x body weight (Hirai 2022 ",
        "Results final-model equation). Weight has no effect on CL/F in this ",
        "model. Cohort mean 57 +/- 15 kg (Table 1)."
      ),
      source_name = "Body weight"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "digoxin",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "digoxin",
      units = "mg",
      specimen = "serum",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 391L,
    n_studies = 1L,
    n_observations = paste0(
      "3465 trough serum digoxin concentrations drawn at least 6 h after ",
      "the last dose and at least 5 days after the start of therapy, all ",
      "treated as steady-state troughs; median 5 [2-12] samples per patient ",
      "(Hirai 2022 Methods and Results)"
    ),
    age_range = "Adults >= 18 years; mean 67 +/- 14 years (Hirai 2022 Table 1)",
    weight_range = "Mean 57 +/- 15 kg (Hirai 2022 Table 1)",
    sex_female_pct = 38,
    race_ethnicity = c(Asian = 100),
    disease_state = paste0(
      "Atrial fibrillation (80% permanent/persistent, 20% paroxysmal) with ",
      "ACC/AHA stage C or D heart failure; NYHA class II/III/IV 291/67/33; ",
      "LVEF 39 +/- 14% (Hirai 2022 Table 1)."
    ),
    dose_range = paste0(
      "Oral digoxin maintenance 0.25 mg/day (13%), 0.125 mg/day (73%), ",
      "0.0625 mg/day (10%) and other regimens (4%); median treatment ",
      "duration 350 [60-1340] days (Hirai 2022 Table 1)"
    ),
    regions = "Japan (single centre; Tokyo Women's Medical University Hospital, 2008-2016)",
    renal_function = paste0(
      "Cockcroft-Gault CLcr median 56.5 [40.7-75.6] mL/min; eGFR 58.1 ",
      "[44.7-71.0] mL/min/1.73 m^2; serum creatinine 0.94 [0.75-1.15] mg/dL ",
      "(Hirai 2022 Table 1)"
    ),
    co_medication = paste0(
      "Amiodarone 63 (16%), diltiazem 33 (8%), verapamil 23 (6%) (Hirai 2022 ",
      "Table 1). Diltiazem was a candidate covariate but only amiodarone was ",
      "retained in the final model."
    ),
    notes = paste0(
      "Retrospective real-world therapeutic-drug-monitoring cohort. Assays: ",
      "COBAS TDM (LLOD 0.3 ng/mL) until November 2016, then Nanopia TDM ",
      "(LLOD 0.05 ng/mL); concentrations below the limit were excluded. ",
      "Fitted in Phoenix NLME 8.1; 1000-replicate bootstrap (100% success, ",
      "Table 3)."
    )
  )

  ini({
    # Hirai 2022 final model (Results):
    #   CL/F = 6.2 x (CLcr/60)^0.41 x (1 - 0.24 x [if amiodarone])
    #   Vd/F = 6.0 x Body weight ; ka = 1.0
    # Table 3 prints the same parameters to three decimals; those are used.
    lka <- fixed(log(1.0))
    label("Absorption rate constant (1/h)") # Hirai 2022 Table 3 'ka (fixed)' = 1.000 1/h; Methods: fixed because only troughs were sampled
    lcl <- log(6.209)
    label("Apparent oral clearance CL/F at CRCL = 60 mL/min without amiodarone (L/h)") # Hirai 2022 Table 3 'CL/F, L/h' = 6.209 (RSE 2.830%)
    lvc <- fixed(log(6.0))
    label("Apparent volume of distribution per kg body weight (L/kg)") # Hirai 2022 Table 3 'Vd/F, L/kg (fixed)' = 6.000; Methods: fixed from literature (ref. 24)

    e_crcl_cl <- 0.409
    label("Power exponent of (CRCL/60) on CL/F (unitless)") # Hirai 2022 Table 3 'CLCR on CL/F' = 0.409 (RSE 9.491%)
    e_amio_cl <- -0.238
    label("Fractional change in CL/F with concurrent amiodarone (unitless)") # Hirai 2022 Table 3 'Amiodarone on CL/F' = -0.238 (RSE 3.158%)

    # Exponential IIV on CL/F (Methods: CL/F = tv CL/F x exp(eta)). Table 3
    # prints omega(CL/F) = 34.4%; taken as the SD of eta, variance 0.344^2.
    # The CV reading log(1 + 0.344^2) = 0.112 differs by 5.6% in variance and
    # cannot be discriminated by Table 4 (see vignette).
    etalcl ~ 0.118336 # Hirai 2022 Table 3 'omega CL/F, %' = 34.4 (RSE 2.9%); 0.344^2

    # Multiplicative residual error (Methods: Cobs = Cpred x (1 + eps)).
    propSd <- 0.366
    label("Proportional residual error (fraction)") # Hirai 2022 Table 3 'Intraindividual variability, Multiplicative, %' = 36.6 (RSE 3.2%)
  })

  model({
    # CLcr above 120 mL/min replaced by 120 mL/min (Hirai 2022 Methods).
    crcl_capped <- min(CRCL, 120)

    ka <- exp(lka)
    cl <- exp(lcl + etalcl) *
      (crcl_capped / 60)^e_crcl_cl *
      (1 + e_amio_cl * CONMED_AMIO)
    vc <- exp(lvc) * WT

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in mg and vc in L give mg/L; x 1000 gives ng/mL (Table 4 units).
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
