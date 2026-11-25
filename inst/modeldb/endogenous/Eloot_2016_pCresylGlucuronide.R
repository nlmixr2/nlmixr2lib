Eloot_2016_pCresylGlucuronide <- function() {
  description <- paste(
    "Two-compartment hemodialysis kinetic model of the uremic solute p-cresylglucuronide (PCG; about 10% protein bound)",
    "in adults on chronic high-flux hemodialysis (Eloot 2016). Plasmatic (V1) and extraplasmatic",
    "(V2) compartments exchange by an intercompartment clearance K21; the solute is generated at a",
    "constant rate G into the plasmatic compartment and removed by the dialyzer clearance K while a",
    "session runs. Net ultrafiltration shrinks both volumes in proportion during dialysis and the",
    "fluid is regained between sessions. Typical values only: the paper reports the between-patient",
    "summary of individual fits, not a population model.",
    sep = " "
  )
  reference <- paste(
    "Eloot S, Schneditz D, Cornelis T, Van Biesen W, Glorieux G, Dhondt A, Kooman J, Vanholder R.",
    "Protein-Bound Uremic Toxin Profiling as a Tool to Optimize Hemodialysis.",
    "PLoS One. 2016;11(1):e0147159. doi:10.1371/journal.pone.0147159.",
    sep = " "
  )
  vignette <- "Eloot_2016_uremicToxins"
  units <- list(time = "min", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    RRT_HEMODIAL_ACTIVE = list(
      description = "Hemodialysis-active indicator (1 while a dialysis session is running, 0 between sessions)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Time-varying. Switches on the dialyzer clearance K and the net ultrafiltration rate, and",
        "switches off the interdialytic fluid regain. The paper's reference schedule is three",
        "240-min sessions per week."
      ),
      source_name = "(dialysis on/off; not a named column in the paper)"
    )
  )

  compartmentData <- list(
    central = list(analyte = "p-cresylglucuronide", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "p-cresylglucuronide", units = "mg", specimen = "tissue", verified = TRUE),
    dialysate = list(analyte = "p-cresylglucuronide", units = "mg", specimen = "dialysate", verified = TRUE),
    fluid_balance = list(analyte = "body water", units = "L", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 10L,
    n_studies = 1L,
    age_range = "69 +/- 12 years (mean +/- SD)",
    weight_range = "75.0 [69.9;82.6] kg (median [IQR])",
    sex_female_pct = 20,
    disease_state = "End-stage kidney disease on chronic thrice-weekly 4-h high-flux hemodialysis",
    dose_range = "No drug administered; endogenous solute generation. One midweek 240-min session at blood flow 300 mL/min and dialysate flow 700 mL/min.",
    regions = "Belgium (Ghent University Hospital)",
    renal_function = "Residual renal function 2.6 [0.4;3.8] mL/min (median [IQR]); not included in the model",
    notes = paste(
      "Ten stable patients (2 women, 8 men), 50 [30;62] months on dialysis, hematocrit",
      "37.2 +/- 4.2%, total plasma protein 61.4 +/- 6.4 g/L, ultrafiltration 0.41 +/- 0.29 L/h,",
      "single-pool Kt/Vurea 1.7 +/- 0.3. Dialyzers: FX800 (n = 6), Evodial, Xenium 210,",
      "Phylter PHF 17G and Phylter HF 17SD (n = 1 each). Nine fistulas, one central venous",
      "catheter. Methods 'Patients'."
    )
  )

  ini({
    lksyn <- log(0.019); label("Solute generation rate G (mg/min)") # Table 2 PCG 'G (mg/min)' = 0.019 [0.016;0.049] (median [IQR])
    lcl_hemodialysis <- log(152 / 1000); label("Dialyzer clearance K while a session runs (L/min)") # Table 2 PCG 'K (mL/min)' = 152 +/- 28 (mean +/- SD), converted to L/min
    lvc <- log(8.2); label("Plasmatic volume V1 at the start of dialysis (L)") # Table 3 PCG 'V1 (L)' = 8.2 +/- 2.4 (mean +/- SD)
    lvp <- log(18.9 - 8.2); label("Extraplasmatic volume V2 = Vtot - V1 at the start of dialysis (L)") # Table 3 PCG 'Vtot (L)' = 18.9 [15.7;25.2] (median) minus V1 = 8.2 (mean)
    lq <- log(103 / 1000); label("Intercompartment clearance K21 (L/min)") # Table 3 PCG 'K21 (mL/min)' = 103 +/- 47 (mean +/- SD), converted to L/min
    lrbase <- log(0.35 * 10); label("Predialysis concentration Cpre, both compartments (mg/L)") # Table 2 PCG 'Cpre (mg/dL)' = 0.35 [0.26;0.79] (median [IQR]), converted to mg/L
    qpfr <- fixed(0.41 / 60); label("Net ultrafiltration rate QUF during dialysis (L/min)") # Methods 'Patients': ultrafiltration rates 'averaged 0.41 +/- 0.29 L/h', converted to L/min; held at the cohort mean
    qgain <- fixed(0.41 / 60 * 720 / (10080 - 720)); label("Interdialytic fluid regain rate (L/min)") # Not printed: weekly ultrafiltrate (0.41 L/h x 12 h) spread over the 156 interdialytic hours of the 3 x 4 h schedule, so the weekly fluid balance is zero (Methods 'Kinetic analysis': same weekly weight loss in every strategy)
  })

  model({
    ksyn <- exp(lksyn)
    cl_hemodialysis <- exp(lcl_hemodialysis)
    q <- exp(lq)
    rbase <- exp(lrbase)

    # Ultrafiltration is distributed over both compartments in proportion to
    # their volumes (Methods 'Calculations and kinetic modeling'), so V1/Vtot is
    # constant and both volumes scale with the current total body volume.
    vtot0 <- exp(lvc) + exp(lvp)
    vscale <- (vtot0 + fluid_balance) / vtot0
    vc <- exp(lvc) * vscale
    vp <- exp(lvp) * vscale

    # Anuric patients: the dialyzer is the only elimination route.
    cl <- RRT_HEMODIAL_ACTIVE * cl_hemodialysis

    Cc <- central / vc
    Cp <- peripheral1 / vp

    # Eq. 1-2 in amount form: d(V1*C1)/dt = G - K*C1 + K21*(C2 - C1) and
    # d(V2*C2)/dt = -K21*(C2 - C1).
    d/dt(central) <- ksyn - cl * Cc + q * (Cp - Cc)
    d/dt(peripheral1) <- q * (Cc - Cp)
    # Cumulative amount collected in spent dialysate (total solute removal).
    d/dt(dialysate) <- cl * Cc
    # Net body-fluid change from the start of the calibrated session (L).
    d/dt(fluid_balance) <- (1 - RRT_HEMODIAL_ACTIVE) * qgain - RRT_HEMODIAL_ACTIVE * qpfr

    # Both compartments start at the measured predialysis plasma concentration
    # (Methods 'Calculations and kinetic modeling').
    central(0) <- rbase * exp(lvc)
    peripheral1(0) <- rbase * exp(lvp)
  })
}
