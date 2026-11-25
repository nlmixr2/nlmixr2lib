Ahlberg_2019_MAB04_mouse_mm <- function() {
  description <- paste(
    "Preclinical (mouse). Two-compartment PK model with parallel linear and Michaelis-Menten",
    "elimination from the central compartment for MAB04, a rat/mouse chimeric IgG2a surrogate",
    "antibody against mouse IL-36R, after intraperitoneal dosing in female C57BL/6 mice",
    "(Ahlberg 2019). First-order absorption from the peritoneal depot with bioavailability",
    "assumed 100%. Fitted to mean data; deterministic (no IIV, no residual error reported).",
    "Its Vc and k10 feed the mouse TMDD model, Ahlberg_2019_MAB04_mouse."
  )
  reference <- paste(
    "Ahlberg J, Giragossian C, Li H, Myzithras M, Raymond E, Caviness G, Grimaldi C,",
    "Brown SE, Perez R, Yang D, Kroe-Barrett R, Joseph D, Pamulapati C, Coble K, Ruus P,",
    "Woska JR, Ganesan R, Hansel S, Mbow ML.",
    "Retrospective analysis of model-based predictivity of human pharmacokinetics for",
    "anti-IL-36R monoclonal antibody MAB92 using a rat anti-mouse IL-36R monoclonal antibody",
    "and RNA expression data (FANTOM5). MAbs. 2019;11(5):956-964.",
    "doi:10.1080/19420862.2019.1615345"
  )
  vignette <- "Ahlberg_2019_spesolimab"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "MAB04", units = "nmol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "MAB04", units = "nmol", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "MAB04", units = "nmol", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "mouse (C57BL/6, female)",
    n_subjects = 9,
    n_studies = 1,
    sex_female_pct = 100,
    disease_state = "Healthy (in-house PK study)",
    dose_range = "Single intraperitoneal dose of MAB04 0.3, 1.5 or 10 mg/kg (N = 3 per dose group)",
    notes = paste(
      "Materials and methods, Mouse studies. Whole-blood microsampling (10 uL) over one week",
      "(0.3 mg/kg) or two weeks (1.5 and 10 mg/kg); serum-equivalent concentrations derived",
      "from whole blood with a hematocrit factor of 0.45. The model was fitted to the mean",
      "concentration-time profile of each dose group simultaneously. Body weight not reported."
    )
  )

  ini({
    # Table 3 prints every rate constant in 'day-1', but the model was run in hours:
    # Vmax is in nmol/h, Figure 4 time axis is in h, and Vc * k10 in 1/h reproduces the
    # Table 2 NCA clearance at 10 mg/kg (0.13 mL/h/kg). The values are used as printed in 1/h.
    lvc <- log(0.000403); label("Central volume of distribution (L)") # Table 3 Vc = 0.000403 L (SE 0.0000308)
    lkel <- log(0.00506); label("Linear elimination rate constant k10 (1/h)") # Table 3 k10 = 0.00506 (SE 0.00147)
    lvmax <- log(0.000533); label("Maximum rate of nonlinear elimination (nmol/h)") # Table 3 Vmax = 0.000533 nmol/h (SE 0.0000626)
    lkm <- log(1.18); label("Michaelis-Menten constant (nmol/L)") # Table 3 Km = 1.18 nM (SE 0.729)
    lk12 <- fixed(log(0.086)); label("Central-to-peripheral rate constant (1/h)") # Table 3 k12 = 0.086 Fixed (in-house typical mouse mAb value)
    lk21 <- fixed(log(0.063)); label("Peripheral-to-central rate constant (1/h)") # Table 3 k21 = 0.063 Fixed (in-house typical mouse mAb value)
    lka <- fixed(log(0.152)); label("Intraperitoneal absorption rate constant (1/h)") # Table 3 ka = 0.152 Fixed (in-house typical mouse mAb value)
    lfdepot <- fixed(log(1)); label("Intraperitoneal bioavailability (fraction)") # Modeling approach text: bioavailability assumed 100%
  })

  model({
    vc <- exp(lvc)
    kel <- exp(lkel)
    vmax <- exp(lvmax)
    km <- exp(lkm)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    ka <- exp(lka)
    fdepot <- exp(lfdepot)

    Cc <- central / vc

    # Figure 2: parallel linear (k10) and Michaelis-Menten (Vmax, Km) elimination from central.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1 -
      vmax * Cc / (km + Cc)
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    f(depot) <- fdepot
  })
}
