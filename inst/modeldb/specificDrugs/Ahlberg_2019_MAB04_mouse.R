Ahlberg_2019_MAB04_mouse <- function() {
  description <- paste(
    "Preclinical (mouse). Two-compartment semi-mechanistic target-mediated drug disposition",
    "(TMDD) model with the rapid-binding approximation for MAB04, a rat/mouse chimeric IgG2a",
    "surrogate antibody against mouse IL-36R, after intraperitoneal dosing in female C57BL/6",
    "mice (Ahlberg 2019). Linear PK from the parallel-elimination fit",
    "(Ahlberg_2019_MAB04_mouse_mm), KD from surface plasmon resonance, kint = kdeg from the",
    "average plasma-membrane turnover, and baseline IL-36R (R0) estimated. Total target is a",
    "turnover state. Fitted to mean data; deterministic (no IIV, no residual error reported)."
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
    central = list(analyte = "MAB04 (free)", units = "nmol", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "MAB04", units = "nmol", specimen = "tissue", verified = TRUE),
    total_target = list(
      analyte = "IL-36R (total, free + MAB04-bound)",
      units = "nmol/L",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "mouse (C57BL/6, female)",
    n_subjects = 9,
    n_studies = 1,
    sex_female_pct = 100,
    disease_state = "Healthy (in-house PK study)",
    dose_range = "Single intraperitoneal dose of MAB04 0.3, 1.5 or 10 mg/kg (N = 3 per dose group)",
    notes = paste(
      "Materials and methods, Mouse studies. Free MAB04 measured by IL-36R-capture ELISA in",
      "whole-blood microsamples, converted to serum-equivalent concentrations (hematocrit",
      "0.45). Fitted to mean profiles (Figure 4). Body weight not reported."
    )
  )

  ini({
    # Table 4 prints the rate constants in 'day-1'; the model runs in hours (Figure 4
    # time axis in h; kdeg = kint = 1.39 is ln(2) / 0.5 h, the plasma-membrane
    # turnover the text cites). Values are used as printed in 1/h.
    lrbase <- log(0.794); label("Baseline total IL-36R concentration R0 (nmol/L)") # Table 4 mouse R0 = 0.794 nM (estimated by model fitting)
    lkd <- fixed(log(0.24)); label("Equilibrium dissociation constant KD, MAB04 to mouse IL-36R (nmol/L)") # Table 4 mouse KD = 0.24 nM (in vitro SPR)
    lkdeg <- fixed(log(1.39)); label("Free target degradation rate constant (1/h)") # Table 4 mouse kdeg = 1.39; text 'average plasma membrane turnover rate of 0.5 h'
    lkint <- fixed(log(1.39)); label("Drug-target complex internalisation rate constant (1/h)") # Table 4 mouse kint = 1.39; text 'applied for both kdeg and kint'
    lvc <- fixed(log(0.000403)); label("Central volume of distribution (L)") # Table 4 mouse Vc = 0.000403 L (from the Table 3 parallel-elimination fit)
    lkel <- fixed(log(0.00506)); label("Linear elimination rate constant (1/h)") # Table 4 mouse kel = 0.00506 (Table 3 k10)
    lk12 <- fixed(log(0.086)); label("Central-to-peripheral rate constant (1/h)") # Table 4 mouse k12 = 0.086
    lk21 <- fixed(log(0.063)); label("Peripheral-to-central rate constant (1/h)") # Table 4 mouse k21 = 0.063
    lka <- fixed(log(0.152)); label("Intraperitoneal absorption rate constant (1/h)") # Table 4 mouse ka = 0.152
    lfdepot <- fixed(log(1)); label("Intraperitoneal bioavailability (fraction)") # Modeling approach text: bioavailability assumed 100%
  })

  model({
    rbase <- exp(lrbase)
    kd <- exp(lkd)
    kdeg <- exp(lkdeg)
    kint <- exp(lkint)
    vc <- exp(lvc)
    kel <- exp(lkel)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    ka <- exp(lka)
    fdepot <- exp(lfdepot)

    # Target synthesis keeps the drug-free target at its baseline: kin = kdeg * R0.
    kin <- kdeg * rbase

    # Free MAB04 concentration in central (nmol/L).
    Cc <- central / vc

    # Rapid-binding TMDD equations as printed in Results (A1 = central, A2 =
    # peripheral1, A0 = depot, Rtot = total_target). The printed A1 equation shows
    # '- A2*k21'; the redistribution term is entered as + k21 * peripheral1, the sign
    # the A2 equation (A1*K12 - A2*K21) requires for mass balance.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1 -
      total_target * kint / (kd + Cc) * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(total_target) <- kin - kdeg * total_target - (kint - kdeg) * total_target / (kd + Cc) * Cc
    total_target(0) <- rbase
    f(depot) <- fdepot
  })
}
