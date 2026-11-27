Ahlberg_2019_spesolimab <- function() {
  description <- paste(
    "Translational two-compartment semi-mechanistic target-mediated drug disposition (TMDD)",
    "model with the rapid-binding approximation for spesolimab (MAB92, BI 655130; humanized",
    "anti-IL-36R IgG1) after intravenous dosing in healthy men, predicted before the",
    "first-in-human study (Ahlberg 2019). Linear PK from Dedrick-scaled cynomolgus monkey",
    "data, KD from surface plasmon resonance, kint = kdeg from the average plasma-membrane",
    "turnover, and baseline IL-36R (R0) = the mouse estimate times 0.013, the human-to-mouse",
    "FANTOM5 IL-36R transcript ratio. Total target is a turnover state. Deterministic",
    "prediction model (no IIV, no residual error)."
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
    central = list(analyte = "spesolimab (free)", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "spesolimab", units = "nmol", specimen = "tissue", verified = TRUE),
    total_target = list(
      analyte = "IL-36R (total, free + spesolimab-bound)",
      units = "nmol/L",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 1,
    sex_female_pct = 0,
    disease_state = "Healthy male volunteers",
    dose_range = "Single intravenous dose 0.03, 0.05, 0.1, 0.3, 1, 3, 6 or 10 mg/kg",
    study_id = "NCT02525679 (Phase 1 single rising dose)",
    notes = paste(
      "No parameter was estimated from human data: the model is the pre-clinical prediction,",
      "checked against first-in-human NCA (Table 5). The paper does not report the number,",
      "age or body weight of the volunteers (Materials and methods, Clinical pharmacokinetic",
      "study)."
    )
  )

  ini({
    # Table 4 'Translated human parameter' column. Rate constants are printed in 'day-1'
    # but run in hours (kdeg = kint = 1.39 is ln(2) / 0.5 h; with 1/h units the model
    # reproduces Table 5 Cmax and AUC). Values are used as printed in 1/h.
    lrbase <- log(0.010); label("Baseline total IL-36R concentration R0 (nmol/L)") # Table 4 translated human R0 = 0.010 nM = mouse 0.794 x 0.013 (FANTOM5 human/mouse TPM ratio, text)
    lkd <- fixed(log(0.020)); label("Equilibrium dissociation constant KD, spesolimab to human IL-36R (nmol/L)") # Table 4 translated human KD = 0.020 nM (in vitro SPR)
    lkdeg <- fixed(log(1.39)); label("Free target degradation rate constant (1/h)") # Table 4 translated human kdeg = 1.39
    lkint <- fixed(log(1.39)); label("Drug-target complex internalisation rate constant (1/h)") # Table 4 translated human kint = 1.39
    lvc <- log(3.0); label("Central volume of distribution (L)") # Table 4 translated human Vc = 3.0 L
    lkel <- log(0.0021); label("Linear elimination rate constant (1/h)") # Table 4 translated human kel = 0.0021
    lk12 <- log(0.0089); label("Central-to-peripheral rate constant (1/h)") # Table 4 translated human k12 = 0.0089
    lk21 <- log(0.0110); label("Peripheral-to-central rate constant (1/h)") # Table 4 translated human k21 = 0.0110
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

    # Target synthesis keeps the drug-free target at its baseline: kin = kdeg * R0.
    kin <- kdeg * rbase

    # Free spesolimab concentration in central (nmol/L).
    Cc <- central / vc

    # Rapid-binding TMDD equations as printed in Results, without the absorption term
    # (human doses are intravenous). The printed A1 equation shows '- A2*k21'; the
    # redistribution term is entered as + k21 * peripheral1, the sign the A2 equation
    # (A1*K12 - A2*K21) requires for mass balance.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 -
      total_target * kint / (kd + Cc) * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(total_target) <- kin - kdeg * total_target - (kint - kdeg) * total_target / (kd + Cc) * Cc
    total_target(0) <- rbase
  })
}
