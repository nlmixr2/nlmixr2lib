Thorsted_2023_neutrophil_invitro <- function() {
  description <- paste(
    "In vitro (human neutrophils with Staphylococcus epidermidis, S. aureus,",
    "Escherichia coli, Streptococcus pyogenes and Corynebacterium in",
    "suspension or fibrin gel). Semi-mechanistic model of neutrophil-mediated",
    "phagocytosis and digestion of bacteria (Thorsted 2023), in vitro arm.",
    "Bacteria follow the Nielsen 2007 growing / resting system with an added",
    "growth-lag state: all bacteria start in bact_lag, move to the growing",
    "state bact_susceptible (klag), replicate (kgrowth) and move to the",
    "resting state bact_resting at kSR = (L + S + R) * (kgrowth - kdeath) /",
    "Bmax, with natural death (kdeath) from every state. A static neutrophil",
    "concentration (neut_tube) phagocytoses bacteria from all three states",
    "into bact_phagocytosed and digests them from there. Both first-order",
    "rates share a maximum (kN,max) and a time-related loss (kN,loss), are",
    "Emax functions of the neutrophil concentration with process-specific",
    "N50, and fall as the phagocytosed-bacteria-to-neutrophil ratio rises",
    "(shared P/N50). Bacteria and neutrophils are both dosed at time 0 as",
    "concentrations: bacteria in CFU/mL into bact_lag, neutrophils in",
    "cells/mL into neut_tube. Three log10 CFU/mL outputs share one additive",
    "residual error: total viable bacteria (log_cfu, extracellular plus",
    "phagocytosed-not-digested), extracellular bacteria",
    "(log_cfu_extracellular) and intracellular bacteria",
    "(log_cfu_intracellular). The in vivo mouse arm of the same simultaneous",
    "fit is Thorsted_2023_neutrophil_mouse."
  )
  reference <- paste(
    "Thorsted A, Pham AD, Friberg LE, Nielsen EI. Model-based assessment of",
    "neutrophil-mediated phagocytosis and digestion of bacteria across in",
    "vitro and in vivo studies. CPT Pharmacometrics Syst Pharmacol.",
    "2023;12(12):1972-1987. doi:10.1002/psp4.13046. PMID: 37700716.",
    "PMCID: PMC10725272. Parameter estimates: Table 2, 'Final model' column.",
    "Equations: Materials and Methods, Equations 1-6. Model code: Data S2",
    "(NONMEM control stream); digitized data: Table S1.",
    sep = " "
  )
  vignette <- "Thorsted_2023_neutrophil_bacterial_killing"
  units <- list(
    time = "h",
    dosing = "CFU/mL (bacterial inoculum); cells/mL (neutrophils)",
    concentration = "log10 CFU/mL (observations); CFU/mL and cells/mL (states)"
  )

  covariateData <- list()

  # The states hold concentrations (the in vitro convention in which the
  # compartment state IS the tube concentration), so dosing amounts are
  # concentrations. The bacterial life-cycle states split one phenotype along
  # a growth axis (lag / growing / resting / phagocytosed), not along the
  # drug-resistance axis of the registered ^bact_(susceptible|...) regex, so
  # they are declared paper-specific as in Nielsen_2007_semimechanistic_antibiotic_pd
  # and Saporta_2026_meropenem.
  paper_specific_compartments <- c(
    "bact_lag",
    "bact_susceptible",
    "bact_resting",
    "bact_phagocytosed",
    "neut_tube",
    "log_cfu_extracellular",
    "log_cfu_intracellular"
  )

  compartmentData <- list(
    neut_tube = list(
      analyte = "neutrophil granulocytes (static in vitro concentration)",
      units = "cells/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    bact_lag = list(
      analyte = "bacteria, growth-lag (L) state",
      units = "CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    bact_susceptible = list(
      analyte = "bacteria, growing (S) state",
      units = "CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    bact_resting = list(
      analyte = "bacteria, resting (R) state",
      units = "CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    bact_phagocytosed = list(
      analyte = "bacteria phagocytosed by neutrophils and not yet digested (PN)",
      units = "CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "in vitro (human neutrophils; S. epidermidis, S. aureus, E. coli, S. pyogenes, Corynebacterium)",
    n_subjects = 115L,
    n_studies = 6L,
    organism = "Staphylococcus epidermidis, Staphylococcus aureus, Escherichia coli, Streptococcus pyogenes, Corynebacterium (pooled; the data did not support species-specific killing parameters)",
    system = "Neutrophil-bacteria co-incubation in suspension (5 studies) or fibrin gel (1 study) with at least 5% human serum for opsonization at 37 C; durations 1-3 h",
    starting_inoculum = "10^3 to 10^10 CFU/mL",
    neutrophil_range = "0 (growth controls) to 1 x 10^7 cells/mL, static within an experiment",
    disease_state = "not applicable (in vitro bactericidal assays)",
    notes = paste(
      "Seven published in vitro studies were digitized (Table 1, In vitro-1 to In vitro-7); In vitro-7 (Malka et al., J Clin Invest; reference 16) was reserved for external evaluation (Figure 3) and not used for estimation.",
      "115 in vitro experiments with 403 CFU observations (Results). n_subjects counts experiments, the unit that NONMEM treated as an individual.",
      "Measurement types: total viable CFU after neutrophil lysis (In vitro-1, -2, -5, -7), extracellular CFU after centrifugal separation (In vitro-3) and intracellular CFU after removal of extracellular bacteria (In vitro-4, -6) (Materials and Methods).",
      "No between-experiment variability was estimated ($OMEGA 0 FIX in Data S2); the model carries typical values plus residual error only."
    )
  )

  ini({
    # ---- Bacterial dynamics, in vitro (Table 2, Final model) ---------------
    lkgrowth <- log(1.84); label("Bacterial first-order growth rate constant, in vitro (1/h)") # Table 2: kgrowth (in vitro) = 1.84 1/h, RSE 15%
    lkdeath <- fixed(log(0.179)); label("Bacterial first-order natural death rate constant (1/h)") # Table 2: kdeath = 0.179 1/h, Fixed (from Nielsen 2007)
    lklag <- log(0.290); label("Bacterial first-order growth-lag rate constant, L to S (1/h)") # Table 2: klag = 0.290 1/h, RSE 40%
    lbmax <- fixed(log(2.5e9)); label("Maximum system bacterial capacity (CFU/mL)") # Table 2: Bmax = 2.5 x 10^9 CFU/mL, Fixed (highest digitized CFU)

    # ---- Neutrophil-mediated killing, in vitro (Table 2, Final model) ------
    # kN,max, kN,loss and P/N50 are each shared between the digestion and the
    # phagocytosis process (Table 2 footnote b; Results: 'Six parameters could
    # be reduced to three'); N50 stays process-specific.
    lkn_max <- log(8.85); label("Neutrophil first-order maximum phagocytosis and digestion rate constant, in vitro (1/h)") # Table 2: kN,max,dig = kN,max,phag = 8.85 1/h, RSE 13%
    lkn_loss <- log(0.779); label("Time-related decrease in neutrophil phagocytosis and digestion, in vitro (1/h)") # Table 2: kN,loss,dig = kN,loss,phag = 0.779 1/h, RSE 13%
    ln50_dig <- log(6.55e6); label("Neutrophil concentration giving 50% of the maximum digestion rate (cells/mL)") # Table 2: N50,dig = 6.55 x 10^6 N/mL, RSE 1.4%
    ln50_phag <- log(1.19e6); label("Neutrophil concentration giving 50% of the maximum phagocytosis rate (cells/mL)") # Table 2: N50,phag = 1.19 x 10^6 N/mL, RSE 4.3%
    lpn50 <- log(3.41); label("Phagocytosed-bacteria-to-neutrophil ratio halving the phagocytosis and digestion rates (CFU/cell)") # Table 2: P/N50,dig = P/N50,phag = 3.41 CFU/N, RSE 22%

    # ---- Residual error ----------------------------------------------------
    # Log10 transform-both-sides with an additive error (Model development);
    # Table 2 reports ERR on the standard-deviation scale. One SIGMA serves
    # every in vitro CFU measurement type (Data S2: W = SQRT(SIGMA(1,1)) for
    # LOCA 1, 2 and 3). nlmixr2 needs one residual parameter per endpoint, so
    # the single estimate is carried three times with the same value.
    addSd_log_cfu <- 0.279; label("Additive residual SD on log10 total viable in vitro CFU (log10 CFU/mL)") # Table 2: ERR (in vitro CFU) = 0.279, RSE 8.0%
    addSd_log_cfu_extracellular <- 0.279; label("Additive residual SD on log10 extracellular in vitro CFU (log10 CFU/mL; same estimate as addSd_log_cfu)") # Table 2: ERR (in vitro CFU) = 0.279, single shared SIGMA
    addSd_log_cfu_intracellular <- 0.279; label("Additive residual SD on log10 intracellular in vitro CFU (log10 CFU/mL; same estimate as addSd_log_cfu)") # Table 2: ERR (in vitro CFU) = 0.279, single shared SIGMA
  })

  model({
    kgrowth <- exp(lkgrowth)
    kdeath <- exp(lkdeath)
    klag <- exp(lklag)
    bmax <- exp(lbmax)
    kn_max <- exp(lkn_max)
    kn_loss <- exp(lkn_loss)
    n50_dig <- exp(ln50_dig)
    n50_phag <- exp(ln50_phag)
    pn50 <- exp(lpn50)

    # Phagocytosed bacteria per neutrophil. Growth controls carry no
    # neutrophils; NONMEM's $ABBR PROTECT keeps the 0/0 finite there, and the
    # ratio is irrelevant because both rates carry the factor N / (N50 + N) = 0.
    pn_ratio <- 0
    if (neut_tube > 0) pn_ratio <- bact_phagocytosed / neut_tube

    # Equations 4 and 5. The exponential loss runs on model time, which is
    # the time since neutrophils and bacteria were combined (Data S2:
    # EXP(-KNA * (T - EXTRT)) with EXTRT = 0 for every in vitro study).
    kn_common <- kn_max * exp(-kn_loss * t) * (1 - pn_ratio / (pn_ratio + pn50))
    kn_dig <- kn_common * neut_tube / (neut_tube + n50_dig)
    kn_phag <- kn_common * neut_tube / (neut_tube + n50_phag)

    # kSR is driven by the extracellular burden L + S + R only (Data S2:
    # CFU = A(2) + A(3) + A(4)).
    bact_ext <- bact_lag + bact_susceptible + bact_resting
    ksr <- bact_ext * (kgrowth - kdeath) / bmax

    # Equations 1-3 extended with phagocytosis from every state, and Equation 6.
    d/dt(neut_tube) <- 0
    d/dt(bact_lag) <- -kdeath * bact_lag - klag * bact_lag - kn_phag * bact_lag
    d/dt(bact_susceptible) <- -kdeath * bact_susceptible + klag * bact_lag - ksr * bact_susceptible + kgrowth * bact_susceptible - kn_phag * bact_susceptible
    d/dt(bact_resting) <- -kdeath * bact_resting + ksr * bact_susceptible - kn_phag * bact_resting
    d/dt(bact_phagocytosed) <- kn_phag * bact_ext - kn_dig * bact_phagocytosed

    # Observations (Data S2 $ERROR): non-positive totals are floored at 1 CFU/mL
    # before the log10 so the prediction stays finite.
    cfu_ext <- bact_ext
    if (cfu_ext <= 0) cfu_ext <- 1
    cfu_total <- bact_ext + bact_phagocytosed
    if (cfu_total <= 0) cfu_total <- 1
    cfu_int <- bact_phagocytosed
    if (cfu_int <= 0) cfu_int <- 1

    log_cfu <- log10(cfu_total)
    log_cfu_extracellular <- log10(cfu_ext)
    log_cfu_intracellular <- log10(cfu_int)

    log_cfu ~ add(addSd_log_cfu)
    log_cfu_extracellular ~ add(addSd_log_cfu_extracellular)
    log_cfu_intracellular ~ add(addSd_log_cfu_intracellular)
  })
}
