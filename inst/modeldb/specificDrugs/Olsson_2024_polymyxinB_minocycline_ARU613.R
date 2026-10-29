Olsson_2024_polymyxinB_minocycline_ARU613 <- function() {
  description <- "In vitro (Klebsiella pneumoniae ARU613, OXA-48-producing clinical isolate; polymyxin B MIC 8 mg/L, minocycline MIC 8 mg/L). Semi-mechanistic PK/PD model of polymyxin B (PMB) and minocycline (MIN) alone and in combination, used for in silico predictions of 72 h dynamic time-kill experiments. Growing drug-susceptible bacteria (bact_susceptible) transfer to drug-insusceptible resting bacteria (bact_resting) as the total count approaches Bmax. Each drug kills the susceptible bacteria through a power function of its concentration, attenuated by a drug-specific irreversible adaptive-resistance system (aroff -> aron at a concentration-driven rate). PMB enhances MIN in two ways: it raises the MIN kill rate through an Emax function of the PMB concentration and it raises the rate at which MIN-induced adaptive resistance develops through a power function of the PMB concentration. Drug concentrations are user-dosed bath states (cpmb, cmino; mg/L). Companion model: Olsson_2024_polymyxinB_minocycline_ARU705 (KPC-3-producing isolate)."
  reference <- paste(
    "Olsson A, Malmberg C, Zhao C, Friberg LE, Nielsen EI, Lagerback P,",
    "Tangden T. (2024). Synergy of polymyxin B and minocycline against KPC-3-",
    "and OXA-48-producing Klebsiella pneumoniae in dynamic time-kill",
    "experiments: agreement with in silico predictions. J Antimicrob Chemother",
    "79(2):391-402.",
    "doi:10.1093/jac/dkad394. PMCID: PMC10832586.",
    "Model equations: Supplementary Table S1 (Eqs 1-10). Parameter values:",
    "Table 3, ARU613 columns, which reproduce the estimates of the original",
    "static time-kill model for this strain (Zhao C et al. Int J Antimicrob",
    "Agents 2020;55:105941, doi:10.1016/j.ijantimicag.2020.105941).",
    sep = " "
  )
  vignette <- "Olsson_2024_polymyxinB_minocycline"
  units <- list(time = "h", dosing = "mg/L", concentration = "mg/L")

  # The "dose" of a time-kill experiment is the medium concentration placed
  # directly into each drug state, so the dosing targets are the bath states.
  dosing <- c("cpmb", "cmino")

  # `bact_resting` (Saporta 2026 precedent), the PMB / MIN bath states and the
  # drug-suffixed adaptive-resistance pairs have no registered canonical:
  # neither `pmb` nor `mino` is a registered drug-suffix token.
  paper_specific_compartments <- c(
    "bact_resting",
    "aroff_pmb",
    "aron_pmb",
    "aroff_mino",
    "aron_mino",
    "cpmb",
    "cmino"
  )

  compartmentData <- list(
    bact_susceptible = list(
      analyte = "Klebsiella pneumoniae ARU613, growing drug-susceptible bacteria (S)",
      units = "CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    bact_resting = list(
      analyte = "Klebsiella pneumoniae ARU613, resting drug-insusceptible bacteria (R)",
      units = "CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    aroff_pmb = list(
      analyte = "polymyxin B adaptive resistance, off state (fraction)",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    ),
    aron_pmb = list(
      analyte = "polymyxin B adaptive resistance, on state (fraction)",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    ),
    aroff_mino = list(
      analyte = "minocycline adaptive resistance, off state (fraction)",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    ),
    aron_mino = list(
      analyte = "minocycline adaptive resistance, on state (fraction)",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    ),
    cpmb = list(analyte = "polymyxin B", units = "mg/L", specimen = "administration site", verified = TRUE),
    cmino = list(analyte = "minocycline", units = "mg/L", specimen = "administration site", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "in vitro (Klebsiella pneumoniae ARU613)",
    n_subjects = 1L,
    n_studies = 1L,
    organism = "Klebsiella pneumoniae ARU613, an OXA-48-producing clinical isolate provided by the Public Health Agency of Sweden. Meropenem MIC >= 16 mg/L (R), minocycline 8 mg/L (I; carries tet(D)), polymyxin B 8 mg/L (R; pmrB T246A, mgrB C39S). Table 1.",
    model_system = "Parameters come from 24 h static time-kill experiments in cation-adjusted Mueller-Hinton broth at 37 C (starting inoculum 6 log10 CFU/mL). Olsson 2024 used them to predict 72 h dynamic time-kill experiments in a 100 mL in vitro PK/PD model mimicking free plasma concentrations in critically ill patients.",
    dose_range = "Dynamic predictions: PMB 2.5 mg/kg loading dose then 1.5 mg/kg q12h (targeted Cmax 5.0 and 4.5 mg/L, Cmin 1.1 mg/L) and MIN 400 mg loading dose then 200 mg q12h (targeted Cmax 3.8 and 4.1 mg/L, Cmin 1.8 mg/L), alone and in combination.",
    disease_state = "not applicable (in vitro; carbapenem-resistant K. pneumoniae)",
    notes = paste(
      "Model fit and simulations were performed in NONMEM 7.4 and mrgsolve.",
      "Predictions started from an inoculum of 6 log10 CFU/mL and assumed that no bacteria could regrow once the total was predicted below 1 CFU in the 100 mL working volume (0.01 CFU/mL).",
      "Olsson 2024 reports typical values only; no between-experiment variability or residual error is given."
    )
  )

  ini({
    # --- Bacteria-related parameters (Table 3; shared by both strains) -------
    kgrowth <- 1.37; label("Rate constant of bacterial growth (1/h)") # Table 3: kgrowth = 1.37 h-1
    kdeath <- 0.179; label("Rate constant of natural bacterial death (1/h)") # Table 3: kdeath = 0.179 h-1
    log10bmax <- 9.53; label("Maximum bacterial concentration in the system, Bmax (log10 CFU/mL)") # Table 3: Bmax = 10^9.53 CFU/mL
    ltlag <- log(0.304); label("Log lag time before bacteria transfer from S to R (h)") # Table 3: Tlag = 0.304 h
    log10cfu0 <- fixed(6); label("Initial inoculum (log10 CFU/mL)") # Methods, PD predictions: 'inoculum sizes in the predictions were set to 6 log10 cfu/mL'

    # --- Minocycline single-drug parameters (Table 3, MIN ARU613) -----------
    slope1_mino <- 0.339; label("Slope of the power function for the MIN kill rate, Slope1 (L/mg/h)") # Table 3: Slope1 MIN ARU613 = 0.339
    slope2_mino <- 0.000179; label("Slope of the linear function for the MIN adaptive-resistance on-rate, Slope2 (L/mg/h)") # Table 3: Slope2 MIN ARU613 = 0.000179
    gamma1_mino <- 0.546; label("Exponent of the power function for the MIN kill rate, gamma1 (unitless)") # Table 3: gamma1 MIN ARU613 = 0.546
    slope3_mino <- -10.2; label("Slope of the linear adaptive-resistance function for MIN, Slope3 (unitless)") # Table 3: Slope3 MIN ARU613 = -10.2

    # --- Polymyxin B single-drug parameters (Table 3, PMB ARU613) -----------
    slope1_pmb <- 0.0690; label("Slope of the power function for the PMB kill rate, Slope1 (L/mg/h)") # Table 3: Slope1 PMB ARU613 = 0.0690
    slope2_pmb <- 0.00402; label("Slope of the linear function for the PMB adaptive-resistance on-rate, Slope2 (L/mg/h)") # Table 3: Slope2 PMB ARU613 = 0.00402
    gamma1_pmb <- 1.20; label("Exponent of the power function for the PMB kill rate, gamma1 (unitless)") # Table 3: gamma1 PMB ARU613 = 1.20
    slope3_pmb <- 1; label("Slope of the linear adaptive-resistance function for PMB, Slope3 (unitless)") # Table 3: Slope3 PMB ARU613 = 1

    # --- PMB-on-MIN interaction (Table 3, combination-specific ARU613) ------
    emax_int <- 2.45; label("Maximum fractional increase of the MIN kill rate by PMB, Emax (unitless)") # Table 3: Emax ARU613 = 2.45
    ec50_int <- 0.285; label("PMB concentration giving half of Emax on the MIN kill rate, EC50 (mg/L)") # Table 3: EC50 ARU613 = 0.285 mg/L
    slope4_int <- 8.32; label("Slope of the PMB power function on the MIN adaptive-resistance on-rate, Slope4 (unitless)") # Table 3: Slope4 ARU613 = 8.32
    gamma4_int <- 0.479; label("Exponent of the PMB power function on the MIN adaptive-resistance on-rate, gamma4 (unitless)") # Table 3: gamma4 ARU613 = 0.479

    # --- Residual error -------------------------------------------------------
    addSd <- fixed(0); label("Additive residual SD on the log10 bacterial count (log10 CFU/mL); not reported") # Not reported in Olsson 2024 (typical-value predictions only)
  })

  model({
    bmax <- 10^log10bmax

    # Eq 3: S -> R transfer grows with the total count and switches on after Tlag.
    btot <- bact_susceptible + bact_resting
    tlag <- exp(ltlag)
    ksr <- btot / bmax * (kgrowth - kdeath) * (t >= tlag)

    # Eqs 4-5: power-function kill rates. PMB raises the MIN kill rate through an
    # Emax function of the PMB concentration (Eq 5 reduces to Eq 4 when cpmb = 0).
    kdrug_pmb <- slope1_pmb * cpmb^gamma1_pmb
    kdrug_mino <- slope1_mino * cmino^gamma1_mino * (1 + emax_int * cpmb / (ec50_int + cpmb))

    # Eq 6: adaptive-resistance inhibition, beta = 1 for PMB and -1 for MIN.
    inh_pmb <- 1 - slope3_pmb * aron_pmb
    inh_mino <- 1 / (1 - slope3_mino * aron_mino)

    # Eqs 9-10: on-rates of adaptive resistance. PMB raises the MIN on-rate
    # through a power function of the PMB concentration (Eq 10).
    kon_pmb <- slope2_pmb * cpmb
    kon_mino <- slope2_mino * cmino * (1 + slope4_int * cpmb^gamma4_int)

    # Eqs 1-2: bacterial states. The growth term is gated off once the total
    # count falls below 0.01 CFU/mL (1 CFU in the 100 mL working volume): the
    # paper assumed no bacteria could regrow below that (Methods, PD
    # predictions), which is what makes the eradicated combination arms stay at
    # the LOD instead of regrowing once the drugs wash out.
    alive <- btot > 0.01
    d/dt(bact_susceptible) <- kgrowth * alive * bact_susceptible -
      (kdeath + kdrug_mino * inh_mino + kdrug_pmb * inh_pmb) * bact_susceptible -
      ksr * bact_susceptible
    d/dt(bact_resting) <- ksr * bact_susceptible - kdeath * bact_resting

    # Eqs 7-8: irreversible adaptive resistance (no return rate in Table S1).
    d/dt(aroff_pmb) <- -kon_pmb * aroff_pmb
    d/dt(aron_pmb) <- kon_pmb * aroff_pmb
    d/dt(aroff_mino) <- -kon_mino * aroff_mino
    d/dt(aron_mino) <- kon_mino * aroff_mino

    # Static time-kill bath concentrations (mg/L), dosed by the user.
    d/dt(cpmb) <- 0
    d/dt(cmino) <- 0

    bact_susceptible(0) <- 10^log10cfu0
    bact_resting(0) <- 0
    aroff_pmb(0) <- 1
    aron_pmb(0) <- 0
    aroff_mino(0) <- 1
    aron_mino(0) <- 0

    log_cfu <- log10(bact_susceptible + bact_resting)
    log_cfu ~ add(addSd)
  })
}
