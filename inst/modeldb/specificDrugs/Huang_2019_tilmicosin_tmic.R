Huang_2019_tilmicosin_tmic <- function() {
  description <- "In vitro (dynamic PK/PD apparatus simulating chicken lung tissue; Mycoplasma gallisepticum strain S6). Inhibitory sigmoid Emax PK/PD-index model for the anti-mycoplasma effect of tilmicosin, driven by %T > MIC (the percentage of a 24 h interval during which the medium concentration exceeds the MIC), fitted to the dynamic time-kill data of Huang 2019. Same equation and inverted naming as the sibling Huang_2019_tilmicosin_aucmic: Huang 2019 fits E = Emax - (Emax - E0) * Ce^N / (EC50^N + Ce^N), where E is the change in log10 CFU/mL over a 24 h interval, the paper's 'Emax' is the drug-free control change and its 'E0' the largest effect. Canonical roles: e0 (zero-exposure effect) = the paper's Emax = 4.00e-6 log10 CFU/mL, emax (maximal effect) = the paper's E0 = -1.07 log10 CFU/mL, EC50 = 46.93 % and Hill N = 9.66 (Huang 2019 Table 4, T>MIC row). This is the poorer of the paper's two index fits (R^2 = 0.49 against 0.87 for AUC24h/MIC); Huang 2019 used the comparison to conclude that tilmicosin acts in a concentration-dependent rather than time-dependent way, and the fit is packaged because the paper reports it in full. The index is supplied as the covariate FTMIC_TILM on the percent scale (0-100), one value per 24 h interval. The bacterial density bact (CFU/mL) is integrated as d/dt(bact) = ln(10) * (E / 24) * bact so that log10(bact) changes by exactly E across each 24 h interval. No variability or residual error magnitude was reported, so no eta parameters are present and addSd is FIXED at 0. Sibling models from the same paper: Huang_2019_tilmicosin_invitro_pk, Huang_2019_tilmicosin_killrate and Huang_2019_tilmicosin_aucmic."
  reference <- paste(
    "Huang Z, Wu Y, Zhou Z, Xia X, Gu X, Cai Q, Shen X, Yang H, Ding H. (2019).",
    "Pharmacokinetic and Pharmacodynamic Integration and Resistance Analysis",
    "of Tilmicosin Against Mycoplasma gallisepticum in an In Vitro Dynamic Model.",
    "Frontiers in Pharmacology 10:670. doi:10.3389/fphar.2019.00670.",
    "PMCID: PMC6598723.",
    "Model equation: Materials and Methods, 'Pharmacokinetic-Pharmacodynamic",
    "Integration and Modeling'. Parameter estimates: Table 4, T>MIC row.",
    "Fitted curve: Figure 6B. Dynamic time-kill data: Figure 4.",
    sep = " "
  )
  vignette <- "Huang_2019_tilmicosin"
  units <- list(
    time = "h",
    dosing = "none (exposure enters as the FTMIC_TILM covariate, in percent)",
    concentration = "log10 CFU/mL (observation)"
  )

  depends <- c("FTMIC_TILM")

  covariateData <- list(
    FTMIC_TILM = list(
      description = "Percentage of the current 24 h interval during which the tilmicosin concentration in the central-chamber medium of the in vitro dynamic model exceeds the MIC against M. gallisepticum S6 (%T > MIC)",
      units = "percent",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying: one value per 24 h interval, on the 0-100 percent scale (EC50 = 46.93 is on the same scale). Huang 2019 does not state which of its four reported MICs defined the threshold; the agar-dilution MIC at the 10^9 CFU inoculum (0.16 ug/mL) is the one consistent with Figure 6B, since it makes the index 0 throughout for the 1 mg group whose Cmax is 0.11 ug/mL. Set to 0 for drug-free control intervals.",
      source_name = "%T > MIC (Huang 2019 Table 4 'T>MIC (%)' and Figure 6B)"
    )
  )

  paper_specific_compartments <- c("bact")

  compartmentData <- list(
    bact = list(
      analyte = "Mycoplasma gallisepticum strain S6",
      units = "CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "in vitro (dynamic PK/PD apparatus; Mycoplasma gallisepticum strain S6 culture)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    organism = "Mycoplasma gallisepticum standard strain S6 (Chinese Veterinary Microorganism Culture Collection Center, Beijing)",
    system = "In vitro dynamic model (absorption chamber, central chamber with a 10 mL dialysis-tube internal compartment, peristaltic pump at 0.13 mL/min) simulating tilmicosin lung-tissue PK in chickens; see Huang_2019_tilmicosin_invitro_pk",
    starting_inoculum = "10 mL of 2 x 10^9 CFU/mL in the internal compartment, incubated 12 h before dosing (time 0)",
    dose_range = "Drug-free control plus seven single doses (nominal clinical 1, 2.5, 5, 7.5, 10, 15 and 20 mg; 0.089 to 1.786 mg injected in vitro)",
    sampling = "100 uL from the internal compartment before dosing and at 6, 12, 24, 36, 48, 72, 96, 120, 144 and 168 h; triplicate experiments",
    disease_state = "not applicable (in vitro)",
    regions = "China (South China Agricultural University, Guangzhou)",
    notes = "Most observations in Figure 6B sit at %T > MIC = 0 or 100, which is why the fit is steep (N = 9.66) and poor (R^2 = 0.49)."
  )

  ini({
    # Huang 2019 'Pharmacokinetic-Pharmacodynamic Integration and Modeling':
    #   E = Emax - (Emax - E0) * Ce^N / (EC50^N + Ce^N)
    # with the paper's Emax = control change and E0 = largest effect,
    # written here with canonical roles (same curve).
    e0 <- 4.00e-06
    label("Change in M. gallisepticum count over 24 h at zero exposure (log10 CFU/mL)") # Huang 2019 Table 4, T>MIC row, 'E max' = 4.00 x 10^-6 (defined in Methods as the control-group change)
    emax <- -1.07
    label("Maximal change in M. gallisepticum count over 24 h (log10 CFU/mL; negative = reduction)") # Huang 2019 Table 4, T>MIC row, 'E 0' = -1.07 (defined in Methods as the largest antimycoplasmal effect)
    lec50 <- log(46.93)
    label("Log %T > MIC producing half the maximal effect EC50 (percent)") # Huang 2019 Table 4, T>MIC row, EC50 = 46.93
    lhill <- log(9.66)
    label("Log Hill coefficient N of the index-effect curve (unitless)") # Huang 2019 Table 4, T>MIC row, Hill's slope = 9.66

    log10_cfu0 <- fixed(9.30103)
    label("log10 starting mycoplasma density in the internal compartment (log10 CFU/mL)") # Huang 2019 'In Vitro Pharmacokinetic and Pharmacodynamic Model and Dosing Regimens': 2 x 10^9 CFU/ml in 10 ml of medium; log10(2e9) = 9.30103

    # Huang 2019 Table 4 reports R^2 only, no residual SD.
    addSd <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (not reported in Huang 2019)") # Huang 2019 Table 4 reports R^2 only
  })

  model({
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    ce <- FTMIC_TILM

    # Change in log10 CFU/mL over a 24 h interval.
    eff_log10 <- e0 - (e0 - emax) * ce^hill / (ce^hill + ec50^hill)

    d/dt(bact) <- log(10) * (eff_log10 / 24) * bact
    bact(0) <- 10^log10_cfu0

    log_cfu <- log10(bact + 1)
    log_cfu ~ add(addSd)
  })
}
