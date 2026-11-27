Huang_2019_tilmicosin_aucmic <- function() {
  description <- "In vitro (dynamic PK/PD apparatus simulating chicken lung tissue; Mycoplasma gallisepticum strain S6). Inhibitory sigmoid Emax PK/PD-index model for the anti-mycoplasma effect of tilmicosin, driven by the AUC24h/MIC index, fitted to the dynamic time-kill data of Huang 2019 (seven single doses into the in vitro dynamic model, counts in the internal dialysis compartment over 168 h). Huang 2019 Materials and Methods ('Pharmacokinetic-Pharmacodynamic Integration and Modeling') fits E = Emax - (Emax - E0) * Ce^N / (EC50^N + Ce^N) in WinNonlin, where E is the change in log10 CFU/mL over a 24 h interval, Ce the PK/PD index, and - in the paper's INVERTED naming - 'Emax' is the change in the drug-free control and 'E0' the largest antimycoplasmal effect. The packaged model uses canonical roles: e0 (effect at zero exposure) = the paper's Emax = 2.98e-4 log10 CFU/mL and emax (maximal effect) = the paper's E0 = -3.35 log10 CFU/mL, with EC50 = 849.68 h and Hill N = 0.33 (Huang 2019 Table 4, AUC24h/MIC row, R^2 = 0.87); only the labels change, not the values. Inverting the fitted curve at E = -1 gives AUC24h/MIC = 63.9 h against the paper's printed 62.58 h (the 2% gap is consistent with the two-decimal rounding of N). The index is supplied as the already-formed ratio covariate AUCMIC_TILM, per 24 h interval, because Huang 2019 reports four different MIC values for strain S6 and never states which one formed the index. The bacterial density bact (CFU/mL) is integrated as d/dt(bact) = ln(10) * (E / 24) * bact so that log10(bact) changes by exactly E across each 24 h interval. The fitted curve is bounded between e0 (about 0) and emax (-3.35), so it cannot describe the regrowth Huang 2019 Figure 4 shows after 120 h. No variability or residual error magnitude was reported, so no eta parameters are present and addSd is FIXED at 0. Sibling models from the same paper: Huang_2019_tilmicosin_invitro_pk, Huang_2019_tilmicosin_killrate and Huang_2019_tilmicosin_tmic (the same structure fitted against %T > MIC)."
  reference <- paste(
    "Huang Z, Wu Y, Zhou Z, Xia X, Gu X, Cai Q, Shen X, Yang H, Ding H. (2019).",
    "Pharmacokinetic and Pharmacodynamic Integration and Resistance Analysis",
    "of Tilmicosin Against Mycoplasma gallisepticum in an In Vitro Dynamic Model.",
    "Frontiers in Pharmacology 10:670. doi:10.3389/fphar.2019.00670.",
    "PMCID: PMC6598723.",
    "Model equation: Materials and Methods, 'Pharmacokinetic-Pharmacodynamic",
    "Integration and Modeling'. Parameter estimates: Table 4, AUC24h/MIC row.",
    "Fitted curve: Figure 6A. Dynamic time-kill data: Figure 4.",
    sep = " "
  )
  vignette <- "Huang_2019_tilmicosin"
  units <- list(
    time = "h",
    dosing = "none (exposure enters as the AUCMIC_TILM covariate, in h)",
    concentration = "log10 CFU/mL (observation)"
  )

  depends <- c("AUCMIC_TILM")

  covariateData <- list(
    AUCMIC_TILM = list(
      description = "Tilmicosin AUC over the current 24 h interval in the central-chamber medium of the in vitro dynamic model, divided by the MIC against M. gallisepticum S6 (AUC24h/MIC)",
      units = "h",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying: one value per 24 h interval (0-24, 24-48, ... h), held constant within the interval. Carried as a ratio because Huang 2019 reports four MICs for strain S6 (0.01 and 0.02 ug/mL by microdilution; 0.02, 0.04 and 0.16 ug/mL by agar dilution) and never states which one formed the index. The agar-dilution MIC at the 10^9 CFU inoculum (0.16 ug/mL), which matches the dynamic model's 2 x 10^9 CFU/mL inoculum, is the one consistent with the paper's figures: with it the per-interval AUC24h/MIC from the Table 2 PK spans 0 to about 340 h (Figure 6A axis 0-350), and %T > MIC is 0 for the 1 mg group whose Cmax (0.11 ug/mL) is below 0.16 (the cluster of zero-index points in Figure 6B). That inference is illustrated in the vignette; it is not built into the model. Set to 0 for drug-free control intervals.",
      source_name = "AUC24h/MIC (Huang 2019 Table 4 and Figure 6A)"
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
    notes = "Huang 2019 identifies AUC24h/MIC as the PK/PD index best correlated with the effect (R^2 = 0.87 versus 0.49 for %T > MIC) and concludes tilmicosin is concentration-dependent against M. gallisepticum. The effect is the change over 24 h intervals of a culture already near its maximal density: the drug-free control stayed at about 9 log10 CFU/mL over 168 h (Figure 4), which is why the zero-exposure effect is close to 0."
  )

  ini({
    # Huang 2019 'Pharmacokinetic-Pharmacodynamic Integration and Modeling':
    #   E = Emax - (Emax - E0) * Ce^N / (EC50^N + Ce^N)
    # with the paper's Emax = control change and E0 = largest effect.
    # Written here with canonical roles (e0 = zero-exposure effect,
    # emax = maximal effect), which is the same curve.
    e0 <- 2.98e-04
    label("Change in M. gallisepticum count over 24 h at zero exposure (log10 CFU/mL)") # Huang 2019 Table 4, AUC24h/MIC row, 'E max' = 2.98 x 10^-4 (defined in Methods as the control-group change)
    emax <- -3.35
    label("Maximal change in M. gallisepticum count over 24 h (log10 CFU/mL; negative = reduction)") # Huang 2019 Table 4, AUC24h/MIC row, 'E 0' = -3.35 (defined in Methods as the largest antimycoplasmal effect)
    lec50 <- log(849.68)
    label("Log AUC24h/MIC producing half the maximal effect EC50 (h)") # Huang 2019 Table 4, AUC24h/MIC row, EC50 = 849.68
    lhill <- log(0.33)
    label("Log Hill coefficient N of the index-effect curve (unitless)") # Huang 2019 Table 4, AUC24h/MIC row, Hill's slope = 0.33

    # Experimental design input: the 0 h state of the internal
    # compartment, 12 h after inoculation with 2 x 10^9 CFU/mL.
    log10_cfu0 <- fixed(9.30103)
    label("log10 starting mycoplasma density in the internal compartment (log10 CFU/mL)") # Huang 2019 'In Vitro Pharmacokinetic and Pharmacodynamic Model and Dosing Regimens': 2 x 10^9 CFU/ml in 10 ml of medium; log10(2e9) = 9.30103

    # Huang 2019 Table 4 reports R^2 only, no residual SD.
    addSd <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (not reported in Huang 2019)") # Huang 2019 Table 4 reports R^2 only
  })

  model({
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    ce <- AUCMIC_TILM

    # Change in log10 CFU/mL over a 24 h interval.
    eff_log10 <- e0 - (e0 - emax) * ce^hill / (ce^hill + ec50^hill)

    # Spread evenly over the interval so log10(bact) moves by exactly
    # eff_log10 between 24 h boundaries.
    d/dt(bact) <- log(10) * (eff_log10 / 24) * bact
    bact(0) <- 10^log10_cfu0

    log_cfu <- log10(bact + 1)
    log_cfu ~ add(addSd)
  })
}
