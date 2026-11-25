Huang_2019_tilmicosin_killrate <- function() {
  description <- "In vitro (Mycoplasma gallisepticum strain S6). Sigmoid Emax kill-rate model for tilmicosin against M. gallisepticum fitted to static time-kill curves at constant concentrations of 0.5 to 64 x MIC (0.01 to 1.28 ug/mL, MIC 0.02 ug/mL at a 10^7 CFU/mL inoculum). Huang 2019 Materials and Methods ('Time-Killing Curve Fitting and Analysis') estimates the kill rate at each concentration as the slope of the log10 count against time over five early windows (0-6, 0-9, 1-6, 1-9 and 3-9 h), averages them, and fits E = E0 + (Emax - E0) * Ce^N / (EC50^N + Ce^N) in WinNonlin, where E is the kill rate in log10 CFU/mL per hour (POSITIVE = killing), E0 the value at zero concentration (negative = net growth) and Emax the saturating kill rate (the paper's 'maximum kill rate of 0.53 h-1'). Parameters from Huang 2019 Table 1: Emax = 0.53, EC50 = 0.12 ug/mL, E0 = -0.03, Hill N = 0.73 (R^2 = 0.99). The bacterial density bact (CFU/mL) is integrated as d/dt(bact) = -ln(10) * E * bact so that log10(bact) changes at exactly -E log10 units per hour, which is the definition of the fitted kill rate. Exposure is static: the assay was a single drug addition to a sealed penicillin bottle, so the `tilmicosin` state holds the bath concentration in ug/mL with d/dt(tilmicosin) = 0. The fit describes the early (0-9 h) decline only: the model has no limit-of-detection floor (100 CFU/mL in the assay) and no carrying capacity, so it should not be used to extrapolate the counts past the fitting window. No variability or residual error magnitude was reported, so no eta parameters are present and addSd is FIXED at 0. Sibling models from the same paper: Huang_2019_tilmicosin_invitro_pk, Huang_2019_tilmicosin_aucmic and Huang_2019_tilmicosin_tmic."
  reference <- paste(
    "Huang Z, Wu Y, Zhou Z, Xia X, Gu X, Cai Q, Shen X, Yang H, Ding H. (2019).",
    "Pharmacokinetic and Pharmacodynamic Integration and Resistance Analysis",
    "of Tilmicosin Against Mycoplasma gallisepticum in an In Vitro Dynamic Model.",
    "Frontiers in Pharmacology 10:670. doi:10.3389/fphar.2019.00670.",
    "PMCID: PMC6598723.",
    "Model equation: Materials and Methods, 'Time-Killing Curve Fitting and",
    "Analysis'. Parameter estimates: Table 1. Fitted curve: Figure 3.",
    "Time-kill design: 'Time-Kill Curve Studies'; time-kill data: Figure 2.",
    sep = " "
  )
  vignette <- "Huang_2019_tilmicosin"
  units <- list(
    time = "h",
    dosing = "ug/mL (tilmicosin bath concentration placed directly into the static `tilmicosin` state at time 0)",
    concentration = "log10 CFU/mL (observation)"
  )

  # The dosing target is a static in vitro bath rather than `depot` or
  # `central`: the dose is the bath concentration itself.
  dosing <- c("tilmicosin")

  paper_specific_compartments <- c("tilmicosin", "bact")

  compartmentData <- list(
    tilmicosin = list(
      analyte = "tilmicosin",
      units = "ug/mL",
      specimen = "administration site",
      verified = TRUE
    ),
    bact = list(
      analyte = "Mycoplasma gallisepticum strain S6",
      units = "CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "in vitro (Mycoplasma gallisepticum strain S6 culture)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    organism = "Mycoplasma gallisepticum standard strain S6 (Chinese Veterinary Microorganism Culture Collection Center, Beijing). Tilmicosin MIC 0.01 and 0.02 ug/mL by broth microdilution at 10^5 and 10^7 CFU/mL inocula, and 0.02, 0.04 and 0.16 ug/mL by agar dilution at 10^5, 10^7 and 10^9 CFU inocula; the 10^7 CFU/mL microdilution value (0.02 ug/mL) sets the x MIC multiples of the time-kill study",
    system = "Static time-kill curves in penicillin bottles (3.86 mL blank medium + 0.04 mL of 10x drug + 0.1 mL logarithmic-phase culture), 37 C, 5% CO2, 48 h",
    starting_inoculum = "10^7 CFU/mL",
    limit_of_detection = "100 CFU/mL (2 log10 CFU/mL)",
    concentration_range = "0 (growth control) and 0.5, 1, 2, 4, 8, 16, 32 and 64 x MIC (0.01 to 1.28 ug/mL)",
    sampling = "100 uL aliquots at 0, 1, 3, 6, 9, 12, 24, 36 and 48 h, ten-fold serial dilution and plating; triplicate experiments",
    disease_state = "not applicable (in vitro)",
    regions = "China (South China Agricultural University, Guangzhou)",
    notes = "Huang 2019 Figure 3 plots the average kill rate at the eight drug-containing concentrations (0.01 to 1.28 ug/mL), rising from about 0.05 to 0.44 log10 CFU/mL per hour; the Table 1 curve reproduces these to the precision the figure can be read."
  )

  ini({
    # Huang 2019 'Time-Killing Curve Fitting and Analysis':
    #   E = E0 + (Emax - E0) * Ce^N / (EC50^N + Ce^N)
    # E is the kill rate in log10 CFU/mL per hour (positive = killing).
    # In this form Emax is the saturating kill rate itself, matching
    # the Results statement 'The maximum kill rate was 0.53 h-1'.
    #
    # e0 stays on the natural scale because it is negative (net growth
    # of the drug-free culture).
    e0 <- -0.03
    label("Kill rate at zero tilmicosin concentration E0 (log10 CFU/mL per h; negative = net growth)") # Huang 2019 Table 1, E 0 = -0.03
    lemax <- log(0.53)
    label("Log maximum kill rate Emax (log10 CFU/mL per h)") # Huang 2019 Table 1, E max = 0.53
    lec50 <- log(0.12)
    label("Log tilmicosin concentration producing half the maximal change in kill rate EC50 (ug/mL)") # Huang 2019 Table 1, EC50 = 0.12; concentration axis of Figure 3 in ug/mL
    lhill <- log(0.73)
    label("Log Hill coefficient N of the concentration-kill-rate curve (unitless)") # Huang 2019 Table 1, Hill's slope = 0.73

    # Experimental design input, not an estimated parameter.
    log10_cfu0 <- fixed(7)
    label("log10 starting mycoplasma density in the time-kill bottle (log10 CFU/mL)") # Huang 2019 'Time-Kill Curve Studies': final inoculation 10^7 CFU/ml

    # Huang 2019 Table 1 reports R^2 = 0.99 only, no residual SD.
    addSd <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (not reported in Huang 2019)") # Huang 2019 Table 1 reports R^2 only
  })

  model({
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # Static in vitro exposure: a bolus of `amt` into `tilmicosin` sets
    # the bath concentration in ug/mL, which then stays constant.
    d/dt(tilmicosin) <- 0

    kill_log10 <- e0 + (emax - e0) * tilmicosin^hill / (ec50^hill + tilmicosin^hill)

    # d(log10 N)/dt = -kill_log10  =>  dN/dt = -ln(10) * kill_log10 * N
    d/dt(bact) <- -log(10) * kill_log10 * bact
    bact(0) <- 10^log10_cfu0

    # 1 CFU/mL floor keeps the log finite if bact falls below 1 CFU/mL.
    log_cfu <- log10(bact + 1)
    log_cfu ~ add(addSd)
  })
}
