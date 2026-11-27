Zhang_2016_doxycycline_killrate <- function() {
  description <- paste(
    "In vitro (Mycoplasma gallisepticum strain S6). Concentration-dependent",
    "kill-rate PK/PD model for doxycycline, coupling the one-compartment",
    "intravenous-bolus PK of Zhang 2016's in-vitro dynamic model to the",
    "sigmoid Emax kill-rate relationship fitted to the static time-kill",
    "curves (Zhang 2016 equation 1, E = Emax * Ce^N / (EC50^N + Ce^N)).",
    "E is the kill rate: the negative slope of log10(CFU/mL) against time",
    "(POSITIVE = killing). The packaged parameters are Table 1, row 0-48 h,",
    "which the paper identifies as the best fit (R^2 = 0.986) and plots in",
    "Figure 2: Emax = 0.11 1/h, EC50 = 0.32 ug/mL, Hill N = 2.11. The",
    "fitted equation has no baseline term, so at zero concentration the",
    "count is constant: the slow growth of the drug-free control",
    "(about 0.02 1/h) and the regrowth seen once the dynamic-model",
    "concentration fell below the MIC are not described. The mycoplasma",
    "density bact (linear CFU/mL) is integrated as",
    "d/dt(bact) = -ln(10) * E * bact, so log10(bact) falls at exactly E",
    "log10 units per hour and the static fit is reproduced at any constant",
    "concentration. The PK defaults to the dynamic-model arm that simulated",
    "the 6.51 h chicken intravenous half-life (Results, 'Pharmacokinetics",
    "in In Vitro Dynamic Model': CL 0.10 L/h, V 0.90 L, measured half-life",
    "6.19 h); the long-acting 12 h arm (CL 0.06 L/h, V 0.99 L) is obtained",
    "by overriding lcl and lvc. A dose of amt mg gives an initial",
    "concentration of amt / 0.90 mg/L (ug/mL); setting the clearance to",
    "a negligible value recovers the static experiment. Zhang 2016 fitted",
    "the kill-rate parameters to static concentrations and ran the dynamic",
    "model as a separate experiment, so this composition is the paper's",
    "PK/PD integration written as an ODE system, not a joint fit, and it",
    "over-predicts the killing measured in the dynamic experiments (Figure",
    "5; see the vignette): for fluctuating exposure prefer",
    "Zhang_2016_doxycycline_pkpd_index, which was fitted to them. No",
    "between-experiment variability or residual error was reported, so",
    "there are no etas and both residual SDs are fixed at zero."
  )
  reference <- paste(
    "Zhang N, Gu X, Ye X, Wu X, Zhang B, Zhang L, Shen X, Jiang H, Ding H.",
    "The PK/PD interactions of doxycycline against Mycoplasma",
    "gallisepticum. Front Microbiol. 2016;7:653.",
    "doi:10.3389/fmicb.2016.00653."
  )
  vignette <- "Zhang_2016_doxycycline"
  units <- list(
    time = "h",
    dosing = "mg (doxycycline into the external compartment; amt / V is the initial concentration in mg/L)",
    concentration = "ug/mL (Cc); log10 CFU/mL (log_cfu)"
  )

  paper_specific_compartments <- c("bact")

  compartmentData <- list(
    central = list(analyte = "doxycycline", units = "mg", specimen = "not applicable", verified = TRUE),
    bact = list(
      analyte = "Mycoplasma gallisepticum strain S6",
      units = "CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "in vitro (Mycoplasma gallisepticum strain S6 culture)",
    n_subjects = 3L,
    n_studies = 1L,
    organism = paste(
      "Mycoplasma gallisepticum standard strain S6 (China Institute of",
      "Veterinary Drugs Control, Beijing), cultured at 37 C and 5% CO2 in",
      "M. gallisepticum artificial medium base with pig serum, NADH and",
      "cysteine. Doxycycline MIC = 0.125 mg/L at inocula of 10^5, 10^6 and",
      "10^7 CFU/mL (modified Tanner and Wu colour-change microdilution)."
    ),
    system = paste(
      "Static time-kill: 7 mL penicillin bottles (0.4 mL exponential-phase",
      "culture, 0.1 mL drug solution, 3.5 mL blank medium) at constant",
      "doxycycline concentrations for 48 h. Dynamic model (PK component):",
      "a three-necked flask holding 320 mL of medium (external compartment)",
      "around a 10 mL dialysis-tube internal compartment carrying the",
      "culture, continuously diluted with drug-free medium so that the",
      "concentration declines monoexponentially (C = C0 * exp(-k * t))."
    ),
    disease_state = "Not applicable (in-vitro culture); starting inoculum approximately 10^7 CFU/mL",
    dose_range = paste(
      "Static: 0, 0.5, 1, 2, 4, 8, 16, 32 and 64 x MIC (0-8 ug/mL).",
      "Dynamic: initial concentrations of 1-64 x MIC with simulated",
      "half-lives of 6.51 h and 12 h, observed for 48 h."
    ),
    sampling = paste(
      "Static: 100 uL at 0, 2, 12, 15, 18, 24, 26, 36, 38 and 48 h, plated in",
      "triplicate, lower limit of detection 100 CFU/mL. Dynamic: external",
      "compartment (drug) and internal compartment (counts) at 2, 12, 15,",
      "18, 21, 24, 36 and 48 h. All experiments in triplicate."
    ),
    regions = "China (South China Agricultural University, Guangzhou)",
    notes = paste(
      "Kill rates were the linear-regression slopes of log10(CFU/mL) over",
      "seven windows (0-24, 0-36, 0-48, 2-24, 2-36, 2-48 and 12-48 h); the",
      "sigmoid Emax was fitted to each window in WinNonlin 6.1 (Table 1,",
      "R^2 0.907-0.986). The 0-48 h window is packaged; the other six are",
      "reproduced in the vignette. The dynamic-model PK was fitted with a",
      "one-compartment intravenous-bolus model in WinNonlin: 6.51 h design",
      "arm half-life 6.19 +/- 0.43 h, MRT 8.29 +/- 0.36 h, V 0.90 +/- 0.10 L,",
      "CL 0.10 +/- 0.01 L/h; 12 h design arm half-life 11.93 +/- 0.74 h,",
      "MRT 17.03 +/- 0.41 h, V 0.99 +/- 0.20 L, CL 0.06 +/- 0.01 L/h. The",
      "chicken intravenous PK (Table 2; 10 mg/kg, 12 Sanhuang broilers) was",
      "analysed non-compartmentally and is not packaged."
    )
  )

  ini({
    # Zhang 2016 Table 1, row 0-48 h; equation (1) E = Emax * Ce^N / (EC50^N + Ce^N),
    # E = kill rate (negative slope of log10 CFU/mL vs time), positive = killing.
    lemax <- log(0.11)
    label("Maximum kill rate Emax (log10 CFU/mL per h)") # Table 1, row '0-48', Emax = 0.11 1/h; Results 'maximum kill rate was 0.11 h-1'
    lec50 <- log(0.32)
    label("Doxycycline concentration producing half the maximum kill rate EC50 (ug/mL)") # Table 1, row '0-48', EC50 = 0.32 ug/mL
    lhill <- log(2.11)
    label("Hill coefficient N of the kill-rate curve (unitless)") # Table 1, row '0-48', Hill's slope = 2.11

    # One-compartment intravenous-bolus PK of the in-vitro dynamic model,
    # 6.51 h design arm. For the 12 h design arm use lcl = log(0.06) and
    # lvc = log(0.99) (same Results paragraph).
    lcl <- log(0.10)
    label("Clearance of the in-vitro dynamic model CL (L/h)") # Results 'Pharmacokinetics in In Vitro Dynamic Model': total body clearance 0.10 +/- 0.01 L/h (6.51 h design)
    lvc <- log(0.90)
    label("Volume of distribution of the in-vitro dynamic model V (L)") # Results 'Pharmacokinetics in In Vitro Dynamic Model': volume of distribution 0.90 +/- 0.10 L (6.51 h design)

    # Starting density, an experimental design input.
    log10_cfu0 <- fixed(7)
    label("Starting mycoplasma density (log10 CFU/mL)") # Methods 'Exposure to Constant Antibiotic Concentrations': inoculum approximately 10^7 CFU/mL

    # No residual error was reported (WinNonlin fits, R^2 only).
    addSd <- fixed(0)
    label("Additive residual SD on doxycycline concentration (ug/mL)") # not reported
    addSd_log_cfu <- fixed(0)
    label("Additive residual SD on log10 CFU/mL") # not reported
  })

  model({
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)
    cl <- exp(lcl)
    vc <- exp(lvc)
    kel <- cl / vc

    # Monoexponential decline of the dynamic model, C = C0 * exp(-k * t).
    d/dt(central) <- -kel * central
    # mg/L = ug/mL
    Cc <- central / vc

    # Equation (1): kill rate in log10 units per hour.
    kill_log10 <- emax * Cc^hill / (ec50^hill + Cc^hill)

    # d(log10 N)/dt = -kill_log10  =>  dN/dt = -ln(10) * kill_log10 * N
    d/dt(bact) <- -log(10) * kill_log10 * bact
    bact(0) <- 10^log10_cfu0

    # 1 CFU/mL floor keeps the log finite.
    log_cfu <- log10(bact + 1)

    Cc ~ add(addSd)
    log_cfu ~ add(addSd_log_cfu)
  })
}
