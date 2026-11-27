Zhang_2016_doxycycline_pkpd_index <- function() {
  description <- paste(
    "In vitro (Mycoplasma gallisepticum strain S6). Inhibitory sigmoid",
    "Imax PK/PD-index model for the 48 h anti-mycoplasma effect of",
    "doxycycline in Zhang 2016's in-vitro dynamic model, driven by the",
    "paper's best-correlating index, the percentage of time the",
    "concentration exceeds the MIC (%T>MIC, R^2 = 0.989 versus 0.953 for",
    "AUC0-48h/MIC and 0.897 for Cmax/MIC). The one-compartment",
    "intravenous-bolus PK of the dynamic model (Results, 'Pharmacokinetics",
    "in In Vitro Dynamic Model') defaults to the arm that simulated the",
    "6.51 h chicken half-life (CL 0.10 L/h, V 0.90 L); the 12 h arm",
    "(CL 0.06 L/h, V 0.99 L) is obtained by overriding lcl and lvc. The",
    "state t_gt_mic accumulates the time above the MIC and pct_tmic =",
    "100 * t_gt_mic / 48 is the index. Zhang 2016 counted the whole",
    "single-dose time above the MIC against the 48 h period WITHOUT",
    "capping it at 48 h (Figure 6C has observed indices of about 120% and",
    "158%, which only the uncapped 12 h-arm values of 124% and 149% at",
    "32 and 64 x MIC can produce), so t_gt_mic keeps accumulating after",
    "48 h and the effect should be read once the concentration has fallen",
    "below the MIC. Equation (2), E = E0 - Imax * Ce^N / (IC50^N + Ce^N),",
    "with Table 3 row 'T > MIC (%)': E0 = 2.25, Imax = 5.93 log10 CFU/mL,",
    "IC50 = 44.39%, N = 1.12. dlog10cfu is the signed change in log10",
    "CFU/mL over the 48 h experiment (positive = net growth; E0 is the",
    "growth of the drug-free control) and the observation log_cfu adds it",
    "to the nominal 7 log10 CFU/mL inoculum. No residual error was",
    "reported, so addSd is fixed at zero."
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

  paper_specific_compartments <- c("t_gt_mic")

  compartmentData <- list(
    central = list(analyte = "doxycycline", units = "mg", specimen = "not applicable", verified = TRUE),
    t_gt_mic = list(analyte = "doxycycline", units = "h", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "in vitro (Mycoplasma gallisepticum strain S6 culture)",
    n_subjects = 3L,
    n_studies = 1L,
    organism = paste(
      "Mycoplasma gallisepticum standard strain S6 (China Institute of",
      "Veterinary Drugs Control, Beijing); doxycycline MIC = 0.125 mg/L at",
      "inocula of 10^5-10^7 CFU/mL."
    ),
    system = paste(
      "In-vitro dynamic model (after Meletiadis 2012): a three-necked flask",
      "holding 320 mL of medium (external compartment) around a 10 mL",
      "dialysis-tube internal compartment inoculated with about 10^7 CFU/mL,",
      "continuously diluted with drug-free medium so that the concentration",
      "declines monoexponentially. Drug was injected into both compartments",
      "at time zero."
    ),
    disease_state = "Not applicable (in-vitro culture); starting inoculum approximately 10^7 CFU/mL",
    dose_range = paste(
      "Single initial concentrations of 1, 2, 4, 8, 16, 32 and 64 x MIC",
      "(0.125-8 ug/mL) in each of two arms simulating half-lives of 6.51 h",
      "(the chicken intravenous value) and 12 h (a long-acting parenteral",
      "formulation), plus a drug-free control; observed for 48 h."
    ),
    sampling = "External (drug) and internal (counts) compartments at 2, 12, 15, 18, 21, 24, 36 and 48 h; triplicate experiments",
    regions = "China (South China Agricultural University, Guangzhou)",
    notes = paste(
      "Indices were computed by the authors from one-compartment WinNonlin",
      "fits. Table 3 also gives the rejected AUC0-48h/MIC (Imax 5.82, IC50",
      "79.12 h, E0 2.30, N 0.55, R^2 0.953) and Cmax/MIC (Imax 6.49, IC50",
      "10.90, E0 2.28, N 0.44, R^2 0.897) fits, reproduced in the vignette.",
      "The paper's stated %T>MIC targets for 0, 2 and 3 log10 CFU/mL",
      "reduction (32.48, 45.68 and 54.36%) are not reproduced by Table 3,",
      "whereas Table 3 does reproduce the fitted curve of Figure 6C; see",
      "the vignette Errata."
    )
  )

  ini({
    # One-compartment intravenous-bolus PK of the in-vitro dynamic model,
    # 6.51 h design arm. For the 12 h design arm use lcl = log(0.06) and
    # lvc = log(0.99) (same Results paragraph).
    lcl <- log(0.10)
    label("Clearance of the in-vitro dynamic model CL (L/h)") # Results 'Pharmacokinetics in In Vitro Dynamic Model': total body clearance 0.10 +/- 0.01 L/h (6.51 h design)
    lvc <- log(0.90)
    label("Volume of distribution of the in-vitro dynamic model V (L)") # Results 'Pharmacokinetics in In Vitro Dynamic Model': volume of distribution 0.90 +/- 0.10 L (6.51 h design)

    # Measured property of the strain, not an estimate.
    mic <- fixed(0.125)
    label("Doxycycline MIC against M. gallisepticum S6 (ug/mL)") # Results 'Minimum Inhibitory Concentration': MIC = 0.125 mg/L

    # Equation (2) E = E0 - Imax * Ce^N / (IC50^N + Ce^N), Table 3 row 'T > MIC (%)'.
    e0 <- 2.25
    label("Change in log10 CFU/mL of the drug-free control over 48 h, E0 (log10 CFU/mL)") # Table 3, row 'T > MIC (%)', E0 = 2.25 log10 CFU/mL
    limax <- log(5.93)
    label("Maximum anti-mycoplasma effect Imax (log10 CFU/mL reduction from E0)") # Table 3, row 'T > MIC (%)', Imax = 5.93 log10 CFU/mL
    lic50 <- log(44.39)
    label("%T>MIC producing half the maximum effect IC50 (percent of 48 h)") # Table 3, row 'T > MIC (%)', IC50 = 44.39 (unit column misprinted as ug/mL)
    lhill <- log(1.12)
    label("Hill coefficient N of the %T>MIC effect curve (unitless)") # Table 3, row 'T > MIC (%)', Hill's slope = 1.12

    # Starting density, an experimental design input.
    log10_cfu0 <- fixed(7)
    label("Starting mycoplasma density in the internal compartment (log10 CFU/mL)") # Methods 'Exposure to Changing Antibiotic Concentrations': 10 mL medium containing approximately 10^7 CFU/mL

    addSd <- fixed(0)
    label("Additive residual SD on log10 CFU/mL") # not reported (R^2 only)
  })

  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    kel <- cl / vc
    imax <- exp(limax)
    ic50 <- exp(lic50)
    hill <- exp(lhill)

    # Monoexponential decline of the dynamic model, C = C0 * exp(-k * t).
    d/dt(central) <- -kel * central
    # mg/L = ug/mL
    Cc <- central / vc

    # Time above the MIC, not capped at the 48 h experiment (see description).
    d/dt(t_gt_mic) <- (Cc > mic)
    pct_tmic <- 100 * t_gt_mic / 48

    # Equation (2): signed change in log10 CFU/mL over 48 h.
    dlog10cfu <- e0 - imax * pct_tmic^hill / (ic50^hill + pct_tmic^hill)
    log_cfu <- log10_cfu0 + dlog10cfu
    log_cfu ~ add(addSd)
  })
}
