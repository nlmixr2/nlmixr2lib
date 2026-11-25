Zhang_2019_doxycycline_cmaxmic <- function() {
  description <- "In vitro (Mycoplasma hyopneumoniae strain ATCC25934). Inhibitory sigmoid Emax PK/PD-index model for the anti-mycoplasma effect of doxycycline against M. hyopneumoniae ATCC25934 in Zhang 2019's in-vitro dynamic model, driven by the Cmax/MIC index. Zhang 2019 Materials and Methods (Integration and modeling of PK/PD) parameterises the effect over the 24 h after each administration as E = Emax - (Emax - E0) * Ce^N / (EC50^N + Ce^N), where E is the SIGNED change in mycoplasma count in log10(CFU/mL) (NEGATIVE = reduction) and Ce is the PK/PD index. THE PAPER'S LABELS ARE REVERSED RELATIVE TO THIS PACKAGE: Zhang 2019 defines its Emax as the change in the untreated control and its E0 as the maximum anti-mycoplasma effect. This file maps the paper's Emax to e0 (the effect at zero exposure) and the paper's E0 to emax (the effect at saturating exposure), so the model equation reads e0 - (e0 - emax) * f, identical in form to the printed one. Parameters from Zhang 2019 Table 2, Cmax/MIC column: paper Emax = -5e-05 log10(CFU/mL) (here e0), paper E0 = -3.14 log10(CFU/mL) (here emax), EC50 = 21, Hill N = 1.01 (R2 = 0.923, the highest of the three indices the paper examined). Substituting these at the paper's headline target Cmax/MIC = 9.89 returns E = -1.00 log10(CFU/mL), confirming both the equation form and the label mapping. The index is formed as the per-interval covariate CMAX_DOXYCYCLINE divided by the parameter mic = 0.125 ug/mL. The bacterial density bacteria (linear CFU/mL) is integrated as d/dt(bacteria) = ln(10) * (E / 24) * bacteria so that log10(bacteria) changes by exactly E across each 24 h interval. There is no PK compartment: Zhang 2019 read the per-interval Cmax directly from the measured concentration-time curve of the dynamic apparatus. Sibling models from the same paper: Zhang_2019_doxycycline_aucmic (the same structure on AUC24h/MIC) and Zhang_2019_doxycycline_killrate (static kill rate composed with the apparatus PK). Zhang 2019 reports no between-replicate variability and no residual error magnitude, so no eta parameters are present and addSd is held at zero for typical-value simulation."
  reference <- paste(
    "Zhang H, Mao C, Li J, Huang Z, Gu X, Shen X, Ding H. (2019).",
    "Pharmacokinetic/pharmacodynamic integration of doxycycline against",
    "Mycoplasma hyopneumoniae in an in vitro model.",
    "Frontiers in Pharmacology 10:1088.",
    "doi:10.3389/fphar.2019.01088. PMCID: PMC6763577.",
    sep = " "
  )
  vignette <- "Zhang_2019_doxycycline"
  units <- list(
    time = "h",
    dosing = "ug/mL (doxycycline peak concentration per dosing interval, supplied as a covariate)",
    concentration = "log10 CFU/mL (observation)"
  )

  depends <- c("CMAX_DOXYCYCLINE")

  compartmentData <- list(
    bacteria = list(
      analyte = "Mycoplasma hyopneumoniae strain ATCC25934",
      units = "CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    CMAX_DOXYCYCLINE = list(
      description = "Doxycycline peak medium concentration after the current dose, held over that 24 h dosing interval",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Per-24-h-interval peak doxycycline concentration in the central chamber of the in-vitro dynamic model, supplied as a piecewise-constant time-varying covariate (one value per 24 h interval). Zhang 2019 read it directly from the measured curve and reports the range 0.16-3.41 ug/mL across all regimens and intervals (Results, Pharmacokinetics in the in vitro dynamic model). Set to 0 for the untreated control so the sigmoid term vanishes and the predicted 24 h change equals e0. The sibling model Zhang_2019_doxycycline_killrate carries the apparatus PK, from which this Cmax can be simulated.",
      source_name = "Cmax (Zhang 2019 Materials and Methods, Intergration and modeling of PK/PD; Table 2)"
    )
  )

  population <- list(
    species = "in vitro (Mycoplasma hyopneumoniae strain ATCC25934 culture)",
    n_subjects = 3L,
    n_studies = 1L,
    organism = "Mycoplasma hyopneumoniae standard strain ATCC25934 (Chinese Veterinary Microorganism Culture Collection Center, Beijing); doxycycline MIC = 0.125 ug/mL by broth microdilution at a 10^6 CFU/mL inoculum",
    system = "In-vitro dynamic model: a 210 mL absorption chamber feeding a central chamber of 300 mL medium (external compartment) holding a 10 mL dialysis tube (internal compartment) with a 10^7 CFU/mL culture, supplied from a fresh-medium reservoir and draining to a waste vessel by peristaltic pumps at 0.58 mL/min, designed to reproduce the porcine plasma doxycycline absorption half-life of 4.15 h and elimination half-life of 5.9 h",
    disease_state = "Not applicable (in-vitro culture)",
    dose_range = "Eight regimens labelled 1, 2.5, 5, 7.5, 10, 15, 20 and 30 mg/kg body weight, once daily for 3 days into the absorption chamber; per-interval AUC24h 1.78-48.4 h*ug/mL and Cmax 0.16-3.41 ug/mL",
    sampling = "Drug: 0 (pre-dose), 1, 2, 4, 6, 12, 24, 30, 36, 48, 54, 60, 72, 84 and 96 h. Counts: 0 (pre-dose), 6, 12, 24, 30, 36, 48, 54, 60, 72, 84 and 96 h. Each regimen in triplicate",
    regions = "China (South China Agricultural University, Guangzhou)",
    notes = "Zhang 2019 Table 2 also fits the same model against AUC24h/MIC (packaged as Zhang_2019_doxycycline_aucmic). %T>MIC was also examined (Figure 6C, R2 = 0.823) but no parameters were reported for it. The paper's headline target for a 1 log10 CFU/mL reduction is Cmax/MIC = 9.89."
  )

  ini({
    # Zhang 2019 equation: E = Emax - (Emax - E0) * Ce^N / (EC50^N + Ce^N),
    # with Emax = change in the untreated control and E0 = maximum effect.
    # Mapped here as e0 = paper Emax, emax = paper E0. Both are signed
    # changes in log10 CFU/mL, so both stay on the natural scale.
    e0 <- -5e-05
    label("Change in mycoplasma count at zero exposure over 24 h e0 (log10 CFU/mL)") # Zhang 2019 Table 2, Cmax/MIC column, row 'E max' = -5x10-05 (paper's label for the control change)
    emax <- -3.14
    label("Maximal change in mycoplasma count over 24 h emax (log10 CFU/mL; negative = reduction)") # Zhang 2019 Table 2, Cmax/MIC column, row 'E 0' = -3.14 (paper's label for the maximum effect)
    lec50 <- log(21)
    label("Cmax/MIC producing half the maximal effect EC50 (unitless)") # Zhang 2019 Table 2, Cmax/MIC column, EC50 = 21
    lhill <- log(1.01)
    label("Hill coefficient N defining the slope of the fitted curve (unitless)") # Zhang 2019 Table 2, Cmax/MIC column, Hill's slope = 1.01

    # Measured property of the challenge strain, not estimated.
    mic <- fixed(0.125)
    label("Doxycycline MIC against M. hyopneumoniae ATCC25934 (ug/mL; broth microdilution)") # Zhang 2019 Results, Minimum inhibitory concentration: 0.125 ug/ml

    # Experimental design input, not an estimated parameter.
    log10_cfu0 <- fixed(7)
    label("log10 starting mycoplasma density in the dialysis tube (log10 CFU/mL)") # Zhang 2019 Exposure to dynamic antibiotic concentrations: inoculum of 10^7 CFU/ml

    addSd <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (not reported in Zhang 2019)") # Zhang 2019 Table 2 reports R2 only, no residual SD
  })

  model({
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # PK/PD index: per-interval peak concentration divided by the MIC.
    cmaxmic <- CMAX_DOXYCYCLINE / mic

    # Signed change in log10(CFU/mL) accrued over the current 24 h interval.
    eff_log10 <- e0 - (e0 - emax) * cmaxmic^hill / (ec50^hill + cmaxmic^hill)

    # Spread the interval change uniformly so log10(bacteria) changes by
    # exactly eff_log10 over each 24 h interval.
    d/dt(bacteria) <- log(10) * (eff_log10 / 24) * bacteria
    bacteria(0) <- 10^log10_cfu0

    log_cfu <- log10(bacteria + 1)
    log_cfu ~ add(addSd)
  })
}
