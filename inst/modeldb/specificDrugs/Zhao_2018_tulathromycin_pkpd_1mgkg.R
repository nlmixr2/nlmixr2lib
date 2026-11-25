Zhao_2018_tulathromycin_pkpd_1mgkg <- function() {
  description <- paste0(
    "Ex vivo (guinea pig serum). Inhibitory sigmoid Emax PK/PD-integration model for the ",
    "antibacterial effect of tulathromycin against Haemophilus parasuis serotype 13 reference ",
    "strain 13R incubated for 24 h in SERUM sampled from H. parasuis-infected neutropenic guinea ",
    "pigs after a single 1 mg/kg intramuscular dose. Zhao 2018 equation (PK/PD modeling and data ",
    "analysis): E = Emax - (Emax - E0) * Ce^N / (EC50^N + Ce^N), where E is the SIGNED change in ",
    "log10 CFU/mL over the 24 h incubation (positive = net growth, negative = net kill) and Ce is ",
    "the PK/PD index AUC0-24h/MIC (h). NOTE THE INVERTED NAMING: Zhao 2018's Emax (+2.65) is the ",
    "change in drug-free control serum and its E0 (-4.01) is the maximal kill; this file maps them ",
    "onto the canonical e0 (effect at zero exposure) and emax (maximal drug effect) without ",
    "changing any value. Parameters from Zhao 2018 Table 6, column '1 mg/kg': e0 = 2.65, ",
    "emax = -4.01 log10 CFU/mL, EC50 = 24.66 h, N = 5.11. The index is formed in model() as ",
    "AUC_TULATHROMYCIN / mic with the strain 13R serum MIC of 0.03 ug/mL. In the ex vivo design ",
    "each serum sample holds one fixed concentration C for 24 h, so AUC_TULATHROMYCIN = 24 * C. ",
    "The parameterisation reproduces the paper's independently printed targets (Table 6: ",
    "bacteriostatic 22.73, bactericidal 34.52, eradication 88.03 h) to within 0.1%. The bacterial ",
    "density bact (CFU/mL) is integrated as d/dt(bact) = ln(10) * (E / 24) * bact so log10(bact) ",
    "changes by exactly E over each 24 h window. No PK component: Zhao 2018 analysed the serum ",
    "concentrations non-compartmentally (WinNonlin model 200), so exposure enters as the ",
    "covariate AUC_TULATHROMYCIN. Neither between-animal variability nor a residual error was ",
    "reported, so there are no etas and addSd is fixed at 0. Sibling: ",
    "Zhao_2018_tulathromycin_pkpd_10mgkg (the separate fit to serum from the 10 mg/kg dose)."
  )
  reference <- paste(
    "Zhao Y, Guo L-L, Fang B, Liu B.",
    "Pharmacokinetic/pharmacodynamic (PK/PD) evaluation of tulathromycin against",
    "Haemophilus parasuis in an experimental neutropenic guinea pig model.",
    "PLoS One. 2018;13(12):e0209177. doi:10.1371/journal.pone.0209177. PMCID: PMC6312216.",
    "Model equation from Materials and methods (PK/PD modeling and data analysis);",
    "parameter values from Table 6, column '1 mg/kg'; serum MIC of strain 13R from Results",
    "(MICs, MBCs, MPC and PAE).",
    sep = " "
  )
  vignette <- "Zhao_2018_tulathromycin"

  units <- list(
    time = "h",
    dosing = "h*ug/mL (tulathromycin serum AUC over a 24 h window, supplied as a covariate)",
    concentration = "log10 CFU/mL (observation)"
  )

  depends <- c("AUC_TULATHROMYCIN")
  paper_specific_compartments <- c("bact")

  compartmentData <- list(
    bact = list(
      analyte = "Haemophilus parasuis strain 13R",
      units = "CFU/mL",
      specimen = "serum",
      verified = TRUE
    )
  )

  covariateData <- list(
    AUC_TULATHROMYCIN = list(
      description = "Tulathromycin serum area under the concentration-time curve over a 24 h exposure window (AUC0-24h); in the ex vivo design this is the static sample concentration times 24 h",
      units = "h*ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Zhao 2018 drives its sigmoid Emax model with Ce = AUC0-24h/MIC (h). Each ex vivo serum ",
        "sample holds one fixed tulathromycin concentration C (ug/mL) for the 24 h incubation, so ",
        "AUC0-24h = 24 * C; Figure 7 confirms the convention, its AUC0-24h/MIC axis being 24 times ",
        "its Cmax/MIC axis point for point. Zhao 2018 reports serum Cmax 1079.5 ng/mL at 0.5 h after ",
        "1 mg/kg (Table 2). Set to 0 for drug-free control serum so the sigmoid term vanishes and the ",
        "predicted 24 h change equals e0."
      ),
      source_name = "AUC0-24h (numerator of the PK/PD index AUC0-24h/MIC serum, Zhao 2018 Methods and Table 6)"
    )
  )

  population <- list(
    species = "ex vivo (serum from neutropenic guinea pig) challenged with Haemophilus parasuis",
    n_subjects = 128L,
    n_studies = 1L,
    age_range = "three weeks old at purchase",
    weight_range = "240-250 g at purchase",
    sex_female_pct = NA_real_,
    disease_state = "Cyclophosphamide-induced neutropenia (100 mg/kg/day intraperitoneally for 3 days; leukocytes < 1000/mm^3) and intraperitoneal Haemophilus parasuis 13R infection (about 10^9 CFU, the ID95)",
    dose_range = "Single intramuscular dose of tulathromycin (Draxxin) 1 mg/kg in the serum donors",
    regions = "China (South China Agricultural University, Guangzhou)",
    organism = "Haemophilus parasuis serotype 13 reference strain 13R; tulathromycin MIC 0.03 ug/mL and MBC 0.06 ug/mL in guinea pig serum (0.5 and 1 ug/mL in CAMHB)",
    notes = paste0(
      "FMMU albino guinea pigs. Eight animals were sacrificed per time point at 16 time points ",
      "(0.083-168 h), so the 1 mg/kg arm used 128 animals; the ex vivo time-kill curves average ",
      "the eight sera of each time point. Ex vivo time-kill: strain 13R incubated in each serum sample ",
      "for 24 h with viable counts at 0, 2, 5, 8 and 24 h (Figure 5); the counts start at about ",
      "10^5.7 CFU/mL. Free fraction in serum 0.56-0.74. Table 6 targets: AUC0-24h/MIC = 22.73 h ",
      "(bacteriostatic, E = 0), 34.52 h (bactericidal, E = -3) and 88.03 h (virtual eradication, ",
      "E = -4)."
    )
  )

  ini({
    # Zhao 2018 equation (Materials and methods, PK/PD modeling and data analysis):
    #   E = Emax - (Emax - E0) * Ce^N / (EC50^N + Ce^N)
    # with Zhao 2018's Emax the drug-free control change (+2.65) and its E0
    # the maximal kill (-4.01): E equals Emax at Ce = 0 and approaches E0 at
    # saturating exposure. The canonical roles used here are
    #   e0   <- Zhao 2018 Emax = +2.65 (effect at zero exposure)
    #   emax <- Zhao 2018 E0   = -4.01 (maximal drug effect)
    # and model() writes the identical E = e0 - (e0 - emax) * f.
    # Solving the packaged form at E = 0, -3 and -4 returns 22.74, 34.54 and
    # 87.98 h against the published 22.73, 34.52 and 88.03 h (Table 6).
    # e0 and emax are signed, so they stay on the natural scale.
    e0 <- 2.65
    label("Change in H. parasuis count over 24 h in drug-free serum (log10 CFU/mL)") # Zhao 2018 Table 6, 1 mg/kg, 'Emax' = 2.65 log10 CFU/mL
    emax <- -4.01
    label("Maximal change in H. parasuis count over 24 h (log10 CFU/mL; negative = reduction)") # Zhao 2018 Table 6, 1 mg/kg, 'E0' = -4.01 log10 CFU/mL
    lec50 <- log(24.66)
    label("AUC0-24h/MIC producing half the maximal effect EC50 (h)") # Zhao 2018 Table 6, 1 mg/kg, EC50 = 24.66
    lhill <- log(5.11)
    label("Hill coefficient N, steepness of the AUC0-24h/MIC effect curve (unitless)") # Zhao 2018 Table 6, 1 mg/kg, Slope (N) = 5.11

    # Measured property of the challenge strain, not an estimate. The serum
    # MIC, not the CAMHB MIC of 0.5 ug/mL, is the index denominator: Table 5
    # AUC168h/MIC serum = 367.3 h equals Table 2 AUC0-168h 11.0196 h*ug/mL
    # divided by 0.03 ug/mL.
    mic <- fixed(0.03)
    label("Tulathromycin MIC against H. parasuis 13R in guinea pig serum (ug/mL)") # Zhao 2018 Results (MICs, MBCs, MPC and PAE): 13R MIC 0.03 ug/mL in guinea pig serum

    # Experimental design input, not an estimate. The 24 h change E does not
    # depend on it.
    log10_cfu0 <- fixed(5.7)
    label("log10 starting bacterial density in the ex vivo incubation (log10 CFU/mL)") # Read by the maintainers from Zhao 2018 Figure 5 (all curves start at about 5.7 log10 CFU/mL)

    # No residual SD is reported for the PK/PD-integration fit.
    addSd <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (0; not reported in Zhao 2018)") # Zhao 2018 Table 6 reports R^2 = 0.9961 only
  })

  model({
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # PK/PD index AUC0-24h/MIC (h)
    aucmic <- AUC_TULATHROMYCIN / mic

    # Zhao 2018 sigmoid Emax equation in canonical roles; signed 24 h change
    # in log10 CFU/mL.
    effect <- e0 - (e0 - emax) * aucmic^hill / (aucmic^hill + ec50^hill)

    # Spread the 24 h change uniformly so log10(bact) moves by exactly
    # `effect` over each 24 h window:
    #   d(log10 N)/dt = effect / 24  =>  dN/dt = ln(10) * (effect / 24) * N
    d/dt(bact) <- log(10) * (effect / 24) * bact
    bact(0) <- 10^log10_cfu0

    # log10 CFU/mL with a 1 CFU/mL floor so the log stays finite.
    log_cfu <- log10(bact + 1)
    log_cfu ~ add(addSd)
  })
}
