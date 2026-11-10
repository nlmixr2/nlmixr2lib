Fratoni_2025_cefiderocol_elf_mouse <- function() {
  description <- "Preclinical (mouse, CD-1 female, neutropenic and renally impaired). One-compartment PK model of cefiderocol in pulmonary EPITHELIAL LINING FLUID (ELF) after subcutaneous injection, with first-order absorption from the injection site directly into the ELF compartment and first-order elimination from it, in the standardized COMBINE murine neutropenic Klebsiella pneumoniae pneumonia model. The ELF concentrations (urea-corrected bronchoalveolar lavage, treated as unbound) were fitted on their own by naive-pooled WinNonlin nonlinear least squares, with no link to plasma, and the fit was used to design and confirm a human-simulated regimen (HSR) reproducing the ELF %fT > MIC of cefiderocol 2 g every 8 h as a 3 h infusion (the human ELF profile is characterized only at 3 h and 5 h). The volume is apparent (V/F) because the subcutaneous bioavailability was not estimated. Companion to the plasma model Fratoni_2025_cefiderocol_plasma_mouse."
  reference <- "Fratoni AJ, Padgett AM, Roenfanz HF, Duffy EM, Nicolau DP. Development and confirmation of humanized plasma and epithelial lining fluid exposures of cefiderocol, cefiderocol and tobramycin in a standardized neutropenic murine pneumonia model. J Antimicrob Chemother. 2025;80(2):478-484. doi:10.1093/jac/dkae432. PMID: 39657022. PMCID: PMC11787888. Parameter estimates: Table 3. Model structure and fitting method: Materials and methods, 'Human-simulated exposure pharmacokinetic studies', and Results ('one-compartment pharmacokinetic models with first order elimination'). ELF by urea correction, unbound: Materials and methods. Dosing regimens: Table 4."
  vignette <- "Fratoni_2025_humanized_exposures_mouse"
  units <- list(time = "h", dosing = "mg/kg", concentration = "mg/L")

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "cefiderocol", units = "mg/kg", specimen = "administration site", verified = TRUE),
    elf = list(analyte = "cefiderocol", units = "mg/kg", specimen = "epithelial lining fluid", verified = TRUE)
  )

  population <- list(
    species = "mouse (CD-1, female, specific pathogen-free, 6-8 weeks old)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    sex_female_pct = 100,
    disease_state = "Neutropenic (cyclophosphamide 150 mg/kg IP on day -4 and 100 mg/kg on day -1) and renally impaired (uranyl nitrate 5 mg/kg IP on day -3) mice with Klebsiella pneumoniae CDC 851 lung infection (intranasal inoculation, 50 uL), per the COMBINE standardized murine pneumonia protocol",
    dose_range = "Subcutaneous cefiderocol ELF HSR: 3.75 mg/kg at 0 h, 5 mg/kg at 1 h, 6.25 mg/kg at 2 h, 1.75 mg/kg at 4 h and 1 mg/kg at 6 h, repeated every 8 h (Table 4)",
    regions = "Hartford Hospital, Hartford, CT, USA",
    notes = paste(
      "Destructive, sparse sampling: groups of six mice per time point, with bronchoalveolar lavage (4 x 0.4 mL saline) at the terminal time point after cardiac puncture.",
      "ELF concentration = BAL concentration * (plasma urea / BAL urea); because BAL contains no albumin, total ELF concentrations were considered unbound.",
      "Across the three drugs, 10 of 108 plasma samples and 7 of 78 BAL samples were excluded as interquartile-range outliers; the number of mice per drug is not reported.",
      "No between-animal variability or residual error model is reported (naive-pooled least-squares fit)."
    )
  )

  ini({
    lka <- log(1.94); label("Apparent first-order absorption rate constant from the subcutaneous site into ELF, Ka-ELF (1/h)") # Table 3, Cefiderocol: Ka-ELF = 1.94 1/h
    lkeff_elf <- log(0.72); label("First-order elimination rate constant from the ELF compartment, Kel-ELF (1/h)") # Table 3, Cefiderocol: Kel-ELF = 0.72 1/h
    lv_elf <- log(0.33); label("Apparent ELF volume of distribution V/F (L/kg)") # Table 3, Cefiderocol: V ELF = 0.33 L/kg

    propSd <- fixed(0); label("Proportional residual SD (fraction); not reported for the least-squares fit") # Materials and methods: WinNonlin nonlinear least squares; weighting and residual model not reported
  })

  model({
    ka <- exp(lka)
    keff_elf <- exp(lkeff_elf)
    v_elf <- exp(lv_elf)

    d/dt(depot) <- -ka * depot
    d/dt(elf) <- ka * depot - keff_elf * elf

    Celf <- elf / v_elf
    Celf ~ prop(propSd)
  })
}
