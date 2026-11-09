Fratoni_2025_meropenem_plasma_mouse <- function() {
  description <- "Preclinical (mouse, CD-1 female, neutropenic and renally impaired). One-compartment PK model of total meropenem in PLASMA after subcutaneous injection, with first-order absorption and first-order elimination, in the standardized COMBINE murine neutropenic Klebsiella pneumoniae pneumonia model. Naive-pooled fit (WinNonlin nonlinear least squares) used to design and confirm a human-simulated regimen (HSR) reproducing the plasma %fT > MIC of meropenem 2 g every 8 h as a 3 h infusion. The volume is apparent (V/F) because the subcutaneous bioavailability was not estimated. The fitted rate constants are flip-flop (Ka = 1.34 1/h much smaller than Kel = 27.13 1/h), so the terminal decline is absorption-limited. Free plasma concentration Cu = fu * Cc with the murine unbound fraction fu = 0.92 (8% binding). Companion to the ELF model Fratoni_2025_meropenem_elf_mouse, which was fitted separately to the epithelial lining fluid concentrations."
  reference <- "Fratoni AJ, Padgett AM, Roenfanz HF, Duffy EM, Nicolau DP. Development and confirmation of humanized plasma and epithelial lining fluid exposures of meropenem, cefiderocol and tobramycin in a standardized neutropenic murine pneumonia model. J Antimicrob Chemother. 2025;80(2):478-484. doi:10.1093/jac/dkae432. PMID: 39657022. PMCID: PMC11787888. Parameter estimates: Table 2. Model structure and fitting method: Materials and methods, 'Human-simulated exposure pharmacokinetic studies', and Results ('one-compartment pharmacokinetic models with first order elimination'). Murine meropenem protein binding 8%: Materials and methods (cited to DeRyke CA et al. Antimicrob Agents Chemother. 2007;51:1481-6). Dosing regimens: Table 4."
  vignette <- "Fratoni_2025_humanized_exposures_mouse"
  units <- list(time = "h", dosing = "mg/kg", concentration = "mg/L")

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "meropenem", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "meropenem", units = "mg/kg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "mouse (CD-1, female, specific pathogen-free, 6-8 weeks old)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    sex_female_pct = 100,
    disease_state = "Neutropenic (cyclophosphamide 150 mg/kg IP on day -4 and 100 mg/kg on day -1) and renally impaired (uranyl nitrate 5 mg/kg IP on day -3) mice with Klebsiella pneumoniae CDC 851 lung infection (intranasal inoculation, 50 uL), per the COMBINE standardized murine pneumonia protocol",
    dose_range = "Subcutaneous meropenem plasma HSR: 65 mg/kg at 0 h, 65 mg/kg at 1.25 h, 45 mg/kg at 3.5 h and 45 mg/kg at 6 h, repeated every 8 h (Table 4)",
    regions = "Hartford Hospital, Hartford, CT, USA",
    notes = paste(
      "Destructive, sparse sampling: groups of six mice per time point, two retro-orbital bleeds plus one terminal cardiac puncture per mouse, with bronchoalveolar lavage at the terminal time point.",
      "Across the three drugs, 10 of 108 plasma samples and 7 of 78 BAL samples were excluded as interquartile-range outliers; the number of mice per drug is not reported.",
      "The model was fitted to TOTAL plasma concentrations: the free concentrations the paper reports (fCmax, fAUC, %fT > MIC, and the 'Free Meropenem Plasma Conc' axis of Figure 1) are fu * Cc.",
      "No between-animal variability or residual error model is reported (naive-pooled least-squares fit)."
    )
  )

  ini({
    lka <- log(1.34); label("Apparent first-order subcutaneous absorption rate constant Ka (1/h)") # Table 2, Meropenem: Ka = 1.34 1/h
    lkel <- log(27.13); label("First-order elimination rate constant from the central compartment Kel (1/h)") # Table 2, Meropenem: Kel = 27.13 1/h
    lvc <- log(0.04); label("Apparent plasma volume of distribution V/F (L/kg)") # Table 2, Meropenem: V plasma = 0.04 L/kg
    fu <- fixed(0.92); label("Fraction of meropenem unbound in mouse plasma (unitless)") # Materials and methods: murine protein binding 8% (ref 18), so fu = 1 - 0.08

    propSd <- fixed(0); label("Proportional residual SD (fraction); not reported for the least-squares fit") # Materials and methods: WinNonlin nonlinear least squares; weighting and residual model not reported
  })

  model({
    ka <- exp(lka)
    kel <- exp(lkel)
    vc <- exp(lvc)

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    Cc <- central / vc
    Cu <- fu * Cc
    Cc ~ prop(propSd)
  })
}
