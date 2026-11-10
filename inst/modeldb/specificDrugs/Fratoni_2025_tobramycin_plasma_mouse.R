Fratoni_2025_tobramycin_plasma_mouse <- function() {
  description <- "Preclinical (mouse, CD-1 female, neutropenic and renally impaired). One-compartment PK model of total tobramycin in PLASMA after subcutaneous injection, with first-order absorption and first-order elimination, in the standardized COMBINE murine neutropenic Klebsiella pneumoniae pneumonia model. Naive-pooled fit (WinNonlin nonlinear least squares) used to design and confirm a human-simulated regimen (HSR) reproducing the plasma fCmax and fAUC0-24 of tobramycin 7 mg/kg as a 30 min infusion. The volume is apparent (V/F) because the subcutaneous bioavailability was not estimated. Free plasma concentration Cu = fu * Cc with the murine unbound fraction fu = 0.776, the mean of ex vivo ultrafiltration binding measured in the same study (22.4%). Companion to the ELF model Fratoni_2025_tobramycin_elf_mouse, which was fitted separately to the epithelial lining fluid concentrations."
  reference <- "Fratoni AJ, Padgett AM, Roenfanz HF, Duffy EM, Nicolau DP. Development and confirmation of humanized plasma and epithelial lining fluid exposures of meropenem, cefiderocol and tobramycin in a standardized neutropenic murine pneumonia model. J Antimicrob Chemother. 2025;80(2):478-484. doi:10.1093/jac/dkae432. PMID: 39657022. PMCID: PMC11787888. Parameter estimates: Table 2. Model structure and fitting method: Materials and methods, 'Human-simulated exposure pharmacokinetic studies', and Results ('one-compartment pharmacokinetic models with first order elimination'). Murine tobramycin protein binding 22.4% (77.6% unbound): Results, 'Ex vivo tobramycin serum protein binding studies'. Dosing regimens: Table 4."
  vignette <- "Fratoni_2025_humanized_exposures_mouse"
  units <- list(time = "h", dosing = "mg/kg", concentration = "mg/L")

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "tobramycin", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tobramycin", units = "mg/kg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "mouse (CD-1, female, specific pathogen-free, 6-8 weeks old)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    sex_female_pct = 100,
    disease_state = "Neutropenic (cyclophosphamide 150 mg/kg IP on day -4 and 100 mg/kg on day -1) and renally impaired (uranyl nitrate 5 mg/kg IP on day -3) mice with Klebsiella pneumoniae CDC 851 lung infection (intranasal inoculation, 50 uL), per the COMBINE standardized murine pneumonia protocol",
    dose_range = "Subcutaneous tobramycin plasma HSR over 24 h: 9 mg/kg at 0 h, 4 mg/kg at 3 h, 2 mg/kg at 6 h and 0.5 mg/kg at 9 h (Table 4); protein-binding sub-study: single subcutaneous doses of 2.5, 5 and 10 mg/kg",
    regions = "Hartford Hospital, Hartford, CT, USA",
    notes = paste(
      "Destructive, sparse sampling: groups of six mice per time point, two retro-orbital bleeds plus one terminal cardiac puncture per mouse, with bronchoalveolar lavage at the terminal time point.",
      "Across the three drugs, 10 of 108 plasma samples and 7 of 78 BAL samples were excluded as interquartile-range outliers; the number of mice per drug is not reported.",
      "Protein binding: triplicate pooled plasma and ultrafiltrate at 1 h after 2.5, 5 and 10 mg/kg (15 mice per dose) gave 27.1%, 20.5% and 19.8% bound, treated as exposure-independent and averaged to 22.4%.",
      "The model was fitted to TOTAL plasma concentrations: the free concentrations the paper reports (fCmax 17.7 mg/L, fAUC0-24 103 mg.h/L, and the 'Free Tobramycin Plasma Conc' axis of Figure 3) are fu * Cc.",
      "No between-animal variability or residual error model is reported (naive-pooled least-squares fit)."
    )
  )

  ini({
    lka <- log(23.68); label("Apparent first-order subcutaneous absorption rate constant Ka (1/h)") # Table 2, Tobramycin: Ka = 23.68 1/h
    lkel <- log(0.31); label("First-order elimination rate constant from the central compartment Kel (1/h)") # Table 2, Tobramycin: Kel = 0.31 1/h
    lvc <- log(0.37); label("Apparent plasma volume of distribution V/F (L/kg)") # Table 2, Tobramycin: V plasma = 0.37 L/kg
    fu <- fixed(0.776); label("Fraction of tobramycin unbound in mouse plasma (unitless)") # Results, ex vivo protein binding: mean 22.4% bound, '77.6% unbound'

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
