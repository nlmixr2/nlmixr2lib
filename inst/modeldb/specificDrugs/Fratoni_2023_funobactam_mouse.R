Fratoni_2023_funobactam_mouse <- function() {
  description <- paste(
    "Preclinical (mouse). One-compartment pharmacokinetic model for the",
    "beta-lactamase inhibitor funobactam (formerly XNW4107) in the neutropenic",
    "murine thigh infection model, with first-order absorption from a",
    "subcutaneous depot and first-order elimination (Fratoni 2023). Funobactam",
    "was given subcutaneously at 1, 10 and 20 mg/kg together with an imipenem",
    "human-simulated regimen in uranyl-nitrate-treated mice (renal impairment",
    "induced to slow clearance), and a single best-fit average parameter set",
    "was estimated in Phoenix WinNonlin from destructive-sampling mean",
    "profiles. Cc is the total plasma concentration that was assayed; Cu is",
    "the unbound concentration (mean mouse plasma protein binding 3.53%) that",
    "the authors used for every %fT>MIC, %fT>CT, fAUC and fCmax exposure",
    "index. The apparent volume is 0.462 L/kg, one tenth of the 4.62 L/kg",
    "printed in the Results: only 0.462 L/kg reproduces the authors' own",
    "simulated exposures in Table S1 and Table 2 (see the vignette). No",
    "between-subject variability or residual error was reported, so the model",
    "has no etas and the residual SD is fixed(0).",
    sep = " "
  )
  reference <- paste(
    "Fratoni AJ, Berry AV, Liu X, Chen X, Wu Y, Nicolau DP, Abdelraouf K.",
    "Imipenem/funobactam (formerly XNW4107) in vivo pharmacodynamics against",
    "serine carbapenemase-producing Gram-negative bacteria: a novel modelling",
    "approach for time-dependent killing. J Antimicrob Chemother.",
    "2023;78(9):2343-2353. doi:10.1093/jac/dkad242.",
    "PK parameters from Results ('PK studies'); protein binding from Results",
    "('Ex vivo plasma protein-binding'); simulated exposures used to pin the",
    "volume from Supplementary Table S1 and Table 2.",
    sep = " "
  )
  vignette <- "Fratoni_2023_funobactam_mouse"
  units <- list(time = "h", dosing = "mg/kg", concentration = "mg/L")

  # Amounts are per kilogram of body weight because the murine doses and the
  # fitted volume are both weight-normalised (mg/kg and L/kg), so the
  # amount / volume quotient is mg/L.
  compartmentData <- list(
    depot = list(
      analyte = "funobactam",
      units = "mg/kg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "funobactam",
      units = "mg/kg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "mouse (ICR, female, neutropenic thigh infection model)",
    n_subjects = 108,
    n_studies = 1,
    weight_range = "20-22 g",
    sex_female_pct = 100,
    disease_state = paste(
      "Neutropenic murine thigh infection. Mice were rendered neutropenic with",
      "cyclophosphamide (150 mg/kg day -4, 100 mg/kg day -1), given uranyl",
      "nitrate 5 mg/kg on day -3 to reduce renal clearance and generate",
      "human-simulated exposures, and inoculated intramuscularly in one thigh",
      "with 0.1 mL of ~10^7 cfu/mL bacterial suspension 2 h before dosing.",
      sep = " "
    ),
    dose_range = paste(
      "funobactam 1, 10 and 20 mg/kg subcutaneously (single dose) with the",
      "imipenem human-simulated regimen (8, 2.5, 1.5 and 0.9 mg/kg at 0, 1.25,",
      "2.5 and 4 h, repeated every 6 h)",
      sep = " "
    ),
    notes = paste(
      "Methods, 'Funobactam single-dose studies': 36 infected mice per",
      "funobactam dose level (108 total), sacrificed in groups of six at",
      "predefined time points for cardiac-puncture plasma; total funobactam",
      "assayed by LC-MS/MS. Dose proportionality was confirmed over 1-20 mg/kg",
      "(Figure S2). No covariates were assessed.",
      sep = " "
    )
  )

  ini({
    # Results, 'PK studies': "PK was best described using a one-compartment
    # model. The best-fit average PK parameters were: apparent volume of
    # distribution, 4.62 L/kg; first-order absorption rate constant, 7.21 h-1;
    # and overall elimination rate constant, 0.41 h-1."
    lka <- log(7.21)
    label("First-order absorption rate constant from the subcutaneous depot (1/h)")
    # Results, 'PK studies': first-order absorption rate constant 7.21 h-1

    lvc <- log(0.462)
    label("Apparent volume of distribution V/F (L/kg)")
    # Results, 'PK studies' prints 4.62 L/kg. Encoded as 0.462 L/kg, a
    # decimal-shift correction: Table S1 (exposures the authors simulated
    # 'using parameters from single dose PK studies') and Table 2 rows 'b'
    # (fAUC0-24 35.68 and 71.36, fCmax 3.48 and 6.95) are reproduced to
    # rounding with 0.462 L/kg and are 10-fold too high for 4.62 L/kg
    # (e.g. 3.6 mg/kg q24h fAUC0-24 18.3 mg.h/L).

    lkel <- log(0.41)
    label("First-order elimination rate constant (1/h)")
    # Results, 'PK studies': overall elimination rate constant 0.41 h-1

    fu <- fixed(0.9647)
    label("Fraction of funobactam unbound in mouse plasma (unitless)")
    # Results, 'Ex vivo plasma protein-binding': overall average binding
    # 3.53% (1, 10 and 20 mg/kg: 3.11, 5.44 and 2.03%); 1 - 0.0353 = 0.9647

    # Residual error was not reported: the single-dose profiles were fitted
    # as naive-pooled means in Phoenix WinNonlin and the paper reports no
    # error model. Fixed at zero; see the vignette Errata.
    propSd <- fixed(0)
    label("Proportional residual error on total plasma concentration (fraction; ZERO - not reported in source)")
  })

  model({
    # 1. Individual parameters. Destructive sampling gives one mean profile
    #    per dose, so no between-subject variability is declared.
    ka <- exp(lka)
    vc <- exp(lvc)
    kel <- exp(lkel)

    # 2. ODE system. The subcutaneous dose enters `depot`; bioavailability is
    #    not identifiable from subcutaneous data and is absorbed into V/F.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # 3. Observation. The assay measured TOTAL plasma funobactam; Cu is the
    #    unbound concentration behind every free-drug index in the paper.
    Cc <- central / vc
    Cu <- fu * Cc
    Cc ~ prop(propSd)
  })
}
