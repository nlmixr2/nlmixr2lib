Abouelhassan_2024_sulbactam_mouse <- function() {
  description <- paste(
    "Preclinical (mouse). One-compartment epithelial lining fluid (ELF)",
    "pharmacokinetic model for sulbactam in the neutropenic murine",
    "Acinetobacter baumannii pneumonia model, with first-order absorption from",
    "a subcutaneous depot and first-order elimination (Abouelhassan 2024).",
    "Sulbactam was given subcutaneously as commercial ampicillin-sulbactam at",
    "1, 10, 25, 100 and 200 mg/kg; ELF concentrations were derived from",
    "bronchoalveolar lavage using the plasma-to-BAL urea ratio, and each dose",
    "level was fitted separately by nonlinear least squares in Phoenix",
    "WinNonlin. The ini() values are the arithmetic means of the four",
    "dose-proportional fits spanning 1-100 mg/kg (Table 2), which is the",
    "parameter set the authors used to simulate sulbactam ELF exposures for",
    "every dose-ranging arm at or below 100 mg/kg. ELF exposure is",
    "dose-proportional across that range (R^2 = 0.997); the 200 mg/kg fit",
    "deviates from proportionality and is recorded in the ini() source-trace",
    "comments rather than carried as a second model. Sulbactam penetrated the",
    "ELF with a mean ELF-to-free-plasma AUC(0-8) ratio of 0.66 (SD 0.05), and",
    "the model reproduces the reported ELF AUC(0-8) at every dose level. No",
    "between-subject variability or residual error was reported because each",
    "dose level was fitted as a single naive-pooled destructive-sampling",
    "profile, so both are encoded as fixed(0).",
    sep = " "
  )
  reference <- paste(
    "Abouelhassan Y, Kuti JL, Nicolau DP, Abdelraouf K.",
    "Pharmacokinetic/pharmacodynamic analysis of sulbactam against",
    "Acinetobacter baumannii pneumonia: establishing in vivo efficacy targets",
    "in the epithelial lining fluid. JAC Antimicrob Resist. 2024;6(6):dlae203.",
    "doi:10.1093/jacamr/dlae203.",
    "Murine ELF PK parameters from Table 2; study design from Methods",
    "('Sulbactam bronchopulmonary PK studies').",
    sep = " "
  )
  vignette <- "Abouelhassan_2024_sulbactam_elf_pneumonia"
  units    <- list(time = "h", dosing = "mg/kg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Amounts are per kilogram of body weight because the
  # murine doses and the fitted volume are both weight-normalised (mg/kg and
  # L/kg), so the amount / volume quotient is still mg/L.
  compartmentData <- list(
    depot = list(
      analyte = "sulbactam", units = "mg/kg",
      specimen = "administration site", verified = TRUE
    ),
    elf = list(
      analyte = "sulbactam", units = "mg/kg",
      specimen = "epithelial lining fluid", verified = TRUE
    )
  )

  population <- list(
    species        = "mouse (ICR/CD-1, female, neutropenic Acinetobacter baumannii pneumonia model)",
    n_subjects     = 180,
    n_studies      = 1,
    weight_range   = "20-22 g",
    disease_state  = paste(
      "Neutropenic murine Acinetobacter baumannii pneumonia. Mice were rendered",
      "neutropenic with cyclophosphamide and uranyl nitrate, then inoculated",
      "intranasally with 10^7 cfu/mL A. baumannii in 3% mucin 2 h before the",
      "first antimicrobial dose.",
      sep = " "
    ),
    dose_range     = "sulbactam 1, 10, 25, 100 and 200 mg/kg subcutaneously (single dose, dosed as ampicillin-sulbactam)",
    notes          = paste(
      "Methods, 'Murine neutropenic pneumonia model' and 'Sulbactam",
      "bronchopulmonary PK studies': groups of 36 mice per dose level (5 dose",
      "levels = 180 mice), euthanised in groups of 6 at 4-6 time points per",
      "dose for destructive blood and bronchoalveolar lavage sampling. ELF",
      "concentrations = BAL sulbactam concentration x (plasma urea / BAL",
      "urea). No covariates were assessed.",
      sep = " "
    )
  )

  ini({
    # Structural parameters -- arithmetic means of the four dose-proportional
    # single-dose fits in Table 2 (1, 10, 25 and 100 mg/kg). Methods,
    # 'Pharmacokinetics/pharmacodynamics analyses': "The PK parameters for the
    # four doses of 1, 10, 25 and 100 mg/kg across the linear range were
    # averaged to determine sulbactam exposures for doses <= 100 mg/kg in the
    # dose-ranging studies". The 200 mg/kg column is deliberately excluded from
    # each mean because the Results state that dose deviates from
    # proportionality; its values are recorded on each line for reference.
    lka   <- log(17.605)
    label("Absorption rate constant Ka from the subcutaneous depot (1/h)")
    # Table 2 Ka row: 8.28 (1), 26.71 (10), 26.79 (25), 8.64 (100) mg/kg
    #   -> mean = (8.28 + 26.71 + 26.79 + 8.64) / 4 = 17.605 1/h
    #   [200 mg/kg fit, excluded from the mean: Ka = 78.44 1/h]

    lvelf <- log(0.5125)
    label("Apparent ELF distribution volume Vc/F (L/kg)")
    # Table 2 Vc/F row: 0.59 (1), 0.52 (10), 0.45 (25), 0.49 (100) mg/kg
    #   -> mean = (0.59 + 0.52 + 0.45 + 0.49) / 4 = 0.5125 L/kg
    #   [200 mg/kg fit, excluded from the mean: Vc/F = 0.57 L/kg]

    lkel  <- log(0.575)
    label("First-order elimination rate constant Kel from ELF (1/h)")
    # Table 2 Kel row: 0.62 (1), 0.48 (10), 0.68 (25), 0.52 (100) mg/kg
    #   -> mean = (0.62 + 0.48 + 0.68 + 0.52) / 4 = 0.575 1/h
    #   [200 mg/kg fit, excluded from the mean: Kel = 0.24 1/h]

    # Residual error was not reported: each dose level was fitted as a single
    # naive-pooled destructive-sampling profile in Phoenix WinNonlin and the
    # paper reports no error model. Fixed at zero per the standing policy for
    # unreported RUV; see the vignette Errata.
    propSd_Celf <- fixed(0)
    label("Proportional residual error on ELF concentration (fraction; ZERO - not reported in source)")
    # Methods, 'Sulbactam bronchopulmonary PK studies': model selection was by
    # visual inspection of the fit and model diagnostics; no error model given.
  })

  model({
    # 1. Individual parameters. No between-subject variability is estimable
    #    from destructive sampling, so no etas are declared.
    ka   <- exp(lka)
    velf <- exp(lvelf)
    kel  <- exp(lkel)

    # 2. ODE system. Subcutaneous dose enters `depot`; the single fitted
    #    disposition compartment is the epithelial lining fluid, because the
    #    paper fitted the ELF concentration-time profiles directly ("ELF
    #    concentrations were best described by a one-compartment model with
    #    first-order elimination"). Vc/F and Kel are therefore ELF-apparent.
    d/dt(depot) <- -ka * depot
    d/dt(elf)   <-  ka * depot - kel * elf

    # 3. Observation. Bioavailability is not separately identifiable from the
    #    subcutaneous data, so it is absorbed into the apparent volume Vc/F
    #    exactly as the paper reports it.
    Celf <- elf / velf
    Celf ~ prop(propSd_Celf)
  })
}
