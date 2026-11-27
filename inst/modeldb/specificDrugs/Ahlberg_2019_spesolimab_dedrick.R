Ahlberg_2019_spesolimab_dedrick <- function() {
  description <- paste(
    "Preclinical (cynomolgus monkey, Dedrick-scaled to human). Linear two-compartment PK",
    "model for spesolimab (MAB92, BI 655130; humanized anti-IL-36R IgG1) fitted to the mean",
    "10 mg/kg intravenous cynomolgus monkey profile after Dedrick transformation to human",
    "time and concentration scales (allometric exponents 1.0 for volume and 0.85 for",
    "clearance), so the parameters are on the human scale (Ahlberg 2019). Spesolimab does not",
    "bind monkey IL-36R, so the fit carries no target-mediated term. Deterministic (no IIV,",
    "no residual error reported)."
  )
  reference <- paste(
    "Ahlberg J, Giragossian C, Li H, Myzithras M, Raymond E, Caviness G, Grimaldi C,",
    "Brown SE, Perez R, Yang D, Kroe-Barrett R, Joseph D, Pamulapati C, Coble K, Ruus P,",
    "Woska JR, Ganesan R, Hansel S, Mbow ML.",
    "Retrospective analysis of model-based predictivity of human pharmacokinetics for",
    "anti-IL-36R monoclonal antibody MAB92 using a rat anti-mouse IL-36R monoclonal antibody",
    "and RNA expression data (FANTOM5). MAbs. 2019;11(5):956-964.",
    "doi:10.1080/19420862.2019.1615345"
  )
  vignette <- "Ahlberg_2019_spesolimab"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  covariateData <- list()

  compartmentData <- list(
    central = list(analyte = "spesolimab", units = "nmol", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "spesolimab", units = "nmol", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "cynomolgus monkey (Mauritius origin, male; Dedrick-scaled to human)",
    n_subjects = 3,
    n_studies = 1,
    weight_range = "2.3-2.7 kg (all 9 monkeys in the study)",
    sex_female_pct = 0,
    disease_state = "Healthy, drug-naive",
    dose_range = "Single 10 mg/kg intravenous dose as a 10-min infusion (the 0.3 and 1.5 mg/kg groups were not used in this fit)",
    notes = paste(
      "Materials and methods, Nonhuman primate studies (Charles River, Reno NV): 9 monkeys,",
      "3 per dose group at 0.3, 1.5 or 10 mg/kg, sampled to 1008 h. Figure 5 shows the",
      "Dedrick-transformed mean 10 mg/kg profile and its fit. ADA in 2 of 3 animals at",
      "10 mg/kg, with no visible PK impact (Results, Monkey PK)."
    )
  )

  ini({
    # Table 4 'Monkey mean' column. Rate constants are printed in 'day-1' but the fit ran
    # in hours (Figure 5 time axis in h; terminal half-life from these values, ~630 h,
    # matches the slope of Figure 5). Values are used as printed in 1/h.
    lvc <- log(3.0); label("Central volume of distribution, human scale (L)") # Table 4 Monkey mean Vc = 3.0 L
    lkel <- log(0.0021); label("Linear elimination rate constant (1/h)") # Table 4 Monkey mean kel = 0.0021
    lk12 <- log(0.0056); label("Central-to-peripheral rate constant (1/h)") # Table 4 Monkey mean k12 = 0.0056
    lk21 <- log(0.0074); label("Peripheral-to-central rate constant (1/h)") # Table 4 Monkey mean k21 = 0.0074
  })

  model({
    vc <- exp(lvc)
    kel <- exp(lkel)
    k12 <- exp(lk12)
    k21 <- exp(lk21)

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
  })
}
