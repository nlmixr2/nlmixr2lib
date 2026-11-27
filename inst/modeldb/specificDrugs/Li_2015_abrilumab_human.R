Li_2015_abrilumab_human <- function() {
  description <- paste(
    "Human projection scaled from cynomolgus monkey (no human data were",
    "fitted). Two-compartment target-mediated drug disposition (TMDD) model",
    "with quasi-equilibrium (QE) binding and constant total receptor for",
    "abrilumab (AMG 181, MEDI7183), a human anti-alpha4beta7 integrin IgG2",
    "monoclonal antibody, used to select first-in-human doses for",
    "inflammatory bowel disease: CL and Vc scaled from a 2.8 kg monkey to a",
    "70 kg human with exponents 0.75 and 1, K12 and K21 with -0.25, Kd set to",
    "the human in vitro value, and KA, F1, Rtot and Kint unscaled (Li 2015)."
  )
  reference <- paste(
    "Li H, Kock K, Wisler JA, Rees WA, Prince PJ, Reynhardt KO, Hsu H, Yu Z,",
    "Borie DC, Salinger DH, Pan WJ. Prediction of clinical pharmacokinetics",
    "of AMG 181, a human anti-alpha4beta7 monoclonal antibody for treating",
    "inflammatory bowel diseases. Pharmacol Res Perspect. 2015;3(1):e00098.",
    "doi:10.1002/prp2.98"
  )
  vignette <- "Li_2015_abrilumab"
  dosing <- c("depot", "central")
  units <- list(time = "day", dosing = "mg", concentration = "nmol/L")

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "abrilumab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "abrilumab (free + alpha4beta7-bound)", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "abrilumab", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 0L,
    weight_median = "70 kg (assumed for the allometric scaling)",
    disease_state = "healthy volunteers and patients with inflammatory bowel disease (projected; the paper fits no human data)",
    dose_range = "0.7-210 mg SC and 70-420 mg IV single doses; 7-210 mg SC Q4W x 3 (Table 4)",
    notes = paste(
      "Projection only. Typical-value parameters of Table 3, allometrically",
      "scaled from the cynomolgus monkey model (Li_2015_abrilumab_cyno; median",
      "monkey weight 2.8 kg). The projection was compared with the",
      "first-in-human study (7 mg SC to 420 mg IV; Figure 7), whose data were",
      "reported separately (Pan et al. 2014) and were not used for fitting."
    )
  )

  ini({
    # Li 2015 Table 3 'Human estimate'. These are projections, not
    # estimates, so every value is fixed.
    lka <- fixed(log(0.846)); label("First-order SC absorption rate constant KA (1/day)") # Table 3: KA = 0.846 /day (not scaled)
    lcl <- fixed(log(144)); label("Linear clearance CL (mL/day)") # Table 3: CL = 144 mL/day (12.9 x (70/2.8)^0.75)
    lvc <- fixed(log(2900)); label("Central volume of distribution Vc (mL)") # Table 3: Vc = 2900 mL (116 x 70/2.8)
    lk12 <- fixed(log(0.276)); label("Central-to-peripheral rate constant K12 (1/day)") # Table 3: K12 = 0.276 /day (0.617 x (70/2.8)^-0.25)
    lk21 <- fixed(log(0.392)); label("Peripheral-to-central rate constant K21 (1/day)") # Table 3: K21 = 0.392 /day (0.877 x (70/2.8)^-0.25)
    lkd <- fixed(log(0.031)); label("Equilibrium dissociation constant Kd (nmol/L)") # Methods 'Allometric scaling and simulation of human PK': Kd fixed at 0.031 nmol/L (in vitro EC50, human primary T cells); Table 3 prints 0.013, see vignette
    lrtot <- fixed(log(9.18)); label("Total alpha4beta7 receptor concentration Rtot (nmol/L)") # Table 3: Rtot = 9.18 nmol/L (not scaled)
    lkint <- fixed(log(0.0215)); label("Internalization rate constant of the drug-receptor complex Kint (1/day)") # Table 3: Kint = 0.0215 /day (not scaled)
    lfdepot <- fixed(log(1)); label("SC bioavailability F1 (fraction)") # Table 3: F1 = 1.0 (not scaled)

    # The human projection carries no variability (typical profiles only).
    propSd <- fixed(0); label("Proportional residual error (fraction; none defined for the projection)")
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    kd <- exp(lkd)
    rtot <- exp(lrtot)
    kint <- exp(lkint)
    fdepot <- exp(lfdepot)

    # Molecular weight about 144 kDa (Methods, 'In vitro pharmacology
    # studies'). Doses and drug amounts are in mg and volumes in mL, so
    # mg/mL x 1e6 / 144 = nmol/L, and nmol/L x mL x 144 / 1e6 = mg.
    mw <- 144
    nm_per_mgml <- 1e6 / mw
    mg_per_nmml <- mw / 1e6

    # QE approximation with constant total receptor (Figure 1): the central
    # state holds total (free + bound) drug; free drug Cf solves
    # Ctot = Cf + Rtot * Cf / (Kd + Cf).
    ctot <- central / vc * nm_per_mgml
    disc <- ctot - rtot - kd
    cfree <- 0.5 * (disc + sqrt(disc * disc + 4 * kd * ctot))
    cbound <- ctot - cfree

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - ((cl / vc + k12) * cfree + kint * cbound) * vc * mg_per_nmml + k21 * peripheral1
    d/dt(peripheral1) <- k12 * cfree * vc * mg_per_nmml - k21 * peripheral1

    f(depot) <- fdepot

    # Free concentration, as in the monkey model; ug/mL = nmol/L x 0.144.
    Cc <- cfree
    Cc ~ prop(propSd)
  })
}
