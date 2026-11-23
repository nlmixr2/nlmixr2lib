Skip_2014_nickel <- function() {
  description <- paste(
    "Ecotoxicology (ground beetle, Pterostichus oblongopunctatus).",
    "One-compartment three-stage toxicokinetic model of nickel (Ni) in adult",
    "beetles fed Ni-contaminated food for 64 days and then uncontaminated",
    "food (Skip 2014 model 3E2, exact solution with separate assimilation and",
    "elimination rate constants for each stage). dCI/dt = kA * CE - kE * CI,",
    "where CI is the whole-body Ni concentration and CE the Ni concentration",
    "in food; kA and kE step at an estimated switch time ts within the uptake",
    "phase (stage I to II) and again, with CE, at the transfer to clean food",
    "tc (stage III; Eqs. 12-14). The food concentrations, the transfer time",
    "and the initial body concentration are exposure-scenario inputs, not",
    "estimates. Fitted by least squares to the geometric-mean concentration",
    "at each sampling time; no between-animal variability or residual error",
    "was reported, and the 95% confidence intervals of the stage-I",
    "parameters span many orders of magnitude (Table 3)."
  )
  reference <- paste(
    "Skip B, Bednarska AJ, Laskowski R. Toxicokinetics of metals in",
    "terrestrial invertebrates: making things straight with the",
    "one-compartment principle. PLoS One. 2014;9(9):e108740.",
    "doi:10.1371/journal.pone.0108740. Data: Laskowski R, Bednarska AJ,",
    "Spurgeon D, Svendsen C, van Gestel CAM. Three-phase metal kinetics in",
    "terrestrial invertebrates exposed to high metal concentrations. Sci",
    "Total Environ. 2010;408:3794-3802."
  )
  vignette <- "Skip_2014_metal_toxicokinetics"
  units <- list(
    time = "day",
    dosing = "none (exposure enters through the food-concentration inputs)",
    concentration = "mg/kg"
  )

  covariateData <- list()

  compartmentData <- list(
    central = list(
      analyte = "nickel",
      units = "mg/kg (whole-body dry weight concentration)",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "ground beetle (Pterostichus oblongopunctatus, adult)",
    n_subjects = 65,
    n_studies = 1,
    disease_state = "none (dietary nickel exposure experiment)",
    dose_range = paste(
      "Artificial food (ground mealworms and apple) spiked with Ni at a nominal",
      "2500 mg/kg dry weight (measured 2525 +/- 76 mg/kg) for 64 days, then",
      "uncontaminated food (measured 5.95 +/- 2.25 mg/kg) to day 96"
    ),
    regions = "Poland (laboratory study)",
    notes = paste(
      "Skip 2014 Case Studies, 'The one-compartment model with three stages',",
      "re-analysing the Laskowski 2010 data. Twenty-one beetles analysed at",
      "day 0 and four per sampling time at days 2, 4, 6, 8, 16, 32 and 64",
      "(uptake) and 66, 68, 80 and 96 (decontamination), 21 + 11 x 4 = 65",
      "beetles. Kept individually at 20 C on moist peat; guts voided 48 h",
      "before Ni analysis by graphite-furnace AAS. The model was fitted to",
      "the geometric-mean concentration at each time in MATLAB R2012b (fit(),",
      "nonlinear least squares, trust-region), cycling the start point and",
      "lower bound of ts over [0, tc]."
    )
  )

  ini({
    # Exposure-scenario inputs ('When estimating model parameters, the
    # following data were used', Case Studies). These are measured quantities
    # that Skip 2014 held fixed in the fit; override them to simulate another
    # exposure.
    cfood_uptake <- fixed(2525)
    label("Ni concentration in food during the uptake phase, CEu (mg/kg dry weight)") # Case Studies text, CEu = 2525 mg/kg
    cfood_depur <- fixed(5.95)
    label("Ni concentration in food during the decontamination phase, CEd (mg/kg dry weight)") # Case Studies text, CEd = 5.95 mg/kg
    tdepur <- fixed(64)
    label("Time of transfer to uncontaminated food, tc (day)") # Case Studies text, tc = 64 days
    lc0 <- fixed(log(0.71))
    label("Whole-body Ni concentration at the start of exposure, CI0 (mg/kg dry weight)") # Case Studies text, CI0 = 0.71 mg/kg

    # Model 3E2 (Table 3, last column).
    tkchange <- 1.978
    label("Switch time from stage I to stage II within the uptake phase, ts (day)") # Table 3, model 3E2: ts = 1.978
    lkassim <- log(0.0751)
    label("Assimilation rate constant, stage I (t <= ts), kA1 (1/day)") # Table 3, model 3E2: kA1 = 0.0751
    lkel <- log(0.6583)
    label("Elimination rate constant, stage I (t <= ts), kE1 (1/day)") # Table 3, model 3E2: kE1 = 0.6583
    lkassim_late <- log(0.0132)
    label("Assimilation rate constant, stage II (ts < t <= tc), kA2 (1/day)") # Table 3, model 3E2: kA2 = 0.0132
    lkel_late <- log(0.9497)
    label("Elimination rate constant, stage II (ts < t <= tc), kE2 (1/day)") # Table 3, model 3E2: kE2 = 0.9497
    lkassim_depur <- log(0.1916)
    label("Assimilation rate constant, stage III (t > tc), kA3 (1/day)") # Table 3, model 3E2: kA3 = 0.1916
    lkel_depur <- log(0.3311)
    label("Elimination rate constant, stage III (t > tc), kE3 (1/day)") # Table 3, model 3E2: kE3 = 0.3311

    addSd <- fixed(0)
    label("Additive residual SD (mg/kg); not reported, least-squares fit to geometric means") # not reported in Skip 2014
  })
  model({
    kassim <- exp(lkassim)
    kel <- exp(lkel)
    kassim_late <- exp(lkassim_late)
    kel_late <- exp(lkel_late)
    kassim_depur <- exp(lkassim_depur)
    kel_depur <- exp(lkel_depur)
    c0 <- exp(lc0)

    # Eq. 12 (t <= ts), Eq. 13 (ts < t <= tc) and Eq. 14 (t > tc). The food
    # concentration changes only at tc; the rate constants change at ts and
    # at tc.
    cfood <- cfood_uptake
    kassim_t <- kassim
    kel_t <- kel
    if (t > tkchange) {
      kassim_t <- kassim_late
      kel_t <- kel_late
    }
    if (t > tdepur) {
      cfood <- cfood_depur
      kassim_t <- kassim_depur
      kel_t <- kel_depur
    }

    # The state is the whole-body Ni concentration CI (mg/kg dry weight),
    # starting from the measured pre-exposure concentration CI0.
    central(0) <- c0
    d/dt(central) <- kassim_t * cfood - kel_t * central

    Cc <- central
    Cc ~ add(addSd)
  })
}
