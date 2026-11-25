Skip_2014_cadmium <- function() {
  description <- paste(
    "Ecotoxicology (ground beetle, Poecilus cupreus). One-compartment",
    "toxicokinetic model of cadmium (Cd) in adult beetles fed Cd-contaminated",
    "housefly larvae for 90 days and then uncontaminated larvae (Skip 2014,",
    "two-stage exact-solution model with separate assimilation and",
    "elimination rate constants for the uptake and decontamination phases).",
    "dCI/dt = kA * CE - kE * CI, where CI is the whole-body Cd concentration",
    "and CE the Cd concentration in food; kA, kE and CE step at the transfer",
    "to clean food (Eqs. 6-7). The food concentrations, the transfer time and",
    "the initial body concentration are exposure-scenario inputs, not",
    "estimates. Fitted by least squares to the mean concentration at each",
    "sampling time; no between-animal variability or residual error was",
    "reported."
  )
  reference <- paste(
    "Skip B, Bednarska AJ, Laskowski R. Toxicokinetics of metals in",
    "terrestrial invertebrates: making things straight with the",
    "one-compartment principle. PLoS One. 2014;9(9):e108740.",
    "doi:10.1371/journal.pone.0108740. Data: Kramarz P. Dynamics of",
    "accumulation and decontamination of cadmium and zinc in carnivorous",
    "invertebrates. 1. The ground beetle, Poecilus cupreus L. Bull Environ",
    "Contam Toxicol. 1999;63:531-537."
  )
  vignette <- "Skip_2014_metal_toxicokinetics"
  units <- list(
    time = "day",
    dosing = "none (exposure enters through the food-concentration inputs)",
    concentration = "mg/kg"
  )

  # There is no dose: the metal enters through the food-concentration
  # inputs cfood_uptake / cfood_depur, and central(0) is the measured
  # pre-exposure body concentration. The model has no dosing compartment.
  dosing <- character(0)

  covariateData <- list()

  compartmentData <- list(
    central = list(
      analyte = "cadmium",
      units = "mg/kg (whole-body dry weight concentration)",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "ground beetle (Poecilus cupreus, adult, laboratory culture)",
    n_subjects = 140,
    n_studies = 1,
    disease_state = "none (dietary cadmium exposure experiment)",
    dose_range = paste(
      "Housefly larvae reared on medium spiked with Cd at 50 mg/kg dry",
      "weight; measured Cd in the larvae 154.6 mg/kg (Figure 1 caption; 155",
      "+/- 2.2 mg/kg in the Case Studies text) for 90 days, then control larvae",
      "(0.92 mg/kg, taken from Maryanski et al.) to day 121"
    ),
    regions = "Poland (laboratory study)",
    notes = paste(
      "Skip 2014 Case Studies, 'The one-compartment model with two stages',",
      "re-analysing Kramarz 1999. Twenty beetles analysed at day 0; six",
      "randomly chosen beetles per sampling time at days 2, 6 and 13 and",
      "weekly to day 90, then at days 92, 96, 100, 107, 114 and 121 (20 + 20 x",
      "6 = 140 beetles). Beetles",
      "were starved 24 h to void the gut before Cd analysis by graphite-furnace",
      "AAS. The model was fitted to the mean concentration at each time in",
      "MATLAB R2012b (fit(), nonlinear least squares, trust-region)."
    )
  )

  ini({
    # Exposure-scenario inputs (Case Studies text and Figure 1 caption). These
    # are measured quantities that Skip 2014 held fixed in the fit; override
    # them to simulate another exposure.
    cfood_uptake <- fixed(154.6)
    label("Cd concentration in food during the uptake phase, CEu (mg/kg dry weight)") # Figure 1 caption (c, d), CEu = 154.6 mg/kg
    cfood_depur <- fixed(0.92)
    label("Cd concentration in food during the decontamination phase, CEd (mg/kg dry weight)") # Figure 1 caption (c, d) and Case Studies text, CEd = 0.92 mg/kg
    tdepur <- fixed(90)
    label("Time of transfer to uncontaminated food, tc (day)") # Case Studies text, 'fed contaminated larvae for 90 days'
    lc0 <- fixed(log(0.15))
    label("Whole-body Cd concentration at the start of exposure, CI0 (mg/kg dry weight)") # Figure 1 caption (a, b), CI0 = 0.15 mg/kg

    # Exact-solution model with kA1, kA2, kE1, kE2 (Table 2, last column).
    lkassim <- log(0.0018)
    label("Assimilation rate constant, uptake phase, kA1 (1/day)") # Table 2, exact solution model kA1, kA2, kE1, kE2: kA1 = 0.0018
    lkel <- log(0.0396)
    label("Elimination rate constant, uptake phase, kE1 (1/day)") # Table 2, exact solution model kA1, kA2, kE1, kE2: kE1 = 0.0396
    lkassim_depur <- log(0.1115)
    label("Assimilation rate constant, decontamination phase, kA2 (1/day)") # Table 2, exact solution model kA1, kA2, kE1, kE2: kA2 = 0.1115
    lkel_depur <- log(0.321)
    label("Elimination rate constant, decontamination phase, kE2 (1/day)") # Table 2, exact solution model kA1, kA2, kE1, kE2: kE2 = 0.321

    addSd <- fixed(0)
    label("Additive residual SD (mg/kg); not reported, least-squares fit to means") # not reported in Skip 2014
  })
  model({
    kassim <- exp(lkassim)
    kel <- exp(lkel)
    kassim_depur <- exp(lkassim_depur)
    kel_depur <- exp(lkel_depur)
    c0 <- exp(lc0)

    # Eq. 6 (t <= tc) and Eq. 7 (t > tc): the food concentration and both
    # rate constants step when the beetles move to uncontaminated food.
    cfood <- cfood_uptake
    kassim_t <- kassim
    kel_t <- kel
    if (t > tdepur) {
      cfood <- cfood_depur
      kassim_t <- kassim_depur
      kel_t <- kel_depur
    }

    # The state is the whole-body Cd concentration CI (mg/kg dry weight),
    # starting from the measured pre-exposure concentration CI0.
    central(0) <- c0
    d/dt(central) <- kassim_t * cfood - kel_t * central

    Cc <- central
    Cc ~ add(addSd)
  })
}
