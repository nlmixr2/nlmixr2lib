Kim_2023_levonorgestrel_biexpRelease <- function() {
  description <- paste(
    "Biexponential in vivo release model for the levonorgestrel (LNG) 150 mg",
    "two-rod subdermal implant (Jadelle), fitted to explant residual-content",
    "data and Jadelle label release rates. Release-only: there is no",
    "levonorgestrel disposition. The loaded dose is split between a fast pool",
    "(fraction frel, the paper's Qb1 / Dose) released first-order at ka_fast",
    "and a slow pool (1 - frel) released first-order at ka_slow, so the whole",
    "load is eventually released. The implant is dosed twice at insertion,",
    "once into depot and once into depot2, each with the full load. The",
    "outputs are the cumulative amount released and the instantaneous release",
    "rate. Kim 2023 used this function as the input to a PK-Sim / MoBi",
    "whole-body PBPK model of levonorgestrel, which is not part of this model;",
    "the first-order alternative the paper carried forward alongside it is",
    "Kim_2023_levonorgestrel_firstOrderRelease."
  )
  reference <- paste(
    "Kim S, Cicali B, Pressly M, Da Silva L, Wendl T, Vozmediano V, Schmidt S,",
    "Cristofoletti R. Model-Based Analysis of In Vivo Release Data of",
    "Levonorgestrel Implants: Projecting Long-Term Systemic Exposure.",
    "Pharmaceutics. 2023;15(5):1393. doi:10.3390/pharmaceutics15051393.",
    "Parameter values from Table 2 (duplicated in Supplementary Table S1);",
    "release function from Table 1 (Biexponential row)."
  )
  vignette <- "Kim_2023_levonorgestrel_implantRelease"
  units <- list(
    time = "day",
    dosing = "mg (levonorgestrel loaded into the implant)",
    concentration = "not applicable (no disposition); released amount in mg, release rate in mg/day"
  )

  # Dose the whole implant load (150 mg for Jadelle) into BOTH depots at
  # insertion; f(depot) and f(depot2) split it between the two pools.
  dosing <- c("depot", "depot2")

  # Bookkeeping sink that integrates the release flux, so the cumulative
  # released amount Q(t) of Table 1 can be read off the solve. It stands in for
  # the systemic circulation the released drug enters; it is not a biological
  # compartment.
  paper_specific_compartments <- c("released")

  compartmentData <- list(
    depot = list(analyte = "levonorgestrel", units = "mg", specimen = "administration site", verified = TRUE),
    depot2 = list(analyte = "levonorgestrel", units = "mg", specimen = "administration site", verified = TRUE),
    released = list(analyte = "levonorgestrel", units = "mg", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 3L,
    studies = c(
      "Callahan 2015 (Contraception 92:350-355): LNG content of explanted implants, mean +/- SE by time since insertion",
      "Fuchs 2020 (Contraception X 2:100039): LNG content of explanted two-rod implants up to 5 years",
      "Jadelle US product label: in vivo release rates at 1, 12 and 24 months, converted to released amount"
    ),
    sex_female_pct = 100,
    disease_state = "Healthy women using the levonorgestrel 150 mg two-rod subdermal implant (Jadelle) for contraception",
    dose_range = "One 150 mg two-rod implant (2 x 75 mg), up to 5 years of use",
    notes = paste(
      "The release model was fitted to group-mean explant data, not to",
      "individual subjects; the paper does not report how many explants",
      "contributed (n_subjects is NA). The abstract's n = 321 is the Jadelle",
      "clinical-trial PK cohort used to check the PBPK model, not the release",
      "fit. Fitted by unweighted nonlinear least squares in MATLAB (Section",
      "2.2); no residual-error variance and no inter-individual variability",
      "were reported, so the model is deterministic."
    )
  )

  ini({
    # Qb1 is the fast-pool amount (mg) with Dose = 150 mg fixed; frel is their
    # ratio, so the dose-driven model reproduces
    # Q(t) = Dose - Qb1 * exp(-kb1 * t) - (Dose - Qb1) * exp(-kb2 * t).
    logitfrel <- log((44.64 / 150) / (1 - 44.64 / 150)); label("Fraction of the loaded dose in the fast-release pool, Qb1 / Dose (fraction)") # Table 2, Biexponential: Qb1 = 44.64 mg (95% CI -134.6 to 223.9), Dose = 150 mg fixed
    lka_fast <- log(0.00127); label("Fast-pool first-order release rate constant k1b (1/day)") # Table 2, Biexponential: k1b = 0.00127 /day (95% CI -0.002461 to 0.005011)
    lka_slow <- log(8.985e-5); label("Slow-pool first-order release rate constant k2b (1/day)") # Table 2, Biexponential: k2b = 8.985e-5 /day (95% CI -0.0005871 to 0.0007668)
  })

  model({
    frel <- expit(logitfrel)
    ka_fast <- exp(lka_fast)
    ka_slow <- exp(lka_slow)

    d/dt(depot) <- -ka_fast * depot
    d/dt(depot2) <- -ka_slow * depot2
    d/dt(released) <- ka_fast * depot + ka_slow * depot2
    f(depot) <- frel
    f(depot2) <- 1 - frel

    # Release rate in mg/day, the quantity of Figure 2 columns 2 and 4.
    releaseRate <- ka_fast * depot + ka_slow * depot2
  })
}
