Kim_2023_levonorgestrel_firstOrderRelease <- function() {
  description <- paste(
    "First-order in vivo release model for the levonorgestrel (LNG) 150 mg",
    "two-rod subdermal implant (Jadelle), fitted to explant residual-content",
    "data and Jadelle label release rates. Release-only: there is no",
    "levonorgestrel disposition. Only a fraction fdepot of the loaded dose is",
    "releasable (the paper's Qmax = (1 - f) * Dose), and that releasable",
    "amount leaves the implant first-order at rate ka. The outputs are the",
    "cumulative amount released and the instantaneous release rate. Kim 2023",
    "used this function as the input to a PK-Sim / MoBi whole-body PBPK model",
    "of levonorgestrel, which is not part of this model; the biexponential",
    "alternative the paper carried forward alongside it is",
    "Kim_2023_levonorgestrel_biexpRelease."
  )
  reference <- paste(
    "Kim S, Cicali B, Pressly M, Da Silva L, Wendl T, Vozmediano V, Schmidt S,",
    "Cristofoletti R. Model-Based Analysis of In Vivo Release Data of",
    "Levonorgestrel Implants: Projecting Long-Term Systemic Exposure.",
    "Pharmaceutics. 2023;15(5):1393. doi:10.3390/pharmaceutics15051393.",
    "Parameter values from Table 2 (duplicated in Supplementary Table S1);",
    "release function from Table 1 (Weibull row, special case l = 1) and",
    "Equation (1)."
  )
  vignette <- "Kim_2023_levonorgestrel_implantRelease"
  units <- list(
    time = "day",
    dosing = "mg (levonorgestrel loaded into the implant)",
    concentration = "not applicable (no disposition); released amount in mg, release rate in mg/day"
  )

  # Dose the whole implant load (150 mg for Jadelle) into depot at insertion;
  # f(depot) keeps only the releasable fraction.
  dosing <- "depot"

  # Bookkeeping sink that integrates the release flux, so the cumulative
  # released amount Q(t) of Table 1 can be read off the solve. It stands in for
  # the systemic circulation the released drug enters; it is not a biological
  # compartment.
  paper_specific_compartments <- c("released")

  compartmentData <- list(
    depot = list(analyte = "levonorgestrel", units = "mg", specimen = "administration site", verified = TRUE),
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
    # Qmax is the releasable amount (mg) with Dose = 150 mg fixed; fdepot is
    # their ratio, the paper's (1 - f), so the dose-driven model reproduces
    # Q(t) = Qmax * (1 - exp(-k1 * t)) for a 150 mg implant.
    lfdepot <- log(67.86 / 150); label("Releasable fraction of the loaded dose, Qmax / Dose (fraction)") # Table 2, First-order: Qmax = 67.86 mg (95% CI 49.28-86.44), Dose = 150 mg fixed; text: 45.3% of original LNG load
    lka <- log(0.0009); label("First-order release rate constant k1 (1/day)") # Table 2, First-order: k1 = 0.0009 /day (95% CI 0.00051-0.0014)
  })

  model({
    fdepot <- exp(lfdepot)
    ka <- exp(lka)

    d/dt(depot) <- -ka * depot
    d/dt(released) <- ka * depot
    f(depot) <- fdepot

    # Release rate in mg/day, the quantity of Figure 2 columns 2 and 4.
    releaseRate <- ka * depot
  })
}
