Han_2023_budesonide_dissolution <- function() {
  description <- paste(
    "In vitro (biorelevant FaSSIF dissolution of Entocort EC pellets).",
    "Weibull dissolution (release) function for the budesonide",
    "controlled-release capsule Entocort EC, fitted by the authors to the",
    "in vitro dissolution profile of the pellets and used as the formulation",
    "input of their Simcyp whole-body PBPK model. Release-only: there is no",
    "budesonide disposition. The undissolved drug in the pellets (depot)",
    "dissolves with the Weibull hazard; the outputs are the cumulative amount",
    "dissolved (released) and the instantaneous dissolution rate. The",
    "Simcyp PBPK model the function fed (multilayer gut wall ADAM absorption,",
    "Rodgers and Rowland tissue distribution, CYP3A4 clearance) is not part",
    "of this model. The paper's six virtual test formulations (Table 2) are",
    "reproduced by overriding lra and lgam1."
  )
  reference <- paste(
    "Han C, Sun T, Chirumamilla SK, Bois FY, Xu M, Rostami-Hodjegan A.",
    "Understanding Discordance between In Vitro Dissolution, Local Gut and",
    "Systemic Bioequivalence of Budesonide in Healthy and Crohn's Disease",
    "Patients through PBPK Modeling. Pharmaceutics. 2023;15(9):2237.",
    "doi:10.3390/pharmaceutics15092237. Weibull form from Equation (4);",
    "Entocort EC parameter values from Supplementary Table S1 (to five",
    "significant figures; Table 2 prints the same fit rounded to 3.12 and",
    "0.94)."
  )
  vignette <- "Han_2023_budesonide_dissolution"
  units <- list(
    time = "h",
    dosing = "mg (budesonide pellet content of the capsule, e.g. 3 mg per Entocort EC capsule)",
    concentration = "not applicable (no disposition); dissolved amount in mg, dissolution rate in mg/h"
  )

  # Dose the capsule content into depot at the start of dissolution (time 0
  # of the in vitro profile).
  dosing <- "depot"

  # Bookkeeping sink that integrates the dissolution flux, so the cumulative
  # amount dissolved F(t) of Equation (4) can be read off the solve. It is
  # the dissolution medium, not a biological compartment.
  paper_specific_compartments <- c("released")

  compartmentData <- list(
    depot = list(analyte = "budesonide", units = "mg", specimen = "administration site", verified = TRUE),
    released = list(analyte = "budesonide", units = "mg", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "in vitro (Entocort EC budesonide pellets in FaSSIF dissolution medium)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    studies = paste(
      "Entocort EC in vitro dissolution profile in fasted-state simulated",
      "intestinal fluid (FaSSIF), taken from Effinger et al. Eur J Pharm Sci",
      "2021;157:105617 (Han 2023 reference 22); sampled every 0.5 h to 6 h,",
      "with 100% dissolution assumed at the end of the profile (Section",
      "2.4.4, Figure 2, Discussion)."
    ),
    disease_state = paste(
      "Not applicable (in vitro). The PBPK model the function fed was applied",
      "to healthy volunteers and to Crohn's disease patients."
    ),
    dose_range = "3 mg budesonide per Entocort EC capsule; the release fraction is dose-independent",
    notes = paste(
      "A deterministic curve fit to one mean dissolution profile, so the model",
      "carries no inter-individual variability and no residual error. Fmax",
      "is 100% for every formulation in Table 2 and is carried as a fixed",
      "lfdepot of log(1). The Weibull lag time Tlag of Equation (4) is not",
      "reported for any formulation and is taken as zero, which is what the",
      "fitted values require to reproduce Figure 2. Release in vivo starts",
      "only once the enteric coating dissolves at pH above the trigger pH",
      "(5.5 for Entocort EC, i.e. after the pellets leave the stomach); that",
      "trigger, the gastric and small-intestinal pellet residence times and",
      "everything downstream of dissolution belong to the Simcyp PBPK model",
      "and are not represented here. Table 2 virtual formulations (Weibull",
      "alpha, beta): +20% (1.53, 0.87), +10% (2.14, 0.89), +5% (2.56, 0.91),",
      "+3% (2.77, 0.92), -5% (3.90, 1.01), -10% (5.01, 1.09); the 'pH",
      "threshold = 5' and 'pH threshold = 6' formulations share the",
      "Entocort EC Weibull parameters and differ only in trigger pH."
    )
  )

  ini({
    # Equation (4) is F = Fmax * (1 - exp(-(t - Tlag)^beta / alpha)), with
    # alpha dividing t^beta. That is the canonical Weibull release input
    # S(t) = exp(-(ra * t)^gam1) with gam1 = beta and ra = alpha^(-1/beta),
    # so neither ini() value below is a number printed in the source.
    lra <- log(3.1196^(-1 / 0.93998)); label("Weibull dissolution rate-scaling parameter, alpha^(-1/beta) (1/h)") # Table S1 'Alpha' = 3.1196 with 'Beta' = 0.93998 (Table 2, Entocort EC: alpha = 3.12, beta = 0.94); ra = 3.1196^(-1/0.93998) = 0.298 1/h
    lgam1 <- log(0.93998); label("Weibull dissolution shape parameter beta (unitless)") # Table S1 'Beta' = 0.93998 (Table 2, Entocort EC: beta = 0.94)
    lfdepot <- fixed(log(1)); label("Maximum fraction dissolved, Fmax / 100 (fraction)") # Table S1 'Fmax' = 100 and Table 2 'Fmax' = 100 for every formulation
  })

  model({
    ra <- exp(lra)
    gam1 <- exp(lgam1)
    fdepot <- exp(lfdepot)

    # Weibull hazard of S(t) = exp(-(ra * t)^gam1). Time after dose is floored
    # at a small positive number because (ra * t)^(gam1 - 1) is not finite at
    # t = 0 for gam1 < 1 and is evaluated on every record.
    tr <- max(tad(depot), 1e-6)
    hz <- gam1 * ra * (ra * tr)^(gam1 - 1)

    d/dt(depot) <- -hz * depot
    d/dt(released) <- hz * depot
    f(depot) <- fdepot

    # Dissolution rate in mg/h; released / dose * 100 is the percent
    # dissolved plotted in Figure 2.
    releaseRate <- hz * depot
  })
}
