Schulthess_2017_directEffect_stim <- function() {
  description <- paste(
    "Theoretical (illustrative; no data fitted). Direct-effect (immediate",
    "Emax) PD model with a STIMULATORY drug function, driven by a",
    "one-compartment i.v. bolus PK model. Supplementary section A.2 and",
    "Figure S1 of Schulthess et al. (2017), included as the degenerate",
    "reference case of the paper's frequency-domain response analysis",
    "(FdRA): because the response has no state of its own, the transfer",
    "function G(s) = +/- Emax/EC50 carries no s at all, so the Bode plot is",
    "a flat line at an amplitude ratio of 4 and every dosing frequency is",
    "passed through unchanged. Contrast this with the four ODE-based PD",
    "classes in the main text, which all filter the input. The inhibitory",
    "flavor differs only in the sign of the drug function and has the same",
    "amplitude ratio, so it is not shipped separately.",
    sep = " "
  )
  reference <- paste(
    "Schulthess P, Post TM, Yates J, van der Graaf PH.",
    "Frequency-domain response analysis for quantitative systems",
    "pharmacology models.",
    "CPT Pharmacometrics Syst Pharmacol. 2018;7(2):111-123.",
    "doi:10.1002/psp4.12266.",
    "Supplementary information section A.2 and Figure S1.",
    sep = " "
  )
  vignette <- "Schulthess_2017_frequency_domain_response_analysis"
  units <- list(
    time = "h",
    dosing = "umol/L",
    concentration = "umol/L"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against the source: the paper doses directly
  # in concentration units (Suppl. Eq. integrates the plasma concentration c
  # itself, and the supplementary R script adds dose = 0.1 to that state), so
  # the implied central volume is 1 and `central` carries a concentration.
  # The response is algebraic (no ODE state) by construction.
  compartmentData <- list(
    central = list(analyte = "drug (generic)", units = "umol/L", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "not applicable (theoretical illustration; no subjects, no data fitted)",
    n_subjects = 0,
    n_studies = 0,
    disease_state = "not applicable (generic direct-effect biomarker)",
    dose_range = paste(
      "illustrative repeated i.v. bolus of 0.1 umol/L given every 4/kel h",
      "(supplementary R script); the Bode-plot analysis sweeps the",
      "elimination rate constant kel over 10^-3 to 10^1 1/h"
    ),
    notes = paste(
      "Schulthess et al. (2017) is a TUTORIAL introducing frequency-domain",
      "response analysis. No clinical or preclinical data were fitted, so",
      "there is no study population, no inter-individual variability and no",
      "residual-error model. Both parameters below are illustrative",
      "constants printed in the Figure S1 caption and are therefore encoded",
      "with fixed(). The response `directEffect` is an algebraic observable",
      "rather than an ODE state, which is precisely why this model has no",
      "frequency-filtering behaviour."
    )
  )

  ini({
    # ---- One-compartment i.v. bolus PK (Suppl. A.2: dc/dt = -kel * c) ----
    # The paper integrates the plasma CONCENTRATION directly and never
    # introduces a volume; the supplementary R script adds the dose to that
    # same state. vc is therefore fixed at 1 so that Cc == central, which
    # keeps the dimensional structure visible without changing any result.
    lvc  <- fixed(log(1))     ; label("Central volume of distribution Vc (L; implied = 1, see notes)")   # Suppl. A.2 - no volume in the source; dose is added directly to the concentration
    # kel is a SWEPT quantity, not an estimate. The value below is the slow
    # elimination rate of the Figure S1a time course; the vignette sweeps
    # kel over 10^-3..10^1 1/h to build the (flat) Bode plot.
    lkel <- fixed(log(1 / 24)); label("Elimination rate constant kel (1/h)")                            # Figure S1a: elimination rates 1/24 and 4 1/h

    # ---- Drug function (Eq. 1: E(c) = 1 +/- Emax * c / (EC50 + c)) -------
    lemax <- fixed(log(1))    ; label("Maximum drug effect Emax (unitless)")                            # Figure S1 caption: Emax = 1
    lec50 <- fixed(log(0.25)) ; label("Concentration producing half-maximal effect EC50 (umol/L)")       # Figure S1 caption: EC50 = 0.25 umol/L
  })

  model({
    vc   <- exp(lvc)
    kel  <- exp(lkel)
    emax <- exp(lemax)
    ec50 <- exp(lec50)

    d/dt(central) <- -kel * central
    Cc <- central / vc

    # Suppl. Eq. in A.2: y = 1 + Emax * c / (EC50 + c). The baseline of 1
    # matches the unforced response of every main-text model at c = 0.
    directEffect <- 1 + emax * Cc / (ec50 + Cc)
  })
}
