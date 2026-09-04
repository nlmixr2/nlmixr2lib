Schulthess_2017_moderator_stim <- function() {
  description <- paste(
    "Theoretical (illustrative; no data fitted). Moderator-mediated feedback",
    "(tolerance) model with a STIMULATORY drug effect on the production of",
    "the response, driven by a one-compartment i.v. bolus PK model. Model",
    "flavor 1 of Figure 5a in Schulthess et al. (2017), used to demonstrate",
    "frequency-domain response analysis (FdRA). The response x1 stimulates a",
    "single endogenous moderator, which in turn divides the production rate",
    "of the response - a negative feedback loop that develops tolerance. The",
    "frequency response combines a low-pass and a band-pass: the linearised",
    "low-frequency gain is 2, the amplitude ratio peaks at 0.1 1/h, and the",
    "cutoff frequency is 0.04 1/h. The structure is case study 18 of",
    "Gabrielsson and Hjorth (2016) - nicotinic acid inhibiting the",
    "production of non-esterified free fatty acids (Isaksson et al. 2009).",
    sep = " "
  )
  reference <- paste(
    "Schulthess P, Post TM, Yates J, van der Graaf PH.",
    "Frequency-domain response analysis for quantitative systems",
    "pharmacology models.",
    "CPT Pharmacometrics Syst Pharmacol. 2018;7(2):111-123.",
    "doi:10.1002/psp4.12266.",
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
  # in concentration units (Eq. 5 integrates the plasma concentration c
  # itself, and the supplementary R script adds dose = 0.1 to that state), so
  # the implied central volume is 1 and `central` carries a concentration.
  compartmentData <- list(
    central    = list(analyte = "drug (generic)", units = "umol/L", specimen = "plasma", verified = TRUE),
    effect     = list(analyte = "biomarker response x1 (the model output)", units = "mL", specimen = "plasma", verified = TRUE),
    moderator1 = list(analyte = "endogenous moderator x2", units = "mL", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "not applicable (theoretical illustration; no subjects, no data fitted)",
    n_subjects = 0,
    n_studies = 0,
    disease_state = "not applicable (generic moderator-mediated turnover process)",
    dose_range = paste(
      "illustrative repeated i.v. bolus of 0.1 umol/L given every 4/kel h",
      "(supplementary R script); the Bode-plot analysis sweeps the",
      "elimination rate constant kel over 10^-3 to 10^1 1/h"
    ),
    notes = paste(
      "Schulthess et al. (2017) is a TUTORIAL introducing frequency-domain",
      "response analysis. No clinical or preclinical data were fitted, so",
      "there is no study population, no inter-individual variability and no",
      "residual-error model. Every parameter below is an illustrative",
      "constant printed in the Figure 5 caption and is therefore encoded",
      "with fixed(). The moderator is a Gabrielsson-Hjorth tolerance chain:",
      "a first-order delay driven by the response with NO mass transfer,",
      "whose value scales the production rate it modulates."
    )
  )

  ini({
    # ---- One-compartment i.v. bolus PK (Eq. 5: dc/dt = -kel * c) ---------
    # The paper integrates the plasma CONCENTRATION directly and never
    # introduces a volume; the supplementary R script adds the dose to that
    # same state. vc is therefore fixed at 1 so that Cc == central, which
    # keeps the dimensional structure visible without changing any result.
    lvc  <- fixed(log(1))     ; label("Central volume of distribution Vc (L; implied = 1, see notes)")   # Eq. 5 - no volume in the source; dose is added directly to the concentration
    # kel is a SWEPT quantity, not an estimate. The value below is the
    # INTERMEDIATE elimination rate of the Figure 5d time course, the one
    # that shows the strongest amplification; the vignette sweeps kel over
    # 10^-3..10^1 1/h to build the Bode plot.
    lkel <- fixed(log(1 / 12)); label("Elimination rate constant kel (1/h)")                            # Figure 5d: elimination rates 1/1440, 1/12 and 4 1/h; 1/12 is the strongly amplifying case

    # ---- Moderator-mediated turnover parameters -------------------------
    lkin  <- fixed(log(1))    ; label("Turnover rate of the response kin (mL/h)")                       # Figure 5 caption: kin = 1 mL/h
    lkout <- fixed(log(1))    ; label("Fractional turnover rate of the response kout (1/h)")            # Figure 5 caption: kout = 1 1/h
    lktol <- fixed(log(0.25)) ; label("Fractional turnover rate of the moderator ktol (1/h)")           # Figure 5 caption: ktol = 0.25 1/h

    # ---- Drug function (Eq. 1: E(c) = 1 +/- Emax * c / (EC50 + c)) -------
    lemax <- fixed(log(1))    ; label("Maximum drug effect Emax (unitless)")                            # Figure 5 caption: Emax = 1
    lec50 <- fixed(log(0.25)) ; label("Concentration producing half-maximal effect EC50 (umol/L)")       # Figure 5 caption: EC50 = 0.25 umol/L
  })

  model({
    vc   <- exp(lvc)
    kel  <- exp(lkel)
    kin  <- exp(lkin)
    kout <- exp(lkout)
    ktol <- exp(lktol)
    emax <- exp(lemax)
    ec50 <- exp(lec50)

    # Unforced (c_SS = 0) stable steady state derived in the source:
    # x1_SS = x2_SS = sqrt(kin / kout). With the published values this is 1.
    effect(0)     <- sqrt(kin / kout)
    moderator1(0) <- sqrt(kin / kout)

    d/dt(central) <- -kel * central
    Cc <- central / vc

    # Eq. 1 with the stimulatory sign.
    drugEffect <- 1 + emax * Cc / (ec50 + Cc)

    # Eq. 15 exactly:
    #   dx1/dt = (kin / x2) * E(c) - kout * x1
    #   dx2/dt = ktol * x1 - ktol * x2
    d/dt(effect)     <- (kin / moderator1) * drugEffect - kout * effect
    d/dt(moderator1) <- ktol * (effect - moderator1)
  })
}
