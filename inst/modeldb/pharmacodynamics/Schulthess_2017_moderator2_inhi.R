Schulthess_2017_moderator2_inhi <- function() {
  description <- paste(
    "Theoretical (illustrative; no data fitted). DOUBLE moderator-mediated",
    "feedback (tolerance) model with an INHIBITORY drug effect on the",
    "production of the response, driven by a one-compartment i.v. bolus PK",
    "model. Model flavor 2 of Figure 6a in Schulthess et al. (2017), used to",
    "demonstrate frequency-domain response analysis (FdRA). Identical to the",
    "single-moderator model except that the tolerance chain has TWO",
    "moderator compartments in series, which the paper shows barely changes",
    "the frequency response: the linearised low-frequency gain is still 2,",
    "the amplitude ratio peaks at 0.06 1/h and the threshold frequency stays",
    "at 0.62 1/h. The structure is case study 19 of Gabrielsson and Hjorth",
    "(2016) - fold mRNA induction of a specific target by an anonymised test",
    "compound.",
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
    effect     = list(analyte = "fold mRNA induction x1 (the model output)", units = "mL", specimen = "not applicable", verified = TRUE),
    moderator1 = list(analyte = "first endogenous moderator x2", units = "mL", specimen = "not applicable", verified = TRUE),
    moderator2 = list(analyte = "second endogenous moderator x3", units = "mL", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "not applicable (theoretical illustration; no subjects, no data fitted)",
    n_subjects = 0,
    n_studies = 0,
    disease_state = "not applicable (generic moderator-mediated gene-induction process)",
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
      "constant printed in the Figure 6 caption and is therefore encoded",
      "with fixed(). Both moderators share the single rate constant ktol,",
      "and it is the SECOND (terminal) moderator that divides the production",
      "rate of the response."
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
    # INTERMEDIATE elimination rate of the Figure 6d time course, the one
    # that shows the fourfold amplification; the vignette sweeps kel over
    # 10^-3..10^1 1/h to build the Bode plot.
    lkel <- fixed(log(1 / 12)); label("Elimination rate constant kel (1/h)")                            # Figure 6d: elimination rates 1/1440, 1/12 and 4 1/h; 1/12 is the strongly amplifying case

    # ---- Moderator-mediated turnover parameters -------------------------
    lkin  <- fixed(log(1))    ; label("Turnover rate of the response kin (mL/h)")                       # Figure 6 caption: kin = 1 mL/h
    lkout <- fixed(log(1))    ; label("Fractional turnover rate of the response kout (1/h)")            # Figure 6 caption: kout = 1 1/h
    lktol <- fixed(log(0.25)) ; label("Fractional turnover rate of both moderators ktol (1/h)")         # Figure 6 caption: ktol = 0.25 1/h

    # ---- Drug function (Eq. 1: E(c) = 1 +/- Emax * c / (EC50 + c)) -------
    lemax <- fixed(log(1))    ; label("Maximum drug effect Emax (unitless)")                            # Figure 6 caption: Emax = 1
    lec50 <- fixed(log(0.25)) ; label("Concentration producing half-maximal effect EC50 (umol/L)")       # Figure 6 caption: EC50 = 0.25 umol/L
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
    # x1_SS = x2_SS = x3_SS = sqrt(kin / kout). With the published values
    # this is 1.
    effect(0)     <- sqrt(kin / kout)
    moderator1(0) <- sqrt(kin / kout)
    moderator2(0) <- sqrt(kin / kout)

    d/dt(central) <- -kel * central
    Cc <- central / vc

    # Eq. 1 with the inhibitory sign.
    drugEffect <- 1 - emax * Cc / (ec50 + Cc)

    # Eq. 19 exactly:
    #   dx1/dt = (kin / x3) * E(c) - kout * x1
    #   dx2/dt = ktol * (x1 - x2)
    #   dx3/dt = ktol * (x2 - x3)
    # Note the TERMINAL moderator (x3) is the one that divides kin.
    d/dt(effect)     <- (kin / moderator2) * drugEffect - kout * effect
    d/dt(moderator1) <- ktol * (effect - moderator1)
    d/dt(moderator2) <- ktol * (moderator1 - moderator2)
  })
}
