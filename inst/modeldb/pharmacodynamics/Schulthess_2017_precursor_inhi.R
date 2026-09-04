Schulthess_2017_precursor_inhi <- function() {
  description <- paste(
    "Theoretical (illustrative; no data fitted). Precursor-pool (tolerance",
    "and rebound) model with an INHIBITORY drug effect on the loss of the",
    "precursor, driven by a one-compartment i.v. bolus PK model. Model",
    "flavor 2 of Figure 4a in Schulthess et al. (2017), used to demonstrate",
    "frequency-domain response analysis (FdRA). This is the one model class",
    "in the paper whose frequency response is a BAND-PASS rather than a",
    "low-pass filter: the linearised low-frequency gain is exactly 0, the",
    "amplitude ratio peaks at 0.16 1/h, and only the frequency window",
    "between the two threshold frequencies 0.07 and 0.38 1/h is amplified.",
    "The structure is case study 16 of Gabrielsson and Hjorth (2016) - the",
    "antilipolytic response of healthy volunteers to an adenosine receptor",
    "agonist (Zannikos et al. 2001).",
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
    precursor1 = list(analyte = "precursor pool x1", units = "mL", specimen = "not applicable", verified = TRUE),
    effect     = list(analyte = "biomarker response x2 (the model output)", units = "mL", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "not applicable (theoretical illustration; no subjects, no data fitted)",
    n_subjects = 0,
    n_studies = 0,
    disease_state = "not applicable (generic precursor-dependent turnover process)",
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
      "constant printed in the Figure 4 caption and is therefore encoded",
      "with fixed(). Note that a SINGLE kout governs both the drug-modulated",
      "loss of the precursor and the loss of the response, which is what",
      "makes the two unforced steady states equal and the linearised",
      "low-frequency gain exactly zero."
    )
  )

  ini({
    # ---- One-compartment i.v. bolus PK (Eq. 5: dc/dt = -kel * c) ---------
    # The paper integrates the plasma CONCENTRATION directly and never
    # introduces a volume; the supplementary R script adds the dose to that
    # same state. vc is therefore fixed at 1 so that Cc == central, which
    # keeps the dimensional structure visible without changing any result.
    lvc  <- fixed(log(1))    ; label("Central volume of distribution Vc (L; implied = 1, see notes)")   # Eq. 5 - no volume in the source; dose is added directly to the concentration
    # kel is a SWEPT quantity, not an estimate. The value below is the
    # INTERMEDIATE elimination rate of the Figure 4d time course, which is
    # the one inside the band-pass passband; the vignette sweeps kel over
    # 10^-3..10^1 1/h to build the Bode plot.
    lkel <- fixed(log(1 / 6)); label("Elimination rate constant kel (1/h)")                            # Figure 4d: elimination rates 1/720, 1/6 and 4 1/h; 1/6 is the amplifying case

    # ---- Precursor-pool turnover parameters -----------------------------
    lkin  <- fixed(log(1))    ; label("Zero-order production rate into the precursor pool kin (mL/h)")  # Figure 4 caption: kin = 1 mL/h
    lkout <- fixed(log(1))    ; label("First-order loss rate of precursor and response kout (1/h)")     # Figure 4 caption: kout = 1 1/h

    # ---- Drug function (Eq. 1: E(c) = 1 +/- Emax * c / (EC50 + c)) -------
    lemax <- fixed(log(1))    ; label("Maximum drug effect Emax (unitless)")                            # Figure 4 caption: Emax = 1
    lec50 <- fixed(log(0.25)) ; label("Concentration producing half-maximal effect EC50 (umol/L)")       # Figure 4 caption: EC50 = 0.25 umol/L
  })

  model({
    vc   <- exp(lvc)
    kel  <- exp(lkel)
    kin  <- exp(lkin)
    kout <- exp(lkout)
    emax <- exp(lemax)
    ec50 <- exp(lec50)

    # Unforced (c_SS = 0) stable steady state derived in the source:
    # x1_SS = x2_SS = kin / kout. Both eigenvalues of the Jacobian are
    # -kout, so the steady state is stable for positive rate constants.
    precursor1(0) <- kin / kout
    effect(0)     <- kin / kout

    d/dt(central) <- -kel * central
    Cc <- central / vc

    # Eq. 1 with the inhibitory sign.
    drugEffect <- 1 - emax * Cc / (ec50 + Cc)

    # Eq. 9 exactly:
    #   dx1/dt = kin - kout * x1 * E(c)
    #   dx2/dt = kout * x1 * E(c) - kout * x2
    # The drug-modulated flux leaving the precursor is the flux entering the
    # response, which is why the low-frequency gain vanishes.
    precursorLoss <- kout * precursor1 * drugEffect

    d/dt(precursor1) <- kin - precursorLoss
    d/dt(effect)     <- precursorLoss - kout * effect
  })
}
