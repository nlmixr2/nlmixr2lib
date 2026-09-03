Schulthess_2017_autoreg_stim_posfb <- function() {
  description <- paste(
    "Theoretical (illustrative; no data fitted). Autoregulation model with",
    "POSITIVE (auto-stimulatory) FEEDBACK on production and a STIMULATORY",
    "drug effect on the loss, driven by a one-compartment i.v. bolus PK",
    "model. Model flavor 1 of Figure 3a in Schulthess et al. (2017), one of",
    "four autoregulation flavors used to demonstrate frequency-domain",
    "response analysis (FdRA). The linearised model is a first-order",
    "low-pass filter with a low-frequency gain of 4 and a threshold",
    "frequency of 0.46 1/h; the positive-feedback flavors amplify slow",
    "dosing frequencies, whereas the negative-feedback flavors attenuate the",
    "input at every frequency.",
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
    central = list(analyte = "drug (generic)", units = "umol/L", specimen = "plasma", verified = TRUE),
    effect  = list(analyte = "autoregulated biomarker response x", units = "umol/L", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "not applicable (theoretical illustration; no subjects, no data fitted)",
    n_subjects = 0,
    n_studies = 0,
    disease_state = "not applicable (generic autoregulated turnover process)",
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
      "constant printed in the Figure 3 caption and is therefore encoded",
      "with fixed(). UNITS CAVEAT: the caption gives kin = 1 mL/h and",
      "K = 0.25 umol/L, but the feedback term K/(K + x) requires K and the",
      "response x to share units, so the two cannot both be right as",
      "printed. Every published value is 1 or 0.25 and the arithmetic is",
      "unit-agnostic as implemented, so no numeric result depends on the",
      "resolution; the response is labelled umol/L here for internal",
      "consistency with K."
    )
  )

  ini({
    # ---- One-compartment i.v. bolus PK (Eq. 5: dc/dt = -kel * c) ---------
    # The paper integrates the plasma CONCENTRATION directly and never
    # introduces a volume; the supplementary R script adds the dose to that
    # same state. vc is therefore fixed at 1 so that Cc == central, which
    # keeps the dimensional structure visible without changing any result.
    lvc  <- fixed(log(1))     ; label("Central volume of distribution Vc (L; implied = 1, see notes)")   # Eq. 5 - no volume in the source; dose is added directly to the concentration
    # kel is a SWEPT quantity, not an estimate. The value below is the slow
    # (amplifying) elimination rate of the Figure 3d time course; the
    # vignette sweeps kel over 10^-3..10^1 1/h to build the Bode plot.
    lkel <- fixed(log(1 / 24)); label("Elimination rate constant kel (1/h)")                            # Figure 3d: elimination rate = 1/24 1/h (paired panel uses 4 1/h)

    # ---- Autoregulated turnover parameters ------------------------------
    lkin  <- fixed(log(1))    ; label("Zero-order production rate of the response kin (mL/h as printed; see units caveat)")  # Figure 3 caption: kin = 1 mL/h
    lkout <- fixed(log(1))    ; label("First-order loss rate of the response kout (1/h)")               # Figure 3 caption: kout = 1 1/h
    # K sets the half-saturation of the autoregulatory feedback term. The
    # source calls it a dissociation constant ("Herein, K could be, for
    # example, a dissociation constant"), which is the registered role of
    # the canonical `kd`.
    lkd   <- fixed(log(0.25)) ; label("Feedback half-saturation (dissociation) constant K (umol/L)")     # Figure 3 caption: K = 0.25 umol/L

    # ---- Drug function (Eq. 1: E(c) = 1 +/- Emax * c / (EC50 + c)) -------
    lemax <- fixed(log(1))    ; label("Maximum drug effect Emax (unitless)")                            # Figure 3 caption: Emax = 1
    lec50 <- fixed(log(0.25)) ; label("Concentration producing half-maximal effect EC50 (umol/L)")       # Figure 3 caption: EC50 = 0.25 umol/L
  })

  model({
    vc   <- exp(lvc)
    kel  <- exp(lkel)
    kin  <- exp(lkin)
    kout <- exp(lkout)
    kd   <- exp(lkd)
    emax <- exp(lemax)
    ec50 <- exp(lec50)

    # Unforced (c_SS = 0) stable steady state for the POSITIVE-feedback
    # flavors, as derived in the source: x_SS = kin/kout - K. With the
    # published values this is 1/1 - 0.25 = 0.75.
    effect(0) <- kin / kout - kd

    d/dt(central) <- -kel * central
    Cc <- central / vc

    # Eq. 1 with the stimulatory sign.
    drugEffect <- 1 + emax * Cc / (ec50 + Cc)

    # Positive (auto-stimulatory) feedback term, defined in the source as
    # phi_1(x) = x / (K + x).
    feedback <- effect / (kd + effect)

    # Eq. 6 exactly: dx/dt = kin * phi(x) - kout * x * E(c). The drug acts on
    # the LOSS in all four autoregulation flavors.
    d/dt(effect) <- kin * feedback - kout * effect * drugEffect
  })
}
