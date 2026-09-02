Schulthess_2017_indirect_stim_kout <- function() {
  description <- paste(
    "Theoretical (illustrative; no data fitted). Indirect response (turnover)",
    "model with STIMULATION OF LOSS, driven by a one-compartment i.v.",
    "bolus PK model. Model flavor 3 of Figure 2a in Schulthess et al. (2017),",
    "the worked example used to derive frequency-domain response analysis",
    "(FdRA) step by step. The linearised model is a first-order low-pass",
    "filter: slow dosing frequencies are amplified fourfold and fast ones are",
    "attenuated. The structure is case study 3 of Gabrielsson and Hjorth",
    "(2016) - a compound acting on the rat urinary bladder sphincter via a",
    "stimulatory alpha-2 adrenergic receptor, with voiding volume as the",
    "biomarker - so the response is a volume in mL.",
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
    effect  = list(analyte = "voiding volume (biomarker response x)", units = "mL", specimen = "urine", verified = TRUE)
  )

  population <- list(
    species = "not applicable (theoretical illustration; no subjects, no data fitted)",
    n_subjects = 0,
    n_studies = 0,
    disease_state = "not applicable (generic turnover biomarker)",
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
      "constant printed in the Figure 2 caption and is therefore encoded",
      "with fixed(). The dosing interval is deliberately tied to the",
      "elimination rate (4/kel) to prevent plasma accumulation and reach a",
      "pseudo steady-state oscillation before amplitudes are measured."
    )
  )

  # NOTE on flavor numbering: the main text and the Figure 2 caption number the
  # four flavors DIFFERENTLY (main text: 2 = inhibit production, 3 = stimulate
  # loss; caption: 2 = stimulate loss, 3 = inhibit production). Both orderings
  # describe the same SET of four models. The numbering used here follows the
  # main text and the E1,2 / E3,4 arrow labels in Figure 2a; the mechanism named
  # in `description` is the unambiguous identifier.
  ini({
    # ---- One-compartment i.v. bolus PK (Eq. 5: dc/dt = -kel * c) ---------
    # The paper integrates the plasma CONCENTRATION directly and never
    # introduces a volume; the supplementary R script adds the dose to that
    # same state. vc is therefore fixed at 1 so that Cc == central, which
    # keeps the dimensional structure visible without changing any result.
    lvc  <- fixed(log(1))     ; label("Central volume of distribution Vc (L; implied = 1, see notes)")   # Eq. 5 - no volume in the source; dose is added directly to the concentration
    # kel is a SWEPT quantity, not an estimate. The value below is the slow
    # (amplifying) elimination rate of the Figure 2d time course; the
    # vignette sweeps kel over 10^-3..10^1 1/h to build the Bode plot.
    lkel <- fixed(log(1 / 24)); label("Elimination rate constant kel (1/h)")                            # Figure 2d: elimination rate = 1/24 1/h (paired panel uses 4 1/h)

    # ---- Turnover (indirect response) parameters -------------------------
    lkin  <- fixed(log(1))    ; label("Zero-order production rate of the response kin (mL/h)")          # Figure 2 caption: kin = 1 mL/h
    lkout <- fixed(log(1))    ; label("First-order loss rate of the response kout (1/h)")               # Figure 2 caption: kout = 1 1/h

    # ---- Drug function (Eq. 1: E(c) = 1 +/- Emax * c / (EC50 + c)) -------
    lemax <- fixed(log(1))    ; label("Maximum drug effect Emax (unitless)")                            # Figure 2 caption: Emax = 1
    lec50 <- fixed(log(0.25)) ; label("Concentration producing half-maximal effect EC50 (umol/L)")       # Figure 2 caption: EC50 = 0.25 umol/L
  })

  model({
    vc   <- exp(lvc)
    kel  <- exp(lkel)
    kin  <- exp(lkin)
    kout <- exp(lkout)
    emax <- exp(lemax)
    ec50 <- exp(lec50)

    # Unforced steady state, Eq. 3 of the source: x_SS = kin / kout.
    # Starting here removes the transient the supplementary R script instead
    # burns off with 250 doses from x = 0.
    effect(0) <- kin / kout

    d/dt(central) <- -kel * central
    Cc <- central / vc

    # Eq. 1 with the stimulatory sign.
    drugEffect <- 1 + emax * Cc / (ec50 + Cc)

    # Drug function moved onto the loss term, per the E3,4 arrow position in
    # Figure 2a: dx/dt = kin - kout*x*E(c). The source prints only Eq. 2
    # (production); the loss-flavor ODE follows the same convention the
    # authors use explicitly in Eq. 6 for the autoregulation models.
    d/dt(effect) <- kin - kout * effect * drugEffect
  })
}
