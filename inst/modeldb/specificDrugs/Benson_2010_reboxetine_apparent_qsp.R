Benson_2010_reboxetine_apparent_qsp <- function() {
  description <- paste(
    "QSP. In vitro (HEK-293 cell membrane homogenate expressing recombinant",
    "human noradrenaline transporter, hNET). Single-competitor competitive",
    "target-binding kinetics of racemic reboxetine against the radioligand",
    "tracer [3H]-desmethylimipramine (DMI), encoded as Scheme 1 / Equations",
    "1-5 of Benson 2010. Two ligands compete for one shared receptor pool:",
    "free hNET (target) binds the tracer, forming the measured drug-target",
    "complex, and binds one unlabelled drug (complex_reboxetine). This is the",
    "estimation model the paper actually fit, and its racemate parameters are",
    "APPARENT: Benson 2010 fits the racemate as if it were a single ligand,",
    "which it is not, and the Discussion is explicit that the parameter",
    "estimates for the racemate should be regarded as observed and not as",
    "molecularly defined properties. The apparent koff,obs of 6 x 10^-6 1/s",
    "(half time 32 h) is slower, and the apparent Kd,obs of 120 pM about",
    "1.5-fold higher, than the corresponding SS-reboxetine values, because",
    "the faster, weaker RR-enantiomer also shapes the observed kinetics.",
    "Association rates are derived inside model() as kon = koff / Kd",
    "(Equation 7) rather than estimated: Benson 2010 Table 1 reports a",
    "standard error, a %CV and a 95% confidence interval for koff and for Kd",
    "but none of the three for kon. Free tracer IS depleted by binding",
    "(Equation 4) but the unlabelled drug concentration is held constant --",
    "Benson 2010 writes no d[D]/dt equation, and the mass balance of Equation",
    "5 spans only the receptor species. Inter-experiment variability on Bmax",
    "was estimated by the paper but never printed, and the proportional",
    "residual error (Equation 17) has no published sigma; both are carried as",
    "fixed(0) so the structure is recorded without inventing a variance.",
    "Sibling model: Benson_2010_reboxetine_qsp, the mechanistic Scheme 2",
    "(Equations 10-15) model in which SS- and RR-reboxetine compete for hNET",
    "separately with enantiomer-specific rate constants.",
    sep = " "
  )
  reference <- paste(
    "Benson N, Snelder N, Ploeger B, Napier C, Sale H, Birdsall NJM, Butt RP,",
    "van der Graaf PH (2010). Estimation of binding rate constants using a",
    "simultaneous mixed-effects method: application to monoamine transporter",
    "reuptake inhibitor reboxetine. British Journal of Pharmacology",
    "160(2):389-398. doi:10.1111/j.1476-5381.2010.00719.x. PMCID PMC2874860.",
    "Structural equations: Scheme 1 and Equations 1-5 (Methods, Parameter",
    "estimation method). Parameter values: Table 1 (mixed-effects analysis",
    "column, rows koff,obs rac / Kd,obs rac / koff DMI / Kd DMI) and the",
    "Results paragraph reporting Bmax.",
    sep = " "
  )
  vignette <- "Benson_2010_reboxetine_hnet_binding"

  # The assay is a fixed-volume in vitro well, so an "amount" dosed into a
  # free-ligand state IS that ligand's molar concentration; there is no
  # volume term anywhere in Equations 1-5. Dosing and concentration units
  # are therefore the same (nM).
  units <- list(time = "h", dosing = "nM", concentration = "nM")

  covariateData <- list()

  # NAMING NOTE. The canonical TMDD pair `target` / `complex` is assigned to
  # the species the assay actually MEASURES: free hNET, and the hNET complex
  # with the labelled ligand [3H]-DMI. DMI is itself a drug (the active
  # metabolite of imipramine) and the bound-radioligand concentration is the
  # fitted observation in every one of Benson 2010's experiments, so
  # `complex` is a drug-target complex in the ordinary canonical sense. The
  # unlabelled competitor of interest -- racemic reboxetine -- carries the
  # explicit `_reboxetine` suffix so the two complexes cannot be confused.
  # The sibling model Benson_2010_reboxetine_qsp uses the same assignment.
  compartmentData <- list(
    target = list(
      analyte = "free (unoccupied) human noradrenaline transporter, hNET",
      units = "nM", specimen = "not applicable", verified = TRUE
    ),
    tracer = list(
      analyte = "free [3H]-desmethylimipramine (radioligand tracer)",
      units = "nM", specimen = "not applicable", verified = TRUE
    ),
    complex = list(
      analyte = "hNET-[3H]-desmethylimipramine complex (RT); the measured species",
      units = "nM", specimen = "not applicable", verified = TRUE
    ),
    reboxetine = list(
      analyte = "free racemic reboxetine, treated as a single apparent ligand",
      units = "nM", specimen = "not applicable", verified = TRUE
    ),
    complex_reboxetine = list(
      analyte = "hNET-reboxetine complex (RD)",
      units = "nM", specimen = "not applicable", verified = TRUE
    )
  )

  # `target` and `complex` are canonical. The free-ligand pools and the
  # reboxetine-specific complex are paper-mechanistic roles with no canonical
  # entry (compartment-names.md registers only the `_csf` / `_isf` location
  # suffixes on `complex`).
  paper_specific_compartments <- c("tracer", "reboxetine", "complex_reboxetine")

  population <- list(
    species        = "in vitro (HEK-293 cell membrane homogenate expressing recombinant human noradrenaline transporter, hNET; assay buffer 20 mM HEPES, 120 mM NaCl, 5 mM KCl, pH 7.4, room temperature ~25 C)",
    n_subjects     = 53L,
    n_studies      = 1L,
    age_range      = NA_character_,
    weight_range   = NA_character_,
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = "In vitro radioligand binding assay; not a clinical population.",
    dose_range     = paste(
      "Experiment 3 (competition kinetics): racemic reboxetine at 0.1, 0.3",
      "and 1 nM against [3H]-DMI tracer at c. 1 nM and hNET at c. 0.1 nM.",
      "Three racemate data sets contributed to the apparent-racemate fit",
      "(Figure 1C bottom row)."
    ),
    regions        = NA_character_,
    notes          = paste(
      "The mixed-effects 'individual' is an EXPERIMENT, not a subject: all",
      "data from Experiments 1-3 were analysed in a single NONMEM VI FOCE-I",
      "step over 1500 data points from 53 experiments (Benson 2010 Results).",
      "The apparent racemate parameters carry the highest %CV in Table 1 (27%",
      "on Kd,obs, 34% on koff,obs) precisely because a single-ligand model is",
      "being fit to a two-enantiomer mixture. Use Benson_2010_reboxetine_qsp",
      "when a mechanistic description of racemic reboxetine is wanted.",
      "All experiments were run at room temperature (25 C); Benson 2010",
      "Discussion cautions that in vivo rate constants at 37 C would differ."
    ),
    n_experiments  = 53L,
    n_observations = 1500L
  )

  ini({
    # ------------------------------------------------------------------
    # Values are carried EXACTLY as published: dissociation rate constants
    # in s^-1 and equilibrium dissociation constants in nM, so that every
    # ini() line round-trips one Benson 2010 Table 1 cell. The per-second
    # to per-hour conversion and the kon = koff / Kd derivation (Equation 7)
    # both happen inside model().
    #
    # koff and Kd are the ESTIMATED quantities -- Table 1 gives each a
    # standard error, a %CV and a 95% confidence interval. kon is reported
    # in Table 1 without any of the three, so it is the derived quantity.
    # ------------------------------------------------------------------

    lbmax <- log(0.0724)
    label("Total hNET receptor concentration, Bmax (nM)")                        # Results: Bmax 72.4 +/- 4.0 pM = 0.0724 nM; equivalently ~28 pmol/mg protein

    # Tracer: [3H]-desmethylimipramine (DMI)
    lkoff_dmi <- log(2.9e-3)
    label("Tracer [3H]-DMI dissociation rate constant, koff (1/s)")               # Table 1, koff DMI: (2.9 +/- 0.2) x 10^-3 s^-1, CV 7%, 95% CI (2.5-3.3) x 10^-3; half time 0.07 h
    lkd_dmi <- log(1.5)
    label("Tracer [3H]-DMI equilibrium dissociation constant, KT (nM)")           # Table 1, Kd DMI: 1.5 +/- 0.15 nM, CV 10%, 95% CI 1.2-1.8; logKd -8.82 (10^-8.82 M = 1.51 nM)

    # Racemic reboxetine, fit as a single apparent ligand
    lkoff_rac <- log(6.0e-6)
    label("Racemic reboxetine apparent dissociation rate constant, koff,obs (1/s)")      # Table 1, koff,obs rac: (0.6 +/- 0.2) x 10^-5 s^-1 = 6.0e-6 s^-1, CV 34%, 95% CI (0.2-1.0) x 10^-5; half time 32 h
    lkd_rac <- log(0.120)
    label("Racemic reboxetine apparent equilibrium dissociation constant, Kd,obs (nM)")  # Table 1, Kd,obs rac: 120 +/- 30 pM = 0.120 nM, CV 27%, 95% CI 55-180 pM; logKd,obs -9.92 (10^-9.92 M = 0.120 nM)

    # Inter-experiment variability. Benson 2010 Results establishes that IIV
    # on Bmax alone improved the fit by 58 MVOF points (P < 0.001) and
    # Equation 16 gives its exponential form, but the omega is never printed
    # in the paper. Carried as fixed(0) so the structure is preserved
    # without inventing a variance; see the vignette Errata.
    etalbmax ~ fixed(0)                                                           # Results + Equation 16: exponential inter-experiment variability on Bmax, magnitude not published

    # Residual error. Equation 17 is a proportional error model on the
    # observed bound radioligand; sigma is never printed. Same treatment.
    propSd <- fixed(0)
    label("Proportional residual error on bound tracer (fraction)")               # Equation 17: Yij = Fij * (1 + eps_PROP,ij); sigma not published
  })

  model({
    # ------------------------------------------------------------------
    # 1. Unit handling. Table 1 rate constants are per SECOND and the
    #    equilibrium constants are molar; the model integrates on an hour
    #    time base (assay runs 2 min to 24 h) with nM concentrations.
    # ------------------------------------------------------------------
    sec_per_h <- 3600

    bmax <- exp(lbmax + etalbmax)

    # ------------------------------------------------------------------
    # 2. Per-ligand rate constants. Equation 7, Kd = koffD / konD, is
    #    inverted to give konD = koffD / Kd. koff is converted from 1/s to
    #    1/h; Kd stays in nM, so kon comes out in 1/(nM*h).
    # ------------------------------------------------------------------
    koff_dmi <- exp(lkoff_dmi) * sec_per_h
    kd_dmi   <- exp(lkd_dmi)
    kon_dmi  <- koff_dmi / kd_dmi

    koff_rac <- exp(lkoff_rac) * sec_per_h
    kd_rac   <- exp(lkd_rac)
    kon_rac  <- koff_rac / kd_rac

    # ------------------------------------------------------------------
    # 3. Scheme 1 ODE system, Benson 2010 Equations 1-5. Equation 5 is the
    #    receptor mass balance Bmax = [R] + [RT] + [RD]; it is not
    #    integrated, it is a consequence of Equations 1-3 summing to zero
    #    and is asserted in the validation vignette instead.
    # ------------------------------------------------------------------

    # Equation 2: free hNET
    d/dt(target) <-
      -(kon_dmi * tracer + kon_rac * reboxetine) * target +
      koff_dmi * complex + koff_rac * complex_reboxetine

    # Equation 4: free tracer IS depleted by binding
    d/dt(tracer) <- -kon_dmi * target * tracer + koff_dmi * complex

    # Equation 1: hNET-tracer complex [RT], the measured species
    d/dt(complex) <- kon_dmi * target * tracer - koff_dmi * complex

    # Benson 2010 writes no d[D]/dt: the unlabelled competitor concentration
    # is an exogenous input held at its nominal assay value for the whole
    # time course. This is the paper's own assumption, not a simplification
    # introduced here, and it is why Equation 5 balances only the receptor
    # species. See the vignette Errata -- at the 0.1 nM competitor
    # concentration of Experiment 3 the free ligand is of the same order as
    # Bmax (0.0724 nM), so a depleting-ligand variant would not give
    # identical curves.
    d/dt(reboxetine) <- 0

    # Equation 3: hNET-reboxetine complex [RD]
    d/dt(complex_reboxetine) <-
      kon_rac * target * reboxetine - koff_rac * complex_reboxetine

    # ------------------------------------------------------------------
    # 4. Initial condition. Benson 2010 Methods: "Depending on the
    #    experiment, the compartment corresponding to the differential
    #    equation of either [R] or [RT] was initialized with the total
    #    receptor concentration (Bmax)." Every experiment reproduced in the
    #    validation vignette starts from free receptor.
    # ------------------------------------------------------------------
    target(0) <- bmax

    # ------------------------------------------------------------------
    # 5. Observables. The fitted observation is bound radioligand, i.e. the
    #    hNET-tracer complex of Equation 1. The remaining outputs are
    #    reported for interpretation and carry no error model.
    # ------------------------------------------------------------------

    # Fractional hNET occupancy by the apparent racemate ligand.
    occupancyReboxetine <- complex_reboxetine / bmax

    # Receptor mass-balance residual; Equation 5 requires this to stay 0.
    massBalance <- target + complex + complex_reboxetine - bmax

    complex ~ prop(propSd)
  })
}
