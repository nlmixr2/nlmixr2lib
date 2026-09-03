Benson_2010_reboxetine_qsp <- function() {
  description <- paste(
    "QSP. In vitro (HEK-293 cell membrane homogenate expressing recombinant",
    "human noradrenaline transporter, hNET). Enantiomer-resolved competitive",
    "target-binding kinetics of SS-reboxetine and RR-reboxetine against the",
    "radioligand tracer [3H]-desmethylimipramine (DMI), encoded as Scheme 2 /",
    "Equations 10-15 of Benson 2010. Three ligands compete for one shared",
    "receptor pool: free hNET (target) binds tracer (complex),",
    "SS-reboxetine (complex_ss) and RR-reboxetine (complex_rr), each with its",
    "own association and dissociation rate constant, under the mutually",
    "exclusive binding assumption of Scheme 2. Association rate constants are",
    "derived inside model() as kon = koff / Kd (Equation 7) rather than",
    "estimated: Benson 2010 Table 1 reports a standard error and a 95%",
    "confidence interval for koff and for Kd but none for kon, so koff and Kd",
    "are the estimated quantities and kon is the derived one. Racemic",
    "reboxetine is simulated by dosing an equimolar SS + RR mixture, which is",
    "how the paper reproduces its Figure 3; setting one enantiomer's",
    "concentration to zero recovers the single-competitor Scheme 1 fit for the",
    "other enantiomer. Free tracer IS depleted by binding (Equation 14), but",
    "the unlabelled enantiomer concentrations are held constant -- Benson 2010",
    "writes no d[SS]/dt or d[RR]/dt equation, and the receptor mass balance of",
    "Equation 15 spans only the receptor species. Inter-experiment variability",
    "on Bmax was estimated by the paper but its magnitude was never printed,",
    "and the proportional residual error (Equation 17) likewise has no",
    "published sigma; both are carried as fixed(0) so the structure is",
    "recorded without inventing a variance. Sibling model:",
    "Benson_2010_reboxetine_apparent_qsp, the Scheme 1 (Equations 1-5) fit in",
    "which racemic reboxetine is treated as a single apparent ligand.",
    sep = " "
  )
  reference <- paste(
    "Benson N, Snelder N, Ploeger B, Napier C, Sale H, Birdsall NJM, Butt RP,",
    "van der Graaf PH (2010). Estimation of binding rate constants using a",
    "simultaneous mixed-effects method: application to monoamine transporter",
    "reuptake inhibitor reboxetine. British Journal of Pharmacology",
    "160(2):389-398. doi:10.1111/j.1476-5381.2010.00719.x. PMCID PMC2874860.",
    "Structural equations: Scheme 2 and Equations 10-15 (Methods, Mutually",
    "exclusive binding of enantiomers). Parameter values: Table 1",
    "(mixed-effects analysis column) and the Results paragraph reporting Bmax.",
    sep = " "
  )
  vignette <- "Benson_2010_reboxetine_hnet_binding"

  # The assay is a fixed-volume in vitro well, so an "amount" dosed into a
  # free-ligand state IS that ligand's molar concentration; there is no
  # volume term anywhere in Equations 10-15. Dosing and concentration units
  # are therefore the same (nM).
  units <- list(time = "h", dosing = "nM", concentration = "nM")

  # Benson 2010 uses no covariates: every experiment shares one parameter
  # set and the only stochastic element is inter-experiment variability on
  # Bmax. Ligand concentrations are experimental design inputs and enter as
  # doses into the free-ligand states, not as covariate columns.
  covariateData <- list()

  # States of Equations 10-15. Every state is a molar concentration in the
  # assay buffer, so `units` is nM for all of them and `specimen` is
  # "not applicable" -- a cell-free membrane-homogenate binding assay is not
  # one of the biological matrices in conventions$specimenVocabulary.
  #
  # NAMING NOTE. The canonical TMDD pair `target` / `complex` is assigned to
  # the species the assay actually MEASURES: free hNET, and the hNET complex
  # with the labelled ligand [3H]-DMI. DMI is itself a drug (the active
  # metabolite of imipramine) and the bound-radioligand concentration is the
  # fitted observation in every one of Benson 2010's experiments, so
  # `complex` is a drug-target complex in the ordinary canonical sense. The
  # unlabelled competitors of interest carry explicit `_ss` / `_rr` suffixes
  # so the three complexes can never be confused.
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
    reboxetine_ss = list(
      analyte = "free SS-reboxetine",
      units = "nM", specimen = "not applicable", verified = TRUE
    ),
    complex_ss = list(
      analyte = "hNET-SS-reboxetine complex (RSS)",
      units = "nM", specimen = "not applicable", verified = TRUE
    ),
    reboxetine_rr = list(
      analyte = "free RR-reboxetine",
      units = "nM", specimen = "not applicable", verified = TRUE
    ),
    complex_rr = list(
      analyte = "hNET-RR-reboxetine complex (RRR)",
      units = "nM", specimen = "not applicable", verified = TRUE
    )
  )

  # `target` and `complex` are canonical. The free-ligand pools and the two
  # enantiomer-specific complexes are paper-mechanistic roles with no
  # canonical entry (compartment-names.md registers only the `_csf` / `_isf`
  # location suffixes on `complex`, and there is no canonical name for a
  # radioligand tracer pool).
  paper_specific_compartments <- c(
    "tracer",
    "reboxetine_ss", "complex_ss",
    "reboxetine_rr", "complex_rr"
  )

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
      "Experiment 3 (competition kinetics, the source of the reboxetine",
      "estimates): unlabelled competitor at 0.1, 0.3 and 1 nM for racemic",
      "reboxetine and SS-reboxetine, and 1, 3, 10 and 30 nM for RR-reboxetine,",
      "against [3H]-DMI tracer at c. 1 nM and hNET at c. 0.1 nM. Experiment 1",
      "(tracer saturation): [3H]-DMI 30 pM to 20 nM against c. 0.5 nM hNET.",
      "Experiment 2 (tracer association): [3H]-DMI typically 1.5 nM against",
      "c. 0.1 nM hNET. [3H]-SS-reboxetine saturation: 0.004 to 16 nM against",
      "c. 0.12 nM hNET."
    ),
    regions        = NA_character_,
    notes          = paste(
      "The mixed-effects 'individual' is an EXPERIMENT, not a subject: all",
      "data from Experiments 1-3 were analysed in a single NONMEM VI FOCE-I",
      "step over 1500 data points from 53 experiments (Benson 2010 Results).",
      "Sampling ran from 2 min to 24 h (Experiment 3 harvest times 2, 4, 6, 8,",
      "10, 12, 14, 16, 18, 20, 25, 30, 35, 40, 45, 50, 55, 60, 90, 120 and",
      "150 min and 4, 6, 8, 20, 22 and 24 h). The paper notes that observed",
      "bound radioligand around 24 h ran below prediction in Experiments 2 and",
      "3, most likely from hNET degradation in the assay, so the model is not",
      "expected to track the extreme tail. All experiments were run at room",
      "temperature (25 C); Benson 2010 Discussion cautions that in vivo rate",
      "constants at 37 C would differ.",
      "The membrane preparation was stored at 5 mg/mL protein with Bmax of",
      "approximately 40 pmol/mg; the fitted Bmax of 72.4 pM corresponds to",
      "approximately 28 pmol/mg."
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
    label("Total hNET receptor concentration, Bmax (nM)")                        # Results: Bmax 72.4 +/- 4.0 pM (CV not tabulated) = 0.0724 nM; equivalently ~28 pmol/mg protein

    # Tracer: [3H]-desmethylimipramine (DMI)
    lkoff_dmi <- log(2.9e-3)
    label("Tracer [3H]-DMI dissociation rate constant, koff (1/s)")               # Table 1, koff DMI: (2.9 +/- 0.2) x 10^-3 s^-1, CV 7%, 95% CI (2.5-3.3) x 10^-3; half time 0.07 h
    lkd_dmi <- log(1.5)
    label("Tracer [3H]-DMI equilibrium dissociation constant, KT (nM)")           # Table 1, Kd DMI: 1.5 +/- 0.15 nM, CV 10%, 95% CI 1.2-1.8; logKd -8.82 (10^-8.82 M = 1.51 nM)

    # SS-reboxetine (the eutomer)
    lkoff_ss <- log(1.05e-5)
    label("SS-reboxetine dissociation rate constant, koff (1/s)")                 # Table 1, koff SS: (1.05 +/- 0.07) x 10^-5 s^-1, CV 7%, 95% CI (0.9-1.2) x 10^-5; half time 18 h
    lkd_ss <- log(0.076)
    label("SS-reboxetine equilibrium dissociation constant, Kd (nM)")             # Table 1, Kd SS: 76 +/- 9 pM, CV 12%, 95% CI 57-94 pM; logKd -10.12 (10^-10.12 M = 0.0759 nM)

    # RR-reboxetine (the distomer)
    lkoff_rr <- log(4.2e-3)
    label("RR-reboxetine dissociation rate constant, koff (1/s)")                 # Table 1, koff RR: (4.2 +/- 0.8) x 10^-3 s^-1, CV 18%, 95% CI (2.7-5.7) x 10^-3; half time 0.05 h
    lkd_rr <- log(9.8)
    label("RR-reboxetine equilibrium dissociation constant, Kd (nM)")             # Table 1, Kd RR: 9.8 +/- 0.8 nM, CV 9%, 95% CI 8.1-11.4; logKd -8.01 (10^-8.01 M = 9.77 nM)

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

    koff_ss <- exp(lkoff_ss) * sec_per_h
    kd_ss   <- exp(lkd_ss)
    kon_ss  <- koff_ss / kd_ss

    koff_rr <- exp(lkoff_rr) * sec_per_h
    kd_rr   <- exp(lkd_rr)
    kon_rr  <- koff_rr / kd_rr

    # ------------------------------------------------------------------
    # 3. Scheme 2 ODE system, Benson 2010 Equations 10-15.
    #
    #    Equation 11 (free receptor) is the sink/source balance for all
    #    three ligands. Equation 15 is the receptor mass balance
    #    Bmax = [R] + [RT] + [RRR] + [RSS]; it is not integrated, it is a
    #    consequence of Equations 10-13 summing to zero and is asserted in
    #    the validation vignette instead.
    #
    #    ERRATUM: Equation 13 as printed reads
    #      d[RRR]/dt = konRR * [R] * [RR] - koffSS * [RRR]
    #    with koffSS on the loss term. That is a typesetting error. Equation
    #    11 carries "+ koffRR * [RRR]" as the matching source term for free
    #    receptor, the Methods variable list defines koffRR as "dissociation
    #    rate constant of RR-reboxetine" and uses it nowhere else, and the
    #    Equation 15 mass balance only closes when the two rates are equal:
    #    summing Equations 10-13 as printed leaves a residual
    #    (koffRR - koffSS) * [RRR] instead of zero. koffRR is used here.
    # ------------------------------------------------------------------

    # Equation 11: free hNET
    d/dt(target) <-
      -(kon_dmi * tracer + kon_ss * reboxetine_ss + kon_rr * reboxetine_rr) * target +
      koff_dmi * complex + koff_ss * complex_ss + koff_rr * complex_rr

    # Equation 14: free tracer IS depleted by binding
    d/dt(tracer) <- -kon_dmi * target * tracer + koff_dmi * complex

    # Equation 10: hNET-tracer complex
    d/dt(complex) <- kon_dmi * target * tracer - koff_dmi * complex

    # Benson 2010 writes no d[SS]/dt or d[RR]/dt: the unlabelled competitor
    # concentrations are exogenous inputs held at their nominal assay value
    # for the whole time course, so their derivatives are zero. This is the
    # paper's own assumption, not a simplification introduced here, and it
    # is why Equation 15 balances only the receptor species. See the
    # vignette Errata -- at the 0.1 nM competitor concentrations of
    # Experiment 3 the free ligand is of the same order as Bmax (0.0724 nM),
    # so a depleting-ligand variant would not give identical curves.
    d/dt(reboxetine_ss) <- 0

    # Equation 12: hNET-SS-reboxetine complex
    d/dt(complex_ss) <- kon_ss * target * reboxetine_ss - koff_ss * complex_ss

    d/dt(reboxetine_rr) <- 0

    # Equation 13 (with koffRR per the erratum note above)
    d/dt(complex_rr) <- kon_rr * target * reboxetine_rr - koff_rr * complex_rr

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
    #    hNET-tracer complex of Equation 10. The remaining outputs are
    #    reported for interpretation and carry no error model.
    # ------------------------------------------------------------------

    # Fractional hNET occupancy by reboxetine, the quantity Equation 9
    # gives in closed form at steady state.
    occupancyReboxetine <- (complex_ss + complex_rr) / bmax

    # Receptor mass-balance residual; Equation 15 requires this to stay 0.
    massBalance <- target + complex + complex_ss + complex_rr - bmax

    complex ~ prop(propSd)
  })
}
