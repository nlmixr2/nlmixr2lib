Zhou_2015_bms911543_hlm <- function() {
  description <- paste(
    "In vitro (pooled human liver microsomes). Michaelis-Menten enzyme-kinetic",
    "model of the oxidative biotransformation of the JAK2 inhibitor BMS-911543",
    "to its metabolite M1, the only drug-related component detected in pooled",
    "human liver microsomes and the primary biotransformation pathway for the",
    "compound. The published fit is a plain two-parameter Michaelis-Menten",
    "velocity, v = Vmax * S / (Km + S), with Vmax expressed per mg of",
    "microsomal protein. Because the microsomal protein concentration of the",
    "incubation IS reported (0.25 mg/mL), this file carries the full",
    "substrate-depletion ODE for the incubation in addition to the fitted",
    "velocity; its three recombinant-enzyme siblings",
    "(Zhou_2015_bms911543_rcyp1a2, Zhou_2015_bms911543_rcyp3a4,",
    "Zhou_2015_bms911543_rcyp2j2) cannot, because the recombinant CYP content",
    "of those incubations is not reported, and are static velocity models.",
    "A fourth sibling, Zhou_2015_bms911543_cyp1a2_tdi, carries the",
    "time-dependent inactivation of CYP1A2 by the same compound.",
    "The Simcyp V12 whole-body PBPK model that consumes these constants is NOT",
    "part of this file and is not reproducible from the published inputs: the",
    "paper reports no tissue volumes, no organ blood flows, no MPPGL, no liver",
    "weight and no hepatic CYP abundances, and the platform's whole-body ODEs",
    "are the vendor's rather than the authors'. The metabolic scaling factor",
    "of 4 that the authors applied to Km (input as Km/4) inside that platform",
    "model to close an in-vitro-to-in-vivo extrapolation gap is a PBPK fitting",
    "device, not an in-vitro measurement, and is deliberately NOT applied to",
    "the Km carried here. See the validation vignette for the full accounting.",
    sep = " "
  )

  reference <- paste(
    "Zhou L, Gan J, Yoshitsugu H, Gu X, Lutz JD, Masson E, Humphreys WG.",
    "Integration of Physiologically-Based Pharmacokinetic Modeling into Early",
    "Clinical Development: An Investigation of the Pharmacokinetic",
    "Nonlinearity.",
    "CPT Pharmacometrics Syst Pharmacol. 2015;4(5):286-294.",
    "doi:10.1002/psp4.35. PMCID: PMC4452934.",
    "Vmax = 48.1 pmol/min/mg protein and Km = 1.9 uM are the annotations on the",
    "'HLM kinetics' panel of Figure 1. Incubation design (0.25 mg/mL human",
    "liver microsomes, 1 mM NADPH, pH 7.4 phosphate buffer, 37 C, 10 min,",
    "substrate 0.1-10 uM, triplicate, nonlinear fit in GraphPad Prism):",
    "Methods, 'BMS-911543 metabolism'. The fraction unbound in incubation",
    "(fumic = 0.78, predicted rather than measured) and the metabolic scaling",
    "factor of 4 applied to Km inside the Simcyp model: Methods, 'PBPK",
    "modeling and simulation'. The statement that M1 was the only drug-related",
    "component detected in pooled HLM and in cDNA-expressed CYP enzymes, and",
    "that CYP1A2 is the primary enzyme with CYP3A4 and CYP2J2 minor:",
    "Results, 'BMS-911543 metabolism'.",
    sep = " "
  )

  vignette <- "Zhou_2015_bms911543_invitro"

  units <- list(
    time = "min",
    dosing = "uM (incubation concentration)",
    concentration = "uM"
  )

  # The "dose" of this static in vitro system is the BMS-911543 concentration
  # spiked into the incubation at time zero, so the dosing target is neither
  # `depot` nor `central`. Same convention as Hyland_2008_maraviroc_hlm.R and,
  # upstream of it, HernandezLozano_2025_apramycin_invitro.R.
  dosing <- c("bms911543")

  # A microsomal incubation has no body compartments, so each chemical species
  # in the incubation is a state named after the species itself. `m1` is the
  # metabolite the paper names only as M1 (structure in Supplementary Figure
  # S1, which is not on disk); it is not a registered metabolite suffix, so
  # the parent + metabolite `<canonical>_<metab>` scheme in
  # compartment-names.md does not apply.
  paper_specific_compartments <- c("bms911543", "m1")

  compartmentData <- list(
    bms911543 = list(
      analyte = "BMS-911543 (a selective small-molecule JAK2 inhibitor; no INN assigned)",
      units = "uM",
      specimen = "administration site",
      verified = TRUE
    ),
    m1 = list(
      analyte = "M1, the primary oxidative metabolite of BMS-911543 and the only major metabolite identified in vitro and in animals",
      units = "uM",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "in vitro (pooled human liver microsomes)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    system = "Pooled human liver microsomes at 0.25 mg/mL, supplemented with 1 mM reduced nicotinamide adenine dinucleotide phosphate in pH 7.4 phosphate buffer",
    temperature = "37 C",
    kinetic_incubation = "10 min; the reaction was terminated with an equal volume of acetonitrile and the supernatant analysed by LC/MS/MS",
    concentration_range = "0.1 to 10 uM BMS-911543",
    replication = "All experiments were performed in triplicate; Km and Vmax were determined by nonlinear fitting in GraphPad Prism and are presented as mean and SE, but the SEs on Km and Vmax are not printed",
    disease_state = "not applicable (in vitro)",
    notes = paste(
      "This is the paper's own bench work and is the reason the extraction exists: the Simcyp V12 whole-body PBPK model built on top of these constants is a vendor platform model whose physiology is not published, and is recorded in the vignette rather than extracted.",
      "Preliminary reaction phenotyping (data not shown in the source) established that CYP1A2, CYP3A4 and CYP2J2 are the only CYP enzymes capable of producing M1. The kinetic data in Figure 1, together with a CYP-inhibition study in HLM using specific chemical inhibitors (also data not shown), identified CYP1A2 as the primary enzyme, with CYP3A4 and CYP2J2 playing a minor role. The Simcyp model predicted CYP1A2 to carry 96 percent of total clearance.",
      "Studies in bile-duct cannulated rats showed direct excretion of BMS-911543 in bile and urine to be minimal, under 2 percent of total clearance, so metabolic clearance dominates.",
      "Three further metabolites (M2, M3, M4) were later found by profiling human plasma samples, and all were judged to be downstream products formed from M1; no kinetic constants are reported for them and they are not carried in this model.",
      "The apparent Km of 1.9 uM is close to the top of the 0.1-10 uM substrate range in units of Km multiples (the highest concentration is 5.3 x Km), so the plateau of the curve is reasonably well determined; the fitted velocity at 10 uM is 40.4 pmol/min/mg against a Vmax of 48.1."
    )
  )

  ini({
    # =====================================================================
    # Michaelis-Menten characterisation of M1 formation in pooled HLM.
    # Figure 1, 'HLM kinetics' panel annotation.
    # The panel prints the two point estimates only; the Methods state the
    # fit is presented as mean and SE, but no SE is printed for either
    # constant, so no uncertainty is encoded.
    # =====================================================================
    km_hlm <- 1.9
    label("Michaelis constant for M1 formation in pooled human liver microsomes (uM)") # Figure 1, 'HLM kinetics' panel: Km = 1.9 uM

    vmax_hlm <- 48.1
    label("Maximum velocity of M1 formation in pooled human liver microsomes (pmol/min/mg protein)") # Figure 1, 'HLM kinetics' panel: Vmax = 48.1 pmol/min/mg protein

    # =====================================================================
    # Measured / declared properties of the incubation.
    # =====================================================================
    prot_inc <- fixed(0.25)
    label("Microsomal protein concentration of the kinetic incubation (mg/mL)") # Methods, 'BMS-911543 metabolism': HLMs 0.25 mg/mL

    # fumic was PREDICTED by the Simcyp built-in method, not measured. It is
    # carried fixed so the unbound Michaelis constant derived in model() is
    # visible and auditable, and so a user who prefers a measured value can
    # override it. Its provenance is flagged in the vignette Errata.
    fumic <- fixed(0.78)
    label("Fraction of BMS-911543 unbound in the incubation at 0.25 mg/mL microsomal protein, predicted by the Simcyp built-in method rather than measured (unitless)") # Methods, 'PBPK modeling and simulation'

    # =====================================================================
    # Residual error is NOT reported. The source fitted the concentration-
    # velocity curve in GraphPad Prism and reports only the point estimates;
    # there is no residual-error model, no assay CV and no goodness-of-fit
    # statistic anywhere in the paper. Per the standing policy on unreported
    # residual error the term is fixed at zero so the model returns the
    # deterministic published curve. Flagged in the vignette Errata.
    # =====================================================================
    addSd <- fixed(0)
    label("Additive residual SD of the M1 formation velocity, ZERO because the source reports no residual-error model (pmol/min/mg protein)")
  })

  model({
    # ===================================================================
    # 1. Published Michaelis-Menten velocity. This is the quantity plotted
    #    against substrate concentration in the 'HLM kinetics' panel of
    #    Figure 1, in pmol M1 formed per min per mg microsomal protein.
    #    Under the initial-rate conditions of the assay (10 min) the
    #    substrate is essentially undepleted, so vM1 at t = 0 reproduces
    #    the figure; the vignette asserts that the depletion over a 10 min
    #    incubation stays inside initial-rate conditions at every
    #    concentration the authors used.
    # ===================================================================
    vM1 <- vmax_hlm * bms911543 / (km_hlm + bms911543)

    # ===================================================================
    # 2. Derived intrinsic clearances. The low-substrate limit of the
    #    Michaelis-Menten route is Vmax/Km, which with Vmax in
    #    pmol/min/mg and Km in uM (= pmol/uL) is a clearance in uL/min/mg
    #    microsomal protein: 48.1 / 1.9 = 25.3 uL/min/mg. Correcting the
    #    apparent Km for microsomal binding gives the unbound Michaelis
    #    constant and the unbound intrinsic clearance. Neither correction
    #    is performed by the authors on the Figure 1 numbers; both are
    #    reproduced here as derived, auditable quantities.
    # ===================================================================
    clint_hlm <- vmax_hlm / km_hlm
    km_u_hlm <- km_hlm * fumic
    clint_u_hlm <- vmax_hlm / km_u_hlm

    # ===================================================================
    # 3. Volumetric rate and the incubation ODE. vmax_hlm * prot_inc has
    #    units pmol/(mL*min) = nmol/L/min, so the factor 1/1000 carries it
    #    to uM/min. Getting this bridge wrong misstates every rate by
    #    1000-fold and nothing else in the model catches it.
    # ===================================================================
    rate_m1 <- vmax_hlm * bms911543 / (km_hlm + bms911543) * prot_inc / 1000

    d/dt(bms911543) <- -rate_m1
    d/dt(m1) <- rate_m1

    # ===================================================================
    # 4. Observations. The fitted endpoint is the formation velocity of
    #    Figure 1. The substrate and metabolite concentrations are
    #    returned as derived quantities without a residual-error model,
    #    because the paper publishes no in-vitro concentration-time data.
    # ===================================================================
    Cbms911543 <- bms911543
    Cm1 <- m1

    vM1 ~ add(addSd)
  })
}
