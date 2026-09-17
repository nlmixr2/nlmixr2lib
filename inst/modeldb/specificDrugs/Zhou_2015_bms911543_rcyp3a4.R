Zhou_2015_bms911543_rcyp3a4 <- function() {
  description <- paste(
    "In vitro (recombinant human CYP3A4). Michaelis-Menten enzyme-kinetic",
    "model of the formation of metabolite M1 from the JAK2 inhibitor",
    "BMS-911543 by heterologously expressed CYP3A4, the least efficient of the",
    "three enzymes that produce M1. The published fit is a plain two-parameter",
    "Michaelis-Menten velocity, v = Vmax * S / (Km + S), with Vmax expressed",
    "per pmol of recombinant enzyme. The recombinant CYP content of the",
    "incubation is NOT reported, so no volumetric depletion ODE can be written",
    "without inventing it; the model is therefore the static velocity",
    "concentration-response that the source actually fitted, driven by the",
    "substrate concentration supplied as the covariate CP_BMS911543_UM.",
    "Siblings: Zhou_2015_bms911543_hlm (pooled human liver microsomes, the one",
    "arm that does carry a depletion ODE), Zhou_2015_bms911543_rcyp1a2 (the",
    "dominant enzyme), Zhou_2015_bms911543_rcyp2j2, and",
    "Zhou_2015_bms911543_cyp1a2_tdi.",
    "The Simcyp V12 whole-body PBPK model that consumes these constants is NOT",
    "part of this file and is not reproducible from the published inputs; the",
    "metabolic scaling factor of 4 the authors applied to Km inside that",
    "platform model is a fitting device and is NOT applied here. See the",
    "validation vignette for the full accounting.",
    sep = " "
  )

  reference <- paste(
    "Zhou L, Gan J, Yoshitsugu H, Gu X, Lutz JD, Masson E, Humphreys WG.",
    "Integration of Physiologically-Based Pharmacokinetic Modeling into Early",
    "Clinical Development: An Investigation of the Pharmacokinetic",
    "Nonlinearity.",
    "CPT Pharmacometrics Syst Pharmacol. 2015;4(5):286-294.",
    "doi:10.1002/psp4.35. PMCID: PMC4452934.",
    "Vmax = 0.048 pmol/min/pmol protein and Km = 1.4 uM are the annotations on",
    "the 'CYP3A4 kinetics' panel of Figure 1. Incubation design (recombinant",
    "CYP3A4, 1 mM NADPH, pH 7.4 phosphate buffer, 37 C, 10 min, substrate",
    "0.1-10 uM, triplicate, nonlinear fit in GraphPad Prism): Methods,",
    "'BMS-911543 metabolism'. The assignment of CYP3A4 as a minor contributor:",
    "Results, 'BMS-911543 metabolism'. The predicted fraction of gut-wall",
    "metabolism (Fg = 1, so intestinal CYP3A4 was not involved in elimination):",
    "Results, 'PBPK modeling and simulation'.",
    sep = " "
  )

  vignette <- "Zhou_2015_bms911543_invitro"

  units <- list(
    time = "min",
    dosing = "(none; static concentration-response model driven by the substrate-concentration covariate CP_BMS911543_UM)",
    concentration = "(the output is an M1 formation velocity in pmol/min/pmol recombinant CYP3A4; the driving covariate CP_BMS911543_UM is the BMS-911543 concentration in the incubation in uM)"
  )

  covariateData <- list(
    CP_BMS911543_UM = list(
      description = "Concentration of BMS-911543 in the recombinant-enzyme incubation, supplied as a covariate. Reused canonical: in this in-vitro model the column carries an incubation-medium substrate concentration rather than a plasma concentration, but the quantity and units (uM) are identical. Same reuse rationale as CP_RIF_UM in Almond_2016_rifampicin_invitro.R, which carries a hepatocyte culture-medium concentration.",
      units = "umol/L (uM)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Zhou 2015 Methods, 'BMS-911543 metabolism': the substrate was incubated over 0.1-10 uM. The individual concentrations are not tabulated; the points plotted in Figure 1 sit at approximately 0.1, 0.2, 0.5, 1, 2, 5 and 10 uM.",
        "The top concentration is 7.1-fold above the fitted Km of 1.4 uM, so the plateau of the CYP3A4 curve is the best determined of the three recombinant systems; the fitted velocity at 10 uM is 0.0421 against a Vmax of 0.048 pmol/min/pmol protein, i.e. 88 percent of maximal.",
        "Set to 0 for a no-substrate condition, at which the model returns a velocity of 0 by construction."
      ),
      source_name = "substrate concentration"
    )
  )

  population <- list(
    species = "in vitro (recombinant human CYP3A4)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    system = "cDNA-expressed recombinant human CYP3A4, supplemented with 1 mM reduced nicotinamide adenine dinucleotide phosphate in pH 7.4 phosphate buffer. The recombinant enzyme concentration of the incubation is NOT reported; the 0.25 mg/mL protein concentration given in the same Methods sentence applies to the human liver microsome arm.",
    temperature = "37 C",
    kinetic_incubation = "10 min; the reaction was terminated with an equal volume of acetonitrile and the supernatant analysed by LC/MS/MS",
    concentration_range = "0.1 to 10 uM BMS-911543",
    replication = "All experiments were performed in triplicate; Km and Vmax were determined by nonlinear fitting in GraphPad Prism and are presented as mean and SE, but the SEs on Km and Vmax are not printed",
    disease_state = "not applicable (in vitro)",
    notes = paste(
      "CYP3A4 is the weakest of the three M1-forming enzymes on a per-pmol-enzyme basis: its Vmax/Km of 0.034 uL/min/pmol is 44-fold below CYP1A2 and 11-fold below CYP2J2. Its Vmax of 0.048 pmol/min/pmol is 67-fold below CYP1A2's, while its Km of 1.4 uM is the lowest of the four systems, so the difference between the isoforms is almost entirely in capacity rather than affinity.",
      "The Simcyp model separately predicted the fraction of drug escaping gut-wall metabolism to be unity, indicating intestinal CYP3A4 played no part in the elimination of BMS-911543; that prediction is a platform output and is recorded in the vignette rather than encoded here.",
      "Only the points at 2 and 10 uM in the 'CYP3A4 kinetics' panel of Figure 1 carry visible error bars, which is consistent with the small absolute velocities being close to the limit of the assay."
    )
  )

  ini({
    # =====================================================================
    # Michaelis-Menten characterisation of M1 formation by recombinant
    # CYP3A4. Figure 1, 'CYP3A4 kinetics' panel annotation.
    # The panel prints the two point estimates only; the Methods state the
    # fit is presented as mean and SE, but no SE is printed for either
    # constant, so no uncertainty is encoded.
    # =====================================================================
    km_cyp3a4 <- 1.4
    label("Michaelis constant for M1 formation by recombinant CYP3A4 (uM)") # Figure 1, 'CYP3A4 kinetics' panel: Km = 1.4 uM

    vmax_cyp3a4 <- 0.048
    label("Maximum velocity of M1 formation by recombinant CYP3A4 (pmol/min/pmol protein)") # Figure 1, 'CYP3A4 kinetics' panel: Vmax = 0.048 pmol/min/pmol protein

    # =====================================================================
    # Residual error is NOT reported; fixed at zero per the standing policy
    # on unreported residual error, as in the sibling files. Flagged in the
    # vignette Errata.
    # =====================================================================
    addSd <- fixed(0)
    label("Additive residual SD of the M1 formation velocity, ZERO because the source reports no residual-error model (pmol/min/pmol protein)")
  })

  model({
    # ===================================================================
    # Published Michaelis-Menten velocity, in pmol M1 formed per min per
    # pmol recombinant CYP3A4. This is exactly the curve drawn through the
    # points of the 'CYP3A4 kinetics' panel of Figure 1.
    # ===================================================================
    vM1 <- vmax_cyp3a4 * CP_BMS911543_UM / (km_cyp3a4 + CP_BMS911543_UM)

    # ===================================================================
    # Derived intrinsic clearance, Vmax/Km in uL/min/pmol enzyme:
    # 0.048 / 1.4 = 0.0343 uL/min/pmol CYP3A4.
    # ===================================================================
    clint_cyp3a4 <- vmax_cyp3a4 / km_cyp3a4

    vM1 ~ add(addSd)
  })
}
