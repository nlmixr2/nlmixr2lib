Zhou_2015_bms911543_rcyp2j2 <- function() {
  description <- paste(
    "In vitro (recombinant human CYP2J2). Michaelis-Menten enzyme-kinetic",
    "model of the formation of metabolite M1 from the JAK2 inhibitor",
    "BMS-911543 by heterologously expressed CYP2J2, the intermediate of the",
    "three enzymes that produce M1. The published fit is a plain two-parameter",
    "Michaelis-Menten velocity, v = Vmax * S / (Km + S), with Vmax expressed",
    "per pmol of recombinant enzyme. The recombinant CYP content of the",
    "incubation is NOT reported, so no volumetric depletion ODE can be written",
    "without inventing it; the model is therefore the static velocity",
    "concentration-response that the source actually fitted, driven by the",
    "substrate concentration supplied as the covariate CP_BMS911543_UM.",
    "Siblings: Zhou_2015_bms911543_hlm (pooled human liver microsomes, the one",
    "arm that does carry a depletion ODE), Zhou_2015_bms911543_rcyp1a2 (the",
    "dominant enzyme), Zhou_2015_bms911543_rcyp3a4, and",
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
    "Vmax = 0.5 pmol/min/pmol protein and Km = 1.3 uM are the annotations on",
    "the 'CYP2J2 kinetics' panel of Figure 1. Incubation design (recombinant",
    "CYP2J2, 1 mM NADPH, pH 7.4 phosphate buffer, 37 C, 10 min, substrate",
    "0.1-10 uM, triplicate, nonlinear fit in GraphPad Prism): Methods,",
    "'BMS-911543 metabolism'. The assignment of CYP2J2 as a minor contributor:",
    "Results, 'BMS-911543 metabolism'.",
    sep = " "
  )

  vignette <- "Zhou_2015_bms911543_invitro"

  units <- list(
    time = "min",
    dosing = "(none; static concentration-response model driven by the substrate-concentration covariate CP_BMS911543_UM)",
    concentration = "(the output is an M1 formation velocity in pmol/min/pmol recombinant CYP2J2; the driving covariate CP_BMS911543_UM is the BMS-911543 concentration in the incubation in uM)"
  )

  covariateData <- list(
    CP_BMS911543_UM = list(
      description = "Concentration of BMS-911543 in the recombinant-enzyme incubation, supplied as a covariate. Reused canonical: in this in-vitro model the column carries an incubation-medium substrate concentration rather than a plasma concentration, but the quantity and units (uM) are identical. Same reuse rationale as CP_RIF_UM in Almond_2016_rifampicin_invitro.R, which carries a hepatocyte culture-medium concentration.",
      units = "umol/L (uM)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Zhou 2015 Methods, 'BMS-911543 metabolism': the substrate was incubated over 0.1-10 uM. The individual concentrations are not tabulated; the points plotted in Figure 1 sit at approximately 0.1, 0.2, 0.5, 1, 2, 5 and 10 uM.",
        "The top concentration is 7.7-fold above the fitted Km of 1.3 uM; the fitted velocity at 10 uM is 0.442 against a Vmax of 0.5 pmol/min/pmol protein, i.e. 88 percent of maximal.",
        "Set to 0 for a no-substrate condition, at which the model returns a velocity of 0 by construction."
      ),
      source_name = "substrate concentration"
    )
  )

  population <- list(
    species = "in vitro (recombinant human CYP2J2)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    system = "cDNA-expressed recombinant human CYP2J2, supplemented with 1 mM reduced nicotinamide adenine dinucleotide phosphate in pH 7.4 phosphate buffer. The recombinant enzyme concentration of the incubation is NOT reported; the 0.25 mg/mL protein concentration given in the same Methods sentence applies to the human liver microsome arm.",
    temperature = "37 C",
    kinetic_incubation = "10 min; the reaction was terminated with an equal volume of acetonitrile and the supernatant analysed by LC/MS/MS",
    concentration_range = "0.1 to 10 uM BMS-911543",
    replication = "All experiments were performed in triplicate; Km and Vmax were determined by nonlinear fitting in GraphPad Prism and are presented as mean and SE, but the SEs on Km and Vmax are not printed",
    disease_state = "not applicable (in vitro)",
    notes = paste(
      "CYP2J2 sits between CYP1A2 and CYP3A4 on a per-pmol-enzyme basis: its Vmax/Km of 0.385 uL/min/pmol is 4.0-fold below CYP1A2 and 11-fold above CYP3A4.",
      "CYP2J2 is expressed at far lower levels in human liver than either CYP1A2 or CYP3A4, so its per-enzyme efficiency substantially overstates its contribution to hepatic clearance; the abundance needed to make that conversion is not reported in this source, and the Simcyp model's own prediction that CYP1A2 carries 96 percent of total clearance is recorded in the vignette as a platform output rather than as a reproducible check.",
      "CYP2J2 is also expressed in cardiac tissue, but the source draws no conclusion about extrahepatic metabolism of BMS-911543 and none is encoded here."
    )
  )

  ini({
    # =====================================================================
    # Michaelis-Menten characterisation of M1 formation by recombinant
    # CYP2J2. Figure 1, 'CYP2J2 kinetics' panel annotation.
    # The panel prints the two point estimates only; the Methods state the
    # fit is presented as mean and SE, but no SE is printed for either
    # constant, so no uncertainty is encoded.
    # =====================================================================
    km_cyp2j2 <- 1.3
    label("Michaelis constant for M1 formation by recombinant CYP2J2 (uM)") # Figure 1, 'CYP2J2 kinetics' panel: Km = 1.3 uM

    vmax_cyp2j2 <- 0.5
    label("Maximum velocity of M1 formation by recombinant CYP2J2 (pmol/min/pmol protein)") # Figure 1, 'CYP2J2 kinetics' panel: Vmax = 0.5 pmol/min/pmol protein

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
    # pmol recombinant CYP2J2. This is exactly the curve drawn through the
    # points of the 'CYP2J2 kinetics' panel of Figure 1.
    # ===================================================================
    vM1 <- vmax_cyp2j2 * CP_BMS911543_UM / (km_cyp2j2 + CP_BMS911543_UM)

    # ===================================================================
    # Derived intrinsic clearance, Vmax/Km in uL/min/pmol enzyme:
    # 0.5 / 1.3 = 0.3846 uL/min/pmol CYP2J2.
    # ===================================================================
    clint_cyp2j2 <- vmax_cyp2j2 / km_cyp2j2

    vM1 ~ add(addSd)
  })
}
