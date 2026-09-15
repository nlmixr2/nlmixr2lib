Zhou_2015_bms911543_rcyp1a2 <- function() {
  description <- paste(
    "In vitro (recombinant human CYP1A2). Michaelis-Menten enzyme-kinetic",
    "model of the formation of metabolite M1 from the JAK2 inhibitor",
    "BMS-911543 by heterologously expressed CYP1A2, the dominant enzyme of the",
    "three that produce M1. The published fit is a plain two-parameter",
    "Michaelis-Menten velocity, v = Vmax * S / (Km + S), with Vmax expressed",
    "per pmol of recombinant enzyme. The recombinant CYP content of the",
    "incubation is NOT reported, so no volumetric depletion ODE can be written",
    "without inventing it; the model is therefore the static velocity",
    "concentration-response that the source actually fitted, driven by the",
    "substrate concentration supplied as the covariate CP_BMS911543_UM. The",
    "sibling Zhou_2015_bms911543_hlm does carry the depletion ODE, because the",
    "microsomal protein concentration of that incubation is reported.",
    "Siblings: Zhou_2015_bms911543_hlm (pooled human liver microsomes),",
    "Zhou_2015_bms911543_rcyp3a4 and Zhou_2015_bms911543_rcyp2j2 (the two minor",
    "enzymes), and Zhou_2015_bms911543_cyp1a2_tdi (time-dependent inactivation",
    "of CYP1A2 by the same compound).",
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
    "Vmax = 3.2 pmol/min/pmol protein and Km = 2.1 uM are the annotations on",
    "the 'CYP1A2 kinetics' panel of Figure 1. Incubation design (recombinant",
    "CYP1A2, 1 mM NADPH, pH 7.4 phosphate buffer, 37 C, 10 min, substrate",
    "0.1-10 uM, triplicate, nonlinear fit in GraphPad Prism): Methods,",
    "'BMS-911543 metabolism'. The assignment of CYP1A2 as the primary enzyme,",
    "with CYP3A4 and CYP2J2 minor, and the statement that these three enzymes",
    "are the only ones capable of producing M1: Results,",
    "'BMS-911543 metabolism'. The predicted 96 percent fraction metabolised by",
    "CYP1A2 and the metabolic scaling factor of 4 applied to Km inside the",
    "Simcyp model: Results and Methods, 'PBPK modeling and simulation'.",
    sep = " "
  )

  vignette <- "Zhou_2015_bms911543_invitro"

  units <- list(
    time = "min",
    dosing = "(none; static concentration-response model driven by the substrate-concentration covariate CP_BMS911543_UM)",
    concentration = "(the output is an M1 formation velocity in pmol/min/pmol recombinant CYP1A2; the driving covariate CP_BMS911543_UM is the BMS-911543 concentration in the incubation in uM)"
  )

  covariateData <- list(
    CP_BMS911543_UM = list(
      description = "Concentration of BMS-911543 in the recombinant-enzyme incubation, supplied as a covariate. Reused canonical: in this in-vitro model the column carries an incubation-medium substrate concentration rather than a plasma concentration, but the quantity and units (uM) are identical. Same reuse rationale as CP_RIF_UM in Almond_2016_rifampicin_invitro.R, which carries a hepatocyte culture-medium concentration.",
      units = "umol/L (uM)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Zhou 2015 Methods, 'BMS-911543 metabolism': the substrate was incubated over 0.1-10 uM. The individual concentrations are not tabulated; the points plotted in Figure 1 sit at approximately 0.1, 0.2, 0.5, 1, 2, 5 and 10 uM.",
        "The top concentration is 4.8-fold above the fitted Km of 2.1 uM, so the plateau of the CYP1A2 curve is only moderately well determined; the fitted velocity at 10 uM is 2.64 against a Vmax of 3.2 pmol/min/pmol protein, i.e. 83 percent of maximal.",
        "Set to 0 for a no-substrate condition, at which the model returns a velocity of 0 by construction."
      ),
      source_name = "substrate concentration"
    )
  )

  population <- list(
    species = "in vitro (recombinant human CYP1A2)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    system = "cDNA-expressed recombinant human CYP1A2, supplemented with 1 mM reduced nicotinamide adenine dinucleotide phosphate in pH 7.4 phosphate buffer. The recombinant enzyme concentration of the incubation is NOT reported; the 0.25 mg/mL protein concentration given in the same Methods sentence applies to the human liver microsome arm.",
    temperature = "37 C",
    kinetic_incubation = "10 min; the reaction was terminated with an equal volume of acetonitrile and the supernatant analysed by LC/MS/MS",
    concentration_range = "0.1 to 10 uM BMS-911543",
    replication = "All experiments were performed in triplicate; Km and Vmax were determined by nonlinear fitting in GraphPad Prism and are presented as mean and SE, but the SEs on Km and Vmax are not printed",
    disease_state = "not applicable (in vitro)",
    notes = paste(
      "Metabolite M1 was the only drug-related component detected in cDNA-expressed CYP enzymes and in pooled human liver microsomes, indicating that its formation is the primary biotransformation pathway for BMS-911543.",
      "Preliminary studies (data not shown in the source) established that CYP1A2, CYP3A4 and CYP2J2 are the only CYP enzymes capable of producing M1.",
      "The three recombinant systems are directly comparable on Vmax/Km because they share the per-pmol-enzyme denominator: CYP1A2 1.52, CYP2J2 0.38 and CYP3A4 0.034 uL/min/pmol enzyme, so per unit of enzyme CYP1A2 is 4.0-fold more efficient than CYP2J2 and 44-fold more efficient than CYP3A4. Converting that per-enzyme ranking into a fraction metabolised requires the hepatic abundance of each isoform, which the source does not report; the Simcyp model's predicted 96 percent for CYP1A2 therefore cannot be reproduced from on-disk values and is recorded in the vignette as a platform output rather than a check.",
      "The Km values of the three recombinant enzymes (2.1, 1.4 and 1.3 uM) and of pooled human liver microsomes (1.9 uM) all sit within a 1.6-fold band, consistent with a single dominant binding mode for the reaction."
    )
  )

  ini({
    # =====================================================================
    # Michaelis-Menten characterisation of M1 formation by recombinant
    # CYP1A2. Figure 1, 'CYP1A2 kinetics' panel annotation.
    # The panel prints the two point estimates only; the Methods state the
    # fit is presented as mean and SE, but no SE is printed for either
    # constant, so no uncertainty is encoded.
    # =====================================================================
    km_cyp1a2 <- 2.1
    label("Michaelis constant for M1 formation by recombinant CYP1A2 (uM)") # Figure 1, 'CYP1A2 kinetics' panel: Km = 2.1 uM

    vmax_cyp1a2 <- 3.2
    label("Maximum velocity of M1 formation by recombinant CYP1A2 (pmol/min/pmol protein)") # Figure 1, 'CYP1A2 kinetics' panel: Vmax = 3.2 pmol/min/pmol protein

    # =====================================================================
    # Residual error is NOT reported. The source fitted the concentration-
    # velocity curve in GraphPad Prism and reports only the point estimates;
    # there is no residual-error model, no assay CV and no goodness-of-fit
    # statistic anywhere in the paper. Per the standing policy on unreported
    # residual error the term is fixed at zero so the model returns the
    # deterministic published curve. Flagged in the vignette Errata.
    # =====================================================================
    addSd <- fixed(0)
    label("Additive residual SD of the M1 formation velocity, ZERO because the source reports no residual-error model (pmol/min/pmol protein)")
  })

  model({
    # ===================================================================
    # Published Michaelis-Menten velocity, in pmol M1 formed per min per
    # pmol recombinant CYP1A2. This is exactly the curve drawn through the
    # points of the 'CYP1A2 kinetics' panel of Figure 1.
    #
    # Structural checks asserted in the validation vignette:
    #   vM1(0)         = 0                (no substrate, no product)
    #   vM1(km_cyp1a2) = vmax_cyp1a2 / 2  (definition of the Michaelis constant)
    #   vM1(Inf)       -> vmax_cyp1a2     (the reported maximum velocity)
    # ===================================================================
    vM1 <- vmax_cyp1a2 * CP_BMS911543_UM / (km_cyp1a2 + CP_BMS911543_UM)

    # ===================================================================
    # Derived intrinsic clearance. The low-substrate limit of the
    # Michaelis-Menten route is Vmax/Km, which with Vmax in pmol/min/pmol
    # enzyme and Km in uM (= pmol/uL) is a clearance in uL/min/pmol
    # enzyme: 3.2 / 2.1 = 1.524 uL/min/pmol CYP1A2. Scaling this to a
    # per-mg-liver or per-body basis needs the hepatic abundance of
    # CYP1A2, an intersystem extrapolation factor and MPPGL, none of
    # which the source reports.
    # ===================================================================
    clint_cyp1a2 <- vmax_cyp1a2 / km_cyp1a2

    vM1 ~ add(addSd)
  })
}
