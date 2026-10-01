He_2022_schizandrolA_cyp3a <- function() {
  description <- paste(
    "In vitro (CYP3A5-genotyped human liver microsomes). Inhibition of CYP3A4",
    "and CYP3A5 by schizandrol A (SZA), a lignan of the Wuzhi capsule",
    "(Schisandra sphenanthera extract), which is co-prescribed with tacrolimus",
    "(FK-506) in China. SZA is the weaker of the two schizandrol lignans and",
    "acts only by time-dependent (mechanism-based) inactivation; it shows",
    "little reversible inhibition of CYP3A (Figure 2A), so NO reversible",
    "constant is carried.",
    "(1) Time-dependent inactivation of CYP3A4 (in CYP3A5*3/*3 microsomes),",
    "measured with testosterone 6-beta-hydroxylation as the probe; the observed",
    "first-order inactivation rate constant is kobs = kinact * I / (KI + I)",
    "(double-reciprocal plot, Figure 4B), encoded as",
    "d/dt(enzyme_3a4) = -kobs * enzyme_3a4 with enzyme_3a4(0) = 1.",
    "(2) SZA produced NO time-dependent inactivation of CYP3A5 (Figure 4C), so",
    "the CYP3A5 activity state enzyme_3a5(0) = 1 has no inactivation term.",
    "The pooled-HLM CYP3A screen (kinact 0.029 /min, KI 15.625 uM, Figure 3B)",
    "is carried as ini() parameters and exposed as a derived pooled kobs.",
    "The inhibitor concentration I is supplied as the covariate CP_SZA_UM.",
    "Sibling from the same paper: He_2022_schizandrolB_cyp3a (the other Wuzhi",
    "schizandrol, the more potent inhibitor). The whole-body Simcyp DDI model",
    "that consumes these constants is a vendor platform model and is NOT part",
    "of this file; see the vignette.",
    sep = " "
  )

  reference <- paste(
    "He Q, Bu F, Wang Q, Li M, Lin J, Tang Z, Mak WY, Zhuang X, Zhu X,",
    "Lin HS, Xiang X. Examination of the Impact of CYP3A4/5 on Drug-Drug",
    "Interaction between Schizandrol A/Schizandrol B and Tacrolimus (FK-506):",
    "A Physiologically Based Pharmacokinetic Modeling Approach.",
    "Int J Mol Sci. 2022;23(9):4485. doi:10.3390/ijms23094485.",
    "PMCID: PMC9103789.",
    "Little reversible inhibition of CYP3A by SZA (Figure 2A, Results 2.2).",
    "Time-dependent inactivation constants: Results 2.3, Discussion, and the",
    "Figure 4B annotation (CYP3A4 kinact = 0.024 /min, KI = 15.38 uM); no TDI",
    "on CYP3A5 (Figure 4C). Pooled-HLM CYP3A screen (kinact = 0.029 /min,",
    "KI = 15.625 uM): Results 2.3 and Figure 3B. IC50 shift 1.57 (IC50 63.46",
    "uM no preincubation, 40.45 uM after NADPH preincubation): Results 2.1.",
    "Michaelis-Menten inactivation form kobs = kinact*I/(KI+I) and its",
    "estimation by the double-reciprocal plot: Materials and Methods 4.4 and",
    "Equation (3). Assay design: Methods 4.3 (RI) and 4.4 (TDI).",
    sep = " "
  )

  vignette <- "He_2022_tacrolimus_wuzhi_schizandrol_cyp3a"

  units <- list(
    time = "min",
    dosing = "(none; the inhibitor concentration is held constant through the preincubation and is supplied as the covariate CP_SZA_UM)",
    concentration = "(the states and outputs enzyme_3a4 and enzyme_3a5 are CYP catalytic activity as a fraction of the no-inhibitor control, dimensionless; the driving covariate CP_SZA_UM is the schizandrol A concentration in the incubation, in uM)"
  )

  compartmentData <- list(
    enzyme_3a4 = list(
      analyte = "Catalytically active CYP3A4, expressed as a fraction of the no-inhibitor control activity and read out as the rate of testosterone 6-beta-hydroxylation",
      units = "fraction of control (dimensionless)",
      specimen = "not applicable",
      verified = TRUE
    ),
    enzyme_3a5 = list(
      analyte = "Catalytically active CYP3A5, expressed as a fraction of the no-inhibitor control activity; SZA does not inactivate CYP3A5 (Figure 4C) so this state stays at its baseline of 1",
      units = "fraction of control (dimensionless)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    CP_SZA_UM = list(
      description = "Concentration of schizandrol A (SZA) in the incubation, supplied as a covariate. This is the quantity the source calls I, the inhibitor concentration, in the mechanism-based inactivation-rate expression. In this in-vitro model the column carries an incubation-buffer concentration rather than a plasma concentration, but the quantity and units (uM) are identical; same reuse rationale as CP_STA_UM in He_2021_schisantherinA_cyp3a.R.",
      units = "umol/L (uM)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "TDI assay, pooled CYP3A (Methods 4.4): SZA at 0, 10, 16, 20, 25, 32 uM.",
        "TDI assay, CYP3A4 (Methods 4.4): SZA at 0, 10, 16, 20, 25, 32 uM; CYP3A5: 0, 10, 20, 30, 40, 50 uM.",
        "Held CONSTANT over the preincubation by the assay design (no depletion correction), so a static covariate value is the faithful encoding.",
        "Set to 0 for the no-inhibitor control, at which kobs = 0 and enzyme_3a4 stays at 1 for all time."
      ),
      source_name = "I (inhibitor concentration)"
    )
  )

  population <- list(
    species = "in vitro (human liver microsomes)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    system = "Pooled human liver microsomes (22 donors) for the CYP3A screen; CYP3A5*3/*3 human liver microsomes (for CYP3A4) and CYP3A5*1/*3 human liver microsomes co-incubated with the selective CYP3A4 inactivator CYP3cide (for CYP3A5). RI assay at 0.2 mg/mL microsomal protein; TDI assay at 0.5 mg/mL microsomal protein, in 0.1 M potassium phosphate buffer pH 7.4.",
    temperature = "37 C",
    probe_reaction = "Tacrolimus disappearance (Dixon plots, reversible inhibition) and testosterone 6-beta-hydroxylation (time-dependent inactivation)",
    disease_state = "not applicable (in vitro)",
    notes = paste(
      "The reported +/- terms on the TDI constants are standard errors of the nonlinear regression, not between-subject variability, so no omega is encoded.",
      "SZA is a weak, purely time-dependent inhibitor of CYP3A4 with no reversible inhibition (Figure 2A) and no CYP3A5 inactivation (Figure 4C).",
      "kinact is reported in per-minute units (Figures 3B, 4B), consistent with the per-minute time axis of the assay; the model's time unit is minutes.",
      "The CYP degradation rate constant kdeg, which governs the in-vivo magnitude of TDI, is a Simcyp default and is not printed; no enzyme-turnover term is carried, so the model describes the inactivation phase of the in-vitro assay only."
    )
  )

  ini({
    # =====================================================================
    # Time-dependent (mechanism-based) inactivation of CYP3A4, from the
    # double-reciprocal (kobs vs 1/[SZA]) plot of Figure 4B.
    # SZA showed little reversible inhibition of CYP3A (Figure 2A), so no
    # reversible constant is carried, and NO time-dependent inactivation of
    # CYP3A5 (Figure 4C), so there is no CYP3A5 inactivation term.
    # No uncertainty is reported, so the constants are point values.
    # =====================================================================
    ki_inact_3a4 <- 15.38
    label("Schizandrol A concentration at half-maximal CYP3A4 inactivation rate KI (uM)") # Figure 4B annotation and Results 2.3 / Discussion: 'KI = 15.38 uM'

    kinact_3a4 <- 0.024
    label("Maximum rate constant of CYP3A4 inactivation by schizandrol A kinact (1/min)") # Figure 4B annotation and Results 2.3 / Discussion: 'kinact = 0.024 min-1'

    # =====================================================================
    # Pooled-HLM CYP3A screen (not isoform-resolved), carried so the
    # paper's initial aggregate characterisation is auditable alongside the
    # CYP3A4-specific values the DDI model consumes. Figure 3B and the
    # Discussion: kinact = 0.029 /min, KI = 15.625 uM.
    # =====================================================================
    ki_inact_3a <- 15.625
    label("Schizandrol A concentration at half-maximal pooled-CYP3A inactivation rate KI (uM)") # Figure 3B annotation and Results 2.3 / Discussion: 'KI ... 15.625 uM'

    kinact_3a <- 0.029
    label("Maximum rate constant of pooled-CYP3A inactivation by schizandrol A kinact (1/min)") # Figure 3B annotation and Results 2.3 / Discussion: 'kinact ... 0.029 min-1'

    # =====================================================================
    # Residual error is NOT reported (the source fits double-reciprocal
    # plots and reports only the resulting constants). Per the standing
    # policy on unreported residual error the term is fixed at zero so the
    # model returns the deterministic published curve. Vignette Errata.
    # =====================================================================
    addSd <- fixed(0)
    label("Additive residual SD of the relative CYP3A4 activity, ZERO because the source reports no residual-error model (fraction of control)")
  })

  model({
    # ===================================================================
    # Pooled-CYP3A observed inactivation rate (Figure 3B), exposed as a
    # derived output. This is the paper's aggregate screen; the CYP3A4
    # inactivation below is the isoform-resolved mechanism.
    # ===================================================================
    kobs_3a <- kinact_3a * CP_SZA_UM / (ki_inact_3a + CP_SZA_UM)

    # ===================================================================
    # Time-dependent inactivation of CYP3A4. The observed first-order rate
    # constant is a hyperbolic (Michaelis-Menten) function of the SZA
    # concentration; its plateau is kinact and it passes through kinact/2
    # at I = KI. enzyme_3a4 is the registered CYP3A4 activity state as a
    # fraction of its untreated baseline, enzyme_3a4(0) = 1.
    # ===================================================================
    kobs_3a4 <- kinact_3a4 * CP_SZA_UM / (ki_inact_3a4 + CP_SZA_UM)

    enzyme_3a4(0) <- 1
    d/dt(enzyme_3a4) <- -kobs_3a4 * enzyme_3a4

    # CYP3A5 has no time-dependent inactivation by SZA (Figure 4C): its
    # activity state holds at baseline of 1.
    enzyme_3a5(0) <- 1
    d/dt(enzyme_3a5) <- 0

    # Figures 3A/4A plot the log of percent activity remaining; exposed as
    # a derived output so the panel can be replicated directly.
    pctActivity_3a4 <- 100 * enzyme_3a4

    enzyme_3a4 ~ add(addSd)
  })
}
