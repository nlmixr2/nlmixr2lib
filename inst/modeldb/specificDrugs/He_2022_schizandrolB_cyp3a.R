He_2022_schizandrolB_cyp3a <- function() {
  description <- paste(
    "In vitro (CYP3A5-genotyped human liver microsomes). Inhibition of CYP3A4",
    "and CYP3A5 by schizandrol B (SZB), a lignan of the Wuzhi capsule",
    "(Schisandra sphenanthera extract), which is co-prescribed with tacrolimus",
    "(FK-506) in China. SZB is the more potent of the two schizandrol lignans",
    "and displays BOTH reversible (RI) and time-dependent (TDI) inhibition of",
    "CYP3A4 and CYP3A5.",
    "(1) Competitive reversible inhibition of both isoforms, characterised by",
    "an inhibition constant Ki from Dixon plots (CYP3A4 Ki = 2.18 uM, Figure",
    "2C; CYP3A5 Ki = 2.03 uM, Figure 2D). The constants scale a victim-drug",
    "clearance as 1/(1 + I/Ki) and are exposed as derived reversible factors.",
    "(2) Time-dependent inactivation of both isoforms (double-reciprocal plots",
    "Figure 5B CYP3A4, Figure 5D CYP3A5): kobs = kinact * I / (KI + I),",
    "encoded as d/dt(enzyme_3a4) = -kobs_3a4 * enzyme_3a4 and",
    "d/dt(enzyme_3a5) = -kobs_3a5 * enzyme_3a5 with both states starting at 1.",
    "The pooled-HLM CYP3A screen (reversible Ki 5.82 uM; TDI kinact 0.044 /min,",
    "KI 0.43 uM) is carried as ini() parameters and exposed as derived pooled",
    "factors. The inhibitor concentration I is supplied as CP_SZB_UM.",
    "Sibling from the same paper: He_2022_schizandrolA_cyp3a (the weaker Wuzhi",
    "schizandrol). The whole-body Simcyp DDI model that consumes these",
    "constants is a vendor platform model and is NOT part of this file; see the",
    "vignette.",
    sep = " "
  )

  reference <- paste(
    "He Q, Bu F, Wang Q, Li M, Lin J, Tang Z, Mak WY, Zhuang X, Zhu X,",
    "Lin HS, Xiang X. Examination of the Impact of CYP3A4/5 on Drug-Drug",
    "Interaction between Schizandrol A/Schizandrol B and Tacrolimus (FK-506):",
    "A Physiologically Based Pharmacokinetic Modeling Approach.",
    "Int J Mol Sci. 2022;23(9):4485. doi:10.3390/ijms23094485.",
    "PMCID: PMC9103789.",
    "Reversible inhibition constants (Ki): Results 2.2 and Figure 2C,D",
    "(SZB on CYP3A4 Ki = 2.18 uM, on CYP3A5 Ki = 2.03 uM); pooled CYP3A",
    "Ki = 5.82 uM (Figure 2B). NOTE: the Discussion text swaps the CYP3A4 and",
    "CYP3A5 Ki values (it prints CYP3A4 Ki = 2.03, CYP3A5 Ki = 2.18); Results",
    "2.2 and the Figure 2C,D annotations are used here (CYP3A4 = 2.18,",
    "CYP3A5 = 2.03). Time-dependent inactivation constants: Results 2.3,",
    "Discussion, and the figure annotations (CYP3A4 kinact = 0.37 /min,",
    "KI = 0.69 uM, Figure 5B; CYP3A5 kinact = 0.009 /min, KI = 0.5 uM, Figure",
    "5D; pooled CYP3A kinact = 0.044 /min, KI = 0.43 uM, Figure 3D). IC50 shift",
    "21.39 (IC50 11.98 uM no preincubation, 0.56 uM after NADPH preincubation):",
    "Results 2.1. Michaelis-Menten inactivation form kobs = kinact*I/(KI+I):",
    "Materials and Methods 4.4 and Equation (3).",
    sep = " "
  )

  vignette <- "He_2022_tacrolimus_wuzhi_schizandrol_cyp3a"

  units <- list(
    time = "min",
    dosing = "(none; the inhibitor concentration is held constant through the preincubation and is supplied as the covariate CP_SZB_UM)",
    concentration = "(the states and outputs enzyme_3a4 and enzyme_3a5 are CYP catalytic activity as a fraction of the no-inhibitor control, dimensionless; the driving covariate CP_SZB_UM is the schizandrol B concentration in the incubation, in uM)"
  )

  compartmentData <- list(
    enzyme_3a4 = list(
      analyte = "Catalytically active CYP3A4, expressed as a fraction of the no-inhibitor control activity and read out as the rate of testosterone 6-beta-hydroxylation",
      units = "fraction of control (dimensionless)",
      specimen = "not applicable",
      verified = TRUE
    ),
    enzyme_3a5 = list(
      analyte = "Catalytically active CYP3A5, expressed as a fraction of the no-inhibitor control activity and read out as the rate of testosterone 6-beta-hydroxylation",
      units = "fraction of control (dimensionless)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    CP_SZB_UM = list(
      description = "Concentration of schizandrol B (SZB) in the incubation, supplied as a covariate. This is the quantity the source calls I, the inhibitor concentration, in both the reversible-inhibition term 1/(1 + CP_SZB_UM/ki) and the mechanism-based inactivation-rate term kinact * CP_SZB_UM / (KI + CP_SZB_UM). In this in-vitro model the column carries an incubation-buffer concentration rather than a plasma concentration, but the quantity and units (uM) are identical; same reuse rationale as CP_STA_UM in He_2021_schisantherinA_cyp3a.R.",
      units = "umol/L (uM)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "RI assay, pooled CYP3A (Methods 4.3): SZB at 0, 1, 2, 4 uM; genotyped CYP3A4/CYP3A5: 0, 1, 2, 4 uM.",
        "TDI assay, pooled CYP3A (Methods 4.4): SZB at 0, 0.1, 0.2, 0.5, 1, 2 uM; CYP3A4: 0, 0.5, 1, 2, 4, 8 uM; CYP3A5: 0, 2, 4, 8, 16 uM.",
        "Held CONSTANT over the preincubation by the assay design (no depletion correction), so a static covariate value is the faithful encoding.",
        "Set to 0 for the no-inhibitor control, at which both kobs values are 0 and the activity states stay at 1 for all time and the reversible factors are 1."
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
      "SZB inhibits BOTH CYP3A4 and CYP3A5 reversibly AND in a time-dependent manner, and is the more potent of the two schizandrol lignans (kinact/KI far exceeds SZA's).",
      "The Discussion text swaps the reversible CYP3A4 and CYP3A5 Ki values relative to Results 2.2 and Figure 2C,D; the figure/Results order is used (CYP3A4 = 2.18 uM, CYP3A5 = 2.03 uM). See the vignette Errata.",
      "kinact is reported in per-minute units (Figures 3D, 5B, 5D), consistent with the per-minute time axis of the assay; the model's time unit is minutes.",
      "The CYP degradation rate constant kdeg, which governs the in-vivo magnitude of TDI, is a Simcyp default and is not printed; no enzyme-turnover term is carried, so the model describes the inactivation phase of the in-vitro assay only."
    )
  )

  ini({
    # =====================================================================
    # Reversible (competitive) inhibition constants, from the Dixon plots.
    # SZB on CYP3A4: Ki = 2.18 uM (Figure 2C); on CYP3A5: Ki = 2.03 uM
    # (Figure 2D). Results 2.2 gives '(2.18 and 2.03 uM, respectively)' for
    # CYP3A4 and CYP3A5 and matches the figure annotations; the Discussion
    # swaps the two, so the figure/Results order is used here.
    # No uncertainty is reported, so they are carried as point values.
    # =====================================================================
    ki_3a4 <- 2.18
    label("Schizandrol B equilibrium dissociation constant Ki for reversible CYP3A4 inhibition (uM)") # Figure 2C annotation and Results 2.2: CYP3A4 Ki = 2.18 uM

    ki_3a5 <- 2.03
    label("Schizandrol B equilibrium dissociation constant Ki for reversible CYP3A5 inhibition (uM)") # Figure 2D annotation and Results 2.2: CYP3A5 Ki = 2.03 uM

    # Pooled-HLM CYP3A reversible constant (not isoform-resolved).
    ki_3a <- 5.82
    label("Schizandrol B equilibrium dissociation constant Ki for reversible pooled-CYP3A inhibition (uM)") # Figure 2B annotation and Results 2.2: 'a Ki value of 5.82 uM'

    # =====================================================================
    # Time-dependent (mechanism-based) inactivation constants, from the
    # double-reciprocal (kobs vs 1/[SZB]) plots. CYP3A4: kinact = 0.37 /min,
    # KI = 0.69 uM (Figure 5B). CYP3A5: kinact = 0.009 /min, KI = 0.5 uM
    # (Figure 5D).
    # =====================================================================
    ki_inact_3a4 <- 0.69
    label("Schizandrol B concentration at half-maximal CYP3A4 inactivation rate KI (uM)") # Figure 5B annotation and Results 2.3 / Discussion: 'KI = 0.69 uM'

    kinact_3a4 <- 0.37
    label("Maximum rate constant of CYP3A4 inactivation by schizandrol B kinact (1/min)") # Figure 5B annotation and Results 2.3 / Discussion: 'kinact = 0.37 min-1'

    ki_inact_3a5 <- 0.5
    label("Schizandrol B concentration at half-maximal CYP3A5 inactivation rate KI (uM)") # Figure 5D annotation and Results 2.3 / Discussion: 'KI = 0.5 uM'

    kinact_3a5 <- 0.009
    label("Maximum rate constant of CYP3A5 inactivation by schizandrol B kinact (1/min)") # Figure 5D annotation and Results 2.3 / Discussion: 'kinact = 0.009 min-1'

    # =====================================================================
    # Pooled-HLM CYP3A time-dependent inactivation screen (Figure 3D).
    # =====================================================================
    ki_inact_3a <- 0.43
    label("Schizandrol B concentration at half-maximal pooled-CYP3A inactivation rate KI (uM)") # Figure 3D annotation and Results 2.3 / Discussion: 'KI = 0.43 uM'

    kinact_3a <- 0.044
    label("Maximum rate constant of pooled-CYP3A inactivation by schizandrol B kinact (1/min)") # Figure 3D annotation and Results 2.3 / Discussion: 'kinact = 0.044 min-1'

    # =====================================================================
    # Residual error is NOT reported. Per the standing policy on unreported
    # residual error the term is fixed at zero so the model returns the
    # deterministic published curve. Flagged in the vignette Errata.
    # =====================================================================
    addSd <- fixed(0)
    label("Additive residual SD of the relative CYP3A4 activity, ZERO because the source reports no residual-error model (fraction of control)")
  })

  model({
    # ===================================================================
    # Reversible-inhibition factors. A victim-drug intrinsic clearance by
    # each isoform is scaled by 1/(1 + I/Ki) under competitive inhibition;
    # exposed as derived outputs (plus the pooled-CYP3A factor).
    # ===================================================================
    riFactor_3a4 <- 1 / (1 + CP_SZB_UM / ki_3a4)
    riFactor_3a5 <- 1 / (1 + CP_SZB_UM / ki_3a5)
    riFactor_3a <- 1 / (1 + CP_SZB_UM / ki_3a)

    # ===================================================================
    # Pooled-CYP3A observed inactivation rate (Figure 3D), derived output.
    # ===================================================================
    kobs_3a <- kinact_3a * CP_SZB_UM / (ki_inact_3a + CP_SZB_UM)

    # ===================================================================
    # Time-dependent inactivation of CYP3A4 and CYP3A5. Each observed
    # first-order rate constant is a hyperbolic (Michaelis-Menten) function
    # of the SZB concentration; each plateau is its kinact and each passes
    # through kinact/2 at I = KI. Both enzyme states are registered CYP
    # activity states as a fraction of untreated baseline, starting at 1.
    # ===================================================================
    kobs_3a4 <- kinact_3a4 * CP_SZB_UM / (ki_inact_3a4 + CP_SZB_UM)
    kobs_3a5 <- kinact_3a5 * CP_SZB_UM / (ki_inact_3a5 + CP_SZB_UM)

    enzyme_3a4(0) <- 1
    d/dt(enzyme_3a4) <- -kobs_3a4 * enzyme_3a4

    enzyme_3a5(0) <- 1
    d/dt(enzyme_3a5) <- -kobs_3a5 * enzyme_3a5

    # Figures 5A/5C plot the log of percent activity remaining; exposed as
    # derived outputs so the panels can be replicated directly.
    pctActivity_3a4 <- 100 * enzyme_3a4
    pctActivity_3a5 <- 100 * enzyme_3a5

    enzyme_3a4 ~ add(addSd)
  })
}
