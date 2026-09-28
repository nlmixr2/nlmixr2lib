He_2021_schisantherinA_cyp3a <- function() {
  description <- paste(
    "In vitro (CYP3A5-genotyped human liver microsomes). Inhibition of CYP3A4",
    "and CYP3A5 by schisantherin A (STA), the most abundant lignan of the Wuzhi",
    "capsule (Schisandra sphenanthera extract), which is co-prescribed with",
    "tacrolimus in China. Two independently fitted mechanisms.",
    "(1) Time-dependent (mechanism-based) inactivation of CYP3A4, measured in",
    "CYP3A5*3/*3 microsomes with testosterone 6-beta-hydroxylation as the probe:",
    "the observed first-order inactivation rate constant is",
    "kobs = kinact * I / (KI + I) (double-reciprocal plot, Figure 3B), encoded",
    "as d/dt(enzyme_3a4) = -kobs * enzyme_3a4 with enzyme_3a4(0) = 1.",
    "(2) Competitive reversible inhibition of both isoforms, characterised by an",
    "inhibition constant Ki from Dixon plots (Figure 2A,B); STA showed no",
    "time-dependent inactivation of CYP3A5 (Figure 3C). The reversible-inhibition",
    "constants are carried as ini() parameters (they scale a victim-drug",
    "clearance as 1/(1 + I/Ki) but are not themselves an ODE) and the CYP3A5",
    "activity state enzyme_3a5(0) = 1 has no inactivation term, by the finding.",
    "The inhibitor concentration I is supplied as the covariate CP_STA_UM.",
    "Sibling from the same paper: He_2021_schisandrinA_cyp3a (the other Wuzhi",
    "lignan). The whole-body Simcyp DDI model that consumes these constants is",
    "a vendor platform model and is NOT part of this file; see the vignette.",
    sep = " "
  )

  reference <- paste(
    "He Q, Bu F, Zhang H, Wang Q, Tang Z, Yuan J, Lin HS, Xiang X.",
    "Investigation of the Impact of CYP3A5 Polymorphism on Drug-Drug",
    "Interaction between Tacrolimus and Schisantherin A/Schisandrin A Based on",
    "Physiologically-Based Pharmacokinetic Modeling.",
    "Pharmaceuticals (Basel). 2021;14(3):198. doi:10.3390/ph14030198.",
    "PMCID: PMC7997453.",
    "Reversible inhibition constants (Ki): Results 2.1 and Figure 2A,B",
    "(STA on CYP3A4 Ki = 0.15 uM, on CYP3A5 Ki = 0.11 uM); Discussion repeats",
    "both. Time-dependent inactivation constants: Results 2.2, Discussion, and",
    "the Figure 3B annotation (kinact = 0.11 /min, KI = 2.45 uM for CYP3A4);",
    "STA produced no TDI on CYP3A5 (Figure 3C). Michaelis-Menten inactivation",
    "form kobs = kinact*I/(KI+I) and its estimation by the double-reciprocal",
    "plot: Materials and Methods 4.4 and Equation (1). Assay design: Methods",
    "4.3 (RI) and 4.4 (TDI); the underlying data points are in Table S1.",
    sep = " "
  )

  vignette <- "He_2021_tacrolimus_wuzhi_cyp3a"

  units <- list(
    time = "min",
    dosing = "(none; the inhibitor concentration is held constant through the preincubation and is supplied as the covariate CP_STA_UM)",
    concentration = "(the states and outputs enzyme_3a4 and enzyme_3a5 are CYP catalytic activity as a fraction of the no-inhibitor control, dimensionless; the driving covariate CP_STA_UM is the schisantherin A concentration in the incubation, in uM)"
  )

  compartmentData <- list(
    enzyme_3a4 = list(
      analyte = "Catalytically active CYP3A4, expressed as a fraction of the no-inhibitor control activity and read out as the rate of testosterone 6-beta-hydroxylation",
      units = "fraction of control (dimensionless)",
      specimen = "not applicable",
      verified = TRUE
    ),
    enzyme_3a5 = list(
      analyte = "Catalytically active CYP3A5, expressed as a fraction of the no-inhibitor control activity; STA does not inactivate CYP3A5 (Figure 3C) so this state stays at its baseline of 1 and only enters the reversible-inhibition term",
      units = "fraction of control (dimensionless)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    CP_STA_UM = list(
      description = "Concentration of schisantherin A (STA) in the incubation, supplied as a covariate. This is the quantity the source calls I, the inhibitor concentration, in the inactivation-rate expression and in the Dixon-plot reversible-inhibition term. In this in-vitro model the column carries an incubation-buffer concentration rather than a plasma concentration, but the quantity and units (uM) are identical; same reuse rationale as CP_RIF_UM in Almond_2016_rifampicin_invitro.R.",
      units = "umol/L (uM)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "TDI assay (Methods 4.4): STA at 0, 0.25, 0.5, 1, 2 uM. RI assay (Methods 4.3): STA at 0, 0.125, 0.25, 0.5 uM.",
        "Held CONSTANT over the preincubation by the assay design (no depletion correction), so a static covariate value is the faithful encoding.",
        "Set to 0 for the no-inhibitor control, at which kobs = 0 and enzyme_3a4 stays at 1 for all time and the reversible term is 1."
      ),
      source_name = "I (inhibitor concentration)"
    )
  )

  population <- list(
    species = "in vitro (human liver microsomes)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    system = "CYP3A5*3/*3 human liver microsomes (for CYP3A4) and CYP3A5*1/*3 human liver microsomes co-incubated with the selective CYP3A4 inactivator CYP3cide (for CYP3A5). RI assay at 0.2 mg/mL microsomal protein; TDI assay at 0.5 mg/mL microsomal protein, in 0.1 M potassium phosphate buffer pH 7.4.",
    temperature = "37 C",
    probe_reaction = "Tacrolimus disappearance (Dixon plots, reversible inhibition) and testosterone 6-beta-hydroxylation (time-dependent inactivation)",
    disease_state = "not applicable (in vitro)",
    notes = paste(
      "The reported +/- terms on the TDI constants are standard errors of the nonlinear regression, not between-subject variability, so no omega is encoded.",
      "STA is a time-dependent AND reversible inhibitor of CYP3A4 but only a reversible inhibitor of CYP3A5 (no TDI, Figure 3C).",
      "kinact is reported in per-minute units (Figure 3B), consistent with the per-minute time axis of the assay; the model's time unit is minutes.",
      "The CYP degradation rate constant kdeg, which governs the in-vivo magnitude of TDI, is a Simcyp default and is not printed; no enzyme-turnover term is carried, so the model describes the inactivation phase of the in-vitro assay only."
    )
  )

  ini({
    # =====================================================================
    # Reversible (competitive) inhibition constants, from the Dixon plots.
    # STA on CYP3A4: Ki = 0.15 uM (Figure 2A); on CYP3A5: Ki = 0.11 uM
    # (Figure 2B). Both are repeated in Results 2.1 and the Discussion.
    # No uncertainty is reported, so they are carried as point values.
    # =====================================================================
    ki_3a4 <- 0.15
    label("Schisantherin A equilibrium dissociation constant Ki for reversible CYP3A4 inhibition (uM)") # Figure 2A and Results 2.1: 'Ki values for the inhibition by STA on CYP3A4 ... of 0.15 uM'

    ki_3a5 <- 0.11
    label("Schisantherin A equilibrium dissociation constant Ki for reversible CYP3A5 inhibition (uM)") # Figure 2B and Results 2.1: 'Ki values for the inhibition by STA on ... CYP3A5 of ... 0.11 uM'

    # =====================================================================
    # Time-dependent (mechanism-based) inactivation of CYP3A4, from the
    # double-reciprocal (kobs vs 1/[STA]) plot of Figure 3B.
    # STA showed NO time-dependent inactivation of CYP3A5 (Figure 3C), so
    # there is no CYP3A5 inactivation term.
    # =====================================================================
    ki_inact_3a4 <- 2.45
    label("Schisantherin A concentration at half-maximal CYP3A4 inactivation rate KI (uM)") # Figure 3B annotation and Discussion: 'KI = 2.45 uM'

    kinact_3a4 <- 0.11
    label("Maximum rate constant of CYP3A4 inactivation by schisantherin A kinact (1/min)") # Figure 3B annotation and Discussion: 'kinact = 0.11 min-1'

    # =====================================================================
    # Residual error is NOT reported (the source fits Dixon and
    # double-reciprocal plots and reports only the resulting constants).
    # Per the standing policy on unreported residual error the term is
    # fixed at zero so the model returns the deterministic published curve.
    # Flagged in the vignette Errata.
    # =====================================================================
    addSd <- fixed(0)
    label("Additive residual SD of the relative CYP3A4 activity, ZERO because the source reports no residual-error model (fraction of control)")
  })

  model({
    # ===================================================================
    # Reversible-inhibition factors. A victim-drug intrinsic clearance by
    # each isoform is scaled by 1/(1 + I/Ki) under competitive inhibition;
    # these are exposed as derived outputs so a downstream user can read
    # the fractional CYP3A4/CYP3A5 activity loss at any STA concentration.
    # ===================================================================
    riFactor_3a4 <- 1 / (1 + CP_STA_UM / ki_3a4)
    riFactor_3a5 <- 1 / (1 + CP_STA_UM / ki_3a5)

    # ===================================================================
    # Time-dependent inactivation of CYP3A4. The observed first-order rate
    # constant is a hyperbolic (Michaelis-Menten) function of the STA
    # concentration; its plateau is kinact and it passes through kinact/2
    # at I = KI. enzyme_3a4 is the registered CYP3A4 activity state as a
    # fraction of its untreated baseline, enzyme_3a4(0) = 1.
    # ===================================================================
    kobs_3a4 <- kinact_3a4 * CP_STA_UM / (ki_inact_3a4 + CP_STA_UM)

    enzyme_3a4(0) <- 1
    d/dt(enzyme_3a4) <- -kobs_3a4 * enzyme_3a4

    # CYP3A5 has no time-dependent inactivation by STA (Figure 3C): its
    # activity state holds at baseline and only the reversible term acts.
    enzyme_3a5(0) <- 1
    d/dt(enzyme_3a5) <- 0

    # Figure 3A plots the log of percent activity remaining; exposed as a
    # derived output so the panel can be replicated directly.
    pctActivity_3a4 <- 100 * enzyme_3a4

    enzyme_3a4 ~ add(addSd)
  })
}
