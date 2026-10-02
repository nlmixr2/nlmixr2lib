He_2021_schisandrinA_cyp3a <- function() {
  description <- paste(
    "In vitro (CYP3A5-genotyped human liver microsomes). Inhibition of CYP3A4",
    "and CYP3A5 by schisandrin A (SIA), a lignan of the Wuzhi capsule",
    "(Schisandra sphenanthera extract), which is co-prescribed with tacrolimus",
    "in China. Two independently fitted mechanisms.",
    "(1) Time-dependent (mechanism-based) inactivation of BOTH CYP3A4 (in",
    "CYP3A5*3/*3 microsomes) and CYP3A5 (in CYP3A5*1/*3 microsomes plus",
    "CYP3cide), measured with testosterone 6-beta-hydroxylation as the probe;",
    "each observed inactivation rate constant is kobs = kinact * I / (KI + I)",
    "(double-reciprocal plots, Figure 4B and 4D), encoded as",
    "d/dt(enzyme_3a4) = -kobs_3a4 * enzyme_3a4 and",
    "d/dt(enzyme_3a5) = -kobs_3a5 * enzyme_3a5 with both states starting at 1.",
    "(2) Reversible inhibition: SIA reversibly inhibits CYP3A5 with Ki = 8.74",
    "uM (Dixon plot, Figure 2D) but shows little reversible inhibition of",
    "CYP3A4 (Figure 2C), so only the CYP3A5 reversible constant is carried.",
    "The inhibitor concentration I is supplied as the covariate CP_SIA_UM.",
    "Sibling from the same paper: He_2021_schisantherinA_cyp3a (the other",
    "Wuzhi lignan, the more potent inhibitor). The whole-body Simcyp DDI model",
    "that consumes these constants is a vendor platform model and is NOT part",
    "of this file; see the vignette.",
    sep = " "
  )

  reference <- paste(
    "He Q, Bu F, Zhang H, Wang Q, Tang Z, Yuan J, Lin HS, Xiang X.",
    "Investigation of the Impact of CYP3A5 Polymorphism on Drug-Drug",
    "Interaction between Tacrolimus and Schisantherin A/Schisandrin A Based on",
    "Physiologically-Based Pharmacokinetic Modeling.",
    "Pharmaceuticals (Basel). 2021;14(3):198. doi:10.3390/ph14030198.",
    "PMCID: PMC7997453.",
    "Reversible inhibition constant on CYP3A5 (Ki = 8.74 uM): Results 2.1 and",
    "Figure 2D; little reversible inhibition of CYP3A4 (Figure 2C).",
    "Time-dependent inactivation constants: Results 2.2, Discussion, and the",
    "Figure 4B / 4D annotations (CYP3A4 kinact = 0.019 /min, KI = 2.54 uM;",
    "CYP3A5 kinact = 0.014 /min, KI = 2.07 uM). Michaelis-Menten inactivation",
    "form kobs = kinact*I/(KI+I) and its estimation by the double-reciprocal",
    "plot: Materials and Methods 4.4 and Equation (1). Assay design: Methods",
    "4.3 (RI) and 4.4 (TDI); the underlying data points are in Table S1.",
    sep = " "
  )

  vignette <- "He_2021_tacrolimus_wuzhi_cyp3a"

  units <- list(
    time = "min",
    dosing = "(none; the inhibitor concentration is held constant through the preincubation and is supplied as the covariate CP_SIA_UM)",
    concentration = "(the states and outputs enzyme_3a4 and enzyme_3a5 are CYP catalytic activity as a fraction of the no-inhibitor control, dimensionless; the driving covariate CP_SIA_UM is the schisandrin A concentration in the incubation, in uM)"
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
    CP_SIA_UM = list(
      description = "Concentration of schisandrin A (SIA) in the incubation, supplied as a covariate. This is the quantity the source calls I, the inhibitor concentration, in the inactivation-rate expression and in the Dixon-plot reversible-inhibition term. In this in-vitro model the column carries an incubation-buffer concentration rather than a plasma concentration, but the quantity and units (uM) are identical; same reuse rationale as CP_RIF_UM in Almond_2016_rifampicin_invitro.R.",
      units = "umol/L (uM)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "TDI assay (Methods 4.4): SIA at 0, 2, 4, 8, 16 uM. RI assay (Methods 4.3): SIA at 0, 2.4, 7.2, 12 uM.",
        "Held CONSTANT over the preincubation by the assay design (no depletion correction), so a static covariate value is the faithful encoding.",
        "Set to 0 for the no-inhibitor control, at which both kobs values are 0 and the activity states stay at 1 for all time and the reversible term is 1."
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
      "SIA inhibits BOTH CYP3A4 and CYP3A5 in a time-dependent manner but reversibly inhibits only CYP3A5 (little reversible inhibition of CYP3A4, Figure 2C). Both inhibition mechanisms are weaker than schisantherin A's.",
      "kinact is reported in per-minute units (Figures 4B, 4D), consistent with the per-minute time axis of the assay; the model's time unit is minutes.",
      "The CYP degradation rate constant kdeg, which governs the in-vivo magnitude of TDI, is a Simcyp default and is not printed; no enzyme-turnover term is carried, so the model describes the inactivation phase of the in-vitro assay only."
    )
  )

  ini({
    # =====================================================================
    # Reversible (competitive) inhibition of CYP3A5, Ki = 8.74 uM, from the
    # Dixon plot of Figure 2D. SIA showed little reversible inhibition of
    # CYP3A4 (Figure 2C), so no CYP3A4 reversible constant is carried.
    # No uncertainty is reported, so it is a point value.
    # =====================================================================
    ki_3a5 <- 8.74
    label("Schisandrin A equilibrium dissociation constant Ki for reversible CYP3A5 inhibition (uM)") # Figure 2D and Results 2.1: 'The value of Ki for the inhibition on CYP3A5 by SIA was 8.74 uM'

    # =====================================================================
    # Time-dependent (mechanism-based) inactivation, from the
    # double-reciprocal plots. CYP3A4: kinact = 0.019 /min, KI = 2.54 uM
    # (Figure 4B). CYP3A5: kinact = 0.014 /min, KI = 2.07 uM (Figure 4D).
    # =====================================================================
    ki_inact_3a4 <- 2.54
    label("Schisandrin A concentration at half-maximal CYP3A4 inactivation rate KI (uM)") # Figure 4B annotation and Discussion: 'KI = 2.54 uM'

    kinact_3a4 <- 0.019
    label("Maximum rate constant of CYP3A4 inactivation by schisandrin A kinact (1/min)") # Figure 4B annotation and Discussion: 'kinact = 0.019 min-1'

    ki_inact_3a5 <- 2.07
    label("Schisandrin A concentration at half-maximal CYP3A5 inactivation rate KI (uM)") # Figure 4D annotation and Discussion: 'KI = 2.07 uM'

    kinact_3a5 <- 0.014
    label("Maximum rate constant of CYP3A5 inactivation by schisandrin A kinact (1/min)") # Figure 4D annotation and Discussion: 'kinact = 0.014 min-1'

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
    # Reversible-inhibition factor on CYP3A5 (competitive, 1/(1 + I/Ki)).
    # Exposed as a derived output. No CYP3A4 reversible term (Figure 2C).
    # ===================================================================
    riFactor_3a5 <- 1 / (1 + CP_SIA_UM / ki_3a5)

    # ===================================================================
    # Time-dependent inactivation of CYP3A4 and CYP3A5. Each observed
    # first-order rate constant is a hyperbolic (Michaelis-Menten) function
    # of the SIA concentration; each plateau is its kinact and each passes
    # through kinact/2 at I = KI. Both enzyme states are registered CYP
    # activity states as a fraction of untreated baseline, starting at 1.
    # ===================================================================
    kobs_3a4 <- kinact_3a4 * CP_SIA_UM / (ki_inact_3a4 + CP_SIA_UM)
    kobs_3a5 <- kinact_3a5 * CP_SIA_UM / (ki_inact_3a5 + CP_SIA_UM)

    enzyme_3a4(0) <- 1
    d/dt(enzyme_3a4) <- -kobs_3a4 * enzyme_3a4

    enzyme_3a5(0) <- 1
    d/dt(enzyme_3a5) <- -kobs_3a5 * enzyme_3a5

    # Figures 4A/4C plot the log of percent activity remaining; exposed as
    # derived outputs so the panels can be replicated directly.
    pctActivity_3a4 <- 100 * enzyme_3a4
    pctActivity_3a5 <- 100 * enzyme_3a5

    enzyme_3a4 ~ add(addSd)
  })
}
