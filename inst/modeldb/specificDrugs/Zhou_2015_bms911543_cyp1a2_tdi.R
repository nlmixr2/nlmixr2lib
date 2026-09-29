Zhou_2015_bms911543_cyp1a2_tdi <- function() {
  description <- paste(
    "In vitro (human liver microsomes). Mechanism-based (time-dependent)",
    "inactivation model of CYP1A2 by the JAK2 inhibitor BMS-911543, measured",
    "with phenacetin O-deethylation to acetaminophen as the CYP1A2 activity",
    "probe. The published model is the standard two-equation inactivation",
    "form: the probe velocity decays exponentially during preincubation,",
    "vI = v0 * exp(-k * t), with a hyperbolic concentration dependence of the",
    "observed inactivation rate constant, k = kinact * I / (KI + I). This file",
    "encodes it as the equivalent first-order ODE on the relative CYP1A2",
    "activity state, d/dt(enzyme_1a2) = -k * enzyme_1a2 with",
    "enzyme_1a2(0) = 1, driven by the preincubation inhibitor concentration",
    "supplied as the covariate CP_BMS911543_UM. This inactivation is the",
    "mechanism the paper identifies as the key driver of the time-dependent",
    "pharmacokinetic nonlinearity BMS-911543 showed in its first-in-human",
    "study. Siblings from the same paper: Zhou_2015_bms911543_hlm,",
    "Zhou_2015_bms911543_rcyp1a2, Zhou_2015_bms911543_rcyp3a4 and",
    "Zhou_2015_bms911543_rcyp2j2, the Michaelis-Menten characterisation of",
    "the metabolite-M1 formation this same enzyme catalyses.",
    "The KI carried here is the EXPERIMENTAL value of 2.9 uM from Figure 2b,",
    "not the 11.2 uM the authors substituted into their Simcyp V12 whole-body",
    "PBPK model to improve the fit to the clinical data. That whole-body model",
    "is NOT part of this file and is not reproducible from the published",
    "inputs; neither is the CYP1A2 autoinduction arm, whose Emax and EC50 were",
    "read from a supplementary table that was not available when this model was built. See the validation",
    "vignette for the full accounting.",
    sep = " "
  )

  reference <- paste(
    "Zhou L, Gan J, Yoshitsugu H, Gu X, Lutz JD, Masson E, Humphreys WG.",
    "Integration of Physiologically-Based Pharmacokinetic Modeling into Early",
    "Clinical Development: An Investigation of the Pharmacokinetic",
    "Nonlinearity.",
    "CPT Pharmacometrics Syst Pharmacol. 2015;4(5):286-294.",
    "doi:10.1002/psp4.35. PMCID: PMC4452934.",
    "The two model equations vI = v0 * exp(-k * t) (Eq. 1) and",
    "k = kinact * I / (KI + I) (Eq. 2), and the definitions of t, vI, v0 and I:",
    "Methods, 'CYP1A2 time-dependent inhibition'. The parameter estimates",
    "KI = 2.9 +/- 0.9 uM and kinact = 1.4 +/- 0.1 per h are given twice, in",
    "Results, 'CYP1A2 TDI by BMS-911543', and as the annotation on Figure 2b.",
    "The assay design (1.0 mg/mL human liver microsomes in 100 mM potassium",
    "phosphate buffer pH 7.4 at 37 C; 5 min preincubation; inactivation",
    "initiated with 1 mM NADPH and run for 3, 10, 20 and 30 min in the absence",
    "of phenacetin; 10-fold dilution into 450 uM phenacetin for a 13.5 min",
    "probe incubation; inhibitor at 0, 0.39, 0.78, 1.56, 3.125, 6.25, 12.5 and",
    "25.0 uM; triplicate; nonlinear regression in GraphPad Prism v.5):",
    "Methods, 'CYP1A2 time-dependent inhibition'. The 11.2 uM KI substituted",
    "into the Simcyp model, and the statement that the other TDI and induction",
    "inputs were used at their experimental values: Methods, 'PBPK modeling",
    "and simulation'. The best-fit clinical KI of 11 +/- 3.4 uM and the",
    "acknowledgement that Km and KI are consequently not identifiable:",
    "Discussion.",
    sep = " "
  )

  vignette <- "Zhou_2015_bms911543_invitro"

  units <- list(
    time = "h",
    dosing = "(none; the inhibitor concentration is held constant through the preincubation and is supplied as the covariate CP_BMS911543_UM)",
    concentration = "(the state and output enzyme_1a2 is CYP1A2 catalytic activity as a fraction of the no-inhibitor control, dimensionless; the driving covariate CP_BMS911543_UM is the BMS-911543 concentration in the preincubation in uM)"
  )

  compartmentData <- list(
    enzyme_1a2 = list(
      analyte = "Catalytically active CYP1A2, expressed as a fraction of the no-inhibitor control activity and read out as the rate of phenacetin O-deethylation to acetaminophen",
      units = "fraction of control (dimensionless)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    CP_BMS911543_UM = list(
      description = "Concentration of BMS-911543 in the preincubation, supplied as a covariate. This is the quantity the source calls I, 'the specific initial inhibitor concentration', in Eq. 2. Reused canonical: in this in-vitro model the column carries an incubation-buffer concentration rather than a plasma concentration, but the quantity and units (uM) are identical. Same reuse rationale as CP_RIF_UM in Almond_2016_rifampicin_invitro.R, which carries a hepatocyte culture-medium concentration.",
      units = "umol/L (uM)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Zhou 2015 Methods, 'CYP1A2 time-dependent inhibition': BMS-911543 was preincubated at 0, 0.39, 0.78, 1.56, 3.125, 6.25, 12.5 and 25.0 uM. The Figure 2a legend rounds these to 0, 0.39, 0.78, 1.6, 3.1, 6.3, 13 and 25 uM.",
        "The concentration is held CONSTANT over the preincubation by the assay design (no depletion correction is applied by the authors), so a static covariate value is the faithful encoding. Applying this model to an in-vivo profile instead requires a time-varying unbound concentration at the site of enzyme interaction, and the source gives no basis for that translation.",
        "Set to 0 for the no-inhibitor control condition, at which k = 0 and the activity state stays at 1 for all time, by construction.",
        "The top concentration of 25 uM is 8.6-fold above the fitted KI of 2.9 uM, so the plateau of the Figure 2b hyperbola is well determined: the fitted k at 25 uM is 1.254 per h, which is 90 percent of kinact."
      ),
      source_name = "I (initial inhibitor concentration)"
    )
  )

  population <- list(
    species = "in vitro (human liver microsomes)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    system = "Human liver microsomes at 1.0 mg/mL in 100 mM potassium phosphate buffer, pH 7.4. Note this is a fourfold higher protein concentration than the 0.25 mg/mL used for the metabolite-M1 formation kinetics in the sibling files, so the two arms' fractions unbound in incubation are not interchangeable.",
    temperature = "37 C",
    kinetic_incubation = "5 min preincubation with inhibitor, then inactivation initiated with 1 mM nicotinamide adenine dinucleotide phosphate and allowed to proceed for 3, 10, 20 or 30 min in the absence of phenacetin; aliquots then diluted 10-fold into 1 mM NADPH plus 450 uM phenacetin and incubated 13.5 min before quenching and LC/MS/MS analysis",
    concentration_range = "0, 0.39, 0.78, 1.56, 3.125, 6.25, 12.5 and 25.0 uM BMS-911543",
    replication = "Values are the mean of three measurements; inhibition constants were determined by nonlinear regression in GraphPad Prism v.5 and are presented as mean and SE",
    probe_reaction = "Phenacetin O-deethylation to acetaminophen, the CYP1A2 activity probe reaction",
    disease_state = "not applicable (in vitro)",
    notes = paste(
      "This experiment is the paper's central piece of reverse translation. CYP1A2 time-dependent inactivation had NOT been characterised before the first-in-human study; only a modest half-maximal inhibitory concentration against tacrine had been seen. When the initial PBPK model built on the pre-clinical package failed to reproduce the accumulation observed after two weeks of twice-daily dosing, a sensitivity analysis pointed at the TDI parameters, and this assay was run to measure them.",
      "The 10-fold dilution step before the probe incubation is what makes the assay a MECHANISM-BASED inactivation measurement rather than a reversible-inhibition measurement: reversible inhibition is largely relieved by the dilution, so residual loss of activity is attributed to irreversible inactivation.",
      "The reported uncertainties (KI +/- 0.9 uM, kinact +/- 0.1 per h) are standard errors of the nonlinear regression, NOT between-subject or between-donor variability, and are therefore not encoded as an omega. No inter-individual variability of any kind is reported for this assay.",
      "The rate constant of CYP1A2 degradation, kdeg, is the third parameter that governs the in-vivo magnitude of TDI. The authors state that they used the Simcyp default value for simplicity and never print it, so no enzyme-turnover term is carried in this file; the model describes the inactivation phase of the in-vitro assay only, over which resynthesis is negligible.",
      "For the same reason the CYP1A2 autoinduction arm of the paper is absent here: its maximal fold induction and EC50 were taken from experimental values reported in Supplementary Table S3, which was not available when this model was built. The simulated interplay of the two effects (62 percent of hepatic CYP1A2 activity remaining at day 15 with induction, 54 percent without, so induction attenuates the inactivation by 17 percent) is a Simcyp platform output and is recorded in the vignette rather than reproduced."
    )
  )

  ini({
    # =====================================================================
    # Time-dependent inactivation constants for CYP1A2, from the nonlinear
    # regression of Figure 2b. The values are printed identically in
    # Results, 'CYP1A2 TDI by BMS-911543', and on the figure panel.
    #
    # The +/- terms printed alongside are standard errors of the
    # regression, not variability, so they are recorded in `population`
    # rather than encoded as an omega.
    #
    # NOTE on which KI this is: the authors deliberately used a DIFFERENT
    # KI of 11.2 uM inside their Simcyp whole-body model, because the
    # experimental 2.9 uM did not reproduce the clinical accumulation
    # (Methods, 'PBPK modeling and simulation'; Discussion gives the
    # best-fit clinical value as 11 +/- 3.4 uM). This file is the IN VITRO
    # model, so it carries the in vitro measurement. The discrepancy is
    # the paper's own headline finding about the difficulty of translating
    # TDI parameters, and is discussed in the vignette.
    # =====================================================================
    ki_1a2 <- 2.9
    label("Inhibitor concentration producing half the maximum rate of CYP1A2 inactivation (uM)") # Figure 2b annotation and Results, 'CYP1A2 TDI by BMS-911543': KI = 2.9 +/- 0.9 uM

    kinact_1a2 <- 1.4
    label("Maximum rate constant of CYP1A2 inactivation (1/h)") # Figure 2b annotation and Results, 'CYP1A2 TDI by BMS-911543': kinact = 1.4 +/- 0.1 per h

    # =====================================================================
    # Residual error is NOT reported. The source fitted Figure 2a by linear
    # regression of log percent activity against time and Figure 2b by
    # nonlinear regression, and reports only the resulting constants; there
    # is no residual-error model and no assay CV. Per the standing policy
    # on unreported residual error the term is fixed at zero so the model
    # returns the deterministic published curve. Flagged in the vignette
    # Errata.
    # =====================================================================
    addSd <- fixed(0)
    label("Additive residual SD of the relative CYP1A2 activity, ZERO because the source reports no residual-error model (fraction of control)")
  })

  model({
    # ===================================================================
    # Eq. 2: the observed first-order rate constant of inactivation is a
    # hyperbolic function of the initial inhibitor concentration. This is
    # the curve drawn through the points of Figure 2b. Its plateau is
    # kinact and it passes through kinact/2 at I = KI, both of which the
    # vignette asserts.
    # ===================================================================
    kobs_1a2 <- kinact_1a2 * CP_BMS911543_UM / (ki_1a2 + CP_BMS911543_UM)

    # Figure 2b plots this constant in per-MINUTE units (its y-axis is
    # lambda in 1/min) while the reported kinact is in per-hour units.
    # The bridge is exposed as a derived output so the figure can be
    # replicated directly without rescaling outside the model.
    lambda_1a2_permin <- kobs_1a2 / 60

    # ===================================================================
    # Eq. 1 as an ODE. vI = v0 * exp(-k * t) is the solution of
    # d(activity)/dt = -k * activity with activity(0) = 1, where activity
    # is vI/v0, the probe velocity as a fraction of the no-inhibitor
    # control. `enzyme_1a2` is the registered canonical for a CYP1A2
    # activity state expressed as a fraction of its untreated baseline,
    # with the same enzyme_1a2(0) = 1 initial condition.
    #
    # Only the inactivation limb is modelled. There is no resynthesis term
    # because the source never reports kdeg, and none is needed for the
    # microsomal assay, in which no enzyme is made.
    # ===================================================================
    enzyme_1a2(0) <- 1
    d/dt(enzyme_1a2) <- -kobs_1a2 * enzyme_1a2

    # The y-axis of Figure 2a is percent activity remaining; exposed as a
    # derived output so that panel can be replicated directly.
    pctActivity_1a2 <- 100 * enzyme_1a2

    enzyme_1a2 ~ add(addSd)
  })
}
