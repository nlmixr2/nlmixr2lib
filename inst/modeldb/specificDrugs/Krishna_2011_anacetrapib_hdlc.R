Krishna_2011_anacetrapib_hdlc <- function() {
  description <- paste0(
    "Proportional Emax exposure-response model for high-density lipoprotein ",
    "cholesterol (HDL-C) raising by the cholesteryl ester transfer protein ",
    "(CETP) inhibitor anacetrapib, in healthy volunteers and patients with ",
    "dyslipidemia (Krishna 2011, Table II integrated dataset and Eq. 8). ",
    "The model is algebraic and has no ODE states: HDL-C at steady state is ",
    "the individual baseline multiplied by (1 + Emax * C / (EC50 + C)), ",
    "where C is the population-predicted anacetrapib trough concentration ",
    "supplied as the covariate column CSS_ANACETRAPIB. Inter-individual ",
    "variability is additive on the natural scale for both the baseline and ",
    "Emax, as fitted in S-PLUS NLME. Neither study population (healthy vs ",
    "dyslipidemic) nor concomitant statin altered the baseline or Emax, so ",
    "the model carries no covariates other than the exposure driver. Drive ",
    "it with the companion popPK model Krishna_2011_anacetrapib.R; the ",
    "LDL-C arm of the same paper is Krishna_2011_anacetrapib_ldlc.R."
  )
  reference <- paste0(
    "Krishna R, Bergman AJ, Green M, Dockendorf MF, Wagner JA, Dykstra K. ",
    "Model-based development of anacetrapib, a novel cholesteryl ester ",
    "transfer protein inhibitor. AAPS J. 2011;13(2):179-190. ",
    "doi:10.1208/s12248-011-9254-0. Structural model: Eq. 8 and Table II, ",
    "'Integrated dataset' column. PK layer: modellib('Krishna_2011_anacetrapib')."
  )
  vignette <- "Krishna_2011_anacetrapib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL", hdlc = "mg/dL")

  covariateData <- list(
    CSS_ANACETRAPIB = list(
      description        = "Population-predicted anacetrapib steady-state trough (24 h post-dose) plasma concentration",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "The sole driver of the model, entering Eq. 8 as ",
        "Emax * CSS_ANACETRAPIB / (EC50 + CSS_ANACETRAPIB). Krishna 2011 ",
        "screened trough and 24 h average CETP inhibition, trough and 24 h ",
        "average CETP activity, anacetrapib trough concentration and ",
        "anacetrapib daily average concentration, and found the trough ",
        "concentration the most predictive of lipid change; the CETP-activity ",
        "cascade was abandoned because CETP activity rebounds with repeated ",
        "dosing while the lipid effects do not. The value is a ",
        "POPULATION-predicted trough, not an individual post hoc one: Table ",
        "III presents both columns and the authors selected the ",
        "population-predicted fit 'since this allows for prediction of ",
        "future lipid effects, given the treatment population and fed state, ",
        "in subjects who may or may not have PK data available'. To ",
        "reproduce that, solve Krishna_2011_anacetrapib.R under ",
        "rxode2::zeroRe() at the regimen, formulation and prandial state of ",
        "interest, take the pre-dose concentration at steady state, and ",
        "supply it here; do NOT propagate PK inter-individual variability ",
        "into this column, because the fitted PD variances below were ",
        "estimated against a population-predicted driver and would ",
        "double-count it. Set to 0 for placebo. Krishna 2011 observed ",
        "troughs above 400 ng/mL only in six normal volunteers given 400 mg ",
        "with a high-fat meal, where the model shows some lack of fit."
      ),
      source_name        = "C_trough (C24h)"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 545L,
    n_observations   = 3089L,
    n_studies        = 5L,
    disease_state    = "Combined phase I and phase IIb dataset of 474 patients with dyslipidemia and 72 normal healthy volunteers (546 subjects); the HDL-C model converged on 545 subjects and 3,089 observations (Table II, integrated dataset).",
    hdlc_baseline    = "Typical baseline HDL-C 50.8 mg/dL (SE 0.5) with 11.7 mg/dL between-subject SD, i.e. a 23% coefficient of variation. Phase I alone gave 53.8 mg/dL and phase IIb alone 50.6 mg/dL, so the pooled estimate is dominated by the much larger phase IIb patient cohort.",
    dose_range       = "Anacetrapib doses across the pooled phase I and phase IIb studies; the model was exercised over simulated arms of placebo and 10, 50, 100, 150, 250 and 300 mg once daily. 100 mg once daily as the hot-melt-extruded tablet was selected for phase III.",
    co_medication    = "Both statin-treated and statin-free participants contributed. No effect of statin administration on baseline HDL-C or on Emax was identified, so atorvastatin does not appear in this model (contrast the LDL-C arm, where it does).",
    notes            = paste0(
      "Per-study designs are in Table IB of the Electronic Supplementary ",
      "Material, which is not on disk (see the vignette Errata). Krishna ",
      "2011 publishes no age, weight, sex or race distribution for the ",
      "exposure-response analysis set in the main text, and no demographic ",
      "covariate is in the model."
    )
  )

  ini({
    # =========================================================================
    # Krishna 2011 Table II, "Integrated dataset" column (phase I + phase IIb),
    # and Eq. 8:
    #
    #   HDL_ij = (BL + eta_BL_i) * (1 + (Emax + eta_Emax_i) * C_i /
    #                                    (EC50 + C_i)) + eps_ij
    #
    # Both random effects are ADDITIVE on the natural scale, which is what
    # Eq. 8 prints and what S-PLUS NLME (the software named in the Methods)
    # parameterises; they are NOT log-normal. The omega and sigma rows are
    # therefore read as STANDARD DEVIATIONS in mg/dL, the nlme reporting
    # convention -- and the arithmetic agrees: 11.7 mg/dL on a 50.8 mg/dL
    # baseline is a 23% CV, whereas reading it as a variance would give
    # sqrt(11.7) = 3.4 mg/dL, a 6.7% CV, far too small for population HDL-C.
    # Note this is the OPPOSITE scale convention from the NONMEM-fitted popPK
    # file, whose Table I omegas are provably variances; the two tables come
    # from two different pieces of software. Both readings are recorded in
    # the vignette Assumptions and deviations section.
    #
    # Parameters are on the natural scale (not log) so that each eta is
    # mu-referenced as theta + eta, matching the printed equation.
    # =========================================================================
    rbase <- 50.8;    label("Typical baseline HDL-C (mg/dL)")                                        # Table II, integrated dataset: Baseline HDL-C = 50.8 mg/dL (SE 0.5); phase I alone 53.8, phase IIb alone 50.6
    emax  <- 1.76;    label("Maximum proportional increase in HDL-C above baseline (unitless)")       # Table II, integrated dataset: Emax = 1.76 (SE 0.07); Results: "an increase of 176% (95% CI, 162% to 190%) over baseline"
    lec50 <- log(135); label("Anacetrapib trough concentration giving half-maximal HDL-C effect (ng/mL)")  # Table II, integrated dataset: EC50 = 135 ng/mL (8% approximate CV); phase I alone 140 (13%), phase IIb alone 108 (10%)

    # ---- IIV, additive on the natural scale (Table II) ----
    etarbase ~ 136.89   # Table II, integrated dataset: w_BL = 11.7 mg/dL (SE 0.4), read as an SD -> variance 11.7^2 = 136.89
    etaemax  ~ 0.16     # Table II, integrated dataset: w_Emax = 0.40 (SE 0.04), read as an SD -> variance 0.40^2 = 0.16; not estimable from the phase I data alone

    # ---- Residual error (Table II) ----
    addSd <- 7.6; label("Additive residual error on HDL-C (mg/dL)")   # Table II, integrated dataset: sigma = 7.6 mg/dL (SE 0.1); Eq. 8 adds eps_ij outside the proportional bracket, so the residual is purely additive
  })

  model({
    # Eq. 8. Individual baseline and Emax are the typical value plus an
    # additive random effect, exactly as printed.
    rbase_i <- rbase + etarbase
    emax_i  <- emax + etaemax
    ec50    <- exp(lec50)

    hdl <- rbase_i * (1 + emax_i * CSS_ANACETRAPIB / (ec50 + CSS_ANACETRAPIB))

    hdl ~ add(addSd)
  })
}
