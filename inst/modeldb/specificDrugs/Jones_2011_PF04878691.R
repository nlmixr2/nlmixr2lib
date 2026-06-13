Jones_2011_PF04878691 <- function() {
  description <- paste(
    "Two-compartment population pharmacokinetic model with first-order",
    "oral absorption and time-varying clearance for the toll-like-receptor-7",
    "(TLR7) agonist PF-04878691 in healthy male and female adult volunteers",
    "(Jones 2011 BJCP, Phase 1 multiple-dose escalation study, twice-weekly",
    "oral doses of 3, 6, or 9 mg over 2 weeks). Observed plasma exposure",
    "increased over the dosing period inconsistently with the 12-16 h",
    "terminal half-life; a standard linear time-invariant two-compartment",
    "model over-estimated Cmax on day 1 and under-estimated exposure on",
    "day 11. The clearance was therefore parameterised with an exponentially",
    "decaying time-dependent component superimposed on a steady-state arm:",
    "CL(t) = CL_SS + CL_TIME * exp(-kdeg * TAFD), where TAFD is the time",
    "after first dose. Reparameterised from the paper's CLF (final = CL_SS)",
    "and CL0 (initial = CL_SS + CL_TIME). The hypothesised mechanism for the",
    "time-varying clearance is IFN-mediated CYP1A2 inhibition by the",
    "TLR7-induced interferon response (Discussion). All disposition",
    "parameters are reported per kilogram body weight (paper: doses were",
    "body-weight-normalised so estimated PK parameters carry per-kg units);",
    "WT is therefore a required covariate. No other covariates retained."
  )
  reference <- paste(
    "Jones HM, Chan PLS, van der Graaf PH, Webster R. Use of modelling and",
    "simulation techniques to support decision making on the progression of",
    "PF-04878691, a TLR7 agonist being developed for hepatitis C.",
    "Br J Clin Pharmacol. 2012;73(1):77-92.",
    "doi:10.1111/j.1365-2125.2011.04047.x"
  )
  vignette <- "Jones_2011_PF04878691_HCV"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required scaling covariate. Jones 2011 estimated all disposition parameters per kilogram body weight (Table 1 unit row 'l h^-1 kg^-1' for clearances and 'l kg^-1' for volumes); per the Methods, 'all doses were normalised for body weight to ensure that any IIV estimated was purely due to variability in the parameter rather than dose'. Population median 79 kg (range 57-97 kg).",
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 24L,
    n_studies      = 1L,
    age_range      = "21-55 years",
    age_median     = "34 years",
    weight_range   = "57-97 kg",
    weight_median  = "79 kg",
    sex_female_pct = 8,
    race_ethnicity = "Not tabulated in Jones 2011.",
    disease_state  = "Healthy adult volunteers (multiple-dose escalation Phase 1 study; ClinicalTrials.gov NCT00810758).",
    dose_range     = "PF-04878691 administered orally as an extemporaneously-prepared solution at 3, 6, or 9 mg twice weekly (days 1, 4, 8, 11) for 2 weeks; n = 6 active per dose cohort plus n = 2 placebo per cohort. Last four subjects withdrew during active treatment after two doses following two SAEs in the 9 mg cohort (study prematurely terminated).",
    regions        = "Not specified.",
    notes          = "Median age and weight from Jones 2011 Methods ('Clinical TLR7 study data'). Two female subjects (8%) and 22 males (92%). Serial PK sampling pre-dose and up to 312 h after the last dose; LLOQ 0.1 ng/mL; HPLC-MS/MS, inter-/intra-assay CV < 4.9%."
  )

  ini({
    # Structural PK parameters (Jones 2011 Table 1). All clearances and
    # volumes were estimated per kilogram body weight; the per-subject
    # values are recovered inside model() as exp(l<param>) * WT.
    #
    # Time-varying clearance reparameterisation. Jones 2011 estimated CLF
    # (final / steady-state CL) and CL0 (initial CL at TAFD = 0); the
    # paper's clearance equation is
    #
    #     CL(TAFD) = CLF + (CL0 - CLF) * exp(-DEG * TAFD)
    #
    # rewritten here in the canonical CL_SS + CL_TIME * exp(-kdeg * t)
    # decomposition used throughout nlmixr2lib (Gibiansky 2014, Lu 2019,
    # Wu 2024, Hussein 1997, ...). The point-estimate mapping is:
    #
    #     cl_ss   = CLF = 1.7  L/h/kg
    #     cl_time = CL0 - CLF = 3.5 - 1.7 = 1.8 L/h/kg (initial offset)
    #     kdeg    = DEG = 0.24 1/h
    #
    # At t = 0 the total clearance is CL_SS + CL_TIME = 3.5 L/h/kg = CL0;
    # at t -> infinity it converges to CL_SS = 1.7 L/h/kg = CLF.
    lcl       <- log(1.7);    label("Steady-state apparent clearance per kg body weight (CL_SS = paper CLF, L/h/kg)")              # Table 1 (CLF = 1.7 L/h/kg, %CV 6.8)
    lcl_time  <- log(1.8);    label("Initial offset of the time-varying clearance component per kg (CL_TIME0 = CL0 - CLF, L/h/kg)") # Derived from Table 1 (CL0 = 3.5, CLF = 1.7)
    lkdeg     <- log(0.24);   label("Exponential decay rate of the time-varying clearance component (paper DEG, 1/h)")              # Table 1 (DEG = 0.24 1/h, %CV 35)
    lvc       <- log(3.3);    label("Apparent central volume of distribution per kg body weight (Vc, L/kg)")                        # Table 1 (Vc = 3.3 L/kg, %CV 47)
    lq        <- log(0.74);   label("Apparent intercompartmental clearance per kg body weight (Q, L/h/kg)")                         # Table 1 (Q = 0.74 L/h/kg, %CV 51)
    lvp       <- log(21);     label("Apparent peripheral volume of distribution per kg body weight (Vp, L/kg)")                     # Table 1 (Vp = 21 L/kg, %CV 17)
    lka       <- log(0.078);  label("First-order absorption rate constant (ka, 1/h)")                                               # Table 1 (ka = 0.078 1/h, %CV 22)

    # Inter-individual variability. Table 1 reports OM1 (IIV on CLF, here
    # IIV on cl_ss) and OM3 (IIV on ka). OM2 is absent in the source
    # table -- no IIV was retained on Vc, Q, Vp, CL0 (the time-varying
    # initial offset), or DEG. The 'l h^-1 kg^-1' units of CLF mean the
    # IIV is on the log-scale variance of the per-kg parameter, which is
    # what etalcl applies to inside model().
    etalcl ~ 0.067                                                                                                                  # Table 1 (IIV CLF = 0.067, %CV 28)
    etalka ~ 0.19                                                                                                                   # Table 1 (IIV ka  = 0.19,  %CV 36)

    # Residual error. Methods state 'intra-individual (residual)
    # variability was generally described using a proportional error
    # model' for PK, OAS and lymphocyte models. Table 1 reports the PK
    # residual SD as 0.046 (%CV 18).
    propSd <- 0.046; label("Proportional residual SD (fraction)")                                                                   # Table 1 (residual error = 0.046, %CV 18)
  })

  model({
    # ----------------------------------------------------------------------
    # Individual disposition parameters (per kg) and per-subject scaling
    # to absolute units via body weight. cl_ss carries the per-subject log-
    # normal IIV (paper IIV on CLF); the time-varying component CL_TIME and
    # the decay rate kdeg are population-typical with no IIV.
    # ----------------------------------------------------------------------
    cl_ss   <- exp(lcl + etalcl) * WT
    cl_time <- exp(lcl_time) * WT
    kdeg    <- exp(lkdeg)
    vc      <- exp(lvc) * WT
    q       <- exp(lq)  * WT
    vp      <- exp(lvp) * WT
    ka      <- exp(lka + etalka)

    # Time-varying clearance: CL(t) = CL_SS + CL_TIME * exp(-kdeg * t).
    # `time` is the rxode2 model time; with dosing initiated at t = 0 this
    # is equivalent to the paper's TAFD (time after first dose).
    cl <- cl_ss + cl_time * exp(-kdeg * time)

    # Microconstants for the two-compartment disposition
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # Plasma concentration in the central compartment. Dose is in mg and
    # vc is in L so central/vc has units mg/L = ug/mL; the multiplier 1000
    # converts to ng/mL because Jones 2011 reports concentrations in ng/mL
    # (LLOQ 0.1 ng/mL) and the OAS / lymphocyte / viral-load downstream
    # PD parameter values (slope, gamma) in Tables 2-4 were estimated on
    # data with Cc in ng/mL.
    Cc <- 1000 * central / vc

    Cc ~ prop(propSd)
  })
}
