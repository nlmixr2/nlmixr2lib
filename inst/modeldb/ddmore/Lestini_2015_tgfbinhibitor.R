Lestini_2015_tgfbinhibitor <- function() {
  description <- "One-compartment first-order absorption PK with indirect-response biomarker turnover (E represents fractional inhibition of TGF-beta signalling) for a small-molecule TGF-beta inhibitor in oncology, simplified by Lestini 2015 from Bueno et al. for use as a population PK/PD test bench in adaptive-design simulations."
  reference <- paste(
    "Lestini G, Dumont C, Mentre F (2015).",
    "Influence of the Size of Cohorts in Adaptive Design for",
    "Nonlinear Mixed Effects Models: An Evaluation by Simulation",
    "for a Pharmacokinetic and Pharmacodynamic Model for a",
    "Biomarker in Oncology. Pharm Res 32(10):3159-69.",
    "doi:10.1007/s11095-015-1693-3.",
    "DDMORE Foundation Model Repository: DDMODEL00000192 (scenario 4).",
    sep = " "
  )
  vignette <- "Lestini_2015_tgfbinhibitor"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L", response = "fraction (0-1) of maximal TGF-beta signalling inhibition")

  ddmore_id    <- "DDMODEL00000192"
  replicate_of <- NULL

  covariateData <- list()

  population <- list(
    n_subjects     = 50L,
    n_studies      = 1L,
    disease_state  = "Simulated oncology cohort receiving a small-molecule TGF-beta inhibitor; the source publication is a simulation study and the bundle ships only a single-dose simulated dataset (Simulated_PKPD.txt).",
    dose_range     = "80 mg single oral dose (DDMODEL00000192 scenario 4 simulated design).",
    notes          = paste0(
      "DDMODEL00000192 is paired with the Lestini 2015 simulation study; ",
      "the bundle does not ship a real-data fit. Population-parameter ",
      "estimates come from the Monolix run reported in ",
      "Output_simulated_PKPD.txt, which re-fits the model on the bundled ",
      "Simulated_PKPD.txt event table (50 subjects, single 80 mg oral dose, ",
      "PK + PD samples at 0.1, 0.5, 1.5, 4, 6, and 12 time-units after ",
      "dose). The bundle does not declare the time unit; this implementation ",
      "labels it 'h' on the basis that the resulting elimination half-life ",
      "(~6.6 h) and biomarker turnover half-life (~2.6 h) are biologically ",
      "consistent with a small-molecule TGF-beta inhibitor. The full ",
      "publication is paywalled and was not available at extraction time, ",
      "so per-cohort demographics, sampling-grid units, and the simulation ",
      "'truth' parameter values are not directly verifiable; see the ",
      "Lestini_2015_tgfbinhibitor vignette Errata for the full deviation list."
    ),
    n_subjects_simulated_dataset = 50L,
    sampling_grid_simulated_dataset = "0.1, 0.5, 1.5, 4, 6, 12 (time units, labelled hours in this implementation)"
  )

  ini({
    # ----------------------------------------------------------------------
    # Structural PK and PD parameter values come from the Monolix MLE refit
    # reported in Output_simulated_PKPD.txt (the only fit listing shipped in
    # the DDMORE bundle). They are estimates from one realisation of the
    # bundle's Simulated_PKPD.txt; the simulation 'truth' values used to
    # generate that dataset live in the (paywalled) Lestini 2015 paper.
    # The .mlxtran initialValues differ materially for Cl and kout
    # (Cl_init = 40, kout_init = 2 vs MLE = 9.91, 0.269), confirming those
    # are optimizer starting points, not the simulation truth or the
    # final-fit estimates.
    # ----------------------------------------------------------------------
    lka   <- log(1.97);  label("First-order absorption rate (ka, 1/h)")            # Output_simulated_PKPD.txt: ka = 1.97
    lcl   <- log(9.91);  label("Apparent clearance (CL/F, L/h)")                   # Output_simulated_PKPD.txt: Cl = 9.91
    lvc   <- log(95.1);  label("Apparent central volume of distribution (Vc/F, L)") # Output_simulated_PKPD.txt: V = 95.1
    lkout <- log(0.269); label("Biomarker turnover rate constant (kout, 1/h)")      # Output_simulated_PKPD.txt: kout = 0.269
    lc50  <- log(0.307); label("Concentration giving 50% of maximal TGF-beta inhibition (C50, mg/L)")  # Output_simulated_PKPD.txt: C50 = 0.307

    # IIV: Output_simulated_PKPD.txt reports omega values as standard deviations
    # on the log scale (Monolix globalSettings withVariance=no). nlmixr2's eta ini
    # value is the variance on the log scale, so use omega^2.
    # ka has no IIV in the .mlxtran INDIVIDUAL block (`ka = {distribution=logNormal, iiv=no}`).
    etalvc   ~ 0.553536  # omega_V    = 0.744 -> 0.744^2  (Output_simulated_PKPD.txt)
    etalcl   ~ 0.670761  # omega_Cl   = 0.819 -> 0.819^2  (Output_simulated_PKPD.txt)
    etalkout ~ 0.498436  # omega_kout = 0.706 -> 0.706^2  (Output_simulated_PKPD.txt)
    etalc50  ~ 0.675684  # omega_C50  = 0.822 -> 0.822^2  (Output_simulated_PKPD.txt)

    # Residual error per the .mlxtran OBSERVATIONS block:
    #   y1 = {prediction=Cc, error=proportional}  -> proportional on Cc (parent: bare propSd)
    #   y2 = {prediction=E,  error=constant}      -> additive on E (non-parent: addSd_E)
    propSd   <- 0.199; label("Proportional residual error on Cc (fraction)")              # Output_simulated_PKPD.txt: b_1 = 0.199
    addSd_E  <- 0.196; label("Additive residual error on E (fraction inhibition; 0-1 scale)") # Output_simulated_PKPD.txt: a_2 = 0.196
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual PK and PD parameters.
    #    ka has no IIV (.mlxtran INDIVIDUAL block); all other lognormal.
    # ------------------------------------------------------------------
    ka   <- exp(lka)
    cl   <- exp(lcl   + etalcl)
    vc   <- exp(lvc   + etalvc)
    kout <- exp(lkout + etalkout)
    c50  <- exp(lc50  + etalc50)

    kel <- cl / vc

    # ------------------------------------------------------------------
    # 2. PK ODE system (one-compartment, first-order absorption).
    #    pkpd_model.txt: Cc = pkmodel(ka, V, Cl)
    # ------------------------------------------------------------------
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central
    Cc <- central / vc

    # ------------------------------------------------------------------
    # 3. Indirect-response biomarker (E = fractional inhibition of
    #    TGF-beta signalling). pkpd_model.txt: E_0 = 0; Imax = 1;
    #    dE/dt = kout * Imax * Cc / (Cc + C50) - kout * E.
    #    With Imax = 1, the production term saturates at kout and the
    #    steady-state under sustained Cc -> infinity is E -> 1.
    # ------------------------------------------------------------------
    d/dt(effect) <- kout * Cc / (Cc + c50) - kout * effect
    effect(0)    <- 0
    E <- effect

    # ------------------------------------------------------------------
    # 4. Observation and error model (.mlxtran OBSERVATIONS block).
    # ------------------------------------------------------------------
    Cc ~ prop(propSd)
    E  ~ add(addSd_E)
  })
}
