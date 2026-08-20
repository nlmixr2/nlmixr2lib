Dodds_2021_covid19_viralkinetic <- function() {
  description <- paste(
    "Mechanistic viral kinetic (target-cell-limited with eclipse phase)",
    "model of SARS-CoV-2 infection dynamics adapted by Dodds et al. (2021)",
    "from the target-cell-limited viral kinetics framework of Goncalves et al.",
    "(2020) that characterized nasopharyngeal viral load in 13 hospitalized",
    "COVID-19 patients. Four ODE states: uninfected epithelial target cells",
    "'target' (T; #/mL), latently (eclipse) infected cells 'eclipse' (I1;",
    "#/mL), productively infected cells 'infected' (I2; #/mL) and free",
    "virions 'virus' (V; copies/mL). Rate constants (all /day) are: k for",
    "eclipse to productive transition, delta for productively-infected cell",
    "death, rho for virion production per productively-infected cell, c for",
    "virion clearance, and beta (mL/#/day) for target-cell infection.",
    "beta is derived from the basic reproductive ratio R0 via",
    "beta = R0 * delta * c / (T0 * (rho - R0 * delta)). All structural",
    "parameters are FIXED to the consensus values reported in Dodds 2021",
    "Section 2.1 (T0 = 1.33e5 #/mL, V0 = 0.1 #/mL, k = 3 /day,",
    "delta = 0.60 /day, rho = 22.7 /day, c = 10 /day, R0 = 8.6).",
    "Simulation-time drug interventions are supported as Heaviside step",
    "modifiers applied at time t_intervention: fractional inhibition on",
    "beta, k, and rho (imax_*, 0-1) and multiplicative stimulation on",
    "delta and c (smax_*, >=0), all defaulting to fixed(0) so the untreated",
    "natural-history dynamics are the baseline behaviour. The model reproduces",
    "the untreated viral trajectory peaking ~9 days post-infection used by",
    "Dodds 2021 to explore single- and combination-target antiviral effects on",
    "viral load AUC, duration of viral shedding (>= 100 copies/mL) and",
    "epithelial cells infected. This is a deterministic simulation model",
    "(no IIV, no residual error, no PK compartment): drug potencies must be",
    "supplied per scenario by overriding imax_*, smax_*, and t_intervention.",
    "The eclipse-phase compartment is declared paper-specific via",
    "paper_specific_compartments; the productive compartment reuses the",
    "existing canonical 'infected' registered by Wang_2018_daclatasvir_asunaprevir.",
    sep = " "
  )

  reference <- paste(
    "Dodds MG, Krishna R, Goncalves A, Rayner CR. (2021).",
    "Model-informed drug repurposing: Viral kinetic modelling to prioritize",
    "rational drug combinations for COVID-19.",
    "Br J Clin Pharmacol 87(9):3439-3450.",
    "doi:10.1111/bcp.14486.",
    "Structural model and consensus parameter values are transcribed from",
    "the upstream Goncalves et al. (2020) analysis cited by Dodds 2021 as",
    "reference 16 (Timing of antiviral treatment initiation is critical to",
    "reduce SARS-CoV-2 viremia. CPT Pharmacometrics Syst Pharmacol",
    "9(9):509-514, doi:10.1002/psp4.12543).",
    sep = " "
  )
  vignette <- "Dodds_2021_covid19_viralkinetic"

  paper_specific_compartments <- c("eclipse")

  units <- list(
    time          = "day",
    dosing        = "no dosing (simulation-only viral kinetic model; drug effects supplied via imax_*, smax_* parameter overrides)",
    concentration = "virus copies/mL; target/eclipse/infected cells #/mL"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    target   = list(analyte = "uninfected epithelial target cells", units = NA_character_, specimen = "blood cell", verified = FALSE),
    eclipse  = list(analyte = "latently (eclipse) infected cells", units = NA_character_, specimen = "blood cell", verified = FALSE),
    infected = list(analyte = "productively infected cells", units = NA_character_, specimen = "blood cell", verified = FALSE),
    virus    = list(analyte = "free virions", units = NA_character_, specimen = "not applicable", verified = FALSE)
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 13L,
    n_studies      = 1L,
    age_range      = "Adults hospitalized with COVID-19 (specific ages not tabulated in Dodds 2021; upstream Goncalves 2020 studied 13 hospitalized patients from Zou et al. 2020 nasopharyngeal-swab dataset).",
    weight_range   = "Not reported.",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported.",
    disease_state  = "Hospitalized COVID-19 patients (upstream Goncalves 2020 dataset). Model is intended for simulation of natural-history and hypothetical antiviral intervention scenarios in the early stage of SARS-CoV-2 infection.",
    dose_range     = "No physical dosing (drug-effect step functions replace an explicit PK layer).",
    regimens       = "N/A -- simulation-only.",
    regions        = "Not reported (upstream Goncalves 2020 patients from China).",
    notes          = paste(
      "The population underlying the parameter estimation is described in",
      "Goncalves et al. (2020) CPT PSP 9(9):509-514, which fitted the",
      "target-cell-limited-with-eclipse model to 13 hospitalized COVID-19",
      "patients from Zou et al. (NEJM 2020). Dodds 2021 does NOT re-fit --",
      "it transcribes the Goncalves consensus (model-averaged) parameter",
      "values from Section 2.1 and uses them for scenario simulation of",
      "drug interventions targeting the five viral cell-cycle rate constants."
    )
  )

  ini({
    # =========================================================================
    # Structural viral kinetic parameters (all FIXED to Dodds 2021 Section 2.1
    # consensus values transcribed from the Goncalves 2020 model-averaging
    # analysis). All rate constants are in /day; state units are #/mL for
    # cells and copies/mL for virions. beta is derived inside model() from
    # R0, delta, c, T0 and rho per the Section 2.1 derivation
    # beta = R0 * delta * c / (T0 * (rho - R0 * delta)); we fix R0 rather
    # than beta so the reproductive-ratio interpretation is transparent.
    # =========================================================================

    lT0    <- fixed(log(1.33e5))
    label("Initial target-cell density T0 (#/mL)")                              # Dodds 2021 Section 2.1: T0 = 1.33E5 #/mL
    lV0    <- fixed(log(0.1))
    label("Initial free-virion concentration V0 (copies/mL)")                   # Dodds 2021 Section 2.1: V0 = 0.1 #/mL
    lk     <- fixed(log(3))
    label("Eclipse-to-productive transition rate k (1/day)")                    # Dodds 2021 Section 2.1: k = 3 1/d
    ldelta <- fixed(log(0.60))
    label("Productively-infected cell death rate delta (1/day)")                # Dodds 2021 Section 2.1: delta = 0.60 1/d
    lrho   <- fixed(log(22.7))
    label("Virion production rate rho (virions/cell/day)")                      # Dodds 2021 Section 2.1: rho = 22.7 1/d
    lc     <- fixed(log(10))
    label("Free-virion clearance rate c (1/day)")                               # Dodds 2021 Section 2.1: c = 10 1/d
    lR0    <- fixed(log(8.6))
    label("Basic reproductive ratio R0 (unitless)")                             # Dodds 2021 Section 2.1: R0 = 8.6 (used to derive beta)

    # =========================================================================
    # Intervention parameters (Dodds 2021 Section 2.2). Interventions are
    # modelled as Heaviside step functions active at t >= t_intervention:
    #   beta_eff  = beta  * (1 - imax_beta  * step)
    #   k_eff     = k     * (1 - imax_k     * step)
    #   rho_eff   = rho   * (1 - imax_rho   * step)
    #   delta_eff = delta * (1 + smax_delta * step)
    #   c_eff     = c     * (1 + smax_c     * step)
    # imax_* are fractional inhibitions in [0, 1]; smax_* are multiplicative
    # stimulations in [0, +Inf). All default to fixed(0) so the untreated
    # natural-history trajectory is the baseline. The user overrides these
    # (e.g. via rxode2 params or a per-subject params data frame) to
    # explore specific intervention scenarios matching Dodds 2021 Figures
    # 2-4. t_intervention defaults to a value beyond any reasonable
    # simulation window so that omitting overrides yields the natural-
    # history dynamics; supply t_intervention = 3 (peak - 6 days),
    # 6 (peak - 3 days), 9 (peak), 12 (peak + 3 days) or 15 (peak + 6 days)
    # per Dodds 2021 Section 2.2 to reproduce the paper's timings.
    # =========================================================================

    imax_beta  <- fixed(0)
    label("Fractional inhibition of beta at t >= t_intervention (0-1); default = 0 (no effect)")  # Dodds 2021 Section 2.2 (inhibitory function on beta)
    imax_k     <- fixed(0)
    label("Fractional inhibition of k at t >= t_intervention (0-1); default = 0 (no effect)")     # Dodds 2021 Section 2.2 (inhibitory function on k)
    imax_rho   <- fixed(0)
    label("Fractional inhibition of rho at t >= t_intervention (0-1); default = 0 (no effect)")   # Dodds 2021 Section 2.2 (inhibitory function on rho)
    smax_delta <- fixed(0)
    label("Multiplicative stimulation of delta at t >= t_intervention (>= 0); default = 0")       # Dodds 2021 Section 2.2 (stimulatory function on delta)
    smax_c     <- fixed(0)
    label("Multiplicative stimulation of c at t >= t_intervention (>= 0); default = 0")           # Dodds 2021 Section 2.2 (stimulatory function on c)
    t_intervention <- fixed(1e6)
    label("Intervention start time (day); default = 1e6 (never during a normal simulation window)")  # Dodds 2021 Section 2.2 (Heaviside onset time)
  })

  model({
    # =========================================================================
    # 1. Back-transform typical values (all FIXED so etas are absent).
    # =========================================================================
    T0    <- exp(lT0)
    V0    <- exp(lV0)
    k     <- exp(lk)
    delta <- exp(ldelta)
    rho   <- exp(lrho)
    c     <- exp(lc)
    R0    <- exp(lR0)

    # =========================================================================
    # 2. Derived infection-rate constant beta (mL/#/day) from R0 (Dodds 2021
    # Section 2.1 derivation: beta = R0 * delta * c / (T0 * (rho - R0 *
    # delta))). Substituting the fixed values above:
    #   beta = 8.6 * 0.60 * 10 / (1.33e5 * (22.7 - 8.6 * 0.60))
    #        = 51.6 / (1.33e5 * 17.54)
    #        = 2.21e-5 mL/#/day
    # which matches the beta value reported in Section 2.1.
    # =========================================================================
    beta <- R0 * delta * c / (T0 * (rho - R0 * delta))

    # =========================================================================
    # 3. Heaviside step function for the intervention onset (Dodds 2021
    # Section 2.2). At t >= t_intervention the intervention effects are
    # switched on; before then, all modifiers are zero. Encoded as an
    # arithmetic (t >= t_intervention) rather than a conditional so the
    # ODE right-hand side remains smooth for the solver.
    # =========================================================================
    step_active <- (t >= t_intervention)

    beta_eff  <- beta  * (1 - imax_beta  * step_active)
    k_eff     <- k     * (1 - imax_k     * step_active)
    rho_eff   <- rho   * (1 - imax_rho   * step_active)
    delta_eff <- delta * (1 + smax_delta * step_active)
    c_eff     <- c     * (1 + smax_c     * step_active)

    # =========================================================================
    # 4. Target-cell-limited-with-eclipse ODE system (Dodds 2021 Section
    # 2.1 prose + Figure 1). Standard Baccam/Goncalves formulation with
    # four states.
    #   dT /dt  = -beta * T * V         (target cells infected by virus)
    #   dI1/dt  =  beta * T * V - k * I1  (eclipse cells acquired then
    #                                      mature into productive cells)
    #   dI2/dt  =  k * I1 - delta * I2   (productive cells produced by
    #                                     eclipse maturation, die at delta)
    #   dV /dt  =  rho * I2 - c * V      (virions produced by productive
    #                                     cells, cleared at rate c)
    # =========================================================================
    d/dt(target)   <- -beta_eff * virus * target
    d/dt(eclipse)  <-  beta_eff * virus * target - k_eff * eclipse
    d/dt(infected) <-  k_eff * eclipse - delta_eff * infected
    d/dt(virus)    <-  rho_eff * infected - c_eff * virus

    # =========================================================================
    # 5. Initial conditions (Dodds 2021 Section 2.1). At t = 0 (infection
    # onset) all epithelial target cells are uninfected, no eclipse or
    # productively-infected cells exist, and the small V0 = 0.1 copies/mL
    # inoculum seeds the infection.
    # =========================================================================
    target(0)   <- T0
    eclipse(0)  <- 0
    infected(0) <- 0
    virus(0)    <- V0

    # =========================================================================
    # 6. Observation outputs. Vlog10 is the log10 free-virion concentration
    # (log10 copies/mL); the 1e-12 floor prevents -Inf when virus -> 0
    # after eradication. cellsInfectedCum is the cumulative loss of
    # target cells from infection (T0 - target), the third metric in
    # Dodds 2021 Table 1 ("epithelial cells infected"; Section 2.4.3).
    # No residual error is defined: this is a deterministic simulation
    # model with all parameters fixed to consensus values.
    # =========================================================================
    Vlog10           <- log10(virus + 1e-12)
    cellsInfectedCum <- T0 - target
  })
}
