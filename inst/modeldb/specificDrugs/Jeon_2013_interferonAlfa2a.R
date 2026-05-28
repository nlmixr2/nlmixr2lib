Jeon_2013_interferonAlfa2a <- function() {
  description <- paste(
    "Joint PK-PD model for a sustained-release subcutaneous formulation of",
    "interferon alfa-2a (SR-IFN-alpha) and the serum neopterin response in",
    "healthy adult male volunteers (Jeon 2013).",
    "Pharmacokinetics: one-compartment with first-order elimination and a",
    "parallel mixture of zero- and first-order absorption. A fraction",
    "Fz = exp(RF)/(1 + exp(RF)) of the dose is absorbed by a zero-order",
    "process with duration D2 entering the central compartment directly;",
    "the remaining 1 - Fz is absorbed by a first-order process (rate Ka)",
    "from a depot compartment with lag time ALAG, accounting for the",
    "second concentration peak observed around 100 h post-injection.",
    "Pharmacodynamics: indirect-response (turnover) model for serum",
    "neopterin (baseline BASE = Kin/Kout) with a single transit",
    "compartment placed between the stimulus and the observed neopterin,",
    "delaying the neopterin response through mean transit time MTT. The",
    "drug stimulates the zero-order production rate of neopterin through",
    "a sigmoid Emax function E(C) = EMAX * C^GA / (EC50(t)^GA + C^GA),",
    "where EC50 is time-dependent and increases monotonically over time",
    "as EC50(t) = ECB * (1 + CA * (1 - exp(-CB * t))) -- an empirical",
    "saturation device that captures the observed loss of the neopterin",
    "dose-response between groups (9, 18, 27, 36 MIU) over the 0-264 h",
    "observation window. No covariate effects were retained in the final",
    "model. Doses are entered in MIU (10^6 IU); the published apparent",
    "clearance (CL/F = 12.2 L/h) and apparent volume of distribution",
    "(V/F = 691 L) match values previously reported for IFN-alpha in",
    "healthy subjects (Reference [19] of Jeon 2013). The model uses an",
    "explicit specific-activity conversion (1 MIU = 4 ug = 4e6 pg, from",
    "the WHO IFN-alpha-2a International Standard at 2.5e8 IU/mg) so the",
    "doses in user data can be entered in MIU and the simulated Cc is",
    "returned in pg/mL. The specific-activity conversion is not stated",
    "in the paper itself; it is documented in the validation vignette's",
    "Assumptions and deviations section."
  )
  reference <- paste(
    "Jeon S, Juhn JH, Han S, Lee J, Hong T, Paek J, Yim DS.",
    "Saturable human neopterin response to interferon-alpha assessed by a",
    "pharmacokinetic-pharmacodynamic model.",
    "J Transl Med. 2013;11:240. doi:10.1186/1479-5876-11-240"
  )
  vignette <- "Jeon_2013_interferonAlfa2a"
  units <- list(time = "h", dosing = "MIU", concentration = "pg/mL")

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 24L,
    n_studies      = 1L,
    age_range      = "18-43 years",
    age_median     = "approximately 22 years (median of group medians 21.5, 22, 22, 24)",
    weight_range   = "53.8-98.8 kg",
    weight_median  = "approximately 75 kg (median of group medians 78.65, 76.15, 74.40, 74.40)",
    height_range   = "168.0-193.0 cm",
    sex_female_pct = 0,
    race_ethnicity = "Not explicitly reported; study conducted at Kendle International BV, Utrecht, Netherlands.",
    disease_state  = "Healthy adult male volunteers, BMI 19-29 kg/m^2, no clinically relevant conditions on medical history, physical examination, lab tests, or ECG; subjects with possible alteration in IFN-alpha metabolism or hypersensitivity to IFN-alpha were excluded.",
    dose_range     = "Single subcutaneous dose of SR-IFN-alpha at 9, 18, 27, or 36 MIU (six subjects per dose group). An 8-subject active-control group received 3 MIU Roferon-A (immediate-release IFN-alpha-2a) and was excluded from the PK-PD model build.",
    regions        = "Netherlands (single-center Phase I; Kendle International BV, Utrecht).",
    notes          = "Randomized, double-blind, active-controlled, dose-escalation Phase I clinical study; PK sampling at pre-dose and 0.75, 1.5, 3, 6, 8, 10, 12, 18, 24, 30, 36, 48, 60, 72, 96, 120, 144, 168, 192 h post-injection; PD (neopterin) sampling at pre-dose and 3, 8, 12, 18, 24, 36, 48, 72, 96, 120, 144, 168, 192, 264 h post-injection. Demographics summarised in Jeon 2013 Table 1; only the four SR-IFN-alpha groups (n = 24) contributed to the PK-PD model."
  )

  ini({
    # ---------------------------------------------------------------
    # PK STRUCTURAL PARAMETERS - Jeon 2013 Table 4, "Pharmacokinetics" rows
    # ---------------------------------------------------------------
    lcl  <- log(12.2)    ; label("Apparent clearance CL/F (L/h)")                                                       # Jeon 2013 Table 4: CL/F = 12.2 L/h
    lvc  <- log(691)     ; label("Apparent central volume of distribution V/F (L)")                                     # Jeon 2013 Table 4: V/F = 691 L
    ld2  <- log(20.2)    ; label("Duration of zero-order absorption D2 (h)")                                            # Jeon 2013 Table 4: D2 = 20.2 h
    lka  <- log(0.00653) ; label("First-order absorption rate constant Ka (1/h)")                                       # Jeon 2013 Table 4: Ka = 0.00653 1/h
    ltlag <- log(85.7)   ; label("Lag time to start of first-order absorption ALAG (h)")                                # Jeon 2013 Table 4: ALAG = 85.7 h
    lrf  <- log(0.185)   ; label("Logit parameter RF for zero-order absorbed fraction; Fz = exp(RF)/(1 + exp(RF))")     # Jeon 2013 Table 4: RF = 0.185 (Fz = e^RF / (1 + e^RF) per Table 4 footnote b)

    # ---------------------------------------------------------------
    # PD STRUCTURAL PARAMETERS - Jeon 2013 Table 4, "Pharmacodynamics" rows
    # ---------------------------------------------------------------
    lbase <- log(5.85)   ; label("Baseline serum neopterin BASE (nmol/L)")                                              # Jeon 2013 Table 4: BASE = 5.85 nmol/L
    lkout <- log(0.0311) ; label("First-order elimination rate of serum neopterin Kout (1/h)")                          # Jeon 2013 Table 4: Kout = 0.0311 1/h
    lemax <- log(16.1)   ; label("Maximum stimulatory effect on neopterin production EMAX (unitless multiplier on Kin)") # Jeon 2013 Table 4: EMAX = 16.1
    lga   <- log(1.24)   ; label("Hill coefficient GA on the sigmoid stimulation function (unitless)")                  # Jeon 2013 Table 4: GA = 1.24
    lca   <- log(405)    ; label("Coefficient CA in EC50(t) = ECB * (1 + CA * (1 - exp(-CB*t))) (unitless)")            # Jeon 2013 Table 4: CA = 405
    lcb   <- log(0.0068) ; label("Coefficient CB in EC50(t) = ECB * (1 + CA * (1 - exp(-CB*t))) (1/h)")                 # Jeon 2013 Table 4: CB = 0.0068 1/h
    lecb  <- log(2.17)   ; label("Baseline of EC50 ECB (pg/mL)")                                                        # Jeon 2013 Table 4: ECB = 2.17 (same units as Cc, pg/mL)
    lmtt  <- log(14.6)   ; label("Mean transit time MTT (h)")                                                           # Jeon 2013 Table 4: MTT = 14.6 h

    # ---------------------------------------------------------------
    # INTER-INDIVIDUAL VARIABILITY
    # ---------------------------------------------------------------
    # Jeon 2013 assumed log-normal IIV with Pj = TVP * exp(eta_j) (Methods,
    # "A log normal distribution was assumed for interindividual variability").
    # IIV is reported in Table 4 as percent CV; the on-log-scale variance
    # is omega^2 = log(CV^2 + 1).
    #
    # PK IIV: CL, V, D2, RF, Ka (no IIV reported on ALAG).
    #   omega^2_CL  = log(0.261^2 + 1) = 0.06598
    #   omega^2_V   = log(0.238^2 + 1) = 0.05504
    #   omega^2_D2  = log(0.357^2 + 1) = 0.12001
    #   omega^2_RF  = log(0.347^2 + 1) = 0.11373
    #   omega^2_Ka  = log(0.770^2 + 1) = 0.45203
    etalcl ~ 0.06598                                                                                                    # Jeon 2013 Table 4: omega_CL = 26.1% -> omega^2 = log(0.261^2 + 1)
    etalvc ~ 0.05504                                                                                                    # Jeon 2013 Table 4: omega_V = 23.8% -> omega^2 = log(0.238^2 + 1)
    etald2 ~ 0.12001                                                                                                    # Jeon 2013 Table 4: omega_D2 = 35.7% -> omega^2 = log(0.357^2 + 1)
    etalrf ~ 0.11373                                                                                                    # Jeon 2013 Table 4: omega_RF = 34.7% -> omega^2 = log(0.347^2 + 1)
    etalka ~ 0.45203                                                                                                    # Jeon 2013 Table 4: omega_Ka = 77.0% -> omega^2 = log(0.770^2 + 1)

    # PD IIV: BASE, CB, GA, ECB, MTT (no IIV reported on Kout, EMAX, CA).
    #   omega^2_BASE = log(0.1385^2 + 1) = 0.01906
    #   omega^2_CB   = log(0.5731^2 + 1) = 0.28367
    #   omega^2_GA   = log(0.1351^2 + 1) = 0.01814
    #   omega^2_ECB  = log(0.2131^2 + 1) = 0.04445
    #   omega^2_MTT  = log(0.1336^2 + 1) = 0.01773
    etalbase ~ 0.01906                                                                                                  # Jeon 2013 Table 4: omega_BASE = 13.85% -> omega^2 = log(0.1385^2 + 1)
    etalcb   ~ 0.28367                                                                                                  # Jeon 2013 Table 4: omega_CB = 57.31% -> omega^2 = log(0.5731^2 + 1)
    etalga   ~ 0.01814                                                                                                  # Jeon 2013 Table 4: omega_GA = 13.51% -> omega^2 = log(0.1351^2 + 1)
    etalecb  ~ 0.04445                                                                                                  # Jeon 2013 Table 4: omega_ECB = 21.31% -> omega^2 = log(0.2131^2 + 1)
    etalmtt  ~ 0.01773                                                                                                  # Jeon 2013 Table 4: omega_MTT = 13.36% -> omega^2 = log(0.1336^2 + 1)

    # ---------------------------------------------------------------
    # RESIDUAL ERROR
    # ---------------------------------------------------------------
    # PK (Cc, IFN-alpha): combined additive + proportional on the linear scale.
    addSd       <- 3.92  ; label("Additive residual SD on Cc (pg/mL)")                                                  # Jeon 2013 Table 4: sigma_add = 3.92 pg/mL
    propSd      <- 0.078 ; label("Proportional residual SD on Cc (fraction)")                                           # Jeon 2013 Table 4: sigma_prop = 7.8%
    # PD (Cneop, neopterin): additive only on the linear scale.
    addSd_Cneop <- 1.14  ; label("Additive residual SD on Cneop (nmol/L)")                                              # Jeon 2013 Table 4: sigma_add (neopterin) = 1.14 nmol/L
  })

  model({
    # ---------------------------------------------------------------
    # IFN-alpha-2a specific-activity conversion (non-paper provenance):
    #   1 MIU = 1e6 IU = (1e6 / 2.5e8) mg = 4 ug = 4e6 pg of IFN-alpha
    # based on the WHO IFN-alpha-2a International Standard nominal
    # specific activity of 2.5e8 IU/mg (Roferon-A package insert; WHO
    # Expert Committee on Biological Standardization). The conversion is
    # required to reconcile the published apparent CL/F = 12.2 L/h and
    # V/F = 691 L (consistent with literature IFN-alpha values; see Jeon
    # 2013 Discussion citing reference [19]) with the published Cmax /
    # AUC values reported in pg/mL, so that the user can enter doses in
    # MIU and receive a Cc in pg/mL. See vignette Assumptions and
    # deviations for the audit trail.
    miu_to_pg <- 4e6   # pg per MIU
    # Cc (pg/mL) = central (MIU) * miu_to_pg (pg/MIU) / vc (L) / 1000 (mL/L)
    #            = central / vc * (miu_to_pg / 1000)
    miu_to_pg_per_mL <- miu_to_pg / 1000  # = 4000 pg/mL per MIU/L

    # ---------------------------------------------------------------
    # Individual PK parameters (log-normal IIV)
    # ---------------------------------------------------------------
    cl   <- exp(lcl   + etalcl)
    vc   <- exp(lvc   + etalvc)
    d2   <- exp(ld2   + etald2)
    ka   <- exp(lka   + etalka)
    alag <- exp(ltlag)
    rf   <- exp(lrf   + etalrf)
    # Fz: zero-order absorbed fraction via logit transform of RF
    # (Jeon 2013 Table 4 footnote b: Fz = e^RF / (1 + e^RF)).
    fz   <- exp(rf) / (1 + exp(rf))

    kel <- cl / vc

    # ---------------------------------------------------------------
    # Individual PD parameters (log-normal IIV where reported; typical-
    # value-only for Kout, EMAX, CA per Table 4)
    # ---------------------------------------------------------------
    base <- exp(lbase + etalbase)
    kout <- exp(lkout)
    emax <- exp(lemax)
    ga   <- exp(lga   + etalga)
    ca   <- exp(lca)
    cb   <- exp(lcb   + etalcb)
    ecb  <- exp(lecb  + etalecb)
    mtt  <- exp(lmtt  + etalmtt)

    # At baseline (no drug) the turnover system holds neopterin = BASE
    # by setting Kin = BASE * Kout. With a single transit compartment
    # placed between Kin and the observed neopterin, the mean transit
    # time MTT is interpreted as 1 / Ktr (residence time in the transit
    # compartment); steady state for the transit state is Kin / Ktr.
    kin <- base * kout
    ktr <- 1 / mtt

    # ---------------------------------------------------------------
    # Time-dependent EC50 (Jeon 2013, EC50(t) equation):
    #   EC50(t) = ECB * (1 + CA * (1 - exp(-CB * t)))
    # CDF-of-exponential shape that raises EC50 from ECB at t = 0 toward
    # an asymptote of ECB * (1 + CA) as t -> infinity.
    # ---------------------------------------------------------------
    ec50t <- ecb * (1 + ca * (1 - exp(-cb * t)))

    # ---------------------------------------------------------------
    # ODE system
    # ---------------------------------------------------------------
    # PK: depot drains first-order to central; central also receives the
    # zero-order absorbed fraction directly via dur(central) below.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # IFN-alpha plasma concentration (pg/mL), reported in Jeon 2013
    # Tables 2 and 4. The specific-activity conversion above maps the
    # internal central-compartment amount (MIU) to pg/mL.
    Cc   <- central / vc * miu_to_pg_per_mL

    # IFN-alpha effect on neopterin production (sigmoid Emax on Kin).
    # ECB is reported on the same pg/mL scale as Cc (Jeon 2013 Table 4
    # baseline of EC50), so the Hill expression uses Cc (pg/mL) directly:
    #   stim = EMAX * Cc^GA / (EC50(t)^GA + Cc^GA)
    stim <- emax * Cc^ga / (ec50t^ga + Cc^ga)

    # Neopterin turnover with one transit compartment delaying the signal.
    # Paper compartments A(3) = transit, A(4) = neopterin observation; the
    # canonical nlmixr2lib names are `transit1` (chain index 1) for the
    # transit state and `effect` for the observed biomarker compartment.
    d/dt(transit1) <-  kin * (1 + stim) - ktr * transit1
    d/dt(effect)   <-  ktr * transit1    - kout * effect

    # Initial conditions: pre-dose steady state of the turnover system.
    # dT/dt = Kin - Ktr*T = 0    -> T_ss   = Kin/Ktr = base * Kout * MTT
    # dE/dt = Ktr*T - Kout*E = 0 -> E_ss   = Kin/Kout = base
    transit1(0) <- base * kout * mtt
    effect(0)   <- base

    # ---------------------------------------------------------------
    # Absorption controls
    # ---------------------------------------------------------------
    # The dose must target both `depot` (first-order arm, fraction
    # 1 - Fz, lag ALAG) and `central` (zero-order arm, fraction Fz,
    # duration D2) -- one dose-event row per administration per
    # compartment (see vignette for an example event table).
    f(depot)     <- 1 - fz
    alag(depot)  <- alag

    f(central)   <- fz
    dur(central) <- d2

    # ---------------------------------------------------------------
    # Observations and residual error
    # ---------------------------------------------------------------
    # Serum neopterin observation derived from the effect compartment
    # state; baseline serves as the pre-dose concentration (nmol/L).
    Cneop <- effect

    Cc    ~ add(addSd) + prop(propSd)
    Cneop ~ add(addSd_Cneop)
  })
}
