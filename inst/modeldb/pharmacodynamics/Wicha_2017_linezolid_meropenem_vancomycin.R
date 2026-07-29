Wicha_2017_linezolid_meropenem_vancomycin <- function() {
  description <- "In vitro (MSSA ATCC 29213). Semimechanistic time-kill pharmacodynamic model of linezolid, meropenem, and vancomycin against methicillin-susceptible Staphylococcus aureus. Bacterial life cycle has three states: growing (gro), replicating (repl), and persisting (pers). LZD inhibits the GRO->REP transition (bacteriostatic via krep) and induces a replication-independent killing rate kdeath_lzd on growing bacteria. MER and VAN, as cell wall-active antibiotics, impair successful doubling at the REP->GRO transition; the joint MER+VAN action is encoded as a modified Bliss-independence term that includes the paradoxical Eagle-effect self-inhibition of MER at high concentrations and the VAN Emax cap. Drug-unsusceptible persisters are generated during replication at rates kper_mer * E_MER and kper_van * E_VAN, then die at kdeath_per. An adaptive-resistance submodel (Tam 2005) inflates the effective EC50 of MER and of VAN over time via fractional ARon states; subinhibitory VAN concentrations inhibit the MER-adaption rate (monodirectional VAN-on-MER PD interaction). MER and VAN solution concentrations decay first-order due to chemical degradation in growth medium (rates fixed from HPLC measurement); LZD is stable. The model is in-vitro PD only -- there is no human PK component; drug exposures are static dosing at t = 0. Random effects (eta) are NOT present: the paper reports replicate-only experimental variability and uses an additive residual error on log10(CFU/mL)."
  reference <- "Wicha SG, Huisinga W, Kloft C. Translational pharmacometric evaluation of typical antibiotic broad-spectrum combination therapies against Staphylococcus aureus exploiting in vitro information. CPT Pharmacometrics Syst Pharmacol. 2017;6(8):512-522. doi:10.1002/psp4.12197."
  vignette <- "Wicha_2017_linezolid_meropenem_vancomycin"
  units <- list(time = "hour", dosing = "mg/L (initial concentration)", concentration = "log10 CFU/mL (observation); mg/L (drug states)")

  # The model has no covariates -- it is a typical-value in-vitro PD model with
  # static drug exposures. Drug initial concentrations are applied via dosing
  # events into the `lzd`, `mer`, and `van` compartments at time 0 with `amt`
  # interpreted as mg/L (consistent with the in-vitro convention where the
  # compartment state IS the bath concentration). Paper-studied ranges per
  # Wicha 2017 Methods: LZD 0.5-32 mg/L, MER 0.015-8 mg/L, VAN 0.06-16 mg/L.
  # MER and VAN decay first-order in the growth medium (kdeg fixed from HPLC
  # measurements); LZD is chemically stable.
  covariateData <- list()

  population <- list(
    species             = "in vitro (Staphylococcus aureus, ATCC 29213 reference strain)",
    n_subjects          = NA_integer_,
    n_studies           = 1L,
    organism            = "Staphylococcus aureus ATCC 29213 (methicillin-susceptible reference strain); model also externally evaluated against clinical MSSA isolates MV13391 and MV13488",
    system              = "Static time-kill experiments, 24-hour exposure, replicates >= 2 per concentration tier, dense sampling (>= 8 timepoints)",
    medium              = "Cation-adjusted Mueller-Hinton broth (per Wicha 2015 upstream dataset, ref 6 in the paper)",
    temperature         = "37 C (standard bacteriology incubation)",
    duration            = "24 h",
    mic_values          = c(LZD = "2.0 mg/L", MER = "0.125 mg/L", VAN = "1.0 mg/L"),
    concentration_range = c(LZD = "0.5-32 mg/L", MER = "0.015-8 mg/L", VAN = "0.06-16 mg/L"),
    regimens            = "Monotherapy of LZD, MER, or VAN; dual combinations of LZD+MER and VAN+MER spanning the clinically relevant concentration range; plus antibiotic-free growth controls.",
    notes               = "In-vitro pharmacodynamic study; no human or animal subjects. n = 1,617 timed CFU data points were available for model building. Final fit reported additive residual variability on log10(CFU/mL) only -- no IIV / between-experiment etas were required. Bootstrap n = 1,198 was used for 95% CIs. See Wicha 2017 Table 1 footnote and Methods (page 2-3)."
  )

  ini({
    # ---- Bacterial life-cycle ----
    cfu0 <- 6.06
    label("Initial bacterial concentration (log10 CFU/mL)")  # Wicha 2017 Table 1, Parameters of the bacterial life-cycle
    cfumax <- 9.43
    label("Maximum bacterial carrying capacity (log10 CFU/mL)")  # Wicha 2017 Table 1
    lklag <- log(0.88)
    label("First-order delay rate constant for lag-phase exit (klag, 1/h)")  # Wicha 2017 Table 1
    lkrep <- log(1.56)
    label("Transit rate constant from gro to repl (krep, 1/h)")  # Wicha 2017 Table 1
    lkdoub <- fixed(log(100))
    label("Replication / doubling rate constant (kdoub, 1/h; FIXED -- not rate-limiting)")  # Wicha 2017 Table 1 (FIXED)
    lkdeath_per <- log(0.23)
    label("Basal death rate of persisters (kdeath_per, 1/h)")  # Wicha 2017 Table 1

    # ---- LZD drug parameters ----
    lec50_lzd <- log(0.68)
    label("LZD half-maximum-effect concentration on krep and kdeath_lzd (EC50_LZD, mg/L)")  # Wicha 2017 Table 1, Drug-related parameters
    lhill_lzd <- log(1.55)
    label("LZD Hill coefficient (H_LZD, unitless)")  # Wicha 2017 Table 1
    lkdeath_lzd <- log(0.10)
    label("Basal death rate of growth-arrested bacteria induced by LZD (kdeath_LZD, 1/h)")  # Wicha 2017 Table 1

    # ---- MER drug parameters (initial-killing arm) ----
    lec50_mer_t0 <- log(0.022)
    label("MER half-maximum-effect concentration on kdoub and kper_mer at t = 0 (EC50_MER,t=0, mg/L)")  # Wicha 2017 Table 1
    lhill_mer <- log(3.23)
    label("MER Hill coefficient for initial killing (H_MER, unitless)")  # Wicha 2017 Table 1

    # ---- MER paradoxical Eagle-effect arm ----
    emax_mer_eagle <- 0.328
    label("Eagle-effect maximum fractional reduction of MER effect at high MER (Emax_MER,Eagle, unitless)")  # Wicha 2017 Table 1 -- "32.8%"
    lec50_mer_eagle <- log(1.35)
    label("MER half-maximum paradoxical Eagle-effect concentration (EC50_MER,Eagle, mg/L)")  # Wicha 2017 Table 1
    lhill_mer_eagle <- fixed(log(4))
    label("MER Eagle-effect Hill coefficient (H_MER,Eagle, unitless; FIXED)")  # Wicha 2017 Table 1 (FIXED)

    # ---- MER adaption submodel ----
    lb_mer <- log(9.53)
    label("MER adaption magnitude factor (b_MER; max adapted EC50 = (1 + b_MER) * EC50_MER,t=0)")  # Wicha 2017 Table 1
    ls_mer <- log(0.47)
    label("MER adaption second-order delay rate constant (s_MER, L/(mg*h))")  # Wicha 2017 Table 1

    # ---- MER persister generation ----
    lkper_mer <- log(0.11)
    label("MER-driven persister development rate (kper_MER, 1/h)")  # Wicha 2017 Table 1

    # ---- MER chemical degradation in growth medium (HPLC-determined; FIXED) ----
    lkdeg_mer <- fixed(log(0.019))
    label("MER first-order degradation rate in growth medium (kdeg_MER, 1/h; FIXED -- HPLC)")  # Wicha 2017 Table 1 (FIXED)

    # ---- VAN drug parameters (initial-killing arm) ----
    emax_van <- 0.743
    label("VAN maximum fractional reduction of successful doubling (Emax_VAN, unitless)")  # Wicha 2017 Table 1 -- "74.3%"
    lec50_van_t0 <- log(0.46)
    label("VAN half-maximum-effect concentration on kdoub and kper_van at t = 0 (EC50_VAN,t=0, mg/L)")  # Wicha 2017 Table 1
    lhill_van <- fixed(log(20))
    label("VAN Hill coefficient (H_VAN, unitless; FIXED -- steep on/off shape)")  # Wicha 2017 Table 1 (FIXED)

    # ---- VAN-on-MER adaption inhibition (PD interaction) ----
    lec50_van_ari <- log(0.39)
    label("VAN half-maximum suppression of MSSA adaption to MER (EC50_VAN,ARI, mg/L)")  # Wicha 2017 Table 1
    lhill_van_ari <- fixed(log(1))
    label("VAN-on-MER-adaption Hill coefficient (H_VAN,ARI, unitless; FIXED)")  # Wicha 2017 Table 1 (FIXED)

    # ---- VAN adaption submodel ----
    lb_van <- log(3.59)
    label("VAN adaption magnitude factor (b_VAN; max adapted EC50 = (1 + b_VAN) * EC50_VAN,t=0)")  # Wicha 2017 Table 1
    ls_van <- log(0.034)
    label("VAN adaption second-order delay rate constant (s_VAN, L/(mg*h))")  # Wicha 2017 Table 1

    # ---- VAN persister generation ----
    lkper_van <- log(0.017)
    label("VAN-driven persister development rate (kper_VAN, 1/h)")  # Wicha 2017 Table 1

    # ---- VAN chemical degradation in growth medium (HPLC-determined; FIXED) ----
    lkdeg_van <- fixed(log(0.0039))
    label("VAN first-order degradation rate in growth medium (kdeg_VAN, 1/h; FIXED -- HPLC)")  # Wicha 2017 Table 1 (FIXED) -- "3.9e-03"

    # ---- Residual error ----
    # Wicha 2017 Table 1, "r [log10 CFU/mL] = 0.63"; residual is additive on
    # log10(CFU/mL). Implemented as an additive error on the Cc observation
    # (named Cc per nlmixr2lib single-output convention; values are log10
    # CFU/mL, not concentration -- see units$concentration).
    addSd <- 0.63
    label("Additive residual SD on log10(bacteria) (log10 CFU/mL)")  # Wicha 2017 Table 1, Residual variability
  })

  model({
    # ---- Per-replicate effective parameters ----
    klag        <- exp(lklag)
    krep        <- exp(lkrep)
    kdoub       <- exp(lkdoub)
    kdeath_per  <- exp(lkdeath_per)

    ec50_lzd    <- exp(lec50_lzd)
    hill_lzd    <- exp(lhill_lzd)
    kdeath_lzd  <- exp(lkdeath_lzd)

    ec50_mer_t0 <- exp(lec50_mer_t0)
    hill_mer    <- exp(lhill_mer)
    ec50_mer_eagle <- exp(lec50_mer_eagle)
    hill_mer_eagle <- exp(lhill_mer_eagle)
    b_mer       <- exp(lb_mer)
    s_mer       <- exp(ls_mer)
    kper_mer    <- exp(lkper_mer)
    kdeg_mer    <- exp(lkdeg_mer)

    ec50_van_t0 <- exp(lec50_van_t0)
    hill_van    <- exp(lhill_van)
    ec50_van_ari <- exp(lec50_van_ari)
    hill_van_ari <- exp(lhill_van_ari)
    b_van       <- exp(lb_van)
    s_van       <- exp(ls_van)
    kper_van    <- exp(lkper_van)
    kdeg_van    <- exp(lkdeg_van)

    # ---- Time-varying effective EC50 via adaption (Wicha 2017 Eqs. 5-7, 12-14) ----
    # ARon is the fractional adapted state (starts at 0); ARoff is the unadapted
    # reservoir (starts at 1). EC50(t) = (1 + b * ARon) * EC50_t=0.
    ec50_mer <- ec50_mer_t0 * (1 + b_mer * aron_mer)
    ec50_van <- ec50_van_t0 * (1 + b_van * aron_van)

    # ---- Hill drug-effect terms (pure 0-to-1 sigmoidal, Eq. 4) ----
    # Emax is held as a separate multiplier when needed (MER Eagle, VAN). For
    # LZD, MER, and the VAN-on-MER adaption-suppression term, Emax is implicit 1.
    eff_lzd      <- lzd  ^ hill_lzd     / (ec50_lzd     ^ hill_lzd     + lzd  ^ hill_lzd)
    eff_mer      <- mer  ^ hill_mer     / (ec50_mer     ^ hill_mer     + mer  ^ hill_mer)
    eff_mer_eagle <- mer ^ hill_mer_eagle / (ec50_mer_eagle ^ hill_mer_eagle + mer ^ hill_mer_eagle)
    eff_van      <- van  ^ hill_van     / (ec50_van     ^ hill_van     + van  ^ hill_van)
    eff_van_ari  <- van  ^ hill_van_ari / (ec50_van_ari ^ hill_van_ari + van  ^ hill_van_ari)

    # ---- MER effect with self-inhibitory Eagle correction (Wicha 2017 Eq. 8) ----
    # Eagle-corrected MER inhibition of doubling; equals E_MER at low MER (no
    # Eagle), and E_MER * (1 - Emax_Eagle) = 0.672 * E_MER at very high MER.
    eff_mer_eagle_corr <- eff_mer * (1 - emax_mer_eagle * eff_mer_eagle)

    # ---- Modified Bliss-Independence success-fraction for the doubling step ----
    # Wicha 2017 Eq. 8 REP -> GRO arm contains the nested success-fraction
    # [1 - E_MER * (1 - Emax_MER,Eagle * E_MER,Eagle) * (1 - E_VAN)] *
    #   (1 - Emax_VAN * E_VAN)
    # The inner (1 - E_VAN) gates the MER contribution: at high VAN it -> 0,
    # so the joint MER+VAN doubling-inhibition collapses to the VAN-only
    # (1 - Emax_VAN * E_VAN) term, matching the paper's "max joint effect was
    # limited to the effect of VAN" claim (cf. M4 vs V4M4 in Figure 2).
    success_doub <-
      (1 - eff_mer_eagle_corr * (1 - eff_van)) *
      (1 - emax_van * eff_van)

    # ---- krep(t, CFU) -- lag-phase + carrying-capacity attenuation (Eq. 11) ----
    cfu_total <- gro + repl + pers
    krep_eff  <- krep * (1 - exp(-klag * t)) * (1 - cfu_total / (10 ^ cfumax))

    # ---- ODE system ----
    # Drug solution concentrations (in vitro static dose + chemical degradation)
    d/dt(lzd) <- 0                            # LZD chemically stable in growth medium
    d/dt(mer) <- -kdeg_mer * mer              # Wicha 2017 page 5 (kdeg fixed from HPLC)
    d/dt(van) <- -kdeg_van * van

    # Adaption states (Tam 2005 formulation; Eqs. 5-6 (MER) and 12-13 (VAN)).
    # VAN inhibits MER adaption monodirectionally via eff_van_ari (Eq. text page 5).
    d/dt(aroff_mer) <- -s_mer * (1 - eff_van_ari) * mer * aroff_mer
    d/dt(aron_mer)  <-  s_mer * (1 - eff_van_ari) * mer * aroff_mer
    d/dt(aroff_van) <- -s_van * van * aroff_van
    d/dt(aron_van)  <-  s_van * van * aroff_van

    # Bacterial life cycle (Wicha 2017 Eqs. 8, 9, 10).
    # The 2 multiplier on the REP->GRO arm reflects that each successful doubling
    # generates two cells from one REP cell.
    d/dt(gro)  <- -kdeath_lzd * eff_lzd * gro -
                   krep_eff * (1 - eff_lzd) * gro +
                   2 * kdoub * success_doub * repl

    d/dt(repl) <-  krep_eff * (1 - eff_lzd) * gro -
                   kdoub * repl -
                   kper_mer * eff_mer * repl -
                   kper_van * eff_van * repl

    d/dt(pers) <-  kper_mer * eff_mer * repl +
                   kper_van * eff_van * repl -
                   kdeath_per * pers

    # Initial conditions: bacteria start in gro at 10^cfu0; adaption reservoir
    # full (ARoff = 1, ARon = 0); drug compartments at 0 unless dosed.
    gro(0)       <- 10 ^ cfu0
    aroff_mer(0) <- 1
    aroff_van(0) <- 1

    # ---- Observation ----
    # The observation Cc is the log10 of the total CFU/mL aggregated across
    # the three bacterial states (growing + replicating + persisting); a 1e-6
    # floor avoids log10(0) when all states are driven to zero. The 'Cc' name
    # is the nlmixr2lib single-output convention; the underlying quantity is
    # log10 CFU/mL, not a drug concentration -- see units$concentration and
    # the model's description for the full semantics. Additive residual on
    # log10(CFU/mL) per Wicha 2017 Table 1 (sigma = 0.63).
    Cc <- log10(gro + repl + pers + 1e-6)
    Cc ~ add(addSd)
  })
}
