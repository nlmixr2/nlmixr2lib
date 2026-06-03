Matsumoto_2005_TF_505 <- function() {
  description <- "Two-compartment first-order-absorption population PK model for the oral 5-alpha-reductase inhibitor TF-505 coupled to an indirect-response PD model for plasma dihydrotestosterone (DHT, expressed as percent of basal) in which the DHT synthesis rate kin is modulated by a 24-h circadian cosine; fit to single- and multiple-dose data from healthy adult male Japanese volunteers (Matsumoto 2005)."
  reference <- "Matsumoto Y, Fujita T, Ishida Y, Shimizu M, Kakuo H, Yamashita K, Majima M, Kumagai Y. Population Pharmacokinetic-Pharmacodynamic Modeling of TF-505 Using Extension of Indirect Response Model by Incorporating a Circadian Rhythm in Healthy Volunteers. Biol Pharm Bull. 2005;28(8):1455-1461. doi:10.1248/bpb.28.1455"
  vignette <- "Matsumoto_2005_TF_505"
  paper_specific_compartments <- c("dht")
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 36L,                                    # Matsumoto 2005 p. 1456 Methods / Table 1: 6 subjects per group x 6 dose groups
    n_studies      = 1L,                                     # Single-centre study, Clinical Investigation Center of Kitasato University East Hospital
    age_range      = "20-64 years",                          # Matsumoto 2005 Table 1 across all six dose groups
    age_median     = NULL,                                   # Not reported; group means range 22.8-52.5 years
    weight_range   = "49.6-72.3 kg",                         # Matsumoto 2005 Table 1 across all six dose groups
    weight_median  = NULL,                                   # Not reported; group means range 61.1-65.7 kg
    sex_female_pct = 0,                                      # All male
    race_ethnicity = c(Japanese = 100),                      # Matsumoto 2005 Methods 'Volunteers and Study Design': healthy male Japanese volunteers
    disease_state  = "Healthy adult male Japanese volunteers",
    dose_range     = "Single 25, 50, 75, 100 mg p.o.; multiple 12.5 or 25 mg p.o. QD x 7 days",
    regions        = "Japan",
    notes          = "Inclusion: males aged >=20 years (low single-dose groups 25 and 50 mg) or >=40 years (high single-dose 75 / 100 mg, and both multiple-dose groups); body weight within 20% of ideal. Single-dose groups dosed at 9 a.m. without breakfast; multiple-dose groups dosed at 9 a.m. after breakfast. Smoking allowed but stopped 1 h pre- through 24 h post-dose. Caffeine, alcohol and grapefruit prohibited. 564 plasma TF-505 and 264 plasma DHT measurements were used in the joint fit."
  )

  ini({
    # ============================================================
    # PK structural parameters (Matsumoto 2005 Table 4 / Results 'Pharmacokinetic and Pharmacodynamic Model').
    # Paper reports micro-rate constants ka, ke, Vc, k12, k21. We
    # reparameterise to canonical lcl, lvc, lq, lk12, lk21 with
    # CL = ke * Vc and Q derived from k12 * Vc; kel = cl / vc is
    # recovered in model(). All values are population means.
    # ============================================================
    lka  <- log(0.197)                  ; label("Absorption rate constant ka (1/h)")                              # Matsumoto 2005 Table 4: ka = 0.197 +/- 0.0361 h^-1
    lcl  <- log(0.0678 * 12.5)          ; label("Apparent clearance CL/F (L/h)")                                  # Matsumoto 2005 Table 4: ke = 0.0678 h^-1 and Vc = 12.5 L; CL = ke * Vc = 0.8475 L/h
    lvc  <- log(12.5)                   ; label("Apparent central volume of distribution Vc/F (L)")               # Matsumoto 2005 Table 4: Vc = 12.5 +/- 2.26 L
    lk12 <- log(0.0645)                 ; label("Central-to-peripheral distribution rate constant k12 (1/h)")     # Matsumoto 2005 Table 4: k12 = 0.0645 +/- 0.0218 h^-1
    lk21 <- log(0.0723)                 ; label("Peripheral-to-central distribution rate constant k21 (1/h)")    # Matsumoto 2005 Table 4: k21 = 0.0723 +/- 0.0130 h^-1

    # ============================================================
    # PD structural parameters (Matsumoto 2005 Table 4 / Eq. (1) and Eq. (2)).
    # Indirect response on DHT (% basal): d/dt(dht) = kin*(1 - Imax*Cp/(IC50+Cp)) - kout*dht;
    # circadian-modulated input: kin(t) = Rm + Ramp*cos(2*pi*(t - Tz)/24).
    # ============================================================
    lic50  <- log(1.01)                 ; label("Plasma TF-505 concentration producing 50% maximum inhibition of DHT synthesis IC50 (ug/mL)")  # Matsumoto 2005 Table 4: IC50 = 1.01 +/- 1.64 ug/mL
    lkout  <- log(0.221)                ; label("First-order DHT elimination rate kout (1/h)")                                                  # Matsumoto 2005 Table 4: kout = 0.221 +/- 0.0486 h^-1
    lrm    <- log(20.4)                 ; label("Mean DHT synthesis rate Rm (% basal / h)")                                                     # Matsumoto 2005 Table 4: Rm = 20.4 +/- 8.08 %/h
    lramp  <- log(5.06)                 ; label("Amplitude of the 24-h circadian variation in DHT synthesis rate Ramp (% basal / h)")           # Matsumoto 2005 Table 4: Ramp = 5.06 +/- 3.43 %/h
    ltacro <- log(5.01)                 ; label("Acrophase Tz: clock time (h after midnight) of maximum DHT synthesis rate")                    # Matsumoto 2005 Table 4: Tz = 5.01 +/- 0.407 h (Results paragraph 'Pharmacokinetic and Pharmacodynamic Model'; peak shifted to 14:00 by the final fit per Discussion)
    limax  <- log(0.706)                ; label("Maximum fractional inhibition of DHT synthesis by TF-505 Imax (unitless, in (0,1))")           # Matsumoto 2005 Table 4: Imax = 0.706 +/- 0.297

    # ============================================================
    # Inter-individual variability (Matsumoto 2005 Table 4 'Inter-individual variability' column;
    # all etas independent on the paper's micro-constant
    # parameterisation per Methods 'Pharmacostatistical Models':
    # 'All of PK and PD parameters were assumed lognormal
    # distribution. All of PK and PD parameters were assumed
    # linearly independence.'). Values reported are omega^2 on the
    # log scale (NONMEM Omega). For our canonical lcl-encoding,
    # etalcl carries the paper's eta_ke plus eta_Vc contribution
    # (since CL = ke * Vc on log scale), so the etalcl / etalvc
    # block is parameterised with the necessary 2x2 covariance to
    # preserve the paper's underlying independence: var(etalcl) =
    # var(eta_ke) + var(eta_Vc); cov(etalcl, etalvc) = var(eta_Vc).
    # See vignette Assumptions and deviations for the algebra.
    # ============================================================
    etalka          ~ 0.0279              # Matsumoto 2005 Table 4: omega^2_ka = 0.0279 (CV 16.70%)
    # etalcl + etalvc 2x2 block: var(etalcl) = var(eta_ke) + var(eta_Vc) = 0.0481 + 0.0180 = 0.0661 (Matsumoto 2005 Table 4 omega^2_ke = 0.0481, omega^2_Vc = 0.0180);
    # cov(etalcl, etalvc) = var(eta_Vc) = 0.0180 (shared eta_Vc from the lcl = lke + lvc reparameterisation); var(etalvc) = 0.0180.
    etalcl + etalvc ~ c(0.0661, 0.0180, 0.0180)
    etalk12         ~ 0.0106              # Matsumoto 2005 Table 4: omega^2_k12 = 0.0106 (CV 10.30%)
    # k21 IIV fixed at zero (Matsumoto 2005 Table 4 footnote a: 'Fixed at zero due to small variance estimates' applies to k21, IC50 and Ramp)

    etalkout        ~ 0.000525            # Matsumoto 2005 Table 4: omega^2_kout = 0.000525 (CV 2.29%)
    etalrm          ~ 0.0117              # Matsumoto 2005 Table 4: omega^2_Rm = 0.0117 (CV 10.82%)
    # IC50, Ramp IIVs fixed at zero in the final fit (paper Table 4 footnote a)
    etaltacro       ~ 0.00700             # Matsumoto 2005 Table 4: omega^2_Tz = 0.00700 (CV 8.37%)
    etalimax        ~ 0.00541             # Matsumoto 2005 Table 4: omega^2_Imax = 0.00541 (CV 7.36%)

    # ============================================================
    # Residual error (Matsumoto 2005 Table 4 'Intra-individual residual variability').
    # Paper reports sigma^2 directly; CV = sqrt(sigma^2) consistent
    # with text 'Intra-individual variabilities were 43.70% for PK
    # and 20.47% for PD' (sqrt(0.191) ~= 0.437, sqrt(0.0419) ~= 0.2047).
    # ============================================================
    propSd     <- sqrt(0.191)             ; label("Proportional residual error on TF-505 plasma concentration (fraction)") # Matsumoto 2005 Table 4: PK sigma^2 = 0.191 +/- 0.0226; sqrt(0.191) = 0.437 (43.7% CV per text)
    propSd_DHT <- sqrt(0.0419)            ; label("Proportional residual error on DHT (% basal)")                          # Matsumoto 2005 Table 4: PD sigma^2 = 0.0419 +/- 0.00664; sqrt(0.0419) = 0.205 (20.47% CV per text)
  })

  model({
    # ------------------------------------------------------------
    # Individual PK parameters
    # ------------------------------------------------------------
    ka  <- exp(lka  + etalka)
    cl  <- exp(lcl  + etalcl)
    vc  <- exp(lvc  + etalvc)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21)                                                  # No IIV (paper fixed at zero)

    kel <- cl / vc

    # ------------------------------------------------------------
    # Individual PD parameters
    # ------------------------------------------------------------
    ic50  <- exp(lic50)                                               # No IIV (paper fixed at zero)
    kout  <- exp(lkout + etalkout)
    rm    <- exp(lrm   + etalrm)
    ramp  <- exp(lramp)                                               # No IIV (paper fixed at zero)
    tacro <- exp(ltacro + etaltacro)
    imax  <- exp(limax  + etalimax)

    # ------------------------------------------------------------
    # Two-compartment PK with first-order absorption (Matsumoto 2005 Methods 'Pharmacokinetic-Pharmacodynamic Model', PK Model 2 in Table 2 selected)
    # ------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-               k12 * central - k21 * peripheral1

    # ------------------------------------------------------------
    # Indirect-response PD with circadian kin (Matsumoto 2005 Eq. (1) and Eq. (2))
    #   kin(t) = Rm + Ramp * cos(2*pi*(t - Tz)/24)
    #   d/dt(dht) = kin * (1 - Imax * Cp / (IC50 + Cp)) - kout * dht
    # `tacro` is the acrophase Tz from the paper (peak time of synthesis).
    # The DHT state is initialised at the no-drug circadian mean
    # rm / kout (per-individual) so the simulation starts in
    # quasi-steady-state; the circadian oscillation is established
    # by the ODE within the first few hours.
    # ------------------------------------------------------------
    Cc  <- central / vc                                               # TF-505 plasma concentration (ug/mL when dose is mg and vc is L)

    kin <- rm + ramp * cos(2 * pi * (t - tacro) / 24)                 # Time-varying DHT synthesis rate (Eq. 2)
    inh <- imax * Cc / (ic50 + Cc)                                    # Sigmoid inhibition of synthesis by drug

    d/dt(dht) <- kin * (1 - inh) - kout * dht                         # Indirect response (Eq. 1)
    dht(0)    <- rm / kout                                            # No-drug circadian mean as initial condition

    DHT <- dht                                                        # Paper-faithful uppercase output alias

    Cc  ~ prop(propSd)
    DHT ~ prop(propSd_DHT)
  })
}
