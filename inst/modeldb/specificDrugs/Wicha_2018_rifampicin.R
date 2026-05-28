Wicha_2018_rifampicin <- function() {
  description <- "Preclinical-to-clinical translational Multistate Tuberculosis Pharmacometric (MTP) framework for high-dose oral rifampicin in adults with pulmonary tuberculosis. The Svensson 2018 HIGHRIF1 plasma PK model (Erlang transit absorption + Michaelis-Menten clearance + enzyme-pool autoinduction + dose-dependent bioavailability anchored at 450 mg) is coupled via the Clewe 2015 epithelial lining fluid (ELF) effect compartment to a new post-antibiotic-effect (PAE) compartment with saturable Michaelis-Menten elimination, driving the Clewe 2016 three-state MTP model (fast-, slow-, and nonmultiplying Mycobacterium tuberculosis substates) at human-specific carrying capacity Bmax = 2.42e8/mL and fast-multiplying growth rate kG = 0.206/day. Time unit is days; all PK rates from Svensson 2018 (reported in 1/h) and the ELF kELF from Clewe 2015 are multiplied by 24 to bring to days. All structural parameters are fixed at the Wicha 2018 Table 1 typical values; only the Svensson 2018 IIV is carried (IOV is omitted because the EBA forward simulation models a single 14-day monotherapy course). The model predicts early bactericidal activity (EBA0-2 / EBA0-5 / EBA0-14) for clinical rifampicin doses 2.5-50 mg/kg without re-estimating any parameter from clinical EBA data."
  reference <- paste(
    "Wicha SG, Clewe O, Svensson RJ, Gillespie SH, Hu Y,",
    "Coates ARM, Simonsson USH. (2018).",
    "Forecasting Clinical Dose-Response From Preclinical Studies in",
    "Tuberculosis Research: Translational Predictions With Rifampicin.",
    "Clin Pharmacol Ther 104(6):1208-1218. doi:10.1002/cpt.1102.",
    "Component-model sources retained inline:",
    "PK structure and parameter values: Svensson RJ et al. (2018)",
    "Clin Pharmacol Ther 103(4):674-683 doi:10.1002/cpt.778 (HIGHRIF1).",
    "ELF effect-compartment structure and ratio: Clewe O et al. (2015)",
    "Eur J Clin Pharmacol 71(3):313-319 doi:10.1007/s00228-014-1798-3.",
    "MTP three-state disease model and growth/transfer rates:",
    "Clewe O et al. (2016) J Antimicrob Chemother 71(4):964-974",
    "doi:10.1093/jac/dkv416.",
    sep = " "
  )
  vignette <- "Wicha_2018_rifampicin"
  units <- list(
    time          = "day",
    dosing        = "mg",
    concentration = "mg/L for plasma Cc and ELF Celf and PAE Cpae; CFU/mL for the bacterial states (the log_cfu observation is on the natural-log scale)"
  )

  covariateData <- list(
    DOSE = list(
      description        = "Per-record administered rifampicin dose (mg) used as the input to the saturable dose-dependent bioavailability function bio = 1 + femax * max(DOSE - 450, 0) / (fed50 + max(DOSE - 450, 0)).",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reference dose 450 mg (the standard adult rifampicin dose) where bio = 1",
        "by construction. The Svensson 2018 HIGHRIF1 calibration spans 600-2100 mg",
        "(10-35 mg/kg in adults); Wicha 2018 extrapolates the EBA simulation up to",
        "50 mg/kg (~3000 mg for a 60 kg adult) and down to 2.5 mg/kg (~150 mg).",
        "Below 450 mg the saturable function is clamped to bio = 1 (max(DOSE-450, 0)",
        "guard); the Svensson 2018 publication notes the function describes the",
        "increase ABOVE the 450 mg reference dose, not behaviour below it.",
        "Wicha 2018 Discussion notes the 50 mg/kg prediction is a slight",
        "extrapolation since Svensson 2018 PK was formally only studied up to",
        "40 mg/kg."
      ),
      source_name        = "DOSE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = 0L,
    age_range      = "adults (TB patients; the Svensson 2018 HIGHRIF1 source cohort enrolled adult pulmonary TB patients)",
    weight_range   = "log-normal with geometric mean 60 kg and geometric SD 10% (Wicha 2018 Methods 'Pharmacokinetics of rifampicin in the different target systems')",
    sex_female_pct = NA_real_,
    race_ethnicity = "not specified; the Wicha 2018 simulation is intended for global TB populations",
    disease_state  = "Adults with active drug-susceptible pulmonary Mycobacterium tuberculosis infection (simulated phase IIa early bactericidal activity cohort)",
    dose_range     = "Oral rifampicin 2.5 to 50 mg/kg once daily for 14 days (Wicha 2018 Methods 'Translational prediction from in vitro to the target systems'). The reported predictions emphasize the 10, 25, 35, and 50 mg/kg dose levels.",
    regions        = "global TB cohorts; the EBA validation in Wicha 2018 Figure 4a pools observed EBA from contemporary clinical trials (Boeree 2015, Jindani 1980, Sirgel 2005, Diacon 2007, Chan 1992, Rustomjee 2008).",
    notes          = paste(
      "Wicha 2018 is a SIMULATION study built on three previously published model",
      "components linked in this work: Svensson 2018 HIGHRIF1 popPK (saturable",
      "Michaelis-Menten clearance + enzyme-pool autoinduction + dose-dependent F),",
      "Clewe 2015 ELF effect compartment, and Clewe 2016 in-vitro MTP three-state",
      "bacterial disease model. The Wicha 2018 paper's novel contributions are the",
      "PAE compartment with saturable Michaelis-Menten elimination and the MIC",
      "scaling for translational potency adjustment. No new subjects were enrolled.",
      "The simulation samples WT from log-normal (geometric mean 60 kg, GSD 10%)",
      "and HT from log-normal (geometric mean 1.75 m, GSD 7.5%) to drive the",
      "Svensson 2018 FFM-allometric scaling on Vmax (exponent 0.75) and Vd",
      "(exponent 1.0); for this model file the parameter values are stored at",
      "their 70 kg FFM reference and the allometric scaling is left to the",
      "vignette covariate-derivation chunk so the model file stays self-contained.",
      "IOV reported in Wicha 2018 Table 1 (IOV km 18.9%, IOV ka 31.4%, IOV MTT",
      "56.4%, IOV F 15.7%) is NOT carried because the EBA simulation models a",
      "single 14-day course; see vignette Assumptions and deviations."
    )
  )

  ini({
    # =========================================================================
    # Svensson 2018 plasma PK (Wicha 2018 Table 1, 'Clinical phase IIa' block).
    # Time unit is DAYS throughout. The Svensson 2018 source publication reports
    # rate constants in 1/h and the time constant MTT in h; values below
    # multiply by 24 for rate constants (or divide by 24 for time constants) to
    # bring to days. Concentrations stay in mg/L.
    #
    # Vmax and Vd are reported standardised to 70 kg FFM in Wicha 2018 Table 1
    # ('70 kg^-1'). The model is stored at the 70 kg reference; the vignette
    # applies the (FFM / 70)^0.75 / (FFM / 70)^1.0 allometric correction on top
    # of these values for the simulated cohort.
    # =========================================================================
    lvmax   <- log(525 * 24)
    label("Michaelis-Menten Vmax of clearance at FFM = 70 kg (mg/day)")
    # Wicha 2018 Table 1 Vmax = 525 mg/h/70kg (Svensson 2018 THETA1). x 24 h/day.

    lkm     <- log(35.3)
    label("Michaelis-Menten Km of clearance (mg/L)")
    # Wicha 2018 Table 1 km = 35.3 mg/L (Svensson 2018 THETA2).

    lvc     <- log(87.2)
    label("Apparent central volume of distribution at FFM = 70 kg (L)")
    # Wicha 2018 Table 1 Vd = 87.2 L/70kg (Svensson 2018 THETA3).

    lka     <- log(1.77 * 24)
    label("First-order absorption rate constant ka (1/day)")
    # Wicha 2018 Table 1 ka = 1.77/h (Svensson 2018 THETA4). x 24 h/day.

    emax    <- 1.16
    label("Autoinduction Emax (unitless multiplier on enzyme synthesis)")
    # Wicha 2018 Table 1 Emax = 1.16 (Svensson 2018 THETA5).

    ec50    <- 0.0699
    label("Autoinduction EC50 on rifampicin Cc (mg/L)")
    # Wicha 2018 Table 1 EC50 = 0.0699 mg/L (Svensson 2018 THETA6).

    kenz    <- 0.00603 * 24
    label("First-order enzyme turnover rate constant (1/day)")
    # Wicha 2018 Table 1 kENZ = 0.00603/h (Svensson 2018 THETA7). x 24 h/day.

    lmtt    <- log(0.513 / 24)
    label("Mean transit time through the absorption chain MTT (day)")
    # Wicha 2018 Table 1 MTT = 0.513 h (Svensson 2018 THETA8). / 24 h/day.

    lnn     <- log(23.8)
    label("Erlang transit-chain shape parameter NN (unitless)")
    # Wicha 2018 Table 1 NN = 23.8 (Svensson 2018 THETA9). Non-integer;
    # rxode2's transit(n, mtt, bio) closed form accepts non-integer n via the
    # gamma-density evaluation.

    femax   <- 0.504
    label("Dose-dependent bioavailability Emax (fraction)")
    # Wicha 2018 Table 1 Fmax = 0.504 (Svensson 2018 THETA10). Saturable Emax of
    # the F-vs-DOSE response anchored at the 450 mg reference dose.

    fed50   <- 67
    label("Dose-dependent bioavailability EC50 offset above 450 mg (mg)")
    # Wicha 2018 Table 1 ED50 = 67.0 mg (Svensson 2018 THETA11).

    # =========================================================================
    # Clewe 2015 epithelial lining fluid (ELF) effect compartment.
    # Wicha 2018 Table 1 'Clinical phase IIa' block carries the Clewe 2015
    # fixed values. kELF reported in 1/h, multiplied by 24 to bring to days.
    # The plasma free fraction fu is reported but is NOT used to scale the ELF
    # input: the Clewe 2015 effect compartment uses TOTAL plasma Cc as the
    # input with R_ELF being the total-plasma ratio (the unbound-plasma ratio
    # R_ELF/unbound = 1.28 is recovered post hoc by dividing by fu = 0.2).
    # =========================================================================
    lkelf   <- fixed(log(41.58 * 24))
    label("Plasma-to-ELF effect-compartment rate constant kELF (1/day)")
    # Wicha 2018 Table 1 kELF = 41.58/h FIX (Clewe 2015). x 24 h/day. Equivalent
    # to a 1-minute equilibration half-life; ELF concentration tracks plasma at
    # the steady-state ratio R_ELF/plasma.

    lrelf   <- fixed(log(0.26))
    label("ELF / total-plasma concentration ratio at pseudo steady-state (unitless)")
    # Wicha 2018 Table 1 R_ELF/plasma = 0.26 FIX (Clewe 2015 Table 1).

    lfu_p   <- fixed(log(0.2))
    label("Rifampicin plasma free fraction fu (unitless)")
    # Wicha 2018 Table 1 fu = 0.2 FIX (Clewe 2015 ref 23). Carried for vignette
    # display; the ELF transit equation uses total plasma Cc (see notes above).

    # =========================================================================
    # PAE (post-antibiotic-effect) compartment. Developed in Wicha 2018 from
    # the Gumbo 2007 PAE dataset. The PAE concentration receives first-order
    # input from the ELF (target-site) concentration at rate ke,in and is
    # eliminated by a saturable Michaelis-Menten kinetics parameterised by
    # ke,out,max (limiting first-order rate at low Cpae) and ke,out,50
    # (Cpae at half the saturation):
    #   dCpae/dt = ke,in * Celf - ke,out,max * Cpae / (1 + Cpae/ke,out,50)
    # which is algebraically identical to a Michaelis-Menten removal with
    # Vmax_pae = ke,out,max * ke,out,50 (in mg/L/day) and Km_pae = ke,out,50.
    # All three rates are reported in /day in Wicha 2018 Table 1.
    # =========================================================================
    lke_in     <- fixed(log(150))
    label("Rate of rifampicin input from ELF into the PAE compartment ke,in (1/day)")
    # Wicha 2018 Table 1 ke,in = 150/day FIX (estimated in this paper from
    # Gumbo 2007 PAE data, AIC vs first-order alternative).

    lke_out_max <- fixed(log(1.091))
    label("Limiting first-order PAE elimination rate at low Cpae, ke,out,max (1/day)")
    # Wicha 2018 Table 1 ke,out,max = 1.091/day FIX (estimated in this paper).

    lke_out_50 <- fixed(log(0.662))
    label("PAE concentration at half saturation, ke,out,50 (mg/L)")
    # Wicha 2018 Table 1 ke,out,50 = 0.662 mg/L FIX (estimated in this paper).

    # =========================================================================
    # MTP (Multistate Tuberculosis Pharmacometric) model parameters (Clewe
    # 2016). Wicha 2018 Table 1 'Pharmacodynamics (MTP model)' block; all
    # rates in /day, all fixed. Carrying capacity Bmax and fast-multiplying
    # growth rate kG are human-specific (Wicha 2018 Table 1 'Value' column
    # 'human' line).
    # =========================================================================
    lkfn     <- fixed(log(0.897e-6))
    label("Transfer rate from fast- to nonmultiplying state kFN (1/day)")
    # Wicha 2018 Table 1 kFN = 0.897e-6/day FIX (Clewe 2016 ref 4).

    lksn     <- fixed(log(0.186))
    label("Transfer rate from slow- to nonmultiplying state kSN (1/day)")
    # Wicha 2018 Table 1 kSN = 0.186/day FIX (Clewe 2016 ref 4).

    lksf     <- fixed(log(0.0145))
    label("Transfer rate from slow- to fast-multiplying state kSF (1/day)")
    # Wicha 2018 Table 1 kSF = 0.0145/day FIX (Clewe 2016 ref 4).

    lkns     <- fixed(log(0.00123))
    label("Transfer rate from non- to slow-multiplying state kNS (1/day)")
    # Wicha 2018 Table 1 kNS = 0.123e-2/day FIX (Clewe 2016 ref 4).

    lkfslin  <- fixed(log(0.00166))
    label("Time-dependent slope of kFS(t) = kFS_lin * t (1/day^2)")
    # Wicha 2018 Table 1 kFS,lin = 0.166e-2/day^2 FIX (Clewe 2016 ref 4).

    lkg      <- fixed(log(0.206))
    label("Fast-multiplying bacterial growth rate kG (1/day; human value)")
    # Wicha 2018 Table 1 kG = 0.206/day for human prediction FIX (the same row
    # also lists kG = 0.150 hollow-fiber estimated and kG = 0.206 mouse;
    # human and mouse share the kG = 0.206 value per the Clewe 2016 reference).

    lf0      <- fixed(log(4.1))
    label("Initial bacterial number of fast-multiplying state F0 (1/mL)")
    # Wicha 2018 Table 1 F0 = 4.1/mL FIX (Clewe 2016 ref 4); the hollow-fiber
    # column lists 4.1 * 50 to account for the larger inoculum, not used here.

    lrbase      <- fixed(log(9770))
    label("Initial bacterial number of slow-multiplying state S0 (1/mL)")
    # Wicha 2018 Table 1 S0 = 9770/mL FIX (Clewe 2016 ref 4); hollow-fiber
    # column lists 9770 * 50 inoculum, not used here.

    lbmax    <- fixed(log(2.42e8))
    label("System carrying capacity Bmax (1/mL; human-specific value)")
    # Wicha 2018 Table 1 Bmax = 2.42e8/mL for human prediction FIX (Clewe 2016
    # ref 4). The same row lists 2.02e9 (hollow-fiber) and 4e6 (mouse).

    lfg_k    <- fixed(log(0.017))
    label("Linear inhibition of fast-multiplying bacterial growth FG,k (L/mg)")
    # Wicha 2018 Table 1 FG,k = 0.017 L/mg FIX (Clewe 2016 ref 4).

    lfd_emax <- fixed(log(2.15))
    label("Maximal fast-multiplying bacterial death rate FD,Emax (1/day)")
    # Wicha 2018 Table 1 FD,Emax = 2.15/day FIX (Clewe 2016 ref 4).

    lfd_ec50 <- fixed(log(0.52))
    label("Rifampicin concentration at half FD,Emax (mg/L)")
    # Wicha 2018 Table 1 FD,EC50 = 0.52 mg/L FIX (Clewe 2016 ref 4).

    lsd_emax <- fixed(log(1.56))
    label("Maximal slow-multiplying bacterial death rate SD,Emax (1/day)")
    # Wicha 2018 Table 1 SD,Emax = 1.56/day FIX (Clewe 2016 ref 4).

    lsd_ec50 <- fixed(log(13.4))
    label("Rifampicin concentration at half SD,Emax (mg/L)")
    # Wicha 2018 Table 1 SD,EC50 = 13.4 mg/L FIX (Clewe 2016 ref 4).

    lnd_k    <- fixed(log(0.24))
    label("Linear nonmultiplying death rate ND,k (L/(mg*day))")
    # Wicha 2018 Table 1 ND,k = 0.24 L/(mg*day) FIX (Clewe 2016 ref 4).
    # Wicha 2018 Table 1 prints the unit string as 'L * mg * days-1'; the
    # canonical Clewe 2016 form is L/(mg*day) for a second-order death rate.

    lfdepot  <- fixed(log(1))
    label("Bioavailability anchor (F = 1 at the 450 mg reference dose)")
    # CL/F and V/F are apparent F-relative in Svensson 2018; the F(DOSE)
    # multiplier is folded into the transit() bio argument inside model().

    # =========================================================================
    # IIV (Wicha 2018 Table 1 IIV rows). Reported as CV%; omega^2 = log(1 +
    # CV^2). The Vmax-km block correlation is 38.9%, reproduced via cov =
    # corr * sqrt(omega^2(vmax) * omega^2(km)). IOV (also in Table 1) is
    # NOT carried -- see file-header description.
    # =========================================================================
    etalvmax + etalkm ~ c(0.08618, 0.03964, 0.12054)
    # log(1 + 0.30^2) = 0.08618 (Wicha 2018 Table 1 IIV Vmax = 30.0% CV);
    # log(1 + 0.358^2) = 0.12054 (Wicha 2018 Table 1 IIV km = 35.8% CV);
    # cov = 0.389 * sqrt(0.08618 * 0.12054) = 0.03964 (Wicha 2018 Table 1
    # 'Correlation Vmax-km [%]' = 38.9%; sign convention positive).

    etalvc  ~ 0.006159
    # log(1 + 0.0786^2) = 0.006159 (Wicha 2018 Table 1 IIV Vd = 7.86% CV).

    etalka  ~ 0.10822
    # log(1 + 0.338^2) = 0.10822 (Wicha 2018 Table 1 IIV ka = 33.8% CV).

    etalmtt ~ 0.13624
    # log(1 + 0.382^2) = 0.13624 (Wicha 2018 Table 1 IIV MTT = 38.2% CV).

    etalnn  ~ 0.47427
    # log(1 + 0.779^2) = 0.47427 (Wicha 2018 Table 1 IIV NN = 77.9% CV).

    # =========================================================================
    # Residual error. Wicha 2018 is a forward simulation, not a popPK / popPD
    # fit, and does not report observation residual variance directly. The PK
    # residual is carried from Svensson 2018 (proportional 23.56% on log Cc,
    # the .lst-derived value used in Svensson_2018_rifampicin.R in this
    # package's ddmore tree). The bacterial residual on log(F+S)/mL is
    # carried from Svensson 2016 (combined 110% + 23.1% additive components
    # on log-CFU/mL -> 1.124 combined SD); see Svensson_2016_rifampicin.R
    # vignette for the L2-nested rationale. These residual SDs make the
    # model usable for VPC-style stochastic simulation in nlmixr2; they are
    # not part of the Wicha 2018 EBA prediction itself.
    # =========================================================================
    propSd            <- 0.2356
    label("Proportional residual error on log-Cc plasma rifampicin (carried from Svensson 2018)")
    # Svensson 2018 .lst SIGMA(1,1) FINAL = 5.55e-2 = 0.2356^2.

    addSd_log_cfu <- 1.124
    label("Additive residual SD on natural-log(F+S)/mL sputum CFU (carried from Svensson 2016 combined estimate)")
    # Svensson 2016 Table 2 sigma_e = 110% CV + sigma_e,repl = 23.1% CV,
    # combined as sqrt(1.10^2 + 0.231^2) = 1.124 per Svensson_2016_rifampicin.R
    # vignette Errata.
  })

  model({
    # -----------------------------------------------------------------------
    # 1. Individual PK parameters with IIVs. All in day-units after the
    #    /h-to-/day conversions in ini().
    # -----------------------------------------------------------------------
    vmax <- exp(lvmax + etalvmax)
    km   <- exp(lkm   + etalkm)
    vc   <- exp(lvc   + etalvc)
    ka   <- exp(lka   + etalka)
    mtt  <- exp(lmtt  + etalmtt)
    nn   <- exp(lnn   + etalnn)

    # Fixed structural and PD parameters back-transformed from log-scale.
    kelf       <- exp(lkelf)
    relf       <- exp(lrelf)
    fu_p       <- exp(lfu_p)
    ke_in      <- exp(lke_in)
    ke_out_max <- exp(lke_out_max)
    ke_out_50  <- exp(lke_out_50)
    kfn        <- exp(lkfn)
    ksn        <- exp(lksn)
    ksf        <- exp(lksf)
    kns        <- exp(lkns)
    kfslin     <- exp(lkfslin)
    kg         <- exp(lkg)
    f0         <- exp(lf0)
    rbase         <- exp(lrbase)
    bmax       <- exp(lbmax)
    fg_k       <- exp(lfg_k)
    fd_emax    <- exp(lfd_emax)
    fd_ec50    <- exp(lfd_ec50)
    sd_emax    <- exp(lsd_emax)
    sd_ec50    <- exp(lsd_ec50)
    nd_k       <- exp(lnd_k)

    # -----------------------------------------------------------------------
    # 2. Saturable dose-dependent bioavailability (Svensson 2018 Table 1
    #    Fmax / ED50). The Wicha 2018 EBA simulation spans 2.5-50 mg/kg
    #    (~150-3000 mg for a 60 kg adult); for DOSE <= 450 mg the saturable
    #    increase is clamped to zero so the function stays at bio = 1 (the
    #    Svensson 2018 formula was calibrated only above 450 mg).
    # -----------------------------------------------------------------------
    dose_excess <- 0
    if (DOSE > 450) dose_excess <- DOSE - 450
    bio <- (1 + femax * dose_excess / (fed50 + dose_excess)) * exp(lfdepot)

    # -----------------------------------------------------------------------
    # 3. Time-dependent fast-to-slow transfer kFS(t) = kFS_lin * t (Clewe
    #    2016 / Wicha 2018). The variable `t` is the elapsed simulation
    #    time in days. Wicha 2018 Methods 'Translational factors' states
    #    that a 150-day preincubation period is assumed for the human EBA
    #    simulation; user event tables should set the rifampicin first-dose
    #    time to t = 150 day so the bacterial states evolve under natural
    #    growth before treatment begins.
    # -----------------------------------------------------------------------
    kfs <- kfslin * t

    # -----------------------------------------------------------------------
    # 4. Plasma concentration Cc (mg/L) and autoinduction driver.
    # -----------------------------------------------------------------------
    Cc  <- central / vc
    eff <- emax * Cc / (ec50 + Cc)

    # -----------------------------------------------------------------------
    # 5. PK ODE system (Svensson 2018 HIGHRIF1). The rxode2 closed-form
    #    transit(n, mtt, bio) returns the gamma-density input rate from the
    #    most recent dose; f(depot) <- 0 below disables normal depot
    #    accumulation so transit() is the only depot input.
    # -----------------------------------------------------------------------
    d/dt(depot)   <- transit(nn, mtt, bio) - ka * depot
    d/dt(central) <- ka * depot - vmax * Cc / (km + Cc) * enzyme
    d/dt(enzyme)  <- kenz * (1 + eff) - kenz * enzyme

    central(0) <- 0.0001
    enzyme(0)  <- 1
    # f(depot) <- 0 disables normal depot accumulation; the transit()
    # closed form is the only depot input. The lfdepot anchor is folded
    # into bio above (multiplicative; lfdepot = log(1) so the effective
    # factor is 1).
    f(depot)   <- 0

    # -----------------------------------------------------------------------
    # 6. Clewe 2015 ELF effect compartment. dCelf/dt = kELF*(R_ELF*Cc - Celf).
    #    At steady state Celf = R_ELF * Cc = 0.26 * Cc; the kELF = 41.58/h
    #    (998/day) input rate gives a ~1-minute equilibration half-life so
    #    Celf tracks plasma rapidly.
    # -----------------------------------------------------------------------
    d/dt(effect1) <- kelf * (relf * Cc - effect1)
    Celf <- effect1

    # -----------------------------------------------------------------------
    # 7. PAE compartment (Wicha 2018 Methods 'PAE model'). The paper
    #    describes a "rapid equilibrium of rifampicin in the effect
    #    compartment with rifampicin concentrations at the target site"
    #    followed by a saturable Michaelis-Menten elimination. Encoded as
    #    an effect compartment with bidirectional first-order exchange at
    #    rate ke_in (fast equilibration; ke_in = 150/day, t_eq half-life
    #    ~6.6 min, so Cpae tracks Celf closely while drug is present) plus
    #    an additional saturable Michaelis-Menten removal. The
    #    Michaelis-Menten form ke_out_max * Cpae / (1 + Cpae/ke_out_50) is
    #    algebraically identical to standard Vmax * Cpae / (Km + Cpae)
    #    with Vmax = ke_out_max * ke_out_50 (in mg/L/day) and Km = ke_out_50
    #    (in mg/L). The multi-day post-antibiotic effect arises primarily
    #    from MTP bacterial-state dynamics (slow N -> S -> F transitions);
    #    the saturable Cpae elimination adds a small additional persistence
    #    contribution. The paper's exact PAE-compartment ODE is not
    #    written out; this interpretation is the most physically sensible
    #    one consistent with the "rapid equilibrium ... followed by ...
    #    Michaelis-Menten" textual description; see vignette Assumptions
    #    and deviations for the alternative one-way-input reading.
    # -----------------------------------------------------------------------
    d/dt(effect2) <- ke_in * (Celf - effect2) -
                     ke_out_max * effect2 / (1 + effect2 / ke_out_50)
    Cpae <- effect2

    # -----------------------------------------------------------------------
    # 8. Drug effects on the MTP bacterial states, driven by Cpae (the
    #    Wicha 2018 translational MTP model uses the PAE concentration as
    #    the universal killing/inhibition driver after the ELF + PAE chain).
    #    MIC scaling = 1 (typical wild-type isolate). For per-subject MIC
    #    variability simulation users can override fd_ec50, sd_ec50, fg_k,
    #    nd_k in the rxode2 model parameter table.
    # -----------------------------------------------------------------------
    fg_inhib <- 1 - fg_k * Cpae
    if (fg_inhib < 0) fg_inhib <- 0
    fd_drug  <- fd_emax * Cpae / (fd_ec50 + Cpae)
    sd_drug  <- sd_emax * Cpae / (sd_ec50 + Cpae)
    nd_drug  <- nd_k * Cpae

    # -----------------------------------------------------------------------
    # 9. Total live bacteria. The 1e-6 floor avoids log(0) when the killing
    #    integrates to ~0 at high doses.
    # -----------------------------------------------------------------------
    total_bact <- fast + slow + nonm + 1e-6

    # -----------------------------------------------------------------------
    # 10. MTP three-state ODE system (Clewe 2016 / Wicha 2018 Eqs. 1-3).
    #     The Gompertz growth term kg * fast * log(bmax / total_bact)
    #     saturates fast multiplication at the carrying capacity bmax. Drug
    #     effects: FG_k inhibits fast growth multiplicatively, FD_drug adds
    #     a death rate on fast bacteria, SD_drug a death rate on slow
    #     bacteria, ND_drug a death rate on nonmultiplying bacteria.
    # -----------------------------------------------------------------------
    d/dt(fast) <- kg * fast * log(bmax / total_bact) * fg_inhib -
                  kfs * fast - kfn * fast + ksf * slow -
                  fd_drug * fast
    d/dt(slow) <- kfs * fast - ksf * slow - ksn * slow + kns * nonm -
                  sd_drug * slow
    d/dt(nonm) <- kfn * fast + ksn * slow - kns * nonm -
                  nd_drug * nonm

    fast(0) <- f0
    slow(0) <- rbase
    nonm(0) <- 0

    # -----------------------------------------------------------------------
    # 11. Observations.
    #     - Cc: plasma rifampicin concentration (mg/L), proportional residual
    #       error carried from Svensson 2018.
    #     - log_cfu: natural log of culturable bacteria (F + S, with
    #       the N state excluded because it is nonculturable per Clewe
    #       2016 / Svensson 2016 Methods), additive residual on the log
    #       scale carried from Svensson 2016 combined estimate.
    # -----------------------------------------------------------------------
    log_cfu <- log(fast + slow + 1e-6)
    Cc           ~ prop(propSd)
    log_cfu ~ add(addSd_log_cfu)
  })
}
