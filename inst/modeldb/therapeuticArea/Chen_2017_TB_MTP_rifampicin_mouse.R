Chen_2017_TB_MTP_rifampicin_mouse <- function() {
  description <- "Preclinical (C57BL/6 mouse). Multistate Tuberculosis Pharmacometric (MTP) model for rifampicin in an acute M. tuberculosis H37Rv lung-infection mouse model: a one-compartment first-order-absorption rifampicin PK model with auto-inducible clearance (fixed from a prior mouse popPK) drives an Emax inhibition of fast-multiplying bacterial growth, a sigmoid Emax stimulation of fast-multiplying bacterial death and a linear stimulation of slow-multiplying bacterial death in a three-state (fast-, slow-, non-multiplying) bacterial system; the CFU observation is the natural log of fast- plus slow-multiplying bacteria in the lungs."
  reference <- paste(
    "Chen C, Ortega F, Rullas J, Alameda L, Angulo-Barturen I, Ferrer S,",
    "Simonsson USH. (2017). The multistate tuberculosis pharmacometric model:",
    "a semi-mechanistic pharmacokinetic-pharmacodynamic model for studying",
    "drug effects in an acute tuberculosis mouse model.",
    "J Pharmacokinet Pharmacodyn 44(2):133-141. doi:10.1007/s10928-017-9508-2.",
    "Final NONMEM control stream: Electronic Supplementary Material S2.",
    "Fixed rifampicin PK from Chen C, Ortega F, Alameda L, Ferrer S,",
    "Simonsson USH. Eur J Pharm Sci 93:319-333 (2016)",
    "doi:10.1016/j.ejps.2016.07.017; fixed MTP transfer rates from Clewe O,",
    "Aulin L, Hu Y, Coates AR, Simonsson USH. J Antimicrob Chemother",
    "71(4):964-974 (2016) doi:10.1093/jac/dkv416."
  )
  vignette <- "Chen_2017_TB_MTP_rifampicin_mouse"
  units <- list(
    time = "h",
    dosing = "mg/kg",
    concentration = "ug/mL (rifampicin blood concentration Cc); log_cfu is the natural log of CFU/lungs"
  )

  covariateData <- list(
    TRT_RIFAMPICIN = list(
      description = "Rifampicin treatment-group indicator: 1 = mouse in a rifampicin-treated group, 0 = untreated (natural-growth) group",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Selects the initial fast-multiplying bacterial number F0: 331000",
        "CFU/lungs in treated mice and 8340 CFU/lungs in untreated mice",
        "(Chen 2017 Table 1). Time-fixed per mouse and set to 1 from the",
        "infection (t = 0), i.e. also before the first dose at 24 h. The",
        "rifampicin drug effects themselves are concentration-driven and need",
        "no indicator."
      ),
      source_name = "DRUG (S2 control stream: IF(DRUG.EQ.0) TVAF = THETA(8))"
    ),
    DOSE_LOW = list(
      description = "Lowest-dose-cohort indicator: 1 = the 1.00 mg/kg/day rifampicin cohort, 0 = any higher dose level (2.83 to 98.0 mg/kg/day) or untreated",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Selects the apparent volume of distribution: 2278.74 mL/kg at the",
        "lowest dose and 1249.14 mL/kg at all higher doses (S2 control stream",
        "IF(DOSE.LE.1020) TVV = THETA(4), DOSE in ug/kg). Only the 1.00 mg/kg",
        "cohort of this study falls at or below the 1.02 mg/kg threshold."
      ),
      source_name = "DOSE (S2 control stream, ug/kg; DOSE_LOW = DOSE <= 1020)"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "rifampicin", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "rifampicin", units = "mg/kg", specimen = "whole blood", verified = TRUE),
    fbugs = list(
      analyte = "fast-multiplying M. tuberculosis",
      units = "CFU/lungs",
      specimen = "not applicable",
      verified = TRUE
    ),
    sbugs = list(
      analyte = "slow-multiplying M. tuberculosis",
      units = "CFU/lungs",
      specimen = "not applicable",
      verified = TRUE
    ),
    nbugs = list(
      analyte = "non-multiplying M. tuberculosis",
      units = "CFU/lungs",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "mouse (C57BL/6)",
    n_subjects = 60L,
    n_studies = 1L,
    disease_state = paste(
      "Acute tuberculosis: intratracheal infection with M. tuberculosis H37Rv",
      "(10^5 CFU per mouse) on Day 0."
    ),
    dose_range = paste(
      "Oral rifampicin once daily from Day 1 (24 h after infection):",
      "1.00, 2.83, 8.88, 26.4 or 98.0 mg/kg/day for 8 days (5 mice per dose,",
      "sacrificed Day 9); 30 mg/kg/day for up to 8 days (5 mice sacrificed on",
      "each of Days 2, 3, 4 and 9); 15 untreated mice sacrificed on Days 1, 9",
      "and 18 (5 per day)."
    ),
    regions = "GlaxoSmithKline Tres Cantos, Madrid, Spain",
    notes = paste(
      "Chen 2017 Methods 'Rifampicin pharmacodynamic mouse study'. Endpoint is",
      "CFU in the homogenised lungs; LLOQ 10 CFU/lungs, 21.6% of observations",
      "below the LLOQ (handled by the NONMEM M3 method). No PK was sampled in",
      "this study: rifampicin blood concentrations come from the earlier mouse",
      "popPK model (population PK parameter approach) with every PK value fixed.",
      "No interindividual variability was estimated ($OMEGA 0 FIX in S2)."
    )
  )

  ini({
    # ---- Rifampicin PK: all FIXED from the upstream mouse popPK (S2 $THETA 1-6) ----
    lcl <- fixed(log(79.321))
    label("Apparent clearance CL/F before the auto-induction step (mL/h/kg)")
    # S2 $THETA(1) = 79.321 FIX 'CL at Days <= 2'; Methods text 79.3 mL/h/kg.

    lcl_late <- fixed(log(131.911))
    label("Apparent clearance CL/F after the auto-induction step (mL/h/kg)")
    # S2 $THETA(2) = 131.911 FIX 'CL at Days >= 3'; Methods text 132 mL/h/kg.

    ltclchange <- fixed(log(48))
    label("Time of the clearance step, from infection (h)")
    # S2 $PK: IF(TIME.GE.48) TVCL = THETA(2). TIME in this dataset is time since
    # infection (F inoculum at TIME 0, first dose at 24 h), so the step falls at
    # the second dose; the $THETA comments 'Days <= 2' / 'Days >= 3' number the
    # upstream PK study's dosing days.

    lvc <- fixed(log(1249.14))
    label("Apparent volume V/F at doses above 1.02 mg/kg (mL/kg)")
    # S2 $THETA(3) = 1249.140 FIX 'V at Dose > 1020 ug/kg'; Methods text 1250 mL/kg.

    lvc_lowdose <- fixed(log(2278.74))
    label("Apparent volume V/F at the lowest dose, <= 1.02 mg/kg (mL/kg)")
    # S2 $THETA(4) = 2278.740 FIX 'V at Dose <= 1020 ug/kg'; Methods text 2280 mL/kg.

    lka_excess <- fixed(log(0.994))
    label("Excess of the absorption rate constant over kel, ka = kel + ka_excess (1/h)")
    # S2 $PK: TVKA = KE + THETA(5); $THETA(5) = 0.994 FIX.

    lfdepot <- fixed(log(0.656))
    label("Oral bioavailability F (fraction)")
    # S2 $THETA(6) = 0.656 FIX; Methods text 'bioavailability was estimated as 65.6%'.

    # ---- MTP natural growth (Table 1) ----
    lf0 <- log(331000)
    label("Initial fast-multiplying bacteria F0 in rifampicin-treated mice (CFU/lungs)")
    # Table 1: F0 treated = 331000, RSE 54.1%.

    lf0_untreated <- log(8340)
    label("Initial fast-multiplying bacteria F0 in untreated mice (CFU/lungs)")
    # Table 1: F0 untreated = 8340, RSE 8.7%.

    lkg <- log(0.0333)
    label("Exponential growth rate of fast-multiplying bacteria kG (1/h)")
    # Table 1: kG = 0.0333 1/h, RSE 4.4%.

    lkfslin <- log(2.96e-5)
    label("Slope of the time-dependent F->S transfer rate, kFS = kFSlin * t (1/h^2)")
    # Table 1: kFSlin = 2.96e-5 1/h^2, RSE 29.8%; Eq. 1.

    lksf <- fixed(log(6.042e-4))
    label("First-order S->F transfer rate kSF (1/h)")
    # Table 1: 6.04e-4 FIX; S2 $THETA(12) = 0.006042 FIX with KSF = THETA(12)/10 (in vitro 0.0145 1/day / 24).

    lkfn <- fixed(log(3.7375e-8))
    label("First-order F->N transfer rate kFN (1/h)")
    # Table 1: 3.74e-8 FIX; S2 $THETA(11) = 0.037375 FIX with KFN = THETA(11)/1e6 (in vitro 0.897e-6 1/day / 24).

    lksn <- fixed(log(7.75e-3))
    label("First-order S->N transfer rate kSN (1/h)")
    # Table 1: 7.75e-3 FIX; S2 $THETA(13) = 0.007750 FIX (in vitro 0.186 1/day / 24).

    lkns <- fixed(log(5.125e-5))
    label("First-order N->S transfer rate kNS (1/h)")
    # Table 1: 5.12e-5 FIX; S2 $THETA(14) = 0.005125 FIX with KNS = THETA(14)/100 (in vitro 0.123e-2 1/day / 24).

    # ---- Rifampicin exposure-response (Table 1; Eqs. 3-5) ----
    lemax_fg_rif <- log(0.716)
    label("FGImax: maximum fractional inhibition of fast-multiplying growth (fraction)")
    # Table 1: FGImax = 0.716, RSE 16.2%.

    lec50_fg_rif <- log(0.0397)
    label("FGIC50: rifampicin concentration at half of FGImax (ug/mL)")
    # Table 1: FGIC50 = 0.0397 ug/mL, RSE 50.2%.

    gamma_fd_r <- 2.19
    label("Sigmoidicity of the fast-multiplying death sigmoid Emax (unitless)")
    # Table 1: gamma = 2.19, RSE 7.2%.

    lemax_fd_rif <- log(352)
    label("FDEmax: maximum fractional stimulation of fast-multiplying death (1/h)")
    # Table 1: FDEmax = 352, RSE 12.8%. (S2 $THETA(18) initial 67.0823;
    # Table 1 reports the final estimate 352, which is authoritative.)

    lec50_fd_rif <- log(212)
    label("FDEC50: rifampicin concentration at half of FDEmax (ug/mL)")
    # Table 1: FDEC50 = 212 ug/mL, RSE 30.4%. (S2 $THETA(19) initial 98.5015;
    # Table 1 reports the final estimate 212, which is authoritative.)

    lksd_rif <- log(4.91e-3)
    label("SDk: second-order slow-multiplying death rate slope (mL/(h*ug))")
    # Table 1: SDk = 4.91e-3, RSE 12.9%; Eq. 5 ESD = SDk * Cc.

    # ---- Residual error ----
    propSd <- fixed(sqrt(0.0235))
    label("Proportional residual SD on the rifampicin blood concentration Cc (fraction)")
    # S2 first $SIGMA 0.0235 FIX (Y1 = IPRED + IPRED*EPS(1), PK records only),
    # carried from the upstream mouse popPK; no PK was sampled in this study.

    addSd <- sqrt(0.217)
    label("Additive residual SD on natural-log CFU/lungs")
    # Table 1: additive residual variance on log scale = 0.217, RSE 12.0%
    # (S2 second $SIGMA 0.216544). SD = sqrt(0.217). Log transformation on both
    # sides in NONMEM, so the additive error is on the natural-log CFU scale.
  })

  model({
    # ---- Rifampicin PK (one-compartment, first-order absorption) ----
    # Auto-induction: CL/F steps up after ltclchange hours since infection.
    cl <- exp(lcl)
    if (t >= exp(ltclchange)) {
      cl <- exp(lcl_late)
    }
    # Volume: the single lowest-dose cohort carries a larger apparent volume.
    vc <- exp(lvc)
    if (DOSE_LOW == 1) {
      vc <- exp(lvc_lowdose)
    }
    kel <- cl / vc
    # Elimination-anchored absorption: ka is pinned a fixed amount above kel
    # (S2 TVKA = KE + THETA(5)), enforcing ka > kel so the fit cannot settle
    # into the flip-flop solution.
    ka <- kel + exp(lka_excess)

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    f(depot) <- exp(lfdepot)
    # Rifampicin blood concentration. Dose is in mg/kg and V/F in mL/kg, so
    # central / vc is mg/mL; * 1000 converts to ug/mL (= mg/L).
    Cc <- 1000 * central / vc

    # ---- MTP rate constants ----
    kg <- exp(lkg)
    kfslin <- exp(lkfslin)
    ksf <- exp(lksf)
    kfn <- exp(lkfn)
    ksn <- exp(lksn)
    kns <- exp(lkns)
    kfs <- kfslin * t # time-dependent F->S transfer (Eq. 1)

    # ---- Rifampicin exposure-response (Eqs. 3-5) ----
    fgimax <- exp(lemax_fg_rif)
    fgic50 <- exp(lec50_fg_rif)
    fdemax <- exp(lemax_fd_rif)
    fdec50 <- exp(lec50_fd_rif)
    sdk <- exp(lksd_rif)

    ifg <- fgimax * Cc / (fgic50 + Cc) # inhibition of F growth (Eq. 3)
    efd <- fdemax * Cc^gamma_fd_r / (fdec50^gamma_fd_r + Cc^gamma_fd_r) # death of F (Eq. 4)
    esd <- sdk * Cc # death of S (Eq. 5)

    # ---- Three-state bacterial system (Eqs. 6-8) ----
    # Initial fast-multiplying count depends on treatment group; slow- and
    # non-multiplying inocula were not supported by the data and are 0.
    fbugs(0) <- exp(lf0) * TRT_RIFAMPICIN + exp(lf0_untreated) * (1 - TRT_RIFAMPICIN)
    sbugs(0) <- 0
    nbugs(0) <- 0

    d/dt(fbugs) <- kg * (1 - ifg) * fbugs - efd * fbugs -
      kfs * fbugs + ksf * sbugs - kfn * fbugs
    d/dt(sbugs) <- kfs * fbugs - ksf * sbugs + kns * nbugs -
      ksn * sbugs - esd * sbugs
    d/dt(nbugs) <- kfn * fbugs + ksn * sbugs - kns * nbugs

    # ---- Observation: natural log of culturable CFU (F + S) ----
    # S2 $ERROR: IPRED = LOG(A(3)+A(4)); non-multiplying N bacteria do not
    # form colonies on solid media. Floor at 1 CFU to keep the log finite
    # when drug effect drives the culturable total below one organism.
    culturable <- fbugs + sbugs
    if (culturable < 1) {
      culturable <- 1
    }
    log_cfu <- log(culturable)
    Cc ~ prop(propSd)
    log_cfu ~ add(addSd)
  })
}
