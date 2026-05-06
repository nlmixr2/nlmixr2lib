Svensson_2018_rifampicin <- function() {
  description <- "One-compartment population PK model for high-dose oral rifampicin in adult pulmonary tuberculosis patients (Svensson 2018, HIGHRIF1), with closed-form transit-compartment absorption (mean transit time and Erlang shape estimated), Michaelis-Menten clearance scaled by an auto-induced enzyme turnover compartment, fat-free-mass allometric scaling on Vmax (0.75) and central volume (1.0), and a saturable dose-dependent bioavailability anchored at the 450 mg reference dose."
  reference <- paste(
    "Svensson RJ, Aarnoutse RE, Diacon AH, Dawson R, Gillespie SH,",
    "Boeree MJ, Simonsson USH. (2018). A population pharmacokinetic",
    "model incorporating saturable pharmacokinetics and autoinduction",
    "for high rifampicin doses. Clin Pharmacol Ther 103(4):674-683.",
    "doi:10.1002/cpt.778.",
    "DDMORE Foundation Model Repository: DDMODEL00000244.",
    sep = " "
  )
  vignette <- "Svensson_2018_rifampicin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")
  ddmore_id    <- "DDMODEL00000244"
  replicate_of <- NULL

  covariateData <- list(
    FFM = list(
      description        = "Fat-free mass (Janmahasatian formula).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on a 70 kg reference fat-free mass: Vmax scales (FFM/70)^0.75, V2 scales (FFM/70)^1. The Svensson 2018 HIGHRIF1 cohort had a typical FFM well below 70 kg (the bundled simulated dataset uses FFM = 34.87 kg corresponding to WT = 46.5 kg / HT = 1.78 m / male). Both scaling exponents are theory-based (Anderson & Holford), not estimated.",
      source_name        = "FFM"
    ),
    DOSE = list(
      description        = "Per-record administered rifampicin dose (mg) used as the input to the saturable bioavailability function f_dose(DOSE) = 1 + femax*(DOSE-450)/(fed50+(DOSE-450)).",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference dose 450 mg (the standard adult rifampicin dose) where f_dose = 1 by construction. Calibrated by Svensson 2018 against HIGHRIF1 cohorts at 10, 20, and 35 mg/kg (~600, 1200, 2100 mg for a 60 kg adult); behaviour below 450 mg or above ~3000 mg is unconstrained and should not be extrapolated. Per record (= per dose event); assumed constant within an occasion.",
      source_name        = "DOSE"
    ),
    OCC = list(
      description        = "Integer-valued dosing-occasion indicator for the IOV multiplexers (1 = day 7, 2 = day 14).",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Two-occasion crossover within subject. Decomposed inside model() into binary indicators oc1 / oc2 that select per-occasion etas on bioavailability, MTT, KM, and KA. The Svensson 2018 design sampled dense profiles on study days 7 and 14 of repeated daily dosing.",
      source_name        = "OCC"
    )
  )

  population <- list(
    n_subjects     = 83L,
    n_studies      = 1L,
    age_range      = "Adults (HIGHRIF1 enrolled adult pulmonary TB patients; per-subject demographics not transcribed from the publication, which is not on disk in this worktree).",
    weight_range   = "Adults (HIGHRIF1; the bundled Simulated_Rif_PK_data.csv carries WT = 46.5 kg as the smoke-test cohort's single weight).",
    disease_state  = "Adult pulmonary tuberculosis on high-dose rifampicin-containing antitubercular therapy.",
    dose_range     = "Oral rifampicin 600, 1200, and 2100 mg once daily (the HIGHRIF1 dose-escalation cohorts at 10, 20, and 35 mg/kg for ~60 kg adults). Dense PK sampling at study days 7 and 14 of repeat daily dosing.",
    regions        = "South Africa and Tanzania (PanACEA HIGHRIF1 trial sites).",
    notes          = "The Svensson 2018 publication itself is not on disk in this worktree, so per-subject demographics here are reproduced from the DDMODEL00000244 RDF abstract and the .mod $INPUT column comments rather than the paper's Table 1. The bundle's Simulated_Rif_PK_data.csv is a single-subject smoke-test cohort (ID 1, WT = 46.5 kg, FFM = 34.87 kg, male, 600 mg QD x ~7 days at occasion 1 then x ~7 days at occasion 2, dense sampling at days 7 and 14). N_subjects = 83 is the total in the listing's 'TOT. NO. OF INDIVIDUALS' field; the RDF model-has-description-long abstract describes the structural model but does not enumerate per-subject demographics."
  )

  ini({
    # Structural parameters - DDMODEL00000244 Output_real_Rif_PK.lst FINAL PARAMETER ESTIMATE
    # block (THETA vector). Reference fat-free mass is 70 kg (.mod $PK lines 51-54: ALLMCL =
    # (FFM/70)^0.75 on Vmax/CL, ALLMV = (FFM/70)^1 on V2). The OFV is OBJV = -1053.189 and the
    # listing reports MINIMIZATION TERMINATED DUE TO ROUNDING ERRORS (NSIG = 2.9 vs 3 required);
    # parameters are stationary across iterations 6 / 8 and the curator has updated the .mod
    # $THETA values to match the .lst final estimates bit-for-bit. See vignette Errata.
    lvmax  <- log(525)    ; label("Michaelis-Menten Vmax of clearance at FFM = 70 kg (mg/h)")              # THETA(1) FINAL = 5.25E+02
    lkm    <- log(35.3)   ; label("Michaelis-Menten Km of clearance (mg/L)")                                # THETA(2) FINAL = 3.53E+01
    lvc    <- log(87.2)   ; label("Apparent central volume of distribution V2 at FFM = 70 kg (L)")          # THETA(3) FINAL = 8.72E+01
    lka    <- log(1.77)   ; label("First-order absorption rate constant from depot to central, ka (1/h)")   # THETA(4) FINAL = 1.77E+00
    emax   <- 1.16        ; label("Auto-induction Emax (unitless multiplier on enzyme synthesis)")          # THETA(5) FINAL = 1.16E+00
    ec50   <- 0.0699      ; label("Auto-induction EC50 on rifampicin Cc (mg/L)")                            # THETA(6) FINAL = 6.99E-02
    kenz   <- 0.00603     ; label("First-order enzyme turnover rate constant kenz (1/h)")                   # THETA(7) FINAL = 6.03E-03
    lmtt   <- log(0.513)  ; label("Mean transit time through the absorption transit chain, MTT (h)")        # THETA(8) FINAL = 5.13E-01
    lnn    <- log(23.8)   ; label("Erlang transit-chain shape parameter NN (number of transit compartments)") # THETA(9) FINAL = 2.38E+01
    femax  <- 0.504       ; label("Dose-dependent bioavailability Emax (fraction; reference dose 450 mg)")  # THETA(10) FINAL = 5.04E-01
    fed50  <- 67          ; label("Dose-dependent bioavailability EC50 offset above 450 mg (mg)")           # THETA(11) FINAL = 6.70E+01

    # Inter-individual variability. NONMEM $OMEGA structure (lst lines 195-213 / 633-682):
    #   $OMEGA BLOCK(2)             — correlated KM (ETA1) - VMAX (ETA2)
    #   $OMEGA  diag                — V2 (ETA3), MTT (ETA4), NN (ETA5), KA (ETA6)
    #   $OMEGA BLOCK(1) + SAME      — IOV in F (BIO) over the 2 occasions  (ETA7, ETA8)
    #   $OMEGA BLOCK(1) + SAME      — IOV in MTT                            (ETA9, ETA10)
    #   $OMEGA BLOCK(1) + SAME      — IOV in KM                             (ETA11, ETA12)
    #   $OMEGA BLOCK(1) + SAME      — IOV in KA                             (ETA13, ETA14)
    # The .mod's KM = TVKM*EXP(ETA(1)+IOVKM) puts ETA(1) on lkm; correspondingly the BLOCK(2)
    # cross-term covar(ETA1,ETA2) carries to corr(etalkm, etalvmax). nlmixr2 syntax keeps the
    # 'l' prefix on the eta name to match the transformed parameter.
    etalkm + etalvmax ~ c(0.128, 0.0418, 0.0901)  # OMEGA(1,1)/(2,1)/(2,2) FINAL = 1.28E-01 / 4.18E-02 / 9.01E-02
    etalvc  ~ 0.00618    # OMEGA(3,3) FINAL = 6.18E-03 — IIV log-V2 variance
    etalmtt ~ 0.146      # OMEGA(4,4) FINAL = 1.46E-01 — IIV log-MTT variance
    etalnn  ~ 0.607      # OMEGA(5,5) FINAL = 6.07E-01 — IIV log-NN variance
    etalka  ~ 0.114      # OMEGA(6,6) FINAL = 1.14E-01 — IIV log-KA variance

    # IOV. NONMEM $OMEGA BLOCK(1) followed by BLOCK(1) SAME re-uses the same single-element
    # variance across the two occasions; the FINAL estimate of each shared variance is one
    # element only. nlmixr2 has no SAME shortcut, so each occasion gets its own eta with the
    # second occasion's variance fix(...)-pinned to the first (matching the Jonsson_2011_ethambutol
    # and Xie_2019_agomelatine patterns).
    etaiov_bio_1 ~ 0.0248        # OMEGA(7,7)  FINAL = 2.48E-02 — estimated occasion-1 IOV in bioavailability
    etaiov_bio_2 ~ fix(0.0248)   # OMEGA(8,8)  fixed equal to OMEGA(7,7) per $OMEGA BLOCK(1) SAME
    etaiov_mtt_1 ~ 0.318         # OMEGA(9,9)  FINAL = 3.18E-01 — estimated occasion-1 IOV in MTT
    etaiov_mtt_2 ~ fix(0.318)    # OMEGA(10,10) fixed equal to OMEGA(9,9) per $OMEGA BLOCK(1) SAME
    etaiov_km_1  ~ 0.0355        # OMEGA(11,11) FINAL = 3.55E-02 — estimated occasion-1 IOV in KM
    etaiov_km_2  ~ fix(0.0355)   # OMEGA(12,12) fixed equal to OMEGA(11,11) per $OMEGA BLOCK(1) SAME
    etaiov_ka_1  ~ 0.0985        # OMEGA(13,13) FINAL = 9.85E-02 — estimated occasion-1 IOV in KA
    etaiov_ka_2  ~ fix(0.0985)   # OMEGA(14,14) fixed equal to OMEGA(13,13) per $OMEGA BLOCK(1) SAME

    # Residual error. The .mod uses log-transformed observations with $ERROR Y = IPRED + EPS(1)
    # where IPRED = LOG(A(2)/S2 + 1e-5); on the back-transformed linear scale this is
    # proportional with proportional-SD = sqrt(SIGMA(1,1)) = sqrt(0.0555) = 0.2356 (NONMEM
    # 'additive on log-scale' ≡ proportional in nlmixr2's linear space — see naming-conventions.md
    # § Residual error). The .mod additionally implements the Beal M3 method for BLOQ data
    # (LLOQ = log(0.13 mg/L), F_FLAG = 1 / Y = PHI((LLOQ - IPRED)/SD)); M3 is an estimation-time
    # construct, not part of the structural model, so it is not carried into the nlmixr2 model.
    propSd <- 0.2356  ; label("Proportional residual error (SD on log-Cc scale ≡ fraction in linear space)") # SIGMA(1,1) FINAL = 5.55E-02 = 0.2356^2
  })

  model({
    # Decompose the integer-valued occasion column into binary indicators for the four IOV
    # multiplexers (matches the .mod IF (OCC.EQ.1)/ELSE assignments to ETA(7..14) in $PK).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)

    iov_bio <- oc1 * etaiov_bio_1 + oc2 * etaiov_bio_2
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2
    iov_km  <- oc1 * etaiov_km_1  + oc2 * etaiov_km_2
    iov_ka  <- oc1 * etaiov_ka_1  + oc2 * etaiov_ka_2

    # Individual PK parameters with allometric scaling on a 70 kg fat-free-mass reference
    # (theory-based exponents: 0.75 on Vmax/CL, 1.0 on V2; .mod $PK lines 51-54).
    vmax <- exp(lvmax + etalvmax)            * (FFM / 70)^0.75
    km   <- exp(lkm   + etalkm  + iov_km)
    vc   <- exp(lvc   + etalvc)              * (FFM / 70)
    ka   <- exp(lka   + etalka  + iov_ka)
    mtt  <- exp(lmtt  + etalmtt + iov_mtt)
    nn   <- exp(lnn   + etalnn)

    # Saturable dose-dependent bioavailability anchored at DOSE = 450 mg (.mod TVBIO line 96):
    #   TVBIO = 1*(1 + FEMAX*(DOSE-450)/(FED50+(DOSE-450))); BIO = TVBIO * EXP(IOVBIO).
    # At DOSE = 450 mg, TVBIO = 1 (reference). The HIGHRIF1 calibration range is 600-2100 mg;
    # do not extrapolate below 450 mg or above ~3000 mg.
    bio <- (1 + femax * (DOSE - 450) / (fed50 + (DOSE - 450))) * exp(iov_bio)

    # Concentration in plasma drives the Michaelis-Menten elimination rate and the
    # auto-induction EFF term. Cc = central / vc, units mg / L.
    Cc  <- central / vc
    eff <- emax * Cc / (ec50 + Cc)

    # Closed-form Erlang transit-compartment input (Wilkins / Savic) replaces the .mod's
    # $DES verbatim FORTRAN
    #   CUMUL = LOG(BIO*PD) + LOG(KTR) - L
    #   DADT(1) = EXP(CUMUL + NN*LOG(KTR*TEMPO) - KTR*TEMPO) - KA*A(1)
    # where L is the Stirling approximation of log(gamma(NN+1)) and PD is the most recent
    # AMT. rxode2's transit(n, mtt, bio) returns precisely the same gamma-density input
    # (using podo() and tad() internally for the most recent dose amount and time-after-dose).
    # The .mod sets F1 = 0 to disable normal dose accumulation in the depot, so transit() is
    # the sole input to depot — preserve that with `f(depot) <- 0` below.
    d/dt(depot)   <-  transit(nn, mtt, bio) - ka * depot
    d/dt(central) <-  ka * depot - vmax * Cc / (km + Cc) * enzyme
    d/dt(enzyme)  <-  kenz * (1 + eff) - kenz * enzyme

    # Initial conditions match the .mod $PK A_0 lines:
    #   A_0(2) = 0.0001  ; numerical-stabilisation seed for M-M denominator (Cp ≈ 1e-6 mg/L)
    #   A_0(3) = 1       ; baseline auto-induction enzyme amount (steady state with no drug)
    central(0) <- 0.0001
    enzyme(0)  <- 1

    # F1 = 0 in the .mod $PK: the dose targets the depot but does not accumulate there
    # directly; the transit() closed form is the only depot input. f(depot) <- 0 enforces this.
    f(depot) <- 0

    # Proportional residual error (NONMEM 'additive on log-scale' Y = LOG(F) + EPS(1) maps to
    # nlmixr2's prop() in linear space).
    Cc ~ prop(propSd)
  })
}
