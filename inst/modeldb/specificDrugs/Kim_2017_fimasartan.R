Kim_2017_fimasartan <- function() {
  description <- "Population PK-PD model for fimasartan (an angiotensin II receptor blocker) in healthy adult Korean men and men with mild or moderate hepatic impairment (Kim 2017). Plasma fimasartan is described by a 2-compartment model with parallel mixed-input absorption: a first-order arm with rate Ka and absorption lag time LAG (fraction F1 = (1 - alpha) * F of the dose) running in parallel with a zero-order arm of virtual duration D2 (fraction F2 = alpha * F of the dose), where the total relative bioavailability F is fixed at 0.18 in healthy subjects (Kim 2008) and incremented to 0.18 + IL1 in mild and 0.18 + IL2 in moderate hepatic impairment to capture the markedly higher Cmax observed in cirrhotic patients via reduced first-pass extraction and intrahepatic shunting. The PD model describes systolic and diastolic blood pressures as indirect-response (turnover) compartments with zero-order synthesis Kin inhibited by fimasartan via a sigmoid-Imax function E(C) = 1 - Emax * Cc / (EC50 + Cc) and first-order loss Kout = Kin / Base; the steady-state baseline rides a fixed cosinor circadian rhythm Bsl(t) = MESOR * (1 + Amp1% * cos(2*pi*(t - AC1)/24) + Amp2% * cos(2*pi*(t - AC2)/12)) with amplitudes and phases inherited from Park 2014 (healthy Korean reference). EC50 is stratified by hepatic-impairment severity: for SBP, healthy versus any-impairment pooled (mild + moderate); for DBP, healthy + mild versus moderate alone, reflecting the contrasting impact of hepatic dysfunction on the two pressure outputs."
  reference   <- "Kim CO, Jeon S, Han S, Hong T, Park MS, Yoon Y-R, Yim D-S. Decreased potency of fimasartan in liver cirrhosis was quantified using mixed-effects analysis. Transl Clin Pharmacol. 2017;25(1):43-49. doi:10.12793/tcp.2017.25.1.43. Bioavailability in healthy subjects (F = 0.18) inherited from Kim TH et al. (Eur J Drug Metab Pharmacokinet 2010). Circadian-rhythm amplitudes and phase shifts (Table 3) inherited from the Park 2014 cosinor model of blood-pressure rhythm in healthy Koreans."
  vignette    <- "Kim_2017_fimasartan"
  units       <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    HEPIMP_MILD = list(
      description        = "Mild hepatic impairment indicator (Child-Pugh Class A, score 5-6)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal hepatic function; mutually exclusive with HEPIMP_MODSEV)",
      notes              = "Time-fixed. Mild = Child-Pugh score 5-6 per Kim 2017 Methods 'Study subjects'. Used in both the PK model (additive increment IL1 on relative bioavailability F: F_mild = 0.18 + IL1) and in the SBP PD model (combined with HEPIMP_MODSEV as the any-impairment indicator for the EC50_AB switch; mild and moderate impairment share a common EC50_AB on SBP). Not used directly in the DBP PD model -- the DBP EC50 switch pools healthy + mild against moderate alone, so HEPIMP_MILD contributes only via the absence of HEPIMP_MODSEV.",
      source_name        = "Mild hepatic impairment"
    ),
    HEPIMP_MODSEV = list(
      description        = "Moderate-or-severe hepatic impairment indicator (Child-Pugh Class B or C). In Kim 2017 the moderate cohort had Child-Pugh score 7-9 (Class B); no severe (Class C) subjects enrolled.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal or mild hepatic function; mutually exclusive with HEPIMP_MILD)",
      notes              = "Time-fixed. Moderate = Child-Pugh score 7-9 per Kim 2017 Methods 'Study subjects'; the column carries the canonical HEPIMP_MODSEV semantics so a downstream user enrolling a severe (Child-Pugh C) cohort would set this to 1 as well -- Kim 2017's reported coefficients IL2 (PK), EC50_AB minus EC50_H (SBP), and EC50_B minus EC50_HA (DBP) are calibrated to the moderate cohort only and should be treated as moderate-specific by users simulating severe disease. Used in the PK model (additive increment IL2 on F: F_moderate = 0.18 + IL2), in the SBP PD model (combined with HEPIMP_MILD as the any-impairment switch for EC50_AB), and in the DBP PD model (alone as the switch for EC50_B versus EC50_HA).",
      source_name        = "Moderate hepatic impairment"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 18L,
    n_studies        = 1L,
    n_observations   = "288 plasma fimasartan concentrations and 180 paired blood-pressure observations (Kim 2017 Results 'Study subjects').",
    age_range        = "26-56 years (means 48.8, 43.2, 48.2 in healthy / mild / moderate cohorts; Kim 2017 Table 1)",
    weight_range     = "57.0-85.0 kg (means 71.8, 70.3, 65.6 in healthy / mild / moderate cohorts; Kim 2017 Table 1)",
    sex_female_pct   = 0,
    race_ethnicity   = c(Korean = 100),
    disease_state    = "Adult Korean men in three balanced groups of six: healthy controls, mild hepatic impairment (Child-Pugh score 5-6), and moderate hepatic impairment (Child-Pugh score 7-9). All moderate-impairment subjects had chronic liver cirrhosis; two of six mild-impairment subjects had cirrhosis (Kim 2017 Discussion).",
    dose_range       = "Single oral 120 mg fimasartan",
    regions          = "Republic of Korea (Yonsei University College of Medicine, Severance Hospital; Catholic University of Korea, Seoul St. Mary's Hospital; Kyungpook National University Hospital, Daegu).",
    hepatic_function = "Six normal, six Child-Pugh Class A (mild), six Child-Pugh Class B (moderate). No Child-Pugh Class C (severe) subjects.",
    notes            = "Open-label, single-dose, parallel-group study. Plasma fimasartan sampled at 0, 0.25, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 10, 12, 16, 24, 32, 48 h post-dose. SBP and DBP sampled at 0, 1, 2, 3, 4, 8, 12, 24, 32, 48 h post-dose, each measurement preceded by a >= 5-minute seated rest. Covariates of age, weight, albumin, bilirubin, AST, ALT, gamma-GT, creatinine, and prothrombin time (INR) were screened but not retained in the final structural PK model (Kim 2017 Discussion); the only covariate effect retained is the hepatic-impairment indicator on relative bioavailability."
  )

  ini({
    # ====================================================================
    # Population PK model -- Kim 2017 Table 2 (final estimates)
    # Two-compartment disposition with mixed zero-and-first-order parallel
    # absorption with a first-order lag. All structural PK THETAs are
    # estimated point values (RSEs and bootstrap CIs reported in Table 2);
    # the bioavailability anchor F = 0.18 in healthy subjects is fixed
    # from Kim 2008 (cited as ref [4] in Kim 2017 'Population pharmaco-
    # kinetic analysis').
    # ====================================================================
    lcl   <- log(27.0);    label("Apparent clearance CL/F at reference covariates (L/h)")                              # Kim 2017 Table 2: CL = 27.0 L/h, RSE 13.3 %
    lvc   <- log(48.7);    label("Apparent central volume of distribution V2/F (L)")                                   # Kim 2017 Table 2: V2 = 48.7 L, RSE 20.4 %
    lvp   <- log(46.5);    label("Apparent peripheral volume of distribution V3/F (L)")                                # Kim 2017 Table 2: V3 = 46.5 L, RSE 11.1 %
    lka   <- log(0.319);   label("First-order absorption rate constant Ka (1/h)")                                      # Kim 2017 Table 2: Ka = 0.319 1/h, RSE 16.3 %
    lq    <- log(3.40);    label("Apparent inter-compartmental clearance Q/F (L/h)")                                   # Kim 2017 Table 2: Q = 3.40 L/h, RSE 12.4 %
    ld2   <- log(0.583);   label("Virtual duration of the zero-order absorption arm D2 (h)")                           # Kim 2017 Table 2: D2 = 0.583 h, RSE 9.8 %
    ltlag  <- log(2.0);     label("Absorption lag time of the first-order arm LAG (h)")                                 # Kim 2017 Table 2: LAG = 2.0 h, bootstrap median 2.0 (1.4-2.5); RSE 0.1 % in the table appears to be a typo given the bootstrap CI -- see vignette Errata
    logitalpha <- qlogis(0.642); label("Logit of the proportionality constant alpha for the zero-order absorption fraction (F2 = alpha * F)")  # Kim 2017 Table 2: alpha = 0.642, RSE 7.4 % -- bounded to [0,1] via logit so the simulated zero-order fraction respects 0 <= F2/F <= 1

    # Bioavailability anchor and hepatic-impairment increments. F (total)
    # = 0.18 (healthy, fixed) + IL1 * HEPIMP_MILD + IL2 * HEPIMP_MODSEV.
    # F1 = (1 - alpha) * F is the first-order absorbed fraction routed via
    # depot; F2 = alpha * F is the zero-order absorbed fraction routed
    # directly into central as a zero-order infusion of duration D2. Note
    # F = 0.18 + IL2 = 1.076 > 1 in the moderate-impairment cohort, which
    # is the paper-stated relative-bioavailability scaling and reflects
    # the ~6x higher absolute exposure observed in that cohort via
    # reduced first-pass extraction / intrahepatic shunting (Kim 2017
    # Discussion). The interpretation is therefore a paper-defined
    # 'relative bioavailability scalar' rather than a strict fraction-
    # bounded-to-1 quantity.
    fhealthy <- fixed(0.18);  label("Reference bioavailability F in healthy subjects (unitless; fixed from Kim 2008)")  # Kim 2017 Methods 'Population pharmacokinetic analysis' citing Kim 2008
    il1      <- 0.0873;       label("Additive increment on F for mild hepatic impairment (unitless)")                   # Kim 2017 Table 2: IL1 = 0.0873, RSE 142.1 % (95 % CI 0.001-0.668)
    il2      <- 0.896;        label("Additive increment on F for moderate hepatic impairment (unitless)")               # Kim 2017 Table 2: IL2 = 0.896, RSE 19.7 % (95 % CI 0.607-1.481)

    # IIV on the PK parameters. Kim 2017 Table 2 reports between-subject
    # variability as CV%; convert to log-scale variance via omega^2 =
    # log(1 + CV^2) per the lognormal IIV convention. None of the IIV
    # entries was reported as fixed.
    etalcl       ~ log(1 + 0.399^2)    # Kim 2017 Table 2: omega_CL/F = 39.9 % CV, RSE 15.8 %
    etalvc       ~ log(1 + 1.214^2)    # Kim 2017 Table 2: omega_V2/F = 121.4 % CV, RSE 41.3 %
    etalka       ~ log(1 + 0.635^2)    # Kim 2017 Table 2: omega_Ka  = 63.5 % CV, RSE 18.0 %
    etalogitalpha ~ log(1 + 0.695^2)   # Kim 2017 Table 2: omega_alpha = 69.5 % CV, RSE 43.8 %

    # Residual error -- Kim 2017 Table 2 (proportional + small fixed
    # additive). The additive component is reported as 'fix' (held at
    # 0.0001 in the source units of ng/mL).
    addSd  <- fixed(0.0001); label("Additive residual error on plasma fimasartan Cc (ng/mL; fixed per Kim 2017 Table 2)")  # Kim 2017 Table 2: sigma_add = 0.0001 (fix)
    propSd <- 0.354;         label("Proportional residual error on plasma fimasartan Cc (fraction)")                        # Kim 2017 Table 2: sigma_prop = 0.354, RSE 8.8 %

    # ====================================================================
    # Systolic blood pressure (SBP) PD model -- Kim 2017 Table 4
    # Indirect-response (turnover) model with Imax-style fimasartan
    # inhibition of synthesis. EC50 differs between healthy and any-
    # hepatic-impairment subjects (mild and moderate pooled per Kim 2017
    # Table 4 row 'EC50,A+B' -- the SBP fit grouped mild + moderate
    # together).
    # ====================================================================
    lkin_sbp    <- log(90.3);  label("Zero-order synthesis rate Kin for SBP (mmHg/h)")                                          # Kim 2017 Table 4: Kin_SBP = 90.3 mmHg/h, RSE 17.7 %
    logitemax_sbp <- qlogis(0.213); label("Logit of maximum fractional inhibition Emax for SBP (unitless; bounded to [0,1])")   # Kim 2017 Table 4: Emax_SBP = 21.3 %, RSE 5.8 % -- logit transform enforces the paper-stated 0 <= Emax <= 1 constraint
    lbase_sbp   <- log(131.0); label("Predose typical SBP value Base (mmHg)")                                                  # Kim 2017 Table 4: Base_SBP = 131.0 mmHg, RSE 2.0 %
    lec50_sbp   <- log(2.28);  label("EC50 of fimasartan for SBP inhibition in healthy subjects (ng/mL)")                       # Kim 2017 Table 4: EC50_H_SBP = 2.28 ng/mL, RSE 20.7 %
    e_hi_any_ec50_sbp <- log(9.19 / 2.28); label("Log-multiplicative effect of any hepatic impairment (mild or moderate, pooled) on EC50 for SBP (unitless)")  # Derived from Kim 2017 Table 4 EC50_A+B_SBP / EC50_H_SBP = 9.19 / 2.28 = 4.03; equivalent to a single coefficient applied jointly to mild and moderate per Kim 2017 'EC50,A+B'

    etalbase_sbp ~ log(1 + 0.053^2)  # Kim 2017 Table 4: omega_Base_SBP = 5.3 % CV, RSE 36.3 % (bootstrap median 5.1 %, consistent with the point estimate)
    propSd_SBP   <- 0.063;     label("Proportional residual error on SBP (fraction)")                                         # Kim 2017 Table 4: sigma_prop_SBP = 0.063, RSE 6.2 %

    # ====================================================================
    # Diastolic blood pressure (DBP) PD model -- Kim 2017 Table 4
    # Same indirect-response structure as SBP. EC50 grouping differs:
    # healthy + mild pooled (EC50_H+A) vs moderate alone (EC50_B), per
    # Kim 2017 Table 4. Encoded as a base EC50 (healthy + mild reference)
    # plus a log-multiplicative covariate effect for the moderate-
    # impairment switch, so that exp(lec50_dbp) reproduces EC50_H+A and
    # exp(lec50_dbp + e_hepmodsev_ec50_dbp) reproduces EC50_B.
    # ====================================================================
    lkin_dbp     <- log(33.1); label("Zero-order synthesis rate Kin for DBP (mmHg/h)")                                          # Kim 2017 Table 4: Kin_DBP = 33.1 mmHg/h, RSE 19.7 %
    logitemax_dbp <- qlogis(0.338); label("Logit of maximum fractional inhibition Emax for DBP (unitless; bounded to [0,1])")   # Kim 2017 Table 4: Emax_DBP = 33.8 %, RSE 8.7 %
    lbase_dbp    <- log(82.3); label("Predose typical DBP value Base (mmHg)")                                                  # Kim 2017 Table 4: Base_DBP = 82.3 mmHg, RSE 3.3 %
    lec50_dbp    <- log(4.82); label("EC50 of fimasartan for DBP inhibition in healthy + mild-impairment subjects (ng/mL)")     # Kim 2017 Table 4: EC50_H+A_DBP = 4.82 ng/mL, RSE 40.5 %
    e_hepmodsev_ec50_dbp <- log(47.3 / 4.82); label("Log-multiplicative effect of moderate hepatic impairment on EC50 for DBP (unitless)")  # Derived from Kim 2017 Table 4 EC50_B_DBP / EC50_H+A_DBP = 47.3 / 4.82 = 9.81; reproduces the moderate-only group's EC50_B when HEPIMP_MODSEV = 1

    # IIV on the DBP PD parameters. The Table 4 point estimate for
    # omega_Base_DBP is reported as 1.25 % CV but the bootstrap median is
    # 10.8 % (95 % CI 6.6-14.5 %); the point estimate appears to be a
    # decimal-placement typo (likely 12.5 % CV) and the bootstrap median
    # is used here as the more reliable value. Documented in the vignette
    # 'Errata and deviations' section. The IIV on EC50 is applied
    # multiplicatively on the (already covariate-modified) DBP EC50 so a
    # single eta describes between-subject variability irrespective of
    # the subject's hepatic-impairment group.
    etalbase_dbp ~ log(1 + 0.108^2)   # Kim 2017 Table 4: omega_Base_DBP bootstrap median 10.8 % CV (point estimate 1.25 % CV in Table 4 appears to be a decimal-place typo; bootstrap value preferred)
    etalec50_dbp ~ log(1 + 0.568^2)   # Kim 2017 Table 4: omega_EC50_DBP = 56.8 % CV, RSE 55.3 % -- shared between EC50_H+A and EC50_B groups since both share a structural EC50 with a covariate switch

    addSd_DBP  <- 6.27;        label("Additive residual error on DBP (mmHg)")                                                  # Kim 2017 Table 4: sigma_add_DBP = 6.27 mmHg, RSE 6.8 %
    propSd_DBP <- fixed(1e-4); label("Proportional residual error on DBP (fraction; fixed at 0.0001 per Kim 2017 Table 4)")  # Kim 2017 Table 4: sigma_prop_DBP = 0.0001 (fix)
  })

  model({
    # ====================================================================
    # Fixed circadian-rhythm constants for the cosinor baseline. Inherited
    # from a prior cosinor analysis of healthy Korean blood pressure
    # (Park 2014, cited as ref [9] in Kim 2017). Amp_i are reported as
    # percent of MESOR; AC_i are phase shifts in hours. Kim 2017 Table 3.
    # The first cosine term is the 24-hour fundamental; the second is
    # the 12-hour harmonic (period 12 h).
    # ====================================================================
    amp1_sbp <- -10.2  # Kim 2017 Table 3: amplitude 1st cosine term SBP -10.2 %
    amp2_sbp <-  4.47  # Kim 2017 Table 3: amplitude 2nd cosine term SBP  4.47 %
    ac1_sbp  <- -3.44  # Kim 2017 Table 3: phase shift 1st cosine term SBP -3.44 h
    ac2_sbp  <-  2.42  # Kim 2017 Table 3: phase shift 2nd cosine term SBP  2.42 h
    amp1_dbp <- -13.8  # Kim 2017 Table 3: amplitude 1st cosine term DBP -13.8 %
    amp2_dbp <-  6.39  # Kim 2017 Table 3: amplitude 2nd cosine term DBP  6.39 %
    ac1_dbp  <- -3.56  # Kim 2017 Table 3: phase shift 1st cosine term DBP -3.56 h
    ac2_dbp  <-  2.28  # Kim 2017 Table 3: phase shift 2nd cosine term DBP  2.28 h

    # ====================================================================
    # Individual PK parameters
    # ====================================================================
    ka    <- exp(lka  + etalka)
    cl    <- exp(lcl  + etalcl)
    vc    <- exp(lvc  + etalvc)
    vp    <- exp(lvp)
    q     <- exp(lq)
    d2    <- exp(ld2)
    lag_t <- exp(ltlag)
    alpha <- expit(logitalpha + etalogitalpha)

    # Total relative bioavailability and the split into first-order
    # (depot) and zero-order (central infusion) fractions.
    f_total <- fhealthy + il1 * HEPIMP_MILD + il2 * HEPIMP_MODSEV
    f1      <- (1 - alpha) * f_total
    f2      <- alpha       * f_total

    # ====================================================================
    # PK ODE system -- 2-compartment with parallel mixed absorption
    # depot   = first-order absorption pool (drains to central with rate ka)
    # central = systemic central compartment (also receives the zero-order
    #           absorption arm directly via a duration-D2 dose record)
    # peripheral1 = first peripheral compartment
    #
    # Each oral administration is encoded by the USER as two parallel
    # dose records in the event table targeting (a) cmt = depot (the
    # first-order arm) and (b) cmt = central (the zero-order arm). The
    # f() multipliers split the dose into F1 and F2 and dur(central)
    # imposes the D2 zero-order infusion duration; lag(depot) imposes
    # the first-order lag.
    # ====================================================================
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - (cl / vc) * central -
                          (q / vc) * central + (q / vp) * peripheral1
    d/dt(peripheral1) <-  (q / vc) * central - (q / vp) * peripheral1

    f(depot)     <- f1
    lag(depot)   <- lag_t
    f(central)   <- f2
    dur(central) <- d2

    # Plasma fimasartan concentration (ng/mL because dose is in mg and
    # vc is in L: mg / L = ug/mL, then * 1000 gives ng/mL).
    Cc <- (central / vc) * 1000

    # ====================================================================
    # Individual PD parameters
    # ====================================================================
    kin_sbp   <- exp(lkin_sbp)
    emax_sbp  <- expit(logitemax_sbp)
    base_sbp  <- exp(lbase_sbp + etalbase_sbp)
    kout_sbp  <- kin_sbp / base_sbp   # Kim 2017 Table 4 footnote b: Kout = Kin / Base

    kin_dbp   <- exp(lkin_dbp)
    emax_dbp  <- expit(logitemax_dbp)
    base_dbp  <- exp(lbase_dbp + etalbase_dbp)
    kout_dbp  <- kin_dbp / base_dbp

    # EC50 stratification by hepatic-impairment severity. SBP encodes the
    # mild + moderate pooled effect as a single log-multiplicative shift
    # on the healthy-reference EC50 (e_hi_any_ec50_sbp). DBP encodes the
    # moderate-only shift relative to the healthy + mild reference
    # (e_hepmodsev_ec50_dbp); the IIV eta on DBP EC50 is shared across
    # both groups via the same multiplicative factor.
    hi_any   <- HEPIMP_MILD + HEPIMP_MODSEV
    ec50_sbp <- exp(lec50_sbp + e_hi_any_ec50_sbp * hi_any)
    ec50_dbp <- exp(lec50_dbp + etalec50_dbp +
                    e_hepmodsev_ec50_dbp * HEPIMP_MODSEV)

    # Sigmoid-Imax drug effect on synthesis (Kim 2017 Methods 'Population
    # pharmacodynamic model'). E(C) approaches 1 - Emax at saturating Cc
    # and equals 1 (no inhibition) at Cc = 0.
    eff_sbp <- 1 - emax_sbp * Cc / (ec50_sbp + Cc)
    eff_dbp <- 1 - emax_dbp * Cc / (ec50_dbp + Cc)

    # ====================================================================
    # Cosinor circadian baseline (Bsl) at the current model time. Bsl(t)
    # rides on a per-subject MESOR_i derived so that Bsl(0) = Base_i:
    #     Bsl(t) / Bsl(0) = circ(t) / circ(0)
    # i.e. Bsl(t) = Base_i * (circ(t) / circ(0)) where
    #     circ(s) = 1 + (Amp1/100) * cos(2 pi (s - AC1) / 24)
    #                 + (Amp2/100) * cos(2 pi (s - AC2) / 12).
    # MESOR_i = Base_i / circ(0). The observed BP is then
    #     BP(t) = Bsl(t) + A(t) - MESOR_i
    # where A(t) is the indirect-response state initialised at MESOR_i
    # so that BP(0) = Base_i. The dimensionless circ(t) factor uses the
    # rxode2 built-in time variable t (hours from dose).
    # ====================================================================
    circ_sbp_t   <- 1 + (amp1_sbp / 100) * cos(2 * pi * (t - ac1_sbp) / 24) +
                        (amp2_sbp / 100) * cos(2 * pi * (t - ac2_sbp) / 12)
    circ_sbp_0   <- 1 + (amp1_sbp / 100) * cos(2 * pi * (0 - ac1_sbp) / 24) +
                        (amp2_sbp / 100) * cos(2 * pi * (0 - ac2_sbp) / 12)
    mesor_sbp    <- base_sbp / circ_sbp_0
    bsl_sbp      <- mesor_sbp * circ_sbp_t

    circ_dbp_t   <- 1 + (amp1_dbp / 100) * cos(2 * pi * (t - ac1_dbp) / 24) +
                        (amp2_dbp / 100) * cos(2 * pi * (t - ac2_dbp) / 12)
    circ_dbp_0   <- 1 + (amp1_dbp / 100) * cos(2 * pi * (0 - ac1_dbp) / 24) +
                        (amp2_dbp / 100) * cos(2 * pi * (0 - ac2_dbp) / 12)
    mesor_dbp    <- base_dbp / circ_dbp_0
    bsl_dbp      <- mesor_dbp * circ_dbp_t

    # ====================================================================
    # Indirect-response ODEs and initial conditions. effect1 and
    # effect2 are the turnover compartments. Without drug, A_ss =
    # Kin / Kout = Base (Kim 2017 Table 4 footnote b). The states are
    # initialised at the per-subject MESOR_i so that the published BP
    # equation BP(t) = Bsl(t) + A(t) - MESOR_i evaluated at t = 0 yields
    # BP(0) = Bsl(0) = Base_i (since Bsl(0) = MESOR_i * circ(0) = Base_i).
    # The minor difference between the IC (MESOR_i) and the long-run
    # steady state (Base_i) generates a small early-time transient that
    # is intrinsic to the paper's stated parameterisation.
    # ====================================================================
    d/dt(effect1) <- kin_sbp * eff_sbp - kout_sbp * effect1
    d/dt(effect2) <- kin_dbp * eff_dbp - kout_dbp * effect2

    effect1(0) <- mesor_sbp
    effect2(0) <- mesor_dbp

    # Observed blood pressures (mmHg). Kim 2017 PD equation:
    #     BP(t) = Bsl(t) + A(t) - MESOR
    SBP <- bsl_sbp + effect1 - mesor_sbp
    DBP <- bsl_dbp + effect2 - mesor_dbp

    Cc  ~ add(addSd) + prop(propSd)
    SBP ~ prop(propSd_SBP)
    DBP ~ add(addSd_DBP) + prop(propSd_DBP)
  })
}
