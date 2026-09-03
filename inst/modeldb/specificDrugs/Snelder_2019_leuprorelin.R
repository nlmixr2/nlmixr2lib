Snelder_2019_leuprorelin <- function() {
  description <- paste(
    "Joint population PK-testosterone-PSA model for the leuprorelin-SR",
    "6-month depot (22.5 or 30 mg single intramuscular dose) in patients",
    "with prostate cancer. PK is a two-compartment disposition model fed by",
    "two parallel release depots: a fast first-order route (Ka1) and a slow",
    "route whose kinetics follow one of two mixture subpopulations -- a",
    "constant first-order rate (Ka2, 40.8% of subjects) or a rate that ramps",
    "linearly with time after the first dose (Ka3_SLP * t, 59.2%).",
    "Testosterone follows the Romero GnRH-agonist model as modified by",
    "Snelder: competitive occupancy of the GnRH receptor by endogenous GnRH",
    "and by leuprorelin drives testosterone synthesis, while the same",
    "occupancy, once it exceeds its baseline value, downregulates receptor",
    "synthesis through a sigmoid Imax term (DR50, nHT). A constant",
    "cyproterone-acetate effect suppresses endogenous GnRH activity during",
    "the flare-prophylaxis window. PSA follows an indirect-response turnover",
    "model whose production rate carries a combined testosterone-dependent",
    "and testosterone-independent (basal) component, Kin_PSA * (1 + EFF)",
    "with EFF a sigmoid Emax function of testosterone."
  )
  reference <- paste(
    "Snelder N, Drenth HJ, Riber Bergmann K, Wood ND, Hibberd M, Scott G.",
    "Population pharmacokinetic-pharmacodynamic modelling of the",
    "relationship between testosterone and prostate specific antigen in",
    "patients with prostate cancer during treatment with leuprorelin.",
    "Br J Clin Pharmacol. 2019;85(6):1247-1259. doi:10.1111/bcp.13891.",
    "Structural equations from Methods 2.5-2.6 (Equations 1, 4, 5, 7);",
    "parameter values from Tables 2-4; the NONMEM control stream in",
    "Supplement 4 resolves the terms the main text leaves implicit.",
    sep = " "
  )
  vignette <- "Snelder_2019_leuprorelin"

  # Doses are given in mg (22.5 / 30 mg depots) and volumes in L, so
  # `central / vc` is mg/L; the *1e6 in model() converts to pg/mL, the unit
  # of the leuprorelin assay (LOQ 16 pg/mL) and of Kd (4.36 pg/mL).
  # Testosterone is carried in ng/dL and PSA in ng/mL, matching Tables 3-4.
  units <- list(time = "h", dosing = "mg", concentration = "pg/mL")

  # Both release depots receive a dose record carrying the FULL nominal dose;
  # f(depot) and f(depot2) then split it (Supplement 4 $PK, F1 = RBIO * FR1 and
  # F3 = RBIO * (1 - FR1)). Declared explicitly because the auto-detection only
  # recognises `depot` and `central`.
  dosing <- c("depot", "depot2")

  # The GnRH-receptor pool RT is a paper-mechanistic state: an apparent total
  # receptor concentration scaled so that its baseline is 1 (Methods 2.5,
  # "The apparent concentration of GnRH receptors at baseline was set to 1 to
  # enable calculation of the fractional receptor occupancy"). It is not a
  # canonical `target` / `complex` TMDD species -- no drug mass is bound out
  # of the central compartment -- so it is declared here rather than renamed.
  paper_specific_compartments <- c("RT")

  compartmentData <- list(
    depot       = list(analyte = "leuprorelin", units = "mg", specimen = "administration site", verified = TRUE),
    depot2      = list(analyte = "leuprorelin", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "leuprorelin", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "leuprorelin", units = "mg", specimen = "serum", verified = TRUE),
    TT          = list(analyte = "testosterone", units = "ng/dL", specimen = "serum", verified = TRUE),
    RT          = list(analyte = "GnRH receptor", units = "unitless", specimen = "not applicable", verified = TRUE),
    PSA         = list(analyte = "prostate-specific antigen", units = "ng/mL", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    CONMED_CYPROTERONE = list(
      description        = "Cyproterone acetate flare-prophylaxis co-medication in effect",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no cyproterone acetate effect in force)",
      notes              = paste(
        "TIME-VARYING. Set to 1 only while the cyproterone acetate (CPA)",
        "effect is assumed to be in force and 0 elsewhere. Study EC403 gave",
        "CPA to at-risk subjects 7 to 3 days before the first leuprorelin",
        "dose (Methods 2.1) and the analysis assumed the effect lasts two",
        "weeks (Methods 2.5: 'Following normal prescribing practice, this CPA",
        "effect was assumed to last for 2 weeks'). Supplement 4 implements",
        "that as a 336 h window opening 168 h before the leuprorelin dose",
        "($DES: START = -FSAMP - 168, END = START + 336, CPA = CPA0 while",
        "START < T < END and ICPA == 1). 90% of the EC403 subjects received",
        "CPA (Results 3.2). The window is supplied as data rather than",
        "hard-coded so that a user can reproduce a different prophylaxis",
        "schedule; set the column to 0 throughout for a CPA-free subject."
      ),
      source_name        = "ICPA / CPA"
    ),
    MIX_RAMP_REL = list(
      description        = "Slow-release mixture class: linearly time-ramping release rate",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (constant first-order slow release, Ka2)",
      notes              = paste(
        "NONMEM $MIXTURE subpopulation indicator (Methods 2.4). 1 = the slow",
        "release depot empties with a rate constant that grows linearly with",
        "time after the first dose (Ka3_SLP * t); 0 = the slow release depot",
        "empties with the constant first-order rate Ka2. Table 2 estimates",
        "the mixture proportion 'prop Pop1' = 0.408 (RSE 18.9%), and Results",
        "3.1 assigns population 1 to first-order release and population 2 to",
        "time-dependent release: '41% and 59% of the subjects were assigned",
        "to population 1 (1st-order release) and population 2 (time",
        "dependent release), respectively.' So a simulated cohort draws",
        "MIX_RAMP_REL ~ Bernoulli(0.592). rxode2 has no $MIXTURE construct,",
        "so the class arrives as data; the proportion is recorded here",
        "rather than in ini() because model() never references it.",
        "NOTE: Supplement 4 numbers the two classes the other way round",
        "(IF(POP.EQ.1) K34 = SLP*TIME/1E6), but that control stream is the",
        "sequential T-PSA run in which the class arrives as a per-subject",
        "data item with run-local numbering; the mechanism-to-fraction",
        "mapping in the main text is explicit and self-consistent."
      ),
      source_name        = "POP / IPOP"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 62L,
    n_studies      = 1L,
    age_range      = "not reported; mean 70.5 years",
    age_median     = "mean 70.5 years",
    sex_female_pct = 0,
    race_ethnicity = c(White = 100),
    disease_state  = paste(
      "Prostate cancer; 83.9% newly diagnosed, the remainder enrolled for",
      "PSA relapse after radical prostatectomy or radiotherapy"
    ),
    dose_range     = "Single 22.5 mg or 30 mg leuprorelin-SR 6-month depot (31 subjects each)",
    co_medication  = paste(
      "Cyproterone acetate (injection or oral) 7 to 3 days before the first",
      "leuprorelin dose in at-risk subjects; 90% of subjects received it"
    ),
    regions        = "not reported",
    notes          = paste(
      "Study EC403: randomised, open-label, multicentre, parallel-group PK/PD",
      "study of two leuprorelin-SR 6-month depot doses (Methods 2.1 and",
      "Table 1). 62 Caucasian men randomised from 71 screened; 59 completed",
      "(28 in the 30 mg and 31 in the 22.5 mg group). EC403 had the densest",
      "sampling scheme and was used for model development; studies EC402,",
      "EC404 Stratum 2 and the Bruchovsky cohort were used only for external",
      "evaluation. Assay LOQs: leuprorelin 16 pg/mL (10.7% BLQ, M4 method),",
      "testosterone 10 ng/dL (21.9% BLQ, M3 method), PSA 0.1 ng/mL (4.5%",
      "BLQ, M3 method). The EC402 1-month and 4-month depot PK models are",
      "packaged separately as Snelder_2019_leuprorelin_1m and",
      "Snelder_2019_leuprorelin_4m."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Leuprorelin PK -- Table 2 (study EC403, 6-month depot).
    # Two-compartment disposition with two parallel release depots.
    # ------------------------------------------------------------------
    lcl       <- log(17.4)      ; label("Clearance (L/h)")                                                    # Table 2: CL = 17.4 (RSE 4.33%)
    lvc       <- log(137)       ; label("Central volume of distribution (L)")                                 # Table 2: Vc = 137 (RSE 6.15%)
    lq        <- log(3.81)      ; label("Inter-compartmental clearance (L/h)")                                # Table 2: Q = 3.81 (RSE 6.67%)
    lvp       <- log(28.1)      ; label("Peripheral volume of distribution (L)")                              # Table 2: Vp = 28.1 (RSE 11.6%)
    lka_fast  <- log(1.57)      ; label("Fast-release first-order rate constant Ka1 (1/h)")                   # Table 2: Ka1 = 1.57 (RSE 9.81%)
    lka_slow  <- log(0.000361)  ; label("Slow-release first-order rate constant Ka2, non-ramping class (1/h)") # Table 2: Ka2 = 0.000361 (RSE 9.72%)
    lka_slope <- log(0.331e-6)  ; label("Slope of the linearly time-ramping slow-release rate Ka3_SLP (1/h^2)") # Table 2: Ka3_SLP = 0.331e-6 (RSE 11.5%)

    # Table 2 reports Fr on the logit scale (footnote a:
    # 1 - exp(Fr)/(1 + exp(Fr)) = 0.55 is the SLOW fraction), so
    # expit(-0.218) = 0.4457 is the fraction released via the fast depot --
    # matching Supplement 4 $PK, `FR1 = IFR1 ; Fr = 0.45`.
    logitfrel <- -0.218         ; label("Logit of the fast-release dose fraction Fr (logit units)")           # Table 2: Fr = -0.218 (RSE 26.5%)

    # RBIO is the relative bioavailability scaling both release depots; its
    # typical value is 1 (Supplement 4 $PK, `RBIO = IRBIO ; 1`) and only its
    # IIV was estimated (Table 2 has no RBIO fixed-effect row).
    lfdepot   <- fixed(log(1))  ; label("Relative bioavailability RBIO (unitless)")                           # Supplement 4 $PK: RBIO = 1

    # ------------------------------------------------------------------
    # Testosterone model -- Table 3 and Methods 2.5, Equation 1.
    # ------------------------------------------------------------------
    # Ratio of the (unknown) endogenous GnRH concentration to its receptor
    # dissociation constant. Methods 2.5: "this ratio was combined into model
    # parameter GnRH, which was arbitrarily fixed to 1 implying 50% activity
    # of the GnRH receptors at baseline". The Discussion warns that DR50, Kd
    # and nHT are conditional on this assumption.
    agonist_kd_ratio <- fixed(1); label("Endogenous GnRH concentration relative to its receptor Kd (unitless)") # Methods 2.5: GnRH fixed to 1

    ldr50     <- log(0.0486)    ; label("Receptor-activation-fraction difference giving half-maximal downregulation DR50 (unitless)") # Table 3: DR50 = 0.0486 (RSE 22.6%)
    lhill_dr  <- log(2.05)      ; label("Hill coefficient of the receptor-downregulation sigmoid nHT (unitless)") # Table 3: nHT = 2.05 (RSE 9.90%)
    lkd       <- log(4.36)      ; label("Leuprorelin GnRH-receptor equilibrium dissociation constant Kd (pg/mL)") # Table 3: Kd = 4.36 (RSE 36.5%)
    lkout_tt  <- log(0.0164)    ; label("Testosterone first-order degradation rate constant KoutT (1/h)")     # Table 3: KoutT = 0.0164 (RSE 4.02%)
    lrbase_tt <- log(449)       ; label("Baseline serum testosterone BSLT (ng/dL)")                           # Table 3: BSLT = 449 (RSE 4.68%)

    # Apparent baseline GnRH-receptor concentration; set to 1 by convention
    # so that receptor occupancy is a fraction (Methods 2.5; Supplement 4
    # $PK, `BSLR = 1`), which also forces KinR = KoutR.
    lrbase_rt <- fixed(log(1))  ; label("Apparent baseline GnRH-receptor concentration BSLR (unitless)")      # Methods 2.5: set to 1

    lcpa      <- log(4.13)      ; label("Cyproterone acetate suppression factor: endogenous GnRH activity is divided by (1 + cpa) (unitless)") # Table 3: CPA = 4.13 (RSE 17.1%)

    # ------------------------------------------------------------------
    # PSA model -- Table 4 and Methods 2.6, Equations 4, 5 and 7 (Model 2,
    # the combined testosterone-dependent and testosterone-independent
    # production form the paper selected, dMVOF = 167 vs Model 1).
    # ------------------------------------------------------------------
    lrbase_psa <- log(9.73)     ; label("Baseline serum PSA BSLP (ng/mL)")                                    # Table 4: BSLP = 9.73 (RSE 4.68%)
    lkout_psa  <- log(0.00241)  ; label("PSA first-order dissipation rate constant KoutP (1/h)")              # Table 4: KoutP = 0.00241 (RSE 5.81%)
    lemax      <- log(57.3)     ; label("Maximum fractional stimulation of PSA production by testosterone Emax (unitless)") # Table 4: Emax = 57.3 (RSE 20.4%)
    lec50      <- log(596)      ; label("Testosterone concentration giving half-maximal PSA stimulation EC50 (ng/dL)") # Table 4: EC50 = 596 (RSE 13.9%)
    lhill_psa  <- log(2.49)     ; label("Hill coefficient of the testosterone-to-PSA sigmoid nHP (unitless)") # Table 4: nHP = 2.49 (RSE 13.8%)

    # Box-Cox shape on the PSA-baseline random effect. Supplement 4 $PK:
    # BSLP = THETA(2)*EXP(((EXP(ETA(5)))**THETA(7) - 1)/THETA(7)), THETA(7) =
    # 0.642 -- the Petersson (2009) form written out in model() below.
    boxcox_lrbase_psa <- 0.642  ; label("Box-Cox shape parameter of the PSA-baseline random effect (unitless)") # Table 4: Shape parameter Box-Cox = 0.642 (RSE 7.29%)

    # ------------------------------------------------------------------
    # Interindividual variability. Methods 2.8: "Random effects were
    # included as exponential terms reflecting lognormal distributions of
    # model parameters", so every eta below is additive on the log scale
    # except etalogitfrel, which is additive on the logit scale.
    # ------------------------------------------------------------------
    # PK etas (Table 2; reported as a diagonal, no covariances given).
    # Table 2: omega^2 RBIO = 0.127 (CV 36.8%, RSE 23.8%)
    etalfdepot   ~ 0.127
    # Table 2: omega^2 CL = 0.0304 (CV 17.6%, RSE 31.9%)
    etalcl       ~ 0.0304
    # Table 2: omega^2 Fr = 0.116 (CV 35.1%, RSE 30.7%); on the logit scale.
    etalogitfrel ~ 0.116
    # Table 2 reports a SINGLE variance for "Ka2 & Ka3_SLP", i.e. NONMEM
    # carried one eta shared by the two slow-release rate parameters (in
    # Supplement 4 they are one symbol, `SLP`, whose typical value switches
    # on the mixture class). A subject belongs to exactly one class, so only
    # one of the two is ever used; reusing this eta on both lines in model()
    # reproduces the published structure exactly.
    # Table 2: omega^2 Ka2 & Ka3_SLP = 0.232 (CV 51.1%, RSE 22.3%)
    etalka_slow  ~ 0.232

    # Testosterone-model etas -- Table 3 reports a full 6x6 block in the
    # order DR50, Kd, KoutT, BSLT, CPA, nHT (6 variances and all 15
    # covariances). Each covariance below reproduces the printed correlation
    # to within 0.0016 and the assembled matrix is positive definite.
    # Row-by-row source trace for the block below (comments must live OUTSIDE
    # the `c(...)`: rxode2 5.1.7 rewrites any comment between `c(` and `)`
    # into a bare `;` and the model then fails to parse).
    #   row 1  omega^2 DR50  = 0.696  (CV 100%,  RSE 28.4%)
    #   row 2  DR50 x Kd     = -0.249 (r -0.164); omega^2 Kd    = 3.30   (CV 511%,  RSE 23.1%)
    #   row 3  DR50 x KoutT  =  0.164 (r  0.642); Kd x KoutT  = -0.163 (r -0.294);
    #          omega^2 KoutT = 0.0933 (CV 31.3%, RSE 19.7%)
    #   row 4  DR50 x BSLT   =  0.159 (r  0.556); Kd x BSLT   = -0.0675 (r -0.108);
    #          KoutT x BSLT  = 0.0192 (r  0.183); omega^2 BSLT = 0.118 (CV 35.4%, RSE 23.1%)
    #   row 5  DR50 x CPA    =  0.626 (r  0.701); Kd x CPA    = -0.153 (r -0.0785);
    #          KoutT x CPA   =  0.179 (r  0.547); BSLT x CPA  =  0.0458 (r 0.125);
    #          omega^2 CPA   = 1.15   (CV 147%,  RSE 25.1%)
    #   row 6  DR50 x nHT    =  0.327 (r  0.912); Kd x nHT    =  0.0814 (r 0.104);
    #          KoutT x nHT   =  0.0680 (r 0.518); BSLT x nHT  =  0.0826 (r 0.560);
    #          CPA x nHT     =  0.272 (r  0.590); omega^2 nHT = 0.185 (CV 45.1%, RSE 15.6%)
    etaldr50 + etalkd + etalkout_tt + etalrbase_tt + etalcpa + etalhill_dr ~ c(
      0.696,
      -0.249,  3.30,
      0.164,  -0.163,  0.0933,
      0.159,  -0.0675, 0.0192,  0.118,
      0.626,  -0.153,  0.179,   0.0458, 1.15,
      0.327,   0.0814, 0.0680,  0.0826, 0.272,  0.185
    )

    # PSA-model etas -- Table 4 BLOCK(3) on KoutP, Emax, EC50 plus two
    # diagonal terms. Supplement 4's $OMEGA prints 0.517 for the EC50
    # variance and -0.00743 for the KoutP-EC50 covariance; Table 4's 0.0517
    # and -0.00734 are the values that reproduce the printed correlation
    # column (0.375 and -0.0834), so Table 4 is used for both.
    # Row-by-row source trace (kept outside the `c(...)`; see the note above).
    #   row 1  omega^2 KoutP = 0.150 (CV 40.2%, RSE 27.7%)
    #   row 2  KoutP x Emax  =  0.0282  (r  0.0469); omega^2 Emax = 2.41   (CV 318%, RSE 22.9%)
    #   row 3  KoutP x EC50  = -0.00734 (r -0.0834); Emax x EC50  = 0.133  (r 0.375);
    #          omega^2 EC50  = 0.0517 (CV 23.0%, RSE 27.9%)
    etalkout_psa + etalemax + etalec50 ~ c(
      0.150,
      0.0282,   2.41,
      -0.00734, 0.133,  0.0517
    )
    # Table 4: omega^2 nHP = 0.0225, reported as held fixed; Supplement 4
    # $OMEGA prints `0.0225 FIX`. (Trace kept off the declaration line: rxode2
    # rewrites a trailing comment on an eta line into label(), and both a
    # double quote and the word "Fixed" break that rewrite or the convention
    # gate.)
    etalhill_psa  ~ fixed(0.0225)
    # Table 4: omega^2 BSLPSA = 0.780 (RSE 16.4%); Box-Cox transformed in model().
    etalrbase_psa ~ 0.780

    # ------------------------------------------------------------------
    # Residual error. Every endpoint used a proportional model; the tables
    # report sigma^2 (a variance), so the SD is its square root.
    # ------------------------------------------------------------------
    propSd     <- 0.5505 ; label("Proportional residual error, leuprorelin (fraction)")  # Table 2: sigma^2 prop = 0.303 (CV 59.5%) -> sqrt(0.303) = 0.5505
    propSd_TT  <- 0.2759 ; label("Proportional residual error, testosterone (fraction)") # Table 3: sigma^2 prop = 0.0761 (CV 28.1%) -> sqrt(0.0761) = 0.2759
    propSd_PSA <- 0.2565 ; label("Proportional residual error, PSA (fraction)")          # Table 4: sigma^2 prop = 0.0658 (CV 26.1%) -> sqrt(0.0658) = 0.2565
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual PK parameters
    # ------------------------------------------------------------------
    fdepot   <- exp(lfdepot + etalfdepot)
    cl       <- exp(lcl + etalcl)
    vc       <- exp(lvc)
    q        <- exp(lq)
    vp       <- exp(lvp)
    ka_fast  <- exp(lka_fast)
    frel     <- expit(logitfrel + etalogitfrel)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    ka_slow  <- exp(lka_slow + etalka_slow)
    # Ka3_SLP carries the SAME random effect as Ka2: Table 2 reports one
    # shared variance for "Ka2 & Ka3_SLP", and Supplement 4 carries the pair
    # as a single symbol `SLP` ($PK: `SLP = ISLP ; 0.331*POP1 + 0.000361*POP2`)
    # with one eta. Multiplying by exp(lka_slope - lka_slow) transfers that
    # eta onto the ramp slope -- algebraically exp(lka_slope + etalka_slow) --
    # while leaving etalka_slow bound to exactly one typical value, which is
    # what rxode2's mu-referencing requires.
    ka_slope <- ka_slow * exp(lka_slope - lka_slow)

    # Slow-depot release rate constant. Supplement 4 $DES:
    #   K34 = SLP * TAD1 / 1E6            (ramping subpopulation)
    #   IF (IPOP .EQ. 2) K34 = SLP        (constant first-order subpopulation)
    # `max(0, tafd())` is 0 before the first dose, so pre-dose baseline
    # observation rows contribute nothing (and do not propagate NA).
    # The ramp is held in its own symbol because rxode2 5.1.7's mu-reference
    # pass fails ("mu-ref err: subscript out of bounds") on the fully inlined
    # `MIX_RAMP_REL * ka_slope * max(0, tafd()) + (1 - MIX_RAMP_REL) * ka_slow`.
    k_ramp <- ka_slope * max(0, tafd())
    k_slow <- MIX_RAMP_REL * k_ramp + (1 - MIX_RAMP_REL) * ka_slow

    # ------------------------------------------------------------------
    # 2. Individual testosterone-model parameters
    # ------------------------------------------------------------------
    dr50     <- exp(ldr50     + etaldr50)
    hill_dr  <- exp(lhill_dr  + etalhill_dr)
    kd       <- exp(lkd       + etalkd)
    kout_tt  <- exp(lkout_tt  + etalkout_tt)
    rbase_tt <- exp(lrbase_tt + etalrbase_tt)
    cpa      <- exp(lcpa      + etalcpa)
    rbase_rt <- exp(lrbase_rt)
    # Supplement 4 $PK: `KoutT = KoutR`, i.e. the receptor pool and
    # testosterone share one degradation rate constant (and one eta).
    kout_rt  <- kout_tt

    # Fraction of activated GnRH receptors at baseline (Equation 1,
    # FRAC_0 = GnRH/(1 + GnRH)); with GnRH fixed to 1 this is 0.5.
    frac0   <- agonist_kd_ratio / (1 + agonist_kd_ratio)
    # Zero-order synthesis rates from the pre-treatment steady state.
    # Receptor: KinR = BSLR * KoutR (Methods 2.5). Testosterone: Supplement 4
    # $PK, KinT = KoutT*BSLT*KoutR/(FRAC0*KinR), which with KinR = BSLR*KoutR
    # reduces to the form below.
    kin_rt  <- kout_rt * rbase_rt
    kin_tt  <- kout_tt * rbase_tt / (frac0 * rbase_rt)

    # ------------------------------------------------------------------
    # 3. Individual PSA-model parameters
    # ------------------------------------------------------------------
    kout_psa <- exp(lkout_psa + etalkout_psa)
    emax     <- exp(lemax     + etalemax)
    ec50     <- exp(lec50     + etalec50)
    hill_psa <- exp(lhill_psa + etalhill_psa)
    # Box-Cox transformed baseline random effect (Petersson 2009 form);
    # Supplement 4 $PK, BSLP = THETA(2)*EXP(((EXP(ETA(5)))**THETA(7)-1)/THETA(7)).
    phi_rbase_psa <- (exp(etalrbase_psa)^boxcox_lrbase_psa - 1) / boxcox_lrbase_psa
    rbase_psa     <- exp(lrbase_psa + phi_rbase_psa)

    # Steady-state PSA production rate, Equation 7 (Model 2): KinP is
    # expressed in terms of BSLP, KoutP and the testosterone effect
    # evaluated at the testosterone baseline.
    eff_bsl <- emax * rbase_tt^hill_psa / (ec50^hill_psa + rbase_tt^hill_psa)
    kin_psa <- kout_psa * rbase_psa / (1 + eff_bsl)

    # ------------------------------------------------------------------
    # 4. Leuprorelin PK ODEs (Figure 1; Supplement 4 $DES compartments 1-5,
    #    whose unused second dose compartment carries F2 = 0 for this study)
    # ------------------------------------------------------------------
    d/dt(depot)       <- -ka_fast * depot
    d/dt(depot2)      <- -k_slow  * depot2
    d/dt(central)     <-  ka_fast * depot + k_slow * depot2 -
      kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot)  <- fdepot * frel
    f(depot2) <- fdepot * (1 - frel)

    # central (mg) / vc (L) = mg/L; * 1e6 -> pg/mL, the assay unit.
    Cc <- central / vc * 1e6

    # ------------------------------------------------------------------
    # 5. Testosterone and GnRH-receptor ODEs -- Equation 1
    # ------------------------------------------------------------------
    # Cyproterone acetate divides the endogenous GnRH activity ratio
    # (Supplement 4 $DES, AGNN = AGN/(1 + CPA)); the covariate column gates
    # the effect to the prophylaxis window.
    agn_cpa <- agonist_kd_ratio / (1 + cpa * CONMED_CYPROTERONE)
    frac    <- (agn_cpa + Cc / kd) / (1 + agn_cpa + Cc / kd)
    rac     <- frac * RT
    # Sigmoid downregulation of receptor synthesis, guarded so that DRR = 1
    # whenever occupancy has not risen above baseline (Supplement 4 $DES:
    # `DR = 1; IF ((FRAC-FRAC0).GT.0) DR = DR50**NH/(DR50**NH+(FRAC-FRAC0)**NH)`).
    drr <- dr50^hill_dr / (dr50^hill_dr + max(0, frac - frac0)^hill_dr)

    d/dt(RT) <- kin_rt * drr - kout_rt * RT
    d/dt(TT) <- kin_tt * rac - kout_tt * TT

    RT(0) <- rbase_rt
    TT(0) <- rbase_tt

    # ------------------------------------------------------------------
    # 6. PSA turnover ODE -- Equations 4 and 5 (Model 2)
    # ------------------------------------------------------------------
    # max(0, TT) guards the fractional power against a solver undershoot
    # into slightly negative testosterone; TT is a concentration and the
    # guard is inactive for every physically meaningful value.
    eff_psa   <- emax * max(0, TT)^hill_psa / (ec50^hill_psa + max(0, TT)^hill_psa)
    d/dt(PSA) <- kin_psa * (1 + eff_psa) - kout_psa * PSA
    PSA(0)    <- rbase_psa

    # ------------------------------------------------------------------
    # 7. Observations
    # ------------------------------------------------------------------
    Cc  ~ prop(propSd)
    TT  ~ prop(propSd_TT)
    PSA ~ prop(propSd_PSA)
  })
}
