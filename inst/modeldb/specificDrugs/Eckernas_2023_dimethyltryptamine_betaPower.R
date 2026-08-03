Eckernas_2023_dimethyltryptamine_betaPower <- function() {
  description <- "Effect-compartment sigmoidal Imax PK/PD model for the suppression of EEG beta power by intravenous N,N-dimethyltryptamine (DMT) in healthy adults (Eckernas 2023). Plasma DMT is described by a two-compartment model with first-order elimination from the central compartment whose parameters are fixed from the upstream population PK analysis (Eckernas 2022, reproduced in Table S1). Beta power is driven by the effect-site (biophase) concentration through a sigmoidal Imax function; unlike alpha power, suppression is only partial, with Imax estimated at 0.70. Between-subject variability is on baseline beta power and on the effect-site IC50 as a correlated block; between-occasion variability is on baseline beta power. Residual variability on beta power is proportional."
  reference <- paste(
    "Eckernas E, Timmermann C, Carhart-Harris R, Roshammar D, Ashton M.",
    "N,N-dimethyltryptamine affects electroencephalography response in a",
    "concentration-dependent manner - A pharmacokinetic/pharmacodynamic analysis.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(4):474-486. doi:10.1002/psp4.12933.",
    "The plasma PK parameters are fixed from Eckernas E, Timmermann C,",
    "Carhart-Harris R, Roshammar D, Ashton M. Population",
    "pharmacokinetic/pharmacodynamic modelling of the psychedelic experience",
    "induced by N,N-dimethyltryptamine - implications for dose considerations.",
    "Clin Transl Sci. 2022;15(12):2928-2937; those fixed values are reproduced in",
    "Table S1 of the 2023 paper and in its Appendix S1 NONMEM control streams.",
    sep = " "
  )
  vignette <- "Eckernas_2023_dimethyltryptamine_EEG"
  units    <- list(time = "min", dosing = "nmol", concentration = "nmol/L")

  covariateData <- list(
    OCC = list(
      description        = "Integer-valued occasion indicator. In the Eckernas 2023 fixed-sequence study design each participant contributed two occasions: OCC = 1 is the placebo visit and OCC = 2 is the DMT visit one week later.",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Constant within an occasion, time-varying within a subject. Decomposed inside model() into the binary indicators oc1 and oc2 that multiplex the two between-occasion etas on baseline beta power. Appendix S1 carries the same construct as the two data columns OCC1 and OCC2 with `$OMEGA BLOCK(1) 0.0436 ; IOV I0` followed by `$OMEGA BLOCK (1) SAME`, i.e. a single shared between-occasion variance across the two occasions.",
      source_name        = "OCC1 / OCC2"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 13L,                                                       # Eckernas 2023 Methods (Clinical study): 13 healthy subjects contributed plasma PK
    n_studies      = 1L,                                                        # single placebo-controlled, single-blind, fixed-sequence pilot study
    age_range      = "22-48 years",                                             # Eckernas 2023 Methods (Clinical study)
    age_median     = "33 years",                                                # Eckernas 2023 Methods (Clinical study)
    sex_female_pct = 46.2,                                                      # Eckernas 2023 Methods: seven men of 13 subjects; (13 - 7) / 13 = 46.2%
    disease_state  = "healthy volunteers",                                      # Eckernas 2023 Methods (Clinical study)
    dose_range     = "Placebo at visit 1 and a single intravenous bolus of 7 mg (n = 3), 14 mg (n = 4), 18 mg (n = 1) or 20 mg (n = 5) DMT fumarate at visit 2, one week later.",  # Eckernas 2023 Methods (Clinical study)
    regions        = "United Kingdom (National Institute of Health Research Imperial Clinical Research Facility, London).",  # Eckernas 2023 Methods (Clinical study)
    n_observations = "252 EEG observations after DMT and 238 after placebo were available across alpha power, beta power, delta power, theta power and LZc score (84, 63, 21 and 84 post-DMT observations for the 7, 14, 18 and 20 mg dose levels); the fixed PK model was built on 93 DMT plasma concentrations.",  # Eckernas 2023 Results, first paragraph
    notes          = "EEG was recorded with a 32-channel Brainproducts system (EasycapMR 32) at 1000 Hz, band-pass filtered at 1-45 Hz, averaged across channels and summarised as mean values per minute for modelling; the beta band is 13-30 Hz. EEG recordings from one participant who received 20 mg were excluded because of excessive movement artifacts, so the PD models were fitted to 12 participants. Only the highest dose produced a clear effect on beta power and Imax was not reached in the study, so the beta-power parameter estimates should be interpreted with caution (Eckernas 2023 Discussion)."
  )

  ini({
    # ---- Plasma PK: fixed from the upstream population PK model --------------
    # Eckernas 2023 used the 'population PK parameters and data' approach of
    # Zhang 2003: the population PK parameters are FIXED and only the individual
    # PK parameters are estimated simultaneously with the PD parameters.
    lcl <- fixed(log(26.0));  label("Clearance CL (L/min)")                                # Table S1: CL = 26.0 (20.6; 33.6) L/min, %RSE 15.1; Appendix S1 beta model `26 FIX ; CL`
    lvc <- fixed(log(221));   label("Central volume of distribution Vc (L)")               # Table S1: Vc = 221 (181-273) L, %RSE 12.1; Appendix S1 `221 FIX ; V1`
    lq  <- fixed(log(2.99));  label("Inter-compartmental clearance Q (L/min)")             # Table S1: Q = 2.99 (1.87; 5.44) L/min, %RSE 36.6; Appendix S1 `2.99 FIX ; Q`
    lvp <- fixed(log(59.0));  label("Peripheral volume of distribution Vp (L)")            # Table S1: Vp = 59.0 (48.0; 82.7) L, %RSE 18.5; Appendix S1 `59 FIX ; V2`

    # ---- PD: beta power (Eckernas 2023 Table 2) -----------------------------
    lrbase <- log(0.0644);    label("Baseline beta power R0 (EEG beta power units)")       # Table 2: R0 = 0.064 (0.052; 0.079), %RSE 13; Appendix S1 THETA(5) = 0.0644
    limax  <- log(0.704);     label("Maximum fractional suppression of beta power Imax (unitless)")  # Table 2: Imax = 0.70 (0.66; 0.72), %RSE 2.4; Appendix S1 THETA(7) = 0.704
    lec50  <- log(137);       label("Effect-site concentration for half-maximal suppression IC50,e (nmol/L)")  # Table 2: IC50,e = 137 (104; 186) nM, %RSE 18; Appendix S1 THETA(8) = 137
    lke0   <- log(1.19);      label("Effect-compartment equilibration rate constant ke0 (1/min)")  # Table 2: Ke0 = 1.2 (0.95; 1.7) 1/min, %RSE 19; Appendix S1 THETA(6) = 1.19
    lhill  <- log(5.19);      label("Hill coefficient gamma (unitless)")                   # Table 2: gamma = 5.2 (4.1; 6.7), %RSE 15; Appendix S1 THETA(9) = 5.19

    # ---- Random effects -----------------------------------------------------
    # Eckernas 2023 reports the variability magnitudes as %CV computed as
    # sqrt(omega^2) * 100 (verifiable against every Appendix S1 $OMEGA value),
    # not as sqrt(exp(omega^2) - 1). The variances below are therefore taken
    # directly from the Appendix S1 $OMEGA blocks.
    etalcl ~ fixed(0.224)              # Appendix S1 `$OMEGA 0.224 FIX ; IIV CL`; Table S1 BSV CL = 47.3 %CV = sqrt(0.224)

    # Correlated block on baseline and IC50,e. Variances from Appendix S1
    # `$OMEGA BLOCK (2) 0.402 ; IIV I0 / 0.21 0.556 ; IIV IC50`, which reproduce
    # Table 2's BSV R0 = 63 %CV (sqrt(0.402)) and BSV IC50,e = 75 %CV
    # (sqrt(0.556)). The covariance is set from Table 2's reported correlation of
    # -46% (95% CI -57; -33): -0.46 * sqrt(0.402 * 0.556) = -0.2175. Appendix S1
    # prints the off-diagonal as `0.21` without a minus sign, which contradicts
    # both Table 2 (CI excludes zero) and the Discussion ("higher baseline values
    # were associated with lower IC50,e values"); the published sign is used.
    etalrbase + etalec50 ~ c(0.402,
                             -0.2175, 0.556)

    etaiov_rbase_1 ~ 0.0436            # Appendix S1 `$OMEGA BLOCK(1) 0.0436 ; IOV I0`; Table 2 BOV R0 = 21 %CV = sqrt(0.0436)
    etaiov_rbase_2 ~ fixed(0.0436)     # Appendix S1 `$OMEGA BLOCK (1) SAME` - the same variance is shared across occasions

    # ---- Residual error -----------------------------------------------------
    # A combined proportional + additive error improved the fit (dOFV = -21) but
    # the additive term was small and degraded parameter precision, so the final
    # model uses a proportional error only (Eckernas 2023 Results).
    propSd           <- fixed(0.499);  label("Proportional residual error on plasma DMT (fraction)")  # Appendix S1 `0.499 FIX ; PK error`; Table S1 residual error DMT = 50.0 (43.9; 58.9) %CV
    propSd_betaPower <- 0.175;         label("Proportional residual error on beta power (fraction)")  # Table 2: proportional error = 18 (17; 19) %CV, %RSE 4.3; Appendix S1 THETA(11) = 0.175
  })

  model({
    # Occasion indicators driving the shared between-occasion variance on the
    # baseline response (Appendix S1 `ETA(4) * OCC1 + ETA(5) * OCC2`).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_rbase <- oc1 * etaiov_rbase_1 + oc2 * etaiov_rbase_2

    # Individual PK parameters. Only CL carries (fixed) between-subject
    # variability; no variability in volume could be estimated from the
    # available PK data (Eckernas 2023 Discussion).
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    # Individual PD parameters.
    rbase <- exp(lrbase + etalrbase + iov_rbase)
    imax  <- exp(limax)
    ec50  <- exp(lec50 + etalec50)
    ke0   <- exp(lke0)
    hill  <- exp(lhill)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment IV disposition (Appendix S1 $DES). Doses are amounts in
    # nmol so that central / vc is in nmol/L (= nM), the unit in which IC50,e is
    # reported.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    Cc <- central / vc

    # Effect compartment (biophase). Appendix S1 carries the effect-site
    # concentration as a state: DADT(3) = K30 * (CONC - A(3)).
    d/dt(effect) <- ke0 * (Cc - effect)

    # Sigmoidal Imax response; imax = 0.70 means beta power plateaus at 30% of
    # baseline rather than being fully suppressed.
    betaPower <- rbase * (1 - imax * effect^hill / (ec50^hill + effect^hill))

    Cc        ~ prop(propSd)
    betaPower ~ prop(propSd_betaPower)
  })
}
