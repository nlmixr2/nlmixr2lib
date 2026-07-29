Eckernas_2023_dimethyltryptamine_LZc <- function() {
  description <- "Effect-compartment sigmoidal Emax PK/PD model for the increase in EEG signal diversity, measured as the Lempel-Ziv complexity (LZc) score, produced by intravenous N,N-dimethyltryptamine (DMT) in healthy adults (Eckernas 2023). Plasma DMT is described by a two-compartment model with first-order elimination from the central compartment whose parameters are fixed from the upstream population PK analysis (Eckernas 2022, reproduced in Table S1). The LZc score is driven by the effect-site (biophase) concentration through a sigmoidal Emax function with a maximum relative increase of about 10% over baseline. Between-subject variability is on baseline LZc, on the effect-site EC50 and on Emax; between-occasion variability is on baseline LZc. Residual variability on the LZc score is additive."
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
  units    <- list(time = "minute", dosing = "nmol", concentration = "nmol/L")

  covariateData <- list(
    OCC = list(
      description        = "Integer-valued occasion indicator. In the Eckernas 2023 fixed-sequence study design each participant contributed two occasions: OCC = 1 is the placebo visit and OCC = 2 is the DMT visit one week later.",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Constant within an occasion, time-varying within a subject. Decomposed inside model() into the binary indicators oc1 and oc2 that multiplex the two between-occasion etas on baseline LZc. Appendix S1 carries the same construct as the two data columns OCC1 and OCC2 with `$OMEGA BLOCK (1) 0.000297` followed by `$OMEGA BLOCK (1) SAME`, i.e. a single shared between-occasion variance across the two occasions.",
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
    notes          = "EEG was recorded with a 32-channel Brainproducts system (EasycapMR 32) at 1000 Hz, band-pass filtered at 1-45 Hz, averaged across channels and summarised as mean values per minute for modelling. Spontaneous signal diversity was computed as a Lempel-Ziv complexity score (see Timmermann 2019, Sci Rep 9:16324, for the EEG analysis details). EEG recordings from one participant who received 20 mg were excluded because of excessive movement artifacts, so the PD models were fitted to 12 participants. The authors cannot be certain that the true maximum response was reached at the doses studied (Eckernas 2023 Discussion)."
  )

  ini({
    # ---- Plasma PK: fixed from the upstream population PK model --------------
    # Eckernas 2023 used the 'population PK parameters and data' approach of
    # Zhang 2003: the population PK parameters are FIXED and only the individual
    # PK parameters are estimated simultaneously with the PD parameters.
    lcl <- fixed(log(26.0));  label("Clearance CL (L/min)")                                # Table S1: CL = 26.0 (20.6; 33.6) L/min, %RSE 15.1; Appendix S1 LZc model `26 FIX ; CL`
    lvc <- fixed(log(221));   label("Central volume of distribution Vc (L)")               # Table S1: Vc = 221 (181-273) L, %RSE 12.1; Appendix S1 `221 FIX ; V1`
    lq  <- fixed(log(2.99));  label("Inter-compartmental clearance Q (L/min)")             # Table S1: Q = 2.99 (1.87; 5.44) L/min, %RSE 36.6; Appendix S1 `2.99 FIX ; Q`
    lvp <- fixed(log(59.0));  label("Peripheral volume of distribution Vp (L)")            # Table S1: Vp = 59.0 (48.0; 82.7) L, %RSE 18.5; Appendix S1 `59 FIX ; V2`

    # ---- PD: Lempel-Ziv complexity score (Eckernas 2023 Table 3) ------------
    lrbase <- log(0.321);     label("Baseline Lempel-Ziv complexity score R0 (unitless)")  # Table 3: R0 = 0.321 (0.315; 0.332), %RSE 1.6; Appendix S1 THETA(5) = 0.321
    lemax  <- log(0.0996);    label("Maximum fractional increase in LZc score Emax (unitless)")  # Table 3: Emax = 0.10 (0.091; 0.11), %RSE 4.8; Appendix S1 THETA(7) = 0.0996
    lec50  <- log(54.1);      label("Effect-site concentration for half-maximal response EC50,e (nmol/L)")  # Table 3: EC50,e = 54 (38; 72) nM, %RSE 19; Appendix S1 THETA(8) = 54.1
    lke0   <- log(0.762);     label("Effect-compartment equilibration rate constant ke0 (1/min)")  # Table 3: Ke0 = 0.76 (0.65; 0.96) 1/min, %RSE 12; Appendix S1 THETA(6) = 0.762
    lhill  <- log(4.79);      label("Hill coefficient gamma (unitless)")                   # Table 3: gamma = 4.8 (3.9; 5.9), %RSE 13; Appendix S1 THETA(9) = 4.79

    # ---- Random effects -----------------------------------------------------
    # Eckernas 2023 reports the variability magnitudes as %CV computed as
    # sqrt(omega^2) * 100 (verifiable against every Appendix S1 $OMEGA value),
    # not as sqrt(exp(omega^2) - 1). The variances below are therefore taken
    # directly from the Appendix S1 $OMEGA blocks. A 92% correlation between the
    # between-subject etas on R0 and Emax was observed but was NOT estimated in
    # the final model because it caused poor precision and ill conditioning
    # (Eckernas 2023 Results), so the etas are kept diagonal here.
    etalcl ~ fixed(0.224)              # Appendix S1 `$OMEGA 0.224 FIX ; IIV CL`; Table S1 BSV CL = 47.3 %CV = sqrt(0.224)
    etalrbase ~ 0.00266                # Appendix S1 `0.00266 ; IIV E0`;   Table 3 BSV R0      = 5.2 %CV = sqrt(0.00266)
    etalec50  ~ 0.599                  # Appendix S1 `0.599 ; IIV EC50`;   Table 3 BSV EC50,e  = 77  %CV = sqrt(0.599)
    etalemax  ~ 0.176                  # Appendix S1 `0.176 ; IIV Emax`;   Table 3 BSV Emax    = 42  %CV = sqrt(0.176)
    etaiov_rbase_1 ~ 0.000297          # Appendix S1 `$OMEGA BLOCK (1) 0.000297`; Table 3 BOV R0 = 1.7 %CV = sqrt(0.000297)
    etaiov_rbase_2 ~ fixed(0.000297)   # Appendix S1 `$OMEGA BLOCK (1) SAME` - the same variance is shared across occasions

    # ---- Residual error -----------------------------------------------------
    propSd    <- fixed(0.499);  label("Proportional residual error on plasma DMT (fraction)")     # Appendix S1 `0.499 FIX ; PK error`; Table S1 residual error DMT = 50.0 (43.9; 58.9) %CV
    addSd_lzc <- 0.0061;        label("Additive residual error on the LZc score (LZc units)")     # Table 3: additive error (SD) = 0.0061 (0.0058; 0.0066), %RSE 4.1; Appendix S1 THETA(11) = 0.0061
  })

  model({
    # Occasion indicators driving the shared between-occasion variance on the
    # baseline response (Appendix S1 `ETA(5) * OCC1 + ETA(6) * OCC2`).
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
    emax  <- exp(lemax + etalemax)
    ec50  <- exp(lec50 + etalec50)
    ke0   <- exp(lke0)
    hill  <- exp(lhill)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment IV disposition (Appendix S1 $DES). Doses are amounts in
    # nmol so that central / vc is in nmol/L (= nM), the unit in which EC50,e is
    # reported.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    Cc <- central / vc

    # Effect compartment (biophase). Appendix S1 carries the effect-site
    # concentration as a state: DADT(3) = K30 * (CONC - A(3)).
    d/dt(effect) <- ke0 * (Cc - effect)

    # Sigmoidal Emax response; emax is a fractional increase over baseline, so
    # the LZc score rises from rbase towards rbase * (1 + emax).
    lzc <- rbase * (1 + emax * effect^hill / (ec50^hill + effect^hill))

    Cc  ~ prop(propSd)
    lzc ~ add(addSd_lzc)
  })
}
