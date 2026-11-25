Astrand_2019_ticagrelor <- function() {
  description <- paste(
    "Joint population PK/PD model for oral ticagrelor, its active metabolite AR-C124910XX",
    "(TAM) and platelet reactivity (VerifyNow P2Y12 reaction units, PRU) in patients with",
    "stable coronary artery disease (ONSET/OFFSET) and with a prior myocardial infarction",
    "(PEGASUS-TIMI 54 platelet function substudy) (Astrand 2019). Ticagrelor: Erlang-type",
    "absorption through eight transit compartments sharing one rate constant KTR, with an",
    "absorption lag time in the prior-MI population only, into a two-compartment disposition.",
    "A fixed fraction Fm = 0.22 of ticagrelor clearance forms TAM, which has its own",
    "two-compartment disposition. PRU falls from a study-specific baseline through a direct",
    "sigmoid Emax function of the ticagrelor plasma concentration (logit-normal Emax).",
    sep = " "
  )
  reference <- paste(
    "Astrand M, Amilon C, Roshammar D, Himmelmann A, Angiolillo DJ, Storey RF, Gurbel PA,",
    "Bonaca MP, Hamren B. Pharmacokinetic-pharmacodynamic modelling of platelet response to ticagrelor",
    "in stable coronary artery disease and prior myocardial infarction patients.",
    "Br J Clin Pharmacol. 2019;85(2):413-421. doi:10.1111/bcp.13812",
    sep = " "
  )
  vignette <- "Astrand_2019_ticagrelor"
  units <- list(time = "h", dosing = "mg", concentration = "nmol/L")

  # Every state carries ticagrelor mass (mg). The metabolite states hold
  # ticagrelor-equivalent mass: Fm is a molar fraction of the metabolised
  # parent, so converting both analytes with the ticagrelor molecular weight
  # gives each concentration in nmol/L without needing the metabolite's own
  # molecular weight.
  compartmentData <- list(
    depot = list(analyte = "ticagrelor", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "ticagrelor", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "ticagrelor", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "ticagrelor", units = "mg", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "ticagrelor", units = "mg", specimen = "administration site", verified = TRUE),
    transit5 = list(analyte = "ticagrelor", units = "mg", specimen = "administration site", verified = TRUE),
    transit6 = list(analyte = "ticagrelor", units = "mg", specimen = "administration site", verified = TRUE),
    transit7 = list(analyte = "ticagrelor", units = "mg", specimen = "administration site", verified = TRUE),
    transit8 = list(analyte = "ticagrelor", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ticagrelor", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ticagrelor", units = "mg", specimen = "tissue", verified = TRUE),
    central_tam = list(
      analyte = "AR-C124910XX (TAM)",
      units = "mg ticagrelor-equivalent",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1_tam = list(
      analyte = "AR-C124910XX (TAM)",
      units = "mg ticagrelor-equivalent",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list(
    STUDY_PEGASUS = list(
      description = "1 = patient from the PEGASUS-TIMI 54 platelet function substudy (prior myocardial infarction 1-3 years before enrolment), 0 = patient from the ONSET/OFFSET study (stable coronary artery disease)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Study and patient population coincide in this pooled analysis, so the paper names the",
        "same split both ways: the absorption lag time is reported for 'prior MI' patients and the",
        "PRU baseline for 'PEGASUS' (Astrand 2019 Table 2). The indicator turns on the 0.48 h",
        "absorption lag and switches the PRU baseline from 311 to 261. No other parameter differed",
        "between the populations (CL/F, EC50, Emax and the Hill coefficient were tested).",
        sep = " "
      ),
      source_name = "study (ONSET/OFFSET vs PEGASUS-TIMI 54)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 303L,
    n_studies = 2L,
    age_range = "not reported in the source",
    weight_range = "not reported in the source",
    sex_female_pct = NA_real_,
    race_ethnicity = "Mainly white: 88% in ONSET/OFFSET and 93% in the PEGASUS-TIMI 54 platelet function substudy",
    disease_state = paste(
      "ONSET/OFFSET: 123 patients with stable coronary artery disease randomised to ticagrelor",
      "(n = 57), clopidogrel (n = 54; clopidogrel data not modelled) or placebo (n = 12).",
      "PEGASUS-TIMI 54 platelet function substudy: 180 patients with a myocardial infarction",
      "1-3 years before enrolment, randomised to placebo (n = 64), ticagrelor 60 mg (n = 58) or",
      "ticagrelor 90 mg twice daily (n = 58). All patients received aspirin.",
      sep = " "
    ),
    dose_range = paste(
      "ONSET/OFFSET: ticagrelor 180 mg oral loading dose then 90 mg twice daily for about",
      "6 weeks. PEGASUS-TIMI 54: ticagrelor 60 or 90 mg twice daily for at least 4 months.",
      sep = " "
    ),
    regions = "ONSET/OFFSET: not stated in the source; PEGASUS-TIMI 54 substudy: three UK centres (Sheffield, Rotherham, Nottingham) and one US centre (Jacksonville, Florida)",
    notes = paste(
      "871 ticagrelor and 872 AR-C124910XX plasma concentrations and 954 PRU measurements",
      "(including 230 from placebo patients) (Astrand 2019 Results; Table 1). NONMEM 7.3.0,",
      "FOCE-I. The PK model was fitted first; the PK/PD model fixed the population PK",
      "parameters and re-estimated the individual PK parameters jointly with the PD data",
      "(PPP&D). Below-LOQ PK samples after the first dose were handled with the M3 method.",
      sep = " "
    )
  )

  ini({
    # ---- Ticagrelor PK (Astrand 2019 Table 2) ----
    lcl <- log(16.6); label("Ticagrelor apparent clearance CL/F (L/h)") # Table 2: CL/F = 16.6 L/h (RSE 3.6%)
    lq <- log(10.4); label("Ticagrelor apparent intercompartmental clearance Q/F (L/h)") # Table 2: Q/F = 10.4 L/h (RSE 4.7%)
    lvc <- log(156); label("Ticagrelor apparent central volume Vc/F (L)") # Table 2: Vc/F = 156 L (RSE 4.2%)
    lvp <- log(55.8); label("Ticagrelor apparent peripheral volume Vp/F (L)") # Table 2: Vp/F = 55.8 L (RSE 8.3%)
    lktr <- log(10.1); label("Erlang absorption transit rate constant KTR (1/h)") # Table 2: KTR = 10.1 1/h (RSE 5.8%)
    ltlag <- log(0.48); label("Absorption lag time, prior-MI (PEGASUS) patients only (h)") # Table 2: 'Absorption lag time prior MI' = 0.48 h (RSE 11%)
    lfdepot <- fixed(log(1)); label("Relative oral bioavailability Frel (unitless)") # Table 2: Frel = 1 FIX

    # ---- AR-C124910XX (TAM) PK (Astrand 2019 Table 2) ----
    lcl_tam <- log(10.2); label("AR-C124910XX apparent clearance CLm/F (L/h)") # Table 2: CLm/F = 10.2 L/h (RSE 3.0%)
    fm <- fixed(0.22); label("Fraction of ticagrelor clearance forming AR-C124910XX, Fm (unitless)") # Table 2: Fm = 0.22 FIX; Methods: 'was fixed to 0.22'
    lq_tam <- log(4.41); label("AR-C124910XX apparent intercompartmental clearance Qm/F (L/h)") # Table 2: Qm/F = 4.41 L/h (RSE 4.7%)
    lvc_tam <- log(7.04); label("AR-C124910XX apparent central volume Vcm/F (L)") # Table 2: Vcm/F = 7.04 L (RSE 6.7%)
    lvp_tam <- log(42.3); label("AR-C124910XX apparent peripheral volume Vpm/F (L)") # Table 2: Vpm/F = 42.3 L (RSE 6.1%)

    # ---- PRU pharmacodynamics (Astrand 2019 Table 2 and the Results equation) ----
    rbase <- 311; label("PRU baseline, ONSET/OFFSET stable CAD patients (PRU)") # Table 2: PRU baseline ONSET/OFFSET = 311 (RSE 1.9%)
    rbase_pegasus <- 261; label("PRU baseline, PEGASUS-TIMI 54 prior-MI patients (PRU)") # Table 2: PRU baseline PEGASUS = 261 (RSE 1.7%)
    lec50 <- log(116); label("Ticagrelor concentration at half-maximal PRU inhibition EC50 (nmol/L)") # Table 2: EC50 = 116 nmol/L (RSE 5.3%)
    logitemax <- qlogis(0.985); label("Logit of the maximal fractional PRU inhibition Emax (logit units)") # Table 2: Emax = 98.5% (RSE 0.1%), back-transformed from the estimated logit(Emax)
    hill <- 1.59; label("Steepness of the exposure-response curve, gamma (unitless)") # Table 2: gamma = 1.59 (RSE 3.4%)

    # ---- Between-subject variability ----
    # Table 2 'BSV (%)' column. For the exponential (log-normal) etas the
    # percentage is read as sqrt(omega^2) * 100, so omega^2 = (BSV / 100)^2;
    # the footnote states the PRU-baseline and Emax entries are standard
    # deviations (additive on PRU, and on the logit scale). The vignette
    # checks this reading against the inhibition quantiles printed in the
    # Results.
    etalcl ~ 0.0576 # Table 2: BSV CL/F 24% -> 0.24^2
    etalq ~ 0.9025 # Table 2: BSV Q/F 95% -> 0.95^2
    etalktr ~ 0.2916 # Table 2: BSV KTR 54% -> 0.54^2
    etalfdepot ~ 0.1024 # Table 2: BSV Frel 32% -> 0.32^2
    etalcl_tam ~ 0.0784 # Table 2: BSV CLm/F 28% -> 0.28^2
    etalvp_tam ~ 0.1369 # Table 2: BSV Vpm/F 37% -> 0.37^2
    # Table 2: PRU baseline SD 35 (additive), EC50 BSV 66%, correlation 0.33,
    # so cov = 0.33 * 35 * 0.66 = 7.623.
    etarbase + etalec50 ~ c(1225, 7.623, 0.4356)
    etalogitemax ~ 4.84 # Table 2: Emax BSV 2.2 (SD on the logit scale) -> 2.2^2

    # ---- Residual error ----
    propSd <- 0.32; label("Proportional residual error, ticagrelor (fraction)") # Table 2: proportional residual error ticagrelor = 32% (RSE 4.5%)
    propSd_tam <- 0.26; label("Proportional residual error, AR-C124910XX (fraction)") # Table 2: proportional residual error metabolite = 26% (RSE 4.9%)
    # Table 2 footnote: PRU residual SD = (PRU/300)^alpha * SD_300, which is
    # rxode2's power error SD = propSd_PRU * PRU^alpha with propSd_PRU = SD_300 / 300^alpha.
    propSd_PRU <- 47.6 / 300^0.48; label("PRU power-error coefficient, SD_300 / 300^alpha (PRU^(1 - alpha))") # Table 2: additive residual error PRU at PRU = 300 is 47.6 (RSE 5.6%)
    powExp_PRU <- 0.48; label("Exponent alpha of the PRU residual error (unitless)") # Table 2: exponent for PRU error alpha = 0.48 (RSE 2.8%)
  })
  model({
    # Ticagrelor molecular weight (g/mol), used to convert mg doses to the
    # nmol/L concentrations the paper reports. The paper's own steady-state
    # averages pin it: 90 mg / (16.6 L/h * 12 h) / 522.57 g/mol = 864 nmol/L
    # and 60 mg gives 576 nmol/L (Astrand 2019 Results).
    mw_ticagrelor <- 522.57

    # ---- Individual parameters ----
    cl <- exp(lcl + etalcl)
    q <- exp(lq + etalq)
    vc <- exp(lvc)
    vp <- exp(lvp)
    ktr <- exp(lktr + etalktr)
    tlag <- exp(ltlag) * STUDY_PEGASUS

    cl_tam <- exp(lcl_tam + etalcl_tam)
    q_tam <- exp(lq_tam)
    vc_tam <- exp(lvc_tam)
    vp_tam <- exp(lvp_tam + etalvp_tam)

    rbase_pop <- rbase * (1 - STUDY_PEGASUS) + rbase_pegasus * STUDY_PEGASUS
    rbase_i <- rbase_pop + etarbase
    ec50 <- exp(lec50 + etalec50)
    logitemax_i <- logitemax + etalogitemax
    emax <- expit(logitemax_i)

    # ---- Micro-constants ----
    # Ticagrelor leaves central at CL/F; the fraction fm of that flux forms
    # AR-C124910XX (Methods: 'The fraction of ticagrelor metabolized to
    # AR-C124910XX (Fm) was fixed to 0.22').
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    kel_tam <- cl_tam / vc_tam
    k12_tam <- q_tam / vc_tam
    k21_tam <- q_tam / vp_tam

    # ---- ODEs ----
    # Erlang absorption: the same KTR moves drug from the dose compartment to
    # the first transit compartment, between the eight transit compartments,
    # and from the last one into central (Results, 'Ticagrelor population PK
    # model').
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(transit4) <- ktr * transit3 - ktr * transit4
    d/dt(transit5) <- ktr * transit4 - ktr * transit5
    d/dt(transit6) <- ktr * transit5 - ktr * transit6
    d/dt(transit7) <- ktr * transit6 - ktr * transit7
    d/dt(transit8) <- ktr * transit7 - ktr * transit8
    d/dt(central) <- ktr * transit8 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(central_tam) <- fm * kel * central - kel_tam * central_tam -
      k12_tam * central_tam + k21_tam * peripheral1_tam
    d/dt(peripheral1_tam) <- k12_tam * central_tam - k21_tam * peripheral1_tam

    f(depot) <- exp(lfdepot + etalfdepot)
    alag(depot) <- tlag

    # ---- Observations ----
    # mg/L * 1e6 / (g/mol) = nmol/L.
    Cc <- central / vc * 1e6 / mw_ticagrelor
    Cc_tam <- central_tam / vc_tam * 1e6 / mw_ticagrelor

    # Direct sigmoid Emax model on the ticagrelor concentration (Results
    # equation): PRU_t = PRU_baseline * (1 - Emax * C^gamma / (EC50^gamma + C^gamma)).
    inhib <- emax * Cc^hill / (ec50^hill + Cc^hill)
    PRU <- rbase_i * (1 - inhib)

    Cc ~ prop(propSd)
    Cc_tam ~ prop(propSd_tam)
    PRU ~ pow(propSd_PRU, powExp_PRU)
  })
}
