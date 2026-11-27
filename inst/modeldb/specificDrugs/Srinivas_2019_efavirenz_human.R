Srinivas_2019_efavirenz_human <- function() {
  description <- paste(
    "Eight-compartment PK model of oral efavirenz in plasma, cerebrospinal",
    "fluid (CSF) and brain tissue of adults living with HIV (Srinivas 2019,",
    "clinical model). Dose enters a depot and passes through three transit",
    "compartments sharing one rate constant (TAU) before first-order absorption",
    "(Ka) into a two-compartment plasma disposition model with linear",
    "elimination; CSF and brain tissue exchange bidirectionally with plasma and",
    "with each other through first-order rate constants on amounts. The",
    "structure was developed in rhesus macaques (see",
    "Srinivas_2019_efavirenz_monkey) and re-estimated on 29 paired plasma/CSF",
    "samples from 24 participants by iterative two-stage (IT2S) Bayesian",
    "estimation in ADAPT 5 with VCSF fixed at 0.15 L; brain tissue was not",
    "sampled in humans, so the brain parameters are model predictions anchored",
    "on the macaque model. Typical values are the geometric means of the",
    "individual estimates (Table 1).",
    sep = " "
  )
  reference <- paste(
    "Srinivas N, Joseph SB, Robertson K, Kincer LP, Menezes P, Adamson L,",
    "Schauer AP, Blake KH, White N, Sykes C, Luciw P, Eron JJ, Forrest A,",
    "Price RW, Spudich S, Swanstrom R, Kashuba ADM. Predicting Efavirenz",
    "Concentrations in the Brain Tissue of HIV-Infected Individuals and",
    "Exploring their Relationship to Neurocognitive Impairment.",
    "Clin Transl Sci. 2019;12(3):302-311. doi:10.1111/cts.12620.",
    sep = " "
  )
  vignette <- "Srinivas_2019_efavirenz"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "efavirenz", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "efavirenz", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "efavirenz", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "efavirenz", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "efavirenz", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "efavirenz", units = "mg", specimen = "plasma", verified = TRUE),
    csf = list(analyte = "efavirenz", units = "mg", specimen = "CSF", verified = TRUE),
    brain = list(analyte = "efavirenz", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 24,
    n_studies = 1,
    age_range = "19-66 years (cohort A median 32, cohort B median 54)",
    weight_range = "45.4-136.1 kg (cohort A median 65.8, cohort B median 71.0)",
    sex_female_pct = 20.8,
    race_ethnicity = c(White = 54.2, Black = 45.8),
    disease_state = "HIV-1 infection on efavirenz-based antiretroviral therapy (3 treatment-naive initiating therapy, 21 treatment-experienced and virally suppressed)",
    dose_range = "efavirenz 600 mg orally once daily (Atripla in 19 of 24)",
    regions = "United States (THINC cohort, PO1 MH094177)",
    notes = paste(
      "Srinivas 2019 Table 2 (cohort A n = 3, cohort B n = 21) and Results.",
      "29 paired plasma and CSF concentrations (one paired sample per visit,",
      "no preselected sampling time) from 24 participants; 2 non-adherent",
      "participants with BLOQ plasma and CSF were excluded. Female 0/3 + 5/21;",
      "race White 0/3 + 13/21, African American 3/3 + 8/21.",
      sep = " "
    )
  )

  ini({
    # Srinivas 2019 Table 1, 'Humans' rows: geometric mean of the individual
    # IT2S estimates (typical value) and CV% (between-subject variability).
    # Parameters are apparent (/F). Rate constants act on amounts (1/h).
    lcl <- log(12.37); label("Apparent clearance CL/F (L/h)") # Table 1 Humans CL/F geometric mean 12.37 L/h (CV 58.5%)
    lvc <- log(155.34); label("Apparent central volume V5/F (L)") # Table 1 Humans V5/F geometric mean 155.34 L (CV 39.1%)
    lka <- log(0.135); label("Absorption rate constant from transit3 to central Ka (1/h)") # Table 1 Humans Ka geometric mean 0.135 1/h (CV 94.6%)
    lvp <- log(405.4); label("Apparent peripheral volume V6/F (L)") # Table 1 Humans V6/F geometric mean 405.4 L (CV 35.8%)
    lq <- log(6.31); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 1 Humans Q/F geometric mean 6.31 L/h (CV 142.0%)
    lktr <- log(1.00); label("Transit rate constant TAU, depot -> transit1 -> transit2 -> transit3 (1/h)") # Table 1 Humans TAU geometric mean 1.00 1/h (CV 27.8%)
    lk_central_csf <- log(4.27e-6); label("Plasma-to-CSF rate constant Kplasma-CSF (1/h)") # Table 1 Humans Kplasma-CSF geometric mean 4.27E-06 1/h (CV 88.8%)
    lk_csf_central <- log(0.244); label("CSF-to-plasma rate constant KCSF-plasma (1/h)") # Table 1 Humans KCSF-plasma geometric mean 0.244 1/h (CV 26.0%)
    lk_central_brain <- log(2.00e-3); label("Plasma-to-brain-tissue rate constant Kplasma-brain (1/h)") # Table 1 Humans Kplasma-brain geometric mean 2.00E-03 1/h (CV 30.2%)
    lk_brain_central <- log(0.175); label("Brain-tissue-to-plasma rate constant Kbrain-plasma (1/h)") # Table 1 Humans Kbrain-plasma geometric mean 0.175 1/h (CV 31.2%)
    lk_csf_brain <- log(3.00e-3); label("CSF-to-brain-tissue rate constant KCSF-brain (1/h)") # Table 1 Humans KCSF-brain geometric mean 3.00E-03 1/h (CV 40.2%)
    lk_brain_csf <- log(2.41e-6); label("Brain-tissue-to-CSF rate constant Kbrain-CSF (1/h)") # Table 1 Humans Kbrain-CSF geometric mean 2.41E-06 1/h (CV 70.5%)
    lv_brain <- log(0.514); label("Apparent brain tissue volume Vbrain/F, conditioned on VCSF (L)") # Table 1 Humans Vbrain/F geometric mean 0.514 L (CV 25.1%)
    lvcsf <- fixed(log(0.15)); label("Apparent CSF volume VCSF/F (L)") # Table 1 Humans VCSF/F = 0.15 L (Fixed); Methods: allometrically scaled from the macaque 0.015 L

    # Between-subject variability: Table 1 CV% converted to log-normal
    # variance omega^2 = log(1 + CV^2). No covariances are reported.
    etalcl ~ 0.29433 # Table 1 Humans CL/F CV 58.5%: log(1 + 0.585^2)
    etalvc ~ 0.14226 # Table 1 Humans V5/F CV 39.1%: log(1 + 0.391^2)
    etalka ~ 0.63917 # Table 1 Humans Ka CV 94.6%: log(1 + 0.946^2)
    etalvp ~ 0.12059 # Table 1 Humans V6/F CV 35.8%: log(1 + 0.358^2)
    etalq ~ 1.10406 # Table 1 Humans Q/F CV 142.0%: log(1 + 1.42^2)
    etalktr ~ 0.07444 # Table 1 Humans TAU CV 27.8%: log(1 + 0.278^2)
    etalk_central_csf ~ 0.58140 # Table 1 Humans Kplasma-CSF CV 88.8%: log(1 + 0.888^2)
    etalk_csf_central ~ 0.06541 # Table 1 Humans KCSF-plasma CV 26.0%: log(1 + 0.26^2)
    etalk_central_brain ~ 0.08728 # Table 1 Humans Kplasma-brain CV 30.2%: log(1 + 0.302^2)
    etalk_brain_central ~ 0.09289 # Table 1 Humans Kbrain-plasma CV 31.2%: log(1 + 0.312^2)
    etalk_csf_brain ~ 0.14980 # Table 1 Humans KCSF-brain CV 40.2%: log(1 + 0.402^2)
    etalk_brain_csf ~ 0.40348 # Table 1 Humans Kbrain-CSF CV 70.5%: log(1 + 0.705^2)
    etalv_brain ~ 0.06110 # Table 1 Humans Vbrain/F CV 25.1%: log(1 + 0.251^2)

    # Residual error: Methods state additive + proportional variance weighting
    # (Results call the error model log-normal); no magnitudes are reported,
    # so both components are encoded and fixed at zero.
    propSd <- fixed(0); label("Proportional residual SD, plasma (fraction); magnitude not reported") # Methods: additive and proportional variance components, values not reported
    addSd <- fixed(0); label("Additive residual SD, plasma (mg/L); magnitude not reported") # Methods: additive and proportional variance components, values not reported
    propSd_Ccsf <- fixed(0); label("Proportional residual SD, CSF (fraction); magnitude not reported") # Methods: additive and proportional variance components, values not reported
    addSd_Ccsf <- fixed(0); label("Additive residual SD, CSF (mg/L); magnitude not reported") # Methods: additive and proportional variance components, values not reported
  })
  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    ka <- exp(lka + etalka)
    vp <- exp(lvp + etalvp)
    q <- exp(lq + etalq)
    ktr <- exp(lktr + etalktr)
    k_central_csf <- exp(lk_central_csf + etalk_central_csf)
    k_csf_central <- exp(lk_csf_central + etalk_csf_central)
    k_central_brain <- exp(lk_central_brain + etalk_central_brain)
    k_brain_central <- exp(lk_brain_central + etalk_brain_central)
    k_csf_brain <- exp(lk_csf_brain + etalk_csf_brain)
    k_brain_csf <- exp(lk_brain_csf + etalk_brain_csf)
    v_brain <- exp(lv_brain + etalv_brain)
    vcsf <- exp(lvcsf)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Srinivas 2019 Figure 1: CMPT 1 depot, 2-4 transit, 5 central,
    # 6 peripheral, 7 CSF, 8 brain tissue.
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ka * transit3
    d/dt(central) <- ka * transit3 + k21 * peripheral1 + k_csf_central * csf +
      k_brain_central * brain -
      (kel + k12 + k_central_csf + k_central_brain) * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(csf) <- k_central_csf * central + k_brain_csf * brain -
      (k_csf_central + k_csf_brain) * csf
    d/dt(brain) <- k_central_brain * central + k_csf_brain * csf -
      (k_brain_central + k_brain_csf) * brain

    Cc <- central / vc
    Ccsf <- csf / vcsf
    Cbrain <- brain / v_brain

    Cc ~ add(addSd) + prop(propSd)
    Ccsf ~ add(addSd_Ccsf) + prop(propSd_Ccsf)
  })
}
