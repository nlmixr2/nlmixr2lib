Srinivas_2019_efavirenz_monkey <- function() {
  description <- paste(
    "Preclinical (rhesus macaque). Eight-compartment PK model of oral",
    "efavirenz in plasma, cerebrospinal fluid (CSF) and brain tissue of adult",
    "rhesus macaques (Srinivas 2019, preclinical model). Dose enters a depot",
    "and passes through three transit compartments sharing one rate constant",
    "(TAU) before first-order absorption (Ka) into a two-compartment plasma",
    "disposition model with linear elimination; CSF and brain tissue exchange",
    "bidirectionally with plasma and with each other through first-order rate",
    "constants on amounts. Apparent clearance steps up by a multiplicative",
    "factor M from 24 h after the first dose (CL/F = CL/F * M for t > 24 h).",
    "Fitted to 79 plasma, CSF and brain tissue concentrations from 8 macaques",
    "dosed 200 mg once daily for 10 days, by iterative two-stage (IT2S)",
    "Bayesian estimation in ADAPT 5 with VCSF fixed at 0.015 L. Typical values",
    "are the geometric means of the individual estimates (Table 1).",
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
    species = "rhesus macaque (Macaca mulatta)",
    n_subjects = 8,
    n_studies = 1,
    disease_state = "healthy (n = 4) or RT-SHIVmac239-infected (n = 4); infection status did not change efavirenz concentrations",
    dose_range = paste(
      "efavirenz 200 mg orally once daily for 10 days, with tenofovir",
      "30 mg/kg SC daily, emtricitabine 16 mg/kg SC daily and raltegravir",
      "100 mg/kg orally twice daily",
      sep = " "
    ),
    regions = "United States (University of California Davis)",
    notes = paste(
      "Srinivas 2019 Methods and Results. Nine adult macaques were dosed;",
      "five were first infected with RT-SHIVmac239 and one infected animal",
      "developed liver failure and was excluded, leaving 8. Plasma at 0, 4, 10",
      "and 24 h after the first dose and 24 h post-dose on days 4, 6, 9 and 10",
      "(necropsy); CSF and brain tissue (median of four regions) at necropsy",
      "only. Age, weight and sex are not reported.",
      sep = " "
    )
  )

  ini({
    # Srinivas 2019 Table 1, 'Rhesus macaques' rows: geometric mean of the
    # individual IT2S estimates (typical value) and CV% (between-animal
    # variability). Parameters are apparent (/F). Rate constants act on
    # amounts (1/h).
    lcl <- log(8.70); label("Apparent clearance CL/F during the first 24 h (L/h)") # Table 1 Rhesus macaques CL/F geometric mean 8.70 L/h (CV 89.8%)
    lvc <- log(6.45); label("Apparent central volume V5/F (L)") # Table 1 Rhesus macaques V5/F geometric mean 6.45 L (CV 51.2%)
    lka <- log(0.183); label("Absorption rate constant from transit3 to central Ka (1/h)") # Table 1 Rhesus macaques Ka geometric mean 0.183 1/h (CV 36.4%)
    lvp <- log(14.6); label("Apparent peripheral volume V6/F (L)") # Table 1 Rhesus macaques V6/F geometric mean 14.6 L (CV 155.0%)
    lq <- log(4.66); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 1 Rhesus macaques Q/F geometric mean 4.66 L/h (CV 89.1%)
    lktr <- log(0.95); label("Transit rate constant TAU, depot -> transit1 -> transit2 -> transit3 (1/h)") # Table 1 Rhesus macaques TAU geometric mean 0.95 1/h (CV 48.5%)
    # The clearance step is encoded as the piecewise-constant canonical pair
    # (ltclchange, lcl_late). Table 1 reports the multiplier M = 2.63, so the
    # late arm is CL/F * M = 8.70 * 2.63 L/h, and it carries both the CL/F and
    # the M random effects (etalcl + etalcl_late), which reproduces
    # CL_late = CL * M exactly.
    ltclchange <- fixed(log(24)); label("Time after the first dose at which CL/F steps up (h)") # Results and Figure 1 legend: 'if time > 24 hours, then CL/F = CL/F x M'
    lcl_late <- log(8.70 * 2.63); label("Apparent clearance CL/F after 24 h, CL/F x M (L/h)") # Table 1 Rhesus macaques CL/F 8.70 L/h x M geometric mean 2.63 (M CV 103.0%)
    lk_central_csf <- log(1.71e-6); label("Plasma-to-CSF rate constant Kplasma-CSF (1/h)") # Table 1 Rhesus macaques Kplasma-CSF geometric mean 1.71E-06 1/h (CV 59.5%)
    lk_csf_central <- log(0.18); label("CSF-to-plasma rate constant KCSF-plasma (1/h)") # Table 1 Rhesus macaques KCSF-plasma geometric mean 0.18 1/h (CV 27.3%)
    lk_central_brain <- log(2.42e-3); label("Plasma-to-brain-tissue rate constant Kplasma-brain (1/h)") # Table 1 Rhesus macaques Kplasma-brain geometric mean 2.42E-03 1/h (CV 31.1%)
    lk_brain_central <- log(0.15); label("Brain-tissue-to-plasma rate constant Kbrain-plasma (1/h)") # Table 1 Rhesus macaques Kbrain-plasma geometric mean 0.15 1/h (CV 26.8%)
    lk_csf_brain <- log(2.56e-3); label("CSF-to-brain-tissue rate constant KCSF-brain (1/h)") # Table 1 Rhesus macaques KCSF-brain geometric mean 2.56E-03 1/h (CV 36.0%)
    lk_brain_csf <- log(5.58e-6); label("Brain-tissue-to-CSF rate constant Kbrain-CSF (1/h)") # Table 1 Rhesus macaques Kbrain-CSF geometric mean 5.58E-06 1/h (CV 68.2%)
    lv_brain <- log(0.037); label("Apparent brain tissue volume Vbrain/F, conditioned on VCSF (L)") # Table 1 Rhesus macaques Vbrain/F geometric mean 0.037 L (CV 21.8%)
    lvcsf <- fixed(log(0.015)); label("Apparent CSF volume VCSF/F (L)") # Table 1 Rhesus macaques VCSF/F = 0.015 L (fixed); Methods: from the literature (ref 17)

    # Between-animal variability: Table 1 CV% converted to log-normal
    # variance omega^2 = log(1 + CV^2). No covariances are reported.
    etalcl ~ 0.59134 # Table 1 Rhesus macaques CL/F CV 89.8%: log(1 + 0.898^2)
    etalvc ~ 0.23281 # Table 1 Rhesus macaques V5/F CV 51.2%: log(1 + 0.512^2)
    etalka ~ 0.12442 # Table 1 Rhesus macaques Ka CV 36.4%: log(1 + 0.364^2)
    etalvp ~ 1.22451 # Table 1 Rhesus macaques V6/F CV 155.0%: log(1 + 1.55^2)
    etalq ~ 0.58438 # Table 1 Rhesus macaques Q/F CV 89.1%: log(1 + 0.891^2)
    etalktr ~ 0.21125 # Table 1 Rhesus macaques TAU CV 48.5%: log(1 + 0.485^2)
    etalcl_late ~ 0.72314 # Table 1 Rhesus macaques M CV 103.0%: log(1 + 1.03^2)
    etalk_central_csf ~ 0.30308 # Table 1 Rhesus macaques Kplasma-CSF CV 59.5%: log(1 + 0.595^2)
    etalk_csf_central ~ 0.07188 # Table 1 Rhesus macaques KCSF-plasma CV 27.3%: log(1 + 0.273^2)
    etalk_central_brain ~ 0.09232 # Table 1 Rhesus macaques Kplasma-brain CV 31.1%: log(1 + 0.311^2)
    etalk_brain_central ~ 0.06936 # Table 1 Rhesus macaques Kbrain-plasma CV 26.8%: log(1 + 0.268^2)
    etalk_csf_brain ~ 0.12186 # Table 1 Rhesus macaques KCSF-brain CV 36.0%: log(1 + 0.36^2)
    etalk_brain_csf ~ 0.38194 # Table 1 Rhesus macaques Kbrain-CSF CV 68.2%: log(1 + 0.682^2)
    etalv_brain ~ 0.04643 # Table 1 Rhesus macaques Vbrain/F CV 21.8%: log(1 + 0.218^2)

    # Residual error: Methods state additive + proportional variance weighting
    # (Results call the error model log-normal); no magnitudes are reported,
    # so both components are encoded and fixed at zero for each matrix.
    propSd <- fixed(0); label("Proportional residual SD, plasma (fraction); magnitude not reported") # Methods: additive and proportional variance components, values not reported
    addSd <- fixed(0); label("Additive residual SD, plasma (mg/L); magnitude not reported") # Methods: additive and proportional variance components, values not reported
    propSd_Ccsf <- fixed(0); label("Proportional residual SD, CSF (fraction); magnitude not reported") # Methods: additive and proportional variance components, values not reported
    addSd_Ccsf <- fixed(0); label("Additive residual SD, CSF (mg/L); magnitude not reported") # Methods: additive and proportional variance components, values not reported
    propSd_Cbrain <- fixed(0); label("Proportional residual SD, brain tissue (fraction); magnitude not reported") # Methods: additive and proportional variance components, values not reported
    addSd_Cbrain <- fixed(0); label("Additive residual SD, brain tissue (mg/L); magnitude not reported") # Methods: additive and proportional variance components, values not reported
  })
  model({
    tclchange <- exp(ltclchange)
    cl_early <- exp(lcl + etalcl)
    # CL_late = CL * M: the M random effect is mu-referenced on lcl_late and
    # the CL random effect carries over multiplicatively.
    cl_late_m <- exp(lcl_late + etalcl_late)
    cl_late <- cl_late_m * exp(etalcl)
    # Time is measured from the first dose (Results: 'if time > 24 hours').
    cl <- cl_early * (t <= tclchange) + cl_late * (t > tclchange)
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
    Cbrain ~ add(addSd_Cbrain) + prop(propSd_Cbrain)
  })
}
