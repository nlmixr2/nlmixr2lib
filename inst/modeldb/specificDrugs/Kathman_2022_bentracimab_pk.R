Kathman_2022_bentracimab_pk <- function() {
  description <- "Two-compartment population PK model for uncomplexed bentracimab (PB2452, the ticagrelor-neutralising monoclonal antibody Fab fragment) given alone as a 30-minute IV infusion to healthy volunteers who did not receive ticagrelor (Kathman 2022, first-in-human cohorts 1-3). Its typical values were carried, fixed, into the PB2452 arm of the joint ticagrelor / active-metabolite / PB2452 PK-PD model Kathman_2022_bentracimab."
  reference <- "Kathman SJ, Wheeler JJ, Bhatt DL, Arnold SE, Lee JS. Population pharmacokinetic-pharmacodynamic modeling of PB2452, a monoclonal antibody fragment being developed as a ticagrelor reversal agent, in healthy volunteers. CPT Pharmacometrics Syst Pharmacol. 2022;11(1):68-81. doi:10.1002/psp4.12734"
  vignette <- "Kathman_2022_bentracimab"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  covariateData <- list()

  compartmentData <- list(
    central = list(analyte = "bentracimab", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "bentracimab", units = "nmol", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = "12 randomised in cohorts 1-3 (9 PB2452 : 3 placebo; Kathman 2022 Table 1); the 9 PB2452-treated subjects contribute PB2452 concentrations",
    n_studies = 1,
    age_range = "18-50 years",
    weight_range = "not reported",
    sex_female_pct = NA_real_,
    disease_state = "Healthy volunteers not pretreated with ticagrelor",
    dose_range = "PB2452 0.1, 0.3 or 1.0 g as a single 30-minute IV infusion (cohorts 1, 2 and 3)",
    regions = "United States (single centre)",
    notes = "Phase I single-ascending-dose trial of PB2452. Concentrations were converted to nmol/L before modelling using the analyte molecular weight (Kathman 2022 Methods 'Data assembly'). Bayesian MCMC estimation in NONMEM 7.4."
  )

  ini({
    lcl <- 0.632; label("Clearance of uncomplexed PB2452 (L/h)") # Kathman 2022 Table 2, CL = EXP(THETA1), THETA1 = 0.632 (1.88 L/h in text)
    lvc <- 1.05; label("Central volume of PB2452 (L)") # Kathman 2022 Table 2, V1 = EXP(THETA2), THETA2 = 1.05 (2.86 L in text)
    lq <- -0.770; label("Intercompartmental clearance of PB2452 (L/h)") # Kathman 2022 Table 2, Q = EXP(THETA3), THETA3 = -0.770
    lvp <- 1.24; label("Peripheral volume of PB2452 (L)") # Kathman 2022 Table 2, V2 = EXP(THETA4), THETA4 = 1.24

    # IIV reported as CV%; variance = CV^2 (the convention the paper uses for its
    # fixed OMEGAs: 0.01 printed as '10%' and 0.0025 as '5%' in Table 3).
    etalcl ~ 0.142884 # Kathman 2022 Table 2, CL IIV 37.8% -> 0.378^2
    etalvc ~ 0.163216 # Kathman 2022 Table 2, V1 IIV 40.4% -> 0.404^2
    etalq ~ 0.183184 # Kathman 2022 Table 2, Q IIV 42.8% -> 0.428^2
    etalvp ~ 0.395641 # Kathman 2022 Table 2, V2 IIV 62.9% -> 0.629^2

    propSd <- 0.0711; label("Proportional residual error (fraction)") # Kathman 2022 Table 2, 'Residual variability: CV = 7.11%'
  })

  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
