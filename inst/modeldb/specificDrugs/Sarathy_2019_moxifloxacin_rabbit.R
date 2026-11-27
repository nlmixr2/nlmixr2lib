Sarathy_2019_moxifloxacin_rabbit <- function() {
  description <- "Preclinical (rabbit, New Zealand White, Mycobacterium tuberculosis HN878-infected). Moxifloxacin plasma and lung-lesion PK: a two-compartment oral plasma model (first-order absorption, linear clearance, per-kg parameters with lognormal between-animal variability) driving three non-depleting site-of-action compartments -- uninvolved lung, cellular areas of caseous granulomas, and caseum of caseous granulomas. Each site compartment holds a concentration and exchanges with plasma through a pair of first-order rate constants, d/dt(C_site) = k_central_site * Cc - k_site_central * C_site, so the steady-state site-to-plasma AUC ratio is k_central_site / k_site_central (7.3 lung, 9.3 cellular, 6.0 caseum). The site compartments are small and do not return mass to plasma (dashed arrows in Figure S3A)."
  reference <- paste(
    "Sarathy JP, Blanc L, Alvarez-Cabrera N, O'Brien P, Dias-Freedman I, Mina M, Zimmerman M, Kaya F, Ho Liang HP, Prideaux B, Dietzold J, Salgame P, Savic RM, Linderman J, Kirschner D, Pienaar E, Dartois V.",
    "Fluoroquinolone efficacy against tuberculosis is driven by penetration into lesions and activity against resident bacterial populations.",
    "Antimicrob Agents Chemother. 2019;63(5):e02516-18. doi:10.1128/AAC.02516-18. PMCID: PMC6496041.",
    "All parameter values: Supplemental Table S1. Model structure: Supplemental Figure S3A and Materials and Methods, 'Computational modeling methods' (iv).",
    "The plasma PK parameters in Table S1 are reprinted from Pienaar E, Sarathy J, Prideaux B, et al. Comparing efficacies of moxifloxacin, levofloxacin and gatifloxacin in tuberculosis granulomas using a multi-scale systems pharmacology approach. PLoS Comput Biol. 2017;13(8):e1005650. doi:10.1371/journal.pcbi.1005650 (Table 3; plasma ODEs in Methods, 'Plasma PK model').",
    sep = " "
  )
  vignette <- "Sarathy_2019_fluoroquinolones_rabbit"
  units <- list(time = "h", dosing = "mg/kg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "moxifloxacin", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "moxifloxacin", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "moxifloxacin", units = "mg/kg", specimen = "plasma", verified = TRUE),
    lung = list(analyte = "moxifloxacin", units = "mg/L", specimen = "tissue", verified = TRUE),
    lesion = list(analyte = "moxifloxacin", units = "mg/L", specimen = "tissue", verified = TRUE),
    caseum = list(analyte = "moxifloxacin", units = "mg/L", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "rabbit (New Zealand White, female)",
    n_subjects = 12L,
    n_studies = 1L,
    weight_range = "2.2-2.6 kg",
    sex_female_pct = 100,
    disease_state = "Pulmonary tuberculosis after nose-only aerosol infection with Mycobacterium tuberculosis HN878, with mature cellular and necrotic (caseous) lung lesions at the time of PK sampling",
    dose_range = "Oral gavage. Plasma PK model fit: 25 or 60 mg/kg daily for 21 doses. Lesion PK: single 100 mg/kg dose. Efficacy-study simulations: 60 mg/kg (human-equivalent of 400 mg) and 25 mg/kg (lower 5th percentile of human exposure) daily.",
    regions = "Public Health Research Institute, Rutgers New Jersey Medical School, Newark, NJ, USA",
    notes = paste(
      "The plasma PK model was fit by nonlinear mixed-effects modeling (SimBiology, stochastic EM, constant error model) to 12 infected rabbits (4 at 25 mg/kg and 8 at 60 mg/kg; Pienaar 2017 Table 6).",
      "The site-of-action rate constants were calibrated by nonlinear least squares, using the plasma PK parameters, to drug concentrations in uninvolved lung, cellular lesion areas and caseum measured by LC-MS/MS of homogenates and laser-capture microdissection after a single 100 mg/kg dose (7 rabbits for homogenate granuloma PK in Pienaar 2017 Table 6; Sarathy 2019 Data Set S1).",
      "Rabbit plasma unbound fraction 69.5 percent (human 71.7 percent; Sarathy 2019 Table S2)."
    )
  )

  ini({
    # ---- Plasma PK (Table S1, 'Rabbit Plasma PK parameters (mean, %CV)') ----
    # Values are the population typical values of a nonlinear mixed-effects
    # fit (Pienaar 2017 Table 3 footnote 3: 'mean and % coefficient of
    # variation determined by non-linear mixed effects modeling'). Per-kg
    # parameters; dose is mg/kg.
    lka <- log(0.9); label("Absorption rate constant (1/h)") # Table S1: ka = 0.9 1/h (CV 57%)
    lcl <- log(1.3); label("Plasma clearance (L/h/kg)") # Table S1: CL = 1.3 L/h/kg (CV 27%)
    lvc <- log(5.35); label("Plasma (central) volume of distribution (L/kg)") # Table S1: Vp = 5.35 L/kg (CV 9%)
    lvp <- log(80.8); label("Peripheral volume of distribution (L/kg)") # Table S1: Vpe = 80.8 L/kg (CV 0.9%)
    lq <- log(2.86); label("Intercompartmental clearance (L/h/kg)") # Table S1: Q = 2.86 L/h/kg (CV 12%)

    # ---- Site-of-action exchange (Table S1, 'Multi-compartment PK model
    # lung tissue PK parameters (baseline (95% CI))'). No between-animal
    # variability was estimated for these. ----
    lk_central_lung <- log(4.4); label("Distribution rate constant, plasma to uninvolved lung (1/h)") # Table S1: k pl-lung = 4.4 (3.9-5) 1/h
    lk_lung_central <- log(0.6); label("Distribution rate constant, uninvolved lung to plasma (1/h)") # Table S1: k lung-pl = 0.6 (0.5-0.7) 1/h
    lk_central_lesion <- log(6.5); label("Distribution rate constant, plasma to cellular areas of granulomas (1/h)") # Table S1: k pl-cell = 6.5 (5.9-7.1) 1/h
    lk_lesion_central <- log(0.7); label("Distribution rate constant, cellular areas of granulomas to plasma (1/h)") # Table S1: k cell-pl = 0.7 (0.6-0.8) 1/h
    lk_central_caseum <- log(1.8); label("Distribution rate constant, plasma to caseous areas of granulomas (1/h)") # Table S1: k pl-cas = 1.8 (1.5-2.1) 1/h
    lk_caseum_central <- log(0.3); label("Distribution rate constant, caseous areas of granulomas to plasma (1/h)") # Table S1: k cas-pl = 0.3 (0.2-0.4) 1/h

    # ---- Between-animal variability: lognormal, omega^2 = log(1 + CV^2) ----
    etalka ~ 0.28134 # Table S1: ka CV 57%
    etalcl ~ 0.070365 # Table S1: CL CV 27%
    etalvc ~ 0.0080674 # Table S1: Vp CV 9%
    etalvp ~ 0.000080997 # Table S1: Vpe CV 0.9%
    etalq ~ 0.014297 # Table S1: Q CV 12%

    # ---- Residual error ----
    # Pienaar 2017 Methods ('Plasma PK calibration') used a constant
    # (additive) error model but does not report its magnitude.
    addSd <- fixed(0); label("Additive residual error on plasma concentration (mg/L); magnitude not reported") # not reported in Sarathy 2019 or Pienaar 2017
  })

  model({
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q <- exp(lq + etalq)

    k_central_lung <- exp(lk_central_lung)
    k_lung_central <- exp(lk_lung_central)
    k_central_lesion <- exp(lk_central_lesion)
    k_lesion_central <- exp(lk_lesion_central)
    k_central_caseum <- exp(lk_central_caseum)
    k_caseum_central <- exp(lk_caseum_central)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Plasma: two-compartment model with first-order absorption (Pienaar 2017
    # 'Plasma PK model'); states are amounts per kg body weight.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc

    # Site-of-action compartments hold CONCENTRATIONS (mg/L) and are driven by
    # the plasma concentration without returning mass to plasma (Figure S3A:
    # 'small physiological compartments that do not significantly affect
    # plasma concentrations').
    d/dt(lung) <- k_central_lung * Cc - k_lung_central * lung
    d/dt(lesion) <- k_central_lesion * Cc - k_lesion_central * lesion
    d/dt(caseum) <- k_central_caseum * Cc - k_caseum_central * caseum

    Clung <- lung
    Clesion <- lesion
    Ccaseum <- caseum

    Cc ~ add(addSd)
  })
}
