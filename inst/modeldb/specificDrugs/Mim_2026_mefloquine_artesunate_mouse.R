Mim_2026_mefloquine_artesunate_mouse <- function() {
  description <- paste(
    "Preclinical (mouse, Plasmodium berghei-infected). Population PK/PD model of",
    "oral artesunate plus mefloquine (ASMQ) given as a single gavage dose. Each drug",
    "has its own two-compartment PK model with first-order absorption and linear",
    "elimination, fitted separately in Monolix (Mim 2026 Tables 1 and 2). The plasma",
    "concentrations of both drugs drive one turnover model of percent parasitemia",
    "in which each drug inhibits parasite production with an Imax = 1, Hill = 1",
    "function and the two inhibitions add (Equation 1, Table 3). The additive sum is",
    "capped at 1 so that production cannot become negative; with the cap the model",
    "reproduces the published combination simulations (Figure 4a), whereas the",
    "uncapped Equation 1 drives parasitemia below zero within 10 h. Between-animal",
    "variability was fixed to a negligible value (sparse destructive sampling) and is",
    "carried as fixed etas. Set the mefloquine dose to zero for artesunate",
    "monotherapy."
  )
  reference <- paste(
    "Mim SR, Vieira Santos V, Pereira LC, Pica K, Miranda ALLdS, Quadros HC,",
    "Moreira DRM, Azeredo FJ. (2026). Pharmacokinetic and pharmacodynamic modeling",
    "of anti-plasmodial drugs mefloquine plus artesunate: insights on translational",
    "application. Antimicrob Agents Chemother 70(3):e01717-25.",
    "doi:10.1128/aac.01717-25."
  )
  vignette <- "Mim_2026_mefloquine_artesunate"

  # Doses in mg/kg and volumes in L/kg (the paper's mL/kg / 1000), so
  # amount / volume is mg/L = ug/mL, the unit of the Table 3 IC50 values.
  units <- list(time = "h", dosing = "mg/kg", concentration = "mg/L")
  dosing <- c("depot_mefloquine", "depot_artesunate")

  # No covariate analysis was performed: Results, 'formal covariate analysis
  # was not performed, recognizing that IIV and residual unexplained
  # variability cannot be reliably separated in this data set'.
  covariateData <- list()

  compartmentData <- list(
    depot_mefloquine = list(analyte = "mefloquine", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central_mefloquine = list(analyte = "mefloquine", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1_mefloquine = list(analyte = "mefloquine", units = "mg/kg", specimen = "tissue", verified = TRUE),
    depot_artesunate = list(analyte = "artesunate", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central_artesunate = list(analyte = "artesunate", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1_artesunate = list(analyte = "artesunate", units = "mg/kg", specimen = "tissue", verified = TRUE),
    parasites = list(
      analyte = "Plasmodium berghei-infected erythrocytes",
      units = "% parasitemia",
      specimen = "whole blood",
      verified = TRUE
    )
  )

  population <- list(
    species = "mouse (Swiss, Plasmodium berghei ANKA gfp/luc-infected)",
    n_subjects = 44L,
    n_studies = 1L,
    age_range = "6-8 weeks",
    weight_range = "18-22 g",
    disease_state = paste(
      "Plasmodium berghei blood-stage malaria after intraperitoneal inoculation of",
      "1e6 parasitized erythrocytes; PK animals dosed once parasitemia reached 5% or",
      "more."
    ),
    dose_range = "Single oral (gavage) dose of artesunate 100 mg/kg plus mefloquine 55 mg/kg",
    regions = "Brazil (Instituto Goncalo Moniz, Fiocruz, Salvador)",
    notes = paste(
      "Male and female Swiss mice (proportions not reported). PK by destructive",
      "sampling, n = 4 per time point at 0.083, 0.25, 0.5, 1, 1.5, 2, 4, 8, 12, 24",
      "and 48 h (Methods, Pharmacokinetic experiments), giving 59 mefloquine and 34",
      "artesunate plasma samples (HPLC-DAD); n_subjects = 11 time points x 4 animals",
      "as described, although 59 mefloquine samples implies some time points had more.",
      "Parasitemia by flow cytometry (Table S1, experiment 2). Monolix 2024R1, SAEM."
    )
  )

  ini({
    # Mefloquine PK, Table 1. Volumes and clearances are per kg; the paper's
    # mL/kg and mL/h/kg values are divided by 1000.
    lka_mefloquine <- log(0.28); label("Mefloquine absorption rate constant Ka (1/h)") # Table 1: Ka = 0.28 1/h (RSE 24.4%)
    lcl_mefloquine <- log(0.28601); label("Mefloquine clearance CL (L/h/kg)") # Table 1: CL = 286.01 mL/h/kg (RSE 5.22%)
    lvc_mefloquine <- log(4.03963); label("Mefloquine central volume V1 (L/kg)") # Table 1: V1 = 4039.63 mL/kg (RSE 11.2%)
    lvp_mefloquine <- log(41.81736); label("Mefloquine peripheral volume V2 (L/kg)") # Table 1: V2 = 41817.36 mL/kg (RSE 25.7%), row mislabelled 'Volume of central compartment'
    lq_mefloquine <- log(20.63648); label("Mefloquine intercompartmental clearance Q (L/h/kg)") # Table 1: Q = 20636.48 mL/h/kg (RSE 36.0%)

    # Artesunate PK, Table 2.
    lka_artesunate <- fixed(log(0.6)); label("Artesunate absorption rate constant Ka (1/h)") # Table 2: Ka = 0.6 1/h (fixed); Methods: fixed from prior information
    lcl_artesunate <- log(3.28552); label("Artesunate clearance CL (L/h/kg)") # Table 2: CL = 3,285.52 mL/h/kg (RSE 11.9%)
    lvc_artesunate <- log(1.9214); label("Artesunate central volume V1 (L/kg)") # Table 2: V1 = 1,921.4 mL/kg (RSE 14.0%)
    lvp_artesunate <- fixed(log(11.3423)); label("Artesunate peripheral volume V2 (L/kg)") # Table 2: V2 = 11,342.3 mL/kg (fixed), row mislabelled 'Volume of central compartment'
    lq_artesunate <- log(19.21083); label("Artesunate intercompartmental clearance Q (L/h/kg)") # Table 2: Q = 19,210.83 mL/h/kg (RSE 6.39%)

    # Parasitemia turnover model, Equation 1 and Table 3.
    lrbase <- fixed(log(17.8)); label("Parasitemia at the start of treatment R0 (%)") # Table 3: R0 = 17.8% (fixed); Results: fixed from Burgert et al.
    lkout <- fixed(log(0.08)); label("First-order parasite elimination rate Kout (1/h)") # Table 3: Kout = 0.08 1/h (fixed); Results: fixed from Burgert et al.
    lic50_mefloquine <- log(0.011); label("Plasma mefloquine concentration for half-maximal inhibition of parasite production IC50 (mg/L)") # Table 3: IC50 (MQ) = 0.011 ug/mL (RSE 11.2%), 29.1 nM
    lic50_artesunate <- log(0.0044); label("Plasma artesunate concentration for half-maximal inhibition of parasite production IC50 (mg/L)") # Table 3: IC50 (AS) = 0.0044 ug/mL (RSE 4.62%), 10.93 nM

    # Between-animal variability. Results: 'inter-individual variability (IIV)
    # was fixed to a minimal value of 0.001'; Methods: fixed to a small value
    # under a log-normal distribution. Monolix parameterises omega as an SD,
    # so the variance is 0.001^2 = 1e-6. The paper does not list which
    # parameters carry it; every parameter is given one, the Monolix default.
    etalka_mefloquine ~ fixed(1e-6) # Results: IIV set to a minimal value of 0.001 (omega SD)
    etalcl_mefloquine ~ fixed(1e-6) # Results: IIV set to a minimal value of 0.001 (omega SD)
    etalvc_mefloquine ~ fixed(1e-6) # Results: IIV set to a minimal value of 0.001 (omega SD)
    etalvp_mefloquine ~ fixed(1e-6) # Results: IIV set to a minimal value of 0.001 (omega SD)
    etalq_mefloquine ~ fixed(1e-6) # Results: IIV set to a minimal value of 0.001 (omega SD)
    etalka_artesunate ~ fixed(1e-6) # Results: IIV set to a minimal value of 0.001 (omega SD)
    etalcl_artesunate ~ fixed(1e-6) # Results: IIV set to a minimal value of 0.001 (omega SD)
    etalvc_artesunate ~ fixed(1e-6) # Results: IIV set to a minimal value of 0.001 (omega SD)
    etalvp_artesunate ~ fixed(1e-6) # Results: IIV set to a minimal value of 0.001 (omega SD)
    etalq_artesunate ~ fixed(1e-6) # Results: IIV set to a minimal value of 0.001 (omega SD)
    etalrbase ~ fixed(1e-6) # Methods: PD IIV set to a small value (taken as 0.001 omega SD)
    etalkout ~ fixed(1e-6) # Methods: PD IIV set to a small value (taken as 0.001 omega SD)
    etalic50_mefloquine ~ fixed(1e-6) # Methods: PD IIV set to a small value (taken as 0.001 omega SD)
    etalic50_artesunate ~ fixed(1e-6) # Methods: PD IIV set to a small value (taken as 0.001 omega SD)

    # Residual error. Results: proportional for PK, additive for PD.
    propSd_mefloquine <- 0.71; label("Mefloquine proportional residual error (fraction)") # Table 1: proportional error model = 0.71 (RSE 18.5%)
    propSd_artesunate <- 0.52; label("Artesunate proportional residual error (fraction)") # Table 2: proportional error model = 0.52 (RSE 15.8%)
    addSd_parasites <- 0.52; label("Parasitemia additive residual error (% parasitemia)") # Table 3: constant error model = 0.52
  })

  model({
    # Individual parameters
    ka_mefloquine <- exp(lka_mefloquine + etalka_mefloquine)
    cl_mefloquine <- exp(lcl_mefloquine + etalcl_mefloquine)
    vc_mefloquine <- exp(lvc_mefloquine + etalvc_mefloquine)
    vp_mefloquine <- exp(lvp_mefloquine + etalvp_mefloquine)
    q_mefloquine <- exp(lq_mefloquine + etalq_mefloquine)

    ka_artesunate <- exp(lka_artesunate + etalka_artesunate)
    cl_artesunate <- exp(lcl_artesunate + etalcl_artesunate)
    vc_artesunate <- exp(lvc_artesunate + etalvc_artesunate)
    vp_artesunate <- exp(lvp_artesunate + etalvp_artesunate)
    q_artesunate <- exp(lq_artesunate + etalq_artesunate)

    rbase <- exp(lrbase + etalrbase)
    kout <- exp(lkout + etalkout)
    ic50_mefloquine <- exp(lic50_mefloquine + etalic50_mefloquine)
    ic50_artesunate <- exp(lic50_artesunate + etalic50_artesunate)
    # Equation 1 legend: Kin = R0 * Kout
    kin <- rbase * kout

    # Micro-constants
    kel_mefloquine <- cl_mefloquine / vc_mefloquine
    k12_mefloquine <- q_mefloquine / vc_mefloquine
    k21_mefloquine <- q_mefloquine / vp_mefloquine
    kel_artesunate <- cl_artesunate / vc_artesunate
    k12_artesunate <- q_artesunate / vc_artesunate
    k21_artesunate <- q_artesunate / vp_artesunate

    # Two-compartment PK with first-order absorption for each drug (Figure 5)
    d/dt(depot_mefloquine) <- -ka_mefloquine * depot_mefloquine
    d/dt(central_mefloquine) <- ka_mefloquine * depot_mefloquine -
      kel_mefloquine * central_mefloquine -
      k12_mefloquine * central_mefloquine +
      k21_mefloquine * peripheral1_mefloquine
    d/dt(peripheral1_mefloquine) <- k12_mefloquine * central_mefloquine -
      k21_mefloquine * peripheral1_mefloquine

    d/dt(depot_artesunate) <- -ka_artesunate * depot_artesunate
    d/dt(central_artesunate) <- ka_artesunate * depot_artesunate -
      kel_artesunate * central_artesunate -
      k12_artesunate * central_artesunate +
      k21_artesunate * peripheral1_artesunate
    d/dt(peripheral1_artesunate) <- k12_artesunate * central_artesunate -
      k21_artesunate * peripheral1_artesunate

    # Plasma concentrations (mg/L = ug/mL)
    Cc_mefloquine <- central_mefloquine / vc_mefloquine
    Cc_artesunate <- central_artesunate / vc_artesunate

    # Equation 1: each drug inhibits parasite production (Imax = 1, Hill = 1)
    # and the two inhibitions add. The sum is capped at 1: both drugs sit far
    # above their IC50 after a dose, so the uncapped sum approaches 2 and
    # production turns negative; the capped form reproduces Figure 4a.
    inh_mefloquine <- Cc_mefloquine / (ic50_mefloquine + Cc_mefloquine)
    inh_artesunate <- Cc_artesunate / (ic50_artesunate + Cc_artesunate)
    inh_total <- min(1, inh_mefloquine + inh_artesunate)

    d/dt(parasites) <- kin * (1 - inh_total) - kout * parasites
    parasites(0) <- rbase

    Cc_mefloquine ~ prop(propSd_mefloquine)
    Cc_artesunate ~ prop(propSd_artesunate)
    parasites ~ add(addSd_parasites)
  })
}
