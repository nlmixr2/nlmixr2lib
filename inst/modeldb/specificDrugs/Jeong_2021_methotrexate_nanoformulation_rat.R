Jeong_2021_methotrexate_nanoformulation_rat <- function() {
  description <- "Preclinical (rat). One-compartment population PK model for methotrexate given orally or intravenously to Sprague-Dawley rats as a free solution or as a nanoformulation (pooled PLGA nanoparticles and olive-oil/Labrasol nanoemulsions) (Jeong 2021, free-solution-vs-nanoformulation model), with first-order absorption, oral bioavailability, log-additive residual error, and a binary nanoformulation indicator acting as a linear deviation on V, CL, ka and F."
  reference <- "Jeong S-H, Jang J-H, Lee Y-B. Pharmacokinetic Comparison between Methotrexate-Loaded Nanoparticles and Nanoemulsions as Hard- and Soft-Type Nanoformulations: A Population Pharmacokinetic Modeling Approach. Pharmaceutics. 2021;13(7):1050. doi:10.3390/pharmaceutics13071050"
  vignette <- "Jeong_2021_methotrexate_nanoformulations"
  units <- list(time = "h", dosing = "mg/kg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "methotrexate", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "methotrexate", units = "mg/kg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    FORM_MTX_NANO = list(
      description = "Methotrexate nanoformulation indicator: 1 = methotrexate given inside a nanoformulation (PLGA nanoparticles or olive-oil/Labrasol nanoemulsions, pooled), 0 = free methotrexate solution",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (free methotrexate solution)",
      notes = "Per-dose-record (per-animal) indicator; each rat received one formulation by one route. Jeong 2021 Equations (1)-(5): 'formulation type Free methotrexate = 0; Nanoformulation = 1'. The nanoparticle and nanoemulsion arms are pooled into the value 1 in this model; the separate nanoparticle-vs-nanoemulsion model is Jeong_2021_methotrexate_nanoemulsion_rat. Doses per arm (Figures 1, 2 and Table 6): free solution 5 mg/kg oral and IV; PLGA nanoparticles 5 mg/kg oral and IV; nanoemulsions 0.06 mg/kg oral and 0.024 mg/kg IV. The authors fitted no dose effect, so the 83- to 208-fold lower nanoemulsion doses are absorbed into this indicator's coefficients.",
      source_name = "formulation type"
    )
  )

  population <- list(
    species = "rat (Sprague-Dawley)",
    n_subjects = 40L,
    n_studies = 2L,
    age_range = "7-9 weeks",
    weight_range = "240-260 g",
    sex_female_pct = 0,
    disease_state = "Healthy male rats (no disease model).",
    dose_range = "Single dose. Free methotrexate solution 5 mg/kg oral or IV; methotrexate-loaded PLGA nanoparticles 5 mg/kg oral or IV; methotrexate-loaded nanoemulsions 0.06 mg/kg oral or 0.024 mg/kg IV (n = 5 per arm).",
    regions = "Republic of Korea (Chonnam National University, Gwangju).",
    formulations = "PLGA nanoparticles 163.7 +/- 10.25 nm, zeta -20.4 mV, encapsulation 93.3%; olive-oil/Labrasol nanoemulsions 173.77 +/- 5.76 nm, zeta -35.63 mV, encapsulation 90.37% (Jeong 2021 Table 5).",
    notes = "Pooled reanalysis of the plasma data of two earlier formulation studies from the same group (Jang 2019 PLGA nanoparticles, Jang 2020 nanoemulsions), each with its own free-solution control arms: 8 arms x 5 rats = 40 rats (Jeong 2021 Section 2.1). Plasma methotrexate by UHPLC-ESI-MS/MS. Phoenix NLME 8.3, FOCE-ELS with eta-epsilon interaction."
  )

  ini({
    # Structural parameters -- Jeong 2021 Table 3 (final free-solution-vs-
    # nanoformulation model); reference = free methotrexate solution
    # (FORM_MTX_NANO = 0). Doses are per kg, so V is in L/kg and CL in L/h/kg.
    lvc <- log(14.889); label("Volume of distribution V for free methotrexate solution (L/kg)") # Table 3 tvV = 14.889 L/kg
    lcl <- log(14.577); label("Clearance CL for free methotrexate solution (L/h/kg)") # Table 3 tvCL = 14.577 L/h/kg
    lka <- log(0.582); label("First-order absorption rate constant ka for free methotrexate solution (1/h)") # Table 3 tvKa = 0.582 1/h
    lfdepot <- log(0.272); label("Oral bioavailability F for free methotrexate solution (fraction)") # Table 3 tvF = 0.272

    # Nanoformulation effects -- Jeong 2021 Table 3 and Equations (1), (2),
    # (4), (5): P = tvP * (1 + dPdFormulation * FORM_MTX_NANO) * exp(eta).
    e_form_mtx_nano_vc <- 0.429; label("Linear-deviation effect of nanoformulation on V (fraction)") # Table 3 dVdFormulation = 0.429
    e_form_mtx_nano_cl <- -0.355; label("Linear-deviation effect of nanoformulation on CL (fraction)") # Table 3 dCLdFormulation = -0.355
    e_form_mtx_nano_ka <- 10.883; label("Linear-deviation effect of nanoformulation on ka (fraction)") # Table 3 dKadFormulation = 10.883
    e_form_mtx_nano_fdepot <- 4.246; label("Linear-deviation effect of nanoformulation on F (fraction)") # Table 3 dFdFormulation = 4.246

    # IIV -- exponential (Jeong 2021 Section 2.3). Table 3 prints omega^2 to
    # three decimals and its IIV (%) column is 100 * sqrt(omega^2) (CL:
    # sqrt(0.338) = 58.14% vs printed 58.130%). Where omega^2 prints as 0.000
    # the variance is back-solved from the IIV (%) column.
    etalvc ~ 6.3504e-6 # Table 3 omega^2 V printed 0.000; IIV 0.252% -> (0.00252)^2
    etalcl ~ 0.338 # Table 3 omega^2 CL = 0.338 (IIV 58.130%)
    etalka ~ 4e-10 # Table 3 omega^2 Ka printed 0.000; IIV 0.002% -> (0.00002)^2
    etalfdepot ~ 0.386 # Table 3 omega^2 F = 0.386 (IIV 62.161%)

    # Residual error -- log-additive (Jeong 2021 Table 1 model 02-04); Phoenix
    # reports sigma as the standard deviation on the log scale.
    expSd <- 1.269; label("Log-additive residual error SD (log scale)") # Table 3 sigma = 1.269
  })

  model({
    # Individual parameters (Jeong 2021 Equations (1), (2), (4), (5)). The
    # absorption lag time (Equation (3)) is omitted: Table 3 prints
    # tvTlag = 0.000 h with SE 0.000, so the lag is below 0.0005 h.
    vc <- exp(lvc + etalvc) * (1 + e_form_mtx_nano_vc * FORM_MTX_NANO)
    cl <- exp(lcl + etalcl) * (1 + e_form_mtx_nano_cl * FORM_MTX_NANO)
    ka <- exp(lka + etalka) * (1 + e_form_mtx_nano_ka * FORM_MTX_NANO)
    fdepot <- exp(lfdepot + etalfdepot) * (1 + e_form_mtx_nano_fdepot * FORM_MTX_NANO)

    kel <- cl / vc

    # One-compartment disposition; oral doses into depot, IV doses into
    # central (Jeong 2021 Figure 7).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    f(depot) <- fdepot

    # Dose in mg/kg and V in L/kg give mg/L; x 1000 for ng/mL (the paper's
    # concentration unit).
    Cc <- 1000 * central / vc
    Cc ~ lnorm(expSd)
  })
}
