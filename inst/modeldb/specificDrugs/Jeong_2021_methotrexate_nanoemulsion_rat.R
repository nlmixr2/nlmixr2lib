Jeong_2021_methotrexate_nanoemulsion_rat <- function() {
  description <- "Preclinical (rat). One-compartment population PK model for methotrexate given orally or intravenously to Sprague-Dawley rats in PLGA nanoparticles (hard-type) or olive-oil/Labrasol nanoemulsions (soft-type) (Jeong 2021, nanoparticle-vs-nanoemulsion model), with first-order absorption, oral bioavailability, log-additive residual error, and a binary nanoemulsion indicator acting as a linear deviation on V, CL, ka and F."
  reference <- "Jeong S-H, Jang J-H, Lee Y-B. Pharmacokinetic Comparison between Methotrexate-Loaded Nanoparticles and Nanoemulsions as Hard- and Soft-Type Nanoformulations: A Population Pharmacokinetic Modeling Approach. Pharmaceutics. 2021;13(7):1050. doi:10.3390/pharmaceutics13071050"
  vignette <- "Jeong_2021_methotrexate_nanoformulations"
  units <- list(time = "h", dosing = "mg/kg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "methotrexate", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "methotrexate", units = "mg/kg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    FORM_MTX_NANOEMULSION = list(
      description = "Methotrexate nanoemulsion indicator: 1 = methotrexate given in olive-oil/Labrasol nanoemulsions (soft-type nanoformulation), 0 = given in PLGA nanoparticles (hard-type nanoformulation)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (methotrexate-loaded PLGA nanoparticles)",
      notes = "Per-dose-record (per-animal) indicator; each rat received one formulation by one route. Jeong 2021 Equations (6)-(10): 'formulation type Nanoparticles = 0; Nanoemulsions = 1'. Free methotrexate solution is outside the scope of this model (see Jeong_2021_methotrexate_nanoformulation_rat). Doses per arm (Table 6): nanoparticles 5 mg/kg oral and IV; nanoemulsions 0.06 mg/kg oral and 0.024 mg/kg IV. The authors fitted no dose effect, so the 83- to 208-fold lower nanoemulsion doses are absorbed into this indicator's coefficients.",
      source_name = "formulation type"
    )
  )

  population <- list(
    species = "rat (Sprague-Dawley)",
    n_subjects = 20L,
    n_studies = 2L,
    age_range = "7-9 weeks",
    weight_range = "240-260 g",
    sex_female_pct = 0,
    disease_state = "Healthy male rats (no disease model).",
    dose_range = "Single dose. Methotrexate-loaded PLGA nanoparticles 5 mg/kg oral or IV; methotrexate-loaded nanoemulsions 0.06 mg/kg oral or 0.024 mg/kg IV (n = 5 per arm).",
    regions = "Republic of Korea (Chonnam National University, Gwangju).",
    formulations = "PLGA nanoparticles 163.7 +/- 10.25 nm, zeta -20.4 mV, encapsulation 93.3%; olive-oil/Labrasol nanoemulsions 173.77 +/- 5.76 nm, zeta -35.63 mV, encapsulation 90.37% (Jeong 2021 Table 5).",
    notes = "Nanoformulation arms of two earlier formulation studies from the same group (Jang 2019 PLGA nanoparticles, Jang 2020 nanoemulsions): 4 arms x 5 rats = 20 rats (Jeong 2021 Sections 2.1 and 3.4, Table 6). Plasma methotrexate by UHPLC-ESI-MS/MS. Phoenix NLME 8.3, FOCE-ELS with eta-epsilon interaction."
  )

  ini({
    # Structural parameters -- Jeong 2021 Table 9 (final nanoparticle-vs-
    # nanoemulsion model); reference = PLGA nanoparticles
    # (FORM_MTX_NANOEMULSION = 0). Doses are per kg, so V is in L/kg and CL in
    # L/h/kg.
    lvc <- log(18.832); label("Volume of distribution V for PLGA nanoparticles (L/kg)") # Table 9 tvV = 18.832 L/kg
    lcl <- log(9.167); label("Clearance CL for PLGA nanoparticles (L/h/kg)") # Table 9 tvCL = 9.167 L/h/kg
    lka <- log(0.714); label("First-order absorption rate constant ka for PLGA nanoparticles (1/h)") # Table 9 tvKa = 0.714 1/h
    lfdepot <- log(0.334); label("Oral bioavailability F for PLGA nanoparticles (fraction)") # Table 9 tvF = 0.334

    # Nanoemulsion effects -- Jeong 2021 Table 9 and Equations (6), (7), (9),
    # (10): P = tvP * (1 + dPdFormulation * FORM_MTX_NANOEMULSION) * exp(eta).
    e_form_mtx_nanoemulsion_vc <- -0.986; label("Linear-deviation effect of nanoemulsion on V (fraction)") # Table 9 dVdFormulation = -0.986
    e_form_mtx_nanoemulsion_cl <- -0.990; label("Linear-deviation effect of nanoemulsion on CL (fraction)") # Table 9 dCLdFormulation = -0.990
    e_form_mtx_nanoemulsion_ka <- 1.552; label("Linear-deviation effect of nanoemulsion on ka (fraction)") # Table 9 dKadFormulation = 1.552
    e_form_mtx_nanoemulsion_fdepot <- 0.193; label("Linear-deviation effect of nanoemulsion on F (fraction)") # Table 9 dFdFormulation = 0.193

    # IIV -- exponential (Jeong 2021 Section 2.3). Table 9 prints omega^2 to
    # three decimals, which leaves at most one significant figure here; the
    # IIV (%) column is 100 * sqrt(omega^2) (Table 3 check: sqrt(0.338) =
    # 58.14% vs 58.130%), so every variance is back-solved from IIV (%).
    etalvc ~ 5.6644e-6 # Table 9 omega^2 V printed 0.000; IIV 0.238% -> (0.00238)^2
    etalcl ~ 0.0084622 # Table 9 omega^2 CL printed 0.008; IIV 9.199% -> (0.09199)^2
    etalka ~ 4e-10 # Table 9 omega^2 Ka printed 0.000; IIV 0.002% -> (0.00002)^2
    etalfdepot ~ 0.0010278 # Table 9 omega^2 F printed 0.001; IIV 3.206% -> (0.03206)^2

    # Residual error -- log-additive (Jeong 2021 Table 7 model 02-04); Phoenix
    # reports sigma as the standard deviation on the log scale.
    expSd <- 0.552; label("Log-additive residual error SD (log scale)") # Table 9 sigma = 0.552
  })

  model({
    # Individual parameters (Jeong 2021 Equations (6), (7), (9), (10)). The
    # absorption lag time (Equation (8)) is omitted: Table 9 prints
    # tvTlag = 0.000 h with SE 0.000, so the lag is below 0.0005 h.
    vc <- exp(lvc + etalvc) * (1 + e_form_mtx_nanoemulsion_vc * FORM_MTX_NANOEMULSION)
    cl <- exp(lcl + etalcl) * (1 + e_form_mtx_nanoemulsion_cl * FORM_MTX_NANOEMULSION)
    ka <- exp(lka + etalka) * (1 + e_form_mtx_nanoemulsion_ka * FORM_MTX_NANOEMULSION)
    fdepot <- exp(lfdepot + etalfdepot) * (1 + e_form_mtx_nanoemulsion_fdepot * FORM_MTX_NANOEMULSION)

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
