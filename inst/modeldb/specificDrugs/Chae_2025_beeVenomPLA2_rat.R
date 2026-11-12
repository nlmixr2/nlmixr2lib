Chae_2025_beeVenomPLA2_rat <- function() {
  description <- paste(
    "Preclinical (rat; Sprague-Dawley). Two-compartment population PK model",
    "with first-order absorption for bee venom phospholipase A2 (bvPLA2) after",
    "a single subcutaneous injection of new composition bee venom (NCBV).",
    "Apparent (per F) parameters are per kg body weight, so doses are entered",
    "in mg/kg. Between-animal variability on CL/F, Q/F and ka; proportional",
    "residual error."
  )
  reference <- paste(
    "Chae SU, Min JS, Jo SJ, Lee CB, Park J, Bae SH, Bae SK.",
    "Prediction of first-in-human dose for new composition bee venom based on",
    "allometric scaling and pharmacokinetic modeling approach.",
    "Transl Clin Pharmacol. 2025;33(1):27-39. doi:10.12793/tcp.2025.33.e4",
    sep = " "
  )
  vignette <- "Chae_2025_beeVenomPLA2"

  units <- list(time = "h", dosing = "mg/kg", concentration = "ng/mL")

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "bvPLA2", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "bvPLA2", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "bvPLA2", units = "mg/kg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "rat (Sprague-Dawley)",
    n_subjects = 32L,
    n_studies = 1L,
    age_range = "7-9 weeks",
    weight_range = "200-330 g",
    sex_female_pct = 50,
    disease_state = "healthy",
    dose_range = "0.2, 0.5, 1 and 2 mg/kg NCBV, single subcutaneous injection",
    regions = "preclinical (Korea)",
    notes = paste(
      "Male and female SD rats, n = 4 per sex per dose group, serial carotid",
      "artery sampling at 0-24 h (Chae 2025 Methods). The Table 2 column header",
      "reads 'Rats (0.2-0.5 mg/kg)' while Methods, Table 1 and Figure 3 cover",
      "0.2-2 mg/kg. Plasma bvPLA2 by ELISA (LLOQ 0.78 ng/mL). Subcutaneous",
      "bioavailability against an IV reference was 11.9-18.4% (Table 1), but",
      "the model is parameterised in apparent (per F) terms."
    )
  )

  ini({
    # Chae 2025 Table 2, 'Rats' column. CL/F, V/F and Q/F are per kg (the
    # Results prose prints Vp/F as '6.89 L'; Table 2 gives L/kg).
    lka <- log(0.803); label("Absorption rate constant (1/h)") # Table 2: ka = 0.803 1/h (RSE 15.3%)
    lcl <- log(1.83); label("Apparent clearance CL/F (L/h/kg)") # Table 2: CL/F = 1.83 L/h/kg (RSE 12.3%)
    lvc <- log(0.712); label("Apparent central volume Vc/F (L/kg)") # Table 2: Vc/F = 0.712 L/kg (RSE 23.5%)
    lvp <- log(6.89); label("Apparent peripheral volume Vp/F (L/kg)") # Table 2: Vp/F = 6.89 L/kg (RSE 29.3%)
    lq <- log(0.862); label("Apparent intercompartmental clearance Q/F (L/h/kg)") # Table 2: Q/F = 0.862 L/h/kg (RSE 21.7%)

    # BSV reported as CV%; log-normal variance = log(CV^2 + 1).
    etalcl ~ 0.4186 # Table 2: BSV CL/F = 72.1% CV (RSE 28.2%)
    etalq ~ 0.3479 # Table 2: BSV Q/F = 64.5% CV (RSE 41.7%)
    etalka ~ 0.1181 # Table 2: BSV ka = 35.4% CV (RSE 39.2%)

    propSd <- 0.322; label("Proportional residual error (fraction)") # Table 2: sigma prop = 32.2% (RSE 6.68%)
  })

  model({
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q <- exp(lq + etalq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # central is in mg/kg and vc in L/kg, so central / vc is mg/L; x 1000 gives ng/mL.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
