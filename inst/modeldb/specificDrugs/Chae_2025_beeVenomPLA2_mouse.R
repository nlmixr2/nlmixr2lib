Chae_2025_beeVenomPLA2_mouse <- function() {
  description <- paste(
    "Preclinical (mouse; ICR). Two-compartment population PK model with",
    "first-order absorption for bee venom phospholipase A2 (bvPLA2) after a",
    "single subcutaneous injection of new composition bee venom (NCBV).",
    "Apparent (per F) parameters are per kg body weight, so doses are entered",
    "in mg/kg. Between-animal variability on CL/F only; proportional residual",
    "error."
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
    species = "mouse (ICR)",
    n_subjects = 96L,
    n_studies = 1L,
    age_range = "7-8 weeks",
    weight_range = "20-35 g",
    sex_female_pct = 50,
    disease_state = "healthy",
    dose_range = "0.5 mg/kg NCBV, single subcutaneous injection",
    regions = "preclinical (Korea)",
    notes = paste(
      "Male and female ICR mice, parallel (terminal) cardiac-puncture sampling",
      "at 0-24 h, n = 4 per time point, 48 per sex (Chae 2025 Methods). Plasma",
      "bvPLA2 by ELISA (LLOQ 0.78 ng/mL). The NCBV product contains up to",
      "76.2% bvPLA2; doses are expressed as NCBV."
    )
  )

  ini({
    # Chae 2025 Table 2, 'Mice (0.5 mg/kg)' column. CL/F, V/F and Q/F are per kg.
    lka <- log(9.08); label("Absorption rate constant (1/h)") # Table 2: ka = 9.08 1/h (RSE 10.9%)
    lcl <- log(0.243); label("Apparent clearance CL/F (L/h/kg)") # Table 2: CL/F = 0.243 L/h/kg (RSE 56.4%)
    lvc <- log(0.912); label("Apparent central volume Vc/F (L/kg)") # Table 2: Vc/F = 0.912 L/kg (RSE 5.77%)
    lvp <- log(5.89); label("Apparent peripheral volume Vp/F (L/kg)") # Table 2: Vp/F = 5.89 L/kg (RSE 194%)
    lq <- log(0.121); label("Apparent intercompartmental clearance Q/F (L/h/kg)") # Table 2: Q/F = 0.121 L/h/kg (RSE 106%)

    # BSV reported as CV%; log-normal variance = log(CV^2 + 1).
    etalcl ~ 0.01806 # Table 2: BSV CL/F = 13.5% CV (RSE 116%)

    propSd <- 0.327; label("Proportional residual error (fraction)") # Table 2: sigma prop = 32.7% (RSE 6.42%)
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q <- exp(lq)

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
