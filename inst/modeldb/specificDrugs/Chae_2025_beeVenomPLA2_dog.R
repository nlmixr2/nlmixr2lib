Chae_2025_beeVenomPLA2_dog <- function() {
  description <- paste(
    "Preclinical (beagle dog). Two-compartment population PK model with",
    "first-order absorption for bee venom phospholipase A2 (bvPLA2) after a",
    "single subcutaneous injection of new composition bee venom (NCBV).",
    "Apparent (per F) parameters are per kg body weight, so doses are entered",
    "in mg/kg. Between-animal variability on CL/F and Vc/F; proportional",
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
    species = "beagle dog",
    n_subjects = 24L,
    n_studies = 1L,
    sex_female_pct = 50,
    disease_state = "healthy",
    dose_range = "0.3, 0.6 and 0.9 mg/kg NCBV, single subcutaneous injection",
    regions = "preclinical (Korea)",
    notes = paste(
      "Beagle dogs, n = 4 per sex per dose group, cephalic-vein sampling at",
      "0-24 h (Chae 2025 Methods); study run by ChemOn. Body weights are not",
      "reported. Plasma bvPLA2 by ELISA."
    )
  )

  ini({
    # Chae 2025 Table 2, 'Beagle dogs (0.3-0.9 mg/kg)' column. CL/F, V/F and Q/F are per kg.
    lka <- log(0.544); label("Absorption rate constant (1/h)") # Table 2: ka = 0.544 1/h (RSE 52.0%)
    lcl <- log(0.150); label("Apparent clearance CL/F (L/h/kg)") # Table 2: CL/F = 0.150 L/h/kg (RSE 19.4%)
    lvc <- log(0.385); label("Apparent central volume Vc/F (L/kg)") # Table 2: Vc/F = 0.385 L/kg (RSE 53.0%)
    lvp <- log(1.47); label("Apparent peripheral volume Vp/F (L/kg)") # Table 2: Vp/F = 1.47 L/kg (RSE 40.5%)
    lq <- log(0.141); label("Apparent intercompartmental clearance Q/F (L/h/kg)") # Table 2: Q/F = 0.141 L/h/kg (RSE 22.1%)

    # BSV reported as CV%; log-normal variance = log(CV^2 + 1).
    etalcl ~ 0.07393 # Table 2: BSV CL/F = 27.7% CV (RSE 41.4%)
    etalvc ~ 0.1547 # Table 2: BSV Vc/F = 40.9% CV (RSE 32.3%)

    propSd <- 0.166; label("Proportional residual error (fraction)") # Table 2: sigma prop = 16.6% (RSE 6.08%)
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
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
