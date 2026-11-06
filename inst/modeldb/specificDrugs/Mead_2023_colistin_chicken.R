Mead_2023_colistin_chicken <- function() {
  description <- "Preclinical (chicken). One-compartment population PK model of colistin in the small-intestine (jejunum + ileum) luminal content of Ross 308 broiler chickens: first-order transit from the crop (depot) into the intestinal lumen (central) and first-order intestinal outflow towards the colon, with between-bird variability on the intestinal clearance only and a multiplicative residual error. Doses by oral gavage (bolus into the crop) or via drinking water (zero-order over 24 h); colistin is not absorbed, so the observation Cc is the luminal content concentration, not a plasma concentration (Mead 2023)"
  reference <- "Mead A, Toutain PL, Richez P, Pelligand L. Targeted dosing for susceptible heteroresistant subpopulations may improve rational dosage regimen prediction for colistin in broiler chickens. Sci Rep. 2023;13:12822. doi:10.1038/s41598-023-39727-w"
  vignette <- "Mead_2023_colistin_chicken"
  # Doses and compartment amounts are mg colistin base per kg body weight
  # (75,000 IU/kg = 2.5 mg/kg at 30,000 IU/mg). Volume and clearance are in
  # g intestinal content per kg body weight, so Cc = 1000 * central / vc is mg
  # colistin per kg intestinal content -- the unit the paper reports.
  units <- list(time = "h", dosing = "mg/kg", concentration = "mg/kg")

  # The specimen vocabulary has no gastrointestinal-content matrix; `faeces`
  # is the vocabulary's luminal-content entry (as in the Echterhof 2026 phage
  # PBPK models), read here as small-intestine luminal digesta.
  compartmentData <- list(
    depot = list(analyte = "colistin", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "colistin", units = "mg/kg", specimen = "faeces", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "chicken (Ross 308 broiler)",
    n_subjects = 148L,
    n_studies = 1L,
    age_range = "13 days at the onset of dosing",
    weight_range = "mean bird weight per group and day 156.0-179.6 g (Supplementary Table S1)",
    sex_female_pct = NA_real_,
    disease_state = "Healthy, uninfected broilers (no E. coli challenge)",
    dose_range = "Colistin sulphate (ColiMeiji) via drinking water at nominal 75,000, 100,000 and 150,000 IU/kg BW/day for 3 days (actual back-calculated 78,023, 97,667 and 141,133 IU/kg BW/day), or a single oral gavage of 75,000 or 150,000 IU/kg BW into the crop",
    regions = "United Kingdom (Royal Veterinary College)",
    sampling = "Destructive sampling, one luminal-content sample per bird. Drinking water: 4 birds at 0, 12, 36, 54, 72, 75, 76, 78, 84 and 96 h. Gavage: 2 birds at 0, 1, 3, 6, 12 and 24 h. Jejunum + ileum content by UHPLC-MS/MS; LOQ 1.1 mg/kg (BLQ handled by M3).",
    notes = "175 birds housed: 38 per drinking-water group (n = 114), 17 per gavage group (n = 34) and 27 untreated controls; n_subjects counts the 148 birds allocated to the treated groups (the number of birds and samples in the final fit is not reported). Fit in Phoenix NLME 8.3 (Laplace). Sex not reported."
  )

  ini({
    lka <- log(0.221); label("Crop-to-small-intestine transit rate constant ka (1/h)") # Table 1, ka = 0.221 1/h
    lvc <- log(5.14); label("Apparent small-intestine content volume V (g intestinal content/kg BW)") # Table 1, V = 5.14 g intestinal content/kg BW
    lcl <- log(29.1); label("Intestinal clearance (luminal outflow) Cl (g intestinal content/kg BW/h)") # Table 1, Cl = 29.1 g intestinal content/kg BW/h
    lfdepot <- fixed(log(1)); label("Fraction of the crop dose reaching the small intestine (unitless)") # Results Eq. 1 text, F = 100% for the gut compartment (no absorption)

    etalcl ~ 0.1519 # Table 1, Cl BSV 40.5% CV; omega^2 = log(1 + 0.405^2) per Methods CV equation

    propSd <- 0.362; label("Multiplicative residual error (fraction)") # Table 1, multiplicative residual error 36.2%
  })
  model({
    ka <- exp(lka)
    vc <- exp(lvc)
    cl <- exp(lcl + etalcl)
    kel <- cl / vc

    d / dt(depot) <- -ka * depot
    d / dt(central) <- ka * depot - kel * central
    f(depot) <- exp(lfdepot)

    # vc is in g content/kg BW; divide by 1000 to give mg colistin per kg content.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
