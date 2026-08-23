Derippe_2024_venetoclax_mouse <- function() {
  description <- "Preclinical (mouse). One-compartment population PK model with first-order absorption for venetoclax in mice, fit to two pooled literature datasets: single 10 mg/kg oral doses in male and female mice (Eisenmann 2020) and venetoclax appearing after a 5 mg/kg intravenous dose of the prodrug ABBV-167 (Salem 2021). Sex is a covariate on volume of distribution and study is a covariate on the absorption rate constant; in the prodrug study the volume term also absorbs the fraction biotransformed and the absorption rate constant is the biotransformation rate constant. This is the PK layer feeding the paper's apoptosis QSP model and minimal agent-based model."
  reference <- paste(
    "Derippe T, Fouliard S, Decleves X, Mager DE.",
    "Quantitative systems pharmacology modeling of tumor heterogeneity in response to",
    "BH3-mimetics using virtual tumors calibrated with cell viability assays.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(7):1252-1263. doi:10.1002/psp4.13158.",
    "Parameter values are from the Supporting Information section 'Mice PK modeling'",
    "(supplementary file PSP4-13-1252-s001.docx).",
    "Source PK data digitized from Eisenmann ED, Jin Y, Weber RH, Sparreboom A, Baker SD,",
    "J Chromatogr B. 2020;1152:122176, and Salem AH, Tao ZF, Bueno OF, et al.,",
    "Mol Cancer Ther. 2021;20(6):999-1008.",
    sep = " "
  )
  vignette <- "Derippe_2024_BH3_mimetics_virtual_tumors"
  units <- list(time = "h", dosing = "mg/kg", concentration = "mg/L")

  covariateData <- list(
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female)",
      notes              = "Female is the REFERENCE here, which is the reverse of the usual nlmixr2lib pattern: the Supplement's 'Mice PK modeling' table lists Female first with the un-flagged V (6.54 L/kg) and reports the male value as a covariate effect (14.6% RSE). Applies only to the Eisenmann 2020 cohort (n = 10 male, n = 5 female); sex was not reported for the Salem 2021 cohort, which the paper treats as its own 'unknown' level. Set SEXF = 0 for Salem rows -- the STUDY_SALEM term neutralises the sex term there. Subject-level (time-fixed).",
      source_name        = "Sex"
    ),
    STUDY_SALEM = list(
      description        = "Salem 2021 ABBV-167 prodrug study indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Eisenmann 2020 oral venetoclax study)",
      notes              = "1 = mouse from Salem 2021 (5 mg/kg ABBV-167 prodrug IV, venetoclax measured as the biotransformation product); 0 = mouse from Eisenmann 2020 (10 mg/kg venetoclax PO). The parameters change meaning across the two studies: for Eisenmann, ka is first-order oral absorption and V is the volume of distribution; for Salem, ka is the prodrug-to-venetoclax biotransformation rate constant and V additionally absorbs the fraction biotransformed. Subject-level (time-fixed).",
      source_name        = "Study"
    )
  )

  covariatesDataExcluded <- list()

  compartmentData <- list(
    depot   = list(analyte = "venetoclax", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "venetoclax", units = "mg/kg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "mouse",
    n_subjects     = 15L,
    n_studies      = 2L,
    sex_female_pct = 33.3,
    weight_range   = "20-30 g (Eisenmann 2020 cohort; not reported for Salem 2021)",
    disease_state  = "Non-tumor-bearing mice",
    dose_range     = "10 mg/kg venetoclax PO single dose (Eisenmann 2020); 5 mg/kg ABBV-167 prodrug IV single dose (Salem 2021), from which only the venetoclax profile was used",
    regions        = "Preclinical (literature-digitized)",
    notes          = "n_subjects = 15 counts the Eisenmann 2020 cohort (10 male, 5 female); the Salem 2021 cohort size is not reported, so sex_female_pct is computed over the Eisenmann cohort only. Three profiles in total were fit jointly. The Eisenmann profiles cover only 0-6 h and were still at a plateau at 6 h, so the terminal phase is informed almost entirely by the Salem profile (0-24 h); the paper notes the terminal half-life is overpredicted. Absence of xenografted tumors means any target-mediated disposition would be invisible in these data (Supplement, 'Mice PK modeling')."
  )

  ini({
    # Reference subject = FEMALE mouse in the Eisenmann 2020 study.
    # Supplement "Mice PK modeling" table:
    #   Female / Eisenmann : ka 0.856 (12.5%), V 6.54 (9.1%),      Cl 0.449 (10.4%)
    #   Male   / Eisenmann : ka 0.856 (same),  V 3.54 (cov 14.6%), Cl 0.449 (same)
    #   Unknown / Salem    : ka 1.65 (cov 37%), V 3.56 (cov 21.6%), Cl 0.449 (same)
    lka <- log(0.856); label("Absorption rate constant in the Eisenmann oral study (1/h)")  # Supplement "Mice PK modeling" table, Female/Eisenmann row (12.5% RSE)
    lvc <- log(6.54);  label("Volume of distribution in female mice (L/kg)")                # Supplement "Mice PK modeling" table, Female row (9.1% RSE)
    lcl <- log(0.449); label("Clearance (L/h/kg)")                                          # Supplement "Mice PK modeling" table (10.4% RSE); the header prints "L/h", see vignette Errata

    # Covariate effects, back-calculated from the reported per-group values as
    # multiplicative factors on the reference (the Supplement prints the
    # resulting parameter per group rather than the coefficient itself).
    e_male_vc  <- fixed(0.5413); label("Ratio of male to female volume of distribution (unitless)")           # 3.54 / 6.54 from the Supplement table (male effect reported with 14.6% RSE)
    e_salem_vc <- fixed(0.5443); label("Ratio of Salem-study to female-reference volume of distribution (unitless)")  # 3.56 / 6.54 from the Supplement table (21.6% RSE)
    e_salem_ka <- fixed(1.9276); label("Ratio of Salem-study to Eisenmann-study absorption rate constant (unitless)") # 1.65 / 0.856 from the Supplement table (37% RSE)

    # Residual error. The Supplement reports the additive and proportional
    # magnitudes without stating variability was estimated on top of them, and
    # reports no IIV at all -- see vignette Errata.
    addSd  <- 0.08536; label("Additive residual SD (mg/L)")                # Supplement "Mice PK modeling": "Additive error = 0.08536 (37%)"
    propSd <- 0.01267; label("Proportional residual SD (fraction)")        # Supplement "Mice PK modeling": "Proportional error = 0.01267 (RSE = 315.8%)"
  })

  model({
    # The sex term is switched off for the Salem cohort, whose own factor
    # replaces it, so the three published rows are reproduced exactly.
    ka <- exp(lka) * e_salem_ka^STUDY_SALEM
    vc <- exp(lvc) * e_male_vc^((1 - SEXF) * (1 - STUDY_SALEM)) * e_salem_vc^STUDY_SALEM
    cl <- exp(lcl)
    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
