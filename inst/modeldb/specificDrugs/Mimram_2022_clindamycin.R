Mimram_2022_clindamycin <- function() {
  description <- "One-compartment population PK model with first-order absorption and elimination for orally administered clindamycin in adults with chronic prosthetic joint infections not receiving rifampicin; no covariates retained (Mimram 2022)."
  reference <- "Mimram L, Magreault S, Kerroumi Y, Salmon D, Kably B, Marmor S, Jannot A-S, Jullien V, Zeller V. Population Pharmacokinetics of Orally Administered Clindamycin to Treat Prosthetic Joint Infections: A Prospective Study. Antibiotics (Basel). 2022;11(11):1462. doi:10.3390/antibiotics11111462"
  vignette <- "Mimram_2022_clindamycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "clindamycin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "clindamycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  # Covariates screened by forward inclusion (Section 4.7). None changed the
  # objective function significantly, so the final model has none (Section 2.2).
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Tested as a categorical covariate (Section 4.7); not retained (Section 2.2).",
      source_name = "sex"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested as a continuous covariate normalised to its median (Section 4.7); not retained (Section 2.2).",
      source_name = "age"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested as a continuous covariate normalised to its median (Section 4.7); not retained (Section 2.2). Patients outside 50-100 kg were excluded.",
      source_name = "weight"
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested as a continuous covariate normalised to its median (Section 4.7); not retained (Section 2.2).",
      source_name = "height"
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested as a continuous covariate normalised to its median (Section 4.7); not retained (Section 2.2).",
      source_name = "bilirubin"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested as a continuous covariate normalised to its median (Section 4.7); not retained (Section 2.2).",
      source_name = "creatinine"
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested as a continuous covariate normalised to its median (Section 4.7); not retained (Section 2.2).",
      source_name = "alanine aminotransferase"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 20L,
    n_studies = 1L,
    n_observations = 140L,
    age_range = "adults >= 18 years; mean 66.5 (SD 15.9) years",
    weight_range = "50-100 kg by exclusion criteria; mean 76.3 (SD 14.4) kg",
    sex_female_pct = 35,
    race_ethnicity = "Not reported",
    disease_state = "Chronic prosthetic joint infection (14 hip, 3 knee, 3 shoulder arthroplasties) caused by clindamycin-susceptible staphylococci, streptococci or anaerobes; hepatocellular insufficiency, cirrhosis, creatinine clearance < 30 mL/min, severe sepsis and CYP3A4/5 inducers or inhibitors excluded.",
    dose_range = "Oral clindamycin 750 mg q8h (< 80 kg) or 900 mg q8h (>= 80 kg), up to 1200 mg q8h (> 95 kg), after 2-5 weeks of IV clindamycin",
    regions = "France (Paris; Diaconesses-Croix Saint-Simon and Cochin hospitals)",
    notes = "Patients not prescribed rifampicin from the prospective NCT02629770 study (December 2015 to November 2019). Seven plasma samples per patient over one steady-state dosing interval (0, 0.5, 1, 2, 4, 6, 8 h) after 2 weeks of oral dosing; concentrations 0.4-13.5 mg/L; LC-MS LLOQ 0.09 mg/L. Demographics from Table 1."
  )

  ini({
    lka <- log(3.53); label("Absorption rate constant (1/h)") # Table 2, Ka = 3.53 /h
    lcl <- log(23); label("Apparent clearance CL/F (L/h)") # Table 2, CL/F = 23.00 L/h
    lvc <- log(103); label("Apparent volume of distribution V/F (L)") # Table 2, V/F = 103.00 L

    # Table 2 prints omega^2 to two decimals (0.14, 0.08, 0.60); the Abstract
    # reports the same values x100 to three figures (14.4%, 8.2%, 59.6%).
    etalcl ~ 0.144 # Table 2 omega2 CL/F = 0.14; Abstract 14.4%
    etalvc ~ 0.082 # Table 2 omega2 V/F = 0.08; Abstract 8.2%
    etalka ~ 0.596 # Table 2 omega2 Ka = 0.60; Abstract 59.6%

    # Table 2 residual-error rows are NONMEM $SIGMA variances: 0.00976 and 0.0801.
    propSd <- 0.0988; label("Proportional residual error (fraction)") # sqrt(0.00976), Table 2 'Proportional error'
    addSd <- 0.283; label("Additive residual error (mg/L)") # sqrt(0.0801), Table 2 'Additive error'
  })
  model({
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
