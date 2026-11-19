Stepniewska_2009_artesunate <- function() {
  description <- paste(
    "One-compartment population PK model of dihydroartemisinin (DHA) after a",
    "3-day oral artesunate-amodiaquine course (fixed-dose combination or loose",
    "tablets) in African children aged 6 months to 5 years with uncomplicated",
    "falciparum malaria (Stepniewska 2009). Artesunate is assumed to convert",
    "completely to DHA (molar conversion); first-order absorption; per-kilogram",
    "apparent clearance with an additive increase after the first dose (acute",
    "illness versus recovery); per-kilogram apparent volume with an additive",
    "linear age effect and IIV; additive residual error."
  )
  reference <- paste(
    "Stepniewska K, Taylor W, Sirima SB, Ouedraogo EB, Ouedraogo A, Gansane A,",
    "Simpson JA, Morgan CC, White NJ, Kiechel JR.",
    "Population pharmacokinetics of artesunate and amodiaquine in African children.",
    "Malar J. 2009;8(Suppl 1):200. doi:10.1186/1475-2875-8-200.",
    sep = " "
  )
  vignette <- "Stepniewska_2009_artesunate_amodiaquine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # The depot receives the administered artesunate dose (mg). The paper assumed
  # complete conversion of artesunate to DHA 'for determination of the dose'
  # (Methods); the published AUCs (4.29 ug*h/mL after a 3.7 mg/kg artesunate
  # dose, Results) equal 3.7 * 284.9 / 384.4 / CL/F, so the conversion is molar.
  # It is applied at the depot-to-central transfer, so central holds DHA mass.
  compartmentData <- list(
    depot = list(analyte = "artesunate", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "dihydroartemisinin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "CL/F and V/F are reported per kilogram (Table 4), so each is multiplied by",
        "WT to give absolute values for a dose in mg. Weight also explained part of",
        "the variation in V/F but was not independent of age; the final model",
        "carries the age effect (Results, 'Artesunate and DHA')."
      ),
      source_name = "weight"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The paper records age in months (Table 1). The model converts AGE (years)",
        "to months and applies the Table 4 'Effect of Age on V/F' (0.063) as an",
        "additive change in V/F (L/kg) per month, centred at the cohort median age",
        "of 32 months: V/F = 2.285 + 0.063 * (AGE_months - 32). The paper prints",
        "neither the equation nor the centring value. This form is the one that",
        "reproduces the published half-life distribution: median 2.5 h equals",
        "log(2) * 2.285 / 0.636, and the 90% range 0.8-4.5 h is additively",
        "symmetric about it, as an age term spanning 7-59 months produces (V/F",
        "0.71-4.0 L/kg). An exponential form per month would span 27-fold and a",
        "proportional linear form per month turns negative at 7 months."
      ),
      source_name = "age"
    ),
    MULTI_DOSE = list(
      description = "Sampling dosing period (0 = after the first dose, 1 = after a later dose)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (first dose)",
      notes = paste(
        "The paper's binary covariate is 0 for samples after dose 1 and 1 for",
        "samples after dose 3; no samples were taken after dose 2, so the value",
        "for the second dosing interval is not informed by the data (this model",
        "treats it as 1). The effect is additive on per-kilogram CL/F: 0.636 +",
        "0.760 = 1.396 L/kg/h. The additive form is confirmed by the published",
        "third-dose AUC (1.96 ug*h/mL = 2.742 / 1.396) and half-life (1.1 h).",
        "The paper attributes the change to recovery from acute illness."
      ),
      source_name = "dosing period"
    )
  )

  covariatesDataExcluded <- list(
    BODYTEMP = list(
      description = "Body temperature at enrolment",
      units = "degC",
      type = "continuous",
      notes = "Screened; not retained."
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Screened; not retained."
    ),
    LNPC = list(
      description = "Natural log of enrolment parasitaemia",
      units = "log(parasites/uL)",
      type = "continuous",
      notes = "Screened; not retained."
    ),
    RESP_RATE = list(
      description = "Respiratory rate at enrolment",
      units = "breaths/min",
      type = "continuous",
      notes = "Screened; not retained."
    ),
    GAMETOCYTE_POS = list(
      description = "Gametocytaemia present at enrolment (1 = yes)",
      units = "(binary)",
      type = "binary",
      notes = "Screened; not retained."
    ),
    FORM_FDC = list(
      description = "Fixed-dose co-formulation (1) versus loose tablets (0)",
      units = "(binary)",
      type = "binary",
      notes = "Screened; no difference in CL/F between arms (p = 0.73). Not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 70L,
    n_subjects_pk = "70 children in the artesunate PK group (35 per arm); 4 excluded for vomiting. Many samples were not reportable because of haemolysis.",
    n_studies = 1L,
    age_range = "7-59 months (artesunate PK group, Table 1)",
    age_median = "37 months (loose tablets) and 32 months (fixed-dose combination)",
    weight_range = "7-19 kg",
    weight_median = "14 kg (loose tablets) and 12 kg (fixed-dose combination)",
    sex_female_pct = 42.9,
    race_ethnicity = "Not reported by race; Burkinabe children",
    disease_state = "Acute uncomplicated Plasmodium falciparum malaria (>1,000 parasites/uL, axillary temperature >= 37.5 degC)",
    dose_range = paste(
      "Artesunate once daily for 3 days, target 4 mg/kg/day (therapeutic range",
      "2-10): fixed-dose combination 25 mg AS per tablet, 1 tablet (<12 months) or",
      "2 tablets (12-60 months); loose Arsumax 50 mg tablets, half or one tablet",
      "by the same age bands."
    ),
    sampling = "Pre-dose 1, then one sample at 0.5, 1, 2 or 4 h (randomised) after dose 1 and after dose 3.",
    regions = "Burkina Faso (Pissy health district, Ouagadougou)",
    notes = paste(
      "Demographics from Stepniewska 2009 Table 1 (artesunate PK group);",
      "sex_female_pct derived from 25/35 and 15/35 males. Analysis in S-PLUS 6.0,",
      "following the approach of Simpson 2006 (see modellib('Simpson_2006_artesunate'))."
    )
  )

  ini({
    lka <- log(4.271)
    label("First-order absorption rate constant (1/h)") # Table 4 DHA: Ka = 4.271 (SE 2.091)
    lcl <- log(0.636)
    label("Apparent DHA clearance per kg after the first dose (CL/F, L/kg/h)") # Table 4 DHA: CL/F = 0.636 (SE 0.075), first dose period
    e_multi_dose_cl <- 0.760
    label("Additive increase in CL/F after later doses (L/kg/h)") # Table 4 DHA: 'Increase in CL/F due to the 3rd dose period' = 0.760 (SE 0.160)
    lvc <- log(2.285)
    label("Apparent DHA volume per kg at age 32 months (V/F, L/kg)") # Table 4 DHA: V/F = 2.285 (SE 0.317)
    e_age_vc <- 0.063
    label("Additive change in V/F per month of age (L/kg per month)") # Table 4 DHA: 'Effect of Age on V/F' = 0.063 (SE 0.015)

    # IIV is log-normal and the reported CV% is the square root of the variance
    # estimate (Methods), so the variance is 0.47^2. Clearance had essentially no
    # between-subject variation (SD < 1e-7, Results), so only V/F carries an eta.
    # Source: Table 4 DHA, V/F inter-subject variation 47%.
    etalvc ~ 0.47^2

    addSd <- 340
    label("Additive residual SD (ng/mL)") # Table 4 DHA: sigma = 340 ng/mL
  })

  model({
    # Molecular weights (g/mol) from Methods ('taking molecular weight as 384.4
    # for artesunate and 284.9 for DHA').
    mw_ars <- 384.4
    mw_dihydroart <- 284.9

    # Per-kilogram CL/F and V/F with additive covariate effects. Age enters in
    # months, centred at the median age of 32 months.
    age_months <- AGE * 12
    clkg <- exp(lcl) + e_multi_dose_cl * MULTI_DOSE
    vckg <- (exp(lvc) + e_age_vc * (age_months - 32)) * exp(etalvc)

    ka <- exp(lka)
    cl <- clkg * WT
    vc <- vckg * WT
    kel <- cl / vc

    # Artesunate (mg) in depot converts mole-for-mole to DHA (mg) on absorption.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot * (mw_dihydroart / mw_ars) - kel * central

    # central (mg) / vc (L) is mg/L; x1000 gives ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd)
  })
}
