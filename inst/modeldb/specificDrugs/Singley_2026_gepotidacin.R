Singley_2026_gepotidacin <- function() {
  description <- paste(
    "Preclinical (mouse, neutropenic male CD-1; thigh infection with",
    "Escherichia coli or Klebsiella pneumoniae). One-compartment population",
    "PK model with first-order absorption for subcutaneous gepotidacin in",
    "infected neutropenic mice (Singley 2026). The dose amount is in mg per",
    "mouse and the observation Cc is the UNBOUND (free) plasma",
    "concentration. The paper converts total to free concentrations with",
    "a fixed 24% plasma protein binding and does not say at which step;",
    "the published parameters predict free concentrations directly, since",
    "they reproduce the free-drug fAUC0-24 and fCmax of Table S2 without",
    "the 0.76 unbound-fraction factor and the predictions axis of Figure",
    "S4 spans the free-drug range of Figure 1. The absorption rate constant falls with the",
    "dose per administration as a power function of DOSE / 2.7 mg (2.7 mg",
    "is 100 mg/kg in a 27 g mouse), which the authors added to describe",
    "the slight dose non-proportionality of the observed exposure (Table",
    "1): absorption is faster at low doses, so the peak rises less than",
    "dose-proportionally.",
    "Between-subject variability on CL and ka; additive residual error."
  )
  reference <- paste(
    "Singley CM, Kurumaddali A, Hoover JL. Pharmacokinetic and",
    "pharmacodynamic characterization of gepotidacin against Escherichia",
    "coli and Klebsiella pneumoniae in a neutropenic mouse thigh infection",
    "model. Antimicrob Agents Chemother. 2026;70(1):e01176-25.",
    "doi:10.1128/aac.01176-25. Parameter values from Supplemental Material",
    "Table S1 and the Figure S1 equation."
  )
  vignette <- "Singley_2026_gepotidacin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL (unbound)")

  covariateData <- list(
    DOSE = list(
      description = "Gepotidacin dose given at each subcutaneous administration",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Use case (a) of the DOSE register entry: the dose per",
        "ADMINISTRATION in mg per mouse (the same unit as amt), entering ka",
        "as (DOSE / 2.7)^-0.437 (Figure S1 equation). The 2.7 mg reference",
        "is 100 mg/kg in a 27 g mouse; to convert a mg/kg dose, multiply",
        "by the body weight in kg (the study mice weighed 27-28 g). For an",
        "every-6-h regimen of equal doses it is the size of each dose, not",
        "the daily total. Supply it on every record and place the event",
        "columns before DOSE in the event table."
      ),
      source_name = "Dose (Figure S1 'Individual Ka = Kapop * (Dose/2.7)^Deff * e^eta')"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "gepotidacin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "gepotidacin (unbound)", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "mouse (CD-1, outbred, male; neutropenic)",
    n_subjects = NA_integer_,
    n_studies = NA_integer_,
    weight_range = "27-28 g",
    sex_female_pct = 0,
    organism = "Escherichia coli and Klebsiella pneumoniae thigh infections",
    disease_state = paste(
      "Neutropenic (cyclophosphamide 150 mg/kg on day -4 and 100 mg/kg on",
      "day -1, intraperitoneally) mice infected intramuscularly in one",
      "thigh, dosed 1 h after infection"
    ),
    dose_range = "Single subcutaneous dose of 6.25-800 mg/kg (PK studies)",
    regions = "United States (GSK)",
    notes = paste(
      "Terminal cardiac-puncture plasma samples at 0.25, 0.5, 1, 2, 4 and 6",
      "h post-dose, N = 3 mice per study, dose and time point, so each",
      "mouse contributed one sample; data pooled across studies (the",
      "number of mice is not reported). Gepotidacin measured by",
      "HPLC-MS/MS, LLOQ 0.05 ug/mL. Fitted in NONMEM 7.3. Plasma protein",
      "binding 24% in the mouse (Bulik 2017), concentration-independent.",
      "High epsilon-shrinkage is reported because of the one-sample-per-",
      "mouse design."
    )
  )

  ini({
    # Table S1 'Final gepotidacin population PK model parameter estimates
    # following subcutaneous administration to mice'. Dose in mg per mouse,
    # CL in L/h and V in L per mouse.
    lka <- log(0.886)
    label("Absorption rate constant ka at the 2.7 mg reference dose (1/h)") # Table S1 'Kapop (1/h)' = 0.886 (RSE 4%)
    lcl <- log(0.101)
    label("Clearance CL (L/h)") # Table S1 'CLpop (L/h)' = 0.101 (RSE 4%)
    lvc <- log(0.018)
    label("Central volume of distribution V (L)") # Table S1 'Vpop (L)' = 0.018 (RSE 13%)

    # Figure S1: Individual Ka = Kapop * (Dose/2.7)^Deff * e^eta.
    e_dose_ka <- -0.437
    label("Power exponent of dose per administration (DOSE/2.7 mg) on ka (unitless)") # Table S1 'Dose effect on Kapop' = -0.437 (RSE 8%); Figure S1 'Deff = -0.437'

    # Exponential between-subject variability. Table S1 gives %CV with
    # %CV = sqrt(exp(variance) - 1) * 100, so variance = log(CV^2 + 1).
    etalcl ~ 0.14090 # Table S1 'IIV on CL' = 38.9%; log(0.389^2 + 1)
    etalka ~ 0.04039 # Table S1 'IIV on Ka' = 20.3%; log(0.203^2 + 1); Figure S1 'eta ~ N(0, 0.04)'

    addSd <- 0.128
    label("Additive residual SD (ug/mL)") # Table S1 'Add err (mg/L)' = 0.128 (RSE 27%)
  })

  model({
    ka <- exp(lka + etalka) * (DOSE / 2.7)^e_dose_ka
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Unbound plasma concentration (mg/L = ug/mL).
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
