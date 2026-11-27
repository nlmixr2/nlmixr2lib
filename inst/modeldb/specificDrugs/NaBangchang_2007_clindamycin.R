NaBangchang_2007_clindamycin <- function() {
  description <- paste(
    "One-compartment pharmacokinetic model with first-order absorption, an",
    "absorption lag time and first-order elimination for oral clindamycin",
    "(600 mg every 12 h for 7 days) given together with oral fosmidomycin",
    "(900 mg every 12 h) to Thai, Karen and Burmese adults with acute",
    "uncomplicated multidrug-resistant Plasmodium falciparum malaria",
    "(Na-Bangchang 2007, combination arm, n = 18). The model was fitted in",
    "ADAPT II and only the median ka, k10, tlag and Vz/F are reported",
    "(Table 3), so this is a typical-value model, not a population model:",
    "every parameter is fixed, there is no between-subject variability and",
    "the residual SD is fixed at zero. Clearance is the derived k10 x Vz/F.",
    "Fosmidomycin from the same patients is",
    "modellib('NaBangchang_2007_fosmidomycin_combination')."
  )
  reference <- paste(
    "Na-Bangchang K, Ruengweerayut R, Karbwang J, Chauemung A, Hutchinson D",
    "(2007). Pharmacokinetics and pharmacodynamics of fosmidomycin",
    "monotherapy and combination therapy with clindamycin in the treatment of",
    "multidrug resistant falciparum malaria. Malaria Journal 6:70.",
    "doi:10.1186/1475-2875-6-70"
  )
  vignette <- "NaBangchang_2007_fosmidomycin_clindamycin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    depot = list(analyte = "clindamycin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "clindamycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 18,
    n_studies = 1,
    age_range = "18-50 years (inclusion criterion)",
    age_median = "29 years (95% CI 18-48)",
    weight_median = "56 kg (95% CI 40-75)",
    sex_female_pct = 27.8,
    race_ethnicity = c(Thai = 55.6, Karen = 16.7, Burmese = 27.8),
    disease_state = paste(
      "Acute uncomplicated Plasmodium falciparum malaria (mono-infection,",
      "parasitaemia 1,000-50,000/uL at screening) in a multidrug-resistant",
      "area on the Thailand-Myanmar border."
    ),
    dose_range = "Clindamycin 600 mg orally every 12 h for 7 days (150 mg film-coated tablets) with fosmidomycin 900 mg every 12 h, taken with 250 mL water and a standard hospital meal.",
    regions = "Thailand (Mae Sot General Hospital, Tak Province)",
    co_medication = "Fosmidomycin 900 mg orally every 12 h for 7 days (part of the regimen).",
    notes = paste(
      "Na-Bangchang 2007 Table 1 and Methods. Open-label study conducted",
      "2000-2002; patients hospitalised for 7 days. Plasma sampled at 0, 24,",
      "48, 156, 157, 158, 160, 162, 164, 168, 174 and 180 h from the first",
      "dose; clindamycin assayed by HPLC-UV (parent only; the N-desmethyl",
      "and sulphoxide metabolites were not measured)."
    )
  )

  ini({
    # Na-Bangchang 2007 Table 3, 'Fosmidomycin + Clindamycin (N = 18)' /
    # 'Clindamycin' column, 'median (%C.V)'. Parameters of the ADAPT II
    # one-compartment model with first-order absorption, absorption lag time
    # and first-order elimination.
    lka <- fixed(log(3.99)); label("Absorption rate constant K01 (1/h)") # Table 3 'K01 (/hr)' = 3.99
    lkel <- fixed(log(0.22)); label("Elimination rate constant K10 (1/h)") # Table 3 'K10 (/hr)' = 0.22
    lvc <- fixed(log(84)); label("Apparent volume of distribution Vz/F (L)") # Table 3 'Vz/F (l)' = 84
    ltlag <- fixed(log(0.38)); label("Absorption lag time tlag (h)") # Table 3 'tlag (hr)' = 0.38

    propSd <- fixed(0); label("Proportional residual SD on clindamycin plasma concentration, not reported (fraction)") # ADAPT II fit; no residual-error model reported
  })

  model({
    ka <- exp(lka)
    kel <- exp(lkel)
    vc <- exp(lvc)
    tlag <- exp(ltlag)
    cl <- kel * vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    alag(depot) <- tlag

    Cc <- central / vc

    Cc ~ prop(propSd)
  })
}
