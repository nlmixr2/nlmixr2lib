NaBangchang_2007_fosmidomycin_monotherapy <- function() {
  description <- paste(
    "One-compartment pharmacokinetic model with first-order absorption, an",
    "absorption lag time and first-order elimination for oral fosmidomycin",
    "given alone (1,200 mg every 8 h for 7 days) to Thai, Karen and Burmese",
    "adults with acute uncomplicated multidrug-resistant Plasmodium",
    "falciparum malaria (Na-Bangchang 2007, monotherapy arm, n = 15). The",
    "model was fitted in ADAPT II and only the median ka, k10, tlag and Vz/F",
    "per arm are reported (Table 3), so this is a typical-value model, not a",
    "population model: every parameter is fixed, there is no between-subject",
    "variability and the residual SD is fixed at zero. Clearance is the",
    "derived k10 x Vz/F. The combination arm (fosmidomycin 900 mg with",
    "clindamycin 600 mg every 12 h) is",
    "modellib('NaBangchang_2007_fosmidomycin_combination') and its clindamycin",
    "is modellib('NaBangchang_2007_clindamycin')."
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
    depot = list(analyte = "fosmidomycin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "fosmidomycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 15,
    n_studies = 1,
    age_range = "18-50 years (inclusion criterion)",
    age_median = "36 years (95% CI 19-45)",
    weight_median = "55 kg (95% CI 42-67)",
    sex_female_pct = 26.7,
    race_ethnicity = c(Thai = 33.3, Karen = 13.3, Burmese = 53.3),
    disease_state = paste(
      "Acute uncomplicated Plasmodium falciparum malaria (mono-infection,",
      "parasitaemia 1,000-50,000/uL at screening) in a multidrug-resistant",
      "area on the Thailand-Myanmar border."
    ),
    dose_range = "Fosmidomycin 1,200 mg orally every 8 h for 7 days (400 mg capsules), taken with 250 mL water and a standard hospital meal.",
    regions = "Thailand (Mae Sot General Hospital, Tak Province)",
    notes = paste(
      "Na-Bangchang 2007 Table 1 and Methods. Open-label study conducted",
      "2000-2002; patients hospitalised for 7 days. Plasma sampled at 0, 1, 2,",
      "4, 8, 16, 24, 48, 72, 96, 120, 144, 168, 192 and 216 h from the first",
      "dose; fosmidomycin assayed by agar disk-diffusion bioassay",
      "(Enterobacter cloacae ATCC 23355). The 28-day cure rate in this arm",
      "was 22% (2/9 evaluable)."
    )
  )

  ini({
    # Na-Bangchang 2007 Table 3, 'Fosmidomycin monotherapy (N = 15)' column,
    # 'median (%C.V)'. Parameters of the ADAPT II one-compartment model with
    # first-order absorption, absorption lag time and first-order elimination.
    lka <- fixed(log(1.27)); label("Absorption rate constant K01 (1/h)") # Table 3 'K01 (/hr)' = 1.27
    lkel <- fixed(log(0.20)); label("Elimination rate constant K10 (1/h)") # Table 3 'K10 (/hr)' = 0.20
    lvc <- fixed(log(214)); label("Apparent volume of distribution Vz/F (L)") # Table 3 'Vz/F (l)' = 214.00
    ltlag <- fixed(log(0.5)); label("Absorption lag time tlag (h)") # Table 3 'tlag (hr)' = 0.5

    propSd <- fixed(0); label("Proportional residual SD on fosmidomycin plasma concentration, not reported (fraction)") # ADAPT II fit; no residual-error model reported
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
