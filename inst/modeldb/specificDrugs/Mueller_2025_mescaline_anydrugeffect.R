Mueller_2025_mescaline_anydrugeffect <- function() {
  description <- paste(
    "Effect-compartment sigmoid-Emax PK-PD model of the any drug effect visual",
    "analog scale (VAS, 0-100 percent) after oral mescaline (dosed as",
    "mescaline hydrochloride capsules) in healthy adults (Mueller 2025).",
    "Mescaline PK is the one-compartment model with first-order absorption,",
    "an absorption lag time and first-order elimination of",
    "Mueller_2025_mescaline, held fixed; the effect-site concentration",
    "follows plasma with the first-order equilibration rate constant ke0, and",
    "the effect is 100 * Ce^gamma / (EC50^gamma + Ce^gamma) with Emax fixed at",
    "100 percent, the top of the scale, and zero baseline.",
    "'Any drug effect' was the primary PD outcome; the same structure was",
    "applied to the exploratory 'good drug effect' and 'ego dissolution'",
    "scales (Mueller_2025_mescaline_gooddrugeffect, _egodissolution).",
    "Mueller 2025",
    "fitted each participant-condition profile individually in Phoenix",
    "WinNonlin with the individual PK frozen (no population model); the PD",
    "parameters are the arithmetic means of 82 individual estimates over all",
    "doses (100-800 mg mescaline HCl; profiles with a predicted maximal",
    "effect below 1 percent excluded) and the PK parameters the geometric",
    "means of the 105 individual PK estimates. The model therefore carries no",
    "between-subject variability; the PK proportional and PD additive",
    "residual errors were used but their magnitudes were not reported and are",
    "fixed at zero. Doses are mg of mescaline base."
  )
  reference <- paste(
    "Mueller L, Klaiber A, Ley L, Becker AM, Thomann J, Luethi D, Schmid Y,",
    "Liechti ME (2025). Pharmacokinetics, Pharmacodynamics, and Urinary",
    "Recovery of Oral Mescaline Hydrochloride in Healthy Participants.",
    "Clin Pharmacokinet 64(10):1495-1506. doi:10.1007/s40262-025-01544-x"
  )
  vignette <- "Mueller_2025_mescaline"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L", anydrugeffectvas = "percent of VAS (0-100)")

  compartmentData <- list(
    depot = list(analyte = "mescaline", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "mescaline", units = "mg", specimen = "plasma", verified = TRUE),
    effect = list(analyte = "mescaline", units = "ug/L", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = 49L,
    n_studies = 2L,
    age_range = "25-65 years (inclusion criterion)",
    age_mean = "30 +/- 7 years",
    weight_mean = "70 +/- 9 kg",
    bmi_mean = "22 +/- 2 kg/m^2 (inclusion 18-29 kg/m^2)",
    sex_female_pct = 51,
    disease_state = "healthy volunteers",
    dose_range = paste(
      "Single oral doses of 100, 200, 300, 400, 500 and 800 mg mescaline HCl;",
      "the analysed doses are the analytically measured mescaline-base",
      "contents 82.5, 164.9, 241.6, 329.9, 404.3 and 658.3 mg. The 800 mg +",
      "40 mg ketanserin condition was excluded from the PK-PD analysis."
    ),
    regions = "Switzerland (University Hospital Basel)",
    n_profiles = "91 PK-PD profiles; 82 contributed to the summary of the primary PD parameters (predicted maximal effect of 1 percent or more)",
    notes = paste(
      "Supplementary Table S1: 49 participants (25 female) pooled from two",
      "randomised double-blind crossover studies (NCT04227756: 300 or 500 mg",
      "mescaline HCl; NCT04849013: 100, 200, 400 and 800 mg). 69 percent had",
      "prior psychedelic experience. The VAS was rated at the plasma sampling",
      "times, to 24 h (Study 1) or 30 h (Study 2). Profiles with emesis within",
      "1 h of dosing or of Cmax were excluded."
    )
  )

  ini({
    # PK layer: Table 2 'All doses' geometric means of the 105 individual PK
    # fits, as in Mueller_2025_mescaline. Individual PK estimates were frozen
    # for the PK-PD fits (Methods 'PK Analyses and PK-PD Modeling').
    lka <- fixed(log(1.24)); label("First-order absorption rate constant (1/h)") # Table 2 All doses: ka = 1.24 (1.07-1.44) 1/h
    ltlag <- fixed(log(0.20)); label("Absorption lag time (h)") # Table 2 All doses: tlag = 0.20 (0.17-0.23) h
    lcl <- fixed(log(43)); label("Apparent clearance Cl/F (L/h)") # Table 2 All doses: Cl/F = 43 (41-44) L/h
    lvc <- fixed(log(217)); label("Apparent volume of distribution Vd/F (L)") # Table 2 All doses: Vd/F = 217 (209-226) L

    # PD layer: Table 3, VAS 'any drug effect' block, 'All doses' row: arithmetic
    # mean (95% CI) of 82 individual estimates.
    emax <- fixed(100); label("Maximal effect on the VAS (percent)") # Methods: sigmoid Emax model 'with an Emax of 100%, according to the range of VASs from 0 to 100%'
    lec50 <- fixed(log(673)); label("Effect-site concentration giving half-maximal effect EC50 (ug/L)") # Table 3 any drug effect All doses: EC50 = 673 (605-741) ug/L
    lhill <- fixed(log(3.5)); label("Sigmoidicity (Hill) coefficient gamma (unitless)") # Table 3 any drug effect All doses: gamma = 3.5 (3.1-3.9)
    lke0 <- fixed(log(1.7)); label("Effect-compartment equilibration rate constant ke0 (1/h)") # Table 3 any drug effect All doses: ke0 = 1.7 (1.4-2.0) 1/h

    propSd <- fixed(0); label("Proportional residual error on plasma mescaline (fraction; not reported)") # Methods: multiplicative residual error model used; magnitude not reported
    addSd_anydrugeffectvas <- fixed(0); label("Additive residual error on the any drug effect VAS (percent; not reported)") # Methods: additive residual error model used for the PD effects; magnitude not reported
  })

  model({
    ka <- exp(lka)
    tlag <- exp(ltlag)
    cl <- exp(lcl)
    vc <- exp(lvc)
    ec50 <- exp(lec50)
    hill <- exp(lhill)
    ke0 <- exp(lke0)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    alag(depot) <- tlag

    # Dose in mg and volume in L give mg/L; x 1000 gives ug/L as reported.
    Cc <- 1000 * central / vc

    # Effect site (ug/L): dCe/dt = ke0 * (Cp - Ce).
    d/dt(effect) <- ke0 * (Cc - effect)

    anydrugeffectvas <- emax * effect^hill / (ec50^hill + effect^hill)

    Cc ~ prop(propSd)
    anydrugeffectvas ~ add(addSd_anydrugeffectvas)
  })
}
