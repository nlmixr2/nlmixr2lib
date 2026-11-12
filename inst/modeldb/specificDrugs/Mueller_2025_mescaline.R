Mueller_2025_mescaline <- function() {
  description <- paste(
    "One-compartment pharmacokinetic model of oral mescaline (dosed as",
    "mescaline hydrochloride capsules) in healthy adults (Mueller 2025), with",
    "first-order absorption after an absorption lag time and first-order",
    "elimination. Clearance and volume are apparent (Cl/F and Vd/F) because",
    "only the oral route was studied; urinary recovery puts the oral",
    "bioavailability at 53 percent or more. Mueller 2025 fitted every",
    "participant-condition profile individually in Phoenix WinNonlin (no",
    "population model), so the parameters are the geometric means of 105",
    "individual estimates over all doses (100-800 mg mescaline HCl) and the",
    "model carries no between-subject variability; the proportional residual",
    "error was used but its magnitude was not reported and is fixed at zero.",
    "Doses are mg of mescaline base. The effect-compartment sigmoid-Emax",
    "PK-PD models for the three subjective-effect scales are",
    "Mueller_2025_mescaline_anydrugeffect, _gooddrugeffect and",
    "_egodissolution."
  )
  reference <- paste(
    "Mueller L, Klaiber A, Ley L, Becker AM, Thomann J, Luethi D, Schmid Y,",
    "Liechti ME (2025). Pharmacokinetics, Pharmacodynamics, and Urinary",
    "Recovery of Oral Mescaline Hydrochloride in Healthy Participants.",
    "Clin Pharmacokinet 64(10):1495-1506. doi:10.1007/s40262-025-01544-x"
  )
  vignette <- "Mueller_2025_mescaline"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  compartmentData <- list(
    depot = list(analyte = "mescaline", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "mescaline", units = "mg", specimen = "plasma", verified = TRUE)
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
      "Single oral doses of 100, 200, 300, 400, 500 and 800 mg mescaline HCl,",
      "and 800 mg mescaline HCl with 40 mg ketanserin; the analysed doses are",
      "the analytically measured mescaline-base contents 82.5, 164.9, 241.6,",
      "329.9, 404.3, 658.3 and 658.9 mg."
    ),
    regions = "Switzerland (University Hospital Basel)",
    n_profiles = "105 PK profiles (113 administrations less 8 with early emesis)",
    notes = paste(
      "Supplementary Table S1: 49 participants (25 female) pooled from two",
      "randomised double-blind crossover studies (NCT04227756: 16 participants",
      "on 300 mg and 16 on 500 mg mescaline HCl alongside LSD, psilocybin and",
      "placebo; NCT04849013: 17 participants on 100, 200, 400 and 800 mg and",
      "800 mg + 40 mg ketanserin). 69 percent had prior psychedelic",
      "experience. Plasma was sampled to 24 h (Study 1) or 30 h (Study 2).",
      "Profiles with emesis within 1 h of dosing or of Cmax were excluded.",
      "The ketanserin condition is included in this PK fit but not in the",
      "PK-PD fits."
    )
  )

  ini({
    # Table 2 'All doses' row: geometric mean (95% CI) of 105 individual
    # Phoenix WinNonlin fits. The secondary parameters in the same row follow
    # from these: ke = 43/217 = 0.198 /h, t1/2 = ln(2)/ke = 3.50 h (Table 2:
    # 3.5 h) and tmax = ln(ka/ke)/(ka - ke) + tlag = 1.96 h (Table 2: 2.0 h).
    lka <- fixed(log(1.24)); label("First-order absorption rate constant (1/h)") # Table 2 All doses: ka = 1.24 (1.07-1.44) 1/h
    ltlag <- fixed(log(0.20)); label("Absorption lag time (h)") # Table 2 All doses: tlag = 0.20 (0.17-0.23) h
    lcl <- fixed(log(43)); label("Apparent clearance Cl/F (L/h)") # Table 2 All doses: Cl/F = 43 (41-44) L/h
    lvc <- fixed(log(217)); label("Apparent volume of distribution Vd/F (L)") # Table 2 All doses: Vd/F = 217 (209-226) L

    propSd <- fixed(0); label("Proportional residual error on plasma mescaline (fraction; not reported)") # Methods: multiplicative residual error model used; magnitude not reported
  })

  model({
    ka <- exp(lka)
    tlag <- exp(ltlag)
    cl <- exp(lcl)
    vc <- exp(lvc)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    alag(depot) <- tlag

    # Dose in mg and volume in L give mg/L; x 1000 gives ug/L as reported.
    Cc <- 1000 * central / vc

    Cc ~ prop(propSd)
  })
}
