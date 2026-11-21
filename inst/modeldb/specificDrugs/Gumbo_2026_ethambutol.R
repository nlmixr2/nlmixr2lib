Gumbo_2026_ethambutol <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order oral absorption and",
    "linear elimination for ethambutol in adults with Mycobacterium kansasii",
    "complex lung disease (Gumbo 2026). This is the plasma (blood-compartment)",
    "PK prior the authors fed to their ADAPT 5 Monte Carlo dose-finding",
    "experiments: clearance fixed at the value for the average 64 kg patient",
    "with normal renal function, no covariates, negligible protein binding, and",
    "log-normal between-subject variability on all five parameters. Gumbo 2026",
    "linked this blood compartment to the centre of lung cavities through a",
    "dynamical sink model whose ethambutol parameters are published elsewhere",
    "and are not part of this file, so the model predicts plasma, not",
    "cavity-centre, concentrations. The companion model",
    "Gumbo_2026_ethambutol_hfs carries the hollow-fiber exposure-response fits",
    "that define the AUC0-24/MIC target."
  )
  reference <- paste(
    "Gumbo T, Boorgula GD, Srivastava S.",
    "Preclinical pharmacokinetics/pharmacodynamics studies defining the role of",
    "ethambutol in Mycobacterium kansasii lung disease. Antimicrob Agents",
    "Chemother. 2026;70(2):e01447-25. doi:10.1128/aac.01447-25.",
    "PK parameters from Materials and Methods, 'Monte Carlo experiments'."
  )
  vignette <- "Gumbo_2026_ethambutol"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "ethambutol", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ethambutol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ethambutol", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Gumbo 2026 Monte Carlo experiments: 'We did not examine for the two",
        "important covariates of weight (including weight banding) and renal",
        "function, but used the clearance for the average 64 kg weight in",
        "patients with MKC lung disease from Poland with normal renal",
        "function'. The 39.7 L/h clearance is therefore the value at 64 kg and",
        "no weight relationship is encoded."
      ),
      source_name = "weight"
    ),
    CRCL = list(
      description = "Creatinine clearance (renal function)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Named by Gumbo 2026 as an important covariate that was deliberately not examined; normal renal function was assumed.",
      source_name = "renal function"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 10000L,
    n_studies = NA_integer_,
    disease_state = "Mycobacterium kansasii complex (MKC) lung disease (virtual patients)",
    weight_range = "64 kg (single typical weight; average weight of MKC lung-disease patients in Poland)",
    renal_function = "normal",
    dose_range = "Oral ethambutol 600, 1,200, 1,800, 2,400 and 3,000 mg once daily",
    regions = "Virtual population; clearance anchored to Polish MKC lung-disease patients",
    notes = paste(
      "Not a fitted cohort. Gumbo 2026 assembled these parameters as the",
      "'subroutine prior' of ADAPT 5 Monte Carlo experiments (10,000 virtual",
      "patients per dose) from previously published ethambutol population PK",
      "and lung-cavity penetration studies (Gumbo 2026 references 11, 18 and",
      "19: Chigutsa 2015, Hall 2012, Dheda 2018). The values are not those of",
      "Hall 2012 evaluated at 64 kg (Hall 2012 gives CL = 42.6 * (WT/45.6)^0.75,",
      "i.e. 54.9 L/h at 64 kg, and larger volumes). Protein binding was assumed",
      "negligible, so total and unbound concentrations coincide."
    )
  )

  ini({
    # Gumbo 2026 Materials and Methods, 'Monte Carlo experiments': 'We assumed
    # a two-compartment model with the following parameters and variances (%)'.
    lka <- log(0.498)
    label("First-order absorption rate constant Ka (1/h)") # Gumbo 2026 Methods, 'absorption constant (Ka) of 0.498 per hour (IIV = 43%)'
    lcl <- log(39.7)
    label("Apparent total clearance CL at 64 kg, normal renal function (L/h)") # Gumbo 2026 Methods, 'total clearance of 39.7 L/h (inter-individual variance [IIV] = 34%)'
    lvc <- log(90)
    label("Apparent central volume of distribution Vc (L)") # Gumbo 2026 Methods, 'central volume of 90 L (IIV = 40%)'
    lq <- log(33.8)
    label("Apparent intercompartmental clearance Q (L/h)") # Gumbo 2026 Methods, 'intercompartmental clearance 33.8 L/h (IIV = 40%)'
    lvp <- log(642)
    label("Apparent peripheral volume of distribution Vp (L)") # Gumbo 2026 Methods, 'peripheral volume of 642 L (IIV = 40%)'

    # Between-subject variability. The paper prints each IIV as a percentage;
    # it is read as a coefficient of variation of a log-normal distribution,
    # omega^2 = log(1 + CV^2). See the vignette Assumptions and deviations.
    etalka ~ 0.16968 # Gumbo 2026 Methods, Ka IIV = 43%: log(1 + 0.43^2) = 0.16968
    etalcl ~ 0.10940 # Gumbo 2026 Methods, CL IIV = 34%: log(1 + 0.34^2) = 0.10940
    etalvc ~ 0.14842 # Gumbo 2026 Methods, Vc IIV = 40%: log(1 + 0.40^2) = 0.14842
    etalq ~ 0.14842 # Gumbo 2026 Methods, Q IIV = 40%: log(1 + 0.40^2) = 0.14842
    etalvp ~ 0.14842 # Gumbo 2026 Methods, Vp IIV = 40%: log(1 + 0.40^2) = 0.14842

    # The Monte Carlo experiments simulate true concentrations; no residual
    # error was reported, so it is held at zero.
    propSd <- fixed(0)
    label("Proportional residual SD (0; not reported in Gumbo 2026)") # Gumbo 2026 reports no residual error for the Monte Carlo PK prior
  })

  model({
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # mg / L = mg/L. Unbound equals total (negligible protein binding assumed).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
