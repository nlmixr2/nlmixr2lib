Grasela_1985_phenobarbital <- function() {
  description <- "One-compartment population PK model for phenobarbital in preterm neonates (Grasela & Donn 1985), derived from routine clinical data via NONMEM."
  reference <- paste(
    "Grasela TH Jr, Donn SM (1985).",
    "Neonatal population pharmacokinetics of phenobarbital derived from",
    "routine clinical data.",
    "Developmental Pharmacology and Therapeutics 8(6):374-383.",
    "doi:10.1159/000457062. PMID:4075936.",
    sep = " "
  )
  vignette <- "Grasela_1985_phenobarbital"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Linear (not allometric) per-kg scaling of CL and V:",
        "CL = 0.0047 L/h/kg * WT and V = 0.96 L/kg * WT.",
        "The abstract reports per-kg typical values, so the implied",
        "weight exponent is 1 (not 0.75); do not retrofit allometric",
        "scaling unless re-fitting on the original data."
      ),
      source_name        = "WT"
    ),
    ASPHYXIA = list(
      description        = "Perinatal asphyxia indicator (1 = 5-minute Apgar score < 5; 0 otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (5-minute Apgar score >= 5)",
      notes              = paste(
        "Time-fixed per subject. Linear-deviation effect on V:",
        "V = V_typical * (1 + 0.13 * ASPHYXIA), i.e., +13% V when the",
        "5-minute Apgar score is < 5. The abstract reports no detected",
        "effect of ASPHYXIA on CL."
      ),
      source_name        = "ASPHYXIA"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 59,
    n_studies      = 1,
    age_range      = "Preterm neonates; gestational age 24-42 weeks (mean 31 weeks). Postnatal age range not reported in the abstract.",
    weight_range   = "Birth weight 0.600-3.620 kg (mean 1.520 kg). Time-varying current weight range not reported in the abstract.",
    sex_female_pct = "Not reported in the abstract.",
    race_ethnicity = "Not reported in the abstract.",
    disease_state  = "Preterm neonates receiving phenobarbital in routine clinical care (typical indication: prevention or treatment of neonatal seizures). The original study was conducted in a single neonatal intensive care unit and analysed retrospectively.",
    dose_range     = "Not reported in the abstract. Real-world neonatal phenobarbital regimens are typically a 15-20 mg/kg IV loading dose followed by 3-5 mg/kg/day maintenance.",
    regions        = "United States (Grasela TH Jr & Donn SM, University of Michigan Medical Center).",
    notes          = paste(
      "Population description is reconstructed from the Grasela & Donn 1985",
      "abstract (PMID:4075936). The full publication was not consulted",
      "(abstract-only source confirmed by the operator at extraction time);",
      "the mean serum half-life of 141 h reported in the abstract is",
      "recovered as a derived check from the typical CL and V."
    )
  )

  ini({
    # Structural typical values (Grasela & Donn 1985 abstract).
    # CL and V are the WT-NORMALIZED per-kg values; the model multiplies by
    # the subject's body weight (linear exponent, NOT 0.75 allometry) inside
    # model(). At WT = 1.52 kg (cohort mean) and ASPHYXIA = 0 the derived
    # half-life is log(2) / (0.0047 / 0.96) = 141.6 h, matching the
    # abstract's reported mean serum half-life of 141 h.
    lcl <- log(0.0047); label("Weight-normalized clearance (L/h/kg)")          # abstract Results
    lvc <- log(0.96);   label("Weight-normalized volume of distribution (L/kg)") # abstract Results

    # Covariate effect: 5-minute Apgar < 5 increases V by 13% (p < 0.05).
    e_asphyxia_vc <- 0.13; label("Effect of perinatal asphyxia (5-min Apgar < 5) on V (fraction)")  # abstract Results

    # Inter-individual variability. Abstract reports CV%; convert to the
    # internal log-normal variance via omega^2 = log(CV^2 + 1).
    etalcl ~ 0.03546   # CV 19% -> log(1 + 0.19^2)
    etalvc ~ 0.02527   # CV 16% -> log(1 + 0.16^2)

    # Residual error: NOT reported in the abstract. Fixed at zero so the
    # model can be loaded and simulated deterministically; see Errata in
    # the vignette. Users who re-fit on real data must replace this with
    # an estimated proportional (or combined) error term.
    propSd <- fixed(0); label("Proportional residual error (fraction; FIXED AT ZERO - not reported in source)")  # abstract: not reported
  })
  model({
    # 1. Individual PK parameters. CL and V scale LINEARLY with body weight
    #    (exponent = 1, per-kg parameterisation in the abstract); ASPHYXIA
    #    is a linear-deviation +13% increment on V.
    cl <- exp(lcl + etalcl) * WT
    vc <- exp(lvc + etalvc) * WT * (1 + e_asphyxia_vc * ASPHYXIA)

    # 2. Micro-constant for the 1-compartment ODE
    kel <- cl / vc

    # 3. ODE system: doses go directly into the central compartment
    #    (the abstract does not report an absorption rate constant).
    d/dt(central) <- -kel * central

    # 4. Observation (mg/L; dose in mg, volume in L) and proportional error.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
