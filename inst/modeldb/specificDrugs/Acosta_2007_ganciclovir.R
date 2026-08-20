Acosta_2007_ganciclovir <- function() {
  description <- paste(
    "One-compartment population PK model for ganciclovir after IV ganciclovir and",
    "an oral liquid valganciclovir formulation in neonates with symptomatic",
    "congenital CMV disease (Acosta 2007), with first-order absorption, estimated",
    "bioavailability, and a steep un-normalized body-weight power effect on",
    "clearance.",
    "Parameters transcribed from the Yang 2023 ganciclovir / valganciclovir",
    "population-PK model repository (Table 3), not from the primary publication;",
    "re-verify against Acosta 2007 when the primary is obtained.",
    sep = " "
  )
  reference <- paste(
    "Acosta EP, Brundage RC, King JR, Sanchez PJ, Sood S, Agrawal V, Homans J,",
    "Jacobs RF, Lang D, Romero JR, et al. Ganciclovir population pharmacokinetics",
    "in neonates following intravenous administration of ganciclovir and oral",
    "administration of a liquid valganciclovir formulation.",
    "Clin Pharmacol Ther. 2007;81(6):867-872. doi:10.1038/sj.clpt.6100150.",
    "Parameters transcribed from Yang W, Mak W, Gwee A, Gu M, Wu Y, Shi Y, He Q,",
    "Xiang X, Han B, Zhu X. Establishment and Evaluation of a Parametric Population",
    "Pharmacokinetic Model Repository for Ganciclovir and Valganciclovir.",
    "Pharmaceutics. 2023;15(7):1801. doi:10.3390/pharmaceutics15071801 (Table 3).",
    sep = " "
  )
  vignette <- "Yang_2023_ganciclovir_model_repository"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    depot   = list(analyte = "ganciclovir", units = "mg", specimen = "administration site", verified = FALSE),
    central = list(analyte = "ganciclovir", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters as RAW (un-normalized) power terms, WT^1.68 on CL and WT^1 on V,",
        "so the structural reference subject is WT = 1 kg rather than a cohort",
        "median. The CL exponent of 1.68 is far above the canonical allometric 0.75",
        "and is the steepest body-size effect in the Yang 2023 repository;",
        "Yang 2023 Section 4.3 notes that in neonates the weight effect on CL is",
        "marked, with clearance at 5 kg being at least 3.34-fold the clearance at",
        "1 kg. Weight was the only covariate retained (Yang 2023 Table 4).",
        "Cohort weights: 2.7 kg (range 2.1-3.4) in study 1.0 and 2.9 kg",
        "(range 1.9-4.4) in study 2.0 (Yang 2023 Table 2).",
        sep = " "
      ),
      source_name        = "BW"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 24L,
    n_studies      = 2L,
    n_observations = 484L,
    age_median     = paste(
      "Postnatal age: study 1.0 30 days (range 11-34); study 2.0 20 days",
      "(range 8-33).",
      sep = " "
    ),
    weight_median  = "study 1.0 2.7 kg (range 2.1-3.4); study 2.0 2.9 kg (range 1.9-4.4)",
    sex_female_pct = 45.8,
    race_ethnicity = "Not reported.",
    disease_state  = "Neonates with symptomatic congenital CMV disease.",
    dose_range     = paste(
      "IV ganciclovir 6 mg/kg/12 h; oral liquid valganciclovir formulation",
      "14 mg/kg/12 h.",
      sep = " "
    ),
    regions        = "United States (prospective; two sequential studies pooled).",
    bioassay       = "LC-MS, LLOQ 0.4 ug/mL.",
    notes          = paste(
      "Demographics and dosing from Yang 2023 Table 2. Sampling was both intensive",
      "(0 pre-dose, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12 and/or 24 h post dose) and",
      "sparse (0 pre-dose, 1, 2, 4, 8, 12 and/or 16 and 24 h post dose).",
      "Covariates tested: weight, sex, postnatal age and BSA; only weight was",
      "retained (Yang 2023 Table 4). This study did NOT convert valganciclovir",
      "doses to ganciclovir equivalents; oral doses are supplied as valganciclovir",
      "milligrams and the estimated bioavailability F = 0.536 absorbs the",
      "molecular-weight ratio. The simulation target in the primary study was an",
      "AUC0-12h of 27 mg*h/L.",
      sep = " "
    )
  )

  ini({
    # Structural PK -- Yang 2023 Table 3, Acosta et al. (2007) row. Body weight
    # enters un-normalized, so the reference subject is WT = 1 kg. Clearance in
    # L/h, volume in L, ka in 1/h. One-compartment model with a depot.
    lcl <- log(0.146); label("Clearance coefficient at WT = 1 kg (CL, L/h)")            # Yang 2023 Table 3 (Acosta 2007): CL = 0.146 * BW^1.68
    lvc <- log(1.15) ; label("Central volume coefficient at WT = 1 kg (V, L)")          # Yang 2023 Table 3 (Acosta 2007): V = 1.15 * BW
    lka <- log(0.591); label("First-order oral absorption rate constant (ka, 1/h)")     # Yang 2023 Table 3 (Acosta 2007): Ka = 0.591

    # Oral bioavailability of ganciclovir from a valganciclovir milligram dose of
    # the liquid formulation (no molecular-weight conversion was applied).
    lfdepot <- log(0.536); label("Oral bioavailability of ganciclovir from liquid valganciclovir (F, fraction)")  # Yang 2023 Table 3 (Acosta 2007): F = 0.536

    # Covariate effects. The CL exponent 1.68 is a non-canonical estimated value;
    # the volume exponent is 1 (linear in weight) as printed and is fixed.
    e_wt_cl <- 1.68    ; label("Power exponent of WT on CL (unitless; raw WT, reference 1 kg)")           # Yang 2023 Table 3 (Acosta 2007): BW^1.68
    e_wt_vc <- fixed(1); label("Power exponent of WT on V (unitless; raw WT, linear, reference 1 kg)")    # Yang 2023 Table 3 (Acosta 2007): 1.15 * BW

    # Between-subject variability. Yang 2023 Methods: %CV = sqrt(omega^2) * 100%,
    # so variance = (BSV% / 100)^2. V and ka carry no BSV in the source table.
    etalcl     ~ 0.080656  # Yang 2023 Table 3 (Acosta 2007): BSV CL = 28.4% -> 0.284^2
    etalfdepot ~ 0.015376  # Yang 2023 Table 3 (Acosta 2007): BSV F  = 12.4% -> 0.124^2

    # Residual unexplained variability. Yang 2023 Table 3 reports 45.4%
    # "exponential error"; an exponential residual on the natural scale is an
    # additive residual on the log scale, which nlmixr2 expresses as a proportional
    # error in linear space.
    propSd <- 0.454; label("Proportional residual error (fraction; exponential residual in the source)")  # Yang 2023 Table 3 (Acosta 2007): 45.4% exponential error
  })

  model({
    cl <- exp(lcl + etalcl) * WT^e_wt_cl
    vc <- exp(lvc)          * WT^e_wt_vc
    ka <- exp(lka)
    fdepot <- exp(lfdepot + etalfdepot)

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Oral valganciclovir doses enter `depot` as valganciclovir milligrams;
    # IV ganciclovir doses go directly to `central`.
    f(depot) <- fdepot

    # Dose in mg, volume in L -> concentration in mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
