Li_2023_RC88_human <- function() {
  description <- paste(
    "Human projection scaled from cynomolgus monkey (no human data were",
    "fitted). Full target-mediated drug disposition (TMDD) PK model of the",
    "mesothelin-targeting antibody-drug conjugate RC88 for patients with",
    "mesothelin-expressing solid tumors: the monkey two-compartment and",
    "binding parameters carried over per kg (allometric exponent 1, target",
    "rate constants 1:1) with baseline mesothelin tripled for the disease",
    "state (Li 2023)."
  )
  reference <- paste(
    "Li Q, Wang L, Zhang J, Zhao G, Liu Z, Ma X, Jiang J.",
    "Translation of the efficacy of antibody-drug conjugates from preclinical",
    "to clinical using a semimechanistic PK/PD model: A case study with RC88.",
    "Clin Transl Sci. 2023;16(7):1232-1242. doi:10.1111/cts.13526"
  )
  vignette <- "Li_2023_RC88"
  units <- list(time = "h", dosing = "mg/kg", concentration = "ug/mL")

  covariateData <- list()

  compartmentData <- list(
    central = list(analyte = "RC88 (antibody-drug conjugate)", units = "mg/kg", specimen = "serum", verified = TRUE),
    peripheral1 = list(
      analyte = "RC88 (antibody-drug conjugate)",
      units = "mg/kg",
      specimen = "tissue",
      verified = TRUE
    ),
    target = list(analyte = "mesothelin (free)", units = "ug/mL", specimen = "serum", verified = TRUE),
    complex = list(analyte = "RC88-mesothelin complex", units = "ug/mL", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 0L,
    disease_state = "patients with mesothelin-expressing malignant solid tumors (projected; the paper fits no human data)",
    dose_range = "0.82-1.96 mg/kg RC88 IV q1w (paper's predicted efficacious range)",
    notes = paste(
      "Projection only. Drug-specific parameters (CL, CLD, V1, V2) were scaled",
      "from monkey with an allometric exponent of 1, so the per-kg values are",
      "unchanged; target-specific rate constants (Kout, Kon, Koff, Ke) were",
      "scaled 1:1; baseline mesothelin was set to three times the monkey value",
      "because disease-state mesothelin is about three-fold higher than in",
      "healthy humans (Methods, 'Clinical PK predictions'; Table 2, 'Human'",
      "column). RC88 was in phase I at publication (NCT05508334)."
    )
  )

  ini({
    # Li 2023 Table 2, 'Human' column ('Scaled value'). These are projections,
    # not estimates, so every value is fixed.
    lcl <- fixed(log(0.556)); label("Linear clearance CL_ADC (mL/h/kg)") # Table 2 Human: CLADC = 0.556 mL/h/kg
    lq <- fixed(log(3.519)); label("Intercompartmental clearance CLD_ADC (mL/h/kg)") # Table 2 Human: CLDADC = 3.519 mL/h/kg
    lvc <- fixed(log(35.8)); label("Central volume of distribution V1 (mL/kg)") # Table 2 Human: V1 = 35.8 mL/kg
    lvp <- fixed(log(16.9)); label("Peripheral volume of distribution V2 (mL/kg)") # Table 2 Human: V2 = 16.9 mL/kg
    lkon <- fixed(log(0.00124)); label("Second-order association rate constant Kon (1/(ug/mL)/h)") # Table 2 Human: Kon = 0.00124
    lkoff <- fixed(log(0.00442)); label("First-order dissociation rate constant Koff (1/h)") # Table 2 Human: Koff = 0.00442
    lkdeg <- fixed(log(0.00496)); label("Fractional turnover rate of free mesothelin Kout (1/h)") # Table 2 Human: Kout = 0.00496 1/h
    lkint <- fixed(log(0.0154)); label("First-order elimination rate of the RC88-mesothelin complex Ke (1/h)") # Table 2 Human: Ke = 0.0154 1/h
    lrbase_target <- fixed(log(109.8)); label("Initial (baseline) mesothelin concentration R in patients (ug/mL)") # Table 2 Human: R = 109.8 ug/mL (3 x monkey 36.6)

    # No variability is defined for the projection.
    propSd <- fixed(0); label("Proportional residual error (fraction; not reported in the source)")
  })
  model({
    cl <- exp(lcl)
    q <- exp(lq)
    vc <- exp(lvc)
    vp <- exp(lvp)
    kon <- exp(lkon)
    koff <- exp(lkoff)
    kdeg <- exp(lkdeg)
    kint <- exp(lkint)
    rbase_target <- exp(lrbase_target)

    # Mesothelin synthesis balances its turnover at the baseline concentration.
    ksyn <- kdeg * rbase_target

    # Doses are mg/kg and volumes mL/kg, so central/vc is mg/mL; x1000 gives
    # ug/mL. Target and complex are carried as ug/mL concentrations in V1,
    # matching the units of R, Kon and Koff in Table 2.
    Cc <- 1000 * central / vc
    binding <- kon * Cc * target - koff * complex

    d/dt(central) <- -(cl / vc) * central - (q / vc) * central + (q / vp) * peripheral1 - binding * vc / 1000
    d/dt(peripheral1) <- (q / vc) * central - (q / vp) * peripheral1
    d/dt(target) <- ksyn - kdeg * target - binding
    d/dt(complex) <- binding - kint * complex

    target(0) <- rbase_target

    Cc ~ prop(propSd)
  })
}
