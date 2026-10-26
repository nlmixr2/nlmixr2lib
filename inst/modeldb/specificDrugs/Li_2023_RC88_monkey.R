Li_2023_RC88_monkey <- function() {
  description <- paste(
    "Preclinical (cynomolgus monkey). Full target-mediated drug disposition",
    "(TMDD) PK model of the mesothelin-targeting antibody-drug conjugate RC88:",
    "two-compartment linear disposition plus reversible binding to mesothelin",
    "in the central compartment, mesothelin turnover, and internalisation of",
    "the drug-target complex, after single IV doses of 0.5, 1.5 and 5 mg/kg",
    "(Li 2023)."
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
    species = "cynomolgus monkey",
    n_subjects = 18L,
    n_studies = 1L,
    sex_female_pct = 50,
    disease_state = "healthy (non-tumor-bearing); RC88 cross-reacts with monkey but not rodent mesothelin",
    dose_range = "0.5, 1.5 and 5 mg/kg RC88, single IV dose (n = 6 per group)",
    regions = "China (preclinical study, RemeGen)",
    notes = paste(
      "Serum sampled pre-dose and at 0.083, 0.25, 2, 6, 24, 48, 72, 96, 168,",
      "336, 504 and 672 h post-dose (Methods, 'Cynomolgus monkey PK studies').",
      "The Results text gives the top dose as 3 mg/kg; Methods, Figure 3b and",
      "Table S1 all give 5 mg/kg."
    )
  )

  ini({
    # Li 2023 Table 2, 'Monkey' column. Values in parentheses are bootstrap %RSE.
    lcl <- log(0.556); label("Linear clearance CL_ADC (mL/h/kg)") # Table 2: CLADC = 0.556 mL/h/kg (7.9%)
    lq <- log(3.519); label("Intercompartmental clearance CLD_ADC (mL/h/kg)") # Table 2: CLDADC = 3.519 mL/h/kg (16.3)
    lvc <- fixed(log(35.8)); label("Central volume of distribution V1 (mL/kg)") # Table 2: V1 = 35.8 mL/kg (fixed)
    lvp <- log(16.9); label("Peripheral volume of distribution V2 (mL/kg)") # Table 2: V2 = 16.9 mL/kg (9.7%)
    lkon <- log(0.00124); label("Second-order association rate constant Kon (1/(ug/mL)/h)") # Table 2: Kon = 0.00124 (10.7%)
    lkoff <- log(0.00442); label("First-order dissociation rate constant Koff (1/h)") # Table 2: Koff = 0.00442 (19.3%); unit printed '1/ug/mL*h'
    lkdeg <- log(0.00496); label("Fractional turnover rate of free mesothelin Kout (1/h)") # Table 2: Kout = 0.00496 1/h (14.9%)
    lkint <- log(0.0154); label("First-order elimination rate of the RC88-mesothelin complex Ke (1/h)") # Table 2: Ke = 0.0154 1/h (16.6%)
    lrbase_target <- log(36.6); label("Initial (baseline) mesothelin concentration R (ug/mL)") # Table 2: R = 36.6 ug/mL (11.8%)

    # Li 2023 reports neither interindividual nor residual variability for
    # this fit. Carried as fixed(0) so the model simulates the published
    # typical-value profile.
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
