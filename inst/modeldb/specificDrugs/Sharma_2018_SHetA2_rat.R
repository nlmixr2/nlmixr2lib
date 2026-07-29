Sharma_2018_SHetA2_rat <- function() {
  description <- paste(
    "Preclinical (rat, Crl:CD Sprague-Dawley).",
    "Two-compartment PK model with first-order absorption for SHetA2",
    "(a flexible heteroarotinoid anti-cancer / chemoprevention drug)",
    "in Crl:CD (SD) rats after intravenous (5 mg/kg single dose) and",
    "oral (100, 500, 2000 mg/kg/day for 28 days) administration.",
    "Disposition (CL, V1, V2, CLD) is reported as absolute total",
    "values for a typical 260-347 g rat, fit by naive-pooled",
    "simultaneous IV+oral least-squares (1/y^2 weighting) in Phoenix",
    "WinNonlin. No IIV was reported. The absorption is slow with",
    "flip-flop kinetics at the higher oral doses; kA was estimated as",
    "a single value across doses (0.0755 1/h) while F is dose-dependent",
    "(1.03% at 100 mg/kg, 1.57% at 500 mg/kg, 0.560% at 2000 mg/kg).",
    "Parameter values from Sharma 2018 Table 3."
  )
  reference <- paste(
    "Sharma A, Benbrook DM, Woo S. (2018).",
    "Pharmacokinetics and interspecies scaling of a novel,",
    "orally-bioavailable anti-cancer drug, SHetA2.",
    "PLoS ONE 13(4):e0194046.",
    "doi:10.1371/journal.pone.0194046.",
    sep = " "
  )
  vignette <- "Sharma_2018_SHetA2"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  covariateData <- list()

  population <- list(
    species        = "rat (Crl:CD Sprague-Dawley)",
    n_subjects     = 12L,
    n_studies      = 1L,
    age_range      = "(not reported in the source publication)",
    weight_range   = "260-347 g",
    sex_female_pct = NA,
    disease_state  = "healthy; not fasted",
    dose_range     = paste(
      "IV: single 5 mg/kg dose in PEG400:Ethanol:Saline",
      "(57.1%:14.3%:28.6%). PO: 100, 500, or 2000 mg/kg/day in 1%",
      "methylcellulose / 0.2% Tween 80 by gavage at 10 mL/kg/day for 28",
      "days. n=3 per dose level. PK profiles were sampled on the IV day",
      "and on week 4/5 for the oral arm. Reference [11] = the underlying",
      "preclinical 28-day toxicology study."
    ),
    regions        = "USA",
    notes          = paste(
      "Crl:CD (SD) rats (260-347 g). Serial blood by retro-orbital",
      "puncture at 0, 0.08, 0.25, 0.5, 1, 2, 4, 6 h after IV dosing and",
      "at 0, 0.5, 1, 2, 4, 6, 8, 24 h after oral dosing on week 4 or 5.",
      "Rat PK data were digitized using GetData Graph Digitizer v2.26.",
      "Plasma SHetA2 quantified by LC/MS/MS (gradient method, 0.01% TFA",
      "in water / acetonitrile); LLOQ 10 ng/mL. Naive-pooled simultaneous",
      "IV+oral fit by least-squares with 1/y^2 weighting. Body weight is",
      "not a model covariate; the absolute parameter values are reported",
      "for the typical 260-347 g animal. The model assumes a single",
      "absorption rate constant across doses; encoding the lowest oral",
      "dose (100 mg/kg) bioavailability value as the typical F. See",
      "Sharma 2018 Materials and methods, Data collection paragraph (rat",
      "section), and Table 3."
    )
  )

  ini({
    # --------------------------------------------------------------
    # Structural PK -- Sharma 2018 Table 3, rat column.
    # Volumes absolute (L), clearances absolute (L/h) at the typical
    # 260-347 g rat. CV% in trailing comments are precisions of the
    # parameter estimates (RSE-like SEs from Gauss-Newton), NOT IIV
    # -- the paper reports a naive-pooled deterministic fit.
    # --------------------------------------------------------------
    lka     <- log(0.0755)  ; label("First-order absorption rate constant kA (1/h)")                  # Table 3: kA = 0.0755 (CV 17%); single value across doses
    lcl     <- log(0.916)   ; label("Clearance CL (L/h) for typical 260-347 g rat")                   # Table 3: CL = 0.916 (CV 10%)
    lvc     <- log(0.969)   ; label("Central volume of distribution V1 (L) for typical rat")          # Table 3: V1 = 0.969 (CV 18%)
    lq      <- log(0.286)   ; label("Distributional clearance CLD (L/h) for typical rat")             # Table 3: CLD = 0.286 (CV 70%)
    lvp     <- log(0.531)   ; label("Peripheral volume V2 (L) for typical rat")                       # Table 3: V2 = 0.531 (CV 48%)

    # --------------------------------------------------------------
    # Bioavailability. The paper estimated F separately at 100, 500,
    # and 2000 mg/kg (1.03%, 1.57%, 0.560%). We encode the 100 mg/kg
    # value (the lowest dose, typically closest to true bioavailability
    # before absorption saturation) as the typical F. The dose-
    # dependence is preserved in the vignette Assumptions section.
    # --------------------------------------------------------------
    lfdepot <- log(0.0103)  ; label("Bioavailability F (fraction) at 100 mg/kg PO; dose-dependent")    # Table 3: F = 1.03% at 100 mg/kg (CV 14%)

    # --------------------------------------------------------------
    # Residual error: 1/y^2 weighting was used in Phoenix WinNonlin
    # but the residual error magnitude was not reported. Encoding
    # a typical 10% proportional SD as a deterministic placeholder;
    # NOT a fitted value. Documented in the vignette Assumptions.
    # --------------------------------------------------------------
    propSd  <- fixed(0.10)  ; label("Proportional residual error (fraction); placeholder, not from paper") # Materials and methods: 1/y^2 weighting (magnitude not reported)
  })

  model({
    # 1. Individual parameters (typical-value model; no eta).
    ka     <- exp(lka)
    cl     <- exp(lcl)
    vc     <- exp(lvc)
    q      <- exp(lq)
    vp     <- exp(lvp)
    fdepot <- exp(lfdepot)

    # 2. Two-compartment disposition with first-order absorption
    #    (Sharma 2018 Eqs 1-3, Fig 2A). State amounts in mg; volumes
    #    in L; concentrations in mg/L * 1000 -> ng/mL.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - (cl + q) / vc * central + q / vp * peripheral1
    d/dt(peripheral1) <-  q / vc * central - q / vp * peripheral1

    # 3. Bioavailability applied to depot.
    f(depot) <- fdepot

    # 4. Observation: plasma concentration in ng/mL.
    Cc <- (central / vc) * 1000

    Cc ~ prop(propSd)
  })
}
