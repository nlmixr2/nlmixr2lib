Sharma_2018_SHetA2_mouse <- function() {
  description <- paste(
    "Preclinical (mouse, CD2F1 female).",
    "Two-compartment PK model with first-order absorption for SHetA2",
    "(a flexible heteroarotinoid anti-cancer / chemoprevention drug)",
    "in non-tumor-bearing CD2F1 female mice after intravenous (20 mg/kg)",
    "and oral (20, 60 mg/kg) administration. Disposition (CL, V1, V2,",
    "CLD) is reported as absolute total values for a typical 20-28 g",
    "mouse, fit by naive-pooled simultaneous IV+oral least-squares (1/y^2",
    "weighting) in Phoenix WinNonlin. No IIV was reported by the authors.",
    "Bioavailability F is dose-dependent (17.7% at 20 mg/kg, 19.5% at 60",
    "mg/kg) but kA is shared across doses. Parameter values from Sharma",
    "2018 Table 3."
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
    species        = "mouse (CD2F1 female)",
    n_subjects     = 9L,
    n_studies      = 1L,
    age_range      = "(not reported in the source publication)",
    weight_range   = "20.0-27.8 g",
    sex_female_pct = 100,
    disease_state  = "healthy, non-tumor-bearing; not fasted",
    dose_range     = paste(
      "IV: 20 mg/kg tail-vein injection in 100 uL of PEG400:Ethanol:Saline",
      "(57.1%:14.3%:28.6%), n=3. PO: 20 and 60 mg/kg by gavage in 100 uL",
      "sesame oil, n=3 per dose. Reference [13] = Liu 2007 (cited as source",
      "of the mouse PK data)."
    ),
    regions        = "USA (University of Oklahoma Health Sciences Center)",
    notes          = paste(
      "CD2F1 female mice (20.0-27.8 g) from a previously-published study",
      "(Sharma 2018 ref [13]). Blood collected by cardiac puncture at",
      "pre-dose and 0.08, 0.15, 0.25, 0.5, 1, 2, 3, 4, 6, 8, 12, 18, 24,",
      "36, 48, 60 h after IV dose; pre-dose and 0.25, 0.5, 1, 2, 3, 4, 6,",
      "8, 12, 18, 24, 36, 48 h after oral dose. Plasma SHetA2 quantified by",
      "SPE-HPLC-UV (341 nm); LLOQ 10 ng/mL. Concentrations were presented",
      "as mean +/- SD; IV and oral data were fit simultaneously by",
      "least-squares with 1/y^2 weighting. Body weight is not a model",
      "covariate; the absolute parameter values are reported for the",
      "typical 20-28 g animal. See Sharma 2018 Materials and methods,",
      "Data collection paragraph (mouse section), and Table 3."
    )
  )

  ini({
    # --------------------------------------------------------------
    # Structural PK -- Sharma 2018 Table 3, mouse column.
    # Volumes are absolute (L) and clearances absolute (L/h) at the
    # typical 20-28 g mouse body weight. CV% values in trailing
    # comments are precisions of the parameter estimates (RSE-like
    # standard errors from the Gauss-Newton fit), NOT inter-individual
    # variability -- the paper reports a naive-pooled deterministic
    # fit with no IIV.
    # --------------------------------------------------------------
    lka     <- log(0.539)   ; label("First-order absorption rate constant kA (1/h)")                 # Table 3: kA = 0.539 (CV 19%)
    lcl     <- log(0.0421)  ; label("Clearance CL (L/h) for typical 20-28 g mouse")                  # Table 3: CL = 0.0421 (CV 8%)
    lvc     <- log(0.121)   ; label("Central volume of distribution V1 (L) for typical mouse")       # Table 3: V1 = 0.121 (CV 16%)
    lq      <- log(0.0323)  ; label("Distributional clearance CLD (L/h) for typical mouse")          # Table 3: CLD = 0.0323 (CV 20%)
    lvp     <- log(0.282)   ; label("Peripheral volume V2 (L) for typical mouse")                    # Table 3: V2 = 0.282 (CV 14%)

    # --------------------------------------------------------------
    # Bioavailability. The paper estimated F separately at 20 and 60
    # mg/kg (17.7% and 19.5%); we encode 18.6% as a single typical
    # value -- this is the paper's own summary of "the maximum extent
    # of absorption of 18.6% at doses < 100 mg/kg" (Fig 4 caption).
    # The dose-dependence (17.7-19.5%) is preserved in the vignette
    # Assumptions and deviations section.
    # --------------------------------------------------------------
    lfdepot <- log(0.186)   ; label("Bioavailability F (fraction); typical value at doses <= 60 mg/kg") # Table 3 / Fig 4 caption: F = 17.7-19.5% (18.6% maximum extent)

    # --------------------------------------------------------------
    # Residual error. The paper fit the model using Phoenix WinNonlin
    # with a 1/y^2 weighting scheme (proportional-error equivalent) and
    # did not report the residual error magnitude. We encode a typical
    # proportional residual SD of 0.10 (10% CV) as a deterministic
    # placeholder for downstream simulations; this is NOT a fitted
    # value. Documented in the vignette Assumptions section.
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

    # 4. Observation: plasma concentration in ng/mL (paper Table 2).
    #    Cc [ng/mL] = central [mg] / vc [L] * 1000 [ng per ug] / 1 [ug per mg]
    #                = central / vc * 1000.
    Cc <- (central / vc) * 1000

    Cc ~ prop(propSd)
  })
}
