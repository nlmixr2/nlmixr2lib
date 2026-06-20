Sharma_2018_SHetA2_dog <- function() {
  description <- paste(
    "Preclinical (beagle dog).",
    "Two-compartment PK model with a 7-compartment gastrointestinal",
    "(GI) transit absorption process for SHetA2 (a flexible",
    "heteroarotinoid anti-cancer / chemoprevention drug) in beagle dogs",
    "after intravenous (5 mg/kg) and oral (100, 400, 1500 mg/kg)",
    "administration. Drug transits through 7 serial GI segments (stomach",
    "G1 = depot, then transit1..transit6 = G2..G7) at a common transit",
    "rate kAT; absorption occurs only from G2 (transit1, rate kA) and",
    "G7 (transit6, rate kA2). Disposition (CL, V1, V2, CLD) is reported",
    "as absolute total values for a typical 6.4-11.2 kg dog, fit by",
    "naive-pooled simultaneous IV+oral least-squares (1/y^2 weighting)",
    "in Phoenix WinNonlin. No IIV was reported. F, kA, kA2, and kAT are",
    "all dose-dependent; this file encodes the 100 mg/kg parameter set",
    "(F=11.2%, kA=1.12/h, kA2=0.929/h, kAT=0.532/h) as the typical",
    "value because 100 mg/kg is the new NOAEL used to derive the",
    "first-in-human dose. Higher-dose parameter sets (400 and 1500",
    "mg/kg) are documented in the vignette. Parameter values from",
    "Sharma 2018 Table 3."
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
    species        = "beagle dog",
    n_subjects     = 14L,
    n_studies      = 1L,
    age_range      = "(not reported in the source publication)",
    weight_range   = "6.4-11.2 kg",
    sex_female_pct = NA,
    disease_state  = "healthy; not fasted",
    dose_range     = paste(
      "IV: 5 mg/kg (n=2; formulation not reported). PO: 100, 400, or",
      "1500 mg/kg in 30% aqueous Solutol HS 15 by gavage at 5 mL/kg/day",
      "(n=4 per oral dose). IV data from NCI PREVENT report",
      "N01-CN-43306 work assignment 21 (pilot PK / tolerability)."
    ),
    regions        = "USA (NCI / University of Oklahoma Health Sciences Center)",
    notes          = paste(
      "Beagle dogs (males 9.0-11.2 kg; females 6.9-8.4 kg for the oral",
      "arms; IV arm 6.4-7.0 kg). Blood from jugular vein at pre-dose and",
      "1, 2, 3, 4, 6, 9, 24 h after oral administration; pre-dose and",
      "0.08, 0.25, 0.5, 1, 2, 4, 8, 24 h after IV. PO data are mean +/-",
      "SD; IV data are mean (n=2 only). Plasma SHetA2 quantified by",
      "LC/MS/MS (gradient method, 0.01% TFA in water / acetonitrile);",
      "LLOQ 10 ng/mL. Naive-pooled simultaneous IV+oral fit by",
      "least-squares with 1/y^2 weighting. The dog data showed",
      "double-peak oral concentration-time profiles at the 400 and 1500",
      "mg/kg doses that required two absorption sites (G2 and G7) to",
      "capture; at 100 mg/kg only the early peak is dominant. Body",
      "weight is not a model covariate; the absolute parameter values",
      "are reported for the typical 6.4-11.2 kg animal. See Sharma 2018",
      "Materials and methods, Data collection paragraph (dog section),",
      "and Table 3."
    )
  )

  ini({
    # --------------------------------------------------------------
    # Structural PK -- Sharma 2018 Table 3, dog column.
    # Volumes absolute (L), clearances absolute (L/h) at the typical
    # 6.4-11.2 kg dog. CV% in trailing comments are precisions of the
    # parameter estimates (RSE-like SEs from Gauss-Newton), NOT IIV.
    # --------------------------------------------------------------
    lcl     <- log(6.39)    ; label("Clearance CL (L/h) for typical 6.4-11.2 kg dog")                 # Table 3: CL = 6.39 (CV 10%)
    lvc     <- log(8.53)    ; label("Central volume of distribution V1 (L) for typical dog")          # Table 3: V1 = 8.53 (CV 18%)
    lq      <- log(3.24)    ; label("Distributional clearance CLD (L/h) for typical dog")             # Table 3: CLD = 3.24 (CV 29%)
    lvp     <- log(22.5)    ; label("Peripheral volume V2 (L) for typical dog")                       # Table 3: V2 = 22.5 (CV 26%)

    # --------------------------------------------------------------
    # GI transit absorption parameters at the 100 mg/kg dose level.
    # The paper estimated dose-dependent F / kA / kA2 / kAT separately
    # for each oral dose (Table 3). We encode the 100 mg/kg set as the
    # typical value because 100 mg/kg is the new NOAEL used to derive
    # the first-in-human dose (Sharma 2018 Discussion). The 400 and
    # 1500 mg/kg parameter sets are recorded in the vignette and may
    # be loaded by zeroing out lka2 / re-fitting locally as needed.
    # --------------------------------------------------------------
    lka     <- log(1.12)    ; label("Absorption rate from G2 (kA, 1/h) at 100 mg/kg PO")              # Table 3: kA = 1.12 at 100 mg/kg (CV 43%)
    lka2    <- log(0.929)   ; label("Absorption rate from G7 (kA2, 1/h) at 100 mg/kg PO")             # Table 3: kA2 = 0.929 at 100 mg/kg (CV 281%)
    lktr    <- log(0.532)   ; label("GI transit rate constant kAT (1/h), common across segments")     # Table 3: kAT = 0.532 at 100 mg/kg (CV 35%)
    lfdepot <- log(0.112)   ; label("Bioavailability F (fraction) at 100 mg/kg PO; dose-dependent")   # Table 3: F = 11.2% at 100 mg/kg (CV 16%)

    # --------------------------------------------------------------
    # Residual error: 1/y^2 weighting was used in Phoenix WinNonlin
    # but the residual error magnitude was not reported. Encoding a
    # typical 10% proportional SD as a deterministic placeholder;
    # NOT a fitted value. Documented in the vignette Assumptions.
    # --------------------------------------------------------------
    propSd  <- fixed(0.10)  ; label("Proportional residual error (fraction); placeholder, not from paper") # Materials and methods: 1/y^2 weighting (magnitude not reported)
  })

  model({
    # 1. Individual parameters (typical-value model; no eta).
    ka     <- exp(lka)
    ka2    <- exp(lka2)
    ktr    <- exp(lktr)
    cl     <- exp(lcl)
    vc     <- exp(lvc)
    q      <- exp(lq)
    vp     <- exp(lvp)
    fdepot <- exp(lfdepot)

    # 2. GI transit absorption + two-compartment disposition
    #    (Sharma 2018 Eqs 4-9, Fig 2B). Seven serial gut compartments:
    #      depot   = G1 stomach (transit out only; no absorption)
    #      transit1 = G2 (transit out at ktr AND absorbs to central at ka)
    #      transit2 = G3, transit3 = G4, transit4 = G5, transit5 = G6
    #                 (transit-only -- no absorption from segments 3-6)
    #      transit6 = G7 (transit out at ktr -- lost from system --
    #                     AND absorbs to central at ka2)
    #    State amounts in mg.
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot     - (ktr + ka)  * transit1
    d/dt(transit2)    <-  ktr * transit1  - ktr * transit2
    d/dt(transit3)    <-  ktr * transit2  - ktr * transit3
    d/dt(transit4)    <-  ktr * transit3  - ktr * transit4
    d/dt(transit5)    <-  ktr * transit4  - ktr * transit5
    d/dt(transit6)    <-  ktr * transit5  - (ktr + ka2) * transit6
    d/dt(central)     <-  ka  * transit1  + ka2 * transit6 -
                          (cl + q) / vc * central + q / vp * peripheral1
    d/dt(peripheral1) <-  q / vc * central - q / vp * peripheral1

    # 3. Bioavailability applied to depot (G1 stomach).
    f(depot) <- fdepot

    # 4. Observation: plasma concentration in ng/mL.
    Cc <- (central / vc) * 1000

    Cc ~ prop(propSd)
  })
}
