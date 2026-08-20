Wang_2024_saxagliptin_rat <- function() {
  description <- paste(
    "Preclinical (rat).",
    "Two-compartment PK model with first-order absorption for orally (intragastric) administered",
    "saxagliptin (SAX), directly linked to a sigmoidal Emax model WITHOUT a baseline term for the",
    "parent-drug contribution to plasma dipeptidyl peptidase-4 (DPP-4) inhibition, in",
    "streptozotocin plus high-fat-diet induced type 2 diabetic male Sprague-Dawley rats.",
    "Wang 2024 Eq 1-3 write the disposition as dA1/dt = F*Ka*Aa + (CL2/V2)*A2 - (CL2/V)*A1 - (CL/V)*A1 and",
    "dA2/dt = (CL2/V)*A1 - (CL2/V2)*A2, i.e. CL2 is a conventional inter-compartmental clearance",
    "(k12 = CL2/V, k21 = CL2/V2) and is encoded here as lq. The PD model is Wang 2024 Eq 8,",
    "E = Emax * C^Gam / (EC50^Gam + C^Gam), applied directly to the plasma SAX concentration with no",
    "effect compartment; the paper states explicitly that no hysteresis was observed.",
    "IMPORTANT: the fitted DPP-4 inhibition ratio for the parent drug is a DERIVED observation, not a",
    "directly measured one. Wang 2024 substituted the measured 5-OH SAX concentrations obtained after",
    "intragastric SAX into the separately fitted metabolite PK/PD model",
    "(see modellib('Wang_2024_5hydroxysaxagliptin_rat')) and subtracted the resulting metabolite",
    "inhibition from the measured total inhibition: In_SAX = In_(SAX+5OHSAX) - In_5OHSAX.",
    "The resulting EC50 ratio EC50_5OHSAX = 0.46 * EC50_SAX is the paper's central finding, confirming",
    "in vivo the roughly two-fold weaker potency of the metabolite reported in vitro.",
    "All volumes and clearances are body-weight normalised (mL/kg, mL/h/kg), so ODE state amounts are",
    "in ng per kg body weight and the 10 mg/kg dose is supplied as amt = 1e7.",
    "Parameters are the arithmetic means of three individually fitted rats (WinNonlin 8.1, naive",
    "individual fitting); no population model and therefore no OMEGA was estimated, so the model",
    "carries no eta terms and is intended for typical-value simulation. The PK residual error was",
    "described as additive but its magnitude was never reported, so addSd is FIXED at 0; the PD",
    "residual SD (stdev0) IS reported and is used.",
    sep = " "
  )
  reference <- paste(
    "Wang T, Tao T, Liu Y, Dong J, Ni S, Liu Y, Li Y, Xu N, Sun Z. (2024).",
    "Pharmacokinetic/Pharmacodynamic modelling of Saxagliptin and its active metabolite,",
    "5-hydroxy Saxagliptin in rats with Type 2 Diabetes Mellitus.",
    "BMC Pharmacology and Toxicology 25:35.",
    "doi:10.1186/s40360-024-00757-3.",
    sep = " "
  )
  vignette <- "Wang_2024_saxagliptin_dpp4_rat"
  units <- list(time = "h", dosing = "ng", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Amounts are body-weight normalised (ng per kg) because
  # Wang 2024 reports V in mL/kg and CL in mL/h/kg (Table 11).
  compartmentData <- list(
    depot       = list(analyte = "saxagliptin", units = "ng", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "saxagliptin", units = "ng", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "saxagliptin", units = "ng", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species       = "rat (Sprague-Dawley)",
    n_subjects    = 3L,
    n_studies     = 1L,
    sex_female_pct = 0,
    weight_range  = "120-140 g at purchase",
    disease_state = paste(
      "Type 2 diabetes mellitus induced by 4 weeks of high-fat diet followed by a single",
      "intraperitoneal injection of streptozocin 40 mg/kg in 0.1 mol/L citrate buffer (pH 4.2-4.5);",
      "rats qualified as T2DM with fasting blood glucose > 7.8 mmol/L and postprandial blood",
      "glucose > 16.7 mmol/L on days 3 and 7 after streptozocin"
    ),
    dose_range    = "10 mg/kg saxagliptin in 0.5% CMC-Na, single intragastric dose",
    regions       = "China",
    notes         = paste(
      "Materials and methods, 'Animals and study design' and 'Pharmacokinetic measurements'.",
      "Male SD rats from the Experimental Animal Center of Nantong University (SCXK (Su):2019-0001),",
      "housed at 22-25 degC on a 12 h light/dark cycle. A matched control group (n = 3) received",
      "0.5% CMC-Na vehicle only. Blood was sampled at 0, 0.08, 0.17, 0.25, 0.5, 0.75, 1, 1.5, 2, 3,",
      "4, 6, 8, 10, 12 and 24 h; only the first 4 h of sampling was used for model fitting",
      "('PK/PD link model' subsection of Results). Per-animal parameter estimates are in",
      "Supplementary Table S22-b; per-animal observations in Supplementary Table S22-a."
    )
  )

  ini({
    # Structural PK. Wang 2024 Table 11 / Supplementary Table S22-b report the arithmetic mean
    # of three individually fitted rats. Volumes are mL/kg and clearances mL/h/kg, so state
    # amounts are ng per kg body weight and Cc = central / vc is ng/mL.
    lka     <- log(9.9021)         ; label("First-order absorption rate constant Ka (1/h)")                 # Table 11: t_1/2,Ka = 0.07 h; Ka = log(2)/0.07 = 9.9021 (Ka itself is not tabulated)
    lcl     <- log(3245.80)        ; label("Apparent systemic clearance CL/F (mL/h/kg)")                    # Table 11: CL = 3245.80 mL/h/kg
    lvc     <- log(2307.24)        ; label("Apparent central volume of distribution V/F (mL/kg)")           # Table 11: V = 2307.24 mL/kg
    lq      <- log(17090.42)       ; label("Apparent inter-compartmental clearance CL2/F (mL/h/kg)")        # Table 11: CL2 = 17090.42 mL/h/kg
    lvp     <- log(18393.30)       ; label("Apparent peripheral volume of distribution V2/F (mL/kg)")       # Table 11: V2 = 18393.30 mL/kg
    lfdepot <- fixed(log(1))       ; label("Oral bioavailability F (fraction); assumed 1, not published")   # F appears in Eq 1 but is never estimated or reported, so CL, V, CL2 and V2 are apparent (CL/F etc.)

    # Direct-effect sigmoidal Emax on the parent-drug DPP-4 inhibition ratio (Wang 2024 Eq 8).
    # There is no E0 term for the parent: the inhibition attributed to SAX is defined as a
    # difference of two inhibition ratios and is therefore 0 at time 0 by construction.
    lemax   <- log(71.47)          ; label("Maximum parent-SAX DPP-4 inhibition ratio Emax (%)")            # Table 12: Emax = 71.47%
    lec50   <- log(544.74)         ; label("SAX concentration producing half of Emax, EC50 (ng/mL)")        # Table 12: EC50 = 544.74 ng/mL
    lhill   <- log(1.38)           ; label("Sigmoidicity (Hill) exponent Gam (unitless)")                   # Table 12: Gam = 1.38

    # Residual error, both additive (Results, 'PK model of parent SAX in T2DM rats').
    addSd          <- fixed(0)     ; label("Additive residual SD for plasma SAX concentration (ng/mL); not published")  # PK error model stated to be additive but no magnitude is reported in Table 11 or Supplementary Table S22-b
    addSd_DPP4inh  <- 2.01         ; label("Additive residual SD for the parent-SAX DPP-4 inhibition ratio (%)")        # Table 12: stdev0 = 2.01
  })

  model({
    # Individual parameters (no IIV: Wang 2024 fitted each rat separately in WinNonlin 8.1
    # and reported means and SDs of the individual estimates rather than a population OMEGA).
    ka   <- exp(lka)
    cl   <- exp(lcl)
    vc   <- exp(lvc)
    q    <- exp(lq)
    vp   <- exp(lvp)
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # Micro-constants. Wang 2024 Eq 1-2 use CL/V, CL2/V and CL2/V2 directly; Eq 3 defines K = CL/V.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot + k21 * peripheral1 - k12 * central - kel * central
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot) <- exp(lfdepot)

    # Observations. Cc is ng/mL (ng per kg divided by mL per kg).
    Cc <- central / vc

    # Wang 2024 Eq 8: sigmoidal Emax with no baseline, driven directly by plasma SAX.
    DPP4inh <- emax * Cc^hill / (ec50^hill + Cc^hill)

    Cc      ~ add(addSd)
    DPP4inh ~ add(addSd_DPP4inh)
  })
}
