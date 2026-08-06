Wang_2024_5hydroxysaxagliptin_rat <- function() {
  description <- paste(
    "Preclinical (rat).",
    "Two-compartment intravenous PK model for the active saxagliptin metabolite 5-hydroxy",
    "saxagliptin (5-OH SAX), directly linked to a sigmoidal Emax model WITH a baseline term E0 for",
    "plasma dipeptidyl peptidase-4 (DPP-4) inhibition, in streptozotocin plus high-fat-diet induced",
    "type 2 diabetic male Sprague-Dawley rats.",
    "Wang 2024 Eq 4-6 write the disposition as dA1/dt = (CL2/V2)*A2 - (CL2/V)*A1 - (CL/V)*A1 and",
    "dA2/dt = (CL2/V)*A1 - (CL2/V2)*A2, i.e. CL2 is a conventional inter-compartmental clearance",
    "(k12 = CL2/V, k21 = CL2/V2) and is encoded here as lq. The PD model is Wang 2024 Eq 7,",
    "E = E0 + Emax * C^Gam / (EC50^Gam + C^Gam), applied directly to the plasma 5-OH SAX",
    "concentration with no effect compartment; the paper states explicitly that 5-OH SAX binds DPP-4",
    "immediately and that no hysteresis was observed.",
    "This model is the metabolite half of Wang 2024 and is the model the authors used to attribute",
    "part of the total DPP-4 inhibition seen after oral saxagliptin to the metabolite; the parent-drug",
    "half is modellib('Wang_2024_saxagliptin_rat'). The fitted potencies satisfy",
    "EC50_5OHSAX = 0.46 * EC50_SAX, the paper's central finding.",
    "All volumes and clearances are body-weight normalised (mL/kg, mL/h/kg), so ODE state amounts are",
    "in ng per kg body weight and the 0.5 mg/kg intravenous dose is supplied as amt = 5e5 into central.",
    "Parameters are the arithmetic means of three individually fitted rats (WinNonlin 8.1, naive",
    "individual fitting); no population model and therefore no OMEGA was estimated, so the model",
    "carries no eta terms and is intended for typical-value simulation. The PK residual error was",
    "described as additive but its magnitude was never reported, so addSd is FIXED at 0; the PD",
    "residual SD (stdev0) IS reported and is used. E0 is kept on the linear scale because the",
    "individual estimates in Supplementary Table S19-b include a negative value (-3.85%).",
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
  # Wang 2024 reports V in mL/kg and CL in mL/h/kg (Table 9).
  compartmentData <- list(
    central     = list(analyte = "5-hydroxy saxagliptin", units = "ng", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "5-hydroxy saxagliptin", units = "ng", specimen = "plasma", verified = TRUE)
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
    dose_range    = "0.5 mg/kg 5-hydroxy saxagliptin in saline, single intravenous tail-vein dose",
    regions       = "China",
    notes         = paste(
      "Materials and methods, 'Animals and study design' and 'Pharmacokinetic measurements'.",
      "Male SD rats from the Experimental Animal Center of Nantong University (SCXK (Su):2019-0001),",
      "housed at 22-25 degC on a 12 h light/dark cycle. A matched control group (n = 3) received",
      "saline only. Blood was sampled at 0, 0.03, 0.08, 0.17, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8,",
      "10, 12 and 24 h; only the first 2 h of sampling was used for model fitting ('PK/PD link model'",
      "subsection of Results). Per-animal parameter estimates are in Supplementary Table S19-b;",
      "per-animal observations in Supplementary Table S19-a."
    )
  )

  ini({
    # Structural PK. Wang 2024 Table 9 / Supplementary Table S19-b report the arithmetic mean
    # of three individually fitted rats. Volumes are mL/kg and clearances mL/h/kg, so state
    # amounts are ng per kg body weight and Cc = central / vc is ng/mL.
    lcl   <- log(882.19)  ; label("Systemic clearance CL (mL/h/kg)")                                # Table 9: CL = 882.19 mL/h/kg
    lvc   <- log(141.49)  ; label("Central volume of distribution V (mL/kg)")                       # Table 9: V = 141.49 mL/kg
    lq    <- log(471.60)  ; label("Inter-compartmental clearance CL2 (mL/h/kg)")                    # Table 9: CL2 = 471.60 mL/h/kg
    lvp   <- log(137.81)  ; label("Peripheral volume of distribution V2 (mL/kg)")                   # Table 9: V2 = 137.81 mL/kg

    # Direct-effect sigmoidal Emax with baseline on the DPP-4 inhibition ratio (Wang 2024 Eq 7).
    lemax <- log(60.88)   ; label("Maximum 5-OH SAX DPP-4 inhibition ratio Emax (%)")               # Table 10: Emax = 60.88%
    lec50 <- log(251.74)  ; label("5-OH SAX concentration producing half of Emax, EC50 (ng/mL)")    # Table 10: EC50 = 251.74 ng/mL
    lhill <- log(1.31)    ; label("Sigmoidicity (Hill) exponent Gam (unitless)")                    # Table 10: Gam = 1.31
    e0    <- 7.03         ; label("Baseline DPP-4 inhibition ratio without drug E0 (%)")            # Table 10: E0 = 7.03; kept linear because individual estimates span -3.85 to 21.15 (Supplementary Table S19-b)

    # Residual error, both additive (Results, 'PK model of 5-OH SAX in T2DM rats').
    addSd         <- fixed(0) ; label("Additive residual SD for plasma 5-OH SAX concentration (ng/mL); not published")  # PK error model stated to be additive but no magnitude is reported in Table 9 or Supplementary Table S19-b
    addSd_DPP4inh <- 1.97     ; label("Additive residual SD for the DPP-4 inhibition ratio (%)")                        # Table 10: stdev0 = 1.97
  })

  model({
    # Individual parameters (no IIV: Wang 2024 fitted each rat separately in WinNonlin 8.1
    # and reported means and SDs of the individual estimates rather than a population OMEGA).
    cl   <- exp(lcl)
    vc   <- exp(lvc)
    q    <- exp(lq)
    vp   <- exp(lvp)
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # Micro-constants. Wang 2024 Eq 4-5 use CL/V, CL2/V and CL2/V2 directly; Eq 6 defines K = CL/V.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central)     <-  k21 * peripheral1 - k12 * central - kel * central
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Observations. Cc is ng/mL (ng per kg divided by mL per kg).
    Cc <- central / vc

    # Wang 2024 Eq 7: sigmoidal Emax with baseline E0, driven directly by plasma 5-OH SAX.
    DPP4inh <- e0 + emax * Cc^hill / (ec50^hill + Cc^hill)

    Cc      ~ add(addSd)
    DPP4inh ~ add(addSd_DPP4inh)
  })
}
