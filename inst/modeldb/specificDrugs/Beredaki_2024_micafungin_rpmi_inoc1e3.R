Beredaki_2024_micafungin_rpmi_inoc1e3 <- function() {
  description <- "In vitro (Candida parapsilosis sensu stricto clinical isolates in serum-free RPMI-1640, starting inoculum 10^3 CFU/mL; two-compartment closed dialysis/diffusion PK/PD model). Micafungin in vitro PK/PD model without serum. Micafungin was added once daily to the 10 mL internal compartment at a target peak (fCmax 0.5-4 mg/L; the protein-free medium leaves all drug unbound) and decayed mono-exponentially with the achieved in vitro half-life of 8 h. The 72 h change in log10 CFU/mL from the starting inoculum was related to the PK/PD index fAUC0-24/MIC by the sigmoidal variable-slope Emax model E = (Emax - Emin) * EI^n / (EI^n + EI50^n) + Emin, with Emax constrained to -2 log10 CFU/mL for this inoculum (Methods, 'PK/PD analysis'). Beredaki 2024 prints no Emin or n for the serum-free curves, so Emin = 4.40, EI50 = 41.3 fAUC0-24/MIC and n = 4.82 were recovered by digitising the fitted curve in Figure 4 (right panel, R^2 >= 0.77); they reproduce the paper's Table 1 stasis (50), EI10 (27) and EI5 (22) fAUC0-24/MIC targets to within 3%, while Table 1's EI50 (46) is 11% above the drawn curve's. The fungal density bact (CFU/mL) is integrated as d/dt(bact) = ln(10) * (change72 / 72) * bact so log10(bact) moves by exactly the fitted 72 h change. Companion files: Beredaki_2024_micafungin_rpmi_inoc1e4 (10^4 CFU/mL without serum) and Beredaki_2024_micafungin_serum_inoc1e4 / _serum_inoc1e3 (10% human serum)."
  reference <- paste(
    "Beredaki MI, Pournaras S, Meletiadis J. (2024).",
    "A new PK/PD target for assessing efficacy of micafungin against Candida",
    "parapsilosis. Journal of Antimicrobial Chemotherapy 79(1):157-165.",
    "doi:10.1093/jac/dkad360.",
    sep = " "
  )
  vignette <- "Beredaki_2024_micafungin"
  units <- list(
    time = "h",
    dosing = "mg (micafungin added to the 10 mL internal compartment)",
    concentration = "mg/L (micafungin in serum-free RPMI, all unbound, Cc); log10 CFU/mL (fungal density, log10cfu)"
  )

  # bact -- linear fungal density (CFU/mL) in the internal compartment.
  paper_specific_compartments <- c("bact")

  population <- list(
    species = "in vitro (Candida parapsilosis sensu stricto clinical isolates; two-compartment closed dialysis/diffusion PK/PD model)",
    n_subjects = 3L,
    n_studies = 1L,
    organism = paste(
      "Candida parapsilosis sensu stricto clinical isolates CP639, CP795 and",
      "CP732 (micafungin EUCAST/CLSI MIC 2, 1 and 0.5 mg/L). The",
      "pan-echinocandin-resistant fks1 F652S isolate (MIC > 8 mg/L) showed no",
      "effect at any exposure. Quality control: C. parapsilosis ATCC 22019 and",
      "Pichia kudriavzevii ATCC 6258."
    ),
    system = paste(
      "Previously described two-compartment closed diffusion/dialysis in vitro",
      "PK/PD model: an external compartment (conical flask on a peristaltic pump)",
      "and an internal compartment consisting of a 10 mL semipermeable cellulose",
      "dialysis tube (Spectra/Por Float-A-Lyzer G2, 300 kDa cut-off) inoculated",
      "with the yeast suspension. Medium: RPMI 1640 with L-glutamine, without",
      "bicarbonate, buffered to pH 7.0 with 0.165 M MOPS, with 100 mg/L",
      "chloramphenicol; no serum."
    ),
    disease_state = "In vitro infection of the internal compartment at a nominal starting inoculum of 10^3 CFU/mL (measured 3.06 +/- 0.12 log10 CFU/mL)",
    dose_range = paste(
      "Micafungin free peak concentrations fCmax 0.5, 1.8 and 4 mg/L q24h for",
      "72 h plus a drug-free control (Figure 2). Target half-life 15 h (human",
      "micafungin); achieved in vitro half-life 8 h (range 6-11 h)."
    ),
    sampling = paste(
      "PK: 100 uL from the internal compartment, assayed by microbiological",
      "diffusion against Aspergillus fumigatus AZN8196 (linear 0.06-8 mg/L in",
      "RPMI). PD: 200 uL at 0, 4, 8, 24, 48 and 72 h, plated on SGC2; limit of",
      "detection 10 CFU/mL."
    ),
    regions = "Greece (Attikon University Hospital, Athens)",
    notes = paste(
      "The drug-free control grew from 3.06 +/- 0.12 to 7.74 +/- 0.04",
      "log10 CFU/mL over 72 h. Exposure-response R^2 >= 0.77. Table 1 reports",
      "stasis at 50 (95% CI 42-65), EI50 46 (38-56), EI10 27 (18-41) and EI5 22",
      "(13-39) fAUC0-24/MIC, equivalent to 20 223, 18 400, 10 800 and 8800",
      "tAUC0-24/MIC under the paper's 99.75% protein binding."
    )
  )

  compartmentData <- list(
    central = list(
      analyte = "micafungin (unbound) in the internal compartment (serum-free RPMI-1640)",
      units = "mg",
      specimen = "not applicable",
      verified = TRUE
    ),
    bact = list(
      analyte = "Candida parapsilosis viable count in the internal compartment",
      units = "CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  ini({
    # ================================================================
    # In vitro pharmacokinetics
    # Beredaki 2024 Methods, "In vitro PK": exposures with Cmax 0.5-64 mg/L
    # and a target half-life of 15 h, "added at the corresponding Cmax
    # values in the in vitro model once daily". Figure 1(a) shows the
    # target profile returning to the same peak at 24 and 48 h, i.e. the
    # compartment is re-spiked to a constant peak.
    # ================================================================
    # Achieved half-life in serum-free RPMI, matching the system that
    # generated the pharmacodynamic data; the 15 h design target can be
    # simulated by overriding lkel.
    lkel <- log(log(2) / 8)
    label("Log first-order rate of micafungin removal from the internal compartment (1/h)") # Beredaki 2024 Results, "In vitro PK": t1/2 'of 8 (6-11) and 11 (7-18) h' in RPMI and RPMI + 10% human serum respectively

    # FIXED: a measured dimension of the apparatus. Converts the amount
    # added to the dialysis tube into a concentration.
    lvc <- fixed(log(0.010))
    label("Log volume of the internal compartment (L; apparatus dimension)") # Beredaki 2024 Methods, "In vitro PK/PD model": 'internal compartment (IC) of a 10 mL semipermeable cellulose dialysis tube' = 0.010 L

    # FIXED design input. Also the integration window of the PK/PD index.
    tau <- fixed(24)
    label("Dosing interval and PK/PD index integration window (h; experimental design input)") # Beredaki 2024 Methods: drug added 'once daily'; the index is AUC0-24/MIC

    # Target peak concentration reached at the start of every 24 h
    # interval. Serum-free RPMI binds no drug, so this is the free peak.
    # FIXED design input; change it to simulate another arm. Default
    # 1.8 mg/L is the middle arm.
    cmax <- fixed(1.8)
    label("Target peak micafungin concentration per 24 h interval (mg/L; experimental design input)") # Beredaki 2024 Figure 2 legend: fCmax 0.5, 1.8 and 4 mg/L q24h

    # ================================================================
    # Sigmoidal variable-slope Emax exposure-response, 72 h endpoint
    # Beredaki 2024 Methods, "PK/PD analysis":
    #   E = (Emax - Emin) * EI^n / (EI^n + EI50^n) + Emin
    # E is the 72 h change in log10 CFU/mL from the starting inoculum,
    # EI the AUC0-24/MIC index, Emin the effect at low exposure, Emax the
    # effect at high exposure. Written with Emax < Emin, the falling curve
    # has a positive exponent n.
    #
    # PROVENANCE: Beredaki 2024 prints Emin and n only for the serum
    # curves. Emin, EI50 and n below were recovered by tracing the fitted
    # 10^3 CFU/mL curve (grey) in Figure 4, right panel, from the 200 dpi
    # publisher PDF figure (249 points spanning fAUC0-24/MIC 0.03-63,
    # calibrated on the labelled decade and y-axis ticks) and
    # least-squares fitting the published equation with Emax held at its
    # constrained value; residual RMSE 0.037 log10 CFU/mL.
    # ================================================================
    e0 <- 4.40
    label("72 h change in log10 CFU/mL at zero exposure, Emin (log10 CFU/mL)") # digitised from Beredaki 2024 Figure 4 (right panel, 10^3 CFU/mL); drug-free control grew 3.06 -> 7.74 log10 CFU/mL (+4.68)
    emax <- fixed(-2)
    label("72 h change in log10 CFU/mL at saturating exposure, Emax (log10 CFU/mL; constrained)") # Beredaki 2024 Methods, "PK/PD analysis": the Emax was constrained '-2 log10 cfu/mL for the 10^3 cfu/mL'
    lec50 <- log(41.3)
    label("Log fAUC0-24/MIC producing 50% of the maximal effect, EI50 (unitless)") # digitised from Beredaki 2024 Figure 4 (right panel, 10^3 CFU/mL); Table 1 reports EI50 46 (95% CI 38-56) fAUC0-24/MIC
    lhill <- log(4.82)
    label("Log Hill coefficient n of the exposure-response relationship (unitless)") # digitised from Beredaki 2024 Figure 4 (right panel, 10^3 CFU/mL)

    # ================================================================
    # Minimum inhibitory concentration. FIXED: a measured property of the
    # isolate. EUCAST and CLSI MICs were identical for every isolate.
    # Default is CP639; change it to apply the model to another MIC.
    # ================================================================
    mic <- fixed(2)
    label("Micafungin EUCAST/CLSI MIC of the simulated isolate (mg/L; measured property of the isolate)") # Beredaki 2024 Methods, "Candida isolates": CP639, CP795, CP732 'with micafungin EUCAST/CLSI MICs of 2, 1 and 0.5 mg/L, respectively'

    # FIXED experimental design input: the measured mean starting density
    # in RPMI (nominal inoculum 10^3 CFU/mL).
    log10_cfu0 <- fixed(3.06)
    label("Log10 starting fungal density in the internal compartment (log10 CFU/mL)") # Beredaki 2024 Results, "In RPMI": 'at 10^3 cfu/mL as starting inoculum, C. parapsilosis grew from a mean +/- SD of 3.06 +/- 0.12 log10 cfu/mL at t = 0 h'

    # ================================================================
    # Residual error. Beredaki 2024 reports no assay CV for the
    # concentrations and only R^2 for the exposure-response fit, so both
    # are held at zero for deterministic typical-value simulation.
    # ================================================================
    propSd <- fixed(0)
    label("Proportional residual error on micafungin concentration (fraction; not reported)") # Beredaki 2024: no assay CV reported
    addSd_log10cfu <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (not reported)") # Beredaki 2024 Results: R^2 >= 0.77 without serum only, no residual SD
  })
  model({
    kel <- exp(lkel)
    vc <- exp(lvc)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # In vitro pharmacokinetics. Dose `central` in mg; the vignette
    # computes the amount that re-spikes the compartment to `cmax` at the
    # start of every interval.
    d/dt(central) <- -kel * central
    Cc <- central / vc

    # PK/PD index: closed-form AUC0-24 of a mono-exponential interval
    # starting at cmax. Serum-free medium, so this AUC is already the free
    # fAUC0-24. Beredaki 2024 fitted a 72 h ENDPOINT against the completed
    # AUC0-24, so the index must be available from t = 0.
    fauc24 <- cmax * (1 - exp(-kel * tau)) / kel
    ei <- fauc24 / mic

    # Beredaki 2024 exposure-response: 72 h change in log10 CFU/mL from
    # the starting inoculum, falling from e0 towards emax.
    change72 <- e0 - (e0 - emax) * ei^hill / (ei^hill + ec50^hill)

    # Spread the fitted 72 h change uniformly over the window so
    # log10(bact) moves by exactly change72 between t = 0 and t = 72 h.
    d/dt(bact) <- log(10) * (change72 / 72) * bact
    bact(0) <- 10^log10_cfu0

    # 1-CFU/mL floor keeps the log finite below 1 CFU/mL.
    log10cfu <- log10(bact + 1)

    Cc ~ prop(propSd)
    log10cfu ~ add(addSd_log10cfu)
  })
}
