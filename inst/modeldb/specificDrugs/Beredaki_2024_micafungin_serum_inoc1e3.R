Beredaki_2024_micafungin_serum_inoc1e3 <- function() {
  description <- "In vitro (Candida parapsilosis sensu stricto clinical isolates in RPMI-1640 with 10% pooled human serum, starting inoculum 10^3 CFU/mL; two-compartment closed dialysis/diffusion PK/PD model). Micafungin in vitro PK/PD model for the low inoculum that yields the paper's proposed new PK/PD target (EI10 = 28 tAUC0-24/MIC). Micafungin was added once daily to the 10 mL internal compartment at a target total peak (tCmax 2-64 mg/L) and decayed mono-exponentially with the achieved in vitro half-life of 11 h. The 72 h change in log10 CFU/mL from the starting inoculum was related to the PK/PD index tAUC0-24/MIC by the sigmoidal variable-slope Emax model E = (Emax - Emin) * EI^n / (EI^n + EI50^n) + Emin, with Emax constrained to -2 log10 CFU/mL for this inoculum (Methods, 'PK/PD analysis'); Emin = 5.23 and EI50 = 124 tAUC0-24/MIC are the published coefficients. The Hill coefficient is 1.45: the Results print n = -1.15, but that value contradicts both the paper's own Table 1 targets (stasis 240, EI10 28, EI5 16 tAUC0-24/MIC, each reproduced by n = 1.45 to within 3%) and the curve drawn in Figure 4 (digitised n = 1.48), so the printed value is treated as a typesetting error for -1.45. The fungal density bact (CFU/mL) is integrated as d/dt(bact) = ln(10) * (change72 / 72) * bact so log10(bact) moves by exactly the fitted 72 h change. Companion files: Beredaki_2024_micafungin_serum_inoc1e4 (10^4 CFU/mL with serum) and Beredaki_2024_micafungin_rpmi_inoc1e4 / _rpmi_inoc1e3 (no serum)."
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
    concentration = "mg/L (total micafungin in RPMI + 10% serum, Cc); log10 CFU/mL (fungal density, log10cfu)"
  )

  # bact -- linear fungal density (CFU/mL) in the internal compartment.
  paper_specific_compartments <- c("bact")

  population <- list(
    species = "in vitro (Candida parapsilosis sensu stricto clinical isolates; two-compartment closed dialysis/diffusion PK/PD model)",
    n_subjects = 2L,
    n_studies = 1L,
    organism = paste(
      "Candida parapsilosis sensu stricto clinical isolates CP639 and CP795",
      "(micafungin EUCAST/CLSI MIC 2 and 1 mg/L). The third susceptible isolate",
      "(CP732, MIC 0.5 mg/L) did not grow in 10% human serum and the",
      "pan-echinocandin-resistant fks1 F652S isolate (MIC > 8 mg/L) showed no",
      "effect at any exposure, so neither contributes to this curve. Quality",
      "control: C. parapsilosis ATCC 22019 and Pichia kudriavzevii ATCC 6258."
    ),
    system = paste(
      "Previously described two-compartment closed diffusion/dialysis in vitro",
      "PK/PD model: an external compartment (conical flask on a peristaltic pump)",
      "and an internal compartment consisting of a 10 mL semipermeable cellulose",
      "dialysis tube (Spectra/Por Float-A-Lyzer G2, 300 kDa cut-off) inoculated",
      "with the yeast suspension. Medium: RPMI 1640 with L-glutamine, without",
      "bicarbonate, buffered to pH 7.0 with 0.165 M MOPS, with 100 mg/L",
      "chloramphenicol and 10% heat-inactivated pooled human serum."
    ),
    disease_state = "In vitro infection of the internal compartment at a nominal starting inoculum of 10^3 CFU/mL (measured 3.2 +/- 0.25 log10 CFU/mL)",
    dose_range = paste(
      "Micafungin total peak concentrations tCmax 2, 8, 32 and 64 mg/L q24h",
      "for 72 h plus a drug-free control (Figure 3). Target half-life 15 h",
      "(human micafungin); achieved in vitro half-life 11 h (range 7-18 h)."
    ),
    sampling = paste(
      "PK: 100 uL from the internal compartment, assayed by microbiological",
      "diffusion against Aspergillus fumigatus AZN8196 (linear 0.125-64 mg/L",
      "with serum). PD: 200 uL at 0, 4, 8, 24, 48 and 72 h, plated on SGC2;",
      "limit of detection 10 CFU/mL."
    ),
    regions = "Greece (Attikon University Hospital, Athens)",
    notes = paste(
      "The drug-free control grew from 3.2 +/- 0.25 to 8.4 +/- 0.02",
      "log10 CFU/mL over 72 h. Exposure-response R^2 >= 0.90. Table 1 reports",
      "stasis at 240 (95% CI 176-352), EI50 124 (88-184), EI10 28 (16-64)",
      "and EI5 16 (8-48) tAUC0-24/MIC; fAUC values use 99.75% protein binding.",
      "EI10 = 28 tAUC0-24/MIC is the new PK/PD target the paper proposes: it",
      "gives >= 95% PTA for susceptible isolates (MIC <= 2 mg/L) at 100 mg q24h."
    )
  )

  compartmentData <- list(
    central = list(
      analyte = "micafungin (total) in the internal compartment (RPMI-1640 + 10% pooled human serum)",
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
    # values in the in vitro model once daily". Figure 1(b) shows the
    # target profile returning to the same peak at 24 and 48 h, i.e. the
    # compartment is re-spiked to a constant peak.
    # ================================================================
    # Achieved half-life in RPMI + 10% serum, matching the system that
    # generated the pharmacodynamic data; the 15 h design target can be
    # simulated by overriding lkel.
    lkel <- log(log(2) / 11)
    label("Log first-order rate of micafungin removal from the internal compartment (1/h)") # Beredaki 2024 Results, "In vitro PK": t1/2 'of 8 (6-11) and 11 (7-18) h' in RPMI and RPMI + 10% human serum respectively

    # FIXED: a measured dimension of the apparatus. Converts the amount
    # added to the dialysis tube into a concentration.
    lvc <- fixed(log(0.010))
    label("Log volume of the internal compartment (L; apparatus dimension)") # Beredaki 2024 Methods, "In vitro PK/PD model": 'internal compartment (IC) of a 10 mL semipermeable cellulose dialysis tube' = 0.010 L

    # FIXED design input. Also the integration window of the PK/PD index.
    tau <- fixed(24)
    label("Dosing interval and PK/PD index integration window (h; experimental design input)") # Beredaki 2024 Methods: drug added 'once daily'; the index is tAUC0-24/MIC

    # Target total peak concentration reached at the start of every 24 h
    # interval. FIXED design input; change it to simulate another arm.
    # Default 32 mg/L is the arm the paper reports as fungistatic for this
    # inoculum against both isolates.
    cmax <- fixed(32)
    label("Target peak total micafungin concentration per 24 h interval (mg/L; experimental design input)") # Beredaki 2024 Figure 3 legend: tCmax 2, 8, 32 and 64 mg/L q24h

    # ================================================================
    # Sigmoidal variable-slope Emax exposure-response, 72 h endpoint
    # Beredaki 2024 Methods, "PK/PD analysis":
    #   E = (Emax - Emin) * EI^n / (EI^n + EI50^n) + Emin
    # E is the 72 h change in log10 CFU/mL from the starting inoculum,
    # EI = tAUC0-24/MIC, Emin the effect at low exposure, Emax the effect
    # at high exposure. The paper reports the GraphPad slope as negative
    # (-1.15 as printed) because its curve falls with exposure; written in the form
    # above with Emax < Emin, the same curve has a positive exponent n.
    # ================================================================
    e0 <- 5.23
    label("72 h change in log10 CFU/mL at zero exposure, Emin (log10 CFU/mL)") # Beredaki 2024 Results, "PK/PD analysis": 'the Emin (4 versus 5.23) and n (-2.9 versus -1.15)' for 10^4 versus 10^3 CFU/mL with serum
    emax <- fixed(-2)
    label("72 h change in log10 CFU/mL at saturating exposure, Emax (log10 CFU/mL; constrained)") # Beredaki 2024 Methods, "PK/PD analysis": the Emax was constrained '-2 log10 cfu/mL for the 10^3 cfu/mL'
    lec50 <- log(124)
    label("Log tAUC0-24/MIC producing 50% of the maximal effect, EI50 (unitless)") # Beredaki 2024 Table 1 and Results: EI50 '124 ... tAUC0-24/MIC for the 10^3 cfu/mL' (95% CI 88-184) with 10% serum

    # Results print n = -1.15 for this curve. That value does not
    # reproduce the paper's own Table 1 targets: with Emin 5.23, Emax -2
    # and EI50 124 it gives stasis 286, EI10 18 and EI5 10 tAUC0-24/MIC
    # against the reported 240, 28 and 16. n = 1.45 gives 241, 27 and 16,
    # and a least-squares trace of the Figure 4 (left panel, grey) curve
    # gives n = 1.48 (Emin 5.24, EI50 0.314 fAUC0-24/MIC = 126
    # tAUC0-24/MIC). The printed -1.15 is therefore read as a typesetting
    # error for -1.45. See the vignette "Assumptions and deviations".
    lhill <- log(1.45)
    label("Log Hill coefficient n of the exposure-response relationship (unitless)") # Beredaki 2024 Table 1 stasis/EI10/EI5 (240, 28, 16 tAUC0-24/MIC) and Figure 4 curve; Results print -1.15, read as -1.45

    # ================================================================
    # Minimum inhibitory concentration. FIXED: a measured property of the
    # isolate. EUCAST and CLSI MICs were identical for every isolate.
    # Default is CP639; change it to apply the model to another MIC.
    # ================================================================
    mic <- fixed(2)
    label("Micafungin EUCAST/CLSI MIC of the simulated isolate (mg/L; measured property of the isolate)") # Beredaki 2024 Methods, "Candida isolates": CP639, CP795, CP732 'with micafungin EUCAST/CLSI MICs of 2, 1 and 0.5 mg/L, respectively'

    # FIXED experimental design input: the measured mean starting density
    # with serum (nominal inoculum 10^3 CFU/mL).
    log10_cfu0 <- fixed(3.2)
    label("Log10 starting fungal density in the internal compartment (log10 CFU/mL)") # Beredaki 2024 Results, "In RPMI with 10% human serum": 'fungal burden increased from 3.2 +/- 0.25 log10 cfu/mL at t = 0 h'

    # ================================================================
    # Residual error. Beredaki 2024 reports no assay CV for the
    # concentrations and only R^2 for the exposure-response fit, so both
    # are held at zero for deterministic typical-value simulation.
    # ================================================================
    propSd <- fixed(0)
    label("Proportional residual error on micafungin concentration (fraction; not reported)") # Beredaki 2024: no assay CV reported
    addSd_log10cfu <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (not reported)") # Beredaki 2024 Results: R^2 >= 0.90 with serum only, no residual SD
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

    # PK/PD index: closed-form total AUC0-24 of a mono-exponential interval
    # starting at cmax. Beredaki 2024 fitted a 72 h ENDPOINT against the
    # completed AUC0-24, so the index must be available from t = 0.
    tauc24 <- cmax * (1 - exp(-kel * tau)) / kel
    ei <- tauc24 / mic

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
