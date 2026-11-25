Zhang_2019_doxycycline_killrate <- function() {
  description <- "In vitro (Mycoplasma hyopneumoniae strain ATCC25934). Concentration-dependent kill-rate PK/PD model for doxycycline, coupling the one-compartment first-order-absorption PK of Zhang 2019's in-vitro dynamic system to the sigmoid Emax kill-rate relationship fitted to the static time-kill curves. Zhang 2019 Materials and Methods (Static time-kill curves fitting and analysis) parameterises the kill rate as E = Emax * Ce^N / (EC50^N + Ce^N), where E is the kill rate in log10(CFU/mL) per hour (POSITIVE = killing) obtained by linear regression of the static time-kill counts over a time window, Ce is the doxycycline concentration and there is no baseline term, so the kill rate is zero in the absence of drug. Parameters are Zhang 2019 Table 1, row 12-48 h, the window the paper identifies as giving the highest correlation (R2 = 0.986): Emax = 0.06 1/h, EC50 = 0.78 ug/mL, Hill N = 2.98. The bacterial density bacteria (linear CFU/mL) is integrated as d/dt(bacteria) = -ln(10) * E * bacteria so that log10(bacteria) falls at exactly E log10(CFU/mL) per hour, reproducing the fitted kill rate at any constant concentration. The PK is the paper's in-vitro dynamic model (Materials and Methods, Description of the in vitro dynamic model): a 210 mL absorption chamber draining into a 300 mL central chamber at 0.58 mL/min, designed to reproduce the porcine plasma absorption half-life of 4.15 h and elimination half-life of 5.9 h (Bousquet 1998), giving ka = log(2)/4.15 = 0.1670 1/h and kel = log(2)/5.9 = 0.1175 1/h. The pump geometry gives 0.58*60/210 = 0.166 1/h and 0.58*60/300 = 0.116 1/h, agreeing with the design half-lives to within 2%. IMPORTANT PROVENANCE NOTE: Zhang 2019 fitted the kill-rate parameters to STATIC concentrations and ran the dynamic apparatus as a separate experiment; the two were not fitted jointly, and the static kill-rate fit does not represent the growth of the untreated control. This model composes them, which is the paper's stated PK/PD integration written as an ODE system. The paper does not state the amount of doxycycline placed in the absorption chamber for each of its mg/kg-labelled regimens, so doses are entered as an amount in ug. Zhang 2019 reports no between-replicate variability and no residual error magnitude for either component, so no eta parameters are present and both residual SDs are held at zero for typical-value simulation."
  reference <- paste(
    "Zhang H, Mao C, Li J, Huang Z, Gu X, Shen X, Ding H. (2019).",
    "Pharmacokinetic/pharmacodynamic integration of doxycycline against",
    "Mycoplasma hyopneumoniae in an in vitro model.",
    "Frontiers in Pharmacology 10:1088.",
    "doi:10.3389/fphar.2019.01088. PMCID: PMC6763577.",
    sep = " "
  )
  vignette <- "Zhang_2019_doxycycline"
  units <- list(
    time = "h",
    dosing = "ug (doxycycline into the 210 mL absorption chamber, depot)",
    concentration = "ug/mL (Cc, central chamber); log10 CFU/mL (log_cfu)"
  )

  compartmentData <- list(
    depot = list(analyte = "doxycycline", units = "ug", specimen = "not applicable", verified = TRUE),
    central = list(analyte = "doxycycline", units = "ug", specimen = "not applicable", verified = TRUE),
    bacteria = list(
      analyte = "Mycoplasma hyopneumoniae strain ATCC25934",
      units = "CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "in vitro (Mycoplasma hyopneumoniae strain ATCC25934 culture)",
    n_subjects = 3L,
    n_studies = 1L,
    organism = "Mycoplasma hyopneumoniae standard strain ATCC25934 (Chinese Veterinary Microorganism Culture Collection Center, Beijing), cultured at 37 C and 5% CO2 in artificial medium base with 20% pig serum, cysteine and NADH. Doxycycline MIC = 0.125 ug/mL by broth microdilution at a 10^6 CFU/mL inoculum",
    system = "Static time-kill curves in 7 mL penicillin bottles (0.4 mL exponential-phase suspension + 0.1 mL doxycycline solution + 3.5 mL medium, about 10^6 CFU/mL) for the kill-rate parameters; the PK component is the separate in-vitro dynamic apparatus, a 210 mL absorption chamber feeding a central chamber of 300 mL medium (external compartment) holding a 10 mL dialysis tube (internal compartment) with the culture, supplied from a fresh-medium reservoir and draining to a waste vessel by peristaltic pumps at 0.58 mL/min",
    disease_state = "Not applicable (in-vitro culture)",
    dose_range = "Static time-kill: 0, 0.5, 1, 2, 4, 6, 8, 16, 32 and 64 times MIC (0 to 8 ug/mL). In-vitro dynamic apparatus: eight regimens labelled 1, 2.5, 5, 7.5, 10, 15, 20 and 30 mg/kg body weight, once daily for 3 days into the absorption chamber",
    sampling = "Static: 100 uL sampled at 0, 1, 3, 6, 9, 12, 24, 36 and 48 h for colony counting, limit of detection 100 CFU/mL, each experiment repeated at least three times",
    regions = "China (South China Agricultural University, Guangzhou)",
    notes = "Zhang 2019 Table 1 reports the sigmoid Emax kill-rate fit over seven time windows (0-24, 0-36, 0-48, 3-24, 3-36, 3-48 and 12-48 h) with R2 from 0.947 to 0.986. This model packages the 12-48 h window, which the paper identifies as the best fit; the other six parameterisations are tabulated in the vignette. The kill-rate model has no growth term, so the untreated culture is predicted to stay at its starting density."
  )

  ini({
    # Zhang 2019 Materials and Methods (Static time-kill curves fitting and
    # analysis): E = Emax * Ce^N / (EC50^N + Ce^N), E = kill rate in
    # log10(CFU/mL) per hour, positive = killing. No baseline term.
    lemax <- log(0.06)
    label("Maximum kill rate Emax (log10 CFU/mL per h)") # Zhang 2019 Table 1, row 12-48 h, Emax = 0.06 h-1
    lec50 <- log(0.78)
    label("Doxycycline concentration producing half the maximum kill rate EC50 (ug/mL)") # Zhang 2019 Table 1, row 12-48 h, EC50 = 0.78 ug/ml
    lhill <- log(2.98)
    label("Hill coefficient N defining the slope of the kill-rate curve (unitless)") # Zhang 2019 Table 1, row 12-48 h, Hill's slope = 2.98

    # In-vitro dynamic-system PK: design constants of the apparatus, set from
    # the porcine plasma doxycycline half-lives the device was built to mimic.
    lka <- fixed(log(log(2) / 4.15))
    label("First-order transfer rate from the absorption chamber to the central chamber (1/h)") # Zhang 2019 Description of the in vitro dynamic model: absorption half-life of 4.15 h
    lkel <- fixed(log(log(2) / 5.9))
    label("First-order elimination rate constant of the central chamber (1/h)") # Zhang 2019 Description of the in vitro dynamic model: elimination time of 5.9 h
    lvc <- fixed(log(300))
    label("Central-chamber medium volume (mL)") # Zhang 2019 Description of the in vitro dynamic model: central chamber comprised 300 ml of sterile medium (external compartment)

    # Experimental design input, not an estimated parameter.
    log10_cfu0 <- fixed(7)
    label("log10 starting mycoplasma density (log10 CFU/mL)") # Zhang 2019 Exposure to dynamic antibiotic concentrations: inoculum of 10^7 CFU/ml in the dialysis tube

    # Zhang 2019 reports R2 only and no residual SD for either component.
    addSd <- fixed(0)
    label("Additive residual SD on doxycycline concentration (ug/mL; not reported in Zhang 2019)") # Zhang 2019 reports R2 only, no residual SD
    addSd_log_cfu <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (not reported in Zhang 2019)") # Zhang 2019 reports R2 only, no residual SD
  })

  model({
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)
    ka <- exp(lka)
    kel <- exp(lkel)
    vc <- exp(lvc)

    # Absorption chamber (depot) drains into the central chamber; the central
    # chamber is diluted by fresh medium at the same flow.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    Cc <- central / vc

    # Zhang 2019 Table 1 sigmoid Emax kill rate (log10 units per hour).
    kill_log10 <- emax * Cc^hill / (ec50^hill + Cc^hill)

    # d(log10 N)/dt = -kill_log10  =>  d(N)/dt = -ln(10)*kill_log10*N
    d/dt(bacteria) <- -log(10) * kill_log10 * bacteria
    bacteria(0) <- 10^log10_cfu0

    # 1-CFU/mL floor keeps the log10 finite if bacteria is driven below 1.
    log_cfu <- log10(bacteria + 1)

    Cc ~ add(addSd)
    log_cfu ~ add(addSd_log_cfu)
  })
}
