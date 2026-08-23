Gao_2025_cefquinome_pkpd_index <- function() {
  description <- paste0(
    "Ex vivo (foal serum). Inhibitory sigmoid Emax PK/PD-integration model for the antibacterial ",
    "effect of cefquinome against Escherichia coli strain HE13, isolated from a septicaemic foal, ",
    "in serum drawn from Ili foals after a single 1 mg/kg intramuscular dose. Gao 2025 Section 2.9 ",
    "parameterises the effect over a 24 h ex vivo incubation as ",
    "E = Emax - (Emax - E0) * Ce^N / (Ce^N + EC50^N), where E is the SIGNED change in log10(CFU/mL) ",
    "between 0 and 24 h of incubation (positive = net growth, negative = net kill) and Ce is the ",
    "PK/PD index. NOTE THE INVERTED NAMING: Gao 2025 defines its Emax as the 24 h change in the ",
    "DRUG-FREE samples and its E0 as the change in the drug-treated samples, which is the opposite ",
    "of the usual convention and of the sibling veterinary models in this library. This file ",
    "therefore maps Gao 2025's Emax = 2.88 log10 CFU/mL onto the canonical `e0` (effect at zero ",
    "exposure) and Gao 2025's E0 = -4.65 log10 CFU/mL onto the canonical `emax` (maximal drug ",
    "effect), and writes the algebraically identical form ",
    "E = e0 - (e0 - emax) * Ce^N / (Ce^N + EC50^N). No value is altered by the remapping. ",
    "Parameters from Gao 2025 Table 4: e0 = 2.88, emax = -4.65 log10 CFU/mL, EC50 = 2.61 ",
    "(dimensionless), Hill N = 4.22. The PK/PD index is AUC0-24h/MIC DIVIDED BY 24 h, which Gao ",
    "2025 adopts so the index is dimensionless; it is formed here as the per-interval covariate ",
    "AUC_CEFQ divided by the parameter mic and by 24. The parameterisation was confirmed against ",
    "the paper's own independently printed targets: solving E = 0 and E = -3 returns index values ",
    "of 2.330 and 3.527 against the published bacteriostatic and bactericidal targets of 2.34 and ",
    "3.53. The published bacterial-elimination target of 4.86 does NOT correspond to E = -4, which ",
    "the fitted curve places at 4.565; the published value instead corresponds to E = -4.14. See ",
    "the vignette Errata. The bacterial density bact (linear CFU/mL) is integrated as ",
    "d/dt(bact) = ln(10) * (E / 24) * bact so that log10(bact) changes by exactly E across each ",
    "24 h window, reproducing the paper's model at the only times bacteria were counted. There is ",
    "no PK component in THIS file: Gao 2025 analysed the serum concentrations non-compartmentally ",
    "in WinNonlin 5.2.1 and defined its PD per 24 h interval on that non-compartmental AUC, so ",
    "exposure enters as the externally supplied covariate AUC_CEFQ. The companion model ",
    "`Gao_2025_cefquinome_foal` supplies a compartmental PK model for the same study if a closed ",
    "PK-PD loop is wanted. Neither between-subject variability nor a residual error magnitude was ",
    "reported, so there are no eta parameters and addSd is FIXED at 0; the model is intended for ",
    "typical-value simulation."
  )
  reference <- paste(
    "Gao T, Liu X, Qiu D, Li Y, Qiu Z, Qi J, Li S, Guo X, Zhang Y, Wang Z, Gao X, Ma Y, Ma T. (2025).",
    "Ex vivo pharmacokinetic/pharmacodynamic integration model of cefquinome",
    "against Escherichia coli in foals.",
    "Veterinary Sciences 12(4):294.",
    "doi:10.3390/vetsci12040294.",
    sep = " "
  )
  vignette <- "Gao_2025_cefquinome"
  units <- list(
    time = "h",
    dosing = "h*ug/mL (cefquinome AUC0-24h per dosing interval, supplied as a covariate)",
    concentration = "log10 CFU/mL (observation)"
  )

  depends <- c("AUC_CEFQ")
  paper_specific_compartments <- c("bact")

  compartmentData <- list(
    bact = list(analyte = "Escherichia coli strain HE13", units = "CFU/mL", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    AUC_CEFQ = list(
      description        = "Cefquinome area under the serum concentration-time curve over the current 24 h dosing interval",
      units              = "h*ug/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-24-h-interval cefquinome serum exposure, supplied as a piecewise-constant time-varying covariate (one value per 24 h interval). Gao 2025 derived it non-compartmentally in WinNonlin 5.2.1 (Section 2.4) and reports AUC0-last as 5.41 +/- 0.81 h*ug/mL after a single 1 mg/kg intramuscular dose and 12.33 +/- 0.69 h*ug/mL after 1 mg/kg intravenously (Table 1). Because the last quantifiable sample is at 12 h (Figure 2) and the terminal half-life is 2.35-4.16 h, AUC0-last and AUC0-24h are numerically interchangeable in this study, and the intramuscular value of 5.41 h*ug/mL is the exposure that generated the ex vivo serum samples. Set to 0 for the drug-free control so the sigmoid term vanishes and the predicted 24 h change equals e0. The companion model Gao_2025_cefquinome_foal can generate this covariate for an arbitrary dose if a closed PK-PD loop is wanted.",
      source_name        = "AUC0-24h (Gao 2025 Sections 2.9-2.11; values in Table 1 as AUC0-last)"
    )
  )

  population <- list(
    species        = "ex vivo (serum from horse, Ili foal) challenged with Escherichia coli",
    n_subjects     = 10L,
    n_studies      = 1L,
    age_range      = "7 months to 1 year (serum donors)",
    weight_range   = "191 +/- 21.7 kg (serum donors, mean +/- SD)",
    organism       = "Escherichia coli strain HE13, isolated from a foal with septicaemia in Heilongjiang Province, China. Cefquinome MIC 0.125 ug/mL in Mueller Hinton broth and 0.062 ug/mL in foal serum; E. coli ATCC 25922 was the quality-control strain (0.054 ug/mL broth, 0.031 ug/mL serum)",
    system         = "Ex vivo time-kill. Serum was drawn from foals before dosing and at 0.083, 0.167, 0.25, 0.5, 0.75, 1, 2, 3, 6, 9, 12 and 24 h after a single 1 mg/kg intramuscular dose; 0.5 mL of each serum sample was mixed with 50 uL of stationary-phase culture to about 1 x 10^6 CFU/mL, and 20 uL aliquots were serially diluted onto trypticase soy agar after 1, 3, 6, 9, 12 and 24 h of incubation (detection limit 200 CFU/mL)",
    disease_state  = "Not applicable (ex vivo challenge of serum from healthy foals); starting inoculum about 1 x 10^6 CFU/mL",
    dose_range     = "Single 1 mg/kg intramuscular dose of cefquinome sulphate injection in the serum donors",
    regions        = "China (Northeast Agricultural University, Harbin; foals in Zhaosu County, Yili)",
    notes          = "Ethics approval XY20230322. Gao 2025 measured MICs for 53 E. coli isolates (Figure 3: 4 strains at 0.016, 9 at 0.031, 20 at 0.062, 9 at 0.125, 7 at 0.25, 2 at 0.5 and 2 at 1 ug/mL in serum), giving MIC50 = 0.062 ug/mL; MICs in Mueller Hinton broth averaged 1.8-fold higher than in serum (Table 2), which is why the serum MIC is used here. Post-antibiotic effects were short and concentration-dependent (Table 3: 0.13-0.58 h over 1x to 4x MIC), so no PAE term is included. Gao 2025 also ran a 10,000-subject Monte Carlo simulation in Crystal Ball, treating AUC0-24h as log-normal and sampling the MIC distribution, and reported target attainment rates of 46.38% for the bactericidal target and 24.48% for the bacterial-elimination target at the recommended 1 mg/kg intramuscular dose; the vignette reproduces both to within about 1 percentage point from this model plus the Figure 3 MIC distribution."
  )

  ini({
    # =================================================================
    # Gao 2025 Table 4 -- inhibitory sigmoid Emax PK/PD-integration model
    # =================================================================
    # Gao 2025 Section 2.9 equation, as printed:
    #   E = Emax - (Emax - E0) * Ce^N / (Ce^N + EC50^N)
    # with Gao 2025's own definitions (Section 2.9 text): its "Emax" is
    # the 24 h change in bacterial count in the samples that received NO
    # drug, and its "E0" is the change in the drug-treated samples. That
    # is the reverse of the usual convention, and the Table 4 signs
    # confirm it: the no-drug value is +2.88 (net growth) and the
    # drug-saturated value is -4.65 (net kill).
    #
    # This file uses the canonical roles, so:
    #   e0   <- Gao 2025's Emax = +2.88  (effect at zero exposure)
    #   emax <- Gao 2025's E0   = -4.65  (maximal drug effect)
    # and model() writes the algebraically identical
    #   E = e0 - (e0 - emax) * Ce^N / (Ce^N + EC50^N)
    # which runs from e0 at Ce = 0 to emax as Ce grows. No value is
    # changed by the remapping -- only the label each number carries.
    #
    # E is the SIGNED change in log10(CFU/mL) over a 24 h incubation, so
    # emax is negative and neither parameter can be log-transformed;
    # both stay on the natural scale.
    #
    # Numerical confirmation against values Gao 2025 printed
    # independently of Table 4: solving E = 0 (bacteriostatic) returns an
    # index of 2.330 against the published 2.34, and solving E = -3
    # (bactericidal) returns 3.527 against the published 3.53.
    e0 <- 2.88
    label("Change in E. coli count over 24 h in drug-free serum (log10 CFU/mL)")  # Gao 2025 Table 4, "Emax" = 2.88 Log10 CFU/mL (defined in Section 2.9 as the no-drug 24 h change)
    emax <- -4.65
    label("Maximal change in E. coli count over 24 h (log10 CFU/mL; negative = reduction)")  # Gao 2025 Table 4, "E0" = -4.65 Log10 CFU/mL (defined in Section 2.9 as the drug-treated 24 h change)
    lec50 <- log(2.61)
    label("PK/PD index producing 50% of the maximal effect EC50 (dimensionless)")  # Gao 2025 Table 4, EC50 = 2.61
    lhill <- log(4.22)
    label("Hill coefficient N defining the steepness of the fitted curve (unitless)")  # Gao 2025 Table 4, N = 4.22

    # =================================================================
    # Minimum inhibitory concentration
    # =================================================================
    # The index driving the sigmoid is AUC0-24h/MIC divided by 24 h, so
    # the MIC is required numerically. FIXED because it is a measured
    # property of the challenge strain, not an estimated parameter. The
    # SERUM MIC is used, not the Mueller Hinton broth MIC of 0.125
    # ug/mL: Gao 2025 built the ex vivo model in serum and argues
    # (Discussion) that broth MICs are artificially elevated, averaging
    # 1.8-fold higher across the 20 strains in Table 2. Change it to
    # apply the model to an isolate with a different susceptibility.
    mic <- fixed(0.062)
    label("Cefquinome MIC against E. coli HE13 in foal serum (ug/mL)")  # Gao 2025 Section 3.3: strain HE13 MIC = 0.062 ug/mL in foal serum (0.125 ug/mL in MHB); also the MIC50 of the 53-strain collection

    # =================================================================
    # Starting bacterial density
    # =================================================================
    # Experimental design input, not an estimated parameter.
    log10_cfu0 <- fixed(6)
    label("log10 starting E. coli density in the ex vivo incubation (log10 CFU/mL)")  # Gao 2025 Section 2.8: serum plus stationary-phase culture to a final concentration of about 1 x 10^6 CFU/mL

    # =================================================================
    # Residual error
    # =================================================================
    # Gao 2025 reported only the correlation coefficient of the sigmoid
    # fit (Figure 6) and gave no residual standard deviation, so the
    # residual SD is held at zero for deterministic typical-value
    # simulation.
    addSd <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (not reported in Gao 2025)")  # Gao 2025 reports R^2 only, no residual SD
  })

  model({
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # PK/PD index: cefquinome serum AUC over the current 24 h interval
    # divided by the MIC and then by 24 h (Gao 2025 Section 2.9, which
    # divides AUC0-24h/MIC by 24 h so the index is dimensionless and can
    # be substituted directly into the paper's dose equation).
    ce <- AUC_CEFQ / mic / 24

    # Gao 2025 inhibitory sigmoid Emax equation, written with canonical
    # parameter roles (see the ini() note on the inverted naming).
    # eff_log10 is the SIGNED change in log10(CFU/mL) accrued over the
    # current 24 h interval (negative = bacterial reduction); it equals
    # e0 at zero exposure and approaches emax at saturating exposure.
    eff_log10 <- e0 - (e0 - emax) * ce^hill / (ce^hill + ec50^hill)

    # Gao 2025 fitted only the per-interval change and counted bacteria
    # at the ends of 24 h incubations. Spreading the interval change
    # uniformly across the interval makes log10(bact) change by exactly
    # eff_log10 over each 24 h window, so the trajectory matches the
    # paper's model at every time the paper actually counted:
    #   d(log10 N)/dt = eff_log10 / 24
    # => d(N)/dt      = ln(10) * (eff_log10 / 24) * N
    d/dt(bact) <- log(10) * (eff_log10 / 24) * bact
    bact(0) <- 10^log10_cfu0

    # log10 CFU/mL observation with a 1-CFU/mL floor so the log10 stays
    # finite if bact is driven below 1 CFU/mL.
    Cc <- log10(bact + 1)
    Cc ~ add(addSd)
  })
}
