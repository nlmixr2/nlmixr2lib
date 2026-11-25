Li_2017_enrofloxacin_chicken_infected <- function() {
  description <- paste0(
    "Preclinical (chicken, Cobb broiler, Salmonella Typhimurium-infected). ",
    "Ex vivo inhibitory sigmoid Emax PK/PD-integration model for the antibacterial effect of oral ",
    "enrofloxacin against Salmonella Typhimurium CVCC541 in intestinal contents collected from ",
    "INFECTED chickens after a single 10 mg/kg oral dose (Li 2017 Table 2, 'Infected' column; ",
    "Supplementary Figure S3b). This is the fit Li 2017 used to design its 4 mg/kg PK/PD dose. The ",
    "effect E is the signed change in log10 CFU/mL over a 24 h ex vivo incubation, driven by the ",
    "PK/PD index AUC24h/MIC: E = e0 - (e0 - emax) * idx^N / (idx^N + EC50^N), with e0 = +2.14 ",
    "(drug-free 24 h growth), emax = -7.01 (maximal kill), EC50 = 684.21 h and N = 4.96. Li 2017 ",
    "labels these two plateaus the other way round ('Log Emax' = 2.14, 'Log E0' = -7.01); the roles ",
    "are fixed by Supplementary Figure S3b and by the paper's own bacteriostatic / bactericidal / ",
    "eradication targets, which the inverted equation reproduces to the printed precision (538.64, ",
    "719.33 and 789.97 h). Exposure enters as the intestinal-content AUC over 24 h (covariate ",
    "AUC_ENRO) divided by the intestinal-content MIC of the challenge strain (0.12 ug/mL), which ",
    "reconstructs the paper's Supplementary Table S1 AUC24h/MIC of 2932 h from the Table 1 AUC0-24 ",
    "of 351.8 h*ug/mL exactly. The bacterial density is integrated as ",
    "d/dt(bact) = ln(10) * (E / 24) * bact so log10(bact) changes by exactly E over each 24 h window. ",
    "There is NO PK component: Li 2017 analysed the concentrations non-compartmentally in WinNonlin ",
    "5.2.1 and printed no structural PK model. No between-animal variability or residual error was ",
    "reported (points are means of five chickens), so the model has no etas and addSd is fixed at 0."
  )
  reference <- paste(
    "Li J, Hao H, Cheng G, Wang X, Ahmed S, Shabbir MAB, Liu Z, Dai M, Yuan Z. (2017).",
    "The effects of different enrofloxacin dosages on clinical efficacy and resistance",
    "development in chickens experimentally infected with Salmonella Typhimurium.",
    "Scientific Reports 7:11676. doi:10.1038/s41598-017-12294-7.",
    sep = " "
  )
  vignette <- "Li_2017_enrofloxacin_chicken"
  units <- list(
    time = "h",
    dosing = "h*ug/mL (enrofloxacin intestinal-content AUC over 24 h, supplied as a covariate)",
    concentration = "log10 CFU/mL (observation)"
  )

  depends <- c("AUC_ENRO")
  paper_specific_compartments <- c("bact")

  compartmentData <- list(
    bact = list(
      analyte = "Salmonella Typhimurium CVCC541",
      units = "CFU/mL",
      specimen = "faeces",
      verified = TRUE
    )
  )

  covariateData <- list(
    AUC_ENRO = list(
      description = "Enrofloxacin area under the intestinal-content concentration-time curve over 24 h (AUC24h)",
      units = "h*ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Divided by the model's mic parameter (0.12 ug/mL) to form the PK/PD index AUC24h/MIC (h) ",
        "that drives the sigmoid. Li 2017 Table 1 reports the intestinal-content AUC0-24 after a ",
        "single 10 mg/kg oral dose as 351.8 h*ug/mL in infected chickens (188.6 in healthy chickens); ",
        "351.8 / 0.12 = 2931.7 h reproduces the Supplementary Table S1 AUC24h/MIC of 2932 h. Li 2017 ",
        "reports intestinal concentrations in ug/mL in Table 1 and ug/g in the Discussion; the two are ",
        "treated as equal. Set to 0 for the drug-free control so the sigmoid term vanishes and the ",
        "predicted 24 h change equals e0."
      ),
      source_name = "AUC24h (Li 2017 Methods 'Design of dosage by PK/PD modeling'; Table 1 row 'AUC0-24'; Supplementary Table S1 row 'AUC24h/MIC')"
    )
  )

  population <- list(
    species = "chicken (Cobb broiler)",
    n_subjects = 65L,
    n_studies = 1L,
    age_range = "30 days at enrolment, inoculated after a 7-day acclimatisation and dosed 24 h later",
    weight_range = "1.0-1.5 kg",
    disease_state = "Experimentally infected with 10^8 CFU Salmonella Typhimurium CVCC541 by gavage 24 h before dosing; the PK/PD group A of Li 2017",
    dose_range = "10 mg/kg enrofloxacin (Baytril 10% oral solution) as a single oral bolus",
    organism = "Salmonella Typhimurium CVCC541 (Chinese Veterinary Culture Collection); enrofloxacin MIC 0.12 ug/mL and MBC 0.25 ug/mL in intestinal contents, 0.06 and 0.12 ug/mL in Mueller-Hinton broth; MPC 0.8 ug/mL in MH broth",
    system = "Ex vivo time-kill assay in intestinal contents collected before dosing and at 0.25-48 h after dosing, counted at 0, 1, 2, 4, 8, 12 and 24 h of incubation (Supplementary Figure S2b); each point is the mean of five chickens",
    regions = "China (Huazhong Agricultural University, Wuhan, Hubei)",
    notes = paste0(
      "130 broilers were split equally into an infected group A (10^8 CFU S. Typhimurium CVCC541 by ",
      "gavage, dosed 24 h later) and a healthy group B; five chickens per group were killed at each ",
      "of 13 sampling times. Infected intestinal-content NCA (Table 1): Cmax 48.99 ug/mL, Tmax 0.5 h, ",
      "T1/2 8.57 h, AUC0-24 351.8 h*ug/mL, ClB/F 0.022. Unbound (feces-free) fraction 0.47. The ",
      "bactericidal target of this fit (719.33 h) gave the 4 mg/kg PK/PD dose tested in the clinical ",
      "efficacy study. The healthy-chicken fit is packaged separately as ",
      "Li_2017_enrofloxacin_chicken_healthy."
    )
  )

  ini({
    # ===============================================================
    # Li 2017 Table 2 -- inhibitory sigmoid Emax PK/PD integration,
    # 'Infected' column (Supplementary Figure S3b)
    # ===============================================================
    # Li 2017 prints 'Log E0' = -7.01 and 'Log Emax' = 2.14. Its footnote
    # calls E0 the drug-free control change, but the drug-free control
    # GREW (Supplementary Figure S2b, red curve, about 7.0 -> 9.2 log10
    # CFU/mL) and Supplementary Figure S3b places +2.14 at AUC24h/MIC = 0
    # and -7.01 on the high-exposure plateau. So the values carry the
    # reverse of the usual labels, as in Gao 2025 and Wei 2026. This file
    # uses canonical roles:
    #   e0   <- Li 2017 'Log Emax' = +2.14 (effect at zero exposure)
    #   emax <- Li 2017 'Log E0'   = -7.01 (maximal drug effect)
    # Only the labels move. Inverting
    #   E = e0 - (e0 - emax) * idx^N / (idx^N + EC50^N)
    # at E = 0, -3 and -4 returns 538.64, 719.33 and 789.97 h, the three
    # targets Table 2 prints, to two decimals. E is signed, so e0 and
    # emax stay on the natural scale.
    e0 <- 2.14
    label("Change in S. Typhimurium count over 24 h in drug-free intestinal contents (log10 CFU/mL)") # Li 2017 Table 2, Infected, 'Log Emax' = 2.14 CFU/mL
    emax <- -7.01
    label("Maximal change in S. Typhimurium count over 24 h (log10 CFU/mL; negative = kill)") # Li 2017 Table 2, Infected, 'Log E0' = -7.01 CFU/mL
    lec50 <- log(684.21)
    label("AUC24h/MIC producing half the maximal effect EC50 (h)") # Li 2017 Table 2, Infected, 'AUC24/MIC EC50' = 684.21 h
    lhill <- log(4.96)
    label("Hill coefficient N, steepness of the AUC24h/MIC effect curve (unitless)") # Li 2017 Table 2, Infected, 'Slope (N)' = 4.96

    # ===============================================================
    # Minimum inhibitory concentration (index denominator)
    # ===============================================================
    # Measured property of the challenge strain, not an estimate. The
    # intestinal-content MIC is used, not the MH-broth MIC (0.06): it is
    # the value that reconstructs Supplementary Table S1 (351.8 / 0.12 =
    # 2932 h; Cmax/MIC 48.99 / 0.12 = 408.3). Change it to apply the
    # model to an isolate of different susceptibility.
    mic <- fixed(0.12)
    label("Enrofloxacin MIC against S. Typhimurium CVCC541 in chicken intestinal contents (ug/mL)") # Li 2017 Results: MIC in intestinal contents 0.12 mg/L (0.06 mg/L in MH broth)

    # ===============================================================
    # Starting bacterial density (experimental design input)
    # ===============================================================
    # Not printed in the text; digitised from the common t = 0 point of
    # every curve in Supplementary Figure S2b. It equals -emax to the
    # digitisation precision: the maximal-kill plateau corresponds to
    # driving the count to 1 CFU/mL (log10 = 0).
    log10_cfu0 <- fixed(7.01)
    label("log10 starting S. Typhimurium density in the ex vivo incubation (log10 CFU/mL)") # Li 2017 Supplementary Figure S2b, t = 0 point, digitised by the maintainers (7.01)

    # ===============================================================
    # Residual error
    # ===============================================================
    # Not reported (the sigmoid was fitted to means of five chickens), so
    # held at 0 for deterministic typical-value simulation.
    addSd <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (not reported in Li 2017)") # Li 2017 reports no residual error for the sigmoid fit
  })

  model({
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # PK/PD index AUC24h/MIC (h), as tabulated in Li 2017 Table 2 and
    # Supplementary Table S1 (no division by 24 h).
    idx <- AUC_ENRO / mic

    # Inhibitory sigmoid Emax equation with canonical roles (see ini()).
    # eff_log10 is the signed change in log10 CFU/mL accrued over a 24 h
    # incubation: e0 at zero exposure, emax at saturating exposure.
    eff_log10 <- e0 - (e0 - emax) * idx^hill / (idx^hill + ec50^hill)

    # Li 2017 fitted the 0-to-24 h difference only. Spreading it uniformly
    # across the window makes log10(bact) change by exactly eff_log10 over
    # 24 h:  d(log10 N)/dt = eff_log10 / 24  =>  dN/dt = ln(10) * (eff_log10 / 24) * N
    d/dt(bact) <- log(10) * (eff_log10 / 24) * bact
    bact(0) <- 10^log10_cfu0

    # log10 CFU/mL with a 1-CFU/mL floor so the log stays finite near the
    # maximal-kill plateau.
    log_cfu <- log10(bact + 1)
    log_cfu ~ add(addSd)
  })
}
