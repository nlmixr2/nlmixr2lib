Zhou_2017_tulathromycin_pkpd_serum <- function() {
  description <- paste0(
    "Veterinary (pig). Ex vivo inhibitory sigmoid Emax PK/PD-integration model for the ",
    "antibacterial effect of tulathromycin against Pasteurella multocida strain CVCC430 incubated ",
    "for 24 h in SERUM sampled from healthy piglets after a single 2.5 mg/kg intramuscular dose. ",
    "Zhou 2017 parameterises the effect as E = Emax - (Emax - E0) * Ce^N / (EC50^N + Ce^N), where ",
    "E is the SIGNED change in log10 CFU/mL over the 24 h incubation (positive = net growth, ",
    "negative = net kill) and Ce is the PK/PD index AUC(0-24 h)/MIC (h). NOTE THE INVERTED NAMING: ",
    "Zhou 2017 defines its Emax as the growth in drug-free control samples and its E0 as the ",
    "maximum growth inhibition, so its Emax = 1.60 maps onto the canonical `e0` (effect at zero ",
    "exposure) and its E0 = -4.91 onto the canonical `emax` (maximal drug effect); the equation ",
    "is algebraically unchanged. Parameters from Zhou 2017 Table 6, column 'Serum': e0 = 1.60, ",
    "emax = -4.91 log10 CFU/mL, EC50 = 58.80 h, N = 4.03. Solving the model for E = 0, -3 and -4 ",
    "returns 44.5, 73.1 and 92.3 h against the published bacteriostatic, bactericidal and ",
    "eradication targets of 44.55, 73.19 and 92.44 h. The index is formed in model() as ",
    "AUC_TULA / mic with the CVCC430 serum MIC of 0.03 ug/mL (Zhou 2017 Table 4); in the ex vivo ",
    "design each serum sample holds one fixed concentration C, so AUC_TULA = 24 h * C. The ",
    "bacterial density bact (CFU/mL) is integrated as d/dt(bact) = ln(10) * (E / 24) * bact so ",
    "log10(bact) changes by exactly E over each 24 h window. No PK component: the paper's PK ",
    "analysis is non-compartmental only, so exposure enters as the covariate AUC_TULA. Neither ",
    "between-animal variability nor a residual error magnitude was reported, so there are no ",
    "etas and addSd is fixed at 0. Siblings: Zhou_2017_tulathromycin_pkpd_transudate and ",
    "Zhou_2017_tulathromycin_pkpd_exudate (the tissue-cage fluid fits of the same experiment)."
  )
  reference <- paste(
    "Zhou Q, Zhang G, Wang Q, Liu W, Huang Y, Yu P, Li Y, Ding H, Fang B.",
    "Pharmacokinetic/Pharmacodynamic Modeling of Tulathromycin against Pasteurella multocida",
    "in a Porcine Tissue Cage Model.",
    "Front Pharmacol. 2017;8:392. doi:10.3389/fphar.2017.00392. PMCID: PMC5487385.",
    "Model equation from Materials and Methods 'PK/PD Integration and Modeling';",
    "parameter values from Table 6, column 'Serum'; MIC from Table 4 (strain CVCC430, serum).",
    sep = " "
  )
  vignette <- "Zhou_2017_tulathromycin"

  units <- list(
    time = "h",
    dosing = "h*ug/mL (tulathromycin AUC over a 24 h window, supplied as a covariate)",
    concentration = "log10 CFU/mL (observation)"
  )

  depends <- c("AUC_TULA")
  paper_specific_compartments <- c("bact")

  compartmentData <- list(
    bact = list(analyte = "Pasteurella multocida strain CVCC430", units = "CFU/mL", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    AUC_TULA = list(
      description = "Tulathromycin area under the concentration-time curve over the 24 h exposure window (AUC0-24h); in the ex vivo design this is the static sample concentration times 24 h",
      units = "h*ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Zhou 2017 'PK/PD Integration and Modeling' estimates the ex vivo AUC(0-24 h) 'by ",
        "multiplying the average drug concentration at each time by the incubation period of ",
        "24 h', so AUC_TULA = 24 * C (C in ug/mL). Figure 2 plots the mean serum concentrations ",
        "after the 2.5 mg/kg intramuscular dose (Cmax 0.74 ug/mL at 0.25 h, Table 1). Set to 0 ",
        "for drug-free control serum so the sigmoid term vanishes and the predicted 24 h change ",
        "reduces to e0."
      ),
      source_name = "AUC(0-24 h) (numerator of the PK/PD index AUC(0-24 h)/MIC, Zhou 2017 'PK/PD Integration and Modeling' and Table 6)"
    )
  )

  population <- list(
    species = "pig (Duroc x Landrace x Yorkshire piglets, ex vivo serum)",
    n_subjects = 6L,
    n_studies = 1L,
    weight_range = "13-15 kg",
    sex_female_pct = 50,
    disease_state = "Healthy piglets with two subcutaneous neck tissue cages (one carrageenan-inflamed); ex vivo challenge with Pasteurella multocida",
    dose_range = "single intramuscular dose of tulathromycin 2.5 mg/kg (Draxxin 100 mg/mL)",
    regions = "China (South China Agricultural University, Guangzhou)",
    organism = "Pasteurella multocida strain CVCC430 (China Veterinary Culture Collection Centre); tulathromycin MIC 0.25 ug/mL in Mueller-Hinton broth and 0.03 ug/mL in serum, MBC 0.06 ug/mL in serum (Table 4)",
    notes = paste0(
      "Twelve piglets (six male, six female) were split into an intravenous group and an ",
      "intramuscular group of six; the ex vivo time-kill curves (Figure 5) used the ",
      "intramuscular-group samples. Each serum sample was inoculated with CVCC430 and counted ",
      "at 0, 2, 5, 8 and 24 h of incubation (limit of detection 200 CFU/mL). The sigmoid Emax ",
      "model was fitted in WinNonlin. Zhou 2017 also reports Cmax/MIC and T>MIC fits for serum ",
      "(Figure 4); the Cmax/MIC fit is the same fit on a rescaled axis (ex vivo AUC(0-24 h) is ",
      "24 h times the sample concentration) and T>MIC fitted poorly (R2 = 0.82), so only the ",
      "AUC(0-24 h)/MIC parameters of Table 6 are packaged."
    )
  )

  ini({
    # =================================================================
    # Zhou 2017 Table 6 -- inhibitory sigmoid Emax PK/PD-integration model
    # column: Serum
    # =================================================================
    # Zhou 2017 'PK/PD Integration and Modeling' equation, as printed:
    #   E = Emax - (Emax - E0) * Ce^N / (EC50^N + Ce^N)
    # Zhou 2017's own definitions: 'Emax is the corresponding bacterial
    # growth in the absence of drug (control samples); E0 is the maximum
    # antibacterial growth inhibition'. That is the reverse of the usual
    # convention, and the Table 6 signs confirm it (Emax > 0, E0 < 0).
    # This file uses the canonical roles:
    #   e0   <- Zhou 2017's Emax = +1.60  (effect at zero exposure)
    #   emax <- Zhou 2017's E0   = -4.91  (maximal drug effect)
    # and model() writes the algebraically identical
    #   E = e0 - (e0 - emax) * Ce^N / (Ce^N + EC50^N).
    # Solving E = 0, -3, -4 gives 44.5, 73.1, 92.3 h versus the Table 6
    # targets 44.55, 73.19, 92.44 h.
    e0 <- 1.60
    label("Change in P. multocida count over 24 h in drug-free serum (log10 CFU/mL)") # Zhou 2017 Table 6, Serum, 'Emax(log10 CFU/mL)' = 1.60 (defined in Methods as the drug-free control growth)
    emax <- -4.91
    label("Maximal change in P. multocida count over 24 h (log10 CFU/mL; negative = reduction)") # Zhou 2017 Table 6, Serum, 'E0 (log10 CFU/mL)' = -4.91 (defined in Methods as the maximum growth inhibition)
    lec50 <- log(58.80)
    label("AUC(0-24 h)/MIC producing 50% of the maximal effect EC50 (h)") # Zhou 2017 Table 6, Serum, 'EC50 (h)' = 58.80
    lhill <- log(4.03)
    label("Hill coefficient N, steepness of the AUC(0-24 h)/MIC effect curve (unitless)") # Zhou 2017 Table 6, Serum, 'Slope (N)' = 4.03

    # =================================================================
    # Minimum inhibitory concentration
    # =================================================================
    # FIXED: a measured property of the challenge strain. The CVCC430
    # serum MIC (0.03 ug/mL) is used, not the 7-strain serum average of
    # 0.04 ug/mL that Table 5 divides by: Methods states the ex vivo
    # index uses the MIC 'determined in biological fluid ... for
    # P. multocida CVCC430', and Figure 4's largest AUC(0-24 h)/MIC
    # (about 595 h) is 24 * 0.74 / 0.03, the Table 1 Cmax over 0.03.
    mic <- fixed(0.03)
    label("Tulathromycin MIC against P. multocida CVCC430 in serum (ug/mL)") # Zhou 2017 Table 4, CVCC430, Serum MIC = 0.03 ug/mL

    # =================================================================
    # Starting bacterial density
    # =================================================================
    # Experimental design input, not an estimated parameter. The text
    # gives no ex vivo inoculum; every curve in Figure 5 starts at about
    # 7.2 log10 CFU/mL, and the serum floor of the fit agrees: the
    # detection limit log10(200) = 2.30 plus 4.91 is 7.21.
    log10_cfu0 <- fixed(7.2)
    label("log10 starting P. multocida density in the ex vivo incubation (log10 CFU/mL)") # Zhou 2017 Figure 5, time-zero counts digitised by the maintainers (about 7.2 in all three panels)

    # =================================================================
    # Residual error
    # =================================================================
    # Zhou 2017 reports only the coefficient of determination of the
    # sigmoid fit (R2 = 0.9969, Figure 4) and no residual SD.
    addSd <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (not reported in Zhou 2017)") # Zhou 2017 reports R2 only, no residual SD
  })

  model({
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # PK/PD index: ex vivo AUC(0-24 h) divided by the MIC (h).
    ce <- AUC_TULA / mic

    # Zhou 2017 inhibitory sigmoid Emax equation, written with canonical
    # parameter roles (see the ini() note on the inverted naming).
    # eff_log10 is the SIGNED change in log10(CFU/mL) over a 24 h window;
    # it equals e0 at zero exposure and approaches emax at saturation.
    eff_log10 <- e0 - (e0 - emax) * ce^hill / (ce^hill + ec50^hill)

    # The paper fitted only the 24 h change. Spreading it uniformly across
    # the window makes log10(bact) change by exactly eff_log10 over 24 h:
    #   d(log10 N)/dt = eff_log10 / 24  =>  d(N)/dt = ln(10) * (eff_log10 / 24) * N
    d/dt(bact) <- log(10) * (eff_log10 / 24) * bact
    bact(0) <- 10^log10_cfu0

    # log10 CFU/mL observation with a 1-CFU/mL floor so the log10 stays
    # finite if bact is driven below 1 CFU/mL.
    log_cfu <- log10(bact + 1)
    log_cfu ~ add(addSd)
  })
}
