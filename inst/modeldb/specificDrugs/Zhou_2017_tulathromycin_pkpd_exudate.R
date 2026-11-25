Zhou_2017_tulathromycin_pkpd_exudate <- function() {
  description <- paste0(
    "Veterinary (pig). Ex vivo inhibitory sigmoid Emax PK/PD-integration model for the ",
    "antibacterial effect of tulathromycin against Pasteurella multocida strain CVCC430 incubated ",
    "for 24 h in EXUDATE (inflammatory fluid from a carrageenan-injected subcutaneous tissue cage) sampled from ",
    "healthy piglets after a single 2.5 mg/kg intramuscular dose. Zhou 2017 parameterises the ",
    "effect as E = Emax - (Emax - E0) * Ce^N / (EC50^N + Ce^N), where E is the SIGNED change in ",
    "log10 CFU/mL over the 24 h incubation (positive = net growth, negative = net kill) and Ce is ",
    "the PK/PD index AUC(0-24 h)/MIC (h). NOTE THE INVERTED NAMING: Zhou 2017 defines its Emax as ",
    "the growth in drug-free control samples and its E0 as the maximum growth inhibition, so its ",
    "Emax = 1.85 maps onto the canonical `e0` (effect at zero exposure) and its E0 = -3.16 onto ",
    "the canonical `emax` (maximal drug effect); the equation is algebraically unchanged. ",
    "Parameters from Zhou 2017 Table 6, column 'Exudate': e0 = 1.85, emax = -3.16 log10 ",
    "CFU/mL, EC50 = 33.55 h, N = 15.55. Solving the model for E = 0 and -3 returns 32.4 and 41.8 h ",
    "against the published bacteriostatic and bactericidal targets of 32.42 and 41.85 h; the ",
    "maximal effect of -3.16 log10 never reaches a 4-log10 kill, which is why Table 6 leaves the ",
    "exudate eradication target blank. The index is formed in model() as AUC_TULA / mic with the ",
    "CVCC430 exudate MIC of 0.03 ug/mL (Zhou 2017 Results, 'MIC and MBC Values'); in the ex vivo ",
    "design each sample holds one fixed concentration C, so AUC_TULA = 24 h * C. The bacterial ",
    "density bact (CFU/mL) is integrated as d/dt(bact) = ln(10) * (E / 24) * bact so log10(bact) ",
    "changes by exactly E over each 24 h window. No PK component: the paper's PK analysis is ",
    "non-compartmental only, so exposure enters as the covariate AUC_TULA. Neither between-animal ",
    "variability nor a residual error magnitude was reported, so there are no etas and addSd is ",
    "fixed at 0. Siblings: Zhou_2017_tulathromycin_pkpd_serum and Zhou_2017_tulathromycin_pkpd_transudate."
  )
  reference <- paste(
    "Zhou Q, Zhang G, Wang Q, Liu W, Huang Y, Yu P, Li Y, Ding H, Fang B.",
    "Pharmacokinetic/Pharmacodynamic Modeling of Tulathromycin against Pasteurella multocida",
    "in a Porcine Tissue Cage Model.",
    "Front Pharmacol. 2017;8:392. doi:10.3389/fphar.2017.00392. PMCID: PMC5487385.",
    "Model equation from Materials and Methods 'PK/PD Integration and Modeling';",
    "parameter values from Table 6, column 'Exudate'; MIC from Results 'MIC and MBC Values'.",
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
    bact = list(
      analyte = "Pasteurella multocida strain CVCC430",
      units = "CFU/mL",
      specimen = "tissue",
      verified = TRUE
    ) # tissue-cage fluid
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
        "24 h', so AUC_TULA = 24 * C (C in ug/mL). Figure 2 plots the mean exudate ",
        "concentrations after the 2.5 mg/kg intramuscular dose (Cmax 0.11 ug/mL at 8 h, ",
        "Table 3). Set to 0 for drug-free control fluid so the sigmoid term vanishes and the ",
        "predicted 24 h change reduces to e0."
      ),
      source_name = "AUC(0-24 h) (numerator of the PK/PD index AUC(0-24 h)/MIC, Zhou 2017 'PK/PD Integration and Modeling' and Table 6)"
    )
  )

  population <- list(
    species = "pig (Duroc x Landrace x Yorkshire piglets, ex vivo tissue-cage exudate)",
    n_subjects = 6L,
    n_studies = 1L,
    weight_range = "13-15 kg",
    sex_female_pct = 50,
    disease_state = "Healthy piglets with two subcutaneous neck tissue cages; exudate is the inflammatory fluid of the cage injected with 0.5 mL of 1% carrageenan; ex vivo challenge with Pasteurella multocida",
    dose_range = "single intramuscular dose of tulathromycin 2.5 mg/kg (Draxxin 100 mg/mL)",
    regions = "China (South China Agricultural University, Guangzhou)",
    organism = "Pasteurella multocida strain CVCC430 (China Veterinary Culture Collection Centre); tulathromycin MIC 0.03 ug/mL and MBC 0.06 ug/mL in exudate",
    notes = paste0(
      "Silicone tissue cages (6.5 cm, 24 holes, 2.3 cm2 exchange surface) were implanted ",
      "subcutaneously on each side of the neck; one received 0.5 mL of 1% carrageenan. After 4-5 ",
      "weeks the carrageenan-injected cage held exudate. Exudate was sampled at 1-360 h after the ",
      "intramuscular dose and each sample was inoculated with CVCC430 and counted at 0, 2, 5, 8 ",
      "and 24 h of incubation (Figure 5C). The sigmoid Emax model was fitted in WinNonlin."
    )
  )

  ini({
    # =================================================================
    # Zhou 2017 Table 6 -- inhibitory sigmoid Emax PK/PD-integration model
    # column: Exudate
    # =================================================================
    # Zhou 2017 'PK/PD Integration and Modeling' equation, as printed:
    #   E = Emax - (Emax - E0) * Ce^N / (EC50^N + Ce^N)
    # with Zhou 2017's 'Emax' = drug-free control growth and 'E0' =
    # maximum growth inhibition (the reverse of the usual convention).
    # This file uses the canonical roles:
    #   e0   <- Zhou 2017's Emax = +1.85  (effect at zero exposure)
    #   emax <- Zhou 2017's E0   = -3.16  (maximal drug effect)
    # and model() writes the algebraically identical
    #   E = e0 - (e0 - emax) * Ce^N / (Ce^N + EC50^N).
    # Solving E = 0, -3 gives 32.4, 41.8 h versus the Table 6 targets
    # 32.42, 41.85 h. N = 15.55 makes the curve almost a step at EC50.
    e0 <- 1.85
    label("Change in P. multocida count over 24 h in drug-free exudate (log10 CFU/mL)") # Zhou 2017 Table 6, Exudate, 'Emax(log10 CFU/mL)' = 1.85 (defined in Methods as the drug-free control growth)
    emax <- -3.16
    label("Maximal change in P. multocida count over 24 h (log10 CFU/mL; negative = reduction)") # Zhou 2017 Table 6, Exudate, 'E0 (log10 CFU/mL)' = -3.16 (defined in Methods as the maximum growth inhibition)
    lec50 <- log(33.55)
    label("AUC(0-24 h)/MIC producing 50% of the maximal effect EC50 (h)") # Zhou 2017 Table 6, Exudate, 'EC50 (h)' = 33.55
    lhill <- log(15.55)
    label("Hill coefficient N, steepness of the AUC(0-24 h)/MIC effect curve (unitless)") # Zhou 2017 Table 6, Exudate, 'Slope (N)' = 15.55

    # =================================================================
    # Minimum inhibitory concentration
    # =================================================================
    # FIXED: a measured property of the challenge strain in this fluid.
    mic <- fixed(0.03)
    label("Tulathromycin MIC against P. multocida CVCC430 in exudate (ug/mL)") # Zhou 2017 Results 'MIC and MBC Values': CVCC430 MIC in transudate and exudate 'were all 0.03 ug/mL'; Table 5 footnote 'MIC (transudate, exudate) = 0.03'

    # =================================================================
    # Starting bacterial density
    # =================================================================
    # Experimental design input, not an estimated parameter. The text
    # gives no ex vivo inoculum; every curve in Figure 5 starts at about
    # 7.2 log10 CFU/mL.
    log10_cfu0 <- fixed(7.2)
    label("log10 starting P. multocida density in the ex vivo incubation (log10 CFU/mL)") # Zhou 2017 Figure 5, time-zero counts digitised by the maintainers (about 7.2 in all three panels)

    # =================================================================
    # Residual error
    # =================================================================
    # Zhou 2017 reports no residual SD for the sigmoid fits.
    addSd <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (not reported in Zhou 2017)") # Zhou 2017 reports no residual SD
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
