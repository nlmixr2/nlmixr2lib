Zhou_2017b_tulathromycin_pkpd_serum <- function() {
  description <- paste0(
    "Veterinary (pig). Ex vivo inhibitory sigmoid Emax PK/PD-integration model for the ",
    "antibacterial effect of tulathromycin against Streptococcus suis ATCC 43765 incubated for ",
    "12 h in SERUM sampled from healthy piglets after a single 2.5 mg/kg intramuscular dose. ",
    "Zhou 2017b parameterises the effect as E = E0 + (Emax - E0) * Ce^N / (EC50^N + Ce^N), where ",
    "E is the SIGNED change in log10 CFU/mL over the 12 h incubation (positive = net growth, ",
    "negative = net kill), E0 is the change in drug-free control serum, Emax the maximal ",
    "antibacterial effect and Ce the PK/PD index AUC(0-12 h)/MIC (h). Parameters from Zhou 2017b ",
    "Table 4: e0 = 1.76, emax = -3.87 log10 CFU/mL, EC50 = 17.9 h, N = 1.18. The index is formed ",
    "in model() as AUC_TULA / mic with the ex vivo serum MIC of 0.031 ug/mL (Figure 2 caption); ",
    "in the ex vivo design each serum sample holds one fixed concentration C, so AUC_TULA = 12 h ",
    "* C. The bacterial density bact (CFU/mL) is integrated as d/dt(bact) = ln(10) * (E / 12) * ",
    "bact so log10(bact) changes by exactly E over each 12 h window. No PK component: the paper's ",
    "PK analysis is non-compartmental only, so exposure enters as the covariate AUC_TULA. Neither ",
    "between-animal variability nor a residual error magnitude was reported, so there are no ",
    "etas and addSd is fixed at 0."
  )
  reference <- paste(
    "Zhou YF, Peng HM, Bu MX, Liu YH, Sun J, Liao XP.",
    "Pharmacodynamic Evaluation and PK/PD-Based Dose Prediction of Tulathromycin:",
    "A Potential New Indication for Streptococcus suis Infection.",
    "Front Pharmacol. 2017;8:684. doi:10.3389/fphar.2017.00684. PMCID: PMC5627010.",
    "Model equation from Materials and Methods 'PK/PD Modeling and Data Analysis';",
    "parameter values from Table 4; ex vivo serum MIC from the Figure 2 caption.",
    sep = " "
  )
  vignette <- "Zhou_2017b_tulathromycin"

  units <- list(
    time = "h",
    dosing = "h*ug/mL (tulathromycin AUC over the 12 h ex vivo incubation, supplied as a covariate)",
    concentration = "log10 CFU/mL (observation)"
  )

  depends <- c("AUC_TULA")
  paper_specific_compartments <- c("bact")

  compartmentData <- list(
    bact = list(analyte = "Streptococcus suis ATCC 43765", units = "CFU/mL", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    AUC_TULA = list(
      description = "Tulathromycin area under the concentration-time curve over the 12 h ex vivo incubation (AUC0-12h); in the ex vivo design this is the static sample concentration times 12 h",
      units = "h*ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Zhou 2017b relates the 12 h change in log10 CFU/mL to the 'ex vivo AUC12h/MIC in ",
        "serum'. Each serum sample holds one fixed concentration C for the 12 h incubation, so ",
        "AUC_TULA = 12 * C (C in ug/mL). The 12 h window is this paper's incubation time: supply ",
        "12 * C, not 24 * C. Figure 1 plots the mean serum concentrations after the 2.5 mg/kg ",
        "intramuscular dose (Cmax 0.82 ug/mL at 0.23 h, Table 1). Set to 0 for drug-free ",
        "control serum so the sigmoid term vanishes and the predicted 12 h change reduces to e0."
      ),
      source_name = "AUC12h (numerator of the PK/PD index AUC12h/MIC, Zhou 2017b 'PK/PD Modeling and Data Analysis' and Table 4)"
    )
  )

  population <- list(
    species = "pig (Duroc x Landrace x Yorkshire piglets, ex vivo serum)",
    n_subjects = 8L,
    n_studies = 1L,
    weight_range = "20.5-23.6 kg",
    sex_female_pct = 50,
    disease_state = "Healthy castrated crossbred piglets; ex vivo challenge of their serum with Streptococcus suis",
    dose_range = "single 2.5 mg/kg tulathromycin (Draxxin) intravenously and intramuscularly in a two-period crossover",
    regions = "China (South China Agricultural University, Guangzhou)",
    organism = "Streptococcus suis ATCC 43765; tulathromycin MIC 1 ug/mL in Mueller-Hinton broth and 0.031 ug/mL in porcine serum (Figure 2 caption; Table 3 serum geometric mean 0.03 ug/mL)",
    notes = paste0(
      "Four male and four female piglets received 2.5 mg/kg intravenously and intramuscularly ",
      "in a two-period crossover. The ex vivo killing curves (Figure 2B) used serum from 0 to ",
      "360 h after the intramuscular dose, inoculated at about 10^6 CFU/mL and counted over 12 h ",
      "(detection limit 40 CFU/mL). The sigmoid Emax model was fitted in WinNonlin 6.1 ",
      "(R2 = 0.9711). Serum protein binding was 36.3% (unbound fraction 0.637)."
    )
  )

  ini({
    # =================================================================
    # Zhou 2017b Table 4 -- sigmoid Emax PK/PD-integration model, serum
    # =================================================================
    # Zhou 2017b 'PK/PD Modeling and Data Analysis' equation, as printed:
    #   E = E0 + (Emax - E0) * Ce^N / (EC50^N + Ce^N)
    # with E0 'the log10 change in visible counts of drug-free control'
    # and Emax 'the maximal antibacterial effect' -- the canonical roles,
    # so no renaming is needed. Table 4 also prints Emax - E0 = -5.63,
    # which these two values reproduce exactly.
    e0 <- 1.76
    label("Change in S. suis count over 12 h in drug-free serum (log10 CFU/mL)") # Zhou 2017b Table 4, 'E0 (log CFU/mL)' = 1.76
    emax <- -3.87
    label("Maximal change in S. suis count over 12 h (log10 CFU/mL; negative = reduction)") # Zhou 2017b Table 4, 'Emax (log CFU/mL)' = -3.87
    lec50 <- log(17.9)
    label("AUC(0-12 h)/MIC producing 50% of the maximal effect EC50 (h)") # Zhou 2017b Table 4, 'EC50 (h)' = 17.9
    lhill <- log(1.18)
    label("Hill coefficient N, steepness of the AUC(0-12 h)/MIC effect curve (unitless)") # Zhou 2017b Table 4, 'Slope (N)' = 1.18

    # =================================================================
    # Minimum inhibitory concentration
    # =================================================================
    # FIXED: a measured property of the challenge strain in the ex vivo
    # matrix. Figure 2B is captioned 'MICserum = 0.031 ug/mL'; Table 3
    # gives the same strain a serum geometric mean of 0.03 ug/mL. The
    # Figure 4 x-axis ends near 315 h, which is 12 * 0.82 / 0.031 = 317 h
    # (the Table 1 intramuscular Cmax over this MIC).
    mic <- fixed(0.031)
    label("Tulathromycin MIC against S. suis ATCC 43765 in porcine serum (ug/mL)") # Zhou 2017b Figure 2 caption, 'MICserum = 0.031 ug/mL'

    # =================================================================
    # Starting bacterial density
    # =================================================================
    # Experimental design input, not an estimated parameter: Methods give
    # 'an initial inoculum of 10^6 CFU/mL' for the ex vivo killing curves.
    # The model predicts a 12 h change, which does not depend on it.
    log10_cfu0 <- fixed(6)
    label("log10 starting S. suis density in the ex vivo incubation (log10 CFU/mL)") # Zhou 2017b Methods, 'Antimicrobial Activities of Tulathromycin and PAEs', initial inoculum 10^6 CFU/mL

    # =================================================================
    # Residual error
    # =================================================================
    # Zhou 2017b reports only the coefficient of determination of the
    # sigmoid fit (R2 = 0.9711) and no residual SD.
    addSd <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (not reported in Zhou 2017b)") # Zhou 2017b reports R2 only, no residual SD
  })

  model({
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # PK/PD index: ex vivo AUC(0-12 h) divided by the MIC (h).
    ce <- AUC_TULA / mic

    # Zhou 2017b sigmoid Emax equation. eff_log10 is the SIGNED change in
    # log10(CFU/mL) over a 12 h window; it equals e0 at zero exposure and
    # approaches emax at saturation.
    eff_log10 <- e0 + (emax - e0) * ce^hill / (ce^hill + ec50^hill)

    # The paper fitted only the 12 h change. Spreading it uniformly across
    # the window makes log10(bact) change by exactly eff_log10 over 12 h:
    #   d(log10 N)/dt = eff_log10 / 12  =>  d(N)/dt = ln(10) * (eff_log10 / 12) * N
    d/dt(bact) <- log(10) * (eff_log10 / 12) * bact
    bact(0) <- 10^log10_cfu0

    # log10 CFU/mL observation with a 1-CFU/mL floor so the log10 stays
    # finite if bact is driven below 1 CFU/mL.
    log_cfu <- log10(bact + 1)
    log_cfu ~ add(addSd)
  })
}
