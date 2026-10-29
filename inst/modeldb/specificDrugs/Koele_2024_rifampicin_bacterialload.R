Koele_2024_rifampicin_bacterialload <- function() {
  description <- "Joint bilinear model for the decline in sputum mycobacterial load during the first 14 days of treatment in adults with sputum-smear-positive pulmonary tuberculosis in the PanACEA HIGHRIF1 trial (rifampicin 10-50 mg/kg once daily as monotherapy for 7 days, then with standard-dose isoniazid, pyrazinamide and ethambutol for 7 days), fitted simultaneously to log10 colony-forming units and log10 time to positivity. Both endpoints share one estimated breakpoint (node, 62.1 h); before and after the node each endpoint changes linearly in time with its own slope. Baselines and slopes carry log-normal between-subject variability, correlated pairwise between the two endpoints, and the residual error is additive on the log10 scale. The fitted model has no exposure-response term (the slopes are pooled across rifampicin doses); it is the base of the hypothetical-drug simulation models Koele_2024_TB_drugA to Koele_2024_TB_drugD. A $PRED regression on time with no differential equations and no dosing events."
  reference <- paste(
    "Koele SE, Dorlo TPC, Upton CM, Aarnoutse RE, Svensson EM. Power to",
    "identify exposure-response relationships in phase IIa pulmonary",
    "tuberculosis trials with multi-dimensional bacterial load modeling.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(3):374-385.",
    "doi:10.1002/psp4.13089.",
    "Typical values, CV% and correlations from Table 2 ('Estimate' column);",
    "structure from the Methods (Equation 1) and the example NONMEM code in",
    "Supplementary File S1.",
    sep = " "
  )
  vignette <- "Koele_2024_tb_bacterialload"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "log10 CFU/mL (log_cfu) and log10 h (log_ttp); no drug concentration"
  )

  # No d/dt() states: log_cfu and log_ttp are algebraic predictions. log_cfu
  # is the registered canonical log-CFU output; log_ttp is declared
  # paper-specific, as in Koele_2025_btz043_bacterialload.
  paper_specific_compartments <- c("log_ttp")

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = 99L,
    n_studies = 1L,
    disease_state = "Adults with sputum-smear-positive pulmonary tuberculosis enrolled in the PanACEA HIGHRIF1 dose-ranging trial of high-dose rifampicin.",
    dose_range = "Rifampicin 10 to 50 mg/kg once daily in separate dose groups, as monotherapy for 7 days, followed by 7 days on which standard-dose isoniazid, pyrazinamide and ethambutol were co-administered; the dose groups are pooled in this model (no exposure-response term).",
    sampling = "Sputum on Days 0, 1, 2, 3, 4, 5, 6, 7, 9 and 14, each sample cultured in duplicate for colony-forming units on solid medium and for time to positivity in liquid culture.",
    n_observations = "1959 CFU and 1968 TTP measurements; 151 (7.7%) CFU and 22 (1.1%) TTP samples culture-negative.",
    notes = "Koele 2024 Methods 'Pharmacodynamic model' and Results 'Pharmacodynamic model'. One further participant was excluded because only one positive culture sample was available. Age, weight, sex and region are not reported in Koele 2024; they are in the primary HIGHRIF1 publications cited there (refs 6 and 13)."
  )

  ini({
    # -----------------------------------------------------------------------
    # Typical values: Koele 2024 Table 2, 'Estimate (SIR 95% CI)' column.
    # Sign convention follows the paper and the Supplementary File S1 code:
    # both slopes are POSITIVE numbers; CFU DECLINES by slope * t and TTP
    # INCREASES by slope * t (hence Table 2's '-log10h*h-1' unit for TTP).
    # -----------------------------------------------------------------------
    le0_cfu <- log(5.92); label("Baseline sputum bacterial load (log10 CFU/mL)") # Table 2 'Baseline CFU (log10 CFU*mL-1) 5.92 (5.75-6.08)'
    le0_ttp <- log(2.02); label("Baseline time to positivity (log10 h)") # Table 2 'Baseline TTP (log10h) 2.02 (1.99-2.04)'
    lslope1_cfu <- log(0.0141); label("First-phase rate of CFU decline (log10 CFU/mL per h)") # Table 2 'Slope 1 CFU (log10 CFU*mL-1*h-1) 0.0141 (0.0121-0.0162)'
    lslope1_ttp <- log(0.00353); label("First-phase rate of TTP increase (log10 h per h)") # Table 2 'Slope 1 TTP (-log10h*h-1) 0.00353 (0.00320-0.00391)'
    lslope2_cfu <- log(0.00662); label("Second-phase rate of CFU decline (log10 CFU/mL per h)") # Table 2 'Slope 2 CFU (log10 CFU*mL-1*h-1) 0.00662 (0.00552-0.00743)'
    lslope2_ttp <- log(0.000910); label("Second-phase rate of TTP increase (log10 h per h)") # Table 2 'Slope 2 TTP (-log10h*h-1) 0.000910 (0.000810-0.00100)'
    lnode_time <- log(62.1); label("Node: time after the start of treatment at which both slopes change (h)") # Table 2 'Node (h) 62.1 (56.3-67.7)'

    # -----------------------------------------------------------------------
    # Between-subject variability. Table 2 reports CV% with the footnote
    # 'CV% was calculated as sqrt(exp(OMEGA2) - 1)', so each variance is
    # log(1 + CV^2) and each covariance is r * sqrt(var1 * var2). The three
    # BLOCK(2) pairs match the $OMEGA BLOCK(2) x 3 of Supplementary File S1.
    # -----------------------------------------------------------------------
    etale0_cfu + etale0_ttp ~ c(0.02344, -0.006016, 0.003020) # Table 2 'IIV baseline CFU 15.4', 'IIV baseline TTP 5.5', 'Correlation baseline CFU-baseline TTP -71.5'
    etalslope1_cfu + etalslope1_ttp ~ c(0.2926, 0.1537, 0.09119) # Table 2 'IIV slope 1 CFU 58.3', 'IIV slope 1 TTP 30.9', 'Correlation slope 1 CFU-slope 1 TTP 94.1'
    etalslope2_cfu + etalslope2_ttp ~ c(0.4377, 0.1881, 0.2011) # Table 2 'IIV slope 2 CFU 74.1', 'IIV slope 2 TTP 47.2', 'Correlation slope 2 CFU-slope 2 TTP 63.4'

    # -----------------------------------------------------------------------
    # Residual error, additive on the log10 scale. Table 2 prints 0.334
    # (CFU) and 0.00667 (TTP); these are read as $SIGMA VARIANCES and the
    # SDs below are their square roots. Evidence: (i) the Figure 1 pcVPC
    # 95% prediction interval at Day 0 is about 0.53 log10 h wide for TTP
    # and 4.4 log10 CFU/mL wide for CFU, which the variance reading
    # reproduces (0.54 and 4.3) and the SD reading does not (0.44 and 3.9);
    # (ii) an SD of 0.00667 log10 h would be a 1.5% replicate CV for MGIT
    # time to positivity, an order of magnitude below the 0.0624 log10 h
    # this group estimated with the same methods in Koele 2025 (BTZ-043),
    # whose $SIGMA variances (0.30-0.32 CFU, 0.0035-0.0043 TTP) are on the
    # same scale as the two Table 2 numbers. See the vignette.
    #
    # Not encoded: the correlation between the CFU and TTP additive errors
    # (Table 2: -47.4%) and the replicate-error correlations (CFU 90.3%,
    # TTP 82.7%) of the $SIGMA BLOCK(4) for duplicate samples; nlmixr2lib
    # has no idiomatic residual correlation across outputs or replicates.
    # -----------------------------------------------------------------------
    addSd_log_cfu <- 0.5779; label("Additive residual error on CFU (log10 CFU/mL)") # Table 2 'Additive error CFU (log10 CFU*mL-1) 0.334 (0.291-0.360)', read as a variance: sqrt(0.334)
    addSd_log_ttp <- 0.08167; label("Additive residual error on TTP (log10 h)") # Table 2 'Additive error TTP (log10h) 0.00667 (0.00606-0.00730)', read as a variance: sqrt(0.00667)
  })

  model({
    # 1. Individual baselines and slopes, log-normal (Supplementary File S1:
    #    THETA * EXP(ETA); Methods: 'Individual baseline bacterial load and
    #    slope were assumed to be log-normally distributed').
    base_cfu <- exp(le0_cfu + etale0_cfu)
    base_ttp <- exp(le0_ttp + etale0_ttp)
    beta1_cfu <- exp(lslope1_cfu + etalslope1_cfu)
    beta1_ttp <- exp(lslope1_ttp + etalslope1_ttp)
    beta2_cfu <- exp(lslope2_cfu + etalslope2_cfu)
    beta2_ttp <- exp(lslope2_ttp + etalslope2_ttp)

    # 2. Bilinear time split about the node (Equation 1). tbefore is the time
    #    spent in the first phase and tafter the time spent in the second, so
    #    the IF (TIME.GT.NODE) branches of the $PRED code are written
    #    branch-free and stay continuous at the node.
    nodetime <- exp(lnode_time)
    tafter <- (time > nodetime) * (time - nodetime)
    tbefore <- time - tafter

    # 3. Predictions; time is time after the start of treatment (h).
    log_cfu <- base_cfu - beta1_cfu * tbefore - beta2_cfu * tafter
    log_ttp <- base_ttp + beta1_ttp * tbefore + beta2_ttp * tafter

    log_cfu ~ add(addSd_log_cfu)
    log_ttp ~ add(addSd_log_ttp)
  })
}
