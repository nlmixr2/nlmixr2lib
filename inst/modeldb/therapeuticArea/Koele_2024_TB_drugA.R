Koele_2024_TB_drugA <- function() {
  description <- "Methodology reference. Simulation (data-generating) PK-PD model for the HYPOTHETICAL antituberculosis 'drug A' of Koele 2024, a clinical-trial-simulation study of the power to identify exposure-response relationships in 14-day phase IIa (early bactericidal activity) trials. Drug A has the highest activity of the four, selected to give at most about a 3 log10 CFU/mL decrease over 14 days. PK: one-compartment model with linear elimination (CL 10 L/h, V 40 L, 30% CV between-subject variability on CL), all author-chosen constants. PD: the joint bilinear CFU + time-to-positivity model fitted to the HIGHRIF1 trial (Koele_2024_rifampicin_bacterialload), with both slopes of both endpoints scaled by an Emax function of the individual daily AUC0-24 sharing one EC50; baselines, between-subject variability and residual error are the HIGHRIF1 estimates. The default EC50 is the paper's low value (5.0 mg*h/L); the high-EC50 scenario is ini(lec50 = log(25)). A $PRED-style regression on time with the exposure computed from the dose; the central compartment exists only to carry the dose."
  reference <- paste(
    "Koele SE, Dorlo TPC, Upton CM, Aarnoutse RE, Svensson EM. Power to",
    "identify exposure-response relationships in phase IIa pulmonary",
    "tuberculosis trials with multi-dimensional bacterial load modeling.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(3):374-385.",
    "doi:10.1002/psp4.13089.",
    "Values from Table 2 ('Simulation value' column); exposure-response",
    "form from Equation 2 and the example NONMEM code in Supplementary",
    "File S1; PK from Methods 'Pharmacokinetic model'.",
    sep = " "
  )
  vignette <- "Koele_2024_tb_bacterialload"
  # Cc is in mg/L; the outputs are log10 CFU/mL (log_cfu) and log10 h
  # (log_ttp).
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # log_ttp is declared paper-specific, as in Koele_2025_btz043_bacterialload.
  paper_specific_compartments <- c("log_ttp")

  compartmentData <- list(
    central = list(
      analyte = "hypothetical drug A",
      units = "mg",
      specimen = "plasma",
      verified = FALSE
    )
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = 1000L,
    n_studies = 0L,
    disease_state = "Virtual adults with pulmonary tuberculosis in a simulated 14-day phase IIa monotherapy trial; the bacterial-load model was fitted to the HIGHRIF1 trial (see Koele_2024_rifampicin_bacterialload).",
    dose_range = "Once daily for 14 days at 1x, 2x, 3x or 4x a base daily dose (250 virtual patients per level). The base dose is not stated in the paper; Figure 2 is reproduced by 100 mg (typical AUC0-24 10 mg*h/L), see the vignette.",
    notes = "Koele 2024 Methods 'Pharmacokinetic model', 'Exposure-response relationship' and 'Clinical trial simulations'. No real patients: the PK constants and the drug-specific Emax values are author-chosen; the remaining PD values are the HIGHRIF1 estimates. No covariate relationships."
  )

  ini({
    # -----------------------------------------------------------------------
    # PK: author-chosen constants, never estimated (the PD analysis model
    # takes AUC0-24 as a data item), so fixed(). Koele 2024 Methods
    # 'Pharmacokinetic model': 'Typical virtual patients were assumed to
    # have a volume of distribution of 40 L and a clearance of 10 L/h with
    # an IIV coefficient of variation of 30%.' The route is not stated and
    # does not affect AUC0-24; doses go into central.
    # -----------------------------------------------------------------------
    lcl <- fixed(log(10)); label("Clearance (L/h)") # Table 2 'Clearance (L/h) 10.0'
    lvc <- fixed(log(40)); label("Volume of distribution (L)") # Table 2 'Volume of distribution (L) 40.0'
    etalcl ~ fixed(0.08618) # Table 2 'IIV clearance (CV%) 30.0'; log(1 + 0.30^2)

    # -----------------------------------------------------------------------
    # PD, Table 2 'Simulation value' column. Baselines and node are the
    # HIGHRIF1 estimates; the four Emax values are drug A's. They are not
    # fixed(): the authors re-estimated this model (full Emax vs reduced
    # constant-slope) on every simulated trial, and Supplementary File S1
    # has no FIX on any THETA. Both slopes are positive numbers; CFU
    # declines and TTP increases.
    # -----------------------------------------------------------------------
    le0_cfu <- log(5.92); label("Baseline sputum bacterial load (log10 CFU/mL)") # Table 2 'Baseline CFU' simulation value 5.92
    le0_ttp <- log(2.02); label("Baseline time to positivity (log10 h)") # Table 2 'Baseline TTP' simulation value 2.02
    lemax1_cfu <- log(0.0329); label("Maximum first-phase rate of CFU decline (log10 CFU/mL per h)") # Table 2 'Slope 1 CFU' simulation value 'Max Emax slope ... A: 0.0329'
    lemax1_ttp <- log(0.0056); label("Maximum first-phase rate of TTP increase (log10 h per h)") # Table 2 'Slope 1 TTP' simulation value 'Max Emax slope ... A: 0.0056'
    lemax2_cfu <- log(0.00550); label("Maximum second-phase rate of CFU decline (log10 CFU/mL per h)") # Table 2 'Slope 2 CFU' simulation value 'Max Emax slope ... A: 0.00550'
    lemax2_ttp <- log(0.00093); label("Maximum second-phase rate of TTP increase (log10 h per h)") # Table 2 'Slope 2 TTP' simulation value 'Max Emax slope ... A: 0.00093'
    lec50 <- log(5.0); label("AUC0-24 producing half the maximum slope, shared by all four slopes (mg*h/L)") # Table 2 'EC50 (mg*h/L) Low:5.0 High:25.0'; low value. Supplementary File S1: EC50TTP = EC50CFU = THETA(5)
    lnode_time <- log(62.1); label("Node: time after the start of treatment at which both slopes change (h)") # Table 2 'Node (h)' simulation value 'A, B, C: 62.1'

    # Between-subject variability: HIGHRIF1 estimates (Table 2 'Simulation
    # value' column, identical to the 'Estimate' column). Variance =
    # log(1 + CV^2) per the Table 2 footnote; covariance = r*sqrt(v1*v2).
    etale0_cfu + etale0_ttp ~ c(0.02344, -0.006016, 0.003020) # Table 2 'IIV baseline CFU 15.4', 'IIV baseline TTP 5.5', 'Correlation -71.5'
    etalemax1_cfu + etalemax1_ttp ~ c(0.2926, 0.1537, 0.09119) # Table 2 'IIV slope 1 CFU 58.3', 'IIV slope 1 TTP 30.9', 'Correlation 94.1'
    etalemax2_cfu + etalemax2_ttp ~ c(0.4377, 0.1881, 0.2011) # Table 2 'IIV slope 2 CFU 74.1', 'IIV slope 2 TTP 47.2', 'Correlation 63.4'

    # Residual error: Table 2 values read as $SIGMA variances (see
    # Koele_2024_rifampicin_bacterialload and the vignette); SD = sqrt.
    addSd_log_cfu <- 0.5779; label("Additive residual error on CFU (log10 CFU/mL)") # Table 2 'Additive error CFU' simulation value 0.334, read as a variance: sqrt(0.334)
    addSd_log_ttp <- 0.08167; label("Additive residual error on TTP (log10 h)") # Table 2 'Additive error TTP' simulation value 0.00667, read as a variance: sqrt(0.00667)
  })

  model({
    # 1. PK. cl and vc give the concentration profile; the PD is driven by
    #    the individual daily AUC0-24, which for once-daily dosing of a
    #    one-compartment linear model at steady state is dose / CL.
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    kel <- cl / vc
    d/dt(central) <- -kel * central
    Cc <- central / vc

    # podo(central) is the most recent dose; it is NA before the first dose,
    # when the exposure is zero.
    dose_amt <- podo(central)
    if (is.na(dose_amt)) {
      dose_amt <- 0
    }
    auc24 <- dose_amt / cl

    # 2. Individual baselines and exposure-dependent slopes (Equation 2:
    #    Slope = Emax * AUC0-24 / (EC50 + AUC0-24) with log-normal
    #    between-subject variability, THETA * EXP(ETA) in Supplementary
    #    File S1).
    base_cfu <- exp(le0_cfu + etale0_cfu)
    base_ttp <- exp(le0_ttp + etale0_ttp)
    ec50 <- exp(lec50)
    feff <- auc24 / (ec50 + auc24)
    beta1_cfu <- exp(lemax1_cfu + etalemax1_cfu) * feff
    beta1_ttp <- exp(lemax1_ttp + etalemax1_ttp) * feff
    beta2_cfu <- exp(lemax2_cfu + etalemax2_cfu) * feff
    beta2_ttp <- exp(lemax2_ttp + etalemax2_ttp) * feff

    # 3. Bilinear time split about the node (Equation 1), branch-free.
    nodetime <- exp(lnode_time)
    tafter <- (time > nodetime) * (time - nodetime)
    tbefore <- time - tafter

    # 4. Predictions; time is time after the start of treatment (h).
    log_cfu <- base_cfu - beta1_cfu * tbefore - beta2_cfu * tafter
    log_ttp <- base_ttp + beta1_ttp * tbefore + beta2_ttp * tafter

    log_cfu ~ add(addSd_log_cfu)
    log_ttp ~ add(addSd_log_ttp)
  })
}
