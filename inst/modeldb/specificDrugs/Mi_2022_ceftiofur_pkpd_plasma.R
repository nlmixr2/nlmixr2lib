Mi_2022_ceftiofur_pkpd_plasma <- function() {
  description <- paste0(
    "Veterinary (pig). Ex vivo inhibitory sigmoid Imax PK/PD-integration model for the ",
    "antibacterial effect of ceftiofur (measured as its active metabolite desfuroylceftiofur, DFC) ",
    "against Pasteurella multocida strain HB13 incubated for 24 h in PLASMA sampled from ",
    "P. multocida-infected pigs after a single 5 mg/kg intramuscular dose. Mi 2022 eq in ",
    "Section 4.6: E = E0 - Imax * INDEX^N / (INDEX^N + INDEX50^N), where E is the SIGNED ",
    "change in log10 CFU/mL over the 24 h incubation (positive = net growth, negative = net ",
    "kill), E0 is that change in drug-free control plasma and INDEX is the PK/PD index ",
    "AUC24h/MIC (h). Parameters from Mi 2022 Table 3, column 'Plasma': Imax = 9.78 log10 ",
    "CFU/mL, E0 = 3.57 log10 CFU/mL, IC50 (INDEX50) = 60.14 h, N = 1.94. The index is formed ",
    "in model() as AUC_CEFTIOFUR / mic with the strain HB13 ex vivo MIC of 0.06 ug/mL (Mi 2022 ",
    "Section 2.1), which reproduces the paper's ex vivo kill at the plasma sample ",
    "concentrations listed in Figure 2A. The bacterial density bact (CFU/mL) is integrated as ",
    "d/dt(bact) = ln(10) * (E / 24) * bact so log10(bact) changes by exactly E over each 24 h ",
    "window. No PK component: exposure enters as the covariate AUC_CEFTIOFUR. Neither ",
    "between-animal variability nor a residual error magnitude was reported (the +/- values ",
    "in Table 3 are SDs of the estimates), so there are no etas and addSd is fixed at 0. ",
    "Siblings: Mi_2022_ceftiofur_pkpd_balf (the bronchoalveolar lavage fluid fit of the same ",
    "experiment) and Mi_2022_ceftiofur_semimech (the paper's semi-mechanistic time-kill PD ",
    "model driven by a one-compartment surrogate of the swine PBPK model)."
  )
  reference <- paste(
    "Mi K, Pu S, Hou Y, Sun L, Zhou K, Ma W, Xu X, Huo M, Liu Z, Xie C, Qu W, Huang L.",
    "Optimization and validation of dosage regimen for ceftiofur against Pasteurella multocida",
    "in swine by physiological based pharmacokinetic-pharmacodynamic model.",
    "Int J Mol Sci. 2022;23(7):3722. doi:10.3390/ijms23073722. PMCID: PMC8998519.",
    "Model equation from Section 4.6; parameter values from Table 3, column 'Plasma';",
    "ex vivo MIC from Section 2.1.",
    sep = " "
  )
  vignette <- "Mi_2022_ceftiofur"

  units <- list(
    time = "h",
    dosing = "h*ug/mL (ceftiofur AUC over a 24 h window, supplied as a covariate)",
    concentration = "log10 CFU/mL (observation)"
  )

  depends <- c("AUC_CEFTIOFUR")
  paper_specific_compartments <- c("bact")

  compartmentData <- list(
    bact = list(analyte = "Pasteurella multocida strain HB13", units = "CFU/mL", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    AUC_CEFTIOFUR = list(
      description = "Ceftiofur (as desfuroylceftiofur) area under the concentration-time curve over the 24 h exposure window (AUC24h); in the ex vivo design this is the static sample concentration times 24 h",
      units = "h*ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Mi 2022 Section 4.6 defines the sigmoid Imax driver INDEX as AUC24h/MIC (h). In the ex ",
        "vivo design each plasma sample holds a fixed concentration C for the 24 h incubation, ",
        "so AUC24h = 24 * C (C in ug/mL). Figure 2A lists the sample concentrations (for example ",
        "0.19 mg/L at 72 h and 0.05 mg/L at 96 h). Set to 0 for drug-free control plasma so the ",
        "sigmoid term vanishes and the predicted 24 h change reduces to E0."
      ),
      source_name = "AUC24h (numerator of the PK/PD index AUC24h/MIC, Mi 2022 Section 4.6 and Table 3)"
    )
  )

  population <- list(
    species = "pig (crossbred, ex vivo plasma)",
    n_subjects = 6L,
    n_studies = 1L,
    age_range = NA_character_,
    weight_range = "20 +/- 2 kg average body weight",
    sex_female_pct = NA_real_,
    disease_state = "Pasteurella multocida (strain HB13) intranasal infection; fever, cough and dyspnea",
    dose_range = "single intramuscular dose of ceftiofur hydrochloride 5 mg/kg",
    regions = "China (Huazhong Agricultural University, Wuhan)",
    organism = "Pasteurella multocida strain HB13, a wild isolate from Hubei (2017); MIC 0.06 ug/mL and MBC 0.125 ug/mL both in vitro (MHB) and ex vivo",
    notes = paste0(
      "Ex vivo time-kill: about 5 x 10^6 CFU/mL of HB13 incubated in plasma sampled at ",
      "0.33-96 h after the dose, viable counts at 0, 3, 6, 9, 12 and 24 h (limit of detection ",
      "10 CFU/mL). Table 3 parameters are reported as mean +/- SD. The bacteriostatic, ",
      "bactericidal and eradication targets (E = 0, -3, -4) are AUC24h/MIC = 44.02, 89.40 and ",
      "119.90 h; with MIC90 = 0.06 ug/mL, Monte Carlo simulation gave doses of 0.22, 0.46 and ",
      "0.64 mg/kg for 90% target attainment."
    )
  )

  ini({
    # Mi 2022 Section 4.6:
    #   E = E0 - Imax * INDEX^N / (INDEX^N + INDEX50^N)
    # E is the SIGNED 24 h change in log10 CFU/mL; it equals E0 at zero
    # exposure and falls to E0 - Imax at saturating exposure.
    le0 <- log(3.57)
    label("Change in bacterial count in drug-free control plasma over 24 h E0 (log10 CFU/mL)") # Mi 2022 Table 3, Plasma, E0 = 3.57 +/- 0.43

    limax <- log(9.78)
    label("Maximum antibacterial effect Imax, the span from control growth to maximum kill (log10 CFU/mL)") # Mi 2022 Table 3, Plasma, Imax = 9.78 +/- 0.37

    lic50 <- log(60.14)
    label("AUC24h/MIC producing half the maximum effect IC50 (h)") # Mi 2022 Table 3, Plasma, IC50 = 60.14 +/- 3.48

    lhill <- log(1.94)
    label("Hill coefficient N, steepness of the AUC24h/MIC effect curve (unitless)") # Mi 2022 Table 3, Plasma, N = 1.94 +/- 0.55

    # Measured property of the challenge strain, not an estimate.
    mic <- fixed(0.06)
    label("Ceftiofur MIC against P. multocida HB13, ex vivo (ug/mL)") # Mi 2022 Section 2.1: in vitro and ex vivo MIC both 0.06 ug/mL

    # Experimental design input, not an estimate.
    log10_cfu0 <- fixed(log10(5e6))
    label("log10 starting bacterial density in the ex vivo incubation (log10 CFU/mL)") # Mi 2022 Section 4.3.2: bacteria (about 5 x 10^6 CFU/mL) cultured in plasma

    # No residual SD is reported for the PK/PD-integration fit.
    addSd <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (0; not reported in Mi 2022)") # Mi 2022 Table 3 reports mean +/- SD of the estimates only
  })

  model({
    e0 <- exp(le0)
    imax <- exp(limax)
    ic50 <- exp(lic50)
    hill <- exp(lhill)

    # PK/PD index AUC24h/MIC (h)
    aucmic <- AUC_CEFTIOFUR / mic

    # Mi 2022 Section 4.6 sigmoid Imax equation; signed 24 h change in
    # log10 CFU/mL.
    effect <- e0 - imax * aucmic^hill / (aucmic^hill + ic50^hill)

    # Spread the 24 h change uniformly so log10(bact) moves by exactly
    # `effect` over each 24 h window:
    #   d(log10 N)/dt = effect / 24  =>  dN/dt = ln(10) * (effect / 24) * N
    d/dt(bact) <- log(10) * (effect / 24) * bact
    bact(0) <- 10^log10_cfu0

    # log10 CFU/mL with a 1 CFU/mL floor so the log stays finite.
    log_cfu <- log10(bact + 1)
    log_cfu ~ add(addSd)
  })
}
