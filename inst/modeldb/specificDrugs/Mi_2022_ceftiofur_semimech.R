Mi_2022_ceftiofur_semimech <- function() {
  description <- paste0(
    "Veterinary (pig) PK driving an in vitro-estimated PD. Semi-mechanistic PK/PD model of ",
    "ceftiofur against Pasteurella multocida strain HB13 as assembled by Mi 2022 for its ",
    "'PBPK/PD' dose validation (Section 4.8.2, last paragraph; Figure 6). PK: a one-compartment ",
    "model with first-order absorption that Mi 2022 fitted in WinNonlin to the unbound plasma ",
    "desfuroylceftiofur (DFC) concentrations predicted by the Lin 2016 swine PBPK model for ",
    "intramuscular doses of 0.22, 0.46 and 0.64 mg/kg in infected pigs (supplement Table S2: ",
    "V/F = 1114 mL/kg, ka = 0.3381 1/h, k = 0.0393 1/h). PD: the Nielsen 2007 two-state ",
    "bacterial model with a growing, drug-susceptible state (bact_s) and a resting, ",
    "drug-insusceptible state (bact_r), estimated in Monolix from in vitro static time-kill ",
    "curves in Mueller-Hinton broth (Mi 2022 Table 5): kgrowth = 0.2 1/h, kdeath = 0.179 1/h ",
    "(fixed), Bmax = 10^8.48 CFU/mL, Emax = 0.11 1/h, EC50 = 0.14 mg/L, gamma = 8.54. The ",
    "drug adds a sigmoid Emax kill rate to the growing state only; transfer from growing to ",
    "resting is kGR = (kgrowth - kdeath) * (G + R) / Bmax. The PD is driven by the unbound ",
    "plasma concentration Cc. Doses are per kg body weight. No between-animal or ",
    "between-experiment variability and no residual error magnitude were reported, so there ",
    "are no etas and addSd is fixed at 0. IMPORTANT: the Table 5 typical values do not ",
    "reproduce the paper's own Figure 1/5 growth controls or its Figure 6 dose predictions ",
    "(see the vignette Assumptions and deviations); the values are shipped as published, not ",
    "tuned. Siblings: Mi_2022_ceftiofur_pkpd_plasma and Mi_2022_ceftiofur_pkpd_balf (the ex ",
    "vivo sigmoid Imax PK/PD-index models of the same paper)."
  )
  reference <- paste(
    "Mi K, Pu S, Hou Y, Sun L, Zhou K, Ma W, Xu X, Huo M, Liu Z, Xie C, Qu W, Huang L.",
    "Optimization and validation of dosage regimen for ceftiofur against Pasteurella multocida",
    "in swine by physiological based pharmacokinetic-pharmacodynamic model.",
    "Int J Mol Sci. 2022;23(7):3722. doi:10.3390/ijms23073722. PMCID: PMC8998519.",
    "PD equations from Section 4.8.2 (and the supplementary Equation section);",
    "PD parameters from Table 5; PK parameters from supplementary Table S2.",
    "The PK layer is a compartmental surrogate of the PBPK model of",
    "Lin Z, Vahl CI, Riviere JE. Human food safety implications of variation in food animal",
    "drug metabolism. Sci Rep. 2016;6:27907. doi:10.1038/srep27907.",
    "The PD structure is that of Nielsen EI, Viberg A, Lowdin E, Cars O, Karlsson MO,",
    "Sandstrom M. Antimicrob Agents Chemother. 2007;51(1):128-136. doi:10.1128/AAC.00604-06.",
    sep = " "
  )
  vignette <- "Mi_2022_ceftiofur"

  units <- list(
    time = "h",
    dosing = "mg/kg",
    concentration = "mg/L (unbound DFC, Cc); log10 CFU/mL (log_cfu observation)"
  )

  # Growing (G) and resting (R) bacterial states of Mi 2022 Figure 7,
  # following the naming of the Nielsen-lineage model
  # Khan_2015_ciprofloxacin.R (bact_s = growing drug-susceptible,
  # bact_r = resting drug-insusceptible).
  paper_specific_compartments <- c("bact_s", "bact_r")

  compartmentData <- list(
    depot = list(analyte = "ceftiofur", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "desfuroylceftiofur (unbound)", units = "mg/kg", specimen = "plasma", verified = TRUE),
    bact_s = list(
      analyte = "Pasteurella multocida HB13, growing drug-susceptible subpopulation (G)",
      units = "CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    bact_r = list(
      analyte = "Pasteurella multocida HB13, resting drug-insusceptible subpopulation (R)",
      units = "CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list()

  population <- list(
    species = "pig (PK, simulated infected swine) + in vitro (PD, P. multocida HB13 in MHB)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    age_range = NA_character_,
    weight_range = "doses are per kg body weight",
    sex_female_pct = NA_real_,
    disease_state = "Pasteurella multocida respiratory infection (PK simulated for infected pigs by lowering the PBPK hepatic metabolic rate)",
    dose_range = "0.22, 0.46 and 0.64 mg/kg single intramuscular doses (bacteriostatic, bactericidal and eradication doses)",
    regions = "China (Huazhong Agricultural University, Wuhan)",
    organism = "Pasteurella multocida strain HB13; MIC 0.06 ug/mL, MBC 0.125 ug/mL, MPC 0.3 ug/mL",
    notes = paste0(
      "PD data: in vitro static time-kill of about 10^6 CFU/mL HB13 in MHB with ceftiofur at ",
      "1/4 to 8 x MIC plus a drug-free control, counts at 0, 3, 6, 9, 12 and 24 h in triplicate ",
      "(limit of detection 10 CFU/mL); fitted by SAEM in Monolix 2018R1. PK: the one-compartment ",
      "surrogate was fitted separately at each dose to the PBPK-predicted unbound plasma ",
      "concentrations (Figure S1); ka and k were identical across doses and V/F differed only ",
      "in the fifth significant figure. The coupled PK/PD was simulated over 72 h in Mlxplore ",
      "(Figure 6)."
    )
  )

  ini({
    # ---- PK: one-compartment surrogate of the Lin 2016 PBPK -------------
    lka <- log(0.3381)
    label("Absorption rate constant after intramuscular dosing ka (1/h)") # Mi 2022 supplement Table S2: ka = 0.3381 1/h at all three doses
    lvc <- log(1.11421)
    label("Apparent volume of the unbound DFC compartment V/F (L/kg)") # Mi 2022 supplement Table S2: V/F = 1114.21 mL/kg at 0.64 mg/kg (1114.11 and 1114.27 at 0.22 and 0.46 mg/kg)
    lkel <- log(0.0393)
    label("Elimination rate constant k (1/h)") # Mi 2022 supplement Table S2: k = 0.0393 1/h at all three doses

    # ---- PD: bacterial system (Table 5) -------------------------------------
    lkgrow <- log(0.2)
    label("Growth rate constant of the growing subpopulation kgrowth (1/h)") # Mi 2022 Table 5: kgrowth = 0.2 1/h (SE 0.28)
    lkdeath <- fixed(log(0.179))
    label("Natural death rate constant of both subpopulations kdeath (1/h)") # Mi 2022 Table 5: kdeath = 0.179 1/h (fixed)
    # Table 5 prints the unit of Bmax as 1/h (copied from the rows above);
    # Bmax is the maximum bacterial concentration and is reported on the
    # log10 CFU/mL scale, as in Nielsen 2007. See the vignette Errata.
    lbmax <- log(10^8.48)
    label("Maximum bacterial concentration in the system Bmax (CFU/mL)") # Mi 2022 Table 5: Bmax = 8.48 (SE 0.28), read as log10 CFU/mL
    linoc <- fixed(log(1e6))
    label("Initial bacterial density, all in the growing state (CFU/mL)") # Mi 2022 Section 4.8.2: bacteria (10^6 CFU/mL) cultured with ceftiofur

    # ---- PD: drug effect (Table 5) ------------------------------------------
    lemax <- log(0.11)
    label("Maximum ceftiofur-induced kill rate Emax (1/h)") # Mi 2022 Table 5: Emax = 0.11 1/h (SE 0.027)
    lec50 <- log(0.14)
    label("Ceftiofur concentration producing half the maximum kill rate EC50 (mg/L)") # Mi 2022 Table 5: EC50 = 0.14 mg/L (SE 0.031)
    lhill <- log(8.54)
    label("Sigmoidicity coefficient gamma (unitless)") # Mi 2022 Table 5: gamma = 8.54 (no SE printed)

    # ---- Residual error -----------------------------------------------------
    addSd <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (0; not reported in Mi 2022)") # Mi 2022 reports no residual error estimate for the Monolix fit
  })

  model({
    ka <- exp(lka)
    vc <- exp(lvc)
    kel <- exp(lkel)
    kgrow <- exp(lkgrow)
    kdeath <- exp(lkdeath)
    bmax <- exp(lbmax)
    inoc <- exp(linoc)
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # ---- PK (supplement Table S2) -------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    # Unbound DFC plasma concentration (mg/L): dose mg/kg over V/F L/kg
    Cc <- central / vc

    # ---- PD (Section 4.8.2 equations) ---------------------------------------
    # EFFECT = Emax * C^gamma / (EC50^gamma + C^gamma)
    effect <- emax * Cc^hill / (ec50^hill + Cc^hill)

    # kGR = (kgrowth - kdeath) * (G + R) / Bmax -- main-text form. The
    # supplement prints (1 - (G + R) / Bmax), which has no carrying
    # capacity; see the vignette Errata.
    btot <- bact_s + bact_r
    kgr <- (kgrow - kdeath) * btot / bmax

    # dG/dt = kgrowth*G - EFFECT*G - kdeath*G - kGR*G
    d/dt(bact_s) <- kgrow * bact_s - effect * bact_s - kdeath * bact_s - kgr * bact_s
    # dR/dt = kGR*G - kdeath*R
    d/dt(bact_r) <- kgr * bact_s - kdeath * bact_r
    bact_s(0) <- inoc

    log_cfu <- log10(btot + 1e-6)
    log_cfu ~ add(addSd)
  })
}
