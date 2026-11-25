Lei_2017b_marbofloxacin_pkpd_index <- function() {
  description <- paste0(
    "Veterinary (pig). Ex vivo inhibitory sigmoid Emax PK/PD-integration model for the ",
    "antibacterial effect of marbofloxacin against Streptococcus suis serotype 2 strain HB2 ",
    "incubated for 24 h in serum collected from pigs after a single 8 mg/kg oral dose. E is the ",
    "SIGNED change in log10 CFU/mL over the 24 h incubation (positive = net growth, negative = ",
    "net kill) and C is the PK/PD index AUC0-24h/MIC (h). Lei 2017b prints equation 1 as ",
    "E = Emax * (1 - C^N / (C^N + EC50^N)), which falls from Emax at C = 0 to 0 at saturation ",
    "and so can never reach the paper's own Table 3 targets E = -3 and -4; the fitted WinNonlin ",
    "inhibitory sigmoid model is E = Emax - (Emax - E0) * C^N / (C^N + EC50^N), whose lower ",
    "plateau E0 is shown in Figure 6 (about -7) but is not printed. Table 3 gives Emax = +1.67 ",
    "log10 CFU/mL (the drug-free change), EC50 = 33.82 h and N = 4.70. The lower plateau is ",
    "back-solved from the three Table 3 targets (E = 0, -3, -4 at 25.23, 35.64 and 39.71 h): ",
    "each target alone gives -6.62, -6.65 and -6.67, and the least-squares value -6.66 ",
    "reproduces all three to within 0.13%. The paper's drug-free Emax is mapped onto the ",
    "canonical `e0` and the maximal-kill plateau onto the canonical `emax`, written as ",
    "E = e0 - (e0 - emax) * C^N / (C^N + EC50^N). The index is formed as AUC_MARBOFLOXACIN / mic ",
    "with the HB2 serum MIC of 2 ug/mL; in the ex vivo design each serum sample holds one ",
    "concentration for the whole 24 h incubation, so AUC_MARBOFLOXACIN = 24 h times that ",
    "concentration. The bacterial density bact (CFU/mL) is integrated as d/dt(bact) = ",
    "ln(10) * (E / 24) * bact so that log10(bact) changes by exactly E across each 24 h window. ",
    "There is no PK component: the serum concentrations were analysed non-compartmentally in ",
    "WinNonlin and no structural PK model was published, so exposure enters as the covariate ",
    "AUC_MARBOFLOXACIN. Neither between-animal variability nor a residual error magnitude was ",
    "reported, so there are no etas and addSd is fixed at 0; the model is intended for ",
    "typical-value simulation."
  )
  reference <- paste(
    "Lei Z, Liu Q, Yang B, Khaliq H, Cao J, He Q. (2017).",
    "PK-PD analysis of marbofloxacin against Streptococcus suis in pigs.",
    "Frontiers in Pharmacology 8:856.",
    "doi:10.3389/fphar.2017.00856. PMCID: PMC5701813.",
    "Model equation 1, Table 3 and Figure 6; serum MIC of strain HB2 from Results.",
    sep = " "
  )
  vignette <- "Lei_2017b_marbofloxacin"
  units <- list(
    time = "h",
    dosing = "h*ug/mL (marbofloxacin ex vivo AUC over a 24 h window in serum, supplied as a covariate)",
    concentration = "log10 CFU/mL (observation)"
  )

  depends <- c("AUC_MARBOFLOXACIN")
  paper_specific_compartments <- c("bact")

  compartmentData <- list(
    bact = list(
      analyte = "Streptococcus suis serotype 2 strain HB2 (ex vivo, in pig serum)",
      units = "CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    AUC_MARBOFLOXACIN = list(
      description = "Marbofloxacin area under the serum concentration-time curve over the 24 h exposure window (AUC0-24h); in the ex vivo design this is the static sample concentration times 24 h",
      units = "h*ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Lei 2017b PK-PD Integration Modeling defines the sigmoid Emax driver C as the ex vivo ",
        "AUC0-24h/MIC (h). In the ex vivo design each serum sample, collected 0.25-24 h after a ",
        "single 8 mg/kg oral dose, holds one concentration for the 24 h incubation, so ",
        "AUC0-24h = 24 * C. The largest index in Figure 6 (about 75 h) matches the mean serum ",
        "Cmax: 6.28 * 24 / 2 = 75.4 h. Set to 0 for drug-free control serum so the sigmoid term ",
        "vanishes and the predicted 24 h change reduces to e0."
      ),
      source_name = "AUC0-24h (numerator of the PK/PD index AUC0-24h/MIC, Lei 2017b equation 1 and Table 3)"
    )
  )

  population <- list(
    species = "pig (healthy; ex vivo serum)",
    n_subjects = 8L,
    n_studies = 1L,
    age_range = "5-6 weeks (Animal Study); 8-10 weeks (Animals and Treatment)",
    weight_range = "16-20 kg (Animal Study); 15-20 kg (Animals and Treatment)",
    sex_female_pct = 50,
    disease_state = "Healthy pigs (serum donors); ex vivo challenge with Streptococcus suis",
    dose_range = "Single oral dose of marbofloxacin 8 mg/kg, then after a 2-week washout a single 8 mg/kg intravenous dose",
    regions = "China (Huazhong Agricultural University, Wuhan)",
    organism = "Streptococcus suis serotype 2 strain HB2, a virulent pig isolate chosen with an MIC close to the MIC90 of 134 clinical isolates; marbofloxacin MIC 1 ug/mL in TSB and 2 ug/mL in pig serum, MBC 4 ug/mL in both, MPC 2.56 ug/mL",
    notes = paste0(
      "Ex vivo time-kill: HB2 at 10^6 CFU/mL co-incubated with serum samples collected ",
      "0.25-24 h after the oral dose (Figure 3B). Serum NCA after the oral dose (Table 1, ",
      "mean +/- SD, n = 8): Cmax 6.28 +/- 0.92 ug/mL, Tmax 0.75 h, AUC0-24h 60.30 +/- 6.47 and ",
      "AUC0-inf 88.94 +/- 6.47 ug*h/mL, t1/2 12.48 h, Clb 0.104 +/- 0.02 L/h (per kg body ",
      "weight), F 94.21%; free fraction 92 +/- 4.3%. MIC distribution of 134 S. suis isolates ",
      "from Chinese pigs, 2016-2017 (Figure 1): MIC50 1 and MIC90 2 ug/mL. Monte Carlo ",
      "simulation (Crystal Ball, 10,000 trials) gave daily doses of 5.88, 8.34 and 9.36 mg/kg ",
      "for 50% and 8.16, 11.31 and 12.35 mg/kg for 90% target attainment of the ",
      "bacteriostatic, bactericidal and elimination targets (Table 4)."
    )
  )

  ini({
    # Lei 2017b equation 1 is printed as E = Emax * (1 - C^N / (C^N + EC50^N)),
    # which tends to 0, not to a maximal kill, and cannot reach the Table 3
    # targets E = -3 and -4. The fitted inhibitory sigmoid model is
    #   E = Emax - (Emax - E0) * C^N / (C^N + EC50^N)
    # with the paper's Emax as the drug-free change (Figure 6 starts near
    # +1.5) and the lower plateau E0 drawn near -7 but never printed.
    #
    # Canonical roles used in this file:
    #   e0   <- paper Emax = +1.67  (effect at zero exposure)
    #   emax <- unprinted plateau, back-solved = -6.66 (maximal drug effect)
    # E is signed, so e0 and emax stay on the natural scale.
    e0 <- 1.67
    label("Change in S. suis count over 24 h in drug-free serum (log10 CFU/mL)") # Lei 2017b Table 3, 'Emax' = 1.67 +/- 0.34 LgCFU/ml (the zero-exposure value of the sigmoid model)
    emax <- fixed(-6.66)
    label("Maximal change in S. suis count over 24 h (log10 CFU/mL; negative = reduction)") # Not printed in Lei 2017b. Back-solved from the Table 3 targets E = 0/-3/-4 at 25.23/35.64/39.71 h with Table 3 Emax/EC50/N: -6.62/-6.65/-6.67 individually, least squares -6.658; Figure 6 drawn plateau about -7.0
    lec50 <- log(33.82)
    label("AUC0-24h/MIC producing half the maximal effect EC50 (h)") # Lei 2017b Table 3, EC50 = 33.82 +/- 2.56 h
    lhill <- log(4.70)
    label("Hill coefficient N, steepness of the AUC0-24h/MIC effect curve (unitless)") # Lei 2017b Table 3, N = 4.70 +/- 0.82

    # Measured property of the challenge strain, not an estimate.
    mic <- fixed(2)
    label("Marbofloxacin MIC against S. suis HB2 in pig serum (ug/mL)") # Lei 2017b Results (MIC, MBC, and MPC Detection of MBF against HB2): MIC 2 ug/mL in pig serum; Table 2 ratios Cmax/MIC 3.14 = 6.28/2 and AUC0-24h/MIC 30.06 ~ 60.30/2

    # Experimental design input, not an estimate.
    log10_cfu0 <- fixed(6)
    label("log10 starting bacterial density in the ex vivo incubation (log10 CFU/mL)") # Lei 2017b Methods (Growth and Killing-Time Curves in Vitro and ex Vivo): bacteria (10^6 CFU/ml) cultivated with serum samples

    # No residual SD is reported for the PK/PD-integration fit.
    addSd <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (0; not reported in Lei 2017b)") # Lei 2017b reports only R2 = 0.989 for the fit (Results, PK-PD Integration Modeling)
  })

  model({
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # PK/PD index AUC0-24h/MIC (h)
    aucmic <- AUC_MARBOFLOXACIN / mic

    # Inhibitory sigmoid Emax model with canonical parameter roles: signed
    # 24 h change in log10 CFU/mL, e0 at zero exposure and emax at saturation.
    effect <- e0 - (e0 - emax) * aucmic^hill / (aucmic^hill + ec50^hill)

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
