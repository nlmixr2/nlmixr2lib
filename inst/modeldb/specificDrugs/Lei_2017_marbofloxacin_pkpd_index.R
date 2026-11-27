Lei_2017_marbofloxacin_pkpd_index <- function() {
  description <- paste0(
    "Veterinary (pig). Ex vivo inhibitory sigmoid Emax PK/PD-integration model for the ",
    "antibacterial effect of marbofloxacin against Escherichia coli strain HB197 incubated for ",
    "24 h in ileum content collected through an ileal cannula from pigs after a single 2 mg/kg ",
    "oral dose. Lei 2017 equation 1 is E = Emax - (Emax - E0) * C^N / (C^N + EC50^N), where E is ",
    "the SIGNED change in log10 CFU/mL over the 24 h incubation (positive = net growth, ",
    "negative = net kill) and C is the PK/PD index AUC0-24h/MIC (h). In that equation E equals ",
    "the paper's Emax at zero exposure and tends to the paper's E0 at saturating exposure, and ",
    "the Table 3 signs agree (Emax = +2.56, E0 = -6.06 log10 CFU/mL), so the paper's Emax is the ",
    "drug-free change and its E0 is the maximal kill, the reverse of the usual naming. This file ",
    "maps the paper's Emax = 2.56 onto the canonical `e0` (effect at zero exposure) and the ",
    "paper's E0 = -6.06 onto the canonical `emax` (maximal drug effect) and writes the ",
    "algebraically identical E = e0 - (e0 - emax) * C^N / (C^N + EC50^N); no value is altered. ",
    "EC50 = 19.93 h and N = 3.88 (Table 3). The parameterisation is confirmed against the ",
    "paper's own targets: solving E = 0, -3 and -4 gives AUC0-24h/MIC = 15.96, 23.25 and 26.86 h ",
    "against the published 16.26, 23.54 and 27.18 h (published values are means of per-animal ",
    "targets, so they need not coincide with the curve at the mean parameters). The index is ",
    "formed as AUC_MARBOFLOXACIN / mic with the HB197 MIC of 2 ug/mL, measured in both LB broth ",
    "and ileum content; in the ex vivo design each ileum-content sample holds one concentration ",
    "for the whole 24 h incubation, so AUC_MARBOFLOXACIN = 24 h times that concentration. The ",
    "bacterial density bact (CFU/mL) is integrated as d/dt(bact) = ln(10) * (E / 24) * bact so ",
    "that log10(bact) changes by exactly E across each 24 h window. There is no PK component: ",
    "Lei 2017 analysed the plasma and ileum-content concentrations non-compartmentally in ",
    "WinNonlin and published no structural PK model, so exposure enters as the covariate ",
    "AUC_MARBOFLOXACIN. Neither between-animal variability nor a residual error magnitude was ",
    "reported, so there are no etas and addSd is fixed at 0; the model is intended for ",
    "typical-value simulation."
  )
  reference <- paste(
    "Lei Z, Liu Q, Xiong J, Yang B, Yang S, Zhu Q, Li K, Zhang S, Cao J, He Q. (2017).",
    "Pharmacokinetic and pharmacodynamic evaluation of marbofloxacin and PK/PD modeling",
    "against Escherichia coli in pigs.",
    "Frontiers in Pharmacology 8:542.",
    "doi:10.3389/fphar.2017.00542. PMCID: PMC5566571.",
    "Model equation 1 and Table 3; MIC of strain HB197 from Results.",
    sep = " "
  )
  vignette <- "Lei_2017_marbofloxacin"
  units <- list(
    time = "h",
    dosing = "h*ug/mL (marbofloxacin ex vivo AUC over a 24 h window in ileum content, supplied as a covariate)",
    concentration = "log10 CFU/mL (observation)"
  )

  depends <- c("AUC_MARBOFLOXACIN")
  paper_specific_compartments <- c("bact")

  compartmentData <- list(
    bact = list(
      analyte = "Escherichia coli strain HB197 (ex vivo, in pig ileum content)",
      units = "CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    AUC_MARBOFLOXACIN = list(
      description = "Marbofloxacin area under the ileum-content concentration-time curve over the 24 h exposure window (AUC0-24h); in the ex vivo design this is the static sample concentration times 24 h",
      units = "h*ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Lei 2017 PK/PD Integration Analysis defines the sigmoid Emax driver C as the ex vivo ",
        "AUC0-24h/MIC (h). In the ex vivo design each ileum-content sample, collected 0-24 h ",
        "after a single 2 mg/kg oral dose, holds one concentration for the 24 h incubation, so ",
        "AUC0-24h = 24 * C. The paper reports ileum-content concentrations per gram (ug/g) and ",
        "the MIC per mL and treats the two as equivalent (density 1 g/mL). The largest index in ",
        "Figure 6 (about 136 h) matches the mean ileum-content Cmax: 11.28 * 24 / 2 = 135.4 h. ",
        "Set to 0 for drug-free control content so the sigmoid term vanishes and the predicted ",
        "24 h change reduces to e0."
      ),
      source_name = "AUC0-24h (numerator of the PK/PD index AUC0-24h/MIC, Lei 2017 equation 1 and Table 3)"
    )
  )

  population <- list(
    species = "pig (healthy, ileal-cannulated; ex vivo ileum content)",
    n_subjects = 8L,
    n_studies = 1L,
    age_range = "4-5 weeks",
    weight_range = "15-20 kg",
    sex_female_pct = 50,
    disease_state = "Healthy pigs fitted with an ileal cannula (content donors); ex vivo challenge with E. coli",
    dose_range = "Single oral (gavage) dose of marbofloxacin 2 mg/kg",
    regions = "China (Huazhong Agricultural University, Wuhan)",
    organism = "Escherichia coli strain HB197, a pathogenic pig isolate chosen with an MIC equal to the MIC90 of the 142 susceptible isolates; marbofloxacin MIC 2 ug/mL and MBC 4 ug/mL in both LB broth and pig ileum content, MPC 8.19 ug/mL",
    notes = paste0(
      "Ex vivo time-kill: HB197 at about 10^6 CFU/mL co-incubated with ileum-content samples ",
      "collected before and 0.25-24 h after the dose (Figure 3B). Ileum-content NCA (Table 1, ",
      "mean +/- SD, n = 8): Cmax 11.28 +/- 1.43 ug/g, AUC0-24h 46.15 +/- 3.62 and AUC0-inf ",
      "77.81 +/- 6.38 ug*h/g, Clb 0.026 kg/h (per kg body weight); plasma Cmax 0.55 ug/mL, ",
      "AUC0-inf 14.67 ug*h/mL, terminal half-life 30.67 h, CL/F 0.14 L/h/kg. MIC distribution of ",
      "218 E. coli isolates from five Chinese provinces (Figure 1): MIC50 0.5 and MIC90 2 ug/mL ",
      "among the 142 susceptible isolates. Table 3 parameters are reported as mean +/- SD. ",
      "Monte Carlo simulation (Crystal Ball, 10,000 trials) gave daily doses of 0.85, 1.22 and ",
      "1.41 mg/kg for 50% and 0.92, 1.33 and 1.53 mg/kg for 90% target attainment of the ",
      "bacteriostatic, bactericidal and elimination targets (Table 4)."
    )
  )

  ini({
    # Lei 2017 equation 1, as printed:
    #   E = Emax - (Emax - E0) * C^N / (C^N + EC50^N)
    # E equals the paper's Emax at C = 0 and tends to the paper's E0 as C
    # grows, and Table 3 gives Emax = +2.56 (net growth) and E0 = -6.06
    # (net kill). The Methods prose defines E0 and Emax together as 'the
    # changes ... in the control samples and those containing MBF samples',
    # and the equation, the Table 3 signs, Figure 6 (curve starts near
    # +2.6 and plateaus near -6) and the three printed targets all agree
    # that the paper's Emax is the drug-free change.
    #
    # Canonical roles used in this file:
    #   e0   <- paper Emax = +2.56  (effect at zero exposure)
    #   emax <- paper E0   = -6.06  (maximal drug effect)
    # and model() writes E = e0 - (e0 - emax) * C^N / (C^N + EC50^N).
    # E is signed, so e0 and emax stay on the natural scale.
    #
    # Check against Table 3's own targets: E = 0, -3, -4 solve to
    # C = 15.96, 23.25, 26.86 h (published 16.26, 23.54, 27.18 h).
    e0 <- 2.56
    label("Change in E. coli count over 24 h in drug-free ileum content (log10 CFU/mL)") # Lei 2017 Table 3, 'Emax' = 2.56 +/- 0.27 LgCFU/ml (the zero-exposure value of equation 1)
    emax <- -6.06
    label("Maximal change in E. coli count over 24 h (log10 CFU/mL; negative = reduction)") # Lei 2017 Table 3, 'E0' = -6.06 +/- 0.34 LgCFU/ml (the saturating-exposure value of equation 1); Table 3 'Emax-E0' = 8.59
    lec50 <- log(19.93)
    label("AUC0-24h/MIC producing half the maximal effect EC50 (h)") # Lei 2017 Table 3, EC50 = 19.93 +/- 2.12 h
    lhill <- log(3.88)
    label("Hill coefficient N, steepness of the AUC0-24h/MIC effect curve (unitless)") # Lei 2017 Table 3, N = 3.88 +/- 0.76

    # Measured property of the challenge strain, not an estimate.
    mic <- fixed(2)
    label("Marbofloxacin MIC against E. coli HB197 (ug/mL)") # Lei 2017 Results (MIC, MBC, and MPC of MBF against E. coli HB197): MIC 2 ug/mL in LB and in pig ileum content; Table 2 ratios Cmax/MIC 5.64 = 11.28/2 and AUC0-24h/MIC 23.08 = 46.15/2

    # Experimental design input, not an estimate.
    log10_cfu0 <- fixed(6)
    label("log10 starting bacterial density in the ex vivo incubation (log10 CFU/mL)") # Lei 2017 Methods (Bacterial Growth Curves ... Ex Vivo): bacteria (10^6 CFU/ml) co-incubated with ileum content samples

    # No residual SD is reported for the PK/PD-integration fit.
    addSd <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (0; not reported in Lei 2017)") # Lei 2017 reports only the observed-vs-predicted correlation 0.995 (Discussion)
  })

  model({
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # PK/PD index AUC0-24h/MIC (h)
    aucmic <- AUC_MARBOFLOXACIN / mic

    # Lei 2017 equation 1 with canonical parameter roles: signed 24 h
    # change in log10 CFU/mL, e0 at zero exposure and emax at saturation.
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
