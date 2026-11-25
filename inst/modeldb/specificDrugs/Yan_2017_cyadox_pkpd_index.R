Yan_2017_cyadox_pkpd_index <- function() {
  description <- paste0(
    "Veterinary (pig). Ex vivo inhibitory sigmoid Emax PK/PD-integration model for the ",
    "antibacterial effect of cyadox, a quinoxaline-1,4-dioxide, against Clostridium perfringens ",
    "strain CPFK122995 incubated for 24 h in free (ultrafiltered) ileum fluid collected from ",
    "pigs after a single 30 mg/kg oral dose. Yan 2017 equation (1) is ",
    "E = Emax - (Emax - E0) * C^N / (C^N + EC50^N), where E is the SIGNED change in log10 ",
    "CFU/mL over the 24 h incubation (positive = net growth, negative = net kill) and C is the ",
    "PK/PD index AUC0-24h/MIC (h). In that equation E equals the paper's Emax at zero exposure ",
    "and tends to the paper's E0 at saturating exposure, and the Table 4 signs agree ",
    "(Emax = +2.71, E0 = -4.61 log10 CFU/mL), so the paper's Emax is the drug-free change and ",
    "its E0 is the maximal kill, the reverse of the usual naming. This file maps the paper's ",
    "Emax = 2.71 onto the canonical `e0` (effect at zero exposure) and the paper's E0 = -4.61 ",
    "onto the canonical `emax` (maximal drug effect) and writes the algebraically identical ",
    "E = e0 - (e0 - emax) * C^N / (C^N + EC50^N); no value is altered. EC50 = 30.07 h and ",
    "N = 4.62 (Table 4). The parameterisation is confirmed against the paper's own targets: ",
    "solving E = 0, -3 and -4 gives AUC0-24h/MIC = 26.80, 39.55 and 50.53 h against the ",
    "published 26.72, 39.54 and 50.69 h. The index is formed as AUC_CYADOX / mic with the ",
    "CPFK122995 MIC of 2 ug/mL, which Yan 2017 measured in both Brucella broth and ileum fluid; ",
    "in the ex vivo design each ileum-fluid sample holds one concentration for the whole 24 h ",
    "incubation, so AUC_CYADOX = 24 h times that concentration. The bacterial density bact ",
    "(CFU/mL) is integrated as d/dt(bact) = ln(10) * (E / 24) * bact so that log10(bact) ",
    "changes by exactly E across each 24 h window. There is no PK component: Yan 2017 analysed ",
    "the ileum-fluid and plasma concentrations non-compartmentally in WinNonlin and published ",
    "no structural PK model, so exposure enters as the covariate AUC_CYADOX. Neither ",
    "between-animal variability nor a residual error magnitude was reported, so there are no ",
    "etas and addSd is fixed at 0; the model is intended for typical-value simulation."
  )
  reference <- paste(
    "Yan L, Xie S, Chen D, Pan Y, Tao Y, Qu W, Liu Z, Yuan Z, Huang L. (2017).",
    "Pharmacokinetic and pharmacodynamic modeling of cyadox against Clostridium perfringens",
    "in swine.",
    "Scientific Reports 7:4064.",
    "doi:10.1038/s41598-017-03970-9. PMCID: PMC5481453.",
    "Model equation (1) and Table 4; MIC of strain CPFK122995 from Results.",
    sep = " "
  )
  vignette <- "Yan_2017_cyadox"
  units <- list(
    time = "h",
    dosing = "h*ug/mL (cyadox ex vivo AUC over a 24 h window, supplied as a covariate)",
    concentration = "log10 CFU/mL (observation)"
  )

  depends <- c("AUC_CYADOX")
  paper_specific_compartments <- c("bact")

  compartmentData <- list(
    bact = list(
      analyte = "Clostridium perfringens strain CPFK122995 (ex vivo, in ileum-fluid ultrafiltrate)",
      units = "CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    AUC_CYADOX = list(
      description = "Cyadox area under the free ileum-fluid concentration-time curve over the 24 h exposure window (AUC0-24h); in the ex vivo design this is the static sample concentration times 24 h",
      units = "h*ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Yan 2017 Materials and Methods (Pharmacokinetic-pharmacodynamic integration and ",
        "modeling) defines the sigmoid Emax driver C as the ex vivo AUC0-24h/MIC (h). In the ex ",
        "vivo design each ileum-fluid sample, collected 0-24 h after a single 30 mg/kg oral dose, ",
        "holds one concentration for the 24 h incubation, so AUC0-24h = 24 * C (C in ug/mL). ",
        "The largest index in Figure 6 (about 285 h) matches the mean ileum-fluid Cmax: ",
        "23.66 * 24 / 2 = 283.9 h. Set to 0 for drug-free control fluid so the sigmoid term ",
        "vanishes and the predicted 24 h change reduces to e0."
      ),
      source_name = "AUC0-24h (numerator of the PK/PD index AUC0-24h/MIC, Yan 2017 equation (1) and Table 4)"
    )
  )

  population <- list(
    species = "pig (crossbred Duroc x Large White x Landrace, castrated; ex vivo ileum fluid)",
    n_subjects = 6L,
    n_studies = 1L,
    age_range = NA_character_,
    weight_range = "15-20 kg",
    sex_female_pct = 0,
    disease_state = "Healthy pigs with a surgically implanted ileum ultrafiltration probe (fluid donors); ex vivo challenge with C. perfringens",
    dose_range = "Single oral dose of cyadox 30 mg/kg (50 mL carboxymethyl cellulose suspension)",
    regions = "China (Huazhong Agricultural University, Wuhan)",
    organism = "Clostridium perfringens strain CPFK122995, a pathogenic pig isolate from Guangdong Province chosen with an MIC near the MIC50; cyadox MIC 2 ug/mL in both Brucella broth and ileum fluid, MBC 4 ug/mL (broth) and 8 ug/mL (ileum fluid), MPC 12 ug/mL",
    notes = paste0(
      "Ex vivo time-kill: CPFK122995 at about 10^6 CFU/mL incubated anaerobically at 37 C in ",
      "ileum-fluid samples collected 0-24 h after the dose, counted at 1, 2, 4, 8, 12 and 24 h ",
      "(limit of detection 10 CFU/mL). Ileum-fluid NCA (Table 2, mean +/- SD, n = 6): ",
      "Cmax 23.66 +/- 0.92 ug/mL, Tmax 1.96 h, terminal half-life 5.86 h, AUC0-24h ",
      "106.40 +/- 12.68 and AUC0-inf 112.35 +/- 14.07 h*ug/mL; CL/F 0.27 L/kg/h (Abstract). ",
      "MIC distribution of 60 clinical pig isolates (Figure 2): 0.5-8 ug/mL, MIC50 2 and MIC90 ",
      "4 ug/mL. Table 4 parameters are reported as mean +/- SD. Monte Carlo simulation ",
      "(Crystal Ball, 100,000 trials) gave daily doses of 13.2, 19.70 and 25.46 mg/kg for 50% ",
      "and 29.30, 42.56 and 54.50 mg/kg for 90% target attainment of the bacteriostatic, ",
      "bactericidal and elimination targets (Table 5)."
    )
  )

  ini({
    # Yan 2017 equation (1), as printed:
    #   E = Emax - (Emax - E0) * C^N / (C^N + EC50^N)
    # E equals the paper's Emax at C = 0 and tends to the paper's E0 as C
    # grows, and Table 4 gives Emax = +2.71 (net growth) and E0 = -4.61
    # (net kill). So the paper's Emax is the drug-free 24 h change and its
    # E0 is the maximal kill. The Methods prose describes the two the
    # other way round; the equation, the Table 4 signs, Figure 6 (curve
    # starts near +2.7 and plateaus near -4.6) and the three printed
    # targets all agree with the reading used here.
    #
    # Canonical roles used in this file:
    #   e0   <- paper Emax = +2.71  (effect at zero exposure)
    #   emax <- paper E0   = -4.61  (maximal drug effect)
    # and model() writes E = e0 - (e0 - emax) * C^N / (C^N + EC50^N).
    # E is signed, so e0 and emax stay on the natural scale.
    #
    # Check against Table 4's own targets: E = 0, -3, -4 solve to
    # C = 26.80, 39.55, 50.53 h (published 26.72, 39.54, 50.69 h).
    e0 <- 2.71
    label("Change in C. perfringens count over 24 h in drug-free ileum fluid (log10 CFU/mL)") # Yan 2017 Table 4, 'Emax' = 2.71 +/- 0.37 Log CFU/mL (the zero-exposure value of equation (1))
    emax <- -4.61
    label("Maximal change in C. perfringens count over 24 h (log10 CFU/mL; negative = reduction)") # Yan 2017 Table 4, 'E0' = -4.61 +/- 0.69 Log CFU/mL (the saturating-exposure value of equation (1)); Table 4 'Emax-E0' = 7.32 = 2.71 + 4.61
    lec50 <- log(30.07)
    label("AUC0-24h/MIC producing half the maximal effect EC50 (h)") # Yan 2017 Table 4, EC50 = 30.07 +/- 6.94 h
    lhill <- log(4.62)
    label("Hill coefficient N, steepness of the AUC0-24h/MIC effect curve (unitless)") # Yan 2017 Table 4, N (slope) = 4.62 +/- 0.72

    # Measured property of the challenge strain, not an estimate.
    mic <- fixed(2)
    label("Cyadox MIC against C. perfringens CPFK122995 (ug/mL)") # Yan 2017 Results: MIC 2 ug/mL in Brucella broth and in pig ileum fluid; also the MIC50 of the 60 isolates

    # Experimental design input, not an estimate.
    log10_cfu0 <- fixed(6)
    label("log10 starting bacterial density in the ex vivo incubation (log10 CFU/mL)") # Yan 2017 Methods (In vitro and ex vivo time-kill study): starting inoculum 10^6 CFU/mL; ex vivo curves 'determined as described above'

    # No residual SD is reported for the PK/PD-integration fit.
    addSd <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (0; not reported in Yan 2017)") # Yan 2017 reports R^2 = 0.999 only (Discussion)
  })

  model({
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # PK/PD index AUC0-24h/MIC (h)
    aucmic <- AUC_CYADOX / mic

    # Yan 2017 equation (1) with canonical parameter roles: signed 24 h
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
