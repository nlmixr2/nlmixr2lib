Lei_2018_florfenicol_pkpd_index <- function() {
  description <- paste0(
    "Ex vivo (pig plasma). Inhibitory sigmoid Emax PK/PD-integration model for the antibacterial ",
    "effect of florfenicol against Streptococcus suis serotype 2 strain WH-2 incubated for 24 h in ",
    "plasma sampled from healthy pigs after a single 30 mg/kg intramuscular dose. E is the SIGNED ",
    "change in log10 CFU/mL between 0 and 24 h of incubation (positive = net growth, negative = net ",
    "kill) and the PK/PD index is AUC0-24h/MIC (h). Lei 2018 Table 3 names the zero-exposure value ",
    "'Emax' (defined in its footnote as the change in the blank sample); it is carried here as the ",
    "canonical e0 = 2.36 log10 CFU/mL, with EC50 = 42.22 h and Hill N = 9.89 from the same table. ",
    "Table 3 does NOT print the maximal (drug-saturated) effect, and the printed Equation 1 is ",
    "garbled and cannot reach the paper's own bactericidal target, so the standard inhibitory form ",
    "E = e0 - (e0 - emax) * C^N / (C^N + EC50^N) is used with emax = -6.00 log10 CFU/mL digitised ",
    "by the maintainers from the lower asymptote of the fitted curve in Figure 6 (a four-parameter ",
    "fit to the digitised curve returns EC50 = 42.23 h and N = 9.97, confirming that the drawn curve ",
    "is the Table 3 fit). These values reproduce the published AUC0-24h/MIC targets of 37.89 ",
    "(bacteriostatic, E = 0), 44.02 (bactericidal, E = -3) and 46.42 h (eradication, E = -4) to ",
    "within 1.4-2.2%. The index is formed in model() as AUC_FLOR / mic with the strain WH-2 plasma ",
    "MIC of 2 ug/mL. The bacterial density bact (CFU/mL) is integrated as ",
    "d/dt(bact) = ln(10) * (E / 24) * bact so that log10(bact) changes by exactly E over each 24 h ",
    "window. There is no PK component: Lei 2018 analysed the plasma concentrations ",
    "non-compartmentally in WinNonlin 5.2.1 and published no compartmental model, so exposure ",
    "enters as the covariate AUC_FLOR (88.85 h*ug/mL after 30 mg/kg intramuscularly, Table 1). ",
    "Neither between-animal variability nor a residual error magnitude was reported, so there are ",
    "no eta parameters and addSd is FIXED at 0; the model is intended for typical-value simulation."
  )
  reference <- paste(
    "Lei Z, Liu Q, Yang S, Yang B, Khaliq H, Li K, Ahmed S, Sajid A, Zhang B, Chen P, Qiu Y,",
    "Cao J, He Q. (2018). PK-PD integration modeling and cutoff value of florfenicol against",
    "Streptococcus suis in pigs. Frontiers in Pharmacology 9:2.",
    "doi:10.3389/fphar.2018.00002. PMCID: PMC5776115.",
    "Model equation from Equation 1 (Methods) in its standard inhibitory sigmoid form;",
    "e0, EC50 and N from Table 3; the drug-saturated effect emax digitised from Figure 6;",
    "plasma MIC of strain WH-2 from Results.",
    sep = " "
  )
  vignette <- "Lei_2018_florfenicol"
  units <- list(
    time = "h",
    dosing = "h*ug/mL (florfenicol AUC over a 24 h window, supplied as a covariate)",
    concentration = "log10 CFU/mL (observation)"
  )

  depends <- c("AUC_FLOR")
  paper_specific_compartments <- c("bact")

  compartmentData <- list(
    bact = list(
      analyte = "Streptococcus suis serotype 2 strain WH-2",
      units = "CFU/mL",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    AUC_FLOR = list(
      description = "Florfenicol total plasma area under the concentration-time curve over a 24 h exposure window (AUC0-24h); in the ex vivo design this is the static sample concentration times 24 h",
      units = "h*ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Lei 2018 forms its PK/PD index as AUC0-24h/MIC (h) from TOTAL plasma florfenicol: Table 2's ",
        "AUC0-24h/MIC of 44.43 h is the Table 1 intramuscular AUC0-24h of 88.85 h*ug/mL divided by ",
        "the 2 ug/mL MIC, with no unbound-fraction correction (fu = 0.85 enters only the dose ",
        "equation). In the ex vivo design each plasma sample holds a fixed concentration C for the ",
        "24 h incubation, so AUC0-24h = 24 * C (C in ug/mL). Set to 0 for drug-free control plasma ",
        "so the sigmoid term vanishes and the predicted 24 h change reduces to e0."
      ),
      source_name = "AUC0-24h (numerator of the PK/PD index AUC0-24h/MIC, Lei 2018 Methods Equation 1, Tables 1-3)"
    )
  )

  population <- list(
    species = "pig (ex vivo plasma)",
    n_subjects = 8L,
    n_studies = 1L,
    age_range = "8-10 weeks (Animals section; the PK design section says 4-5 weeks)",
    weight_range = "15-20 kg",
    sex_female_pct = 50,
    disease_state = "Healthy pigs (serum donors); ex vivo challenge with Streptococcus suis",
    dose_range = "Single 30 mg/kg intramuscular dose of a 30% florfenicol injection (Wuhan Huisheng Biotechnology); the same pigs received 30 mg/kg intravenously after a 2-week washout for the bioavailability arm",
    regions = "China (Huazhong Agricultural University, Wuhan)",
    organism = "Streptococcus suis serotype 2 strain WH-2, chosen from 124 Chinese clinical isolates (MIC50 = 1, MIC90 = 2 ug/mL) because its MIC equals the MIC90; florfenicol MIC 2 ug/mL and MBC 4 ug/mL in both tryptone soy broth and pig plasma, mutant prevention concentration 3.2 ug/mL",
    notes = paste0(
      "Ex vivo time-kill: about 10^6 CFU/mL of WH-2 incubated for 24 h in plasma drawn at 0, 0.25, ",
      "0.5, 1, 2, 4, 6, 8, 10, 12 and 24 h after the intramuscular dose (Methods), counted on TSA ",
      "(limit of detection 10 CFU/mL). Plasma protein binding: free fraction 85 +/- 3.3%. ",
      "Non-compartmental PK (Table 1, mean +/- SD, n = 8): intramuscular Cmax 4.44 ug/mL at Tmax 4 h, ",
      "AUC0-24h 88.85 and AUC0-inf 158.56 h*ug/mL, CL/F 0.185 L/h/kg, terminal half-life 14.46 h, ",
      "F 99.14%; intravenous AUC0-inf 159.93 h*ug/mL, CL 0.187 L/h/kg, half-life 6.54 h. Targets ",
      "(Table 3): AUC0-24h/MIC 37.89 h (E = 0), 44.02 h (E = -3) and 46.42 h (E = -4). Predicted ",
      "doses (Table 4) for 50% / 90% target attainment: 16.50 / 21.55 (bacteriostatic), ",
      "19.17 / 25.02 (bactericidal) and 20.14 / 26.85 mg/kg (eradication)."
    )
  )

  ini({
    # Lei 2018 Equation 1 is garbled in the published text
    # ('E=Emax=(1-C^N/(C^N+EC50^N))'); read literally, E = Emax * (1 - f)
    # with Emax = 2.36 can never fall below zero and so cannot reach the
    # paper's own bactericidal (E = -3) and eradication (E = -4) targets.
    # The standard WinNonlin inhibitory sigmoid form is used instead:
    #   E = e0 - (e0 - emax) * C^N / (C^N + EC50^N)
    # with C = AUC0-24h/MIC. E is the SIGNED 24 h change in log10 CFU/mL,
    # so emax is negative and neither e0 nor emax is log-transformed.
    #
    # Lei 2018 Table 3 calls the zero-exposure value 'Emax' and its
    # footnote defines it as the change in the blank (drug-free) sample;
    # Figure 3B's blank curve rises from about 6.4 to 8.75 log10 CFU/mL,
    # matching it. It is the canonical e0 here; no value is changed.
    e0 <- 2.36
    label("Change in S. suis count over 24 h in drug-free plasma (log10 CFU/mL)") # Lei 2018 Table 3, 'Emax' = 2.36 +/- 0.48 Lg CFU/ml (footnote: change in bacterial counts of blank sample)

    # NOT PRINTED in Lei 2018. Digitised by the maintainers from the lower
    # asymptote of the fitted ('Predicted') curve in Figure 6. A free
    # four-parameter sigmoid fitted to the digitised curve returns
    # top = 2.55, bottom = -6.00, EC50 = 42.23 h, N = 9.97; the EC50 and
    # N agree with Table 3 to under 1%, so the drawn curve is the Table 3
    # fit and its bottom is the missing maximal effect.
    emax <- -6.00
    label("Maximal change in S. suis count over 24 h (log10 CFU/mL; negative = reduction)") # Digitised from Lei 2018 Figure 6, lower asymptote of the 'Predicted' curve (not printed in Table 3)

    lec50 <- log(42.22)
    label("AUC0-24h/MIC producing 50% of the maximal effect EC50 (h)") # Lei 2018 Table 3, EC50 = 42.22 +/- 1.20 h
    lhill <- log(9.89)
    label("Hill coefficient N defining the steepness of the AUC0-24h/MIC effect curve (unitless)") # Lei 2018 Table 3, N = 9.89 +/- 1.15

    # Measured property of the challenge strain, not an estimate. Change it
    # to apply the model to an isolate of different susceptibility.
    mic <- fixed(2)
    label("Florfenicol MIC against S. suis WH-2 in pig plasma (ug/mL)") # Lei 2018 Results, 'MIC, MBC, and MPC of WH-2': MIC 2 ug/mL in TSB and plasma (equal to the MIC90 of 124 isolates)

    # Experimental design input, not an estimate.
    log10_cfu0 <- fixed(6)
    label("log10 starting S. suis density in the ex vivo incubation (log10 CFU/mL)") # Lei 2018 Methods, ex vivo time-killing curves: bacteria at 10^6 CFU/ml co-incubated with serum samples

    # Lei 2018 reports only R^2 = 0.999 for the sigmoid fit (Discussion)
    # and no residual standard deviation.
    addSd <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (not reported in Lei 2018)") # Lei 2018 reports R^2 only, no residual SD
  })

  model({
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # PK/PD index AUC0-24h/MIC (h), on total plasma florfenicol.
    aucmic <- AUC_FLOR / mic

    # Inhibitory sigmoid Emax; signed 24 h change in log10 CFU/mL. Equals
    # e0 at zero exposure and approaches emax at saturating exposure.
    eff_log10 <- e0 - (e0 - emax) * aucmic^hill / (aucmic^hill + ec50^hill)

    # Spread the 24 h change uniformly so log10(bact) moves by exactly
    # eff_log10 over each 24 h window:
    #   d(log10 N)/dt = eff_log10 / 24  =>  dN/dt = ln(10) * (eff_log10 / 24) * N
    d/dt(bact) <- log(10) * (eff_log10 / 24) * bact
    bact(0) <- 10^log10_cfu0

    # log10 CFU/mL with a 1 CFU/mL floor so the log stays finite.
    log_cfu <- log10(bact + 1)
    log_cfu ~ add(addSd)
  })
}
