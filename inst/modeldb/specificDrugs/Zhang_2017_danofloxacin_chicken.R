Zhang_2017_danofloxacin_chicken <- function() {
  description <- paste0(
    "Preclinical (chicken, Mycoplasma gallisepticum-infected). Sigmoid Emax PK/PD-index model for ",
    "the in vivo antimycoplasmal effect of oral danofloxacin against M. gallisepticum strain S6 in ",
    "an intratracheal chicken infection model. The effect E is the reduction in lung M. gallisepticum ",
    "load (log10 CFU/lung) relative to the untreated control group, measured 24 h after each of three ",
    "once-daily oral doses (POSITIVE = reduction), and the driver is the lung-tissue AUC24/MIC index. ",
    "Zhang 2017 writes E = Emax - (Emax - E0) * Ce^N / (EC50^N + Ce^N) with INVERTED names: its Emax ",
    "(0.72) is the effect in the untreated control and its E0 (3.70) is the maximal reduction. This file ",
    "maps them to the canonical e0 (effect at zero exposure) and emax (maximal effect); no value is ",
    "changed. Exposure enters as the lung AUC over the 24 h dosing interval (AUC_DANOFLOXACIN_LUNG) divided by the ",
    "fixed strain-S6 MIC of 0.15 mg/L, so the model can be re-targeted to another isolate by changing mic. ",
    "Zhang 2017 reports pharmacokinetics by non-compartmental analysis only (WinNonlin), so there is no ",
    "PK model; supply AUC_DANOFLOXACIN_LUNG from the paper's Table 5 or another danofloxacin PK model."
  )
  reference <- paste(
    "Zhang N, Wu Y, Huang Z, Yao L, Zhang L, Cai Q, Shen X, Jiang H, Ding H. (2017).",
    "The PK-PD Relationship and Resistance Development of Danofloxacin against",
    "Mycoplasma gallisepticum in An In Vivo Infection Model.",
    "Frontiers in Microbiology 8:926. doi:10.3389/fmicb.2017.00926. PMCID: PMC5447713.",
    "Model equation: Materials and Methods, 'PK-PD Analysis'. Parameter estimates: Table 6.",
    "MIC: Table 2 (strain S6). Lung AUC24: Table 5. Target indices and dose: Results, 'PK and PD",
    "Analysis' and 'Dosage Calculation'.",
    sep = " "
  )
  vignette <- "Zhang_2017_danofloxacin_chicken"

  units <- list(
    time = "h",
    dosing = "ug*h/g (lung danofloxacin AUC over the 24 h dosing interval, supplied as a covariate)",
    concentration = "log10 CFU/lung (observation: reduction relative to the untreated control)"
  )

  covariateData <- list(
    AUC_DANOFLOXACIN_LUNG = list(
      description = "Danofloxacin area under the lung-tissue concentration-time curve over the 24 h dosing interval preceding the effect measurement (AUC24)",
      units = "ug*h/g",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Zhang 2017 derived lung AUC24 non-compartmentally (WinNonlin 6.1) from pooled lung ",
        "samples (n = 8 chickens per time point) after once-daily oral danofloxacin at 1, 2.5, 5, ",
        "10 and 20 mg/kg (Table 5: 2.72, 4.56, 10.26, 21.99 and 58.58; the table's unit header ",
        "reads ug*h/ml although the matrix is lung tissue). Divided inside model() by the strain ",
        "MIC (mg/L) to form the AUC24/MIC index of Table 6 and Figure 5, whose unit Zhang 2017 ",
        "prints as h (Figure 5, Table 6) or L*h/kg (Results text). Set to 0 for the untreated ",
        "control so the sigmoid term vanishes and the prediction reduces to e0."
      ),
      source_name = "AUC24 (Zhang 2017 Table 5, lung tissues; Materials and Methods, 'PK-PD Analysis')"
    )
  )

  population <- list(
    species = "chicken (Mycoplasma gallisepticum-infected; breed not reported)",
    n_subjects = 210L,
    n_studies = 1L,
    age_range = "obtained at one day old (Guangdong Academy of Agricultural Sciences); age at infection not reported",
    organism = "Mycoplasma gallisepticum strain S6 (Chinese Veterinary Microorganism Culture Collection Center); danofloxacin MIC 0.15 mg/L by agar dilution (Table 2)",
    system = "Intratracheal infection: 0.2 mL of 10^9 CFU/mL M. gallisepticum twice daily for 3 days. Air sacs and lungs harvested 24 h after each daily dose (n = 10 per group per time point) and plated by 10-fold serial dilution",
    disease_state = "experimentally induced M. gallisepticum chronic respiratory disease",
    dose_range = "danofloxacin 1, 2.5, 5, 7.5, 10, 15 or 20 mg/kg by oral gavage once daily for 3 days (PD cohort); 1, 2.5, 5, 10 or 20 mg/kg in the PK cohort",
    regions = "China (South China Agricultural University, Guangzhou)",
    notes = "The 210 chickens are the pharmacodynamic cohort; a separate 840-chicken cohort supplied the non-compartmental plasma (Table 4) and lung (Table 5) pharmacokinetics (n = 8 per sampling time). 40 infected chickens given saline formed the growth control. The sigmoid was fitted in WinNonlin to the lung effect only (R^2 = 0.8901); air-sac counts (Figure 3) were not modelled."
  )

  ini({
    # =================================================================
    # Zhang 2017 Table 6 -- sigmoid Emax PK/PD-index model (lung)
    # =================================================================
    # Equation as printed (Materials and Methods, 'PK-PD Analysis'):
    #   E = Emax - (Emax - E0) * Ce^N / (EC50^N + Ce^N)
    # with Zhang 2017's own definitions, which INVERT the usual roles:
    #   Emax: 'the reduction in log10 CFU/lung for the untreated control'
    #   E0:   'the maximum reduction after administration'
    # Figure 5 confirms the inversion: the fitted curve starts near 0.7
    # at AUC24/MIC = 0 and plateaus near 3.5.
    #
    # Canonical roles used here:
    #   e0   <- Zhang 2017 Emax = 0.72  (effect at zero exposure)
    #   emax <- Zhang 2017 E0   = 3.70  (maximal effect)
    # model() writes the algebraically identical
    #   E = e0 - (e0 - emax) * Ce^N / (Ce^N + EC50^N)
    #
    # E is a POSITIVE reduction relative to the untreated control, so
    # e0 and emax stay on the natural scale; EC50 and N are strictly
    # positive and carried as logs.
    e0 <- 0.72
    label("Reduction in lung M. gallisepticum load versus untreated control at zero exposure (log10 CFU/lung)") # Zhang 2017 Table 6, 'Emax (log10 CFU/ml)' = 0.72 +/- 0.50 (defined in 'PK-PD Analysis' as the untreated-control reduction)
    emax <- 3.70
    label("Maximal reduction in lung M. gallisepticum load versus untreated control (log10 CFU/lung)") # Zhang 2017 Table 6, 'E0 (log10 CFU/ml)' = 3.70 +/- 0.84 (defined in 'PK-PD Analysis' as the maximum reduction)
    lec50 <- log(115.02)
    label("Log AUC24/MIC producing 50% of the maximal effect EC50 (h)") # Zhang 2017 Table 6, 'EC50 (h)' = 115.02 +/- 39.22
    lhill <- log(1.77)
    label("Log Hill coefficient N (unitless)") # Zhang 2017 Table 6, 'Slope (N)' = 1.77 +/- 1.34

    # =================================================================
    # Minimum inhibitory concentration
    # =================================================================
    # Measured property of the challenge strain, not an estimate. The
    # paper forms AUC24/MIC from the in vitro MIC of strain S6, and its
    # dose equation uses the same 0.15 mg/L as the projected MIC90.
    mic <- fixed(0.15)
    label("Danofloxacin MIC against M. gallisepticum strain S6 (mg/L)") # Zhang 2017 Table 2, S6 / Danofloxacin = 0.15 mg/L; 'Dosage Calculation' uses it as the projected MIC90

    # =================================================================
    # Residual error
    # =================================================================
    # Zhang 2017 reports the parameter estimates with SDs and R^2 =
    # 0.8901 but no residual SD, so it is held at zero.
    addSd <- fixed(0)
    label("Additive residual SD on the reduction in log10 CFU/lung (not reported)") # Zhang 2017 reports R^2 = 0.8901 only ('PK and PD Analysis'), no residual SD
  })

  model({
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # Lung AUC24/MIC index (h); AUC_DANOFLOXACIN_LUNG in ug*h/g, mic in mg/L.
    ce <- AUC_DANOFLOXACIN_LUNG / mic

    # Zhang 2017 sigmoid Emax with canonical roles. Cc is the reduction
    # in lung M. gallisepticum load relative to the untreated control,
    # 24 h after a dose (log10 CFU/lung; positive = reduction).
    Cc <- e0 - (e0 - emax) * ce^hill / (ce^hill + ec50^hill)
    Cc ~ add(addSd)
  })
}
