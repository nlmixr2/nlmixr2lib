Qi_2026_nemonoxacin_camhb_cmaxmic_aubkc <- function() {
  description <- paste0(
    "In vitro (Escherichia coli, one-compartment dynamic urinary tract infection model). Sigmoid ",
    "Emax PK/PD-index model for nemonoxacin against E. coli in cation-adjusted Mueller-Hinton ",
    "broth (CAMHB): Qi 2026 Table 2 fit of the area under the log10 CFU/mL bacterial kill curve ",
    "over 0-24 h (AUBKC) against Cmax/MIC, the peak nemonoxacin concentration in the medium ",
    "divided by the standard CAMHB broth-microdilution MIC of the isolate. E = E0 - Emax * X^H / ",
    "(X^H + EC50^H) with E0 = 208 log10 CFU/mL*h, Emax = 120 log10 CFU/mL*h, EC50 = 2.45 and H = ",
    "1.89 (R^2 = 0.934). The model is algebraic: the output auc_log_cfu is the predicted 0-24 h ",
    "AUBKC itself. There is no PK component: exposure enters as the covariates CMAX_NEMONOXACIN ",
    "and MIC, which the model divides to form the index. Neither between-isolate variability nor ",
    "a residual error magnitude was reported. One of 17 sibling fits from Qi 2026 sharing the ",
    "vignette Qi_2026_nemonoxacin."
  )
  reference <- paste0(
    "Qi C, Zhang R, Qu X, Liu X, Bian X, Zhu F, Chen J, Li X, Zhang J. Nemonoxacin ",
    "pharmacokinetics/pharmacodynamics against Escherichia coli in an in vitro dynamic urinary ",
    "tract infection model. Antimicrob Agents Chemother. 2026;70(7):e01368-25. ",
    "doi:10.1128/aac.01368-25"
  )
  vignette <- "Qi_2026_nemonoxacin"
  units <- list(
    time = "h",
    dosing = "none (nemonoxacin exposure is supplied as covariates)",
    concentration = "log10 CFU/mL*h (AUBKC observation)"
  )

  depends <- c("CMAX_NEMONOXACIN", "MIC")

  covariateData <- list(
    CMAX_NEMONOXACIN = list(
      description = "Nemonoxacin peak concentration in the medium of the in vitro dynamic urinary tract infection model, used as the numerator of the Cmax/MIC index",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Qi 2026 read Cmax non-compartmentally (Phoenix WinNonlin) from the central-flask concentrations measured by LC-MS/MS. The flask profile reproduced day-10 urine concentrations after oral 500 mg once daily (Table 1), whose 0-2 h interval concentration is 150.76 mg/L (minimum, low exposure) and 535.34 mg/L (median, high exposure). The paper does not tabulate the per-experiment Cmax values. Set to 0 for drug-free growth controls.",
      source_name = "Cmax (Qi 2026 Methods, PK/PD index and targets; numerator of the Cmax/MIC and Cmax/MICSHU indices of Table 2)"
    ),
    MIC = list(
      description = "Nemonoxacin MIC against the challenge Escherichia coli isolate by standard CLSI broth microdilution in CAMHB (median of triplicate)",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Qi 2026 Table S3, column CAMHB / NEN: 0.5-128 mg/L across the clinical isolates and below 0.03 mg/L for ATCC 25922. This is the denominator of the Cmax/MIC and AUC0-24/MIC indices in BOTH media; the MIC measured in synthetic human urine is the separate column MIC_SYNTHETIC_URINE. Must be positive.",
      source_name = "MIC (Qi 2026 Methods, MIC determination in CAMHB and SHU; Table S3 column CAMHB / NEN)"
    )
  )

  population <- list(
    species = "in vitro (Escherichia coli clinical isolates)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    organism = "Escherichia coli: 12 clinical isolates (5 producing extended-spectrum beta-lactamase) plus ATCC 25922 for susceptibility testing (Table S3; nemonoxacin CAMHB MIC 0.5-128 mg/L for the clinical isolates, below 0.03 mg/L for ATCC 25922), and NTU-001 and NTU-002 added to the high-exposure SHU arm",
    system = "One-compartment open dilution model after Satta et al. (Figure 1): a central flask holding the bacterial suspension (inoculated at 1e6 CFU/mL and preincubated 1 h at 35 C), a drug reservoir, a drug-free medium reservoir and a waste flask, driven by a Masterflex L/S peristaltic pump at the Table S2 flow rates; drug concentrations and viable counts sampled at 0, 0.5, 1, 2, 4, 6, 8 and 24 h",
    medium = "cation-adjusted Mueller-Hinton broth (CAMHB; BD, USA)",
    experimental_arms = "low urine exposure only (Figure 5a strains W035-079, W049-027, W042-007 and W049-099) plus drug-free growth controls",
    disease_state = "in vitro model of uncomplicated urinary tract infection",
    dose_range = "simulated day-10 urine exposure after oral nemonoxacin 500 mg once daily in healthy Chinese volunteers: low exposure = Table 1 minimum and high exposure = Table 1 median interval urine concentrations",
    regions = "China (Huashan Hospital, Fudan University, Shanghai)",
    notes = "Qi 2026 fitted the sigmoid Emax model separately for each combination of medium (CAMHB, SHU), PK/PD index (Cmax/MIC, AUC0-24/MIC, %T>MIC and, in SHU only, the same three against the SHU MIC) and endpoint (Delta log10 CFU, AUBKC); Table 2 reports 17 fits (the AUC0-24/MICSHU x Delta log10 CFU fit is absent from both Table 2 and Figure 6), packaged as 17 files. AUC0-24/MIC was the index best correlated with effect (Delta log10 CFU R^2 0.951 in CAMHB, 0.982 in SHU). Parameters were fitted with naive pooling of all isolates; the number of experiments behind each fit is not stated."
  )

  ini({
    # Qi 2026 Table 2, CAMHB, Cmax/MIC, AUBKC; R^2 = 0.934.
    # Methods eq: E = E0 - Emax * X^H / (X^H + EC50^H); the Table 2 column
    # headed gamma is the Hill coefficient H.
    le0 <- log(208)
    label("AUBKC over 0-24 h at zero exposure, the drug-free E0 (log10 CFU/mL*h)") # Qi 2026 Table 2, CAMHB, Cmax/MIC, AUBKC, E0 = 208

    lemax <- log(120)
    label("Maximum drug effect Emax, the largest attainable reduction of AUBKC (log10 CFU/mL*h)") # Qi 2026 Table 2, CAMHB, Cmax/MIC, AUBKC, Emax = 120

    lec50 <- log(2.45)
    label("Cmax/MIC producing half the maximum effect EC50 (unitless)") # Qi 2026 Table 2, CAMHB, Cmax/MIC, AUBKC, EC50 = 2.45

    lhill <- log(1.89)
    label("Hill coefficient H of the index-effect curve (unitless)") # Qi 2026 Table 2, CAMHB, Cmax/MIC, AUBKC, gamma = 1.89

    # Table 2 reports R^2 only, no residual SD and no standard errors, so the
    # residual SD is held at zero for deterministic typical-value simulation.
    addSd <- fixed(0)
    label("Additive residual SD on AUBKC (log10 CFU/mL*h)") # not reported in Qi 2026 (Table 2 gives R^2 only)
  })

  model({
    e0 <- exp(le0)
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # PK/PD index Cmax/MIC: exposure divided by the isolate MIC.
    pkpd_index <- CMAX_NEMONOXACIN / MIC

    # Qi 2026 Methods sigmoid Emax equation; the effect equals e0 at zero
    # exposure and approaches e0 - emax at saturating exposure.
    effect <- e0 - emax * pkpd_index^hill / (pkpd_index^hill + ec50^hill)

    # AUBKC is the area under the log10 CFU/mL curve over 0-24 h, a single
    # per-experiment value, so the output is the effect itself.
    auc_log_cfu <- effect
    auc_log_cfu ~ add(addSd)
  })
}
