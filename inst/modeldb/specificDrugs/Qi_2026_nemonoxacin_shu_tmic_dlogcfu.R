Qi_2026_nemonoxacin_shu_tmic_dlogcfu <- function() {
  description <- paste0(
    "In vitro (Escherichia coli, one-compartment dynamic urinary tract infection model). Sigmoid ",
    "Emax PK/PD-index model for nemonoxacin against E. coli in synthetic human urine (SHU, pH ",
    "6.5): Qi 2026 Table 2 fit of the change in log10 CFU/mL between 0 and 24 h (Delta log10 CFU) ",
    "against %T>MIC, the percentage of the 24 h interval during which the nemonoxacin ",
    "concentration exceeds the standard CAMHB broth-microdilution MIC of the isolate. E = E0 - ",
    "Emax * X^H / (X^H + EC50^H) with E0 = 2.70 log10 CFU/mL, Emax = 66.6 log10 CFU/mL, EC50 = ",
    "129 percent (Table 2 prints 1.29 on the 0-1 fraction scale; multiplied by 100 to the percent ",
    "scale of FTMIC_NEMONOXACIN) and H = 9.95 (R^2 = 0.735). The fitted EC50 exceeds 100 percent, ",
    "so only part of the sigmoid is reachable and Emax is an extrapolation far beyond the data ",
    "(R^2 = 0.735, the weakest fit in Table 2). The 24 h change is spread uniformly across the 24 ",
    "h window from the nominal 1e6 CFU/mL inoculum, d/dt(bact) = ln(10) * (E / 24) * bact, so ",
    "log10(bact) moves by exactly E at 24 h; the intermediate time course is a construction, not ",
    "the observed kill curve (which falls within 2-4 h and then regrows, Figure 5). There is no ",
    "PK component: exposure enters as the covariate FTMIC_NEMONOXACIN. Neither between-isolate ",
    "variability nor a residual error magnitude was reported. One of 17 sibling fits from Qi 2026 ",
    "sharing the vignette Qi_2026_nemonoxacin."
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
    concentration = "log10 CFU/mL (observation)"
  )

  depends <- c("FTMIC_NEMONOXACIN")

  paper_specific_compartments <- c("bact")

  compartmentData <- list(
    bact = list(analyte = "Escherichia coli", units = "CFU/mL", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    FTMIC_NEMONOXACIN = list(
      description = "Percentage of the 24 h interval during which the nemonoxacin concentration in the medium exceeds the CAMHB broth-microdilution MIC of the isolate (%T>MIC)",
      units = "percent (0-100)",
      type = "continuous",
      reference_category = NULL,
      notes = "Qi 2026 fits this index on the 0-1 fraction scale (Figure 6 x-axes run 0-1 and the Table 2 EC50 values are fractions); the packaged EC50 is the Table 2 value multiplied by 100 so this column is supplied in percent, matching the FTMIC_<DRUG> family convention. The rescaling is exact because the sigmoid term depends only on the index divided by EC50. Set to 0 for drug-free growth controls.",
      source_name = "%T > MIC (Qi 2026 Methods, PK/PD index and targets; Table 2; Figure 6)"
    )
  )

  population <- list(
    species = "in vitro (Escherichia coli clinical isolates)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    organism = "Escherichia coli: 12 clinical isolates (5 producing extended-spectrum beta-lactamase) plus ATCC 25922 for susceptibility testing (Table S3; nemonoxacin CAMHB MIC 0.5-128 mg/L for the clinical isolates, below 0.03 mg/L for ATCC 25922), and NTU-001 and NTU-002 added to the high-exposure SHU arm",
    system = "One-compartment open dilution model after Satta et al. (Figure 1): a central flask holding the bacterial suspension (inoculated at 1e6 CFU/mL and preincubated 1 h at 35 C), a drug reservoir, a drug-free medium reservoir and a waste flask, driven by a Masterflex L/S peristaltic pump at the Table S2 flow rates; drug concentrations and viable counts sampled at 0, 0.5, 1, 2, 4, 6, 8 and 24 h",
    medium = "synthetic human urine (SHU) at pH 6.5, modified with 0.1% (vol/vol) casamino acids and 0.1% (vol/vol) yeast extract (Table S1 recipe)",
    experimental_arms = "low urine exposure (Figure 5b strains W049-080, W031-023, W035-079 and W042-007) and high urine exposure (Figure 5c strains W035-079, W042-007, NTU-001 and NTU-002) plus drug-free growth controls",
    disease_state = "in vitro model of uncomplicated urinary tract infection",
    dose_range = "simulated day-10 urine exposure after oral nemonoxacin 500 mg once daily in healthy Chinese volunteers: low exposure = Table 1 minimum and high exposure = Table 1 median interval urine concentrations",
    regions = "China (Huashan Hospital, Fudan University, Shanghai)",
    notes = "Qi 2026 fitted the sigmoid Emax model separately for each combination of medium (CAMHB, SHU), PK/PD index (Cmax/MIC, AUC0-24/MIC, %T>MIC and, in SHU only, the same three against the SHU MIC) and endpoint (Delta log10 CFU, AUBKC); Table 2 reports 17 fits (the AUC0-24/MICSHU x Delta log10 CFU fit is absent from both Table 2 and Figure 6), packaged as 17 files. AUC0-24/MIC was the index best correlated with effect (Delta log10 CFU R^2 0.951 in CAMHB, 0.982 in SHU). Parameters were fitted with naive pooling of all isolates; the number of experiments behind each fit is not stated."
  )

  ini({
    # Qi 2026 Table 2, SHU, %T > MIC, Delta log10 CFU; R^2 = 0.735.
    # Methods eq: E = E0 - Emax * X^H / (X^H + EC50^H); the Table 2 column
    # headed gamma is the Hill coefficient H.
    le0 <- log(2.70)
    label("Change in log10 CFU/mL over 0-24 h at zero exposure, the drug-free growth E0 (log10 CFU/mL)") # Qi 2026 Table 2, SHU, %T > MIC, Delta log10 CFU, E0 = 2.70

    lemax <- log(66.6)
    label("Maximum drug effect Emax, the span from drug-free growth to maximal kill (log10 CFU/mL)") # Qi 2026 Table 2, SHU, %T > MIC, Delta log10 CFU, Emax = 66.6

    lec50 <- log(129)
    label("%T>MIC producing half the maximum effect EC50 (percent)") # Qi 2026 Table 2, SHU, %T > MIC, Delta log10 CFU, EC50 = 1.29 on the 0-1 fraction scale, x 100 = 129 percent

    lhill <- log(9.95)
    label("Hill coefficient H of the index-effect curve (unitless)") # Qi 2026 Table 2, SHU, %T > MIC, Delta log10 CFU, gamma = 9.95

    # Starting density: an experimental design input, not an estimate.
    log10_cfu0 <- fixed(6)
    label("log10 starting bacterial density in the central flask (log10 CFU/mL)") # Qi 2026 Methods, In vitro UTIs model: inoculated at 1 x 10^6 CFU/mL

    # Table 2 reports R^2 only, no residual SD and no standard errors, so the
    # residual SD is held at zero for deterministic typical-value simulation.
    addSd <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (log10 CFU/mL)") # not reported in Qi 2026 (Table 2 gives R^2 only)
  })

  model({
    e0 <- exp(le0)
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # PK/PD index %T>MIC, supplied in percent (0-100).
    pkpd_index <- FTMIC_NEMONOXACIN

    # Qi 2026 Methods sigmoid Emax equation; the effect equals e0 at zero
    # exposure and approaches e0 - emax at saturating exposure.
    effect <- e0 - emax * pkpd_index^hill / (pkpd_index^hill + ec50^hill)

    # Delta log10 CFU is the change over 0-24 h. Spreading it uniformly over
    # the 24 h window makes log10(bact) move by exactly `effect` by 24 h:
    #   d(log10 N)/dt = effect / 24  =>  dN/dt = ln(10) * (effect / 24) * N
    d/dt(bact) <- log(10) * (effect / 24) * bact
    bact(0) <- 10^log10_cfu0

    # log10 CFU/mL with a 1 CFU/mL floor so the log stays finite.
    log_cfu <- log10(bact + 1)
    log_cfu ~ add(addSd)
  })
}
