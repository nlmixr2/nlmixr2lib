Chen_2023_tilmicosin <- function() {
  description <- "Preclinical (piglet, Duroc x Landrace x Yorkshire crossbred). Sigmoidal Emax PK/PD-integration model for the antibacterial effect of orally administered tilmicosin against Pasteurella multocida serovar D:7 (strain C44-15, MIC = 0.25 ug/mL) in subcutaneously implanted tissue-cage fluid (TCF). Chen 2023 Section 2.6 parameterises the antibacterial effect over a 24 h dosing interval as E = E0 + (Emax - E0) * Ce^N / (EC50^N + Ce^N), where E is the MAGNITUDE of the log10(CFU/mL) reduction accrued over that interval (positive = kill), E0 is the corresponding change in the untreated control, and Ce is a PK/PD index. The packaged model uses the paper's best-correlating index AUC24h/MIC (R^2 = 0.92, versus 0.90 for Cmax/MIC and 0.83 for %T>MIC), formed as the per-interval covariate AUC_TILM divided by the parameter mic. Parameters from Chen 2023 Table 2: Emax = 1.09 log10(CFU/mL), E0 = 0.003 log10(CFU/mL), EC50 = 26.66 h, Hill N = 2.69. The bacterial density bact (linear CFU/mL) is integrated as d/dt(bact) = -ln(10) * (E / 24) * bact so that log10(bact) falls by exactly E across each 24 h interval, reproducing the paper's per-interval model exactly at every time the paper actually counted bacteria (24 h boundaries only). There is NO PK component: tilmicosin exposure enters as the externally supplied per-interval AUC_TILM, because Chen 2023 analysed the tissue-cage-fluid concentrations non-compartmentally in WinNonlin (Table 1 reports NCA AUC0-24h, Cmax and T>MIC per dose group per dosing day) and reported no structural PK model. Neither between-subject variability nor a residual error magnitude was reported, so no eta parameters are present and addSd is FIXED at 0; the model is intended for typical-value simulation."
  reference <- paste(
    "Chen Y, Ji X, Zhang S, Wang W, Zhang H, Ding H. (2023).",
    "Pharmacokinetic/pharmacodynamic integration of tilmicosin against",
    "Pasteurella multocida in a piglet tissue cage model.",
    "Frontiers in Veterinary Science 10:1260990.",
    "doi:10.3389/fvets.2023.1260990.",
    sep = " "
  )
  vignette <- "Chen_2023_tilmicosin"
  units <- list(
    time = "h",
    dosing = "ug*h/mL (tilmicosin AUC0-24h per dosing interval, supplied as a covariate)",
    concentration = "log10 CFU/mL (observation)"
  )

  depends <- c("AUC_TILM")
  paper_specific_compartments <- c("bact")

  covariateData <- list(
    AUC_TILM = list(
      description        = "Tilmicosin area under the tissue-cage-fluid concentration-time curve over the current 24 h dosing interval",
      units              = "ug*h/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-24-h-interval tilmicosin exposure in infected tissue-cage fluid, supplied as a piecewise-constant time-varying covariate (one value per 24 h interval). Chen 2023 obtained these by non-compartmental analysis in WinNonlin 6.1 (Section 2.5) and reports them in Table 1 for each dose group and each of the three dosing days, plus in Section 3.2 for the 24-48 h window after the final dose. Set to 0 for the untreated control group so the sigmoid term vanishes and the predicted change equals E0. Because the paper published no structural PK model, this covariate is the model's only route for drug exposure.",
      source_name        = "AUC0-24h (Chen 2023 Table 1; and Section 3.2 for the 24-48 h post-final-dose interval)"
    )
  )

  population <- list(
    species             = "pig (crossbred piglet, Duroc x Landrace x Yorkshire)",
    n_subjects          = 10L,
    n_studies           = 1L,
    weight_range        = "25-30 kg",
    sex_female_pct      = 50,
    organism            = "Pasteurella multocida strain C44-15, serovar D:7 (China Institute of Veterinary Drug Control, Beijing); tilmicosin MIC = 0.25 ug/mL in BOTH tissue-cage fluid and tryptic soy broth by CLSI (2013) microdilution",
    system              = "Subcutaneous tissue-cage infection model. Food-grade silicone tissue cages (65 mm long, 13 mm inner / 18 mm outer diameter, twelve holes at each end) implanted subcutaneously on both sides of the neck, equidistant between the jugular vein and the spinal cord, under pentobarbital-sodium general and procainamide-hydrochloride local anaesthesia; 4-week recovery before use",
    disease_state       = "Healthy piglets with an experimentally induced localised P. multocida infection: each sterile tissue cage inoculated with 1 mL of ~1.5 x 10^8 CFU/mL suspension, incubated 48 h, and only cages reaching ~10^7 CFU/mL entered the experiment",
    dose_range          = "30, 40, 50 and 60 mg/kg body weight tilmicosin (tilmicosin phosphate, 80.4% purity) orally, once daily for 3 days; plus an infected untreated control group given the same volume of normal saline",
    design              = "Ten piglets randomised to four treatment groups and one control group (two piglets and four tissue cages per group, one male and one female per group)",
    sampling            = "PK: ~0.5 mL tissue-cage fluid at 1, 3, 6, 9, 12 and 24 h after each administration, plus 48 and 72 h after the final dose. PD: ~0.1 mL for bacterial counting at 24 h after each administration, plus 48 and 72 h after the final dose",
    regions             = "China (South China Agricultural University, Guangzhou)",
    notes               = "Ethics approval 2018A001 (Committee on the Ethics of Animals of South China Agricultural University). The PK/PD index magnitudes reported in Chen 2023 Table 2 apply to a single 24 h dosing interval. The paper's headline total bacterial reductions after the three-day course (1.48, 2.82, 3.39 and 3.52 log10 CFU/mL at 30, 40, 50 and 60 mg/kg; Section 3.3) accumulate over FOUR consecutive 24 h intervals, i.e. through 96 h -- the three dosing intervals plus the 24-48 h window after the final dose (Discussion: 'the number of bacteria could be reduced by 3 log10CFU/mL within 96 h after administration'). See the vignette for the numerical reproduction."
  )

  ini({
    # =============================================================
    # Chen 2023 Table 2 -- sigmoidal Emax PK/PD-integration model
    # =============================================================
    # Chen 2023 Section 2.6 equation:
    #   E = E0 + (Emax - E0) * Ce^N / (EC50^N + Ce^N)
    # E is the change in bacterial count (log10 CFU/mL) over a 24 h
    # dosing interval, reported as a positive REDUCTION magnitude: the
    # printed form makes E rise from E0 at zero exposure to Emax at
    # infinite exposure, and both Table 2 values are positive. This sign
    # convention is confirmed numerically by the paper's own reported
    # thresholds -- substituting Table 2 into the equation returns
    # E = 0.335 at Ce = 19.64 h (the reported 1/3-log threshold) and
    # E = 0.749 at Ce = 35.64 h (the reported 3/4-log threshold).
    e0 <- 0.003
    label("Change in bacterial count in the untreated control over a 24 h interval E0 (log10 CFU/mL)")  # Chen 2023 Table 2, E0 = 0.003; kept on the natural scale because E0 is a near-zero net change that is not sign-constrained
    lemax <- log(1.09)
    label("Maximum antibacterial effect over a 24 h interval Emax (log10 CFU/mL reduction)")  # Chen 2023 Table 2, Emax = 1.09
    lec50 <- log(26.66)
    label("AUC24h/MIC producing 50% of the maximal effect EC50 (h)")  # Chen 2023 Table 2, EC50 = 26.66 h
    lhill <- log(2.69)
    label("Hill coefficient N defining the slope of the effect curve (unitless)")  # Chen 2023 Table 2, Slope (N) = 2.69

    # =============================================================
    # Minimum inhibitory concentration
    # =============================================================
    # The PK/PD index driving the sigmoid is AUC24h/MIC, so the MIC is
    # required numerically. FIXED because it is a measured property of
    # the challenge strain, not an estimated parameter. Change it to
    # apply the model to an isolate with a different susceptibility.
    mic <- fixed(0.25)
    label("Tilmicosin MIC against P. multocida D:7 (ug/mL;, measured)")  # Chen 2023 Section 3.1: MIC = 0.25 ug/mL in both TCF and TSB broth

    # =============================================================
    # Starting bacterial density
    # =============================================================
    # FIXED experimental design input, not an estimated parameter.
    log10_cfu0 <- fixed(7)
    label("log10 starting bacterial density in the tissue cage (log10 CFU/mL)")  # Chen 2023 Section 2.2 (cages with ~10^7 CFU/mL selected) and Section 3.3 (untreated control remained around 10^7 CFU/mL)

    # =============================================================
    # Residual error
    # =============================================================
    # Chen 2023 reported only the coefficient of determination of the
    # sigmoid fit (R^2 = 0.92 for AUC24h/MIC; Section 3.4) and gave no
    # residual standard deviation, so the density-scale residual SD is
    # held at zero for deterministic typical-value simulation. See the
    # vignette Assumptions and deviations section.
    addSd <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (0; not reported in Chen 2023)")  # Chen 2023 reported R^2 only, no residual SD
  })

  model({
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # PK/PD index: tilmicosin AUC over the current 24 h dosing interval
    # divided by the MIC (Chen 2023 Section 2.6). Units: h.
    aucmic <- AUC_TILM / mic

    # Chen 2023 Section 2.6 sigmoidal Emax equation. kill_log10 is the
    # MAGNITUDE of the log10(CFU/mL) reduction accrued over the current
    # 24 h dosing interval (positive = bacterial kill); it equals e0 at
    # zero exposure and approaches emax at saturating exposure.
    kill_log10 <- e0 + (emax - e0) * aucmic^hill / (ec50^hill + aucmic^hill)

    # Chen 2023 fitted only the per-interval reduction, and counted
    # bacteria only at 24 h interval boundaries. Spreading the interval
    # reduction uniformly across the interval makes log10(bact) fall by
    # exactly kill_log10 over each 24 h interval, so the trajectory
    # matches the paper's model at every counted time point:
    #   d(log10 N)/dt = -kill_log10 / 24
    # => d(N)/dt      = -ln(10) * (kill_log10 / 24) * N
    d/dt(bact) <- -log(10) * (kill_log10 / 24) * bact
    bact(0) <- 10^log10_cfu0

    # log10 CFU/mL observation with a 1-CFU/mL floor (matches the
    # Wen 2016 / Bulitta 2009 / Yadav 2017 in-vitro PD convention so the
    # log10 stays finite if bact is driven below 1 CFU/mL).
    Cc <- log10(bact + 1)
    Cc ~ add(addSd)
  })
}
