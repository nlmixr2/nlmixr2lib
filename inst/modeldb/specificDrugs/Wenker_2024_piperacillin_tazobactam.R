Wenker_2024_piperacillin_tazobactam <- function() {
  description <- paste(
    "In vitro (hollow-fibre infection model, HFIM; six Gram-negative strains).",
    "Sigmoidal Emax PK/PD-index model for piperacillin/tazobactam in its clinical 8:1 ratio.",
    "Wenker 2024 ran a 21-experiment dose-fractionation study (7 growth controls plus 14 treated",
    "arms) in an HFIM and fit, by non-linear least squares on the 24 h timepoint, the change in",
    "bacterial density over 24 h as a function of each candidate PK/PD index. fT>MIC (the",
    "percentage of the dosing interval during which the free piperacillin/tazobactam",
    "concentration exceeds the isolate MIC) was the best-fitting index (AIC 94.1, R^2 0.691)",
    "versus fCmax/MIC (AIC 102, R^2 0.526) and fAUC/MIC (AIC 96.9, R^2 0.62).",
    "The packaged model is the fT>MIC model of Wenker 2024 Figure 1:",
    "effect = e0 - emax * FTMIC_TZP^hill / (ec50^hill + FTMIC_TZP^hill),",
    "where effect is the change in log10(CFU/mL) accrued over the 24 h experiment (positive =",
    "net growth, negative = net kill), e0 = 4.155 (FIXED) is the change in the untreated growth",
    "controls, emax = 8.851 (FIXED) is the maximum achievable reduction of that change,",
    "ec50 = 51.01 %fT>MIC is the index value giving half-maximal effect (the paper calls it",
    "EI50), and hill = 2.11 is the sigmoidicity coefficient (the paper calls it gam).",
    "There is NO pharmacokinetic component: exposure enters solely through the externally",
    "supplied per-experiment PK/PD index FTMIC_TZP, and the model has no ODE states because",
    "Wenker 2024 fit only the 24 h endpoint and never reported the starting inoculum numerically.",
    "The corresponding PK/PD targets recovered from these coefficients are 48% fT>MIC for",
    "bacteriostasis, 60% fT>MIC for 1 log10 kill and 75% fT>MIC for 2 log10 kill; the paper",
    "reports 48%, 60% and 77% respectively (see the vignette Errata for the 2 log10 discrepancy).",
    "Neither between-experiment variability nor a residual error magnitude was reported, so no",
    "eta parameters are present and addSd is FIXED at 0; the model is intended for typical-value",
    "and target-derivation simulation.",
    sep = " "
  )
  reference <- paste(
    "Wenker SAM, Alabdulkarim N, Readman JB, Slob EMA, Satta G, Ali S, Gadher N,",
    "Shulman R, Standing JF. (2024).",
    "Defining the pharmacokinetic/pharmacodynamic index of piperacillin/tazobactam within a",
    "hollow-fibre infection model to determine target attainment in intensive care patients.",
    "JAC-Antimicrobial Resistance 6(2):dlae036.",
    "doi:10.1093/jacamr/dlae036.",
    "Structural model and all four parameter estimates: main-text Figure 1 (panel annotation).",
    "HFIM system design (half-life, circulating volumes, dilution rates, dose levels):",
    "Supplementary data, Supplementary Methods.",
    sep = " "
  )
  vignette <- "Wenker_2024_piperacillin_tazobactam"
  units <- list(
    time = "hour",
    dosing = "% of the dosing interval above the MIC (PK/PD index supplied as a covariate)",
    concentration = "log10(CFU/mL) change over 24 h (observation)"
  )

  depends <- c("FTMIC_TZP")

  covariateData <- list(
    FTMIC_TZP = list(
      description        = "Percentage of the dosing interval during which the free piperacillin/tazobactam concentration exceeds the MIC of the challenge isolate",
      units              = "% (0-100)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-experiment PK/PD index, supplied as a covariate because Wenker 2024 published no",
        "structural PK model for the HFIM (Supplementary Methods states only that a model was",
        "fit in nlmixr2 to the bioassay concentrations, with the fit shown graphically in",
        "Figure S2 and no parameter values reported). Set to 0 for the untreated growth-control",
        "experiments, at which the sigmoid term vanishes and the predicted 24 h change reduces",
        "to e0. The scale is percent (0-100), NOT a fraction: the fitted ec50 of 51.01 is on the",
        "percent scale (Wenker 2024 Figure 1 x-axis, labelled 'Percent of time > MIC').",
        "The HFIM circulates protein-free Mueller-Hinton broth, so the free fraction is 1 and",
        "fT>MIC equals T>MIC in this system.",
        "Not previously in inst/references/covariate-columns.md; registered with this extraction",
        "as the first member of the FTMIC_<DRUG> family, parallel to the AUC_<DRUG> family",
        "(AUC_TILM and siblings) for interval-integrated exposure driving a PK/PD-index model.",
        sep = " "
      ),
      source_name        = "fT>MIC / 'Percent of time > MIC' (Wenker 2024 Figure 1 x-axis and Results, Dose-fractionation study; per-patient medians in Table 1)"
    )
  )

  population <- list(
    species            = "in vitro (hollow-fibre infection model; Escherichia coli, Klebsiella pneumoniae, Pseudomonas aeruginosa)",
    n_subjects         = 21L,
    n_studies          = 1L,
    organism           = paste(
      "Six Gram-negative strains, one laboratory reference strain and five clinical isolates",
      "(Supplementary Methods). Piperacillin/tazobactam MICs by CLSI broth microdilution",
      "(Results, Dose-fractionation study): Escherichia coli ATCC 25922 = 2 mg/L,",
      "Escherichia coli DWEC107 = 32 mg/L, Klebsiella pneumoniae DWKC01 = 256 mg/L,",
      "Klebsiella pneumoniae JRKC01 = 32 mg/L, Pseudomonas aeruginosa SWPC02 = 64 mg/L,",
      "Pseudomonas aeruginosa SWPC04 = 4 mg/L.",
      sep = " "
    ),
    system             = paste(
      "Hollow-fibre infection model (Supplementary Methods). Cartridges primed 24 h with PBS",
      "then 24 h with Mueller-Hinton broth. Pumps set to mimic a piperacillin half-life of 2 h.",
      "Cartridge volume 28 mL, tubing 39 mL, central reservoir 10 mL (run 1) or 30 mL (run 2),",
      "giving a total circulating volume of 77 mL (run 1) or 97 mL (run 2); dilution rate",
      "26.6 mL/h (run 1) or 33.6 mL/h (run 2). All experiments at 37 C.",
      sep = " "
    ),
    medium             = "Mueller-Hinton broth (protein-free, so free fraction = 1)",
    duration           = "24 h per experiment, with piperacillin/tazobactam infused every 6 h",
    regimens           = paste(
      "Dose levels were set from the predicted bolus Cmax of piperacillin/tazobactam",
      "(Supplementary Methods): high 128/16 mg/L, middle 32/8 mg/L, low 8/1 mg/L.",
      "Infusion duration 30 min or 3 h, every 6 h over 24 h.",
      "Piperacillin and tazobactam were used in an 8:1 ratio throughout.",
      sep = " "
    ),
    starting_inoculum  = "Not reported numerically; the inoculum titre was estimated by UV spectrophotometer OD before each experiment (Methods, Method for dose-fractionation assays) and the resulting time courses are shown graphically in Figure S1",
    n_experiments      = "21 total: 7 growth control and 14 dose-fractionation experiments (Methods, Method for dose-fractionation assays)",
    disease_state      = "Not applicable (in vitro infection model)",
    notes              = paste(
      "The clinical arm of Wenker 2024 (a nine-patient ICU case series, Table 1) is NOT part of",
      "this packaged model. That target-attainment analysis simulated plasma concentrations",
      "using the primary PK parameters of a separate publication (Lonsdale et al. 2020 ABDose,",
      "doi:10.1093/jac/dkaa363) plus a dialysis clearance parameter taken from an ertapenem",
      "model (Eyler et al. 2014, doi:10.1128/AAC.02090-12; packaged separately as",
      "Eyler_2014_ertapenem). Wenker 2024 reports no parameter values from either source, so no",
      "PK layer is extractable from this paper.",
      sep = " "
    )
  )

  ini({
    # ==================================================================
    # Wenker 2024 Figure 1 fT>MIC PK/PD-index model. All four parameter
    # estimates are printed in the Figure 1 panel annotation:
    #   E0 = 4.155 FIX, EI50 = 51.01, Emax = 8.851 FIX, gam = 2.11
    # (with AIC = 94.1 and R^2 = 0.68 also printed there; the Results
    # text gives R^2 = 0.691).
    #
    # The fitted function is the change in log10(CFU/mL) over the 24 h
    # experiment as a function of the index I = %fT>MIC:
    #   effect(I) = E0 - Emax * I^gam / (EI50^gam + I^gam)
    # so effect(0) = E0 = +4.155 (untreated growth) and the asymptote as
    # I -> Inf is E0 - Emax = -4.696 log10(CFU/mL).
    # ==================================================================

    le0 <- fixed(log(4.155))
    label("Change in log10(CFU/mL) over 24 h without drug, E0 (log10 CFU/mL)")
    # Wenker 2024 Figure 1 panel annotation: E0 = 4.155 FIX. Marked FIX on
    # the figure; the paper does not state how the value was derived, but
    # it is the y-intercept of the fitted curve and corresponds to the
    # seven untreated growth-control experiments (Methods).

    lemax <- fixed(log(8.851))
    label("Maximum drug-induced reduction of the 24 h change, Emax (log10 CFU/mL)")
    # Wenker 2024 Figure 1 panel annotation: Emax = 8.851 FIX. Marked FIX
    # on the figure. Subtractive: the fitted curve falls from E0 towards
    # the asymptote E0 - Emax = -4.696 log10(CFU/mL).

    lec50 <- log(51.01)
    label("Index value giving half-maximal effect, EI50 (% fT>MIC)")
    # Wenker 2024 Figure 1 panel annotation: EI50 = 51.01. Canonical
    # sigmoid-Emax half-maximal-driver parameter; the driver here is the
    # PK/PD index in percent, not a concentration, so the paper names it
    # EI50 rather than EC50.

    lhill <- log(2.11)
    label("Sigmoidicity (Hill) coefficient on the PK/PD index (unitless)")
    # Wenker 2024 Figure 1 panel annotation: gam = 2.11. Canonical PD Hill
    # exponent (parameter-names.md: lhill is the canonical for the
    # sigmoidicity coefficient of a sigmoidal Emax function, whatever the
    # source paper calls it).

    # ==================================================================
    # Residual error. Wenker 2024 reports only AIC and R^2 for the
    # non-linear least-squares fit and gives no residual standard
    # deviation, so addSd is FIXED at 0 rather than invented. See the
    # vignette Errata. Same convention as Wen_2016_enrofloxacin_* and
    # Chen_2023_tilmicosin.
    # ==================================================================
    addSd <- fixed(0)
    label("Additive residual SD on the 24 h change in log10(CFU/mL); not reported by the source")
    # Wenker 2024 Results, Dose-fractionation study: only AIC = 94.1 and
    # R^2 = 0.691 are reported for the fT>MIC model.
  })

  model({
    # ==================================================================
    # Back-transform from the log scale.
    # ==================================================================
    e0   <- exp(le0)
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # ==================================================================
    # Wenker 2024 Figure 1 sigmoidal Emax PK/PD-index model.
    #
    # effect is the change in log10(CFU/mL) accrued over the 24 h
    # experiment: positive = net growth, negative = net kill. The three
    # published PK/PD targets are the index values at which effect equals
    # 0 (bacteriostasis), -1 (1 log10 kill) and -2 (2 log10 kill).
    #
    # FTMIC_TZP is on the percent scale (0-100) to match the fitted
    # ec50 of 51.01; it is 0 for the untreated growth-control arms, at
    # which effect reduces to e0.
    # ==================================================================
    effect <- e0 - emax * FTMIC_TZP^hill / (ec50^hill + FTMIC_TZP^hill)

    effect ~ add(addSd)
  })
}
