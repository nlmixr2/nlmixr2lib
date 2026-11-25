Srivastava_2017_linezolid_hfs_adr <- function() {
  description <- paste(
    "In vitro (hollow fiber system model of tuberculosis). Quadratic",
    "exposure-response model of the linezolid-resistant Mycobacterium",
    "tuberculosis H37Rv subpopulation (log10 CFU/mL on agar with 3x MIC",
    "linezolid) against the percentage of time the concentration exceeds the",
    "MIC (%T>MIC), one quadratic per sampling day from day 14 to day 49",
    "(Srivastava 2017 Figure 2C, the 'system of changing U-shaped curves').",
    "The paper prints no coefficients; every coefficient was obtained by",
    "digitizing the drawn curves of Figure 2C and refitting the quadratic.",
    "Algebraic model with no ODE: model time is the study day and selects the",
    "day-specific curve; %T>MIC enters as a covariate. Predictions are floored",
    "at 0 log10 CFU/mL, where Figure 2C clips the curves. No residual error",
    "was reported, so addSd is fixed at 0. Siblings:",
    "Srivastava_2017_linezolid_hfs_aucmic (total burden) and",
    "Srivastava_2017_linezolid (clinical Monte Carlo pharmacokinetics)."
  )
  reference <- paste(
    "Srivastava S, Magombedze G, Koeuth T, Sherman C, Pasipanodya JG,",
    "Raj P, Wakeland E, Deshpande D, Gumbo T. Linezolid dose that maximizes",
    "sterilizing effect while minimizing toxicity and resistance emergence for",
    "tuberculosis. Antimicrob Agents Chemother. 2017;61(8):e00751-17.",
    "doi:10.1128/AAC.00751-17. PMCID: PMC5527615. Quadratic curves of",
    "Figure 2C (coefficients digitized by the nlmixr2lib maintainers)."
  )
  vignette <- "Srivastava_2017_linezolid"
  units <- list(
    time = "day (study day; selects the day-specific curve)",
    dosing = "percent (%T>MIC, supplied as a covariate)",
    concentration = "log10 CFU/mL (log_cfu_res, linezolid-resistant subpopulation)"
  )

  depends <- c("FTMIC_LINEZOLID")

  # log_cfu_res (linezolid-resistant subpopulation) is an algebraic
  # observation output of this paper.
  paper_specific_compartments <- c("log_cfu_res")

  covariateData <- list(
    FTMIC_LINEZOLID = list(
      description = "Percentage of the 24-h dosing period during which the linezolid concentration exceeds the MIC of the challenge strain (%T>MIC)",
      units = "percent (0-100)",
      type = "continuous",
      reference_category = NULL,
      notes = "Srivastava 2017 Figure 2C x-axis. The hollow-fiber arms took %T>MIC values of 0, about 29, 33, 54.5, 83 and 100 (the x positions of the Figure 2C data points). Hollow-fiber concentrations are free drug, against an MIC of 1 mg/L. Set to 0 for non-treated controls.",
      source_name = "%TMIC (Figure 2C)"
    )
  )

  population <- list(
    species = "in vitro (hollow fiber system model of tuberculosis, HFS-TB; Mycobacterium tuberculosis H37Rv)",
    n_subjects = 15L,
    n_studies = 1L,
    organism = "Mycobacterium tuberculosis H37Rv, semidormant at pH 5.8; linezolid MIC 1 mg/L; resistant subpopulation counted on Middlebrook agar with 3x MIC linezolid after up to 6 weeks of incubation",
    system = paste(
      "Fifteen HFS-TB units treated with linezolid once or twice daily to",
      "free AUC0-24/MIC values of 5 to 326 (see",
      "Srivastava_2017_linezolid_hfs_aucmic), peripheral compartment cultured",
      "weekly to day 56."
    ),
    disease_state = "in vitro model of tuberculosis sterilizing effect (semidormant, acidic-pH bacilli)",
    dose_range = "%T>MIC 0 to 100",
    notes = paste(
      "Figure 2C plots seven sampling days; the day-7 curve and a second",
      "flat curve both lie at about 0.1 log10 CFU/mL along the x-axis and",
      "cannot be told apart, so day 7 is not encoded. The resistant",
      "subpopulation was below the limit of detection on day 42 (Table 1",
      "footnote), although Figure 2C draws a day-42 curve; that curve is",
      "encoded as drawn."
    )
  )

  ini({
    # Every coefficient below is digitized (Figure 2C curve pixels, unweighted
    # quadratic refit; RMSE of the refit to the drawn curve <= 0.03 log10 CFU/mL).
    # Form for each day d: log_cfu_res = res0_dd + res1_dd * FTMIC + res2_dd * FTMIC^2.
    res0_d14 <- 6.0285
    label("Day-14 resistant burden at 0% T>MIC (log10 CFU/mL)") # Figure 2C day-14 curve, digitized (refit RMSE 0.007)
    res1_d14 <- 0.01015
    label("Day-14 linear coefficient (log10 CFU/mL per %)") # Figure 2C day-14 curve, digitized
    res2_d14 <- -0.00035631
    label("Day-14 quadratic coefficient (log10 CFU/mL per %^2)") # Figure 2C day-14 curve, digitized

    res0_d21 <- 5.5116
    label("Day-21 resistant burden at 0% T>MIC (log10 CFU/mL)") # Figure 2C day-21 curve, digitized (refit RMSE 0.011)
    res1_d21 <- 0.054815
    label("Day-21 linear coefficient (log10 CFU/mL per %)") # Figure 2C day-21 curve, digitized
    res2_d21 <- -0.0010528
    label("Day-21 quadratic coefficient (log10 CFU/mL per %^2)") # Figure 2C day-21 curve, digitized

    res0_d28 <- 6.4024
    label("Day-28 resistant burden at 0% T>MIC (log10 CFU/mL)") # Figure 2C day-28 curve, digitized (refit RMSE 0.026)
    res1_d28 <- 0.028077
    label("Day-28 linear coefficient (log10 CFU/mL per %)") # Figure 2C day-28 curve, digitized
    res2_d28 <- -0.00085407
    label("Day-28 quadratic coefficient (log10 CFU/mL per %^2)") # Figure 2C day-28 curve, digitized

    res0_d35 <- 5.8067
    label("Day-35 resistant burden at 0% T>MIC (log10 CFU/mL)") # Figure 2C day-35 curve, digitized (refit RMSE 0.006)
    res1_d35 <- 0.0073879
    label("Day-35 linear coefficient (log10 CFU/mL per %)") # Figure 2C day-35 curve, digitized
    res2_d35 <- -0.00060824
    label("Day-35 quadratic coefficient (log10 CFU/mL per %^2)") # Figure 2C day-35 curve, digitized

    res0_d42 <- 1.1852
    label("Day-42 resistant burden at 0% T>MIC (log10 CFU/mL)") # Figure 2C day-42 curve, digitized over 0-46% (refit RMSE 0.008)
    res1_d42 <- -0.033664
    label("Day-42 linear coefficient (log10 CFU/mL per %)") # Figure 2C day-42 curve, digitized
    res2_d42 <- 0.00023075
    label("Day-42 quadratic coefficient (log10 CFU/mL per %^2)") # Figure 2C day-42 curve, digitized

    res0_d49 <- 2.0459
    label("Day-49 resistant burden at 0% T>MIC (log10 CFU/mL)") # Figure 2C day-49 curve, digitized over 0-25% (refit RMSE 0.009)
    res1_d49 <- -0.072401
    label("Day-49 linear coefficient (log10 CFU/mL per %)") # Figure 2C day-49 curve, digitized
    res2_d49 <- 0.00053139
    label("Day-49 quadratic coefficient (log10 CFU/mL per %^2)") # Figure 2C day-49 curve, digitized

    # No residual error was reported (r^2 and AIC only), so the model is for typical-value simulation.
    addSd <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (not reported)") # Table 1 reports AIC only; no residual SD or standard errors
  })

  model({
    # Day windows centred on the weekly sampling days; times before day 17.5
    # use the day-14 curve and times from day 45.5 on use the day-49 curve.
    w14 <- 1 * (t < 17.5)
    w21 <- 1 * (t >= 17.5) * (t < 24.5)
    w28 <- 1 * (t >= 24.5) * (t < 31.5)
    w35 <- 1 * (t >= 31.5) * (t < 38.5)
    w42 <- 1 * (t >= 38.5) * (t < 45.5)
    w49 <- 1 * (t >= 45.5)

    res0 <- w14 * res0_d14 + w21 * res0_d21 + w28 * res0_d28 +
      w35 * res0_d35 + w42 * res0_d42 + w49 * res0_d49
    res1 <- w14 * res1_d14 + w21 * res1_d21 + w28 * res1_d28 +
      w35 * res1_d35 + w42 * res1_d42 + w49 * res1_d49
    res2 <- w14 * res2_d14 + w21 * res2_d21 + w28 * res2_d28 +
      w35 * res2_d35 + w42 * res2_d42 + w49 * res2_d49

    # Quadratic function of Gumbo et al. for the size of the resistant
    # subpopulation, floored at 0 log10 CFU/mL where Figure 2C clips it.
    res_raw <- res0 + res1 * FTMIC_LINEZOLID + res2 * FTMIC_LINEZOLID^2
    log_cfu_res <- max(res_raw, 0)
    log_cfu_res ~ add(addSd)
  })
}
