Dorey_2017_florfenicol_Pmultocida_serum <- function() {
  description <- paste0(
    "In vitro (porcine field isolates of P. multocida, serum). ",
    "Static-concentration time-kill sigmoid Emax PK/PD-index model for the antibacterial effect of ",
    "florfenicol against six EU pig-pneumonia field isolates of Pasteurella multocida grown in porcine serum (Invitrogen Gibco Porcine Serum). ",
    "Dorey 2017 Figure 2 parameterises the effect as E = E0 + Emax * Ce^N / (EC50^N + Ce^N), where E is ",
    "the SIGNED change in log10 CFU/mL between 0 and 24 h of incubation (positive = net growth, ",
    "negative = net kill), E0 is that change in the drug-free control and Ce is the PK/PD index ",
    "AUC24h/MIC. The coefficient multiplying the sigmoid is the Table 4 row 'Log Emax - Log E0', so the ",
    "model is written with canonical roles as E = e0 + (emax - e0) * Ce^N / (Ce^N + EC50^N), where ",
    "e0 is the Table 4 row 'Log E0' and emax is the row 'Log Emax' (the 24 h change at saturating ",
    "exposure). Parameters (means of six per-isolate fits, Table 4 column 'Serum'): e0 = 2.43, emax = -6.39 log10 CFU/mL, ",
    "EC50 = 32.3 h, N = 3.63. Solving the mean-parameter curve for E = 0, -3 and -4 returns AUC24h/MIC = ",
    "24.7, 36.8 and 42.4 h against the published means of the per-isolate breakpoints, 24.2, 37.3 and 43.9 h. ",
    "There is NO PK component: Dorey 2017 analysed the plasma florfenicol concentrations of 34 pigs ",
    "non-compartmentally (WinNonlin 6.4; Table 1) and published no structural PK model, so exposure ",
    "enters as the externally supplied covariate AUCMIC_FLOR, carried as the AUC24h/MIC ratio because ",
    "each isolate was dosed at multiples of its own MIC. The SDs in Table 4 describe variability across ",
    "the six isolates, not estimated variance components, so there are no eta parameters and addSd ",
    "is FIXED at 0; the model is intended for typical-value simulation. Sibling models from the same ",
    "paper: Dorey_2017_florfenicol_Pmultocida_broth, Dorey_2017_florfenicol_Apleuropneumoniae_broth, Dorey_2017_florfenicol_Apleuropneumoniae_serum."
  )
  reference <- paste(
    "Dorey L, Pelligand L, Cheng Z, Lees P. (2017).",
    "Pharmacokinetic/pharmacodynamic integration and modelling of florfenicol for the pig",
    "pneumonia pathogens Actinobacillus pleuropneumoniae and Pasteurella multocida.",
    "PLoS One 12(5):e0177568.",
    "doi:10.1371/journal.pone.0177568. PMCID: PMC5446118.",
    "Model equation: Figure 2. Parameter estimates: Table 4, Serum column.",
    sep = " "
  )
  vignette <- "Dorey_2017_florfenicol"

  units <- list(
    time = "h",
    dosing = "h (florfenicol AUC24h/MIC index, supplied as a covariate)",
    concentration = "log10 CFU/mL (observation: signed 24 h change in bacterial count)"
  )

  depends <- c("AUCMIC_FLOR")

  covariateData <- list(
    AUCMIC_FLOR = list(
      description = "Florfenicol PK/PD index: area under the concentration-time curve over 24 h divided by the MIC of the test isolate (AUC24h/MIC)",
      units = "h",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Dorey 2017 exposed each isolate to eight static florfenicol concentrations set as multiples ",
        "of that isolate's own MIC (0.25 to 8 x MIC) for 24 h, so the index equals 24 h times the MIC ",
        "multiple (6 to 192 h). It is carried as the ratio rather than as an absolute AUC divided by a ",
        "model `mic` parameter because the six isolates have different MICs and only their geometric ",
        "mean (0.30 ug/mL in serum, SD 0.16) is reported; the fitted EC50 and breakpoints are means of ",
        "per-isolate values each normalised by its own MIC. To drive the model from a plasma exposure, ",
        "form the ratio from the AUC0-24 and the MIC of the isolate of interest; Dorey 2017 corrected ",
        "broth MIC distributions by the serum:broth MIC ratio before using the serum breakpoints. Set ",
        "to 0 for drug-free control records so the sigmoid term vanishes and the predicted 24 h change ",
        "reduces to e0."
      ),
      source_name = "AUC24h/MIC (Dorey 2017 Figure 2 Ce; Table 4 column 'Serum')"
    )
  )

  population <- list(
    species = "in vitro (porcine field isolates of Pasteurella multocida)",
    n_subjects = 6L,
    n_studies = 1L,
    organism = "Six EU (UK, France, Germany) field isolates of Pasteurella multocida from cases of pig pneumonia, selected for logarithmic growth in both broth and pig serum, susceptibility to florfenicol by disk diffusion, and a spread of broth MICs (highest, lowest and four intermediate). Geometric mean florfenicol MIC 0.30 ug/mL (SD 0.16) in serum",
    system = "Static-concentration time-kill assay in porcine serum (Invitrogen Gibco Porcine Serum): eight florfenicol concentrations of 0.25 to 8 x the isolate MIC, 24 h incubation, each isolate tested in triplicate. Maximum count set to 1 x 10^10 CFU/mL; lower limit of quantification 33 CFU/mL. The sigmoid Emax model was fitted to AUC24h/MIC versus the 24 h change in log10 CFU/mL per isolate in WinNonlin, and the six per-isolate estimates averaged",
    disease_state = "in vitro; isolates from pig pneumonia",
    notes = paste0(
      "The plasma PK used for PK/PD integration came from 34 healthy Landrace Cross male pigs ",
      "(about 2 months old, 32.5-52.0 kg) given 15 mg/kg florfenicol intramuscularly as one of two ",
      "bioequivalent products (Nuflor or Norfenicol, n = 17 each), analysed non-compartmentally: ",
      "geometric mean Cmax 3.04 mg/L, AUC 53.7 h*mg/L, Cl/F 0.23 L/h/kg (Table 1). Dorey 2017 fitted ",
      "the sigmoid Emax model separately for each organism in broth and in serum and reports the four ",
      "parameter sets side by side in Tables 4 and 5, so all four are packaged."
    )
  )

  ini({
    # =================================================================
    # Dorey 2017 Table 4, Serum column -- sigmoid Emax PK/PD-index model
    # =================================================================
    # Figure 2 equation as printed:
    #   E = E0 + Emax * Ce^N / (EC50^N + Ce^N) x 100
    # The trailing 'x 100' is a typesetting error: with it the curve
    # would span hundreds of log10 units, while Figure 3 spans about
    # +2 to -6. The coefficient multiplying the sigmoid is the span
    # 'Log Emax - Log E0' (-8.83 here), not the 'Log Emax' row: only the
    # span brings the mean-parameter curve close to the published
    # breakpoints (see the vignette). The file uses canonical roles, e0 = effect at zero
    # exposure and emax = effect at saturating exposure, and writes
    #   E = e0 + (emax - e0) * Ce^N / (Ce^N + EC50^N)
    # e0 and emax are signed log10 differences and stay on the natural
    # scale; EC50 and N are positive and carried as logs. Values are
    # means of six per-isolate fits (SD across isolates in comments).
    e0 <- 2.43
    label("Change in log10 CFU/mL over 24 h in the drug-free control E0 (log10 CFU/mL)")  # Dorey 2017 Table 4, Serum, row 'Log E0 (CFU/mL)' = 2.43 (SD 0.50)
    emax <- -6.39
    label("Change in log10 CFU/mL over 24 h at saturating florfenicol exposure Emax (log10 CFU/mL; negative = kill)")  # Dorey 2017 Table 4, Serum, row 'Log Emax (CFU/mL)' = -6.39 (SD 1.51); row 'Log Emax - Log E0' = -8.83 (SD 1.00)
    lec50 <- log(32.3)
    label("Log AUC24h/MIC producing half the maximal effect EC50 (h)")  # Dorey 2017 Table 4, Serum, row 'AUC24h/MIC50 (h)' = 32.3 (SD 3.13)
    lhill <- log(3.63)
    label("Log Hill coefficient N, steepness of the AUC24h/MIC-effect curve (unitless)")  # Dorey 2017 Table 4, Serum, row 'Gamma' = 3.63 (SD 1.23)

    # =================================================================
    # Residual error
    # =================================================================
    # Dorey 2017 reports no residual SD for the sigmoid fits, so the
    # residual SD is held at zero for typical-value simulation.
    addSd <- fixed(0)
    label("Additive residual SD on the 24 h change in log10 CFU/mL (0; not reported in Dorey 2017)")  # Dorey 2017 Methods, PK/PD modelling: no residual error reported
  })

  model({
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # PK/PD index driving the sigmoid, supplied as an already-formed
    # ratio (see the covariateData notes).
    ce <- AUCMIC_FLOR

    # Dorey 2017 Figure 2 sigmoid Emax equation with canonical roles.
    # Cc is the SIGNED change in log10 CFU/mL between 0 and 24 h of
    # static-concentration incubation (positive = net growth, negative =
    # net kill); it equals e0 at zero exposure and approaches emax at
    # saturating exposure.
    Cc <- e0 + (emax - e0) * ce^hill / (ce^hill + ec50^hill)
    Cc ~ add(addSd)
  })
}
