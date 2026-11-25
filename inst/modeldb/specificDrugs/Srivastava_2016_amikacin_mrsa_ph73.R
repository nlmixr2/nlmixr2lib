Srivastava_2016_amikacin_mrsa_ph73 <- function() {
  description <- paste(
    "In vitro (static test-tube study, methicillin-resistant Staphylococcus",
    "aureus ATCC 33591, pH 7.3). Inhibitory sigmoid Emax",
    "concentration-effect model of the 4-h MRSA burden against a static",
    "amikacin concentration expressed as a multiple of the MIC (Srivastava",
    "2016 Figure 2A). Srivastava 2016 used MRSA as a fast-growing surrogate",
    "to test whether aminoglycoside activity falls at acidic pH. Emax, the",
    "Hill coefficient and EC50 are from the Results; the paper does not",
    "print E0, which was digitised by the maintainers from the fitted curve",
    "in Figure 2A. Algebraic model with no ODE: the concentration enters as",
    "a covariate. No residual error was reported, so addSd is fixed at 0.",
    "Sibling: Srivastava_2016_amikacin_mrsa_ph50."
  )
  reference <- paste(
    "Srivastava S, Modongo C, Siyambalapitiyage Dona CW, Pasipanodya JG,",
    "Deshpande D, Gumbo T. Amikacin optimal exposure targets in the hollow-fiber",
    "system model of tuberculosis. Antimicrob Agents Chemother.",
    "2016;60(10):5922-5927. doi:10.1128/AAC.00961-16. PMCID: PMC5038304.",
    "Emax, H and EC50 from Results (Figure 2A); E0 digitised from Figure 2A;",
    "MIC from Results."
  )
  vignette <- "Srivastava_2016_amikacin"
  units <- list(
    time = "h (4-h readout; the model has no time dependence)",
    dosing = "mg/L (static amikacin concentration, supplied as a covariate)",
    concentration = "log10 CFU/mL (log_cfu)"
  )

  depends <- c("CONC_AMK_MGL")

  covariateData <- list(
    CONC_AMK_MGL = list(
      description = "Static in-vitro amikacin concentration in the growth medium",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Srivastava 2016 Figure 2A plots concentration as a multiple of the MRSA MIC (0 to 8 x MIC); the model divides this absolute concentration by mic = 0.25 mg/L to form that ratio. Tested concentrations 0-2 mg/L. Set to 0 for drug-free controls.",
      source_name = "amikacin concentration (Methods), 'Amikacin concentration (xMIC)' (Figure 2A)"
    )
  )

  population <- list(
    species = "in vitro (static test-tube concentration-effect study; methicillin-resistant Staphylococcus aureus, ATCC 33591)",
    n_subjects = 48L,
    n_studies = 1L,
    organism = "methicillin-resistant Staphylococcus aureus (MRSA) ATCC 33591, log-phase growth; amikacin MIC 0.25 mg/L",
    system = paste(
      "Log-phase MRSA co-incubated with amikacin for 4 h at 35 C with rapid shaking in",
      "Mueller-Hinton broth at neutral pH (7.3), then washed twice in saline, serially diluted",
      "and plated overnight on Mueller-Hinton agar."
    ),
    disease_state = "in vitro concentration-effect study at neutral and acidic pH",
    dose_range = "amikacin 0, 0.03, 0.06, 0.125, 0.25, 0.5, 1 and 2 mg/L, static",
    notes = "n_subjects counts the 8 concentrations x 3 replicates x 2 runs of the experiment. The paper used MRSA rather than M. tuberculosis for a faster readout when testing whether aminoglycoside activity falls at acidic pH (Figure 2)."
  )

  ini({
    # Amikacin MIC of the challenge strain; a measured design constant, not estimated.
    mic <- fixed(0.25)
    label("Amikacin MIC of MRSA ATCC 33591 (mg/L)") # Results, 'The amikacin MIC against MRSA was 0.25 mg/liter'

    # Inhibitory sigmoid Emax model, Methods: 'The inhibitory sigmoid Emax model, where Emax is maximal kill, was employed'
    le0 <- log(6.07)
    label("4-h MRSA burden at zero concentration E0 (log10 CFU/mL)") # digitised by the maintainers from the fitted curve in Figure 2A (E0 not printed); least-squares E0 with the printed Emax, H and EC50 held fixed, rmse 0.02 log10
    lemax <- log(2.34)
    label("Maximal reduction in bacterial burden Emax (log10 CFU/mL)") # Results, Figure 2A, Emax = 2.34 log10 CFU/mL
    lhill <- log(0.94)
    label("Hill coefficient H (unitless)") # Results, Figure 2A, H = 0.94
    lec50 <- log(0.12)
    label("Concentration/MIC at half-maximal effect EC50 (unitless)") # Results, Figure 2A, EC50 = 0.12 times the MIC

    # No residual error was reported (only r^2, plus AIC for the HFS fits), so the model is for
    # typical-value simulation.
    addSd <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (not reported)") # Figure 2A reports r^2 only; no residual SD or standard errors
  })

  model({
    e0 <- exp(le0)
    emax <- exp(lemax)
    hill <- exp(lhill)
    ec50 <- exp(lec50)

    # Concentration as a multiple of the MIC (unitless), the Figure 2 x-axis.
    concmic <- CONC_AMK_MGL / mic

    # Inhibitory sigmoid Emax model of the bacterial burden; e0 - emax is
    # the burden at saturating exposure.
    log_cfu <- e0 - emax * concmic^hill / (concmic^hill + ec50^hill)
    log_cfu ~ add(addSd)
  })
}
