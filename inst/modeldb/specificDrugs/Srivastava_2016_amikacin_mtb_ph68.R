Srivastava_2016_amikacin_mtb_ph68 <- function() {
  description <- paste(
    "In vitro (static test-tube study, M. tuberculosis H37Ra, pH 6.8).",
    "Inhibitory sigmoid Emax concentration-effect model of the day-7 burden",
    "of log-phase-growth Mycobacterium tuberculosis against a static",
    "amikacin concentration expressed as a multiple of the MIC (Srivastava",
    "2016 Figure 1A). Emax, the Hill coefficient and EC50 are from the",
    "Results; the paper does not print E0, which was digitised by the",
    "maintainers from the fitted curve in Figure 1A. With E0 below Emax the",
    "fitted floor is slightly negative (-0.12); the figure clips it at 0.",
    "Algebraic model with no ODE: the concentration enters as a covariate.",
    "No residual error was reported, so addSd is fixed at 0. Sibling:",
    "Srivastava_2016_amikacin_mtb_ph58."
  )
  reference <- paste(
    "Srivastava S, Modongo C, Siyambalapitiyage Dona CW, Pasipanodya JG,",
    "Deshpande D, Gumbo T. Amikacin optimal exposure targets in the hollow-fiber",
    "system model of tuberculosis. Antimicrob Agents Chemother.",
    "2016;60(10):5922-5927. doi:10.1128/AAC.00961-16. PMCID: PMC5038304.",
    "Emax, H and EC50 from Results (Figure 1A); E0 digitised from Figure 1A;",
    "MIC from Results."
  )
  vignette <- "Srivastava_2016_amikacin"
  units <- list(
    time = "day (day-7 readout; the model has no time dependence)",
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
      notes = "Srivastava 2016 Figure 1 plots concentration as a multiple of the MIC (0.25 to 16 x MIC on a log2 axis); the model divides this absolute concentration by mic = 0.5 mg/L to form that ratio. Tested concentrations 0-8 mg/L. Set to 0 for drug-free controls.",
      source_name = "amikacin concentration (Methods), 'Amikacin (xMIC)' (Figure 1)"
    )
  )

  population <- list(
    species = "in vitro (static test-tube concentration-effect study; Mycobacterium tuberculosis H37Ra, ATCC 25177)",
    n_subjects = 24L,
    n_studies = 1L,
    organism = "Mycobacterium tuberculosis H37Ra (ATCC 25177); amikacin MIC 0.5 mg/L",
    system = paste(
      "Static test tubes: 9 mL of day-4 log-phase M. tuberculosis in Middlebrook",
      "7H9 broth with 10% OADC (measured pH 6.8) at 1.5 x 10^5 CFU/mL plus 1 mL",
      "of amikacin solution, in triplicate, incubated 7 days at 37 C under 5%",
      "CO2, then washed twice and plated on Middlebrook 7H10 agar."
    ),
    disease_state = "in vitro tuberculosis concentration-effect study",
    dose_range = "amikacin 0, 0.125, 0.25, 0.5, 1, 2, 4 and 8 mg/L (0 to 16 x MIC), static",
    notes = "n_subjects counts the 8 concentrations x triplicate tubes of the experiment. The paper fitted the inhibitory sigmoid Emax model to the day-7 bacterial burden against concentration expressed as a multiple of the MIC (Figure 1)."
  )

  ini({
    # Amikacin MIC of the challenge strain; a measured design constant, not estimated.
    mic <- fixed(0.5)
    label("Amikacin MIC of M. tuberculosis H37Ra (mg/L)") # Results, 'The amikacin MIC against M. tuberculosis was 0.5 mg/liter on two separate occasions, using two separate methods'

    # Inhibitory sigmoid Emax model, Methods: 'The inhibitory sigmoid Emax model, where Emax is maximal kill, was employed'
    le0 <- log(5.27)
    label("Day-7 bacterial burden at zero concentration E0 (log10 CFU/mL)") # digitised by the maintainers from the fitted curve in Figure 1A (E0 not printed); with Emax, H and EC50 left free the same digitisation returns 5.39, 2.97 and 0.99
    lemax <- log(5.39)
    label("Maximal reduction in bacterial burden Emax (log10 CFU/mL)") # Results, Figure 1A, Emax = 5.39 log10 CFU/mL
    lhill <- log(2.96)
    label("Hill coefficient H (unitless)") # Results, Figure 1A, H = 2.96
    lec50 <- log(0.99)
    label("Concentration/MIC at half-maximal effect EC50 (unitless)") # Results, Figure 1A, EC50 = 0.99 times the MIC

    # No residual error was reported (only r^2, plus AIC for the HFS fits), so the model is for
    # typical-value simulation.
    addSd <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (not reported)") # Figure 1A reports r^2 only; no residual SD or standard errors
  })

  model({
    e0 <- exp(le0)
    emax <- exp(lemax)
    hill <- exp(lhill)
    ec50 <- exp(lec50)

    # Concentration as a multiple of the MIC (unitless), the Figure 1 x-axis.
    concmic <- CONC_AMK_MGL / mic

    # Inhibitory sigmoid Emax model of the bacterial burden; e0 - emax is
    # the burden at saturating exposure.
    log_cfu <- e0 - emax * concmic^hill / (concmic^hill + ec50^hill)
    log_cfu ~ add(addSd)
  })
}
