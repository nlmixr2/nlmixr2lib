Srivastava_2016_streptomycin_mrsa_ph50 <- function() {
  description <- paste(
    "In vitro (static test-tube study, methicillin-resistant Staphylococcus",
    "aureus ATCC 33591, pH 5.0). Inhibitory sigmoid Emax",
    "concentration-effect model of the 4-h MRSA burden against a static",
    "streptomycin concentration expressed as a multiple of the MIC",
    "(Srivastava 2016 Figure 2C). Srivastava 2016 used MRSA as a",
    "fast-growing surrogate to test whether aminoglycoside activity falls at",
    "acidic pH. Emax, the Hill coefficient and EC50 are from the Results;",
    "the paper does not print E0, which was digitised by the maintainers",
    "from the fitted curve in Figure 2C. Algebraic model with no ODE: the",
    "concentration enters as a covariate. No residual error was reported, so",
    "addSd is fixed at 0. Sibling: Srivastava_2016_streptomycin_mrsa_ph73."
  )
  reference <- paste(
    "Srivastava S, Modongo C, Siyambalapitiyage Dona CW, Pasipanodya JG,",
    "Deshpande D, Gumbo T. Amikacin optimal exposure targets in the hollow-fiber",
    "system model of tuberculosis. Antimicrob Agents Chemother.",
    "2016;60(10):5922-5927. doi:10.1128/AAC.00961-16. PMCID: PMC5038304.",
    "Emax, H and EC50 from Results (Figure 2C); E0 digitised from Figure 2C;",
    "MIC from Results."
  )
  vignette <- "Srivastava_2016_amikacin"
  units <- list(
    time = "h (4-h readout; the model has no time dependence)",
    dosing = "mg/L (static streptomycin concentration, supplied as a covariate)",
    concentration = "log10 CFU/mL (log_cfu)"
  )

  depends <- c("CONC_STR_MGL")

  covariateData <- list(
    CONC_STR_MGL = list(
      description = "Static in-vitro streptomycin concentration in the growth medium",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Srivastava 2016 Figure 2C plots concentration as a multiple of the MRSA MIC (0 to 1 x MIC); the model divides this absolute concentration by mic = 2 mg/L to form that ratio. Tested concentrations 0-2 mg/L. Set to 0 for drug-free controls.",
      source_name = "streptomycin concentration (Methods), 'Streptomycin concentration (xMIC)' (Figure 2C)"
    )
  )

  population <- list(
    species = "in vitro (static test-tube concentration-effect study; methicillin-resistant Staphylococcus aureus, ATCC 33591)",
    n_subjects = 48L,
    n_studies = 1L,
    organism = "methicillin-resistant Staphylococcus aureus (MRSA) ATCC 33591, log-phase growth; streptomycin MIC 2 mg/L",
    system = paste(
      "Log-phase MRSA co-incubated with streptomycin for 4 h at 35 C with rapid shaking in",
      "Mueller-Hinton broth acidified to pH 5.0 with citric acid, then washed twice in saline, serially diluted",
      "and plated overnight on Mueller-Hinton agar."
    ),
    disease_state = "in vitro concentration-effect study at neutral and acidic pH",
    dose_range = "streptomycin 0, 0.03, 0.06, 0.125, 0.25, 0.5, 1 and 2 mg/L, static",
    notes = "n_subjects counts the 8 concentrations x 3 replicates x 2 runs of the experiment. The paper used MRSA rather than M. tuberculosis for a faster readout when testing whether aminoglycoside activity falls at acidic pH (Figure 2)."
  )

  ini({
    # Streptomycin MIC of the challenge strain; a measured design constant, not estimated.
    mic <- fixed(2)
    label("Streptomycin MIC of MRSA ATCC 33591 (mg/L)") # Results, 'streptomycin MIC = 2 mg/liter'

    # Inhibitory sigmoid Emax model, Methods: 'The inhibitory sigmoid Emax model, where Emax is maximal kill, was employed'
    le0 <- log(6.99)
    label("4-h MRSA burden at zero concentration E0 (log10 CFU/mL)") # digitised by the maintainers from the fitted curve in Figure 2C (E0 not printed); plateau of the dashed curve (6.67) plus 0.97 x Emax
    lemax <- log(0.32)
    label("Maximal reduction in bacterial burden Emax (log10 CFU/mL)") # Results, Figure 2C, pH 5.0, 'the streptomycin Emax fell to 0.32' log10 CFU/mL
    lhill <- log(1.89)
    label("Hill coefficient H (unitless)") # Results, Figure 2C, H = 1.89
    lec50 <- log(0.03)
    label("Concentration/MIC at half-maximal effect EC50 (unitless)") # Results, Figure 2C, EC50 = 0.03 times the MIC

    # No residual error was reported (only r^2, plus AIC for the HFS fits), so the model is for
    # typical-value simulation.
    addSd <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (not reported)") # Figure 2C reports r^2 only; no residual SD or standard errors
  })

  model({
    e0 <- exp(le0)
    emax <- exp(lemax)
    hill <- exp(lhill)
    ec50 <- exp(lec50)

    # Concentration as a multiple of the MIC (unitless), the Figure 2 x-axis.
    concmic <- CONC_STR_MGL / mic

    # Inhibitory sigmoid Emax model of the bacterial burden; e0 - emax is
    # the burden at saturating exposure.
    log_cfu <- e0 - emax * concmic^hill / (concmic^hill + ec50^hill)
    log_cfu ~ add(addSd)
  })
}
