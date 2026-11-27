Srivastava_2017_linezolid_hfs_aucmic <- function() {
  description <- paste(
    "In vitro (hollow fiber system model of tuberculosis). Inhibitory",
    "sigmoid Emax exposure-response model of linezolid sterilizing effect",
    "against semidormant Mycobacterium tuberculosis H37Rv (pH 5.8) driven by",
    "the free-drug AUC0-24/MIC ratio, day-28 total bacterial burden",
    "(Srivastava 2017 Results; r^2 = 0.95, EC80 = AUC0-24/MIC 119). Algebraic",
    "model with no ODE and no time dependence: exposure enters as a",
    "covariate, and predictions refer to day 28 of once- or twice-daily",
    "linezolid. The hollow-fiber pharmacokinetics (one-compartment, elimination",
    "rate constant 0.17/h) has no reported volume and is not included. No",
    "residual error was reported, so addSd is fixed at 0. Siblings:",
    "Srivastava_2017_linezolid_hfs_adr (linezolid-resistant subpopulation) and",
    "Srivastava_2017_linezolid (the clinical Monte Carlo pharmacokinetic model",
    "used to translate the EC80 into a dose)."
  )
  reference <- paste(
    "Srivastava S, Magombedze G, Koeuth T, Sherman C, Pasipanodya JG,",
    "Raj P, Wakeland E, Deshpande D, Gumbo T. Linezolid dose that maximizes",
    "sterilizing effect while minimizing toxicity and resistance emergence for",
    "tuberculosis. Antimicrob Agents Chemother. 2017;61(8):e00751-17.",
    "doi:10.1128/AAC.00751-17. PMCID: PMC5527615. Day-28 equation and EC80",
    "printed in Results; MIC from Results."
  )
  vignette <- "Srivastava_2017_linezolid"
  units <- list(
    time = "day (day-28 readout; the model has no time dependence)",
    dosing = "mg*h/L (free linezolid AUC0-24, supplied as a covariate)",
    concentration = "log10 CFU/mL (log_cfu)"
  )

  depends <- c("AUC_LINEZOLID")

  covariateData <- list(
    AUC_LINEZOLID = list(
      description = "Free (unbound) linezolid area under the concentration-time curve over 0-24 h (AUC0-24), numerator of the AUC0-24/MIC PK/PD index",
      units = "mg*h/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Srivastava 2017 reports exposures only as AUC0-24/MIC in the hollow fiber system, which contains no protein, so the hollow-fiber AUC is the free AUC. The ratio is split into the absolute AUC0-24 and the model parameter mic = 1 mg/L, the MIC of the H37Rv strain, so the covariate and the index are numerically equal here. Arms: AUC0-24/MIC 5.08, 10.92, 16.24, 73.60 (once daily) and 5.25, 9.85, 20.55, 43.47, 73.60, 111.20, 157.33, 326.13 (twice daily) per Figure 2. Set to 0 for non-treated controls. A patient exposure must be converted to the free epithelial-lining-fluid AUC before use; Srivastava_2017_linezolid does this with the paper's ELF penetration ratio 3.29 and fraction unbound 0.7.",
      source_name = "AUC0-24/MIC (Results, Figure 2) x MIC"
    )
  )

  population <- list(
    species = "in vitro (hollow fiber system model of tuberculosis, HFS-TB; Mycobacterium tuberculosis H37Rv)",
    n_subjects = 15L,
    n_studies = 1L,
    organism = "Mycobacterium tuberculosis H37Rv, semidormant at pH 5.8; linezolid MIC 1 mg/L by Etest, broth macrodilution and MGIT (each done twice)",
    system = paste(
      "Fifteen HFS-TB units inoculated with semidormant M. tuberculosis at pH",
      "5.8 (5.41 log10 CFU/mL at the start of therapy). Linezolid was dosed",
      "into the central compartment to mimic the 4 h human half-life and the",
      "3.3-fold epithelial-lining-fluid AUC0-24 penetration, at free AUC0-24",
      "values of 5, 10, 20, 40 and 80 mg*h/L once or twice daily, and 111, 157",
      "and 326 mg*h/L twice daily, plus non-treated controls. The peripheral",
      "compartment was cultured on days 0, 3, 7 and then weekly to day 56."
    ),
    disease_state = "in vitro model of tuberculosis sterilizing effect (semidormant, acidic-pH bacilli)",
    dose_range = "free linezolid AUC0-24/MIC 0 to 326, once or twice daily",
    notes = paste(
      "Day-28 fit only: the paper prints the inhibitory sigmoid Emax equation",
      "for day 28, when AUC0-24/MIC was the index with the lowest AIC before",
      "the resistant subpopulation replaced most of the burden (Table 1).",
      "Analyses in ADAPT 5."
    )
  )

  ini({
    # Linezolid MIC of the challenge strain; a measured design constant, not estimated.
    mic <- fixed(1)
    label("Linezolid MIC of M. tuberculosis H37Rv (mg/L)") # Results, 'The linezolid MIC was 1 mg/liter, based on the 3 methods, on two different occasions each'

    # Results: 'effect (log10 CFU/ml) = (6.52 - 5.62) x AUC/MIC^3.41/(AUC/MIC^3.41 + 79.46^3.41), r2 = 0.95'
    # read as E = 6.52 - 5.62 x ... (see the vignette for why the bracketed form is a typesetting error)
    le0 <- log(6.52)
    label("Day-28 bacterial burden at zero exposure E0 (log10 CFU/mL)") # Results day-28 equation, E0 = 6.52
    lemax <- log(5.62)
    label("Maximal reduction in bacterial burden Emax (log10 CFU/mL)") # Results day-28 equation, Emax = 5.62
    lhill <- log(3.41)
    label("Hill coefficient H (unitless)") # Results day-28 equation, H = 3.41
    lec50 <- log(79.46)
    label("AUC0-24/MIC at half-maximal effect EC50 (h)") # Results day-28 equation, EC50 = 79.46 (EC80 = 79.46 x 4^(1/3.41) = 119, Results print 119)

    # No residual error was reported (r^2 and AIC only), so the model is for typical-value simulation.
    addSd <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (not reported)") # Results report r^2 = 0.95 only; no residual SD or standard errors
  })

  model({
    e0 <- exp(le0)
    emax <- exp(lemax)
    hill <- exp(lhill)
    ec50 <- exp(lec50)

    # PK/PD index: free AUC0-24/MIC (h).
    aucmic <- AUC_LINEZOLID / mic

    # Inhibitory sigmoid Emax model of the total bacterial burden; e0 - emax
    # is the burden at saturating exposure.
    log_cfu <- e0 - emax * aucmic^hill / (aucmic^hill + ec50^hill)
    log_cfu ~ add(addSd)
  })
}
