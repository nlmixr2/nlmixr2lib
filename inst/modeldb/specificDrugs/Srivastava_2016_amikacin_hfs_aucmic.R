Srivastava_2016_amikacin_hfs_aucmic <- function() {
  description <- paste(
    "In vitro (hollow fiber system model of tuberculosis). Inhibitory",
    "sigmoid Emax exposure-response model of amikacin against Mycobacterium",
    "tuberculosis driven by the AUC0-24/MIC ratio (Srivastava 2016 Figure 5;",
    "AIC -9.448, r^2 = 0.98, EC90 = AUC0-24/MIC 102.74). The paper ranked",
    "this index second to Cmax/MIC. Day-14 burden in the peripheral",
    "compartment of the hollow fiber system model of tuberculosis (HFS-TB),",
    "M. tuberculosis H37Ra (amikacin MIC 0.5 mg/L). Algebraic model with no",
    "ODE and no time dependence: exposure enters as a covariate, and",
    "predictions refer to day 14 of once-daily (or every-48-h / twice-daily)",
    "amikacin. The paper's two-compartment model of the HFS",
    "central-compartment concentrations has no reported parameters and is",
    "not included. No residual error was reported, so addSd is fixed at 0.",
    "Siblings for the other indices: Srivastava_2016_amikacin_hfs_cmaxmic,",
    "Srivastava_2016_amikacin_hfs_tmic,",
    "Srivastava_2016_amikacin_hfs_cminmic."
  )
  reference <- paste(
    "Srivastava S, Modongo C, Siyambalapitiyage Dona CW, Pasipanodya JG,",
    "Deshpande D, Gumbo T. Amikacin optimal exposure targets in the hollow-fiber",
    "system model of tuberculosis. Antimicrob Agents Chemother.",
    "2016;60(10):5922-5927. doi:10.1128/AAC.00961-16. PMCID: PMC5038304.",
    "Equation and estimates printed in Figure 5; MIC from Results."
  )
  vignette <- "Srivastava_2016_amikacin"
  units <- list(
    time = "day (day-14 readout; the model has no time dependence)",
    dosing = "mg*h/L (amikacin AUC0-24, supplied as a covariate)",
    concentration = "log10 CFU/mL (log_cfu)"
  )

  depends <- c("AUC_AMIKACIN")

  covariateData <- list(
    AUC_AMIKACIN = list(
      description = "Amikacin area under the concentration-time curve over 0-24 h (AUC0-24), numerator of the AUC0-24/MIC PK/PD index",
      units = "mg*h/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Srivastava 2016 calculated AUC0-24 from the fitted two-compartment model of the hollow-fiber central-compartment concentrations and reports exposures only as AUC0-24/MIC (Figure 5 x-axis, 0 to about 800). The ratio is split into the absolute AUC0-24 and the model parameter mic = 0.5 mg/L, the MIC of the H37Ra strain. For the every-48-h schedule the paper still used the AUC over 0-24 h. Set to 0 for non-treated controls.",
      source_name = "AUC0-24/MIC (Figure 5) x MIC"
    )
  )

  population <- list(
    species = "in vitro (hollow fiber system model of tuberculosis, HFS-TB; Mycobacterium tuberculosis H37Ra, ATCC 25177)",
    n_subjects = 12L,
    n_studies = 1L,
    organism = "Mycobacterium tuberculosis H37Ra (ATCC 25177), log-phase growth; amikacin MIC 0.5 mg/L by broth macrodilution and by Etest on Middlebrook 7H10 agar (each done twice); amikacin mutation frequency 2 x 10^-7 on 6x MIC agar",
    system = paste(
      "Hollow-fiber cartridges (FiberCell) with 20 mL of 10^6 CFU/mL log-phase",
      "M. tuberculosis in the peripheral compartment, in Middlebrook 7H9 broth",
      "with 10% OADC at pH 6.8. Amikacin was infused over 30 min into the",
      "central compartment to mimic the concentration-time profiles in the lungs",
      "of tuberculosis patients, assuming negligible protein binding. Exposures",
      "matched human doses of 0 to 30 mg/kg/day given once daily, six of them",
      "(1.25 to 30 mg/kg) also given every 48 h, and 30 mg/kg also split twice",
      "daily. Central-compartment drug concentrations were sampled 0.5, 2.5, 5,",
      "7.5, 10, 15 and 23.5 h after the first dose and fitted with a",
      "two-compartment model whose parameter estimates are not reported; the",
      "peripheral compartment was cultured on days 3, 7, 14, 21 and 28."
    ),
    disease_state = "in vitro model of pulmonary tuberculosis (log-phase-growth bacilli)",
    dose_range = "human-equivalent amikacin exposures of 0, 0.35, 0.75, 1.25, 2.5, 5, 10, 20 and 30 mg/kg/day; every-48-h and twice-daily schedules as above",
    notes = paste(
      "Day-14 exposure-response fits only: the paper fitted the inhibitory",
      "sigmoid Emax model on each sampling day and kept day 14, which had the",
      "lowest AIC and the highest r^2 for every PK/PD index. n_subjects counts",
      "the 12 day-14 points plotted in Figures 4-6 (Methods describes 16 HFS-TB",
      "units in total). No amikacin-resistant isolates arose in any unit, so",
      "the model has no resistant subpopulation. Analyses in ADAPT 5."
    )
  )

  ini({
    # Amikacin MIC of the challenge strain; a measured design constant, not estimated.
    mic <- fixed(0.5)
    label("Amikacin MIC of M. tuberculosis H37Ra (mg/L)") # Results, 'The amikacin MIC against M. tuberculosis was 0.5 mg/liter on two separate occasions, using two separate methods'

    # Figure 5: 'E (log10 CFU/mL)=6.89-6.68*(AUC0-24/MIC)^1.70/(AUC0-24/MIC^1.70+28.21^1.70); r2>0.98'
    le0 <- log(6.89)
    label("Day-14 bacterial burden at zero exposure E0 (log10 CFU/mL)") # Figure 5 equation, E0 = 6.89
    lemax <- log(6.68)
    label("Maximal reduction in bacterial burden Emax (log10 CFU/mL)") # Figure 5 equation, Emax = 6.68
    lhill <- log(1.7)
    label("Hill coefficient H (unitless)") # Figure 5 equation, H = 1.70
    lec50 <- log(28.21)
    label("AUC0-24/MIC at half-maximal effect EC50 (h)") # Figure 5 equation, EC50 = 28.21 (EC90 = 28.21 x 9^(1/1.70) = 102.7, Results print 102.74)

    # No residual error was reported (only r^2, plus AIC for the HFS fits), so the model is for
    # typical-value simulation.
    addSd <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (not reported)") # Figure 5 reports r^2 only; no residual SD or standard errors
  })

  model({
    e0 <- exp(le0)
    emax <- exp(lemax)
    hill <- exp(lhill)
    ec50 <- exp(lec50)

    # PK/PD index AUC0-24/MIC (h).
    aucmic <- AUC_AMIKACIN / mic

    # Inhibitory sigmoid Emax model of the bacterial burden; e0 - emax is
    # the burden at saturating exposure.
    log_cfu <- e0 - emax * aucmic^hill / (aucmic^hill + ec50^hill)
    log_cfu ~ add(addSd)
  })
}
