Busso_2017_exerciseTraining <- function() {
  description <- "Secondary-signal (indirect-response) model TI of the dose-response relationship between daily exercise training and physical performance in healthy men (Busso 2017). Each day's training dose (arbitrary training units, tu) enters a secondary signal pool that dissipates first-order; the signal stimulates production of performance (5-min maximal cycling power, W), and that production is inhibited in proportion to the same day's training dose. Performance is lost first-order towards its baseline. Typical values are the medians of the six individual fits in Table 2; no between-subject variability was estimated. The paper fitted a once-daily discrete recursion of the ODEs; this file reproduces that recursion exactly at whole days."
  reference <- paste(
    "Busso T. From an indirect response pharmacodynamic model towards a",
    "secondary signal model of dose-response relationship between exercise",
    "training and physical performance. Sci Rep. 2017;7:40422.",
    "doi:10.1038/srep40422.",
    "Experimental data: Busso T. Variable dose-response relationship between",
    "exercise training and performance. Med Sci Sports Exerc.",
    "2003;35(7):1188-1195. doi:10.1249/01.MSS.0000074465.13621.37."
  )
  vignette <- "Busso_2017_exerciseTraining"
  units <- list(
    time = "day",
    dosing = "tu (training units, arbitrary; one daily training dose per bolus into signal)",
    concentration = "W (performance, mean power of a 5-min maximal cycle-ergometer test)"
  )

  covariateData <- list()

  # `signal` is the latent secondary signal of the training stimulus, not a
  # drug-disposition compartment. `perf` (physical performance, the modelled
  # and observed endpoint) is a registered canonical.
  paper_specific_compartments <- "signal"

  compartmentData <- list(
    signal = list(
      analyte = "secondary signal of the training stimulus (latent)",
      units = "tu",
      specimen = "not applicable",
      verified = TRUE
    ),
    perf = list(
      analyte = "physical performance (Plim5', mean power of a 5-min maximal cycle-ergometer test)",
      units = "W",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 6L,
    n_studies = 1L,
    sex_female_pct = 0,
    disease_state = "Healthy male volunteers undertaking a 15-week cycle-ergometer training programme",
    dose_range = "Daily training doses of 0 to about 460 tu: 100 tu for each performance test; training sessions of 4 or 5 repetitions of 5 min at 85% of the last measured Plim5' (Methods, Equation 17)",
    regions = "France (Saint-Etienne)",
    notes = paste(
      "Methods 'Experimental data': data from Busso 2003 (the paper's",
      "reference 2). 8 weeks of 3 sessions/week, 1 week without training,",
      "4 weeks of 5 sessions/week, then 2 weeks without training; 40 to 46",
      "performance tests per subject. Each subject was fitted separately",
      "(maximum likelihood, R nlm); Table 2 gives the six individual",
      "estimates and the Results text their medians. The full data set is",
      "Supplementary Table S1."
    )
  )

  ini({
    # Perf0 is not estimated: the paper fixes each subject's baseline to their
    # first measured performance (Discretization section, before Equation 8),
    # and its simulations with the median parameters use 250 W.
    lrbase <- fixed(log(250)); label("Baseline performance Perf0 (W)") # Figure 6 caption 'baseline performance equal to 250 W'
    # Typical values are the medians of the six Table 2 estimates, the values
    # the paper simulates with (Figures 5 and 6). With six subjects each
    # median is the mean of the 3rd and 4th ranked values; the Results text
    # prints them rounded (0.0305, 0.0074, 0.292, 0.0021).
    lkout <- log(0.03055); label("First-order rate constant for loss of performance koff (1/day)") # Table 2 column 'koff', median of 0.0194, 0.0267, 0.0285, 0.0326, 0.0328, 0.0341 = (0.0285 + 0.0326) / 2; Results text 0.0305
    lkon_signal <- log(0.00739); label("Rate constant for transformation of the secondary signal into performance kon_s (1/day)") # Table 2 column 2 (kon_s), median of 0.00071, 0.00462, 0.00695, 0.00783, 0.00845, 0.03751 = (0.00695 + 0.00783) / 2; Results text 0.0074
    lkout_signal <- log(0.2915); label("First-order rate constant for dissipation of the secondary signal kout_s (1/day)") # Table 2 column 3 (kout_s), median of 0.053, 0.189, 0.286, 0.297, 0.405, 1.066 = (0.286 + 0.297) / 2; Results text 0.292
    lkinhib <- log(0.00206); label("Inhibition of performance production per unit of same-day training kin_i (1/tu)") # Table 2 column 4 (kin_i), median of 0.00148, 0.00199, 0.00201, 0.00211, 0.00213, 0.00246 = (0.00201 + 0.00211) / 2; Results text 0.0021

    # The paper maximises a Gaussian likelihood with the residual variance
    # profiled out (Methods, Equation 18), so sigma = sqrt(RSS / n) for each
    # subject and is not printed. This value is the median over the six
    # subjects of that sigma, recomputed by the maintainers from Model TI with
    # the Table 2 estimates on the Supplementary Table S1 data.
    addSd <- 6.02; label("Additive residual error on performance (W)") # derived: median of sqrt(RSS/n) over subjects 1-6 = 5.57, 6.09, 6.57, 5.94, 5.00, 7.17 W
  })

  model({
    rbase <- exp(lrbase)
    kout <- exp(lkout)
    kon_signal <- exp(lkon_signal)
    kout_signal <- exp(lkout_signal)
    kinhib <- exp(lkinhib)

    # Inhibition on a given day is proportional to that day's training dose
    # (Equations 5 and 12). Each daily dose is a bolus into `signal` at the
    # start of the day, so podo(signal) is today's dose for one day after it
    # and zero on a day without training.
    tad_signal <- tad(signal)
    podo_signal <- podo(signal)
    w_today <- 0
    if (!is.na(tad_signal)) {
      if (tad_signal < 1) {
        w_today <- podo_signal
      }
    }
    inhib <- kinhib * w_today

    # Signal dissipation plus inhibited transformation into performance
    # (Equation 13 as implemented in the authors' spreadsheet, which applies
    # day i's inhibition to the decay from day i to day i + 1).
    kloss_signal <- kout_signal + kon_signal * (1 - inhib)

    # The paper solves the model as a once-daily recursion (Equations 7-14):
    #   Perf(i + 1) = Perf(i) exp(-koff) + Prod(i)
    #   Prod(i) = kon0 + kon_s (1 - Inhib(i)) Signal(i), kon0 = Perf0 (1 - exp(-koff))
    # with Signal(i) the value just after day i's dose. Inside each day the
    # signal decays as exp(-kloss_signal t), so the production coefficient below
    # makes the within-day integral deliver exactly kon_s (1 - Inhib) Signal(i)
    # by the end of the day, and the zero-order term koff * Perf0 delivers
    # exactly kon0. Performance at whole days then equals the recursion.
    dk <- kloss_signal - kout
    fprod <- kon_signal * (1 - inhib) * exp(kout)
    if (abs(dk) > 1e-10) {
      fprod <- kon_signal * (1 - inhib) * dk / (exp(-kout) - exp(-kloss_signal))
    }

    d/dt(signal) <- -kloss_signal * signal
    perf(0) <- rbase
    d/dt(perf) <- kout * rbase + fprod * signal - kout * perf

    perf ~ add(addSd)
  })
}
