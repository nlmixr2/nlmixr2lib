Beredaki_2024_amphotericinB_liposomal_cauris24h <- function() {
  description <- "In vitro (four Candida auris clinical isolates; one-compartment dilution PK/PD model). Liposomal amphotericin B (L-AMB) in vitro PK/PD model for the 24 h exposure-response relationship. Beredaki 2024 simulated human L-AMB pharmacokinetics in the internal compartment of a dilution model, re-spiking to a constant target peak once daily for 48 h with an average half-life of 10 h (range 5-12 h; target 9 h). The 24 h change in log10 CFU/mL from the starting inoculum was then related to the PK/PD index Cmax/MIC with the sigmoidal variable-slope Emax model E = Emax * EI^n / (EI50^n + EI^n) (Methods, 'PK/PD analysis'), where E is the REDUCTION in log10 CFU/mL relative to the drug-free control and the MIC is the CLSI M27 amphotericin B deoxycholate MIC. Beredaki 2024 printed the equation and the fitted curves (Figure 5, R^2 = 0.86) but no coefficients, so e0, lemax, lec50 and lhill were recovered by digitising the 24 h curve in Figure 5; the recovered coefficients place the stasis exposure at 5.02 Cmax/MIC against the 5 (3-7) reported in the Figure 5 legend. The fungal density bact (linear CFU/mL) is integrated as d/dt(bact) = ln(10) * ((e0 - kill24) / 24) * bact so that log10(bact) changes by exactly (e0 - kill24) across the paper's 24 h observation window, reproducing the endpoint model exactly at the time the paper actually fitted it. The companion 48 h parameterisation -- the one the paper carries into its Monte Carlo target-attainment analysis -- is Beredaki_2024_amphotericinB_liposomal_cauris48h, and the Candida albicans validation arm is Beredaki_2024_amphotericinB_liposomal_calbicans."
  reference <- paste(
    "Beredaki MI, Sanidopoulos I, Pournaras S, Meletiadis J. (2024).",
    "Defining optimal doses of liposomal amphotericin B against Candida auris:",
    "data from an in vitro pharmacokinetic/pharmacodynamic model.",
    "The Journal of Infectious Diseases 229(2):599-607.",
    "doi:10.1093/infdis/jiad583.",
    sep = " "
  )
  vignette <- "Beredaki_2024_amphotericinB_liposomal"
  units <- list(
    time = "h",
    dosing = "mg (L-AMB added to the internal compartment of the dilution model)",
    concentration = "mg/L (L-AMB, Cc); log10 CFU/mL (fungal density, log10cfu)"
  )

  # bact -- linear fungal density (CFU/mL) in the internal compartment.
  # Follows the in-vitro PD accumulator precedent in Chen_2023_tilmicosin.
  paper_specific_compartments <- c("bact")

  # specimen is "not applicable" for both states: the matrix is the RPMI-1640
  # growth medium in the internal compartment of the in vitro dilution model,
  # which has no entry in the biological-specimen vocabulary. The matrix is
  # recorded in the analyte strings and in population$system.
  compartmentData <- list(
    central = list(analyte = "liposomal amphotericin B in the internal compartment (RPMI-1640 growth medium)", units = "mg", specimen = "not applicable", verified = TRUE),
    bact    = list(analyte = "Candida auris viable count in the internal compartment (RPMI-1640 growth medium)", units = "CFU/mL", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species        = "in vitro (four Candida auris clinical isolates; one-compartment dilution PK/PD model)",
    n_subjects     = 4L,
    n_studies      = 1L,
    organism       = paste(
      "Four Candida auris isolates (Beredaki 2024 Table 1, median MIC in mg/L,",
      "CLSI M27-A3, amphotericin B deoxycholate / liposomal amphotericin B):",
      "C. auris 51 clade I 1 / 1; C. auris 52 clade I 2 (1-2) / 4;",
      "C. auris 55 clade I 0.5 / 0.5; C. auris 60 clade II 0.5 / 0.125 (0.06-0.125).",
      "Isolates provided by J. Meis, Canisius Wilhelmina Hospital, Nijmegen.",
      "The exposure index is formed on the CLSI amphotericin B deoxycholate MIC",
      "(Figure 5 x axis 'Cmax L-AMB / CLSI MIC'; Figure 4 panel titles quote",
      "'CLSI AMB/L-AMB' pairs and the Monte Carlo analysis is indexed on",
      "'CLSI AMB MICs')."
    ),
    system         = paste(
      "Previously validated one-compartment in vitro PK/PD dilution model:",
      "a 250 mL conical glass culture vessel (internal compartment) holding an",
      "initial 5 mL of fresh RPMI-1640, connected to a peristaltic pump",
      "(Minipuls Evolution, Gilson) that adds fresh medium at a rate matching",
      "the clearance of L-AMB in human plasma. The vessel is covered with",
      "aluminium foil to minimise light exposure and held at 37 degrees C on a",
      "magnetic stirrer; its volume rises to approximately 100 mL by 48 h.",
      "Growth medium was RPMI 1640 with L-glutamine, without bicarbonate,",
      "buffered to pH 7.0 with 0.165 M MOPS and supplemented with 100 mg/L",
      "chloramphenicol. L-AMB was AmBisome (Gilead Sciences) reconstituted to",
      "4 mg/mL."
    ),
    disease_state  = "In vitro infection of the internal compartment at a starting inoculum of 10^4 CFU/mL",
    dose_range     = paste(
      "L-AMB target peak concentrations Cmax 0.25-64 mg/L, added once daily for",
      "48 h. Target average half-life 9 h (the longer second distribution",
      "half-life of human L-AMB, which covers most of the 24 h dosing",
      "interval); achieved in vitro half-life 10 h (range 5-12 h) with peak",
      "concentrations within 10% of target."
    ),
    sampling       = paste(
      "PK: samples from the internal compartment assayed by microbiological",
      "diffusion against a Paecilomyces variotii strain; lowest limit of",
      "detection 0.03 mg/L for the 80% partial growth inhibition zone, linear",
      "over 0.03-1 mg/L (Figure 3A, r^2 = 0.9373), so samples expected above",
      "1 mg/L were pre-diluted into the linear range.",
      "PD: 500 uL sampled at regular intervals up to 48 h, 10-fold serially",
      "diluted in normal saline, 20 uL subcultured on Sabouraud dextrose agar,",
      "incubated at 30 degrees C for 24 h; dilutions yielding 10-50 colonies",
      "were counted (Figure 4 shows 0, 3, 6, 24 and 48 h)."
    ),
    regions        = "Greece (Clinical Microbiology Laboratory, Attikon University Hospital, Athens)",
    notes          = paste(
      "All experiments were repeated twice. Each drug-free control grew by",
      "more than 2.5 log10 CFU/mL, from 4.04 +/- 0.24 log10 CFU/mL at t = 0 h",
      "to 7.52 +/- 0.56 log10 CFU/mL at t = 48 h across all isolates. Because",
      "the internal-compartment volume rises about 20-fold by 48 h, repeating",
      "the analysis on absolute log10 CFU (CFU/mL multiplied by the 48 h",
      "volume) shifted the time-kill curves upward by roughly 1 log10 CFU but",
      "gave similar PK/PD targets. No carryover effect was detected at L-AMB",
      "concentrations at or above 16 times the MIC, and liposomes without",
      "amphotericin B had no antifungal activity. Data were analysed in",
      "GraphPad Prism 5.0."
    )
  )

  ini({
    # ================================================================
    # In vitro pharmacokinetics
    # Beredaki 2024 Methods, "Pharmacokinetic analysis": L-AMB exposures
    # with peak concentrations Cmax 0.25-64 mg/L and an average half-life
    # of 9 hours were simulated, with drug added at the corresponding
    # Cmax values once daily for 48 hours. Figure 3B shows the target
    # profile returning to the SAME peak at t = 24 h (target 64 mg/L at
    # both t = 0 and t = 24 h), i.e. the compartment is re-spiked to a
    # constant peak rather than dosed additively.
    # ================================================================
    # Beredaki 2024 Results, "Pharmacokinetics" (C. auris): "The
    # pharmacokinetic parameters of L-AMB were well simulated in the in
    # vitro model with an average half-life of 10 (range, 5-12) hours".
    # The ACHIEVED half-life is used, matching the system that generated
    # the pharmacodynamic data; the 9 h design target is one third of a
    # half-life away and is reproduced in the vignette by overriding lkel.
    lkel <- log(log(2) / 10)
    label("Log first-order rate of L-AMB removal from the internal compartment (1/h)")  # Beredaki 2024 Results: achieved in vitro t1/2 = 10 h (range 5-12 h); Methods design target 9 h

    # Volume of the internal compartment at the time the peak occurs.
    # FIXED: a measured dimension of the apparatus, not an estimated
    # parameter. It converts the L-AMB amount added to the culture vessel
    # into a concentration. Beredaki 2024 reports the compartment volume
    # rising to about 100 mL by 48 h, and abstracts the resulting dilution
    # into the reported mono-exponential half-life above, so a constant
    # volume anchored at the initial 5 mL reproduces the paper's own
    # concentration-time description exactly. See the vignette
    # "Assumptions and deviations".
    lvc <- fixed(log(0.005))
    label("Log volume of the internal compartment at dosing (L; apparatus dimension)")  # Beredaki 2024 Methods, "In vitro PK/PD model": culture vessel "containing fresh RPMI-1640 medium to an initial volume of 5 mL" = 0.005 L

    # Target peak L-AMB concentration reached at the start of every 24 h
    # dosing interval, and the numerator of the paper's PK/PD index.
    # FIXED experimental design input; change it to simulate another arm.
    # The default is the mean clinical peak of the standard 3 mg/kg q24h
    # i.v. dose used in the paper's Monte Carlo analysis.
    cmax <- fixed(21.87)
    label("Target peak L-AMB concentration per 24 h interval (mg/L; experimental design input)")  # Beredaki 2024 Methods, "Prediction of PK/PD target attainment": 3 mg/kg q24h i.v. gives L-AMB Cmax 21.87 +/- 12.47 mg/L; in vitro arms spanned Cmax 0.25-64 mg/L

    # ================================================================
    # Sigmoidal variable-slope Emax exposure-response, 24 h endpoint
    # Beredaki 2024 Methods, "PK/PD analysis":
    #   E = Emax * EI^n / (EI50^n + EI^n)
    # where EI is the exposure index Cmax/MIC, Emax is the maximum growth
    # rate, EI50 is the index giving 50% of Emax and n is the Hill
    # coefficient. Read with E as the REDUCTION relative to the drug-free
    # control, the quantity plotted in Figure 5 is (e0 - E) -- the printed
    # form alone gives E = 0 at zero exposure whereas the 24 h curve of
    # Figure 5 starts near +2.4 log10 CFU/mL.
    #
    # PROVENANCE: Beredaki 2024 printed this equation and the fitted
    # curves (Figure 5, R^2 = 0.86) but no coefficients. The four values
    # below were recovered by digitising the 24 h curve of Figure 5 from
    # the 600 dpi page render (1913 traced points spanning Cmax/MIC
    # 0.011-136) and least-squares fitting the four-parameter
    # log-logistic form; residual RMSE was 0.023 log10 CFU/mL, i.e. the
    # recovered coefficients reproduce the drawn curve to within its line
    # width. See the vignette "Source trace" and "Assumptions and
    # deviations".
    # ================================================================
    # Top asymptote: the 24 h change in log10 CFU/mL at zero exposure.
    # Kept on the natural scale because it is a signed net change.
    # Beredaki 2024 reports no measured 24 h drug-free control growth, so
    # unlike the 48 h fit this asymptote has no direct cross-check
    # against a quoted control value.
    e0 <- 2.417
    label("24 h change in log10 CFU/mL at zero exposure (fitted top asymptote)")  # digitised from Beredaki 2024 Figure 5, 24 h curve
    lemax <- log(5.491)
    label("Log maximum reduction Emax relative to the drug-free control (log10 CFU/mL)")  # digitised from Beredaki 2024 Figure 5, 24 h curve; implies a bottom asymptote of 2.417 - 5.491 = -3.07 log10 CFU/mL
    lec50 <- log(6.006)
    label("Log Cmax/MIC producing 50% of Emax (EI50, unitless)")  # digitised from Beredaki 2024 Figure 5, 24 h curve
    lhill <- log(1.339)
    label("Log Hill coefficient n of the exposure-response relationship (unitless)")  # digitised from Beredaki 2024 Figure 5, 24 h curve

    # ================================================================
    # Minimum inhibitory concentration
    # ================================================================
    # The PK/PD index driving the sigmoid is Cmax/MIC, so the MIC is
    # required numerically. FIXED because it is a measured property of the
    # isolate, not an estimated parameter. Beredaki 2024 indexes on the
    # CLSI M27 amphotericin B DEOXYCHOLATE MIC, not the L-AMB MIC.
    # Default is C. auris 51; change it to apply the model to another
    # isolate or to an arbitrary MIC from the CLSI distribution.
    mic <- fixed(1)
    label("Amphotericin B CLSI M27 MIC of the simulated isolate (mg/L; measured property of the isolate)")  # Beredaki 2024 Table 1: C. auris 51 (clade I) CLSI amphotericin B deoxycholate median MIC 1 mg/L (range 1)

    # ================================================================
    # Starting fungal density
    # ================================================================
    # FIXED experimental design input: the measured mean density at
    # t = 0 h, so the simulated trajectory starts where Figure 4 does.
    # The nominal inoculum was 10^4 CFU/mL.
    log10_cfu0 <- fixed(4.04)
    label("Log10 starting fungal density in the internal compartment (log10 CFU/mL)")  # Beredaki 2024 Results, "Pharmacodynamics": "from 4.04 +/- 0.24 log10 CFU/mL at t = 0 hours to 7.52 +/- 0.56 log10 CFU/mL at t = 48 hours for all isolates"; nominal inoculum 10^4 CFU/mL

    # ================================================================
    # Residual error
    # ================================================================
    # L-AMB concentration: Beredaki 2024 reports only that the achieved
    # peaks were "within 10% of target values" (Results,
    # "Pharmacokinetics") and gives no assay CV point estimate, so the
    # proportional residual SD is FIXED at that reported bound.
    propSd <- fixed(0.10)
    label("Proportional residual error on L-AMB concentration (fraction; the reported bound on peak accuracy)")  # Beredaki 2024 Results: "L-AMB Cmax concentrations within 10% of target values"

    # Fungal density: Beredaki 2024 reported only the coefficient of
    # determination of the Emax fit (R^2 = 0.86) and no residual standard
    # deviation, so the log10-density residual SD is held at zero for
    # deterministic typical-value simulation. See the vignette
    # "Assumptions and deviations".
    addSd_log10cfu <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (0; not reported in Beredaki 2024)")  # Beredaki 2024 Figure 5 reports R^2 = 0.86 only, no residual SD
  })

  model({
    kel <- exp(lkel)
    vc <- exp(lvc)
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # ---------------------------------------------------------------
    # In vitro pharmacokinetics (Beredaki 2024 Methods, "Pharmacokinetic
    # analysis"). L-AMB is added to the internal compartment once daily
    # and decays mono-exponentially. Dose the `central` compartment in
    # mg; the vignette computes the amount that re-spikes the compartment
    # to `cmax` at the start of every interval, matching the constant
    # target peaks in Figure 3B.
    # ---------------------------------------------------------------
    d/dt(central) <- -kel * central
    Cc <- central / vc

    # ---------------------------------------------------------------
    # PK/PD index and exposure-response (Beredaki 2024 Methods, "PK/PD
    # analysis"). Cmax is a design input rather than a state because the
    # paper fitted a 24 h ENDPOINT against the per-interval target peak,
    # so the index must be available from t = 0 for the trajectory below
    # to land on the fitted endpoint. kill24 is the paper's E: the
    # reduction in log10 CFU/mL over the 24 h observation window relative
    # to the drug-free control, rising from 0 at zero exposure towards
    # emax at saturating exposure.
    # ---------------------------------------------------------------
    ei <- cmax / mic
    kill24 <- emax * ei^hill / (ei^hill + ec50^hill)

    # Beredaki 2024 fitted the 24 h change in log10 CFU/mL for this curve.
    # Spreading that change uniformly over the window makes log10(bact)
    # move by exactly (e0 - kill24) between t = 0 and t = 24 h, so the
    # trajectory matches the fitted endpoint model at the time the paper
    # fitted it:
    #   d(log10 N)/dt = (e0 - kill24) / 24
    #   => d(N)/dt    = ln(10) * ((e0 - kill24) / 24) * N
    d/dt(bact) <- log(10) * ((e0 - kill24) / 24) * bact
    bact(0) <- 10^log10_cfu0

    # log10 CFU/mL with a 1-CFU/mL floor so the log stays finite if bact
    # is driven below 1 CFU/mL (matches the in-vitro PD convention used
    # by Chen_2023_tilmicosin).
    log10cfu <- log10(bact + 1)

    Cc ~ prop(propSd)
    log10cfu ~ add(addSd_log10cfu)
  })
}
