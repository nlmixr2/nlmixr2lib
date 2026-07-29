Marier_2002_tobramycin_rat_conventional <- function() {
  description <- paste(
    "Preclinical (rat).",
    "Two-compartment population PK model for the conventional",
    "(non-liposomal) formulation of tobramycin (Tobi inhalation solution,",
    "PathoGenesis) after a single 1,200 ug intratracheal dose to male",
    "Sprague-Dawley rats with chronic Burkholderia cepacia (strain BC",
    "1368) pulmonary infection. NONMEM ADVAN4 (depot, central,",
    "peripheral) parameterised in rate-constant form: first-order",
    "absorption ka into a lung central compartment carrying drug amount",
    "(not concentration -- volumes of distribution were not fitted because",
    "the dependent variable was the amount of tobramycin recovered from",
    "homogenised lung tissue, calculated as the measured tissue",
    "concentration times the lung volume per animal), inter-compartmental",
    "rate constants k12 and k21 between central and peripheral, first-",
    "order elimination kel from central, and a fitted lung",
    "bioavailability FL accounting for the fraction of the intratracheal",
    "dose actually reaching the lung tissue compartment. Comparator arm",
    "for Marier_2002_tobramycin_rat_liposomal; the conventional",
    "formulation shows faster absorption, faster elimination, and ~8-fold",
    "lower lung AUC than the liposomal formulation in the source paper",
    "(Table 1, Results).",
    sep = " ")
  reference <- paste(
    "Marier JF, Lavigne J, Ducharme MP.",
    "Pharmacokinetics and efficacies of liposomal and conventional",
    "formulations of tobramycin after intratracheal administration in rats",
    "with pulmonary Burkholderia cepacia infection.",
    "Antimicrob Agents Chemother. 2002;46(12):3776-3781.",
    "doi:10.1128/aac.46.12.3776-3781.2002.",
    sep = " ")
  vignette <- "Marier_2002_tobramycin_rat_conventional"
  units <- list(
    time          = "h",
    dosing        = "ug",
    concentration = "ug (amount in lung tissue; Vc not fitted)"
  )

  covariateData <- list()

  population <- list(
    species        = "rat (Sprague-Dawley, male)",
    n_subjects     = 36L,
    n_studies      = 1L,
    age_range      = "Adult (specific age not reported)",
    weight_range   = "175-225 g at study entry",
    sex_female_pct = 0,
    disease_state  = paste(
      "Chronic pulmonary Burkholderia cepacia infection (strain BC 1368;",
      "genomovar III, tobramycin MIC 128 ug/mL). Infection established 6",
      "days pre-dose by intratracheal inoculation of agar beads carrying",
      "10^6 CFU; confirmed by throat-swab culture 3-5 days post-inoculation."
    ),
    dose_range     = paste(
      "Single 1,200 ug intratracheal bolus of the conventional Tobi",
      "(PathoGenesis Canada Ltd.) tobramycin inhalation solution",
      "delivered as 100 uL via calibrated pipette followed by an air",
      "bolus to ensure complete intratracheal delivery."
    ),
    regions        = "Canada (MDS Pharma Services, St-Laurent, Montreal)",
    notes          = paste(
      "Sparse sampling: 3 rats per time point at 0.5, 1, 2, 3, 5, 7, 9, 11,",
      "13, 15, 16, 18, and 24 h post-dose (cardiac-puncture exsanguination",
      "under ketamine + xylazine anaesthesia). HPLC-UV quantitation of",
      "tobramycin in homogenised lung tissue. Two of the 39 dosed",
      "animals died during the surgical procedure (3 h post-dose) and",
      "the infection level could not be confirmed in a third animal,",
      "leaving n = 36 evaluable. See Marier 2002 Methods and Table 1",
      "(Conventional column) for the design and fitted model output."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters -- Marier 2002 Table 1 (Conventional column).
    # NONMEM v5 ADVAN4 in rate-constant form (see the liposomal sibling
    # file for the full provenance discussion). Paper's k_23 / k_32 /
    # k_20 translate to canonical k12 / k21 / kel.
    # ------------------------------------------------------------------
    lka     <- log(4.68)  ; label("Absorption rate constant ka into lung tissue (1/h)")                       # Marier 2002 Table 1: k_a = 4.68 /h (CV 38.7%)
    lfdepot <- log(0.77)  ; label("Lung bioavailability F_L (fraction of intratracheal dose reaching lung)")  # Marier 2002 Table 1: F_L = 0.77 (CV 17.2%)
    lk12    <- log(5.58)  ; label("Rate constant central -> peripheral (paper's k_23, 1/h)")                  # Marier 2002 Table 1: k_23 = 5.58 /h (CV 43.9%)
    lk21    <- log(0.330) ; label("Rate constant peripheral -> central (paper's k_32, 1/h)")                  # Marier 2002 Table 1: k_32 = 0.330 /h (CV 37.8%)
    lkel    <- log(1.10)  ; label("Elimination rate from lung tissue (paper's k_20, 1/h)")                    # Marier 2002 Table 1: k_20 = 1.10 /h (CV 19.8%)

    # ------------------------------------------------------------------
    # IIV (BSV in paper) -- omega^2 = log(CV^2 + 1) conversion of the
    # Marier 2002 Table 1 CV%. Table 1 reports diagonal CV% only.
    #   ka:    CV 38.7% -> log(1 + 0.387^2)  = 0.13958
    #   F_L:   CV 17.2% -> log(1 + 0.172^2)  = 0.02915
    #   k_23:  CV 43.9% -> log(1 + 0.439^2)  = 0.17626
    #   k_32:  CV 37.8% -> log(1 + 0.378^2)  = 0.13357
    #   k_20:  CV 19.8% -> log(1 + 0.198^2)  = 0.03845
    # ------------------------------------------------------------------
    etalka     ~ 0.13958  # Marier 2002 Table 1: CV(k_a)   = 38.7%
    etalfdepot ~ 0.02915  # Marier 2002 Table 1: CV(F_L)   = 17.2%
    etalk12    ~ 0.17626  # Marier 2002 Table 1: CV(k_23)  = 43.9%
    etalk21    ~ 0.13357  # Marier 2002 Table 1: CV(k_32)  = 37.8%
    etalkel    ~ 0.03845  # Marier 2002 Table 1: CV(k_20)  = 19.8%

    # ------------------------------------------------------------------
    # Residual error. Methods declare a combined proportional +
    # additive error model; Results report a single combined population
    # residual variability of 18.4% for the conventional arm. The
    # individual proportional and additive SDs are not separately
    # reported, so the combined residual is encoded as proportional only.
    # See vignette Assumptions and deviations.
    # ------------------------------------------------------------------
    propSd <- 0.184 ; label("Proportional residual variability on lung-tissue amount (fraction)")  # Marier 2002 Results paragraph after Table 1: 18.4% residual variability (conventional)
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual parameters (log-normal IIV per Methods).
    # ------------------------------------------------------------------
    ka     <- exp(lka     + etalka)
    fdepot <- exp(lfdepot + etalfdepot)
    k12    <- exp(lk12    + etalk12)
    k21    <- exp(lk21    + etalk21)
    kel    <- exp(lkel    + etalkel)

    # ------------------------------------------------------------------
    # 2. ODE system -- NONMEM ADVAN4 with central / peripheral
    # parameterised in amount units (ug). See liposomal sibling file
    # for the full rationale on the rate-constant-only parameterisation.
    # ------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                 k12 * central - k21 * peripheral1

    f(depot) <- fdepot

    # ------------------------------------------------------------------
    # 3. Observation. Vc was not fitted (see Methods 'Pharmacokinetic
    # analysis'); the central state holds the amount of tobramycin in
    # lung tissue (ug). Cc carries amount units, not concentration. See
    # vignette Assumptions and deviations for the amount-vs-concentration
    # mapping.
    # ------------------------------------------------------------------
    Cc <- central
    Cc ~ prop(propSd)
  })
}
