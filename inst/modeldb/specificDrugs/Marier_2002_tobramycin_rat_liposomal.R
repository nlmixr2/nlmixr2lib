Marier_2002_tobramycin_rat_liposomal <- function() {
  description <- paste(
    "Preclinical (rat).",
    "Two-compartment population PK model for the liposomal formulation of",
    "tobramycin (DPPC:DMPG 10:1 phospholipids, 230-400 nm extruded) after a",
    "single 1,200 ug intratracheal dose to male Sprague-Dawley rats with",
    "chronic Burkholderia cepacia (strain BC 1368) pulmonary infection.",
    "NONMEM ADVAN4 (depot, central, peripheral) parameterised in rate-",
    "constant form: first-order absorption ka into a lung central",
    "compartment carrying drug amount (not concentration -- volumes of",
    "distribution were not fitted because the dependent variable was the",
    "amount of tobramycin recovered from homogenised lung tissue,",
    "calculated as the measured tissue concentration times the lung",
    "volume per animal), inter-compartmental rate constants k12 and k21",
    "between central and peripheral, first-order elimination kel from",
    "central, and a fitted lung bioavailability FL accounting for the",
    "fraction of the intratracheal dose actually reaching the lung tissue",
    "compartment.",
    sep = " ")
  reference <- paste(
    "Marier JF, Lavigne J, Ducharme MP.",
    "Pharmacokinetics and efficacies of liposomal and conventional",
    "formulations of tobramycin after intratracheal administration in rats",
    "with pulmonary Burkholderia cepacia infection.",
    "Antimicrob Agents Chemother. 2002;46(12):3776-3781.",
    "doi:10.1128/aac.46.12.3776-3781.2002.",
    sep = " ")
  vignette <- "Marier_2002_tobramycin_rat_liposomal"
  units <- list(
    time          = "h",
    dosing        = "ug",
    concentration = "ug (amount in lung tissue; Vc not fitted)"
  )

  covariateData <- list()

  population <- list(
    species        = "rat (Sprague-Dawley, male)",
    n_subjects     = 39L,
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
      "Single 1,200 ug intratracheal bolus of liposomal tobramycin (Tobi",
      "rehydrated into DPPC:DMPG 10:1 phospholipid lyophilisate, extruded",
      "to 230-400 nm) delivered as 100 uL via calibrated pipette followed",
      "by an air bolus to ensure complete intratracheal delivery."
    ),
    regions        = "Canada (MDS Pharma Services, St-Laurent, Montreal)",
    notes          = paste(
      "Sparse sampling: 3 rats per time point at 0.5, 1, 2, 3, 5, 7, 9, 11,",
      "13, 15, 16, 18, and 24 h post-dose (cardiac-puncture exsanguination",
      "under ketamine + xylazine anaesthesia). HPLC-UV quantitation of",
      "tobramycin in homogenised lung tissue (Symmetry C18 150 x 4.6 mm, 5",
      "um; acetonitrile : 0.1 N acetic acid 90:10; UV 350 nm; LOQ ~6.25",
      "ug/mL ~ 10 ug/lungs). Assay measures total (free + encapsulated)",
      "tobramycin in the homogenate. See Marier 2002 Methods (Animal",
      "housing; Bacterial strain; Experimental infection and antibiotic",
      "treatment; Sample collection and bacterial growth; Analytical assay;",
      "Pharmacokinetic analysis) and Table 1 (Liposomal column) for the",
      "design and fitted model output."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters -- Marier 2002 Table 1 (Liposomal column).
    # NONMEM v5 ADVAN4 (depot + central + peripheral) parameterised in
    # rate-constant form because amounts in lung tissue (ug) -- not
    # concentrations -- were the dependent variable; consequently Vc is
    # not identifiable from the data and is not estimated. The paper's
    # k_23 / k_32 / k_20 labels follow NONMEM ADVAN4 numbering (depot=1,
    # central=2, peripheral=3); the canonical nlmixr2 micro-constant
    # names are k12 (central -> peripheral, paper's k_23), k21
    # (peripheral -> central, paper's k_32), and kel (elimination from
    # central, paper's k_20). The lung-tissue bioavailability F_L is
    # encoded as the depot bioavailability lfdepot.
    #
    # See vignette Assumptions and deviations for the amount-vs-
    # concentration encoding (Cc carries amount units, ug, in this
    # model).
    # ------------------------------------------------------------------
    lka     <- log(2.30)  ; label("Absorption rate constant ka into lung tissue (1/h)")                       # Marier 2002 Table 1: k_a = 2.30 /h (CV 30.8%)
    lfdepot <- log(0.85)  ; label("Lung bioavailability F_L (fraction of intratracheal dose reaching lung)")  # Marier 2002 Table 1: F_L = 0.85 (CV 12.2%)
    lk12    <- log(1.19)  ; label("Rate constant central -> peripheral (paper's k_23, 1/h)")                  # Marier 2002 Table 1: k_23 = 1.19 /h (CV 38.5%)
    lk21    <- log(0.407) ; label("Rate constant peripheral -> central (paper's k_32, 1/h)")                  # Marier 2002 Table 1: k_32 = 0.407 /h (CV 40.2%)
    lkel    <- log(0.155) ; label("Elimination rate from lung tissue (paper's k_20, 1/h)")                    # Marier 2002 Table 1: k_20 = 0.155 /h (CV 30.2%)

    # ------------------------------------------------------------------
    # IIV (BSV in paper) -- omega^2 = log(CV^2 + 1) conversion of the
    # Marier 2002 Table 1 CV%. Table 1 reports diagonal CV% only; no
    # between-parameter correlations are given, so all etas are
    # uncorrelated.
    #   ka:    CV 30.8% -> log(1 + 0.308^2)  = 0.09063
    #   F_L:   CV 12.2% -> log(1 + 0.122^2)  = 0.01478
    #   k_23:  CV 38.5% -> log(1 + 0.385^2)  = 0.13824
    #   k_32:  CV 40.2% -> log(1 + 0.402^2)  = 0.14978
    #   k_20:  CV 30.2% -> log(1 + 0.302^2)  = 0.08728
    # ------------------------------------------------------------------
    etalka     ~ 0.09063  # Marier 2002 Table 1: CV(k_a)   = 30.8%
    etalfdepot ~ 0.01478  # Marier 2002 Table 1: CV(F_L)   = 12.2%
    etalk12    ~ 0.13824  # Marier 2002 Table 1: CV(k_23)  = 38.5%
    etalk21    ~ 0.14978  # Marier 2002 Table 1: CV(k_32)  = 40.2%
    etalkel    ~ 0.08728  # Marier 2002 Table 1: CV(k_20)  = 30.2%

    # ------------------------------------------------------------------
    # Residual error. Methods declare a combined proportional + additive
    # error model; Results report a single combined population residual
    # variability of 18.1% for the liposomal arm. The individual
    # proportional and additive SDs are not separately reported, so the
    # combined residual is encoded as proportional only (the additive
    # component is documented in the vignette Assumptions and deviations
    # as a structural fact of the source paper). See vignette Errata.
    # ------------------------------------------------------------------
    propSd <- 0.181 ; label("Proportional residual variability on lung-tissue amount (fraction)")  # Marier 2002 Results paragraph after Table 1: 18.1% residual variability (liposomal)
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
    # parameterised in amount units (ug). The depot receives the
    # intratracheal bolus; bioavailability fdepot scales the fraction
    # that reaches the central (lung-tissue) compartment.
    # ------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                 k12 * central - k21 * peripheral1

    f(depot) <- fdepot

    # ------------------------------------------------------------------
    # 3. Observation. Vc was not fitted (see Methods 'Pharmacokinetic
    # analysis'); the central state holds the amount of tobramycin in
    # lung tissue (ug). Cc is conventionally the central observation
    # variable in nlmixr2lib; here it carries amount units, not
    # concentration. See vignette Assumptions and deviations for the
    # amount-vs-concentration mapping.
    # ------------------------------------------------------------------
    Cc <- central
    Cc ~ prop(propSd)
  })
}
