Liu_2018_anidulafungin_mouse_qsp <- function() {
  description <- "QSP. Preclinical (mouse, C3H/HeN). Integrated quantitative systems pharmacology model of anidulafungin against Pneumocystis murina. A three-compartment PK module (administration compartment -> plasma -> peripheral tissue, with first-order decay from both plasma and tissue and the two dimensionless scaling factors RAP and RTP that Liu 2018 uses instead of volumes) drives the two-stage Pneumocystis life-cycle PD module of trophic forms and asci. Echinocandins block cell-wall construction in the ascus, so anidulafungin acts only on the asci: it adds a sigmoidal Emax increment to the asci death rate and inhibits the trophic-to-asci transformation with the same EC50, leaving the trophic burden essentially untouched. All states are concentrations or counts, not amounts -- the PK module carries no volume term, so a dose is administered as the initial plasma or administration-compartment concentration increment tabulated in Liu 2018 Table 3. PK rate constants are reported in 1/h and PD rate constants in 1/day, so the PK right-hand side is scaled by 24 to put the whole model on a day time base. Reproduces the Figure 4a dose ranking in which anidulafungin reduces the day-56 asci burden at 1, 0.5 and 0.1 mg/kg while the trophic form is unchanged (Figure 4b). Drug-free, this file reduces exactly to modellib('Liu_2018_pneumocystis_mouse_qsp')."
  reference <- paste(
    "Liu GS, Ballweg R, Ashbaugh A, Zhang Y, Facciolo J, Cushion MT, Zhang T.",
    "(2018). A quantitative systems pharmacology (QSP) model for Pneumocystis",
    "treatment in mice. BMC Syst Biol 12(1):77.",
    "doi:10.1186/s12918-018-0603-9.",
    "PK construction data (Liu 2018 Table 2 / reference 35) from",
    "Gumbo T, Drusano GL, Liu W, Ma L, Deziel MR, Drusano MF, Louie A. (2006).",
    "Anidulafungin pharmacokinetics and microbial response in neutropenic mice",
    "with disseminated candidiasis. Antimicrob Agents Chemother 50(11):3695-3700.",
    "doi:10.1128/AAC.00507-06.",
    "PK validation data (Liu 2018 Table 2 / reference 42) from",
    "Andes D, Diekema DJ, Pfaller MA, Prince RA, Marchillo K, Ashbeck J, Hou J.",
    "(2008). In vivo pharmacodynamic characterization of anidulafungin in a",
    "neutropenic murine candidiasis model.",
    "Antimicrob Agents Chemother 52(2):539-550. doi:10.1128/AAC.01061-07.",
    "Organism-burden data constraining the drug effect (Liu 2018 reference 33) from",
    "Cushion MT, Linke MJ, Ashbaugh A, Sesterhenn T, Collins MS, Lynch K,",
    "Brubaker R, Walzer PD. (2010). Echinocandin treatment of pneumocystis",
    "pneumonia in rodent models depletes cysts leaving trophic burdens that",
    "cannot transmit the infection. PLoS One 5(1):e8524.",
    "doi:10.1371/journal.pone.0008524."
  )
  vignette <- "Liu_2018_pneumocystis_qsp"

  # The two Pneumocystis life-cycle stages are paper-mechanistic states with no
  # canonical analogue in inst/references/compartment-names.md: `trophic` is the
  # vegetative, actively proliferating trophic form and `asci` the
  # non-proliferating ascus (cyst) form the echinocandins target.
  paper_specific_compartments <- c("trophic", "asci")

  units <- list(
    time = "day",
    dosing = "ug/mL (an initial concentration increment, not a mass -- see description)",
    concentration = "ug/mL for the plasma anidulafungin observation Cc; organisms per mouse lung for trophic and asci"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE -- checked against Liu 2018 Methods
  # 'Construction of the PK module in mice', Figure 2 and Table 1.
  compartmentData <- list(
    depot = list(
      analyte = "anidulafungin",
      units = "ug/mL (the AC state is a concentration; the module has no volume term)",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "anidulafungin",
      units = "ug/mL",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "anidulafungin",
      units = "ug/mL",
      specimen = "tissue",
      verified = TRUE
    ),
    trophic = list(
      analyte = "Pneumocystis murina trophic form",
      units = "organisms per mouse lung",
      specimen = "not applicable",
      verified = TRUE
    ),
    asci = list(
      analyte = "Pneumocystis murina asci (cyst form)",
      units = "organisms per mouse lung",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  # No subject-level covariates. Liu 2018 induced variability by resampling the
  # PD rate constants from a uniform distribution spanning 70-130% of the
  # Table 4 basal values, not by covariate effects.
  covariateData <- list()

  population <- list(
    species = "mouse (C3H/HeN)",
    n_subjects = "PK module built on published murine plasma-concentration profiles digitised from Gumbo 2006 (construction) and Andes 2008 (validation); PD and drug-effect layers constrained by the Pneumocystis burdens of Cushion 2010 and by the authors' own 6-week-old male C3H/HeN mice",
    n_studies = 3L,
    age_range = "6 weeks old at the start of the authors' own experiments",
    sex_female_pct = 0,
    disease_state = "Pneumocystis murina lung infection in immunosuppressed mice (the model of Pneumocystis pneumonia, PCP)",
    dose_range = "0.1, 0.5, 1, 2.5, 5 and 10 mg/kg intraperitoneally (Liu 2018 Table 3); the Figure 4a treatment arms are 0.1, 0.5 and 1 mg/kg three times a week for three weeks",
    regions = "United States (University of Cincinnati; mice supplied by Charles River)",
    notes = "Liu 2018 Experimental methods: 6-week-old male C3H/HeN mice. Organism burden quantified by RT-qPCR of Pneumocystis mitochondrial large-subunit rRNA (total nuclei) or by microscopic quantification with cresyl echt violet (asci) and a rapid Wright-Giemsa stain (all stages)."
  )

  ini({
    # ==================================================================
    # PK module -- Liu 2018 Table 1, anidulafungin column. All rate
    # constants are reported in 1/h; the value inside each log() below is
    # the Table 1 entry verbatim in the paper's own 1/h units, and model()
    # scales the whole PK right-hand side by 24 h/day.
    #
    # Every value was chosen manually by the authors ('all parameters sets
    # and initial conditions were derived manually using a trial and error
    # method to find a plausible parameter set that visually recaptures the
    # experimentally observed data'), so none carries an uncertainty
    # estimate and all are encoded as fixed(). Goodness of fit for this
    # parameter set was SSE = 459.26 (Table 2).
    # ==================================================================
    lka <- fixed(log(5))
    label("Absorption rate constant Kabs, administration compartment to plasma (1/h)")  # Table 1
    lk12 <- fixed(log(1.5))
    label("Plasma to peripheral tissue transfer rate constant KPT (1/h)")  # Table 1
    lk21 <- fixed(log(5))
    label("Peripheral tissue to plasma transfer rate constant KTP (1/h)")  # Table 1
    lkel <- fixed(log(0.035))
    label("Plasma decay rate constant KdP (1/h)")  # Table 1
    lkdt <- fixed(log(0.035))
    label("Peripheral tissue decay rate constant KdT (1/h)")  # Table 1
    rap <- fixed(3)
    label("Administration-compartment scaling factor RAP (dimensionless)")  # Table 1
    rtp <- fixed(0.2)
    label("Peripheral tissue scaling factor RTP (dimensionless)")  # Table 1

    # ==================================================================
    # Anidulafungin effect on the asci -- Liu 2018 Table 5, echinocandin
    # block. ME is ADDED to the asci death rate (which is in 1/day), so ME
    # carries units of 1/day; the same EC50 and Hill coefficient also gate
    # the inhibition of asci formation, which has no separate maximal
    # effect (complete inhibition at saturating concentration).
    # ==================================================================
    ec50 <- fixed(0.039)
    label("Plasma anidulafungin concentration giving half the maximal asci effects Ec50 (ug/mL)")  # Table 5
    emax <- fixed(0.42)
    label("Maximal anidulafungin-induced increment in the asci death rate ME (1/day)")  # Table 5
    hill <- fixed(1)
    label("Hill coefficient of the anidulafungin asci effects n (dimensionless)")  # Table 5

    # ==================================================================
    # PD module -- Liu 2018 Table 4 'Basal Parameter Values', shared by
    # all four integrated models.
    # ==================================================================
    kstro <- fixed(1)
    label("Trophic form proliferation rate KsTro (1/day)")  # Table 4
    kdtro <- fixed(1e-7)
    label("Trophic form second-order death rate KdTro (1/day/organism)")  # Table 4
    kta <- fixed(0.1)
    label("Trophic form to asci transformation rate KTA (1/day)")  # Table 4
    kat <- fixed(0.1)
    label("Asci to trophic form transformation rate KAT (1/day)")  # Table 4
    kdasci <- fixed(2e-12)
    label("Asci death rate KdAsci (1/day)")  # Table 4

    # Initial organism burden. NOT REPORTED BY LIU 2018; back-solved from
    # the paper's own description of the untreated time course (slow
    # accumulation over two weeks, exponential growth from week three,
    # 'about 35 days, the levels of both trophic forms and asci reached a
    # steady state of about 10^7'). A unit inoculum reproduces that
    # timeline. See the vignette 'Assumptions and deviations'.
    trophic0 <- fixed(1)
    label("Initial trophic form burden (organisms per lung)")  # back-solved; not reported
    asci0 <- fixed(1)
    label("Initial asci burden (organisms per lung)")  # back-solved; not reported
  })

  model({
    # ==================================================================
    # Liu 2018 reports every PK rate constant in 1/h (Table 1) but every PD
    # rate constant in 1/day (Table 4). The integrated QSP model is run on
    # a day time base -- the readouts are at day 35 and day 56 -- so the
    # entire PK right-hand side, each term of which is first order in a 1/h
    # rate constant, is multiplied by the 24 h in a day. This is a pure
    # unit conversion; it changes no reported value.
    # ==================================================================
    hrday <- 24

    # ==================================================================
    # Three-compartment PK module -- Liu 2018 Table 1 / Equations (a)-(c).
    #
    #   dDrugAC/dt = -RAP*Kabs*DrugAC
    #   dDrugP/dt  =  Kabs*DrugAC - (KPT + KdP)*DrugP + KTP*DrugT
    #   dDrugT/dt  =  RTP*(-KTP*DrugT + KPT*DrugP) - KdT*DrugT
    #
    # Note the deliberate asymmetries, which are reproduced verbatim rather
    # than "corrected" into a mass-conserving system: RAP multiplies only
    # the loss from the administration compartment and not the matching
    # gain in plasma (so a fraction 1/RAP of the administered
    # concentration reaches plasma), and RTP multiplies both transfer terms
    # of the tissue equation but neither of the plasma equation (so RTP
    # plays the role a plasma-to-tissue volume ratio would). Liu 2018 calls
    # both 'non-dimensional scaling factor[s]' whose 'values ... are
    # estimated from the observed data for each drug'.
    # ==================================================================
    # The five PK rate constants are log-parameterised in ini() (the
    # library convention for strictly-positive fixed-effect PK
    # parameters); the value inside each log() is the Table 1 entry in
    # the paper's own 1/h units.
    ka <- exp(lka)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    kel <- exp(lkel)
    kdt <- exp(lkdt)

    d/dt(depot) <- (-rap * ka * depot) * hrday
    d/dt(central) <- (ka * depot - (k12 + kel) * central +
      k21 * peripheral1) * hrday
    d/dt(peripheral1) <- (rtp * (-k21 * peripheral1 + k12 * central) -
      kdt * peripheral1) * hrday

    # The plasma state IS the plasma concentration: the module carries no
    # volume, and a dose is given as an initial concentration increment
    # (Table 3). An i.v. dose is administered into `central` directly ('in
    # order to mimic i.v. injection, we elevated the initial level of the
    # drug in the plasma compartment'); an i.p. or p.o. dose goes into
    # `depot`, the paper's administration compartment (AC).
    Cc <- central

    # ==================================================================
    # Echinocandin effect -- Liu 2018 Table 5. 'Echinocandins ... block the
    # construction of the cellular wall of the asci. Therefore, this family
    # of drugs were assumed to reduce the level of asci by promoting their
    # death as well as inhibiting their formation.'
    #
    #   vdAsci = kdAsci + ME*Echipla^n/(Echipla^n + Ec50Echi^n)
    #   vTA    = kTA*(1 - Echipla^n/(Echipla^n + Ec50Echi^n))
    #
    # Echipla is 'the current plasma level of echinocandin', i.e. Cc; there
    # is no effect on the trophic proliferation or trophic death rates.
    # ==================================================================
    echieff <- Cc^hill / (Cc^hill + ec50^hill)
    vdasci <- kdasci + emax * echieff
    vta <- kta * (1 - echieff)

    # ==================================================================
    # Two-stage Pneumocystis life cycle -- Liu 2018 Table 4 / Equations (d)
    # and (e), with kTA and kdAsci replaced by their drug-dependent
    # counterparts vTA and vdAsci. The trophic loss term is deliberately
    # second order (logistic crowding); the asci loss term is first order.
    # ==================================================================
    d/dt(trophic) <- kstro * trophic - kdtro * trophic * trophic -
      vta * trophic + kat * asci
    d/dt(asci) <- vta * trophic - kat * asci - vdasci * asci

    trophic(0) <- trophic0
    asci(0) <- asci0

    # Liu 2018 plots both stages on a log10 scale (Figures 4a-4d).
    log10_trophic <- log10(trophic)
    log10_asci <- log10(asci)
  })
}
