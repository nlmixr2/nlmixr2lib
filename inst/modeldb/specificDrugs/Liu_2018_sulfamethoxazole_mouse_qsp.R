Liu_2018_sulfamethoxazole_mouse_qsp <- function() {
  description <- "QSP. Preclinical (mouse, C3H/HeN). Integrated quantitative systems pharmacology model of trimethoprim-sulfamethoxazole (TMP-SMX) against Pneumocystis murina. A three-compartment PK module (administration compartment -> plasma -> peripheral tissue, with first-order decay from both plasma and tissue and the two dimensionless scaling factors RAP and RTP that Liu 2018 uses instead of volumes) drives the two-stage Pneumocystis life-cycle PD module of trophic forms and asci. Because the TMP:SMX ratio was fixed at 1:5 in the constraining data, Liu 2018 uses the sulfamethoxazole level alone as the proxy for the combination, so this model carries a single drug. TMP-SMX represses folate synthesis and therefore genome replication, so unlike the echinocandins it acts on all three PD rates: it inhibits trophic proliferation (complete inhibition at saturating concentration), multiplies the second-order trophic death rate by up to 1 + 650, and adds a sigmoidal Emax increment to the asci death rate. The effect is driven not by the current plasma level but by the plasma level seven days earlier (SMXeff = SMXpla(t - tau), tau = 7 days), a genuine delay differential equation encoded with rxode2's delay(); this reproduces the delayed onset the authors observed in Figure 4d. All states are concentrations or counts, not amounts -- the PK module carries no volume term, so a dose is administered as an initial concentration increment (500 ug/mL for 200 mg/kg, 550 ug/mL for 250 mg/kg; Liu 2018 Table 3 footnote). PK rate constants are reported in 1/h and PD rate constants in 1/day, so the PK right-hand side is scaled by 24 to put the whole model on a day time base. Drug-free, this file reduces exactly to modellib('Liu_2018_pneumocystis_mouse_qsp')."
  reference <- paste(
    "Liu GS, Ballweg R, Ashbaugh A, Zhang Y, Facciolo J, Cushion MT, Zhang T.",
    "(2018). A quantitative systems pharmacology (QSP) model for Pneumocystis",
    "treatment in mice. BMC Syst Biol 12(1):77.",
    "doi:10.1186/s12918-018-0603-9.",
    "PK construction and validation data (Liu 2018 Table 2 / reference 41) from",
    "Misiek M, Buck RE, Pursiano TA, Chisholm DR, Tsai YH, Price KE,",
    "Leitner F. (1985). Antibacterial activity of phosphanilic acid, alone and",
    "in combination with trimethoprim.",
    "Antimicrob Agents Chemother 28(6):761-765. doi:10.1128/AAC.28.6.761.",
    "Organism-burden data constraining the drug effect (Liu 2018 reference 33) from",
    "Cushion MT, Linke MJ, Ashbaugh A, Sesterhenn T, Collins MS, Lynch K,",
    "Brubaker R, Walzer PD. (2010). Echinocandin treatment of pneumocystis",
    "pneumonia in rodent models depletes cysts leaving trophic burdens that",
    "cannot transmit the infection. PLoS One 5(1):e8524.",
    "doi:10.1371/journal.pone.0008524.",
    "Folate-synthesis mechanism (Liu 2018 reference 13) from",
    "Huang L, Crothers K, Atzori C, Benfield T, Miller R, Rabodonirina M,",
    "Helweg-Larsen J. (2004). Dihydropteroate synthase gene mutations in",
    "pneumocystis and sulfa resistance. Emerg Infect Dis 10(10):1721-1728.",
    "doi:10.3201/eid1010.030994."
  )
  vignette <- "Liu_2018_pneumocystis_qsp"

  # The two Pneumocystis life-cycle stages are paper-mechanistic states with no
  # canonical analogue in inst/references/compartment-names.md: `trophic` is the
  # vegetative, actively proliferating trophic form and `asci` the
  # non-proliferating ascus (cyst) form.
  paper_specific_compartments <- c("trophic", "asci")

  units <- list(
    time = "day",
    dosing = "ug/mL (an initial concentration increment, not a mass -- see description)",
    concentration = "ug/mL for the plasma sulfamethoxazole observation Cc; organisms per mouse lung for trophic and asci"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE -- checked against Liu 2018 Methods
  # 'Construction of the PK module in mice', Figure 2 and Table 1.
  compartmentData <- list(
    depot = list(
      analyte = "sulfamethoxazole",
      units = "ug/mL (the AC state is a concentration; the module has no volume term)",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "sulfamethoxazole",
      units = "ug/mL",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "sulfamethoxazole",
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
    n_subjects = "PK module built on published murine plasma sulfamethoxazole profiles digitised from Misiek 1985 (construction and validation); PD and drug-effect layers constrained by the Pneumocystis burdens of Cushion 2010 and by the authors' own 6-week-old male C3H/HeN mice",
    n_studies = 2L,
    age_range = "6 weeks old at the start of the authors' own experiments",
    sex_female_pct = 0,
    disease_state = "Pneumocystis murina lung infection in immunosuppressed mice (the model of Pneumocystis pneumonia, PCP)",
    dose_range = "200 and 250 mg/kg sulfamethoxazole orally, giving initial administration-compartment concentrations of 500 and 550 ug/mL respectively (Liu 2018 Table 3 footnote); the PK validation profile of Figure 5d is an oral dose of 50 mg/kg, and the Figure 4 treatment arm is 200 mg/kg three times a week for three weeks",
    co_medication = "trimethoprim, co-administered at a fixed TMP:SMX ratio of 1:5 in the constraining data; Liu 2018 therefore models the sulfamethoxazole level alone as a proxy for the combination",
    regions = "United States (University of Cincinnati; mice supplied by Charles River)",
    notes = "Liu 2018 Experimental methods: 6-week-old male C3H/HeN mice. Organism burden quantified by RT-qPCR of Pneumocystis mitochondrial large-subunit rRNA (total nuclei) or by microscopic quantification with cresyl echt violet (asci) and a rapid Wright-Giemsa stain (all stages)."
  )

  ini({
    # ==================================================================
    # PK module -- Liu 2018 Table 1, TMP/SMX column. All rate constants
    # are reported in 1/h; the value inside each log() below is the
    # Table 1 entry verbatim in the paper's own 1/h units, and model()
    # scales the whole PK right-hand side by 24 h/day.
    #
    # Every value was chosen manually by the authors ('all parameters sets
    # and initial conditions were derived manually using a trial and error
    # method to find a plausible parameter set that visually recaptures the
    # experimentally observed data'), so none carries an uncertainty
    # estimate and all are encoded as fixed(). Goodness of fit for this
    # parameter set was SSE = 131.15 (Table 2).
    # ==================================================================
    lka <- fixed(log(5))
    label("Absorption rate constant Kabs, administration compartment to plasma (1/h)")  # Table 1
    lk12 <- fixed(log(0.17))
    label("Plasma to peripheral tissue transfer rate constant KPT (1/h)")  # Table 1
    lk21 <- fixed(log(5))
    label("Peripheral tissue to plasma transfer rate constant KTP (1/h)")  # Table 1
    lkel <- fixed(log(0.2))
    label("Plasma decay rate constant KdP (1/h)")  # Table 1
    lkdt <- fixed(log(0.2))
    label("Peripheral tissue decay rate constant KdT (1/h)")  # Table 1
    rap <- fixed(3)
    label("Administration-compartment scaling factor RAP (dimensionless)")  # Table 1
    rtp <- fixed(0.01)
    label("Peripheral tissue scaling factor RTP (dimensionless)")  # Table 1

    # ==================================================================
    # TMP/SMX effects -- Liu 2018 Table 5, TMP/SMX block. One EC50 and one
    # Hill coefficient gate all three effects. MEAsci is ADDED to the asci
    # death rate (in 1/day) and so carries units of 1/day; METro MULTIPLIES
    # the second-order trophic death rate and is dimensionless. The
    # inhibition of trophic proliferation has no maximal-effect parameter
    # (complete inhibition at saturating concentration).
    # ==================================================================
    ec50 <- fixed(0.2)
    label("Delayed plasma sulfamethoxazole concentration giving half the maximal effects Ec50SMX (ug/mL)")  # Table 5
    hill <- fixed(2)
    label("Hill coefficient of the sulfamethoxazole effects n (dimensionless)")  # Table 5
    emax_asci <- fixed(0.75)
    label("Maximal SMX-induced increment in the asci death rate MEAsci (1/day)")  # Table 5
    emax_trophic <- fixed(650)
    label("Maximal fractional SMX-induced increase in the trophic death rate METro (dimensionless)")  # Table 5
    tau <- fixed(7)
    label("Delay between the plasma SMX level and its pharmacodynamic effect tau (day)")  # Table 5

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
    # rate constant in 1/day (Table 4) and the effect delay tau in days
    # (Table 5). The integrated QSP model is run on a day time base -- the
    # readouts are at day 35 and day 56 -- so the entire PK right-hand
    # side, each term of which is first order in a 1/h rate constant, is
    # multiplied by the 24 h in a day. This is a pure unit conversion; it
    # changes no reported value.
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
    # (Table 3 footnote). TMP-SMX was administered orally, so the dose goes
    # into `depot`, the paper's administration compartment (AC).
    Cc <- central

    # ==================================================================
    # Delay in the SMX effect -- Liu 2018 Table 5, 'Delay in SMX effect':
    #
    #   SMXeff(t) = SMXpla(t - tau),  tau = 7 days
    #
    # 'in comparison with the rapid antifungal effect of anidulafungin, the
    # experimental evidence suggests that the effect of TMP-SMX was
    # delayed. This time delay was incorporated into our QSP model.'
    #
    # delay() interpolates the plasma state at t - tau from the solver's
    # dense output, so this is a genuine delay differential equation rather
    # than a transit-chain approximation. Before t = tau the constant
    # initial condition central(0) = 0 is returned, which is exactly right:
    # no drug had been given.
    # ==================================================================
    smxeff <- delay(central, tau)

    # ==================================================================
    # TMP/SMX effects -- Liu 2018 Table 5. 'TMP-SMX represses folate
    # synthesis which is essential for genome replication in the organism.
    # Therefore, in our simplified model, TMP-SMX was assumed to inhibit
    # the proliferation rate of the trophic forms and increase the death
    # rates of both the trophic forms and asci.'
    #
    #   vsTro  = ks*(1 - SMXeff^n/(SMXeff^n + Ec50SMX^n))
    #   vdTro  = kdTro*(1 + METro*SMXeff^n/(SMXeff^n + Ec50SMX^n))
    #   vdAsci = kdAsci + MEAsci*SMXeff^n/(SMXeff^n + Ec50SMX^n)
    # ==================================================================
    smxfrac <- smxeff^hill / (smxeff^hill + ec50^hill)
    vstro <- kstro * (1 - smxfrac)
    vdtro <- kdtro * (1 + emax_trophic * smxfrac)
    vdasci <- kdasci + emax_asci * smxfrac

    # ==================================================================
    # Two-stage Pneumocystis life cycle -- Liu 2018 Table 4 / Equations (d)
    # and (e), with ks, kdTro and kdAsci replaced by their drug-dependent
    # counterparts vsTro, vdTro and vdAsci. TMP-SMX does not alter the
    # transformation rates kTA and kAT. The trophic loss term is
    # deliberately second order (logistic crowding); the asci loss term is
    # first order.
    # ==================================================================
    d/dt(trophic) <- vstro * trophic - vdtro * trophic * trophic -
      kta * trophic + kat * asci
    d/dt(asci) <- kta * trophic - kat * asci - vdasci * asci

    trophic(0) <- trophic0
    asci(0) <- asci0

    # Liu 2018 plots both stages on a log10 scale (Figures 4a-4d).
    log10_trophic <- log10(trophic)
    log10_asci <- log10(asci)
  })
}
