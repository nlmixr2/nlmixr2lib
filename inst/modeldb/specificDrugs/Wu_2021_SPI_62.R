Wu_2021_SPI_62 <- function() {
  description <- paste(
    "Population target-mediated drug disposition (TMDD) PK model for the small-molecule",
    "11-beta-hydroxysteroid dehydrogenase type 1 (HSD-1) inhibitor SPI-62 in healthy adults",
    "(Wu 2021), fit to plasma concentrations alone. Oral absorption through a chain of four",
    "identical first-order transit steps (depot -> transit1 -> transit2 -> transit3 -> central,",
    "all governed by a single Ktr), two-compartment linear disposition (central + peripheral1)",
    "with clearance CL and distribution flow Q, and explicit second-order binding of SPI-62 in",
    "the central compartment to a high-affinity, low-capacity target (association rate constant",
    "Kon, dissociation rate constant Koff) to form a drug-target complex. The total target",
    "amount Rtotal is constant and the complex is not internalised (kint = 0), so free target",
    "is Rtotal minus the complex amount. The saturable binding explains SPI-62's nonlinear PK:",
    "very low plasma exposure after single low doses, dose-proportional PK at steady state, and",
    "unusually high accumulation ratios at low doses. Fit in NONMEM 7.4.3 (FOCEI, ADVAN13) to",
    "the <= 10 mg cohorts of the SPI-62 single-ascending-dose (1-10 mg) and low-dose",
    "multiple-ascending-dose (0.2-2 mg once daily) phase 1 trials. All volumes and clearances",
    "are apparent (per unit bioavailability) because F is unknown. Exponential inter-individual",
    "variability on Vcentral, CL, Ktr, Koff, and Rtotal; proportional residual error. No",
    "covariates retained (age, sex, body weight, and race were tested). Amounts are carried in",
    "nmol so that A / V is in nmol/L (= nM), the unit of Kon; convert an mg dose to nmol by",
    "multiplying by 1e6 / 424.4 (SPI-62 molecular weight 424.4 g/mol). The same group's later",
    "joint PK/PD refit of this structure is Wu_2023_SPI_62.",
    sep = " "
  )
  reference <- paste(
    "Wu N, Katz DA, An G. A target-mediated drug disposition model to explain nonlinear",
    "pharmacokinetics of the 11beta-hydroxysteroid dehydrogenase type 1 inhibitor SPI-62 in",
    "healthy adults. J Clin Pharmacol. 2021;61(11):1442-1453. doi:10.1002/jcph.1925.",
    "PMCID:PMC8596879.",
    sep = " "
  )
  vignette <- "Wu_2021_SPI_62"

  units <- list(
    time = "h",
    dosing = "nmol",
    concentration = "nmol/L",
    dosing_notes = paste(
      "Amounts are carried in nmol and volumes in L, so C = A / V is in nmol/L, identical to",
      "nM -- the unit in which Wu 2021 Table 2 reports Kon (nM^-1 h^-1). Rtotal is an amount",
      "(nmol), matching Table 2. Wu 2021 does not print a molecular weight; convert an mg dose",
      "to nmol with nmol = mg * 1e6 / 424.4 and a concentration from nM to ng/mL by",
      "multiplying by 0.4244. The 424.4 g/mol value is the one Wu 2023 (the same group's",
      "follow-up paper on the same drug) implies through its conversion 0.0787 nM = 0.0334",
      "ng/mL, and it agrees with Wu 2021's own statement that Rtotal = 6070 nmol 'corresponds",
      "to approximately 2.5 mg of SPI-62' (6070 nmol * 424.4 g/mol = 2.58 mg).",
      sep = " "
    )
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = FALSE means NOT checked against the source
  # paper.
  compartmentData <- list(
    depot = list(analyte = "SPI-62", units = "nmol", specimen = "administration site", verified = FALSE),
    transit1 = list(analyte = "SPI-62", units = "nmol", specimen = "administration site", verified = FALSE),
    transit2 = list(analyte = "SPI-62", units = "nmol", specimen = "administration site", verified = FALSE),
    transit3 = list(analyte = "SPI-62", units = "nmol", specimen = "administration site", verified = FALSE),
    central = list(analyte = "SPI-62", units = "nmol", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "SPI-62", units = "nmol", specimen = "plasma", verified = FALSE),
    complex = list(analyte = "drug-target complex", units = "nmol", specimen = "plasma", verified = FALSE)
  )

  # No covariates are used by this model: Wu 2021 Results ('Parameter
  # estimation') reports that age, body weight, race, and sex were tested on
  # the SPI-62 PK parameters and none showed a significant impact.
  covariateData <- list()

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at study entry.",
      units = "years",
      type = "continuous",
      notes = paste(
        "Tested by forward addition / backward elimination (Wu 2021 Methods, 'Covariate",
        "model') but not retained ('none of them showed any significant impact'). Analysis",
        "cohort age range 20-54 years.",
        sep = " "
      ),
      source_name = "age"
    ),
    WT = list(
      description = "Body weight at study entry.",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Tested but not retained (Wu 2021 Results). No allometric scaling is applied.",
        "Analysis cohort mean +/- SD 76.5 +/- 12.1 kg.",
        sep = " "
      ),
      source_name = "body weight"
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = male",
      notes = "Tested as 'sex' but not retained (Wu 2021 Results). 33 male and 11 female subjects.",
      source_name = "sex"
    ),
    RACE_WHITE = list(
      description = "White race indicator (1 = White, 0 = non-White).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = non-White",
      notes = paste(
        "Race was tested but not retained (Wu 2021 Results). Wu 2021 does not report the race",
        "composition or how race levels were grouped; no race term enters model().",
        sep = " "
      ),
      source_name = "race"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 44L,
    n_studies = 2L,
    n_observations = 996L,
    study_names = c(
      "SAD/FE -- SPI-62 first-in-human single-ascending-dose and food-effect trial; the fasted 1, 3, 6, and 10 mg cohorts (n = 6 active per cohort) only",
      "MAD part B -- SPI-62 low-dose multiple-ascending-dose cohorts (0.2, 0.4, 0.7, and 2 mg once daily)"
    ),
    age_range = "20-54 years",
    weight_range = "mean +/- SD 76.5 +/- 12.1 kg (range not reported)",
    sex_female_pct = 25.0,
    disease_state = "Healthy adult volunteers.",
    dose_range = paste(
      "SAD: 1, 3, 6, and 10 mg single oral doses, fasted. MAD part B: 3 mg loading dose on day",
      "1 then 0.2 mg once daily on days 2-14 (n = 4); 0.4 mg once daily on days 1-14 (n = 4);",
      "0.7 mg or 2 mg single dose on day 1, 6-day washout, then once daily on days 7-20 (n = 6",
      "each). Higher-dose cohorts (> 10 mg) were excluded from the analysis.",
      sep = " "
    ),
    regions = "Not reported.",
    notes = paste(
      "996 plasma concentrations (774 above and 222 below the LLOQ; BLQ replaced with LLOQ/2",
      "per the Discussion). LLOQ 0.1 ng/mL in the SAD/FE trial and 0.004 ng/mL in MAD part B.",
      "Wu 2021 Methods, 'Data Source', and Table 1.",
      sep = " "
    )
  )

  ini({
    # Absorption: four identical sequential first-order transit steps (depot ->
    # transit1 -> transit2 -> transit3 -> central) sharing one rate constant
    # Ktr. Wu 2021 Eqs. 1-5 and Figure 1.
    lktr <- log(8.52); label("Transit absorption rate constant Ktr (1/h)") # Wu 2021 Table 2 (Ktr = 8.52 1/h, RSE 11%)

    # Two-compartment linear disposition. Table 2 footnote a: apparent
    # parameters because bioavailability F is unknown.
    lcl <- log(10.1); label("Apparent clearance CL/F (L/h)") # Wu 2021 Table 2 (CL = 10.1 L/h, RSE 6%)
    lvc <- log(141); label("Apparent central volume of distribution Vcentral/F (L)") # Wu 2021 Table 2 (Vcentral = 141 L, RSE 16%)
    lq <- log(2.31); label("Apparent distribution flow Q/F (L/h)") # Wu 2021 Table 2 (Q = 2.31 L/h, RSE 12%)
    lvp <- log(114); label("Apparent peripheral volume of distribution Vperipheral/F (L)") # Wu 2021 Table 2 (Vperipheral = 114 L, RSE 7%)

    # Target binding (TMDD). Concentrations are nM and Rtotal is an amount in
    # nmol, so Kon * C * (Rtotal - RC) is in nmol/h. Wu 2021 Eqs. 5 and 7.
    lkon <- log(7.1); label("Second-order association rate constant Kon (1/(nM*h))") # Wu 2021 Table 2 (Kon = 7.1 nM^-1 h^-1, RSE 7%)
    lkoff <- log(0.249); label("First-order dissociation rate constant Koff (1/h)") # Wu 2021 Table 2 (Koff = 0.249 1/h, RSE 24%); Discussion: Kd = Koff/Kon = 35.1 pM
    lrtot <- log(6070); label("Total amount of target binding sites Rtotal (nmol)") # Wu 2021 Table 2 (Rtotal = 6070 nmol, RSE 9%)

    # Inter-individual variability, exponential model (Wu 2021 Eq. 8). Table 2
    # reports each IIV as a percent; it is read as a coefficient of variation
    # and converted with omega^2 = log(CV^2 + 1). Diagonal only -- no
    # covariances are reported.
    #   Vcentral CV 54.3% -> log(1 + 0.543^2) = 0.258394
    #   CL       CV 20.6% -> log(1 + 0.206^2) = 0.041560
    #   Ktr      CV 50.6% -> log(1 + 0.506^2) = 0.227961
    #   Koff     CV  116% -> log(1 + 1.16^2)  = 0.852541
    #   Rtotal   CV 35.8% -> log(1 + 0.358^2) = 0.120592
    etalvc ~ 0.258394 # Wu 2021 Table 2 (IIV Vcentral = 54.3%, RSE 27%, shrinkage 14%)
    etalcl ~ 0.041560 # Wu 2021 Table 2 (IIV CL = 20.6%, RSE 59%, shrinkage 22%)
    etalktr ~ 0.227961 # Wu 2021 Table 2 (IIV Ktr = 50.6%, RSE 29%, shrinkage 8%)
    etalkoff ~ 0.852541 # Wu 2021 Table 2 (IIV Koff = 116%, RSE 50%, shrinkage 15%)
    etalrtot ~ 0.120592 # Wu 2021 Table 2 (IIV Rtotal = 35.8%, RSE 38%, shrinkage 12%)

    # Proportional residual error, Wu 2021 Eq. 9 (C = Cpred * (1 + eps)).
    propSd <- 0.268; label("Proportional residual error (fraction)") # Wu 2021 Table 2 (proportional residual variability = 26.8%, RSE 2%, shrinkage 7%)
  })

  model({
    # Units: amounts nmol, volumes L, time h, so A / V is nmol/L == nM,
    # matching Kon. Doses passed to rxode2::et() must be in nmol.

    # Individual parameters. Exponential IIV on Vcentral, CL, Ktr, Koff,
    # and Rtotal; none on Q, Vperipheral, or Kon (Wu 2021 Table 2).
    ktr <- exp(lktr + etalktr)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)

    kon <- exp(lkon)
    koff <- exp(lkoff + etalkoff)
    rtot <- exp(lrtot + etalrtot)

    # Central-compartment SPI-62 concentration (nM).
    Cc <- central / vc

    # Target binding flux (nmol/h). Free target is Rtotal - RC because the
    # total number of binding sites is constant and kint = 0 (Wu 2021
    # Methods, 'Structural model').
    bind_rate <- kon * Cc * (rtot - complex)
    unbind_rate <- koff * complex

    # ODE system, Wu 2021 Eqs. 1-7. All initial conditions are 0 except the
    # depot, which receives the dose (Adepot(0) = Dose). Eq. 1 is printed as
    # -ktr * Dose; the depot amount Adepot is used here, which is the
    # standard first-order form and the only one consistent with Eq. 2.
    d/dt(depot) <- -ktr * depot # Eq. 1
    d/dt(transit1) <- ktr * depot - ktr * transit1 # Eq. 2
    d/dt(transit2) <- ktr * transit1 - ktr * transit2 # Eq. 3
    d/dt(transit3) <- ktr * transit2 - ktr * transit3 # Eq. 4
    d/dt(central) <- ktr * transit3 - bind_rate + unbind_rate -
      (cl / vc) * central -
      (q / vc) * central + (q / vp) * peripheral1 # Eq. 5
    d/dt(peripheral1) <- (q / vc) * central - (q / vp) * peripheral1 # Eq. 6
    d/dt(complex) <- bind_rate - unbind_rate # Eq. 7

    Cc ~ prop(propSd)
  })
}
