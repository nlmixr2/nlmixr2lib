Li_2015_abrilumab_cyno <- function() {
  description <- paste(
    "Preclinical (cynomolgus monkey). Two-compartment target-mediated drug",
    "disposition (TMDD) model with quasi-equilibrium (QE) binding and",
    "constant total receptor for abrilumab (AMG 181, MEDI7183), a human",
    "anti-alpha4beta7 integrin IgG2 monoclonal antibody, after single and",
    "repeated IV or SC doses of 0.01-80 mg/kg, with a direct inhibitory Emax",
    "model of free alpha4beta7 receptor on CD4+ T cells (Li 2015)."
  )
  reference <- paste(
    "Li H, Kock K, Wisler JA, Rees WA, Prince PJ, Reynhardt KO, Hsu H, Yu Z,",
    "Borie DC, Salinger DH, Pan WJ. Prediction of clinical pharmacokinetics",
    "of AMG 181, a human anti-alpha4beta7 monoclonal antibody for treating",
    "inflammatory bowel diseases. Pharmacol Res Perspect. 2015;3(1):e00098.",
    "doi:10.1002/prp2.98"
  )
  vignette <- "Li_2015_abrilumab"
  units <- list(time = "day", dosing = "mg", concentration = "nmol/L")

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "abrilumab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "abrilumab (free + alpha4beta7-bound)", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "abrilumab", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "cynomolgus monkey (Macaca fascicularis)",
    n_subjects = 119L,
    n_studies = 5L,
    weight_range = "2.1-6.4 kg",
    weight_median = "2.8 kg (mean 3.43 kg)",
    sex_female_pct = 38.7,
    disease_state = "healthy (pharmacology and toxicology studies)",
    dose_range = paste(
      "Single IV 0.01, 0.1, 0.3, 3 and 9 mg/kg; single SC 3 and 9 mg/kg;",
      "SC QW x 2 at 0.5 and 80 mg/kg; SC QW x 13 at 5, 20 and 80 mg/kg;",
      "IV QW x 13 at 80 mg/kg; SC QW x 24 at 20 and 80 mg/kg (Table 1)."
    ),
    regions = "Amgen preclinical studies A-E (two PK/PD, one non-GLP and two GLP toxicology studies)",
    notes = paste(
      "119 actively treated animals (73 male, 46 female; 86 SC, 33 IV),",
      "followed up to 420 days after the first dose. 120 of 1743",
      "concentrations (6.88%) that were anti-drug-antibody positive and",
      "visibly affected were excluded, leaving 1623 for modelling (Methods,",
      "'AMG 181 PK data for modeling'). The Emax PD model was fitted only to",
      "Study A (single IV 0.01-3 mg/kg or SC 3 mg/kg, 3 males per group):",
      "free alpha4beta7 mean fluorescence intensity on total CD4+ T cells",
      "against time-matched serum concentration."
    )
  )

  ini({
    # Li 2015 Table 2 (final two-compartment TMDD QE model, NONMEM 7.2
    # FOCE-I). Concentrations were modelled in nmol/L.
    lka <- log(0.846); label("First-order SC absorption rate constant KA (1/day)") # Table 2: KA = 0.846 /day (3.8 %SE)
    lcl <- log(12.9); label("Linear clearance CL (mL/day)") # Table 2: CL = 12.9 mL/day (5.9 %SE)
    lvc <- log(116); label("Central volume of distribution Vc (mL)") # Table 2: Vc = 116 mL (7.9 %SE)
    lk12 <- log(0.617); label("Central-to-peripheral rate constant K12 (1/day)") # Table 2: K12 = 0.617 /day (13 %SE)
    lk21 <- log(0.877); label("Peripheral-to-central rate constant K21 (1/day)") # Table 2: K21 = 0.877 /day (8.8 %SE)
    lkd <- fixed(log(0.013)); label("Equilibrium dissociation constant Kd (nmol/L)") # Table 2: Kd fixed to 0.013 nmol/L (in vitro EC50, monkey primary T cells)
    lrtot <- log(9.18); label("Total alpha4beta7 receptor concentration Rtot (nmol/L)") # Table 2: Rtot = 9.18 nmol/L (20 %SE)
    lkint <- log(0.0215); label("Internalization rate constant of the drug-receptor complex Kint (1/day)") # Table 2: Kint = 0.0215 /day (6.0 %SE)
    lfdepot <- fixed(log(1)); label("SC bioavailability F1 (fraction)") # Results: F1 'set to 1.0 (or 100%) to stabilize the model'

    # Emax model of free alpha4beta7 on total CD4+ T cells, Study A
    # (Results, 'Modeling of AMG 181 PD in cynomolgus monkeys'; equation
    # E = E0 * (1 - Emax * C / (EC50 + C)), Methods).
    le0 <- log(531); label("Baseline free alpha4beta7 receptor signal E0 (MFI)") # Results: E0 = 531 MFI
    lemax <- log(0.922); label("Maximum fractional reduction of free alpha4beta7 Emax (fraction)") # Results: Emax = 0.922
    lec50 <- log(0.0140); label("Abrilumab concentration at half-maximal effect EC50 (ug/mL)") # Results: EC50 = 0.0140 +/- 0.0082 ug/mL

    # IIV: Table 2 'IIV (%SE)' column, read as CV% of a log-normal eta,
    # omega^2 = log(1 + CV^2). KA, K12, K21, Kd and Kint have no IIV.
    etalcl + etalvc ~ c(
      log(1 + 0.441^2),
      0.642 * sqrt(log(1 + 0.441^2) * log(1 + 0.461^2)),
      log(1 + 0.461^2)
    ) # Table 2: IIV CL 44.1%, Vc 46.1%, correlation CL-Vc 0.642
    etalrtot ~ log(1 + 1.80^2) # Table 2: IIV Rtot 180%

    propSd <- 0.201; label("Proportional residual error (fraction)") # Table 2: proportional residual error 20.1% CV
    addSd <- 1.11; label("Additive residual error (nmol/L)") # Table 2: additive residual error SD 1.11 nmol/L
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    kd <- exp(lkd)
    rtot <- exp(lrtot + etalrtot)
    kint <- exp(lkint)
    fdepot <- exp(lfdepot)
    e0 <- exp(le0)
    emax <- exp(lemax)
    ec50 <- exp(lec50)

    # Molecular weight about 144 kDa (Methods, 'In vitro pharmacology
    # studies'). Doses and drug amounts are in mg and volumes in mL, so
    # mg/mL x 1e6 / 144 = nmol/L, and nmol/L x mL x 144 / 1e6 = mg.
    mw <- 144
    nm_per_mgml <- 1e6 / mw
    mg_per_nmml <- mw / 1e6

    # QE approximation with constant total receptor (Figure 1): the central
    # state holds total (free + bound) drug; free drug Cf solves
    # Ctot = Cf + Rtot * Cf / (Kd + Cf). Table 2 carries no receptor
    # degradation rate, so total receptor is held at Rtot (kdeg = kint).
    ctot <- central / vc * nm_per_mgml
    disc <- ctot - rtot - kd
    cfree <- 0.5 * (disc + sqrt(disc * disc + 4 * kd * ctot))
    cbound <- ctot - cfree

    # Linear elimination and distribution act on free drug; the
    # drug-receptor complex is internalized at Kint.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - ((cl / vc + k12) * cfree + kint * cbound) * vc * mg_per_nmml + k21 * peripheral1
    d/dt(peripheral1) <- k12 * cfree * vc * mg_per_nmml - k21 * peripheral1

    f(depot) <- fdepot

    # Free alpha4beta7 receptor MFI on total CD4+ T cells, driven by the
    # serum concentration in ug/mL (Methods Emax equation). Algebraic only:
    # the paper reports no PD variability.
    conc_ugml <- cfree * mw / 1000
    mfi_a4b7 <- e0 * (1 - emax * conc_ugml / (ec50 + conc_ugml))

    # The dependent variable is the free concentration Cf (Methods, 'PK
    # structural models'), in nmol/L. The human projection in Table 4
    # reproduces only on Cf, not on total drug.
    Cc <- cfree
    Cc ~ add(addSd) + prop(propSd)
  })
}
