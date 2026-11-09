Kovalenko_2025_alirocumab_adult <- function() {
  description <- "Full target-mediated drug disposition (TMDD) population PK model for total alirocumab and free PCSK9 in healthy adults from two densely sampled phase I studies (Kovalenko 2025, adult model). Two-compartment disposition of free alirocumab, free PCSK9 and both alirocumab:PCSK9 complexes (1:1 and 1:2) sharing the same intercompartmental rate constants; SC absorption through one depot-transit compartment; two-step binding of alirocumab to one and then a second PCSK9 molecule in the central compartment; zero-order PCSK9 synthesis. Weight on Vc (power, 75 kg reference); kon and koff fixed to Biacore values. Its absorption, distribution and bioavailability estimates were carried into the pediatric model (Kovalenko_2025_alirocumab_pediatric)."
  reference <- "Kovalenko P, Harnisch L, Mendell J, Wang Y, Davis JD, DiCioccio AT. Target-Mediated Modeling of Alirocumab in Adolescents and Children >=8 to <12 Years of Age Using Phase II and III Data. Clin Pharmacol Drug Dev. 2025;14(5):347. doi:10.1002/cpdd.1523"
  vignette <- "Kovalenko_2025_alirocumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  # complex_db is the 1:2 alirocumab:PCSK9:PCSK9 complex (the paper's L:R:R)
  # formed when the second arm of the bivalent antibody binds a second PCSK9;
  # complex_db_peripheral1 is its peripheral-compartment counterpart. The 1:1
  # complex is the canonical `complex`.
  paper_specific_compartments <- c("complex_db", "complex_db_peripheral1")

  covariateData <- list(
    WT = list(
      description = "Body weight at baseline",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on Vc; Table 2 footnote states the adult volumes are estimated at 75 kg.",
      source_name = "weight"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "alirocumab", units = "nmol", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "alirocumab", units = "nmol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "alirocumab (free)", units = "nmol", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "alirocumab (free)", units = "nmol", specimen = "not applicable", verified = TRUE),
    target = list(analyte = "PCSK9 (free)", units = "nmol", specimen = "serum", verified = TRUE),
    target_peripheral1 = list(analyte = "PCSK9 (free)", units = "nmol", specimen = "not applicable", verified = TRUE),
    complex = list(analyte = "alirocumab:PCSK9 (1:1 complex)", units = "nmol", specimen = "serum", verified = TRUE),
    complex_peripheral1 = list(
      analyte = "alirocumab:PCSK9 (1:1 complex)",
      units = "nmol",
      specimen = "not applicable",
      verified = TRUE
    ),
    complex_db = list(
      analyte = "alirocumab:PCSK9:PCSK9 (1:2 complex)",
      units = "nmol",
      specimen = "serum",
      verified = TRUE
    ),
    complex_db_peripheral1 = list(
      analyte = "alirocumab:PCSK9:PCSK9 (1:2 complex)",
      units = "nmol",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_studies = 2L,
    disease_state = "Healthy adults",
    dose_range = "Single IV and single SC doses of alirocumab with dense PK sampling",
    studies = "Phase I NCT01026597 (R727-CL-0902) and NCT01074372 (R727-CL-0904)",
    co_medication = "No high-dose statins (healthy volunteers), so no statin covariate was tested",
    notes = "Kovalenko 2025 reports the adult parameter estimates (Table 2, 'Adults' column) but not the adult subject count or demographics; the studies are described in its references 7-9. The adult fit supplied the fixed absorption, distribution and bioavailability values of the pediatric model."
  )

  ini({
    # ---- Absorption and distribution (Table 2, adults) ----
    lktr <- log(22.5); label("Transit rate constant depot -> transit1, ktr (1/day)") # Table 2 'ktr' adults 22.5
    lka <- log(0.333); label("Absorption rate constant transit1 -> central, ka (1/day)") # Table 2 'ka' adults 0.333
    lfdepot <- log(0.655); label("SC bioavailability F (fraction)") # Table 2 'F' adults 0.655
    lk12 <- log(0.194); label("Central-to-peripheral rate constant kcp (1/day)") # Table 2 'kcp' adults 0.194
    lk21 <- log(0.280); label("Peripheral-to-central rate constant kpc (1/day)") # Table 2 'kpc' adults 0.280

    # ---- Linear PK and PCSK9 turnover (Table 2, adults) ----
    lvc <- log(3.01); label("Central volume Vc at 75 kg (L)") # Table 2 'Vc (L)' adults 3.01
    lkel <- log(0.0698); label("Elimination rate constant of free alirocumab keL (1/day)") # Table 2 'keL (1/day)' adults 0.0698
    lkint <- log(0.320); label("Elimination rate constant of alirocumab:PCSK9 complexes keLR (1/day)") # Table 2 'keLR (1/day)' adults 0.320
    lksyn <- log(4.97); label("PCSK9 synthesis rate ksyn (nM/day)") # Table 2 'ksyn (nM/day)' adults 4.97
    lr0 <- log(2.28); label("Baseline free PCSK9 concentration R0 (nM)") # Table 2 'R0 (nM)' adults 2.28

    # ---- Binding: fixed to Biacore in-vivo values (Methods; Table S11) ----
    lkon <- fixed(log(57.9)); label("Association rate constant kon (1/nM/day)") # Table S11 'kon (1/day/nM)' 57.9 (fixed, Biacore)
    lkoff <- fixed(log(106)); label("Dissociation rate constant koff (1/day)") # Table S11 'koff (1/day)' 106 (fixed, Biacore)

    # ---- Covariate effect (Table 2) ----
    e_wt_vc <- 0.787; label("Power exponent of body weight (WT / 75 kg) on Vc (unitless)") # Table 2 'Vc ~ weight' adults 0.787

    # ---- Between-subject variability (Table 2 SD of log-parameter, squared) ----
    etalvc ~ 0.058081 # Table 2 'SD(ln[Vc])' adults 0.241 -> variance 0.241^2
    etalkint ~ 0.073441 # Table 2 'SD(ln[keLR])' adults 0.271 -> variance 0.271^2
    etalksyn ~ 0.125316 # Table 2 'SD(ln[ksyn])' adults 0.354 -> variance 0.354^2
    etalr0 ~ 0.045796 # Table 2 'SD(ln[R0])' adults 0.214 -> variance 0.214^2
    etalkel ~ 0.058564 # Table 2 'SD(ln[keL])' adults 0.242 -> variance 0.242^2
    etalka ~ 0.334084 # Table 2 'SD(ln[ka])' adults 0.578 -> variance 0.578^2
    etalktr ~ 1.0201 # Table 2 'SD(ln[ktr])' adults 1.01 -> variance 1.01^2

    # ---- Residual error: proportional per analyte (Table 2) ----
    propSd <- 0.147; label("Proportional residual error, total alirocumab (fraction)") # Table 2 'SD(alirocumab)' adults 0.147
    propSd_Ctarget <- 0.249; label("Proportional residual error, free PCSK9 (fraction)") # Table 2 'SD(PCSK9)' adults 0.249
  })

  model({
    # Alirocumab molecular weight. Kovalenko 2025 states parameters in nM
    # without giving the mass conversion; approximately 146 kDa is the value
    # on the Praluent (alirocumab) label.
    mw_alirocumab <- 146000 # g/mol
    nmol_per_mg <- 1e6 / mw_alirocumab

    # ---- Individual parameters ----
    ktr <- exp(lktr + etalktr)
    ka <- exp(lka + etalka)
    fdepot <- exp(lfdepot)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    vc <- exp(lvc + etalvc) * (WT / 75)^e_wt_vc
    kel <- exp(lkel + etalkel)
    kint <- exp(lkint + etalkint)
    ksyn <- exp(lksyn + etalksyn)
    r0 <- exp(lr0 + etalr0)
    kon <- exp(lkon)
    koff <- exp(lkoff)
    # Methods: keR = ksyn / R0, so free PCSK9 starts at its steady state
    kdeg <- ksyn / r0

    # ---- Binding fluxes (nmol/day) in the central compartment (Figure 1) ----
    # L + R <-> L:R and L:R + R <-> L:R:R, both with kon and koff.
    bind1 <- kon * central * target / vc - koff * complex
    bind2 <- kon * complex * target / vc - koff * complex_db

    # ---- ODEs: amounts in nmol; every species shares kcp / kpc (Figure 1) ----
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ka * transit1
    d/dt(central) <- ka * transit1 - kel * central - k12 * central + k21 * peripheral1 - bind1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(target) <- ksyn * vc - kdeg * target - k12 * target + k21 * target_peripheral1 - bind1 - bind2
    d/dt(target_peripheral1) <- k12 * target - k21 * target_peripheral1
    d/dt(complex) <- bind1 - bind2 - kint * complex - k12 * complex + k21 * complex_peripheral1
    d/dt(complex_peripheral1) <- k12 * complex - k21 * complex_peripheral1
    d/dt(complex_db) <- bind2 - kint * complex_db - k12 * complex_db + k21 * complex_db_peripheral1
    d/dt(complex_db_peripheral1) <- k12 * complex_db - k21 * complex_db_peripheral1

    # Pre-dose steady state: R0 in the central compartment and the matching
    # distribution equilibrium in the peripheral compartment.
    target(0) <- r0 * vc
    target_peripheral1(0) <- r0 * vc * k12 / k21

    # Doses are in mg; SC doses enter the depot (bioavailability F), IV doses
    # the central compartment. Both are converted to nmol here, so give IV
    # infusions with `dur` rather than `rate`.
    f(depot) <- fdepot * nmol_per_mg
    f(central) <- nmol_per_mg

    # ---- Observations ----
    # Total alirocumab assay measures free plus 1- and 2-PCSK9-bound drug
    # (Analytical Methods); reported in mg/L. Free PCSK9 is reported in nM.
    Cc <- (central + complex + complex_db) / vc / nmol_per_mg
    Ctarget <- target / vc
    Cc ~ prop(propSd)
    Ctarget ~ prop(propSd_Ctarget)
  })
}
