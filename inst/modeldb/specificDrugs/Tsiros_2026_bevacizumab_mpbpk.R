Tsiros_2026_bevacizumab_mpbpk <- function() {
  description <- "mPBPK. Minimal physiologically-based PK model (Cao 2013 second-generation structure; plasma, tight-tissue ISF, leaky-tissue ISF and lymph) for TAB008, a bevacizumab biosimilar, in healthy Chinese adult men after a single 1 mg/kg 90-min IV infusion. Fitted by distributional physics-informed neural networks (D-PINNs) to digitised mean and SD plasma concentrations; plasma clearance carries log-normal interindividual variability, physiological volumes scale with body weight, and the residual error is log-normal with its SD fixed at 0.1."
  reference <- paste(
    "Tsiros P, Minadakis V, Sarimveis H. A physics-informed neural network",
    "approach for estimating population-level pharmacokinetic parameters from",
    "aggregated concentration data. J Pharmacokinet Pharmacodyn.",
    "2026;53(2):11. doi:10.1007/s10928-026-10019-w. PMID 41699348. PMCID",
    "PMC12909361. Aggregated data from Wang J et al. Front Pharmacol.",
    "2019;10:905. doi:10.3389/fphar.2019.00905 (TAB008 phase I). Model",
    "structure and fixed physiology from Cao Y, Balthasar JP, Jusko WJ.",
    "J Pharmacokinet Pharmacodyn. 2013;40(5):597-607.",
    "doi:10.1007/s10928-013-9332-2.",
    sep = " "
  )
  vignette <- "Tsiros_2026_bevacizumab_mpbpk"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  dosing <- "plasma"

  covariateData <- list(
    WT = list(
      description = "Body weight. Scales the blood, plasma, lymph and interstitial-fluid volumes (Tsiros 2026 Appendix C Eqs. 67-72) and, in simulation, the 1 mg/kg dose.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Tsiros 2026 sampled WT from Uniform(50, 75) kg, the range reported for the TAB008 phase I cohort, because individual weights were not available. No reference weight: volumes are proportional to WT.",
      source_name = "BW"
    )
  )

  compartmentData <- list(
    plasma = list(analyte = "TAB008 (bevacizumab biosimilar)", units = "mg", specimen = "plasma", verified = TRUE),
    tight = list(analyte = "TAB008 (bevacizumab biosimilar)", units = "mg", specimen = "tissue", verified = TRUE),
    leaky = list(analyte = "TAB008 (bevacizumab biosimilar)", units = "mg", specimen = "tissue", verified = TRUE),
    lymph = list(analyte = "TAB008 (bevacizumab biosimilar)", units = "mg", specimen = "lymph", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 49L,
    n_studies = 1L,
    age_range = "18-45 years",
    weight_range = "50-75 kg",
    sex_female_pct = 0,
    race_ethnicity = "Chinese",
    disease_state = "Healthy volunteers",
    dose_range = "1 mg/kg single IV infusion over 90 min",
    regions = "China",
    notes = "Tsiros 2026 Methods 'Monoclonal antibody data': the PK arm of the TAB008 phase I biosimilarity study (Wang 2019). The model was fitted to 15 mean and SD plasma concentrations digitised from the published figure; no individual data were used. Body weight was sampled from Uniform(50, 75) kg."
  )

  ini({
    lcl <- log(0.010832); label("Typical (median) plasma clearance CLp (L/h)") # Tsiros 2026 Table 3 D-PINNs: mu_CLp 0.0115, sigma_CLp 0.0041 (arithmetic mean and SD of a log-normal CLp); median = 0.0115 / sqrt(1 + (0.0041/0.0115)^2) = 0.010832
    sigma_tight <- 0.8378; label("Vascular reflection coefficient for tight tissues sigma1 (unitless)") # Tsiros 2026 Table 3 D-PINNs
    sigma_leaky <- 0.7529; label("Vascular reflection coefficient for leaky tissues sigma2 (unitless)") # Tsiros 2026 Table 3 D-PINNs

    sigmal <- fixed(0.2); label("Lymphatic capillary reflection coefficient sigmaL (unitless)") # Tsiros 2026 Appendix C text (from Cao 2013)
    kp <- fixed(0.8); label("Fraction of interstitial fluid available to an IgG1 antibody Kp (unitless)") # Tsiros 2026 Appendix C text (native IgG1, from Cao 2013)
    ltot <- fixed(0.125); label("Total lymph flow L (L/h)") # Tsiros 2026 Appendix C text: 3 L/day = 0.125 L/h
    hct <- fixed(0.43); label("Hematocrit (fraction)") # Tsiros 2026 Appendix C text: 0.43 for male volunteers

    etalcl ~ 0.119654 # Tsiros 2026 Table 3 D-PINNs: CV = 0.0041/0.0115 = 0.3565; omega^2 = log(1 + CV^2) = 0.119654

    expSd <- fixed(0.1); label("Log-normal residual error SD (log scale)") # Tsiros 2026 Results 'Case study 2': noise SD held at 0.1 (Eq. 14 log-normal error)
  })
  model({
    # Physiological volumes (Tsiros 2026 Appendix C Eqs. 67-72)
    vblood <- 0.073 * WT # Eq. 69
    visf <- 0.15 * WT # Eq. 70
    vplasma <- vblood * (1 - hct) # Eq. 71
    vlymph <- vblood # Eq. 72
    vtight <- 0.65 * visf * kp # Eq. 67
    vleaky <- 0.35 * visf * kp # Eq. 68

    # Lymph flow split (Eqs. 65-66); not weight-scaled in Appendix C
    l1 <- 0.33 * ltot
    l2 <- 0.67 * ltot

    cl <- exp(lcl + etalcl)

    cp <- plasma / vplasma
    ctight <- tight / vtight
    cleaky <- leaky / vleaky
    clymph <- lymph / vlymph

    # Eqs. 28-31 written in amounts (each concentration ODE multiplied by its volume)
    d/dt(plasma) <- clymph * ltot - cp * l1 * (1 - sigma_tight) - cp * l2 * (1 - sigma_leaky) - cl * cp
    d/dt(tight) <- l1 * (1 - sigma_tight) * cp - l1 * (1 - sigmal) * ctight
    d/dt(leaky) <- l2 * (1 - sigma_leaky) * cp - l2 * (1 - sigmal) * cleaky
    d/dt(lymph) <- l1 * (1 - sigmal) * ctight + l2 * (1 - sigmal) * cleaky - clymph * ltot

    Cc <- cp
    Cc ~ lnorm(expSd)
  })
}
