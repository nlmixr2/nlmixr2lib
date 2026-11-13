Duong_2026_cobimetinib <- function() {
  description <- "QSP. Single-tumour-cell stochastic PK-PD model of cobimetinib (Duong 2026): one-compartment first-order-absorption plasma PK for a representative patient on 60 mg once daily 21/7, free plasma drug driving uptake into the cell at a rate proportional to the cell's influx-transporter count (log-normal across cells, geometric mean 500, sigma 0.5), first-order intracellular metabolism, and a sigmoid-Emax intracellular-concentration-driven death hazard integrated into a cumulative hazard. Each simulated subject is one tumour cell; the per-cell random effect is the transporter heterogeneity. All parameter values are fixed literature/calibration values, not fitted."
  reference <- "Duong TT, et al. Cellular Heterogeneity in Drug Uptake Amplifies Pharmacodynamic Variability: A Stochastic PK-PD Analysis. CPT Pharmacometrics Syst Pharmacol. 2026;15:e70282. doi:10.1002/psp4.70282"
  vignette <- "Duong_2026_cobimetinib"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "mg/L (the source labels plasma, free and intracellular concentrations uM, but its code computes them as dose in mg over V/F in L with no molecular-weight conversion; EC50 = 2.0 is on that same numerical scale)"
  )

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "cobimetinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "cobimetinib", units = "mg", specimen = "plasma", verified = TRUE),
    int_tumor = list(
      analyte = "cobimetinib",
      units = "mg/L",
      specimen = "tumor",
      verified = TRUE,
      notes = "Intracellular drug concentration of a single tumour cell (the paper's C_i); a driven state that does not deplete plasma."
    ),
    cumhaz = list(
      analyte = "hazard",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE,
      notes = "Cumulative death hazard H_i(t) of the cell; the cell dies when H_i reaches -log(U_i), U_i ~ Uniform(0, 1)."
    )
  )

  population <- list(
    species = "human (in silico)",
    n_subjects = NA_integer_,
    n_studies = 0L,
    disease_state = "Hypothetical solid tumour treated with cobimetinib; proof-of-principle simulation study, no patient data were fitted.",
    dose_range = "60 mg orally once daily for 21 days followed by 7 days off (21/7), repeated; sensitivity analyses at 0.25x, 0.5x, 1x and 2x that dose.",
    notes = "Each simulated subject is one tumour cell of a representative patient. Plasma PK is deterministic and identical for all cells; cells differ only in influx-transporter count T_i ~ LogNormal(log(500), 0.5^2) (Methods, Cellular Pharmacokinetics). Populations of N = 100-2000 cells were simulated (Experiments 2 and 4)."
  )

  ini({
    # Systemic PK (Methods, Systemic Pharmacokinetics; supplement pk.jl)
    lka <- fixed(log(2.0)); label("Absorption rate constant ka (1/h)") # Methods: 'ka=2.0 h-1'
    lvc <- fixed(log(806)); label("Apparent volume of distribution V/F (L)") # Methods: 'V/F=806 L'
    lkel <- fixed(log(0.0171)); label("Elimination rate constant ke (1/h)") # Methods: 'ke=0.0171 h-1'
    fu <- fixed(0.05); label("Fraction unbound in plasma (unitless)") # Methods: '95% protein binding (fu=0.05)'

    # Cellular PK (Methods, Cellular Pharmacokinetics; supplement cellular_dynamics.jl)
    # Uptake rate of the median cell = k_in x T_median = 0.038 1/h per transporter x 500 transporters.
    lkin_int_tumor <- fixed(log(0.038 * 500)); label("Uptake rate constant of the median cell, k_in * 500 transporters (1/h)") # Methods: 'kin=0.038 h-1' per transporter; 'muT=log500'
    lkout_int_tumor <- fixed(log(0.1)); label("First-order intracellular metabolic clearance rate k_met (1/h)") # Methods: 'kmet=0.1 h-1'

    # Cellular PD (Methods, Cellular Pharmacodynamics and Death Process)
    lemax <- fixed(log(1.5)); label("Maximum death rate mu_max (1/day)") # Methods: 'mumax=1.5 day-1'
    lec50 <- fixed(log(2.0)); label("Intracellular concentration at half-maximal death rate EC50 (mg/L; uM in the source)") # Methods: 'EC50=2.0 uM'
    lhill <- fixed(log(2.5)); label("Hill coefficient n (unitless)") # Methods: 'n=2.5'

    # Cell-to-cell transporter heterogeneity: log(T_i) ~ N(log(500), sigma_T^2), sigma_T = 0.5,
    # so the per-cell uptake rate k_in * T_i carries an eta of variance 0.25.
    etalkin_int_tumor ~ fixed(0.25) # Methods: 'sigmaT=0.5 yields a coefficient of variation of 53%'
  })

  model({
    # Systemic one-compartment PK with first-order absorption
    ka <- exp(lka)
    vc <- exp(lvc)
    kel <- exp(lkel)

    # Cellular parameters; exp(etalkin_int_tumor) = T_i / 500
    kin_int_tumor <- exp(lkin_int_tumor + etalkin_int_tumor)
    kout_int_tumor <- exp(lkout_int_tumor)
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    Cc <- central / vc
    # Free plasma concentration drives transmembrane uptake
    Cu <- fu * Cc

    # Deterministic drift of the paper's SDE dC_i = (k_in T_i C_p,u - k_met C_i) dt + sigma_C dW_i;
    # the sigma_C = 0.01 diffusion term is omitted (see vignette).
    d/dt(int_tumor) <- kin_int_tumor * Cu - kout_int_tumor * int_tumor

    # Sigmoid-Emax death rate in 1/day; the hazard integrates on the hour scale
    hazard <- emax * int_tumor^hill / (ec50^hill + int_tumor^hill)
    d/dt(cumhaz) <- hazard / 24

    # Probability that the cell is still alive at time t
    sur <- exp(-cumhaz)
  })
}
