Kim_2024_niclosamide_hamster_linear <- function() {
  description <- paste(
    "Preclinical (hamster). Linear partition-coefficient population PK",
    "model for a long-acting intramuscular niclosamide formulation in male",
    "Golden Syrian hamsters, fit jointly to plasma and lung-tissue",
    "concentrations; the alternative the authors compared against their",
    "final power-law model (Kim_2024_niclosamide_hamster), which it shares",
    "in structure except for the plasma-lung exchange. The dose is split",
    "between a fast first-order depot (fraction 1 - F, rate ka) and a slow",
    "depot (fraction F) whose first-order rate is switched on by a Weibull",
    "cumulative-distribution function of time, ka2(t) = kamax * (1 -",
    "exp(-(t / lambda)^gamma)). Disposition is two compartments (central,",
    "peripheral1) with linear clearance plus a perfusion-limited lung",
    "compartment whose volume and blood flow were held at allometrically",
    "interpolated hamster values, exchanging as dAlung/dt = Qlung *",
    "(Cplasma - Clung / Kp). The power-law model fit better (OFV -256.33",
    "versus -200.84).",
    sep = " "
  )
  reference <- paste(
    "Kim T, Jung W, Cho S, Kim G, Yun H-y, Chae J-w. Predicting lung",
    "exposure of intramuscular niclosamide as an antiviral agent:",
    "Power-law based pharmacokinetic modeling. Clin Transl Sci.",
    "2024;17(5):e13833. doi:10.1111/cts.13833. PMCID: PMC11128490.",
    "Parameter estimates: Appendix S1 Table S3. Structure: Figure 1 and",
    "Equation 1.",
    sep = " "
  )
  vignette <- "Kim_2024_niclosamide_hamster"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "niclosamide", units = "mg", specimen = "plasma", verified = TRUE),
    depot2 = list(analyte = "niclosamide", units = "mg", specimen = "plasma", verified = TRUE),
    central = list(analyte = "niclosamide", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "niclosamide", units = "mg", specimen = "plasma", verified = TRUE),
    lung = list(analyte = "niclosamide", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "hamster (male Golden Syrian, SLC Japan)",
    n_subjects = 25L,
    n_studies = 1L,
    age_range = "7 weeks at purchase, dosed after a 7-day acclimatisation",
    weight_range = "Not reported",
    sex_female_pct = 0,
    disease_state = "Healthy (uninfected) animals",
    dose_range = paste(
      "Single intramuscular injection of a long-acting niclosamide",
      "formulation at 55, 96, 128 or 240 mg/kg (plasma, n = 5 per dose);",
      "lung tissue from a separate group dosed at 240 mg/kg (n = 5 per",
      "sampling time)."
    ),
    regions = "Republic of Korea (author affiliations; the study site is not stated)",
    notes = paste(
      "Plasma and lung sampled at 0, 0.25, 0.5, 1, 2, 4, 8, 24, 48, 96, 168,",
      "240, 336, 504 and 672 h. NONMEM 7.5.1, FOCE-I, log-transformed",
      "observations. Lung volume and blood flow were interpolated from an",
      "allometric regression of log physiological value on log body weight",
      "across human, pig, monkey, rat and mouse. Body weights of the",
      "hamsters are not reported."
    )
  )

  ini({
    # Disposition: Appendix S1 Table S3 (linear model), final estimates,
    # printed to two decimals.
    lcl <- log(0.49); label("Clearance (L/h)") # Table S3: CL = 0.49 L/h
    lvc <- log(0.35); label("Central volume of distribution (L)") # Table S3: Vc = 0.35 L
    lq <- log(0.34); label("Intercompartmental clearance, central to peripheral (L/h)") # Table S3: Q inter = 0.34 L/h
    lvp <- log(1.58); label("Peripheral volume of distribution (L)") # Table S3: Vp = 1.58 L

    # Physiological lung parameters, not estimated (Table S3 asterisk).
    lq_lung <- fixed(log(0.0613)); label("Lung blood flow (L/h)") # Table S3: Q lung = 0.0613 L/h, asterisk
    lv_lung <- fixed(log(0.000919)); label("Lung volume (L)") # Table S3: Vl = 0.000919 L, asterisk
    lkp_lung <- log(0.27); label("Lung-to-plasma partition coefficient (unitless)") # Table S3: Kp = 0.27; Equation 1

    # Absorption: Figure 1.
    lka <- log(0.04); label("Fast first-order absorption rate constant (1/h)") # Table S3: KA1 = 0.04 1/h
    lkamax <- log(0.002); label("Slow absorption rate constant reached once the Weibull switch is fully on (1/h)") # Table S3: KA2 final printed as 0.00 at two decimals; bootstrap median 0.002 (0.001, 0.002) used
    lra <- log(1 / 164.60); label("Weibull time-scale of the slow-absorption switch, 1/lambda (1/h)") # Table S3: Lambda = 164.60 h; ra = 1/lambda
    lgam1 <- log(22.08); label("Weibull shape of the slow-absorption switch (unitless)") # Table S3: Gamma = 22.08
    lfdepot2 <- log(0.42); label("Fraction of the dose entering the slow (Weibull) depot (fraction)") # Table S3: F = 0.42

    # IIV: Table S3 'Inter-individual variability', read as exponential-model
    # omega^2 variances (same reading as Table 2 of the power-law model).
    etalvc ~ 0.53 # Table S3 IIV Vc = 0.53
    etalcl ~ 0.16 # Table S3 IIV CL = 0.16
    etalfdepot2 ~ 0.33 # Table S3 IIV F = 0.33

    # Residual error: log-transformed observations, proportional error
    # (Table S3 footnote), i.e. additive on the log scale; values read as SDs.
    expSd <- 0.39; label("Log-scale residual SD, plasma (log ng/mL)") # Table S3: Residuals Plasma = 0.39
    expSd_Clung <- 0.50; label("Log-scale residual SD, lung tissue (log ng/mL)") # Table S3: Residuals Lung = 0.50
  })

  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)
    q_lung <- exp(lq_lung)
    v_lung <- exp(lv_lung)
    kp_lung <- exp(lkp_lung)
    ka <- exp(lka)
    kamax <- exp(lkamax)
    ra <- exp(lra)
    gam1 <- exp(lgam1)
    # Exponential IIV on a fraction can exceed 1 (about 7% of animals at
    # omega^2 = 0.33); capped at 1 so the fast-depot fraction never goes
    # negative.
    fdepot2 <- min(exp(lfdepot2 + etalfdepot2), 1)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Weibull switch on the slow depot, a function of time since the dose
    # (Figure 1, WB = 1 - exp(-(time / lambda)^gamma)).
    ka2 <- kamax * (1 - exp(-(ra * t)^gam1))

    # Concentrations in ng/mL (= ug/L) from amounts in mg and volumes in L.
    Cc <- 1000 * central / vc
    Clung <- 1000 * lung / v_lung

    # Equation 1: dAlung/dt = Qlung * (Cplasma - Clung / Kp); divided by 1000
    # to give mg/h. The same flux leaves central (Figure 1 two-way arrow).
    flux_lung <- q_lung * (Cc - Clung / kp_lung) / 1000

    d/dt(depot) <- -ka * depot
    d/dt(depot2) <- -ka2 * depot2
    d/dt(central) <- ka * depot + ka2 * depot2 - kel * central - k12 * central + k21 * peripheral1 - flux_lung
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(lung) <- flux_lung

    # Each intramuscular dose is given as two records with the full amount,
    # one to depot and one to depot2; f() performs the split.
    f(depot) <- 1 - fdepot2
    f(depot2) <- fdepot2

    Cc ~ lnorm(expSd)
    Clung ~ lnorm(expSd_Clung)
  })
}
