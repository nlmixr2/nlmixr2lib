Kim_2024_niclosamide_hamster <- function() {
  description <- paste(
    "Preclinical (hamster). Final power-law population PK model for a",
    "long-acting intramuscular niclosamide formulation in male Golden",
    "Syrian hamsters, fit jointly to plasma and lung-tissue concentrations.",
    "The dose is split between a fast first-order depot (fraction 1 - F,",
    "rate ka) and a slow depot (fraction F) whose first-order rate is",
    "switched on by a Weibull cumulative-distribution function of time,",
    "ka2(t) = kamax * (1 - exp(-(t / lambda)^gamma)). Disposition is two",
    "compartments (central, peripheral1) with linear clearance plus a",
    "perfusion-limited lung compartment whose volume and blood flow were",
    "held at allometrically interpolated hamster values. Exchange with the",
    "lung follows the power law dAlung/dt = Qlung * (Cplasma^Xcl - Clung)",
    "(Kopelman anomalous reaction order Xcl), so lung concentration at",
    "equilibrium is Cplasma^Xcl rather than proportional to plasma. The",
    "power law is tied to the paper's concentration units: Cplasma and",
    "Clung are in ng/mL. A linear partition-coefficient alternative is",
    "Kim_2024_niclosamide_hamster_linear.",
    sep = " "
  )
  reference <- paste(
    "Kim T, Jung W, Cho S, Kim G, Yun H-y, Chae J-w. Predicting lung",
    "exposure of intramuscular niclosamide as an antiviral agent:",
    "Power-law based pharmacokinetic modeling. Clin Transl Sci.",
    "2024;17(5):e13833. doi:10.1111/cts.13833. PMCID: PMC11128490.",
    "Parameter estimates: Table 2. Structure: Figure 1 and Equation 2.",
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
    # Disposition: Table 2 (power-law model), final estimate (RSE%).
    lcl <- log(0.485); label("Clearance (L/h)") # Table 2: CL = 0.485 L/h (RSE 7%)
    lvc <- log(0.371); label("Central volume of distribution (L)") # Table 2: Vc = 0.371 L (RSE 12.3%)
    lq <- log(0.365); label("Intercompartmental clearance, central to peripheral (L/h)") # Table 2: Q inter = 0.365 L/h (RSE 14.0%)
    lvp <- log(1.647); label("Peripheral volume of distribution (L)") # Table 2: Vp = 1.647 L (RSE 18.4%)

    # Physiological lung parameters, not estimated (Table 2 asterisk): from
    # allometric interpolation across species (Methods, PK Model development).
    lq_lung <- fixed(log(0.0613)); label("Lung blood flow (L/h)") # Table 2: Q lung = 0.0613 L/h, asterisk
    lv_lung <- fixed(log(0.000919)); label("Lung volume (L)") # Table 2: V l = 0.000919 L, asterisk
    lpow_lung <- log(0.753); label("Anomalous reaction order of plasma concentration in the plasma-to-lung exchange (unitless)") # Table 2: X cl = 0.753 (RSE 1%); Equation 2

    # Absorption: Figure 1 (fast depot 1 - F at Ka1; slow depot F at Ka2 * WB,
    # WB = 1 - exp(-(time / lambda)^gamma)).
    lka <- log(0.042); label("Fast first-order absorption rate constant (1/h)") # Table 2: K A1 = 0.042 1/h (RSE 2.89%)
    lkamax <- log(0.002); label("Slow absorption rate constant reached once the Weibull switch is fully on (1/h)") # Table 2: K A2 = 0.002 1/h (RSE 16.0%)
    lra <- log(1 / 155.80); label("Weibull time-scale of the slow-absorption switch, 1/lambda (1/h)") # Table 2: Lambda = 155.80 h (RSE 4.93%); ra = 1/lambda
    lgam1 <- log(7.199); label("Weibull shape of the slow-absorption switch (unitless)") # Table 2: Gamma = 7.199 (RSE 5.70%)
    lfdepot2 <- log(0.430); label("Fraction of the dose entering the slow (Weibull) depot (fraction)") # Table 2: F = 0.430 (RSE 39%)

    # IIV: Table 2 'Interindividual variability', read as exponential-model
    # omega^2 variances.
    etalvc ~ 0.299 # Table 2 IIV Vc = 0.299 (RSE 27.65%, shrinkage 21.67%)
    etalcl ~ 0.0233 # Table 2 IIV CL = 0.0233 (RSE 25.56%, shrinkage 19.94%)
    etalfdepot2 ~ 0.169 # Table 2 IIV F = 0.169 (RSE 19.55%, shrinkage 10.54%)

    # Residual error: observations were log-transformed (Methods) and the
    # Table 2 footnote calls the error proportional, i.e. additive on the log
    # scale; Table 2 values read as SDs.
    expSd <- 0.373; label("Log-scale residual SD, plasma (log ng/mL)") # Table 2: Residuals Plasma = 0.373 (RSE 4%)
    expSd_Clung <- 0.361; label("Log-scale residual SD, lung tissue (log ng/mL)") # Table 2: Residuals Lung = 0.361 (RSE 9%)
  })

  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)
    q_lung <- exp(lq_lung)
    v_lung <- exp(lv_lung)
    pow_lung <- exp(lpow_lung)
    ka <- exp(lka)
    kamax <- exp(lkamax)
    ra <- exp(lra)
    gam1 <- exp(lgam1)
    # Exponential IIV on a fraction can exceed 1 (about 2% of animals at
    # omega^2 = 0.169); capped at 1 so the fast-depot fraction never goes
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

    # Equation 2: dAlung/dt = Qlung * (Cplasma^Xcl - Clung), with the
    # concentrations in ng/mL; divided by 1000 to give mg/h. The same flux
    # leaves central (Figure 1 two-way central-lung arrow). Cc is floored at 0
    # inside the power because a non-integer power of an integrator undershoot
    # below zero is NaN (plasma sits near 0 for days when all of the dose is in
    # the slow depot).
    flux_lung <- q_lung * (max(Cc, 0)^pow_lung - Clung) / 1000

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
