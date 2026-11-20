Oliveira_2026_islatravir_rabbit <- function() {
  description <- "Preclinical (rabbit). Two-compartment PK model for plasma islatravir in New Zealand White rabbits after a single 5 mg/kg IV bolus, fitted in Phoenix NLME in the macro-constant (A, alpha, B, beta) parameterisation and used as the unit impulse response for deconvolving absorption from a subcutaneous reservoir-style poly(caprolactone) islatravir implant. The implant input is the published double-Weibull model of the median fraction absorbed over one year (Finf, F1, MDT1, b1, MDT2, b2), applied to the loaded implant mass through two parallel release depots. Dosing is per kg of body weight."
  reference <- paste(
    "Oliveira D, Cruz D, Levin LPK, Li L, Norton CA, Dobek G, Wang X,",
    "Veazey R, Watkins M, Schauer AP, Dumond JB, Johnson LM, Cottrell ML.",
    "Harnessing Pharmacokinetic Modeling to Develop a Long-Acting Subcutaneous",
    "HIV Treatment Platform for Young Children. Pharmaceutics.",
    "2026;18(5):522. doi:10.3390/pharmaceutics18050522",
    sep = " "
  )
  vignette <- "Oliveira_2026_antiretroviral_implant_rabbit"
  # Doses were given per kg of body weight (5 mg/kg IV), and no rabbit body
  # weights are reported, so the model is written per kg: dose in mg/kg and
  # volume in L/kg. An implant loaded with M mg in a rabbit weighing W kg is
  # recorded as M / W mg/kg on BOTH depot and depot2 (the release fractions
  # are applied as bioavailabilities of the two parallel depots).
  units <- list(time = "h", dosing = "mg/kg", concentration = "ng/mL")
  dosing <- c("depot", "depot2", "central")

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "islatravir", units = "mg/kg", specimen = "administration site", verified = TRUE),
    depot2 = list(analyte = "islatravir", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "islatravir", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "islatravir", units = "mg/kg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "rabbit (New Zealand White)",
    n_subjects = 16,
    n_studies = 1,
    sex_female_pct = 50,
    disease_state = "healthy rabbits",
    dose_range = paste(
      "5 mg/kg single IV bolus (marginal ear vein, with emtricitabine 30",
      "mg/kg, with or without bictegravir 0.75 mg/kg); after a 7-week washout,",
      "one subcutaneous islatravir-eluting implant (mean 62 mg loaded) placed",
      "between the scapulae and followed for one year"
    ),
    regions = "USA (Tulane University vivarium)",
    notes = paste(
      "Oliveira 2026 Methods 'In Vivo Absorption': 16 rabbits received",
      "islatravir IV (ISL/FTC and BIC/ISL/FTC groups, 8 each). Sparse design:",
      "each rabbit gave 2 samples between 0.5 and 30 h. Each implant group had",
      "equal numbers of males and females. Body weights are not reported.",
      "Absorption was deconvolved for 16 implanted rabbits; the double-Weibull",
      "model was fitted to the median fraction absorbed of the rabbits",
      "followed for the whole year (Table 4 gives n = 8 at 6 and 12 months)."
    )
  )

  ini({
    # Phoenix NLME two-compartment macro model (Table 3; the %CV printed
    # beside each value is its relative standard error). A and B are
    # concentrations (ng/mL) after the 5 mg/kg IV dose, so the dose-normalised
    # impulse response is (A / Dref) exp(-alpha t) + (B / Dref) exp(-beta t)
    # with Dref = 5 mg/kg, written with the two-compartment identities
    #   vc  = Dref / (A + B)
    #   k21 = (A * beta + B * alpha) / (A + B)
    #   kel = alpha * beta / k21
    #   k12 = alpha + beta - k21 - kel
    #
    # Table 3 prints the ISL final Alpha and Beta to one decimal place only
    # (0.5 and 0.1). The bootstrap medians in the same table, 0.45 and 0.06,
    # lie inside those rounding intervals and, with the final A and B,
    # reproduce the mean IV profile of Figure 3B, which the one-decimal values
    # miss by up to 1.2 log units at 30 h; the vignette shows the comparison.
    # The bootstrap medians are used for the two rate constants.
    lalpha <- log(0.45); label("Distribution-phase hybrid rate constant alpha (1/h)") # Table 3 ISL Alpha bootstrap estimate = 0.45 1/h (95% CI 0.41-0.49); final estimate printed as 0.5 (%CV 6)
    lbeta <- log(0.06); label("Terminal-phase hybrid rate constant beta (1/h)") # Table 3 ISL Beta bootstrap estimate = 0.06 1/h (95% CI 0.04-0.09); final estimate printed as 0.1 (%CV 40)
    lvc <- log(5000 / (5275.9 + 6.3)); label("Central volume of distribution per kg, Dref / (A + B) (L/kg)") # derived: Dref = 5 mg/kg (Methods); Table 3 ISL A = 5275.9 ng/mL, B = 6.3 ng/mL; 5000 / 5282.2 = 0.947 L/kg
    lk21 <- log((5275.9 * 0.06 + 6.3 * 0.45) / (5275.9 + 6.3)); label("Peripheral to central rate constant k21, (A beta + B alpha) / (A + B) (1/h)") # derived: Table 3 ISL A, B (final) and Alpha, Beta (bootstrap); 319.39 / 5282.2 = 0.0605 1/h

    # Implant input: Phoenix double-Weibull dissolution model (Model 603)
    # fitted to the median fraction absorbed over time (Supplementary Table
    # S1, ISL 'Double Weibull', the best fit, r2 = 0.9999, AIC = -891):
    #   F(t) = Finf * (F1 * (1 - exp(-(t / MDT1)^b1))
    #                  + (1 - F1) * (1 - exp(-(t / MDT2)^b2))),  t in hours.
    # Written with the registered stems ra = 1 / MDT1, gam1 = b1,
    # ra2 = 1 / MDT2, gam2 = b2 and frel = F1.
    lfdepot <- log(0.544194); label("Asymptotic fraction of the loaded implant mass absorbed, Finf (unitless)") # Table S1 ISL Double Weibull FINF = 0.544194 (SE 0.000568)
    logitfrel <- log(0.656434 / (1 - 0.656434)); label("Logit of the fraction of Finf absorbed by the first Weibull process, F1 (unitless)") # Table S1 ISL Double Weibull F1 = 0.656434 (SE 0.040172)
    lra <- log(1 / 2288.246); label("Weibull rate-scaling parameter of the first process, ra = 1/MDT1 (1/h)") # Table S1 ISL Double Weibull MDT1 = 2288.246 h (SE 95.21243)
    lgam1 <- log(1.426693); label("Weibull shape of the first process, b1 (unitless)") # Table S1 ISL Double Weibull B1 = 1.426693 (SE 0.066116)
    lra2 <- log(1 / 589.4632); label("Weibull rate-scaling parameter of the second process, ra2 = 1/MDT2 (1/h)") # Table S1 ISL Double Weibull MDT2 = 589.4632 h (SE 28.99826)
    lgam2 <- log(1.400658); label("Weibull shape of the second process, b2 (unitless)") # Table S1 ISL Double Weibull B2 = 1.400658 (SE 0.04603)

    # Results: 'BIC and ISL included random effects on all four parameters'
    # (A, Alpha, B, Beta), but no variance is reported. The etas on alpha and
    # beta are kept at zero; the etas on A and B have no one-to-one
    # counterpart on lvc and lk21 and are not carried.
    etalalpha ~ fixed(0) # Results: random effect on Alpha, variance not reported
    etalbeta ~ fixed(0) # Results: random effect on Beta, variance not reported

    # Results: 'log-additive error'; magnitude not reported.
    expSd <- fixed(0); label("Log-scale residual error SD (not reported)") # Results: log-additive error, value not reported
  })

  model({
    # 1. Disposition
    alpha <- exp(lalpha + etalalpha)
    beta <- exp(lbeta + etalbeta)
    vc <- exp(lvc)
    k21 <- exp(lk21)

    # 2. Micro-constants from alpha * beta = kel * k21 and
    #    alpha + beta = kel + k12 + k21
    kel <- alpha * beta / k21
    k12 <- alpha + beta - k21 - kel

    # Clearance-parameterisation equivalents (derived outputs, per kg)
    cl <- kel * vc
    q <- k12 * vc
    vp <- q / k21

    # 3. Implant absorption. Each depot empties with the hazard of its own
    #    Weibull survival function S(t) = exp(-(ra * t)^gam),
    #    h(t) = gam * ra * (ra * t)^(gam - 1); with the bioavailabilities
    #    Finf * F1 and Finf * (1 - F1) the two depots together release
    #    Dose * dF/dt. Time after the implant dose is floored at a small
    #    positive number and capped where the cumulative hazard reaches 40
    #    (depot below exp(-40) of its dose), beyond which a shape above 1
    #    would keep growing the hazard of an empty depot.
    fdepot <- exp(lfdepot)
    frel <- expit(logitfrel)
    ra <- exp(lra)
    gam1 <- exp(lgam1)
    ra2 <- exp(lra2)
    gam2 <- exp(lgam2)
    tr1 <- min(max(tad0(depot), 1e-6), 40^(1 / gam1) / ra)
    tr2 <- min(max(tad0(depot2), 1e-6), 40^(1 / gam2) / ra2)
    h1 <- gam1 * ra * (ra * tr1)^(gam1 - 1)
    h2 <- gam2 * ra2 * (ra2 * tr2)^(gam2 - 1)

    # 4. ODE system
    d/dt(depot) <- -h1 * depot
    d/dt(depot2) <- -h2 * depot2
    d/dt(central) <- h1 * depot + h2 * depot2 - (kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    f(depot) <- fdepot * frel
    f(depot2) <- fdepot * (1 - frel)

    # 5. Observation: mg/kg over L/kg is mg/L; x 1000 gives ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ lnorm(expSd)
  })
}
