Oliveira_2026_bictegravir_rabbit <- function() {
  description <- "Preclinical (rabbit). Two-compartment PK model for plasma bictegravir in New Zealand White rabbits after a single 0.75 mg/kg IV bolus, fitted in Phoenix NLME in the macro-constant (A, alpha, B, beta) parameterisation and used as the unit impulse response for deconvolving absorption from a subcutaneous reservoir-style poly(caprolactone) bictegravir implant. The implant input is the published single-Weibull model of the median fraction absorbed over one year (Finf, MDT, b), applied to the loaded implant mass. Dosing is per kg of body weight."
  reference <- paste(
    "Oliveira D, Cruz D, Levin LPK, Li L, Norton CA, Dobek G, Wang X,",
    "Veazey R, Watkins M, Schauer AP, Dumond JB, Johnson LM, Cottrell ML.",
    "Harnessing Pharmacokinetic Modeling to Develop a Long-Acting Subcutaneous",
    "HIV Treatment Platform for Young Children. Pharmaceutics.",
    "2026;18(5):522. doi:10.3390/pharmaceutics18050522",
    sep = " "
  )
  vignette <- "Oliveira_2026_antiretroviral_implant_rabbit"
  # Doses were given per kg of body weight (0.75 mg/kg IV), and no rabbit body
  # weights are reported, so the model is written per kg: dose in mg/kg and
  # volume in L/kg. An implant loaded with M mg in a rabbit weighing W kg is a
  # depot dose of M / W mg/kg.
  units <- list(time = "h", dosing = "mg/kg", concentration = "ng/mL")
  # IV bolus doses go to central; an implant dose goes to depot.
  dosing <- c("depot", "central")

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "bictegravir", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "bictegravir", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "bictegravir", units = "mg/kg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "rabbit (New Zealand White)",
    n_subjects = 16,
    n_studies = 1,
    sex_female_pct = 50,
    disease_state = "healthy rabbits",
    dose_range = paste(
      "0.75 mg/kg single IV bolus (marginal ear vein); after a 7-week washout,",
      "subcutaneous bictegravir-eluting implants (mean 165 mg loaded per",
      "implant) placed between the scapulae and followed for up to one year"
    ),
    regions = "USA (Tulane University vivarium)",
    notes = paste(
      "Oliveira 2026 Methods 'In Vivo Absorption': 24 rabbits, 8 per IV group",
      "(BIC 0.75 mg/kg only; ISL/FTC 5/30 mg/kg; BIC/ISL/FTC 0.75/5/30 mg/kg),",
      "so 16 rabbits received bictegravir IV. Sparse design: each rabbit gave",
      "2 samples between 0.5 and 30 h (4 observations per group per time).",
      "Each implant group had equal numbers of males and females. Body",
      "weights are not reported. The Weibull absorption model was fitted to",
      "the median deconvolved fraction absorbed of the rabbits followed for",
      "the whole year (Regression Analyses; Table 4 gives n = 4 at 6 and 12",
      "months)."
    )
  )

  ini({
    # Phoenix NLME two-compartment macro model (Table 3, 'Final Estimate'
    # column; the %CV printed beside each value is its relative standard
    # error). The coefficients A and B are concentrations (ng/mL) after the
    # 0.75 mg/kg IV dose, so the dose-normalised impulse response is
    # (A / Dref) exp(-alpha t) + (B / Dref) exp(-beta t) with Dref = 0.75
    # mg/kg. It is written with the two-compartment identities
    #   vc  = Dref / (A + B)
    #   k21 = (A * beta + B * alpha) / (A + B)
    #   kel = alpha * beta / k21
    #   k12 = alpha + beta - k21 - kel
    # with A = 1115.6 ng/mL and B = 102.0 ng/mL. Dref in ug/kg over
    # (A + B) in ug/L gives vc in L/kg.
    lalpha <- log(0.438); label("Distribution-phase hybrid rate constant alpha (1/h)") # Table 3 BIC Alpha = 0.438 1/h (%CV 20)
    lbeta <- log(0.170); label("Terminal-phase hybrid rate constant beta (1/h)") # Table 3 BIC Beta = 0.170 1/h (%CV 23)
    lvc <- log(750 / (1115.6 + 102.0)); label("Central volume of distribution per kg, Dref / (A + B) (L/kg)") # derived: Dref = 0.75 mg/kg (Methods); Table 3 BIC A = 1115.6 ng/mL, B = 102.0 ng/mL; 750 / 1217.6 = 0.616 L/kg
    lk21 <- log((1115.6 * 0.170 + 102.0 * 0.438) / (1115.6 + 102.0)); label("Peripheral to central rate constant k21, (A beta + B alpha) / (A + B) (1/h)") # derived: Table 3 BIC A, B, Alpha, Beta; 234.33 / 1217.6 = 0.192 1/h

    # Implant input: Phoenix single-Weibull dissolution model (Model 602)
    # fitted to the median fraction absorbed over time (Supplementary Table
    # S1, BIC 'Weibull', the best fit, r2 = 0.9999, AIC = -981):
    #   F(t) = Finf * (1 - exp(-(t / MDT)^b)),  t in hours, TLAG = 0.
    # Written with the registered Weibull stems ra = 1 / MDT and gam1 = b.
    # Finf is the asymptotic fraction of the loaded mass absorbed and is
    # greater than 1 here: the fit only spans the first year, over which 35%
    # of the load is absorbed (Results: 65% remains at 12 months).
    lfdepot <- log(1.666046); label("Asymptotic fraction of the loaded implant mass absorbed, Finf (unitless)") # Table S1 BIC Weibull FINF = 1.666046 (SE 0.178462)
    lra <- log(1 / 39994.39); label("Weibull rate-scaling parameter of implant absorption, ra = 1/MDT (1/h)") # Table S1 BIC Weibull MDT = 39994.39 h (SE 5537.618); ra = 1/MDT
    lgam1 <- log(0.949174); label("Weibull shape of implant absorption, b (unitless)") # Table S1 BIC Weibull B = 0.949174 (SE 0.007761)

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

    # 3. Implant absorption. The depot empties with the hazard of the
    #    Weibull survival function S(t) = exp(-(ra * t)^gam1),
    #    h(t) = gam1 * ra * (ra * t)^(gam1 - 1), so a depot dosed with the
    #    loaded mass and scaled by f(depot) = Finf releases exactly
    #    Dose * dF/dt. Time after the implant dose is floored at a small
    #    positive number because (ra * t)^(gam1 - 1) is not finite at t = 0
    #    for a shape below 1, and capped where the cumulative hazard reaches
    #    40 (depot below exp(-40) of its dose).
    fdepot <- exp(lfdepot)
    ra <- exp(lra)
    gam1 <- exp(lgam1)
    tr1 <- min(max(tad0(depot), 1e-6), 40^(1 / gam1) / ra)
    h1 <- gam1 * ra * (ra * tr1)^(gam1 - 1)

    # 4. ODE system
    d/dt(depot) <- -h1 * depot
    d/dt(central) <- h1 * depot - (kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    f(depot) <- fdepot

    # 5. Observation: mg/kg over L/kg is mg/L; x 1000 gives ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ lnorm(expSd)
  })
}
