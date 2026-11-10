Protzenko_2025_busulfan <- function() {
  description <- "One-compartment IV-infusion population PK prior for high-dose once-daily busulfan in adults undergoing haematopoietic cell transplant conditioning, parameterised on the inverse volume of distribution (1/V) and the elimination rate constant k with correlated additive-normal between-subject variability, from a standard two-stage analysis of 16 richly sampled patients and used as the Bayesian prior for three-sample adaptive dosing (Protzenko 2025)"
  reference <- "Protzenko D, Bouchacourt B, Carriat L, Maroselli P, Boeri C, Devillier R, Ciccolini J. Adaptive dosing of high-dose busulfan in real-world adult patients undergoing haematopoietic cell transplant conditioning. Br J Clin Pharmacol. 2025;91(4):1171. doi:10.1111/bcp.16343"
  vignette <- "Protzenko_2025_busulfan"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The population distribution is Gaussian on the linear scale of the two
  # estimated parameters, 1/V and k (Methods: 'assuming Gaussian distribution
  # for interindividual variability'; standard two-stage mean vector and
  # covariance matrix). The etas are therefore additive on 1/V and on k, not
  # exponential on V and k.
  paper_specific_etas <- c("etavcinv", "etakel")

  covariateData <- list()

  compartmentData <- list(
    central = list(analyte = "busulfan", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 16L,
    n_studies = 1L,
    age_range = "adults (demographics of the 16 prior patients not reported)",
    disease_state = "Adults undergoing haematopoietic cell transplantation and treated with high-dose IV busulfan, enrolled in the HGB-207 protocol (EudraCT 2015-004122-33)",
    dose_range = "High-dose busulfan as a 3-h IV infusion (prior cohort doses not reported)",
    regions = "France (Institut Paoli Calmettes and Assistance Publique Hopitaux de Marseille)",
    sampling = "Rich sampling, eight samples per patient from 2 to 12 h after the start of infusion (Methods); the Discussion describes the same protocol as six samples. Busulfan in plasma by LC-MS/MS, 50-3000 ng/mL.",
    application_cohort = "71 adults (Table 1; January 2016 - September 2022): age 53 +/- 14 years (range 19-71), 49% female, weight 73 +/- 17 kg (range 40-135), BSA 1.87 +/- 0.26 m2; acute myeloid leukaemia 31%, acute lymphoblastic leukaemia 21%, myelofibrosis 14%. Busulfan-fludarabine(-thiotepa) conditioning, 3 mg/kg once daily as a 3-h infusion for 4 (n = 60), 5 (n = 5) or 6 (n = 6) doses, targeting a cumulative AUC of 66.0, 82.6 or 87.6 mg*h/L; three samples (3.25, 9 and 12 h after the start of the first infusion) for Bayesian estimation.",
    notes = "The population prior was built by the standard two-stage method: weighted-least-squares one-compartment fits (proportional error variance) to the 16 rich profiles, then the mean and covariance of the individual (1/V, k) estimates. The 71-patient real-world cohort is the application set for Bayesian adaptive dosing in Kinetic Pro (MIPPS) and was not used to estimate the prior."
  )

  ini({
    # Population mean vector of the standard two-stage analysis (Methods):
    #   [1/V, k] = [2.277e-2 L^-1, 2.448e-1 h^-1]
    # The mean of 1/V is the typical value on the scale the paper estimates,
    # so the typical V is its reciprocal, 43.92 L.
    lvc <- log(1 / 2.277e-2); label("Typical volume of distribution, the reciprocal of the population mean of 1/V (L)") # Methods, mean population parameters: 1/V = 2.277e-2 L^-1
    lkel <- log(2.448e-1); label("Typical elimination rate constant k (1/h)") # Methods, mean population parameters: k = 2.448e-1 h^-1

    # Covariance matrix of [1/V, k] (Methods), on the linear scale:
    #   var(1/V) = 2.018e-5 L^-2, cov(1/V, k) = 3.104e-5 L^-1 h^-1,
    #   var(k) = 5.338e-4 h^-2
    etavcinv + etakel ~ c(2.018e-5, 3.104e-5, 5.338e-4) # Methods, associated covariance matrix of the two-stage mean vector

    # Residual error: the first-stage fits used a proportional error variance
    # model (Methods) but its magnitude, and the error model of the Bayesian
    # step, are not reported. Held at zero so simulations return individual
    # predictions; see the vignette.
    propSd <- fixed(0); label("Proportional residual error (fraction); magnitude not reported") # Methods: proportional error variance model, value not reported
  })
  model({
    # Individual parameters on the paper's Gaussian scale
    vcinv <- 1 / exp(lvc) + etavcinv
    kel <- exp(lkel) + etakel

    vc <- 1 / vcinv
    cl <- kel * vc

    d/dt(central) <- -kel * central

    # Dose in mg, vc in L: central / vc is mg/L
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
