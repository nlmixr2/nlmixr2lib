Zhang_2023_covid19_tcle_immune_age <- function() {
  description <- paste(
    "QSP (within-host viral dynamics). Final simulation model of Zhang et al.",
    "(2023): the target-cell-limited model with an eclipse phase and an",
    "immune response that increases viral clearance (M7), with age as a power",
    "covariate on the loss rate of productively infected cells delta,",
    "re-fitted to the unvaccinated patients of the UK ATACCC community cohort",
    "and used for the paper's simulation of the tixagevimab-cilgavimab phase",
    "III trial (Figure 4).",
    "Five ODE states: uninfected target cells 'target', eclipse-phase cells",
    "'eclipse', productively infected cells 'infected', free virus 'virus'",
    "(copies/mL) and immune response 'immune'. Time is days since contact with",
    "the infectious case. The antiviral blockade of infection used in the",
    "paper's simulations, beta*(1 - eta*H(t)) with H(t) = 1 for t >= t*, is",
    "available through imax_beta (eta) and t_intervention (t*), both",
    "defaulting to no treatment. Parameter values are the estimates the",
    "authors coded into their published simulation script (Supplement S2);",
    "the IIV variances (0.01) and residual SD (3 on the natural-log scale) are",
    "the values used in that script, not tabulated final estimates. The",
    "all-patient fit without the age effect is Zhang_2023_covid19_tcle_immune."
  )
  reference <- paste(
    "Zhang S, Agyeman AA, Hadjichrysanthou C, Standing JF.",
    "SARS-CoV-2 viral dynamic modeling to inform model selection and timing",
    "and efficacy of antiviral therapy.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(10):1450-1460.",
    "doi:10.1002/psp4.13022"
  )
  vignette <- "Zhang_2023_covid19_viral_dynamics"

  paper_specific_compartments <- c("eclipse", "immune")

  units <- list(
    time = "day",
    dosing = "no dosing (antiviral effect enters as a fractional inhibition of beta)",
    concentration = "copies/mL (observed as log10 copies/mL); cells/mL for target, eclipse and infected cells; arbitrary units for the immune response"
  )

  covariateData <- list(
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power covariate on delta, (AGE/35)^e_age_delta. Supplement S2 model code: log_age = log(Age/35). Fitted ages spanned 7-78 years; the paper's trial simulation used ages 18-78.",
      source_name = "Age"
    )
  )

  compartmentData <- list(
    target = list(
      analyte = "uninfected upper-respiratory-tract epithelial target cells",
      units = "cells/mL",
      specimen = "tissue",
      verified = TRUE
    ),
    eclipse = list(
      analyte = "eclipse-phase (infected, not yet virus-producing) epithelial cells",
      units = "cells/mL",
      specimen = "tissue",
      verified = TRUE
    ),
    infected = list(
      analyte = "productively infected epithelial cells",
      units = "cells/mL",
      specimen = "tissue",
      verified = TRUE
    ),
    virus = list(
      analyte = "SARS-CoV-2 viral RNA (free virus)",
      units = "copies/mL",
      specimen = "tissue",
      verified = TRUE
    ),
    immune = list(
      analyte = "immune response enhancing viral clearance (latent, unitless level)",
      units = "arbitrary units",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 114L,
    n_studies = 1L,
    age_range = "Not reported for the unvaccinated subset; whole cohort 7-78 years (median 36).",
    weight_range = "Not reported.",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported.",
    disease_state = paste(
      "Community SARS-CoV-2 infection (pre-alpha, alpha and delta variants) in",
      "unvaccinated household and non-household contacts of index cases."
    ),
    dose_range = "No antiviral treatment.",
    regions = "United Kingdom (ATACCC study).",
    notes = paste(
      "Methods 'The tixagevimab-cilgavimab clinical trial and simulation': the",
      "simulation parameters came from the model fitted to unvaccinated",
      "patients (114 of the 162 ATACCC patients in Table 1). Age, vaccination",
      "status, contact type, variant and symptomatic status were screened;",
      "only age was tested numerically (on beta, delta, V0, phi and theta;",
      "Supplement S5) and age on delta gave the lowest BIC. The other",
      "covariates showed no trend on graphical inspection and were not",
      "retained."
    )
  )

  ini({
    lbeta <- -14.65386721
    label("Rate constant for virus infection beta, log scale (mL/(copies*day))") # Supplement S2 simulation code: tbeta = -14.65386721 '## beta from NM'
    ldelta <- -0.02423612
    label("Loss rate of productively infected cells delta at age 35 years, log scale (1/day)") # Supplement S2 simulation code: tdelta = -0.02423612 '## delta from NM'
    lphi <- -1.54922962
    label("Maximal immune enhancement of viral clearance phi, log scale (1/day)") # Supplement S2 simulation code: tphi = -1.54922962
    ltheta <- 6.28098060
    label("Immune-response level for half-maximal effect theta, log scale (arbitrary units)") # Supplement S2 simulation code: ttheta = 6.28098060
    lV0 <- 3.19332606
    label("Viral load at contact V(0), log scale (copies/mL)") # Supplement S2 simulation code: tv0 = 3.19332606 '## V0 from NM'
    lrho <- fixed(log(6000))
    label("Viral production rate per productively infected cell rho (copies/cell/day)") # Table 3: rho = 6000 (fixed); Supplement S2: trho = fix(log(6000))
    lc <- fixed(log(20))
    label("Free-virus clearance rate c (1/day)") # Table 3: c = 20 (fixed); Supplement S2: tc = fix(log(20))
    lk <- fixed(log(3))
    label("Eclipse-to-productive transition rate k (1/day)") # Table 3: k = 3 (fixed); Supplement S2: tk = fix(log(3))
    ldF <- fixed(log(0.4))
    label("Elimination rate of the immune response dF (1/day)") # Table 3: dF = 0.4 (fixed); Supplement S2: tdF = fix(log(0.4))
    lT0 <- fixed(log(1.33e5))
    label("Target cells at contact T(0) (cells/mL)") # Table 3: T0 = 1.33 x 10^5 (fixed); Supplement S2: A_T(0) = 1.33 * 10^5
    e_age_delta <- -0.11
    label("Power exponent of age on delta, (AGE/35)^e_age_delta (unitless)") # Supplement S2 simulation code: delta = exp(tdelta + edelta + -0.11 * log(Sims$age[i]/35))

    imax_beta <- fixed(0)
    label("Antiviral efficacy eta: fractional inhibition of beta at t >= t_intervention (0-1)") # Methods 'Simulations for antiviral therapies', Eq. 3-6: (1 - eta*H(t)) * beta; 0 = untreated
    t_intervention <- fixed(1e6)
    label("Treatment initiation time t* (day since contact); default never") # Methods 'Simulations for antiviral therapies': H(t) = 0 if t < t*, otherwise 1

    # Supplement S2 simulation code: ombeta, omdelta, omv0, omphi, omtheta =
    # 0.01 (variances, drawn as rnorm(sd = sqrt(om))). Same as the estimation
    # initial values; final estimates are not reported in the paper.
    etalbeta ~ 0.01
    etaldelta ~ 0.01
    etalV0 ~ 0.01
    etalphi ~ 0.01
    etaltheta ~ 0.01

    addSd <- 1.3029
    label("Additive residual SD on log10 viral load (log10 copies/mL)") # Supplement S2 simulation code: exp(rnorm(mean = log(pred), sd = 3)), i.e. SD 3 on natural-log viral load = 3/log(10) on log10
  })

  model({
    beta <- exp(lbeta + etalbeta)
    delta <- exp(ldelta + etaldelta + e_age_delta * log(AGE / 35))
    phi <- exp(lphi + etalphi)
    theta <- exp(ltheta + etaltheta)
    V0 <- exp(lV0 + etalV0)
    rho <- exp(lrho)
    c <- exp(lc)
    k <- exp(lk)
    dF <- exp(ldF)
    T0 <- exp(lT0)

    # Eq. 3-7: blockade of de novo infection, (1 - eta*H(t)) * beta.
    treated <- (t >= t_intervention)
    beta_eff <- beta * (1 - imax_beta * treated)

    # Supplement S1 M7 (TCLE + increasing viral clearance) and Eq. 3-7.
    d/dt(target) <- -beta_eff * virus * target
    d/dt(eclipse) <- beta_eff * virus * target - k * eclipse
    d/dt(infected) <- k * eclipse - delta * infected
    d/dt(virus) <- rho * infected - c * virus -
      phi * immune / (immune + theta) * virus -
      beta_eff * virus * target
    d/dt(immune) <- infected - dF * immune
    # Methods 'Assumptions on parameters': T = T0, I1 = 0, I2 = 1/30, V = V0.
    target(0) <- T0
    eclipse(0) <- 0
    infected(0) <- 1 / 30
    virus(0) <- V0
    immune(0) <- 0

    # Basic reproductive number with F(0) = 0 (Results 'Parameter estimation').
    R0 <- rho * beta * T0 / (delta * (c + beta * T0))
    # Fraction of target cells still uninfected (Figure 3d-f, 'Target cell ratio').
    f_target <- target / T0

    log10_viral_load <- log10(virus)
    log10_viral_load ~ add(addSd)
  })
}
