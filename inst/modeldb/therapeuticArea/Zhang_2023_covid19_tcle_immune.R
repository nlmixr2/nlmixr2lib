Zhang_2023_covid19_tcle_immune <- function() {
  description <- paste(
    "QSP (within-host viral dynamics). Target-cell-limited model with an",
    "eclipse phase and an immune-response compartment that increases viral",
    "clearance (model M7, TCLE + increasing viral clearance) for SARS-CoV-2",
    "upper-respiratory-tract viral load, fitted by Zhang et al. (2023) to 2155",
    "measurements from 162 household and non-household contacts of the UK",
    "ATACCC community cohort (no antiviral treatment). Five ODE states:",
    "uninfected target cells 'target' (T, cells/mL) are infected by free virus",
    "at rate beta and enter the eclipse phase 'eclipse' (I1), convert at rate k",
    "to productively infected cells 'infected' (I2) that die at rate delta and",
    "shed virus 'virus' (V, copies/mL) at rate rho; free virus is cleared at",
    "rate c, lost by cell entry (beta*T*V) and additionally cleared by the",
    "immune response 'immune' (F, produced by I2 and eliminated at rate dF)",
    "through phi*F/(F+theta). Time is days since contact with the infectious",
    "case. The paper reports the basic reproductive number R0 rather than beta;",
    "beta is derived in model() from R0 = rho*beta*T0/(delta*(c + beta*T0)).",
    "The antiviral blockade of infection used in the paper's simulations,",
    "beta*(1 - eta*H(t)) with H(t) = 1 for t >= t*, is available through",
    "imax_beta (eta) and t_intervention (t*), both defaulting to no treatment.",
    "The paper does not report the IIV or residual-error estimates; the IIV",
    "variances (0.01 on beta, delta, V0, phi and theta) and residual SD (3 on",
    "the natural-log scale) are the values in the authors' published",
    "simulation script (Supplement S2). With them this parameter set",
    "reproduces Figure 3. The age-on-delta model fitted to unvaccinated",
    "patients is Zhang_2023_covid19_tcle_immune_age."
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

  covariateData <- list()

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
    n_subjects = 162L,
    n_observations = 2155L,
    n_studies = 1L,
    age_range = "7-78 years (median 36); 14.8% < 18 years, 3.1% >= 65 years (Table 1).",
    weight_range = "Not reported.",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported.",
    disease_state = paste(
      "Community SARS-CoV-2 infection in household (124, 76.5%) and",
      "non-household (38, 23.5%) contacts of index cases; 90 symptomatic,",
      "56 asymptomatic, 16 unknown. Variants: pre-alpha 50, alpha 42, delta 70.",
      "Vaccination: unvaccinated 114, partially vaccinated 10, fully",
      "vaccinated 38."
    ),
    dose_range = "No antiviral treatment.",
    regions = "United Kingdom (ATACCC study).",
    notes = paste(
      "Assessment of Transmission and Contagiousness of COVID-19 in Contacts",
      "(ATACCC) longitudinal cohort; daily upper-respiratory-tract swabs from",
      "enrolment (within 5 days of index-case symptom onset). ORF1ab Ct values",
      "converted to copies/mL by 133.3333 * exp((37.933 - Ct)/1.418) (Eq. 1);",
      "Ct = 40 (31.0 copies/mL) is the limit of detection and lower values were",
      "censored in the fit. Time zero is contact with the infectious case,",
      "imputed for household contacts as 7 days before the viral-load peak.",
      "Fitted with nlmixr 2.0.7 (SAEM)."
    )
  )

  ini({
    lR0 <- log(13.74)
    label("Basic reproductive number R0 (unitless)") # Table 3, c and k fixed column: R0 = 13.74 (RSE 0.78%)
    ldelta <- log(1.04)
    label("Loss rate of productively infected cells delta (1/day)") # Table 3, c and k fixed column: delta = 1.04 1/day
    lphi <- log(0.24)
    label("Maximal immune enhancement of viral clearance phi (1/day)") # Table 3, c and k fixed column: phi = 0.24 (RSE 5.44%)
    ltheta <- log(542.46)
    label("Immune-response level for half-maximal effect theta (arbitrary units)") # Table 3, c and k fixed column: theta = 542.46 (RSE 8.01%)
    lV0 <- log(23.33)
    label("Viral load at contact V(0) (copies/mL)") # Table 3, c and k fixed column: V0 = 23.33 copies/mL (RSE 7.70%)
    lrho <- fixed(log(6000))
    label("Viral production rate per productively infected cell rho (copies/cell/day)") # Table 3: rho = 6000 (fixed)
    lc <- fixed(log(20))
    label("Free-virus clearance rate c (1/day)") # Table 3, c and k fixed column: c = 20 (fixed)
    lk <- fixed(log(3))
    label("Eclipse-to-productive transition rate k (1/day)") # Table 3, c and k fixed column: k = 3 (fixed)
    ldF <- fixed(log(0.4))
    label("Elimination rate of the immune response dF (1/day)") # Table 3: dF = 0.4 (fixed); Methods 'Assumptions on parameters'
    lT0 <- fixed(log(1.33e5))
    label("Target cells at contact T(0) (cells/mL)") # Table 3: T0 = 1.33 x 10^5 (fixed); Methods 'Assumptions on parameters'

    imax_beta <- fixed(0)
    label("Antiviral efficacy eta: fractional inhibition of beta at t >= t_intervention (0-1)") # Methods 'Simulations for antiviral therapies', Eq. 3-6: (1 - eta*H(t)) * beta; 0 = untreated
    t_intervention <- fixed(1e6)
    label("Treatment initiation time t* (day since contact); default never") # Methods 'Simulations for antiviral therapies': H(t) = 0 if t < t*, otherwise 1

    # Supplement S2 simulation code: ombeta, omdelta, omv0, omphi, omtheta =
    # 0.01 (variances, drawn as rnorm(sd = sqrt(om))), equal to the estimation
    # initial values; final IIV estimates are not reported in the paper.
    etalbeta ~ 0.01
    etaldelta ~ 0.01
    etalV0 ~ 0.01
    etalphi ~ 0.01
    etaltheta ~ 0.01

    addSd <- 1.3029
    label("Additive residual SD on log10 viral load (log10 copies/mL)") # Supplement S2 simulation code: exp(rnorm(mean = log(pred), sd = 3)), i.e. SD 3 on natural-log viral load = 3/log(10) on log10
  })

  model({
    R0 <- exp(lR0)
    delta_typ <- exp(ldelta)
    delta <- delta_typ * exp(etaldelta)
    phi <- exp(lphi + etalphi)
    theta <- exp(ltheta + etaltheta)
    V0 <- exp(lV0 + etalV0)
    rho <- exp(lrho)
    c <- exp(lc)
    k <- exp(lk)
    dF <- exp(ldF)
    T0 <- exp(lT0)

    # Results 'Parameter estimation': R0 = rho*beta*T0 / (delta*(c + beta*T0 +
    # phi*F0/(F0 + theta))) with F0 = F(0) = 0, solved for the typical beta.
    # IIV then acts on beta (Supplement S2: beta = exp(tbeta + ebeta)).
    beta_typ <- R0 * delta_typ * c / (T0 * (rho - R0 * delta_typ))
    beta <- beta_typ * exp(etalbeta)

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

    # Fraction of target cells still uninfected (Figure 3d-f, 'Target cell ratio').
    f_target <- target / T0

    log10_viral_load <- log10(virus)
    log10_viral_load ~ add(addSd)
  })
}
