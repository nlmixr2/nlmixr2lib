Agyeman_2022_covid19_tcle_gastine <- function() {
  description <- paste(
    "Target-cell-limited model with eclipse phase (TCLE) for SARS-CoV-2",
    "upper-respiratory viral load, fitted by Agyeman et al. (2022) to dataset",
    "A: 252 untreated COVID-19 patients (747 samples within 14 days of symptom",
    "onset) compiled by the Gastine et al. (2021) patient-level meta-analysis.",
    "Four ODE states: uninfected target cells 'target' (T, cells/mL) are",
    "infected by free virus at rate beta and become latently (eclipse-phase)",
    "infected cells 'eclipse' (I1), which convert at rate k to productively",
    "infected cells 'infected' (I2) that die at rate delta and produce free",
    "virus 'virus' (V, copies/mL) at rate rho; free virus is cleared at rate c",
    "and is also lost by entry into target cells (beta*T*V). Time is days since",
    "INFECTION (T(0) = 1.3e5 cells/mL, I1(0) = 1/30 cells/mL, I2(0) = 0,",
    "V(0) = 0.1 copies/mL); symptom onset occurs at the fixed incubation",
    "period tinc = 0.5 day, so an observation at s days after symptom onset is",
    "simulated at time s + tinc (output tsymptom = time - tinc). Typical",
    "values only: the paper estimated log-normal IIV on each estimated",
    "parameter and an additive residual error on the natural-log viral load",
    "but did not report their magnitudes, so the etas and residual SD are",
    "carried as fixed(0). The eclipse-phase state is declared paper-specific",
    "via paper_specific_compartments, as in Dodds_2021_covid19_viralkinetic."
  )
  reference <- paste(
    "Agyeman AA, You T, Chan PLS, Lonsdale DO, Hadjichrysanthou C, Mahungu T,",
    "Wey EQ, Lowe DM, Lipman MCI, Breuer J, Kloprogge F, Standing JF.",
    "Comparative assessment of viral dynamic models for SARS-CoV-2 for",
    "pharmacodynamic assessment in early treatment trials.",
    "Br J Clin Pharmacol. 2022;88(12):5428-5433. doi:10.1111/bcp.15518"
  )
  vignette <- "Agyeman_2022_covid19_viral_dynamics"

  paper_specific_compartments <- c("eclipse")

  units <- list(
    time = "day",
    dosing = "no dosing (untreated natural-history viral dynamics)",
    concentration = "copies/mL (observed as log10 copies/mL); cells/mL for target, eclipse and infected cells"
  )

  covariateData <- list()

  compartmentData <- list(
    target = list(
      analyte = "uninfected respiratory epithelial target cells",
      units = "cells/mL",
      specimen = "tissue",
      verified = TRUE
    ),
    eclipse = list(
      analyte = "latently (eclipse-phase) infected respiratory epithelial cells",
      units = "cells/mL",
      specimen = "tissue",
      verified = TRUE
    ),
    infected = list(
      analyte = "productively infected respiratory epithelial cells",
      units = "cells/mL",
      specimen = "tissue",
      verified = TRUE
    ),
    virus = list(
      analyte = "SARS-CoV-2 viral RNA (free virus)",
      units = "copies/mL",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 252L,
    n_observations = 747L,
    n_studies = NA_integer_,
    age_range = "Not reported in Agyeman 2022 (see Gastine et al. 2021).",
    weight_range = "Not reported.",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported.",
    disease_state = "COVID-19 (SARS-CoV-2 infection), predominantly mild disease with one reported death; untreated patients only.",
    dose_range = "No antiviral treatment (untreated patients selected from the dataset).",
    regions = "Multinational (patient-level data from published case series).",
    notes = paste(
      "Dataset A of Agyeman 2022 (Section 2.1): patient-level viral load data",
      "compiled by the systematic review and meta-analysis of Gastine et al.",
      "(Clin Pharmacol Ther 2021;110:321-333). Restricted to untreated",
      "patients, upper-respiratory-tract sampling sites and 0-14 days after",
      "symptom onset. Viral loads below the limit of detection were censored",
      "in the fit. Patient demographics are given in the Gastine 2021 source",
      "publication, not in Agyeman 2022."
    )
  )

  ini({
    lbeta <- log(7.76e-5)
    label("Rate constant for virus infection beta ((copies/mL)^-1 day^-1)") # Table S2, TCLE column: beta = 7.76 x 10^-5 (RSE 6.2%)
    ldelta <- log(0.69)
    label("Death rate of productively infected cells delta (1/day)") # Table S2, TCLE column: delta = 0.69 1/day (RSE 12.6%)
    lrho <- log(2.96e3)
    label("Viral production rate rho (copies/mL per infected cell/mL per day)") # Table S2, TCLE column: rho = 2.96 x 10^3 copies/ml.day^-1 (RSE 4.3%)
    lc <- log(11.6)
    label("Viral clearance rate c (1/day)") # Table S2, TCLE column: c = 11.6 1/day (RSE 9.4%)
    lk <- log(9.23)
    label("Conversion rate from eclipse-phase to productively infected cells k (1/day)") # Table S2, TCLE column: k = 9.23 1/day (RSE 10.6%)

    lT0 <- fixed(log(1.3e5))
    label("Target cells at infection T(0) (cells/mL)") # Table S2 footnote: T(0) = 1.3 x 10^5 cells/ml (for TCLE)
    lI10 <- fixed(log(1 / 30))
    label("Eclipse-phase cells at infection I1(0) (cells/mL)") # Table S2 footnote: I1(0) = 1/30 cells/ml (for TCLE)
    lV0 <- fixed(log(0.1))
    label("Free virus at infection V(0) (copies/mL)") # Table S2 footnote: V(0) = 0.1 copies/ml (for TCLE)
    tinc <- fixed(0.5)
    label("Incubation period from infection to symptom onset (day)") # Results and Figure S1: 0.5 day gave the lowest BIC for TCLE with dataset A

    # Section 2.3: log-normal IIV on every estimated parameter; variances not
    # reported in the paper or Supporting Information, so held at zero.
    etalbeta ~ fixed(0)
    etaldelta ~ fixed(0)
    etalrho ~ fixed(0)
    etalc ~ fixed(0)
    etalk ~ fixed(0)

    # Section 2.3: additive (normal) residual error on the log-transformed viral
    # load; magnitude not reported, so held at zero.
    addSd <- fixed(0)
    label("Additive residual SD on log10 viral load (log10 copies/mL); not reported")
  })

  model({
    beta <- exp(lbeta + etalbeta)
    delta <- exp(ldelta + etaldelta)
    rho <- exp(lrho + etalrho)
    c <- exp(lc + etalc)
    k <- exp(lk + etalk)
    T0 <- exp(lT0)
    I10 <- exp(lI10)
    V0 <- exp(lV0)

    # Equation 4. The printed dI1/dt carries a leading minus on beta*T*V; the
    # infection flux is a gain for I1 (Figure 1D, and Equation 7 reproduces the
    # tabulated R0 only with this sign).
    d/dt(target) <- -beta * target * virus
    d/dt(eclipse) <- beta * target * virus - k * eclipse
    d/dt(infected) <- k * eclipse - delta * infected
    d/dt(virus) <- rho * infected - c * virus - beta * target * virus
    target(0) <- T0
    eclipse(0) <- I10
    infected(0) <- 0
    virus(0) <- V0

    tsymptom <- t - tinc

    log10_viral_load <- log10(virus)
    log10_viral_load ~ add(addSd)
  })
}
