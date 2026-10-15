Agyeman_2022_covid19_tcle_neant <- function() {
  description <- paste(
    "Target-cell-limited model with eclipse phase (TCLE) for SARS-CoV-2",
    "upper-respiratory viral load, fitted by Agyeman et al. (2022) to dataset",
    "B: 321 hospitalised COVID-19 patients (563 nasopharyngeal samples within",
    "14 days of symptom onset) of the French COVID cohort analysed by Neant et",
    "al. (2021).",
    "Four ODE states: uninfected target cells 'target' (T, cells/mL) are",
    "infected by free virus at rate beta and become latently (eclipse-phase)",
    "infected cells 'eclipse' (I1), which convert at rate k to productively",
    "infected cells 'infected' (I2) that die at rate delta and produce free",
    "virus 'virus' (V, copies/mL) at rate rho; free virus is cleared at rate c",
    "and is also lost by entry into target cells (beta*T*V). Time is days since",
    "INFECTION (T(0) = 1.3e5 cells/mL, I1(0) = 1/30 cells/mL, I2(0) = 0,",
    "V(0) = 0.1 copies/mL); symptom onset occurs at the fixed incubation",
    "period tinc = 1.5 day, so an observation at s days after symptom onset is",
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
    n_subjects = 321L,
    n_observations = 563L,
    n_studies = NA_integer_,
    age_range = "Not reported in Agyeman 2022 (see Neant et al. 2021).",
    weight_range = "Not reported.",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported.",
    disease_state = "COVID-19 (SARS-CoV-2 infection) requiring hospitalisation in conventional wards or intensive care units; 78 deaths during follow-up.",
    dose_range = "No study drug; patients received routine care including antiviral, antibiotic, antifungal or corticosteroid treatment, which the model does not represent.",
    regions = "France (French COVID cohort).",
    notes = paste(
      "Dataset B of Agyeman 2022 (Section 2.1): nasopharyngeal viral load",
      "data from the prospective French COVID cohort of hospitalised patients",
      "analysed by Neant et al. (Proc Natl Acad Sci U S A 2021;118(8)).",
      "Restricted to 0-14 days after symptom onset. Viral loads below the",
      "limit of detection were censored in the fit. Patient demographics are",
      "given in the Neant 2021 source publication, not in Agyeman 2022."
    )
  )

  ini({
    lbeta <- log(9.79e-5)
    label("Rate constant for virus infection beta ((copies/mL)^-1 day^-1)") # Table S3, TCLE column: beta = 9.79 x 10^-5 (RSE 7.8%)
    ldelta <- log(0.81)
    label("Death rate of productively infected cells delta (1/day)") # Table S3, TCLE column: delta = 0.81 1/day (RSE 12.4%)
    lrho <- log(2.84e3)
    label("Viral production rate rho (copies/mL per infected cell/mL per day)") # Table S3, TCLE column: rho = 2.84 x 10^3 copies/ml.day^-1 (RSE 3.5%)
    lc <- log(2.62)
    label("Viral clearance rate c (1/day)") # Table S3, TCLE column: c = 2.62 1/day (RSE 10%)
    lk <- log(8.87)
    label("Conversion rate from eclipse-phase to productively infected cells k (1/day)") # Table S3, TCLE column: k = 8.87 1/day (RSE 9.9%)

    lT0 <- fixed(log(1.3e5))
    label("Target cells at infection T(0) (cells/mL)") # Table S3 footnote: T(0) = 1.3 x 10^5 cells/ml (for TCLE)
    lI10 <- fixed(log(1 / 30))
    label("Eclipse-phase cells at infection I1(0) (cells/mL)") # Table S3 footnote: I1(0) = 1/30 cells/ml (for TCLE)
    lV0 <- fixed(log(0.1))
    label("Free virus at infection V(0) (copies/mL)") # Table S3 footnote: V(0) = 0.1 copies/ml (for TCLE)
    tinc <- fixed(1.5)
    label("Incubation period from infection to symptom onset (day)") # Results and Figure S2: 1.5 day gave the lowest BIC for TCLE with dataset B

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
