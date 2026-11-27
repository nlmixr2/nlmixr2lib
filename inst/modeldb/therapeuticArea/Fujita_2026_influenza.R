Fujita_2026_influenza <- function() {
  description <- paste(
    "Target-cell-limited viral dynamics model of nasal-wash influenza A",
    "(A/Texas/91 H1N1) and influenza B (B/Yamagata/16/88) viral titre in",
    "healthy adults experimentally inoculated in three oseltamivir dose-ranging",
    "challenge trials (Fujita 2026; Roche protocols PV15616, NP15717 and",
    "NP15827).",
    "Two ODE states under the quasi-steady-state approximation for free virus:",
    "the fraction of uninfected target cells 'target' (f, 1 at inoculation) is",
    "depleted at rate beta*f*V, and the viral titre 'virus' (V, TCID50/mL)",
    "grows at gamma*f*V and declines at the infected-cell death rate delta.",
    "Time is days since inoculation. Influenza type B (INFLUENZA_B = 1) shifts",
    "all four parameters (beta, gamma, delta, V(0)) on the log scale.",
    "Oseltamivir is NOT represented: the authors fitted the model to viral",
    "titres from day 1 onward (after treatment started) and read any drug",
    "effect from differences in the individual estimates, so the typical",
    "values describe the pooled placebo and oseltamivir cohort.",
    "Residual error was not reported and is carried as fixed(0)."
  )
  reference <- paste(
    "Fujita Y, Akao M, Tatematsu D, Iwami S, Nakamura N, Iwanami S.",
    "Modeling viral shedding and symptom outcomes in oseltamivir-treated",
    "experimental influenza infection.",
    "PLoS One. 2026;21(2):e0342676. doi:10.1371/journal.pone.0342676.",
    "Population parameters are from S8 Table; the model is Equation 1 of",
    "Materials and methods."
  )
  vignette <- "Fujita_2026_influenza"

  units <- list(
    time = "day",
    dosing = "no dosing (the model has no drug input; the inoculum is the initial viral titre V(0))",
    concentration = "TCID50/mL of nasal wash for virus (observed as log10 TCID50/mL); target is a unitless fraction of its value at inoculation"
  )

  covariateData <- list(
    INFLUENZA_B = list(
      description = "Infecting influenza virus type indicator: 1 = influenza B, 0 = influenza A.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = influenza A (A/Texas/91 (H1N1), Flu A study PV15616)",
      notes = paste(
        "1 = influenza B/Yamagata/16/88 (Flu B studies NP15717 and NP15827).",
        "The paper calls this covariate 'genotype' (S8 Table column",
        "'Covariate for genotype B'); it is perfectly collinear with study, so",
        "it also absorbs any protocol difference between the Flu A and Flu B",
        "trials (inoculum 10^6 vs 10^7 TCID50, treatment start 28 h vs 24 h)."
      ),
      source_name = "genotype"
    )
  )

  compartmentData <- list(
    target = list(
      analyte = "uninfected respiratory epithelial target cells (fraction of the value at inoculation)",
      units = "(fraction)",
      specimen = "tissue",
      verified = TRUE
    ),
    virus = list(
      analyte = "influenza virus infectious titre",
      units = "TCID50/mL",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 103L,
    n_studies = 3L,
    age_range = "17-64 years; mean 25.9 years (22.0 years Flu A, 28.8 years Flu B) among the 103 modelled participants (S7 Table)",
    weight_range = "Not reported.",
    sex_female_pct = 56.3,
    race_ethnicity = "Not reported.",
    disease_state = paste(
      "Healthy adult volunteers intranasally inoculated with influenza",
      "A/Texas/91 (H1N1) (10^6 TCID50) or influenza B/Yamagata/16/88",
      "(10^7 TCID50)."
    ),
    dose_range = paste(
      "Oral oseltamivir for 5 days from 28 h (Flu A) or 24 h (Flu B) after",
      "inoculation. Modelled participants: Flu A placebo 9, 20 mg BID 8,",
      "100 mg BID 8, 200 mg BID 10, 200 mg QD 9; Flu B placebo 28,",
      "75 mg BID 26, 150 mg BID 5 (S7 and S9 Tables)."
    ),
    regions = "Not reported.",
    notes = paste(
      "Of 257 inoculated participants (Flu A 80, Flu B study A 60, Flu B",
      "study B 117), the 103 with at least 3 viral-titre values above the",
      "detection limit were used for the fit (44 influenza A, 59 influenza B;",
      "Fig 1A and S7 Table). Daily nasal-wash titres were digitised from",
      "Hooker and Ganusov's published dataset; the day-0 sample was excluded",
      "and the model was fitted in Monolix (SAEM) with log-normal individual",
      "parameters and virus type as a covariate."
    )
  )

  ini({
    lbeta <- log(5.55e-4)
    label("Scaled rate constant for infection beta, influenza A ((TCID50/mL)^-1 day^-1)") # S8 Table: beta population value 5.55 x 10^-4 (unit printed as (RNA copies/mL)^-1 day^-1)
    lgamma <- log(14.7)
    label("Maximum rate constant for viral replication gamma, influenza A (1/day)") # S8 Table: gamma population value 14.7 day^-1
    ldelta <- log(5.18)
    label("Death rate of infected cells delta, influenza A (1/day)") # S8 Table: delta population value 5.18 day^-1
    lV0 <- log(0.492)
    label("Viral titre at inoculation V(0), influenza A (TCID50/mL)") # S8 Table: V(0) population value 0.492 TCID50/mL

    e_influenza_b_beta <- 1.65
    label("Influenza B shift in log(beta) (unitless)") # S8 Table: covariate for genotype B on beta = 1.65 (Wald p = 1.37 x 10^-2)
    e_influenza_b_gamma <- -0.739
    label("Influenza B shift in log(gamma) (unitless)") # S8 Table: covariate for genotype B on gamma = -0.739 (Wald p = 2.47 x 10^-13)
    e_influenza_b_delta <- -0.595
    label("Influenza B shift in log(delta) (unitless)") # S8 Table: covariate for genotype B on delta = -0.595 (Wald p = 2.98 x 10^-3)
    e_influenza_b_V0 <- 0.959
    label("Influenza B shift in log(V(0)) (unitless)") # S8 Table: covariate for genotype B on V(0) = 0.959 (Wald p = 0.659)

    # S8 Table column 'Variance of individual parameters'. Read as the Monolix
    # standard deviation omega and squared here; see the vignette for the
    # evidence from the S9 Table individual estimates.
    etalbeta ~ 7.1824 # S8 Table: beta 2.68; 2.68^2
    etalgamma ~ 0.00912025 # S8 Table: gamma 9.55 x 10^-2; 0.0955^2
    etaldelta ~ 0.168921 # S8 Table: delta 0.411; 0.411^2
    etalV0 ~ 3.8416 # S8 Table: V(0) 1.96; 1.96^2

    addSd <- fixed(0)
    label("Additive residual SD on log10 viral titre (log10 TCID50/mL); not reported")
  })

  model({
    beta <- exp(lbeta + e_influenza_b_beta * INFLUENZA_B + etalbeta)
    gamma <- exp(lgamma + e_influenza_b_gamma * INFLUENZA_B + etalgamma)
    delta <- exp(ldelta + e_influenza_b_delta * INFLUENZA_B + etaldelta)
    V0 <- exp(lV0 + e_influenza_b_V0 * INFLUENZA_B + etalV0)

    # Equation 1: target-cell fraction f and viral titre V, with free virus
    # in quasi-steady state with infected cells.
    d/dt(target) <- -beta * target * virus
    d/dt(virus) <- gamma * target * virus - delta * virus
    target(0) <- 1
    virus(0) <- V0

    log10_viral_load <- log10(virus)
    log10_viral_load ~ add(addSd)
  })
}
