Neant_2023_covid19_news2_soc <- function() {
  description <- paste(
    "Joint SARS-CoV-2 viral-kinetic / clinical-score model of Neant et al.",
    "(2023), standard-of-care arm fit: the final model structure re-estimated",
    "on the 329 DisCoVeRy patients randomised to standard of care (SoC) alone",
    "(Table S3), used by the authors to compare remdesivir and SoC empirical",
    "Bayes estimates. No remdesivir effect. Viral load follows a target-cell",
    "limited model with an eclipse phase (states target, eclipse, infected,",
    "infectious 'virus' and noninfectious 'virus_ni'; beta derived from R0 =",
    "beta * p * T0 * mu / (c * delta)); time is days since INFECTION and the",
    "incubation period tinc maps it to days since symptom onset (output",
    "tsymptom). The normalised viral load drives an effect compartment with",
    "k1e = ke0, and NEWS-2 is a Hill function of log10(effect) between the",
    "pre-infection baseline N0 (logit-normal, bounded 0-3) and the score",
    "maximum 20. delta and ke0 take separate typical values for patients",
    "aged >= 65 years (AGE_GE65). The main final model fitted to both arms",
    "is Neant_2023_covid19_news2."
  )
  reference <- paste(
    "Neant N, Lingas G, Gaymard A, Belhadi D, Hites M, Staub T, Greil R,",
    "Paiva JA, Poissy J, Peiffer-Smadja N, Costagliola D, Yazdanpanah Y,",
    "Bouscambert-Duchamp M, Gagneux-Brunon A, Ader F, Mentre F, Wallet F,",
    "Burdet C, Guedj J; the DisCoVeRy study group. Association between",
    "SARS-CoV-2 viral kinetics and clinical score evolution in hospitalized",
    "patients. CPT Pharmacometrics Syst Pharmacol. 2023;12(12):2027.",
    "doi:10.1002/psp4.13051"
  )
  vignette <- "Neant_2023_covid19_news2"

  # Eclipse-phase cells and the noninfectious virion pool are states of this
  # viral-kinetic framework with no register entry.
  paper_specific_compartments <- c("eclipse", "virus_ni")

  units <- list(
    time = "day",
    dosing = "no dosing (standard-of-care arm; untreated natural-history viral dynamics)",
    concentration = paste(
      "cells (target, eclipse, infected) and virions (virus, virus_ni) per",
      "whole nasopharyngeal tract; viral load observed as log10 copies per",
      "10^4 cells; effect compartment in copies per 10^4 cells; NEWS-2 in",
      "score points (0-20)"
    )
  )

  covariateData <- list(
    AGE_GE65 = list(
      description = "Age 65 years or older at randomisation (1 = age >= 65 y, 0 = age < 65 y)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (age < 65 years)",
      notes = paste(
        "Selects the age-specific typical values of the infected-cell loss",
        "rate delta and the effect-compartment rate constant ke0 (Table 2",
        "rows delta<65 / delta>=65 and ke0<65 / ke0>=65). The two strata share",
        "one random effect on each parameter."
      ),
      source_name = "age >= 65 years"
    )
  )

  compartmentData <- list(
    target = list(
      analyte = "uninfected nasopharyngeal epithelial target cells",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    eclipse = list(
      analyte = "eclipse-phase (not yet productive) infected epithelial cells",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    infected = list(
      analyte = "productively infected epithelial cells",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    virus = list(
      analyte = "infectious SARS-CoV-2 virions",
      units = "virions",
      specimen = "tissue",
      verified = TRUE
    ),
    virus_ni = list(
      analyte = "noninfectious SARS-CoV-2 virions",
      units = "virions",
      specimen = "tissue",
      verified = TRUE
    ),
    effect = list(
      analyte = "hypothetical effect compartment driven by the normalised viral load",
      units = "copies/10^4 cells",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 329L,
    n_studies = 1L,
    age_range = "Adults; median 64 years (IQR 53-72)",
    age_median = "64 years",
    weight_range = "Not reported",
    sex_female_pct = 32.5,
    race_ethnicity = "Not reported",
    disease_state = paste(
      "Hospitalised adults with laboratory-confirmed SARS-CoV-2 infection",
      "(COVID-19) randomised to standard of care alone, at most 20 days after",
      "symptom onset (median 9 days, IQR 7-11); NEWS-2 >= 7 at randomisation",
      "in 74.5%."
    ),
    dose_range = "Standard of care alone (no antiviral study drug)",
    regions = "France, Belgium, Portugal, Austria, Luxembourg (48 sites)",
    notes = paste(
      "Standard-of-care arm of the DisCoVeRy analysis set (NCT04315948):",
      "222 (67.5%) male, 160 (48.6%) aged >= 65 years (Table 1). This arm was",
      "also used to build the covariate model. Parameter estimates: Table S3."
    )
  )

  ini({
    # Viral-kinetic parameters (Table S3; log-normal IIV)
    lR0 <- log(6.17)
    label("Basic reproductive number R0 (unitless)") # Table S3: R0 = 6.17 (RSE 7.74%)
    ldelta_agelt65 <- log(0.89)
    label("Loss rate of productively infected cells, age < 65 years (1/day)") # Table S3: delta<65 = 0.89 1/day (RSE 6.02%)
    ldelta_agege65 <- log(0.80)
    label("Loss rate of productively infected cells, age >= 65 years (1/day)") # Table S3: delta>=65 = 0.80 1/day (RSE 55.60%)
    lp <- log(1.42e5)
    label("Viral production rate per productively infected cell (virions/cell/day)") # Table S3: p = 1.42 x 10^5 virus/cell/day (RSE 35.00%)

    # Viral-kinetic constants held at literature values (Methods, 'Assumptions on parameter values'; the same constants are used for the supplementary fits)
    lk <- fixed(log(4))
    label("Eclipse-to-productive transition rate k (1/day)") # Methods: k = 4 days^-1
    lc <- fixed(log(10))
    label("Virion clearance rate c (1/day)") # Methods: c = 10 days^-1
    mu <- fixed(1e-4)
    label("Fraction of produced virions that are infectious (fraction)") # Methods: mu = 10^-4
    T0 <- fixed(4e6)
    label("Initial number of target cells in the nasopharyngeal tract (cells)") # Methods: T0 = 4 x 10^6 cells
    fscale <- fixed(2.5e-6)
    label("Scaling factor from virions in the tract to copies per 10^4 sampled cells (copies/10^4 cells per virion)") # Methods: f = 10^-3 x 10^4 / (4 x 10^6) = 2.5 x 10^-6
    ltinc <- fixed(log(5))
    label("Incubation period from infection to symptom onset (day)") # Methods: prior log-normal incubation period, median 5 days

    # NEWS-2 dynamics (Table S3)
    lke0_agelt65 <- log(1)
    label("Effect-compartment rate constant ke0 = k1e, age < 65 years (1/day)") # Table S3: ke0<65 = 1 1/day (RSE 26.20%)
    lke0_agege65 <- log(0.29)
    label("Effect-compartment rate constant ke0 = k1e, age >= 65 years (1/day)") # Table S3: ke0>=65 = 0.29 1/day (RSE 24.40%)
    lhill <- log(0.98)
    label("Hill coefficient of the log10 viral load / NEWS-2 relationship (unitless)") # Table S3: n = 0.98 (RSE 10.40%)
    llvl50 <- log(4.43)
    label("Effect-compartment log10 viral load giving half-maximal NEWS-2 increase LVL50 (log10 copies/10^4 cells)") # Table S3: LVL50 = 4.43 log10 RNA copies/10^4 cells (RSE 5.89%)
    logitn0 <- log(0.8 / (3 - 0.8))
    label("Logit of the pre-infection NEWS-2 baseline N0 on its 0-3 range (logit)") # Table S3: N0 = 0.8 (RSE 43.50%), logit-normal bounded 0-3; log(0.8/(3-0.8)) = -1.0116

    # Between-subject variability. Table S3 reports SDs of the random effects;
    # the variances below are their squares.
    etalR0 ~ fixed(0.25) # Table S3: SD of the R0 random effect 0.50, held at 50% per Methods; 0.50^2
    etaldelta ~ 0.1024 # Table S3: SD 0.32 (RSE 14.01%), shared by delta<65 and delta>=65; 0.32^2
    etalp ~ 2.89 # Table S3: SD 1.70 (RSE 13.60%); 1.70^2
    etalke0 ~ 2.5921 # Table S3: SD 1.61 (RSE 13.00%), shared by ke0<65 and ke0>=65; 1.61^2
    etalhill ~ 0.0784 # Table S3: SD 0.28 (RSE 33.70%); 0.28^2
    etallvl50 ~ 0.5184 # Table S3: SD 0.72 (RSE 10.70%); 0.72^2
    etalogitn0 ~ 11.0889 # Table S3: SD 3.33 (RSE 32.00%) on the logit scale; 3.33^2
    etaltinc ~ fixed(0.015625) # Methods: incubation-period prior SD 0.125 (90% of patients 4-6 days); 0.125^2

    # Residual error (Table S3)
    addSd_log10_viral_load <- 1.34
    label("Additive residual SD on log10 viral load (log10 copies/10^4 cells)") # Table S3: sigma = 1.34 log10 RNA copies/10^4 cells (RSE 3.87%)
    addSd_news2 <- 1.56
    label("Additive residual SD on NEWS-2 (score points)") # Table S3: additive error (NEWS-2) 1.56 (RSE 3.97%)
    propSd_news2 <- 0.04
    label("Proportional residual SD on NEWS-2 (fraction)") # Table S3: proportional error (NEWS-2) 0.04 (RSE 18.40%)
  })

  model({
    # Individual parameters; age >= 65 years selects the second typical value
    # of delta and ke0 (Table S3), with one random effect shared by both strata.
    R0 <- exp(lR0 + etalR0)
    delta <- exp(ldelta_agelt65 * (1 - AGE_GE65) + ldelta_agege65 * AGE_GE65 + etaldelta)
    p <- exp(lp + etalp)
    k <- exp(lk)
    c <- exp(lc)
    ke0 <- exp(lke0_agelt65 * (1 - AGE_GE65) + lke0_agege65 * AGE_GE65 + etalke0)
    hill <- exp(lhill + etalhill)
    lvl50 <- exp(llvl50 + etallvl50)
    # Logit-normal N0 bounded to (0, 3), as in the main model (Table 2 footnote)
    logit_n0_i <- logitn0 + etalogitn0
    n0 <- 3 * expit(logit_n0_i)
    tinc <- exp(ltinc + etaltinc)

    # Infection rate constant from R0 = beta * p * T0 * mu / (c * delta) (Methods)
    beta <- R0 * c * delta / (p * T0 * mu)

    # Equation 1 (no antiviral treatment in the SoC arm)
    d/dt(target) <- -beta * virus * target
    d/dt(eclipse) <- beta * virus * target - k * eclipse
    d/dt(infected) <- k * eclipse - delta * infected
    d/dt(virus) <- mu * p * infected - c * virus
    d/dt(virus_ni) <- (1 - mu) * p * infected - c * virus_ni

    # Normalised (observed-scale) viral load, Vobs = f * (Vi + Vni)
    vl_obs <- fscale * (virus + virus_ni)

    # Equation 2 with k1e = ke0. The effect compartment is driven by the
    # normalised viral load so that log10(effect) is on the scale of LVL50
    # (log10 copies/10^4 cells); see the vignette for the evidence.
    d/dt(effect) <- ke0 * vl_obs - ke0 * effect

    # At infection: T = T0, I1 = 0, I2 = 1, Vi = Vni = 0 (Methods)
    target(0) <- T0
    infected(0) <- 1

    tsymptom <- t - tinc

    # Equation 3. log10(effect) is floored at 0 (effect below 1 copy/10^4
    # cells) so the Hill term stays in [0, 1] once the virus has cleared.
    leffect <- log10(max(effect, 1))
    news2 <- n0 + (20 - n0) * leffect^hill / (leffect^hill + lvl50^hill)

    # The 1e-10 floor only keeps the log finite at the infection time.
    log10_viral_load <- log10(vl_obs + 1e-10)

    log10_viral_load ~ add(addSd_log10_viral_load)
    news2 ~ add(addSd_news2) + prop(propSd_news2) + combined1()
  })
}
