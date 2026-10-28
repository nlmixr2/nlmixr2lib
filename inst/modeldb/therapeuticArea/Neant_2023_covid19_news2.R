Neant_2023_covid19_news2 <- function() {
  description <- paste(
    "Joint SARS-CoV-2 viral-kinetic / clinical-score model of Neant et al.",
    "(2023), fitted to nasopharyngeal viral load and the National Early",
    "Warning Score 2 (NEWS-2) of 664 hospitalised COVID-19 patients of the",
    "DisCoVeRy trial randomised to standard of care (SoC) or SoC plus",
    "remdesivir (final model, Table 2). Viral load follows a target-cell",
    "limited model with an eclipse phase: target cells 'target' are infected",
    "by infectious virus 'virus' at rate beta into eclipse-phase cells",
    "'eclipse', which become productively infected cells 'infected' at rate k",
    "and die at rate delta; infected cells produce virions at rate p, of",
    "which a fraction mu is infectious ('virus') and the rest noninfectious",
    "('virus_ni'), both cleared at rate c. beta is derived from the basic",
    "reproductive number R0 = beta * p * T0 * mu / (c * delta). Time is days",
    "since INFECTION (one productively infected cell, T0 = 4e6 target cells);",
    "the individual incubation period tinc maps it to days since symptom",
    "onset (output tsymptom = time - tinc). Remdesivir reduces viral",
    "production by the fraction eps while TRT_REMDESIVIR = 1. The normalised",
    "viral load (f * (virus + virus_ni), copies per 10^4 cells) drives an",
    "effect compartment 'effect' with k1e = ke0, and NEWS-2 is a Hill",
    "function of log10(effect) between the pre-infection baseline N0 (logit-",
    "normal, bounded 0-3) and the score maximum 20. delta and ke0 take",
    "separate typical values for patients aged >= 65 years (AGE_GE65)."
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
    dosing = "no dosing; remdesivir enters through the TRT_REMDESIVIR on-treatment indicator",
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
    ),
    TRT_REMDESIVIR = list(
      description = "Remdesivir on-treatment indicator (1 = receiving remdesivir, 0 = not)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no remdesivir; standard of care)",
      notes = paste(
        "Time-varying in model time (days since infection): 0 before the",
        "first remdesivir infusion, 1 from it onward while treatment is",
        "active. Multiplies viral production p by (1 - eps) (Equation 4).",
        "In DisCoVeRy remdesivir started at randomisation (200 mg on Day 1,",
        "then 100 mg once daily for up to 10 days). Set to 0 throughout for",
        "standard-of-care-only subjects."
      ),
      source_name = "SoC + remdesivir arm"
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
    n_subjects = 664L,
    n_studies = 1L,
    age_range = "Adults; median 64 years (IQR 53-72) in SoC and 63 years (IQR 55-73) in SoC + remdesivir",
    age_median = "64 years (SoC), 63 years (SoC + remdesivir)",
    weight_range = "Not reported",
    sex_female_pct = 31.3,
    race_ethnicity = "Not reported",
    disease_state = paste(
      "Hospitalised adults with laboratory-confirmed SARS-CoV-2 infection",
      "(COVID-19), randomised at most 20 days after symptom onset (median 9",
      "days, IQR 7-11); NEWS-2 >= 7 at randomisation in 74.5% (SoC) and 67.7%",
      "(remdesivir) of patients. Data collected before variants of concern",
      "and before vaccination."
    ),
    dose_range = "SoC alone (N = 329) or SoC + IV remdesivir 200 mg on Day 1 then 100 mg once daily for up to 10 days (N = 335)",
    regions = "France, Belgium, Portugal, Austria, Luxembourg (48 sites)",
    notes = paste(
      "DisCoVeRy trial (NCT04315948; EudraCT 2020-000936-23), March 2020 to",
      "January 2021. Of 833 patients in the primary modified intention-to-",
      "treat analysis, 684 had at least one nasopharyngeal viral load; 17",
      "randomised more than 20 days after symptom onset and 3 without NEWS-2",
      "data were excluded (Figure S1). 456 (68.7%) male; 349 (52.5%) younger",
      "than 65 years. Viral load (limit of detection 1 log10 copies/10^4",
      "cells; censored in the fit) measured at randomisation and Days 3, 5, 8,",
      "11, 15 and 29; NEWS-2 daily from randomisation to discharge and on",
      "Days 15 and 29. Baseline characteristics: Table 1."
    )
  )

  ini({
    # Viral-kinetic parameters (Table 2; log-normal IIV)
    lR0 <- log(6.71)
    label("Basic reproductive number R0 (unitless)") # Table 2: R0 = 6.71 (RSE 5.6%)
    ldelta_agelt65 <- log(0.84)
    label("Loss rate of productively infected cells, age < 65 years (1/day)") # Table 2: delta<65 = 0.84 1/day (RSE 4.3%)
    ldelta_agege65 <- log(0.74)
    label("Loss rate of productively infected cells, age >= 65 years (1/day)") # Table 2: delta>=65 = 0.74 1/day (RSE 31.6%)
    lp <- log(1.56e5)
    label("Viral production rate per productively infected cell (virions/cell/day)") # Table 2: p = 1.56 x 10^5 virus/cell/day (RSE 22.4%)

    # Viral-kinetic constants held at literature values (Methods, 'Assumptions on parameter values')
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

    # NEWS-2 dynamics (Table 2)
    lke0_agelt65 <- log(0.82)
    label("Effect-compartment rate constant ke0 = k1e, age < 65 years (1/day)") # Table 2: ke0<65 = 0.82 1/day (RSE 17.7%)
    lke0_agege65 <- log(0.37)
    label("Effect-compartment rate constant ke0 = k1e, age >= 65 years (1/day)") # Table 2: ke0>=65 = 0.37 1/day (RSE 21.4%)
    lhill <- log(1)
    label("Hill coefficient of the log10 viral load / NEWS-2 relationship (unitless)") # Table 2: n = 1 (RSE 8.4%)
    llvl50 <- log(4.81)
    label("Effect-compartment log10 viral load giving half-maximal NEWS-2 increase LVL50 (log10 copies/10^4 cells)") # Table 2: LVL50 = 4.81 log10 RNA copies/10^4 cells (RSE 4.0%)
    logitn0 <- log(0.86 / (3 - 0.86))
    label("Logit of the pre-infection NEWS-2 baseline N0 on its 0-3 range (logit)") # Table 2: N0 = 0.86 (RSE 26.1%), logit-normal bounded 0-3; log(0.86/(3-0.86)) = -0.9116

    # Remdesivir antiviral effect (Table 2; logit-normal)
    logiteps <- log(0.22 / (1 - 0.22))
    label("Logit of remdesivir fractional inhibition of viral production eps (logit)") # Table 2: eps = 22% (RSE 67.1%); log(0.22/0.78) = -1.2657

    # Between-subject variability. Table 2 reports SDs of the random effects;
    # the variances below are their squares.
    etalR0 ~ fixed(0.25) # Table 2: SD of the R0 random effect 0.50, held at 50% per Methods; 0.50^2
    etaldelta ~ 0.09 # Table 2: SD 0.30 (RSE 9.9%), shared by delta<65 and delta>=65; 0.30^2
    etalp ~ 2.4336 # Table 2: SD 1.56 (RSE 9.4%); 1.56^2
    etalke0 ~ 2.0164 # Table 2: SD 1.42 (RSE 9.9%), shared by ke0<65 and ke0>=65; 1.42^2
    etalhill ~ 0.0529 # Table 2: SD 0.23 (RSE 34.6%); 0.23^2
    etallvl50 ~ 0.5184 # Table 2: SD 0.72 (RSE 8.6%); 0.72^2
    etalogitn0 ~ 11.2225 # Table 2: SD 3.35 (RSE 25.1%) on the logit scale; 3.35^2
    etalogiteps ~ 5.9049 # Table 2: SD 2.43 (RSE 29.8%) on the logit scale; 2.43^2
    etaltinc ~ fixed(0.015625) # Methods: incubation-period prior SD 0.125 (90% of patients 4-6 days); 0.125^2

    # Residual error (Table 2)
    addSd_log10_viral_load <- 1.23
    label("Additive residual SD on log10 viral load (log10 copies/10^4 cells)") # Table 2: additive error 1.23 log10 RNA copies/10^4 cells (RSE 2.5%)
    addSd_news2 <- 1.53
    label("Additive residual SD on NEWS-2 (score points)") # Table 2: additive error (NEWS-2) 1.53 (RSE 2.8%)
    propSd_news2 <- 0.04
    label("Proportional residual SD on NEWS-2 (fraction)") # Table 2: proportional error (NEWS-2) 0.04 (RSE 10.5%)
  })

  model({
    # Individual parameters; age >= 65 years selects the second typical value
    # of delta and ke0 (Table 2), with one random effect shared by both strata.
    R0 <- exp(lR0 + etalR0)
    delta <- exp(ldelta_agelt65 * (1 - AGE_GE65) + ldelta_agege65 * AGE_GE65 + etaldelta)
    p <- exp(lp + etalp)
    k <- exp(lk)
    c <- exp(lc)
    ke0 <- exp(lke0_agelt65 * (1 - AGE_GE65) + lke0_agege65 * AGE_GE65 + etalke0)
    hill <- exp(lhill + etalhill)
    lvl50 <- exp(llvl50 + etallvl50)
    # Logit-normal parameters: N0 bounded to (0, 3), eps to (0, 1) (Table 2 footnote)
    logit_n0_i <- logitn0 + etalogitn0
    n0 <- 3 * expit(logit_n0_i)
    logit_eps_i <- logiteps + etalogiteps
    eps <- expit(logit_eps_i)
    tinc <- exp(ltinc + etaltinc)

    # Infection rate constant from R0 = beta * p * T0 * mu / (c * delta) (Methods)
    beta <- R0 * c * delta / (p * T0 * mu)

    # Remdesivir reduces viral production by the fraction eps (Equation 4)
    p_eff <- p * (1 - eps * TRT_REMDESIVIR)

    # Equation 1 with Equation 4 production terms
    d/dt(target) <- -beta * virus * target
    d/dt(eclipse) <- beta * virus * target - k * eclipse
    d/dt(infected) <- k * eclipse - delta * infected
    d/dt(virus) <- mu * p_eff * infected - c * virus
    d/dt(virus_ni) <- (1 - mu) * p_eff * infected - c * virus_ni

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
