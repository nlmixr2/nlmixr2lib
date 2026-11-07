Neant_2023_covid19_news2_lowvl <- function() {
  description <- paste(
    "Joint SARS-CoV-2 viral-kinetic / clinical-score model of Neant et al.",
    "(2023), subgroup fit: the final model re-estimated on the DisCoVeRy",
    "patients with a nasopharyngeal viral load < 3.5 log10 copies/10^4",
    "cells at admission",
    "(Table S2), standard of care with or without remdesivir. Viral load",
    "follows a target-cell limited model with an eclipse phase (states",
    "target, eclipse, infected, infectious 'virus' and noninfectious",
    "'virus_ni'; beta derived from R0 = beta * p * T0 * mu / (c * delta));",
    "time is days since INFECTION and the incubation period tinc maps it to",
    "days since symptom onset (output tsymptom). Remdesivir reduces viral",
    "production by the fraction eps while TRT_REMDESIVIR = 1. The normalised",
    "viral load drives an effect compartment with k1e = ke0, and NEWS-2 is a",
    "Hill function of log10(effect) between the pre-infection baseline N0",
    "(logit-normal, bounded 0-3) and the score maximum 20. delta and ke0 take",
    "separate typical values for patients aged >= 65 years (AGE_GE65). The",
    "main final model fitted to all patients is Neant_2023_covid19_news2."
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
    n_subjects = NA_integer_,
    n_studies = 1L,
    age_range = "Adults (whole-cohort median 64 years in SoC, 63 years in SoC + remdesivir)",
    weight_range = "Not reported",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported",
    disease_state = paste(
      "Hospitalised adults with laboratory-confirmed SARS-CoV-2 infection",
      "(COVID-19) whose nasopharyngeal viral load at admission was below",
      "3.5 log10 copies/10^4 cells, the threshold the authors associate with",
      "a positive viral culture."
    ),
    dose_range = "SoC alone or SoC + IV remdesivir 200 mg on Day 1 then 100 mg once daily for up to 10 days",
    regions = "France, Belgium, Portugal, Austria, Luxembourg (48 sites)",
    notes = paste(
      "Subgroup of the 664-patient DisCoVeRy analysis set (NCT04315948).",
      "The subgroup size is not reported; Table 1 lists 236 (107 SoC, 129 remdesivir)",
      "patients with an observed admission viral load in this band, and",
      "39.8% (SoC) / 34.0% (remdesivir) of patients had no admission viral",
      "load. Parameter estimates: Table S2."
    )
  )

  ini({
    # Viral-kinetic parameters (Table S2; log-normal IIV)
    lR0 <- log(7.60)
    label("Basic reproductive number R0 (unitless)") # Table S2: R0 = 7.60 (RSE 13.40%)
    ldelta_agelt65 <- log(0.63)
    label("Loss rate of productively infected cells, age < 65 years (1/day)") # Table S2: delta<65 = 0.63 1/day (RSE 7.18%)
    ldelta_agege65 <- log(0.62)
    label("Loss rate of productively infected cells, age >= 65 years (1/day)") # Table S2: delta>=65 = 0.62 1/day (RSE 100%)
    lp <- log(2.55e3)
    label("Viral production rate per productively infected cell (virions/cell/day)") # Table S2: p = 2.55 x 10^3 virus/cell/day (RSE 29.70%)

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

    # NEWS-2 dynamics (Table S2)
    lke0_agelt65 <- log(1.1)
    label("Effect-compartment rate constant ke0 = k1e, age < 65 years (1/day)") # Table S2: ke0<65 = 1.1 1/day (RSE 29.60%); printed as a Word auto-numbered list label, see the vignette
    lke0_agege65 <- log(0.29)
    label("Effect-compartment rate constant ke0 = k1e, age >= 65 years (1/day)") # Table S2: ke0>=65 = 0.29 1/day (RSE 25.30%)
    lhill <- log(0.94)
    label("Hill coefficient of the log10 viral load / NEWS-2 relationship (unitless)") # Table S2: n = 0.94 (RSE 12.70%)
    llvl50 <- log(3.33)
    label("Effect-compartment log10 viral load giving half-maximal NEWS-2 increase LVL50 (log10 copies/10^4 cells)") # Table S2: LVL50 = 3.33 log10 RNA copies/10^4 cells (RSE 8.21%)
    logitn0 <- log(1.15 / (3 - 1.15))
    label("Logit of the pre-infection NEWS-2 baseline N0 on its 0-3 range (logit)") # Table S2: N0 = 1.15 (RSE 56%), logit-normal bounded 0-3; log(1.15/(3-1.15)) = -0.4754

    # Remdesivir antiviral effect (Table S2; logit-normal)
    logiteps <- log(0.15 / (1 - 0.15))
    label("Logit of remdesivir fractional inhibition of viral production eps (logit)") # Table S2: eps = 15% (RSE 197%); log(0.15/0.85) = -1.7346

    # Between-subject variability. Table S2 reports SDs of the random effects;
    # the variances below are their squares.
    etalR0 ~ fixed(0.25) # Table S2: SD of the R0 random effect 0.50, held at 50% per Methods; 0.50^2
    etaldelta ~ 0.0324 # Table S2: SD 0.18 (RSE 41.90%), shared by delta<65 and delta>=65; 0.18^2
    etalp ~ 1.8225 # Table S2: SD 1.35 (RSE 19.30%); 1.35^2
    etalke0 ~ 2.7225 # Table S2: SD 1.65 (RSE 13.80%), shared by ke0<65 and ke0>=65; 1.65^2
    etalhill ~ 1.3225 # Table S2: SD 1.15 (RSE 224%); 1.15^2
    etallvl50 ~ 0.5929 # Table S2: SD 0.77 (RSE 15.70%); 0.77^2
    etalogitn0 ~ 14.2884 # Table S2: SD 3.78 (RSE 29.80%) on the logit scale; 3.78^2
    etalogiteps ~ 0.64 # Table S2: SD 0.8 (RSE 128%) on the logit scale; 0.8^2
    etaltinc ~ fixed(0.015625) # Methods: incubation-period prior SD 0.125 (90% of patients 4-6 days); 0.125^2

    # Residual error (Table S2)
    addSd_log10_viral_load <- 1.02
    label("Additive residual SD on log10 viral load (log10 copies/10^4 cells)") # Table S2: sigma = 1.02 log10 RNA copies/10^4 cells (RSE 4.08%)
    addSd_news2 <- 1.59
    label("Additive residual SD on NEWS-2 (score points)") # Table S2: additive error (NEWS-2) 1.59 (RSE 4.29%)
    propSd_news2 <- 0.03
    label("Proportional residual SD on NEWS-2 (fraction)") # Table S2: proportional error (NEWS-2) 0.03 (RSE 25%)
  })

  model({
    # Individual parameters; age >= 65 years selects the second typical value
    # of delta and ke0 (Table S2), with one random effect shared by both strata.
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
