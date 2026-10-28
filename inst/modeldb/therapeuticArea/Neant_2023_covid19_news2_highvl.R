Neant_2023_covid19_news2_highvl <- function() {
  description <- paste(
    "Joint SARS-CoV-2 viral-kinetic / clinical-score model of Neant et al.",
    "(2023), subgroup fit: the final model re-estimated on the DisCoVeRy",
    "patients with a nasopharyngeal viral load >= 3.5 log10 copies/10^4",
    "cells at admission",
    "(Table S1), standard of care with or without remdesivir. Viral load",
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
      "(COVID-19) whose nasopharyngeal viral load at admission was at or above",
      "3.5 log10 copies/10^4 cells, the threshold the authors associate with",
      "a positive viral culture."
    ),
    dose_range = "SoC alone or SoC + IV remdesivir 200 mg on Day 1 then 100 mg once daily for up to 10 days",
    regions = "France, Belgium, Portugal, Austria, Luxembourg (48 sites)",
    notes = paste(
      "Subgroup of the 664-patient DisCoVeRy analysis set (NCT04315948).",
      "The subgroup size is not reported; Table 1 lists 183 (91 SoC, 92 remdesivir)",
      "patients with an observed admission viral load in this band, and",
      "39.8% (SoC) / 34.0% (remdesivir) of patients had no admission viral",
      "load. Parameter estimates: Table S1."
    )
  )

  ini({
    # Viral-kinetic parameters (Table S1; log-normal IIV)
    lR0 <- log(7.31)
    label("Basic reproductive number R0 (unitless)") # Table S1: R0 = 7.31 (RSE 16.9%)
    ldelta_agelt65 <- log(0.75)
    label("Loss rate of productively infected cells, age < 65 years (1/day)") # Table S1: delta<65 = 0.75 1/day (RSE 7.56%)
    ldelta_agege65 <- log(0.70)
    label("Loss rate of productively infected cells, age >= 65 years (1/day)") # Table S1: delta>=65 = 0.70 1/day (RSE 92%)
    lp <- log(3.06e6)
    label("Viral production rate per productively infected cell (virions/cell/day)") # Table S1: p = 3.06 x 10^6 virus/cell/day (RSE 16.10%)

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

    # NEWS-2 dynamics (Table S1)
    lke0_agelt65 <- log(0.52)
    label("Effect-compartment rate constant ke0 = k1e, age < 65 years (1/day)") # Table S1: ke0<65 = 0.52 1/day (RSE 26.50%)
    lke0_agege65 <- log(0.33)
    label("Effect-compartment rate constant ke0 = k1e, age >= 65 years (1/day)") # Table S1: ke0>=65 = 0.33 1/day (RSE 64.50%)
    lhill <- log(1.44)
    label("Hill coefficient of the log10 viral load / NEWS-2 relationship (unitless)") # Table S1: n = 1.44 (RSE 17.90%)
    llvl50 <- log(6.07)
    label("Effect-compartment log10 viral load giving half-maximal NEWS-2 increase LVL50 (log10 copies/10^4 cells)") # Table S1: LVL50 = 6.07 log10 RNA copies/10^4 cells (RSE 5.72%)
    logitn0 <- log(1.5 / (3 - 1.5))
    label("Logit of the pre-infection NEWS-2 baseline N0 on its 0-3 range (logit)") # Table S1: N0 = 1.5 (RSE 27.50%), logit-normal bounded 0-3; log(1.5/(3-1.5)) = 0.0000

    # Remdesivir antiviral effect (Table S1; logit-normal)
    logiteps <- log(0.78 / (1 - 0.78))
    label("Logit of remdesivir fractional inhibition of viral production eps (logit)") # Table S1: eps = 78% (RSE 22.10%); log(0.78/0.22) = 1.2657

    # Between-subject variability. Table S1 reports SDs of the random effects;
    # the variances below are their squares.
    etalR0 ~ fixed(0.25) # Table S1: SD of the R0 random effect 0.50, held at 50% per Methods; 0.50^2
    etaldelta ~ 0.0289 # Table S1: SD 0.17 (RSE 24.70%), shared by delta<65 and delta>=65; 0.17^2
    etalp ~ 0.0196 # Table S1: SD 0.14 (RSE 131%); 0.14^2
    etalke0 ~ 1.7956 # Table S1: SD 1.34 (RSE 17.50%), shared by ke0<65 and ke0>=65; 1.34^2
    etalhill ~ 0.1156 # Table S1: SD 0.34 (RSE 44.20%); 0.34^2
    etallvl50 ~ 0.2704 # Table S1: SD 0.52 (RSE 18.10%); 0.52^2
    etalogitn0 ~ 2.6569 # Table S1: SD 1.63 (RSE 87.60%) on the logit scale; 1.63^2
    etalogiteps ~ 3.2761 # Table S1: SD 1.81 (RSE 113%) on the logit scale; 1.81^2
    etaltinc ~ fixed(0.015625) # Methods: incubation-period prior SD 0.125 (90% of patients 4-6 days); 0.125^2

    # Residual error (Table S1)
    addSd_log10_viral_load <- 1.19
    label("Additive residual SD on log10 viral load (log10 copies/10^4 cells)") # Table S1: sigma = 1.19 log10 RNA copies/10^4 cells (RSE 5.35%)
    addSd_news2 <- 1.47
    label("Additive residual SD on NEWS-2 (score points)") # Table S1: additive error (NEWS-2) 1.47 (RSE 6.77%)
    propSd_news2 <- 0.05
    label("Proportional residual SD on NEWS-2 (fraction)") # Table S1: proportional error (NEWS-2) 0.05 (RSE 18.20%)
  })

  model({
    # Individual parameters; age >= 65 years selects the second typical value
    # of delta and ke0 (Table S1), with one random effect shared by both strata.
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
