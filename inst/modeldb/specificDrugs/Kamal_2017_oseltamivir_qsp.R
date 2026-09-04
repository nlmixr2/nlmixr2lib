Kamal_2017_oseltamivir_qsp <- function() {
  description <- paste(
    "QSP. Population-level influenza transmission model (stochastic SEIR)",
    "linked to oseltamivir pharmacology, forming the pharmacology-",
    "epidemiology half of the interdisciplinary platform of Kamal 2017.",
    "Four ODE states (susceptible, exposed, infected, recovered) over a",
    "closed population of 100 000 individuals followed for a single",
    "influenza season; no births, deaths or waning immunity. Oseltamivir",
    "acts only by shortening the duration of viral shedding Tshed, and the",
    "SEIR recovery rate is gamma = 1/Tshed, so the infected compartment",
    "empties through three parallel routes weighted by the fraction of the",
    "population that is untreated (Tshed = 6 days), treated with an",
    "oseltamivir-carboxylate daily AUC at or below the 14 180 ng.h/mL",
    "PK/PD breakpoint (Tshed = 3 days), or treated above that breakpoint",
    "(Tshed = 1.9 days). The pharmacology module is upstream and enters",
    "only through two scalars: f_uptake, the fraction of infected",
    "individuals receiving antiviral therapy, and f_auchigh, the fraction",
    "of treated patients whose OC AUC exceeds the breakpoint (0.326 for",
    "75 mg twice daily and 0.795 for 150 mg twice daily for 5 days,",
    "simulated by the authors in 5000 virtual 70-kg adults from the",
    "Kamal 2013 oseltamivir population PK model; see",
    "modellib('Kamal_2013_oseltamivir')). All five parameters carry the",
    "between-simulation variability the authors sampled over in their",
    "Monte Carlo pandemics: log-normal on each Tshed and normal on",
    "f_auchigh, so one rxode2 subject is one simulated pandemic rather",
    "than one patient. There are no dosing events; treatment scenarios are",
    "set by overriding f_uptake, f_auchigh and the infectivity rate",
    "lbeta_trans (0.21 /day moderate, 0.41 /day high transmissibility).",
    "The paper's third module, a health-economic cost-utility decision",
    "tree, is not part of this model file - it has no ODE or time",
    "dimension and its branch topology is published only as the bitmap",
    "Figure 3; the vignette reconstructs its mortality arm from Table 2.",
    sep = " "
  )
  reference <- paste(
    "Kamal MA, Smith PF, Chaiyakunapruk N, Wu DBC, Pratoomsoot C,",
    "Lee KKC, Chong HY, Nelson RE, Nieforth K, Dall G, Toovey S,",
    "Kong DCM, Kamauu A, Kirkpatrick CM, Rayner CR.",
    "Interdisciplinary pharmacometrics linking oseltamivir pharmacology,",
    "influenza epidemiology and health economics to inform antiviral use",
    "in pandemics. Br J Clin Pharmacol. 2017;83(7):1580-1594.",
    "doi:10.1111/bcp.13229.",
    "The oseltamivir carboxylate exposures behind f_auchigh come from",
    "Kamal MA, Van Wart SA, Rayner CR, Subramoney V, Reynolds DK,",
    "Bulik CC, et al. Population pharmacokinetics of oseltamivir:",
    "pediatrics through geriatrics. Antimicrob Agents Chemother.",
    "2013;57(7):3470-3477; see modellib('Kamal_2013_oseltamivir').",
    "The AUC breakpoint and the Tshed distributions come from",
    "Rayner CR, Bulik CC, Kamal MA, Reynolds DK, Toovey S, Hammel JP,",
    "et al. Pharmacokinetic-pharmacodynamic determinants of oseltamivir",
    "efficacy using data from phase 2 inoculation studies.",
    "Antimicrob Agents Chemother. 2013;57(7):3478-3487.",
    sep = " "
  )
  vignette <- "Kamal_2017_oseltamivir"

  # The four SEIR states are epidemiological population counts, not
  # pharmacokinetic compartments, so none of the canonical compartment
  # roles apply.
  paper_specific_compartments <-
    c("susceptible", "exposed", "infected", "recovered")

  units <- list(
    time          = "day",
    dosing        = paste(
      "no dosing events (population-level transmission model); the",
      "oseltamivir regimen enters through the f_uptake and f_auchigh",
      "parameters, which the authors simulated for 75 mg and 150 mg",
      "twice daily for 5 days"
    ),
    concentration = "individuals / population of 100 000"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. These states hold people, not an analyte in a
  # biological matrix, so specimen is "not applicable" throughout.
  compartmentData <- list(
    susceptible = list(
      analyte = "susceptible individuals",
      units = "individuals", specimen = "not applicable", verified = TRUE
    ),
    exposed = list(
      analyte = "exposed (latently infected, not yet infectious) individuals",
      units = "individuals", specimen = "not applicable", verified = TRUE
    ),
    infected = list(
      analyte = "infected (infectious) individuals",
      units = "individuals", specimen = "not applicable", verified = TRUE
    ),
    recovered = list(
      analyte = "recovered and immune individuals",
      units = "individuals", specimen = "not applicable", verified = TRUE
    )
  )

  # The published SEIR model carries no subject-level covariates: the
  # intervention scenario is set through the f_uptake / f_auchigh /
  # lbeta_trans parameters rather than through a data column. The
  # oseltamivir dose is documented here because f_auchigh is
  # dose-specific and a user comparing the two published regimens must
  # know which value belongs to which dose.
  covariatesDataExcluded <- list(
    DOSE_OSELTAMIVIR_MG = list(
      description = paste(
        "Oseltamivir dose per administration (twice daily for 5 days) in",
        "the simulated treatment arm"
      ),
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Not a model covariate. The dose acts only by selecting the",
        "value of f_auchigh: Kamal 2017 Table 1 gives f_auchigh = 0.326",
        "(SD 0.048) for 75 mg twice daily and 0.795 (SD 0.095) for",
        "150 mg twice daily, both for 5 days. Set f_auchigh (and its",
        "variance) to the value matching the regimen being simulated."
      ),
      source_name        = "oseltamivir dosage regimen"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 100000,
    n_studies      = 2,
    age_range      = "18-65 years (the simulated adult population)",
    weight_median  = "70 kg",
    disease_state  = paste(
      "A simulated closed population of 100 000 individuals susceptible",
      "to pandemic influenza, seeded with a single infected index case",
      "and followed for one influenza season (365 days)."
    ),
    dose_range     = paste(
      "no treatment, oseltamivir 75 mg twice daily, or 150 mg twice",
      "daily, each for 5 days; antiviral uptake 25%, 50% or 80% of",
      "infected individuals"
    ),
    regions        = "USA (health-economic and epidemic-validation data)",
    notes          = paste(
      "This is a virtual population, not an observed cohort. The",
      "pharmacodynamic inputs (the viral shedding durations Tshed and the",
      "14 180 ng.h/mL oseltamivir carboxylate AUC breakpoint) were",
      "derived from 140 subjects in two phase 2 influenza inoculation",
      "studies (Rayner 2013). The exposure inputs (f_auchigh) were",
      "simulated in 5000 virtual 70-kg adults aged 18-65 years with",
      "normal renal function using the Kamal 2013 oseltamivir population",
      "PK model, which was itself built on 390 subjects aged 1-78 years",
      "over a 20-1000 mg dose range. The SEIR structure was validated",
      "against the 2007-2008 Midwestern USA influenza season (20 263",
      "patients tested, 4970 confirmed), for which the authors report a",
      "fit of beta = 0.73 /day, 1/gamma = 4.1 days and 1/kappa = 1 day",
      "(Kamal 2017 Figure 2B)."
    )
  )

  ini({
    # ---------------------------------------------------------------
    # Transmission parameters (Kamal 2017 Table 1). Every value in this
    # model is a published simulation input, not a maximum-likelihood
    # estimate, so all are fixed().
    # ---------------------------------------------------------------

    # Table 1 gives two infectivity scenarios: 0.21 /day ("moderate
    # infectivity", the default here) and 0.41 /day ("high infectivity").
    # Override lbeta_trans with log(0.41) to run the high-transmissibility
    # scenario. Named _trans to keep it distinct from the within-host
    # target-cell infection rate lbeta_inf of Kamal_2015_oseltamivir.
    lbeta_trans <- fixed(log(0.21))
    label("Population infectivity rate beta, moderate-transmissibility scenario (1/day)")  # Table 1: beta, moderate infectivity = 0.21 /day; high infectivity = 0.41 /day

    # Table 1 reports the latency period 1/kappa = 1 day, so kappa = 1 /day.
    lkappa <- fixed(log(1))
    label("Rate of progression from exposed to infectious kappa (1/day)")  # Table 1: latency period 1/kappa = 1 day

    # ---------------------------------------------------------------
    # Duration of viral shedding by oseltamivir exposure group
    # (Kamal 2017 Table 1). The SEIR recovery rate is gamma = 1/Tshed
    # (Methods, Epidemiology module), so these three durations are what
    # oseltamivir acts on.
    # ---------------------------------------------------------------
    ltshed0 <- fixed(log(6))
    label("Duration of viral shedding without antiviral treatment (days)")  # Table 1: Tshed(0), no treatment = 6 (SD 2.5) days

    ltshed_low <- fixed(log(3))
    label("Duration of viral shedding, oseltamivir carboxylate AUC at or below 14180 ng.h/mL (days)")  # Table 1: Tshed(low) = 3 (SD 0.58) days

    ltshed_high <- fixed(log(1.9))
    label("Duration of viral shedding, oseltamivir carboxylate AUC above 14180 ng.h/mL (days)")  # Table 1: Tshed(high) = 1.9 (SD 0.51) days

    # ---------------------------------------------------------------
    # Intervention scenario inputs. Defaults reproduce the untreated
    # arm; override per scenario (Kamal 2017 Table 3 varies uptake over
    # 25%, 50% and 80% and dose over 75 mg and 150 mg twice daily).
    # ---------------------------------------------------------------
    f_uptake <- fixed(0)
    label("Fraction of infected individuals receiving oseltamivir (unitless)")  # Methods, Epidemiology module: drug uptake simulated at 25%, 50% and 80%; 0 = the no-treatment arm

    # Dose-specific; see covariatesDataExcluded$DOSE_OSELTAMIVIR_MG.
    f_auchigh <- fixed(0.326)
    label("Fraction of treated patients with oseltamivir carboxylate AUC above 14180 ng.h/mL, 75 mg twice daily (unitless)")  # Table 1: FAUChigh (75 mg BID) = 0.326 (SD 0.048); FAUChigh (150 mg BID) = 0.795 (SD 0.095)

    # ---------------------------------------------------------------
    # Between-simulation variability. The authors ran 1000 Monte Carlo
    # pandemics per scenario, drawing Tshed from a log-normal
    # distribution (Methods, Pharmacology module) and f_auchigh from the
    # simulated density of the population target-attainment fraction
    # (Methods, Pharmacology module). One rxode2 subject is therefore one
    # simulated pandemic, not one patient.
    #
    # Table 1 reports these distributions as mean (SD). For the
    # log-normal Tshed terms the variance follows the standard
    # omega^2 = log(CV^2 + 1) identity; the paper does not state whether
    # its log-normal was matched on the mean or the median, and this
    # encoding treats the tabulated value as the median (see the
    # vignette's Assumptions and deviations section).
    # ---------------------------------------------------------------
    etaltshed0     ~ fixed(0.160085)  # Table 1: Tshed(0) 6 (2.5) days -> CV 41.67% -> omega^2 = log(1 + 0.4167^2)
    etaltshed_low  ~ fixed(0.036696)  # Table 1: Tshed(low) 3 (0.58) days -> CV 19.33% -> omega^2 = log(1 + 0.1933^2)
    etaltshed_high ~ fixed(0.069573)  # Table 1: Tshed(high) 1.9 (0.51) days -> CV 26.84% -> omega^2 = log(1 + 0.2684^2)

    # Additive (not log-normal) because f_auchigh is a bounded population
    # fraction reported as mean (SD) on its natural scale. Set to
    # 0.095^2 = 0.009025 for the 150 mg twice daily regimen.
    etaf_auchigh ~ fixed(0.002304)  # Table 1: FAUChigh (75 mg BID) 0.326 (0.048) -> variance 0.048^2
  })

  model({
    # Total population N (Kamal 2017 Methods, Epidemiology module:
    # "N equals the total population (assumed to be 100 000)").
    npop <- 100000

    # ---- Per-simulation (per simulated pandemic) parameters ----
    beta_trans <- exp(lbeta_trans)
    kappa      <- exp(lkappa)
    tshed0     <- exp(ltshed0 + etaltshed0)
    tshed_low  <- exp(ltshed_low + etaltshed_low)
    tshed_high <- exp(ltshed_high + etaltshed_high)

    # f_auchigh is a fraction; the clamp is a numerical guard so a draw
    # in the upper tail of the 150 mg distribution (mean 0.795, SD 0.095)
    # cannot make the low-exposure fraction negative. It is not part of
    # the published model.
    f_auchigh_i <- max(0, min(1, f_auchigh + etaf_auchigh))

    # Recovery rate is the reciprocal of the shedding duration
    # (Methods: "the recovery rate gamma in the SEIR model is inversely
    # related to Tshed (Tshed = 1/gamma)").
    gamma0     <- 1 / tshed0
    gamma_low  <- 1 / tshed_low
    gamma_high <- 1 / tshed_high

    # Population fractions entering Equations 3 and 4. F0 is the fraction
    # not receiving therapy; FAUClow and FAUChigh split the treated
    # fraction at the 14 180 ng.h/mL breakpoint.
    frac0    <- 1 - f_uptake
    frachigh <- f_uptake * f_auchigh_i
    fraclow  <- f_uptake * (1 - f_auchigh_i)

    # ---- SEIR system (Kamal 2017 Equations 1-4) ----
    d/dt(susceptible) <- -(beta_trans / npop) * susceptible * infected
    d/dt(exposed) <-
      (beta_trans / npop) * susceptible * infected - exposed * kappa
    d/dt(infected) <-
      exposed * kappa -
      frac0 * gamma0 * infected -
      frachigh * gamma_high * infected -
      fraclow * gamma_low * infected
    d/dt(recovered) <-
      frac0 * gamma0 * infected +
      frachigh * gamma_high * infected +
      fraclow * gamma_low * infected

    # Initial conditions (Methods, Epidemiology module: "Initial
    # conditions for each compartment were as follows: S = 100 000;
    # E = 0; I = 1; R = 0").
    susceptible(0) <- 100000
    exposed(0)     <- 0
    infected(0)    <- 1
    recovered(0)   <- 0

    # ---- Derived outputs ----
    # Population-average recovery rate and the resulting reproductive
    # number, R0 ~ beta / gamma ~ beta * Tshed (Methods, Epidemiology
    # module).
    gammaeff    <- frac0 * gamma0 + frachigh * gamma_high + fraclow * gamma_low
    tshedeff    <- 1 / gammaeff
    reproNumber <- beta_trans * tshedeff

    # Cumulative number of individuals infected since the start of the
    # season - the quantity tabulated in Kamal 2017 Table 3 and passed to
    # the health-economics module.
    cumInfected <- exposed + infected + recovered
    attackRate  <- cumInfected / npop
  })
}
