Duflot_2022_morphine <- function() {
  description <- "Three-compartment population parent-metabolite PK model for morphine, morphine-3-glucuronide (M3G) and morphine-6-glucuronide (M6G) in healthy adult volunteers given a single intravenous or nebulized morphine dose, with a Savic transit-compartment absorption for the nebulized route, Savic-kernel delayed formation of each glucuronide through one transit state, and first-order glucuronide elimination (Duflot 2022)"
  reference <- "Duflot T, Pereira T, Tavolacci MP, Joannides R, Aubrun F, Lamoureux F, Lvovschi VE. Pharmacokinetic modeling of morphine and its glucuronides: Comparison of nebulization versus intravenous route in healthy volunteers. CPT Pharmacometrics Syst Pharmacol. 2022;11(1):82-93. doi:10.1002/psp4.12735"
  vignette <- "Duflot_2022_morphine"
  units <- list(time = "min", dosing = "nmol", concentration = "nmol/L")

  covariateData <- list(
    ROUTE_IV = list(
      description = "Administration route indicator: 1 = intravenous morphine bolus, 0 = nebulized morphine",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (nebulized); the source model's own reference is intravenous, so the paper's NEB = 1 indicator enters every covariate term as (1 - ROUTE_IV)",
      notes = "Per-subject (parallel-group) route. The route also selects the dose compartment: intravenous doses go to `central`, nebulized doses to `depot` (the depot's Savic transit input carries the nebulized bioavailability, and f(depot) = 0 suppresses the ordinary bolus). Duflot 2022 Methods 'Pharmacokinetic modeling' codes the route 0 = i.v. (reference) and 1 = NEB and applies it as log(theta) = log(theta_pop) + beta; covariate building (Supplementary Material S3) retained it on k13 and km3.",
      source_name = "ROUTE (IV / NEB)"
    ),
    SEXF = list(
      description = "Biological sex indicator: 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male) in the register; the source model's reference is female, so the paper's male = 1 indicator enters the km6 term as (1 - SEXF)",
      notes = "Duflot 2022 Methods 'Pharmacokinetic modeling' codes sex 0 = women (reference) and 1 = men; beta_km6_Sex(=Male) = -0.483 (Table 3) therefore multiplies (1 - SEXF).",
      source_name = "SEX (F / M)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested as a mean-centred, log-transformed continuous covariate (Methods 'Pharmacokinetic modeling') but not retained in the final model (Supplementary Material S3 covariate building)."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Tested as a mean-centred, log-transformed continuous covariate (Methods 'Pharmacokinetic modeling') but not retained in the final model (Supplementary Material S3 covariate building)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Tested as a mean-centred, log-transformed continuous covariate (Methods 'Pharmacokinetic modeling') but not retained in the final model (Supplementary Material S3 covariate building)."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "morphine", units = "nmol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "morphine", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "morphine", units = "nmol", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "morphine", units = "nmol", specimen = "tissue", verified = TRUE),
    transit1_m3g = list(
      analyte = "morphine-3-glucuronide",
      units = "nmol",
      specimen = "not applicable",
      verified = TRUE
    ),
    central_m3g = list(analyte = "morphine-3-glucuronide", units = "nmol", specimen = "plasma", verified = TRUE),
    transit1_m6g = list(
      analyte = "morphine-6-glucuronide",
      units = "nmol",
      specimen = "not applicable",
      verified = TRUE
    ),
    central_m6g = list(analyte = "morphine-6-glucuronide", units = "nmol", specimen = "plasma", verified = TRUE)
  )

  # The nebulized dose enters `depot` and the intravenous dose enters
  # `central`; declare both so the registry records the two dose targets.
  dosing <- c("depot", "central")

  population <- list(
    species = "human",
    n_subjects = 27L,
    n_studies = 1L,
    age_range = "18-60 years (inclusion criterion); median 25 [IQR 24-34] years i.v., 27 [25-50] years NEB",
    weight_range = "median 71 [IQR 62-76] kg i.v., 68 [63-75] kg NEB",
    sex_female_pct = 48.1,
    disease_state = "Healthy adult volunteers (BMI 19-29 kg/m^2) in an experimental RIII-reflex pain model",
    dose_range = "Single morphine hydrochloride dose by Dixon up-and-down titration: 1-5 mg i.v. bolus (n = 14) or 3-8 mg nebulized over 5 min (n = 13)",
    regions = "France (Rouen University Hospital)",
    notes = "Parallel-group randomized phase I trial NCT01975753 (Duflot 2022 Methods 'Study design'; Table 1 demographics). 14 i.v. and 13 NEB participants; 7 men in each arm. Morphine, M3G and M6G quantified by LC-MS/MS and modelled in nmol/L (morphine 285.34 and glucuronides 461.46 g/mol)."
  )

  ini({
    # Structural parameters. Time in minutes, amounts in nmol, volumes in L.
    # Duflot 2022 Table 3 'Population pharmacokinetic parameter estimates';
    # '(-)' in the RSE column marks a value the authors fixed (Table 3
    # footnote and Results: fixed where the Fisher-information RSE exceeded
    # 50 percent). Mlxtran structure from Supplementary Material S4.
    lfdepot <- fixed(log(0.035))
    label("Bioavailability of the nebulized dose F (fraction)") # Table 3 row F = 0.035 (-); logit-normal in Monolix, no BSV
    lka <- fixed(log(0.046))
    label("First-order absorption rate from depot to central after nebulization ka (1/min)") # Table 3 row ka = 0.046 (-)
    lktr <- fixed(log(1.23))
    label("Transit rate constant of the nebulized-dose absorption chain Ktr (1/min)") # Table 3 row ktr = 1.23 (-)
    lmtt <- log(2.35)
    label("Mean transit time of the nebulized-dose absorption chain MTT (min)") # Table 3 row MTT = 2.35 (RSE 6.83)
    lvc <- log(1.75)
    label("Morphine central volume of distribution Vc, shared by M3G and M6G (L)") # Table 3 row Vc = 1.75 (RSE 10.6)
    lk12 <- log(0.188)
    label("Transfer rate constant central to peripheral1 k12 (1/min)") # Table 3 row k12 = 0.188 (RSE 19.2)
    lk21 <- fixed(log(0.143))
    label("Transfer rate constant peripheral1 to central k21 (1/min)") # Table 3 row k21 = 0.143 (-)
    lk13 <- fixed(log(0.306))
    label("Transfer rate constant central to peripheral2 k13, intravenous reference (1/min)") # Table 3 row k13 = 0.306 (-)
    lk31 <- fixed(log(0.010))
    label("Transfer rate constant peripheral2 to central k31 (1/min)") # Table 3 row k31 = 0.010 (-)
    lktr_m3g <- log(0.642)
    label("Morphine-to-M3G transit rate constant ktr1 (1/min)") # Table 3 row ktr1 = 0.642 (RSE 3.77)
    lmtt_m3g <- fixed(log(8.16))
    label("Mean transit time of delayed M3G formation MTT1 (min)") # Table 3 row MTT 1 = 8.16 (-)
    lka_m3g <- log(0.172)
    label("Transfer rate from the M3G transit state to the M3G central compartment kam3 (1/min)") # Table 3 row kam3 = 0.172 (RSE 19.7)
    lkel_m3g <- log(0.0038)
    label("M3G first-order elimination rate constant km3, intravenous reference (1/min)") # Table 3 row km3 = 0.0038 (RSE 1.78)
    lktr_m6g <- fixed(log(0.040))
    label("Morphine-to-M6G transit rate constant ktr2 (1/min)") # Table 3 row ktr2 = 0.040 (-)
    lmtt_m6g <- log(57.5)
    label("Mean transit time of delayed M6G formation MTT2 (min)") # Table 3 row MTT 2 = 57.5 (RSE 2.72)
    lka_m6g <- fixed(log(0.186))
    label("Transfer rate from the M6G transit state to the M6G central compartment kam6 (1/min)") # Table 3 row kam6 = 0.186 (-)
    lkel_m6g <- log(0.0081)
    label("M6G first-order elimination rate constant km6, female reference (1/min)") # Table 3 row km6 = 0.0081 (RSE 1.79)

    # Covariate effects: log(theta) = log(theta_pop) + beta * indicator
    # (Methods 'Pharmacokinetic modeling'). All three are fixed (RSE '(-)').
    e_neb_k13 <- fixed(1.39)
    label("Nebulized-route effect on log k13 (unitless)") # Table 3 row k13, beta = 1.39 (-); printed as 'beta_km3_Route (=NEB)' in the k13 row, a typo for beta_k13_Route (covariate building S3: 'Route on k13')
    e_neb_kel_m3g <- fixed(-1.33)
    label("Nebulized-route effect on log km3 (unitless)") # Table 3 row km3, beta_km3_Route (=NEB) = -1.33 (-)
    e_male_kel_m6g <- fixed(-0.483)
    label("Male-sex effect on log km6 (unitless)") # Table 3 row km6, beta_km6_Sex (=Male) = -0.483 (-)

    # Between-subject variability. Monolix reports omega as the SD of the
    # log-normal random effect, so each variance below is omega^2.
    # Random effects of F, ka, Ktr, k21, k31, MTT1, ktr2, km6 and kam6 were
    # removed (Results, 'values below 0.05').
    etalmtt ~ 0.0441 # Table 3 row MTT, BSV 0.21 (RSE 27.5); 0.21^2
    etalvc ~ 0.2401 # Table 3 row Vc, BSV 0.49 (RSE 15.8); 0.49^2
    etalk12 ~ 0.3721 # Table 3 row k12, BSV 0.61 (RSE 32.9); 0.61^2
    etalktr_m3g ~ 0.012996 # Table 3 row ktr1, BSV 0.114 (RSE 20.3); 0.114^2
    etalmtt_m6g ~ 0.0081 # Table 3 row MTT 2, BSV 0.090 (RSE 22.5); 0.090^2
    etalka_m3g ~ 0.887364 # Table 3 row kam3, BSV 0.942 (RSE 16); 0.942^2
    # k13 and km3 random effects: omega_k13 = 0.403 (fixed in the source,
    # RSE '(-)'), omega_km3 = 0.18 (RSE 8.35), corr_km3_k13 = -1 (RSE 0.81).
    # A correlation of exactly -1 makes the block singular, so the
    # covariance is scaled by 0.99 (correlation -0.99); the two published
    # variances are kept exactly.
    etalk13 + etalkel_m3g ~ c(
      0.162409,
      -0.99 * 0.403 * 0.18, 0.0324
    ) # Table 3 rows k13 BSV 0.403 (-), km3 BSV 0.18, corr_km3_k13 = -1

    # Residual error: proportional, Monolix 'b' coefficient per analyte.
    propSd <- 0.29
    label("Proportional residual error, morphine (fraction)") # Table 3 row 'Residual error for morphine' = 0.29 (RSE 5.01)
    propSd_m3g <- 0.18
    label("Proportional residual error, M3G (fraction)") # Table 3 row 'Residual error for M3G' = 0.18 (RSE 6.86)
    propSd_m6g <- 0.27
    label("Proportional residual error, M6G (fraction)") # Table 3 row 'Residual error for M6G' = 0.27 (RSE 9.07)
  })

  model({
    # Route and sex indicators in the source's own coding (0 = reference).
    neb <- 1 - ROUTE_IV
    male <- 1 - SEXF

    # Individual parameters
    fdepot <- exp(lfdepot)
    ka <- exp(lka)
    ktr <- exp(lktr)
    mtt <- exp(lmtt + etalmtt)
    vc <- exp(lvc + etalvc)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21)
    k13 <- exp(lk13 + etalk13 + e_neb_k13 * neb)
    k31 <- exp(lk31)
    ktr_m3g <- exp(lktr_m3g + etalktr_m3g)
    mtt_m3g <- exp(lmtt_m3g)
    ka_m3g <- exp(lka_m3g + etalka_m3g)
    kel_m3g <- exp(lkel_m3g + etalkel_m3g + e_neb_kel_m3g * neb)
    ktr_m6g <- exp(lktr_m6g)
    mtt_m6g <- exp(lmtt_m6g + etalmtt_m6g)
    ka_m6g <- exp(lka_m6g)
    kel_m6g <- exp(lkel_m6g + e_male_kel_m6g * male)

    # Number of transit compartments, n = MTT * Ktr - 1 (Monolix depot()
    # convention and Supplementary Material S4 N1 / N2). All three are
    # non-integer, so each chain is the analytical Savic gamma kernel.
    ntr <- mtt * ktr - 1
    ntr_m3g <- ktr_m3g * mtt_m3g - 1
    ntr_m6g <- ktr_m6g * mtt_m6g - 1

    # Nebulized absorption: Monolix depot(type = 2, p = F, ka, Mtt, Ktr)
    # feeds the dose through the Savic transit kernel into an absorption
    # state emptied at rate ka. podo(depot) / tad(depot) give the most recent
    # nebulized dose; f(depot) = 0 below removes the ordinary bolus so the
    # kernel alone delivers F * dose. Before any depot dose (every
    # intravenous subject) tad(depot) is undefined and the kernel would be
    # NaN, so the input is held at 0 until the first nebulized dose.
    input_depot <- 0
    if (tad(depot) > 0) {
      input_depot <- fdepot * exp(log(podo(depot)) + log(ktr) + ntr * log(ktr * tad(depot)) - ktr * tad(depot) - lgamma(ntr + 1))
    }
    d/dt(depot) <- input_depot - ka * depot

    # Morphine disposition (Supplementary Material S4 ddt_Ac, ddt_Ap1, ddt_Ap2).
    # Morphine is eliminated only through the two glucuronidation outflows
    # ktr1 * Ac and ktr2 * Ac.
    d/dt(central) <- ka * depot - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2 - ktr_m3g * central - ktr_m6g * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # Glucuronide formation, transcribed from Supplementary Material S4:
    #   ddt_AM3Gtr = exp(log(Ac) + log(ktr1) + N1*log(ktr1*t) - ktr1*t - LN1fac) - kam3*AM3Gtr
    # i.e. the current morphine amount times a Savic gamma kernel evaluated
    # on the model time t (time since the start of the record, which is the
    # single dose time in the source data), not on time after dose.
    d/dt(transit1_m3g) <- central * exp(log(ktr_m3g) + ntr_m3g * log(ktr_m3g * t) - ktr_m3g * t - lgamma(ntr_m3g + 1)) - ka_m3g * transit1_m3g
    d/dt(central_m3g) <- ka_m3g * transit1_m3g - kel_m3g * central_m3g
    d/dt(transit1_m6g) <- central * exp(log(ktr_m6g) + ntr_m6g * log(ktr_m6g * t) - ktr_m6g * t - lgamma(ntr_m6g + 1)) - ka_m6g * transit1_m6g
    d/dt(central_m6g) <- ka_m6g * transit1_m6g - kel_m6g * central_m6g

    f(depot) <- 0

    # Concentrations (nmol/L). The glucuronide volumes equal Vc (Table 3
    # footnote; S4 M3G = AM3G/V, M6G = AM6G/V).
    Cc <- central / vc
    Cc_m3g <- central_m3g / vc
    Cc_m6g <- central_m6g / vc

    Cc ~ prop(propSd)
    Cc_m3g ~ prop(propSd_m3g)
    Cc_m6g ~ prop(propSd_m6g)
  })
}
