Wang_2026_oseltamivir <- function() {
  description <- paste(
    "Integrated parent-metabolite population PK and influenza viral-dynamic",
    "(drug-disease) model for oral oseltamivir in infants and toddlers",
    "(birth through 23 months) with laboratory-confirmed influenza (CASG 114).",
    "Oseltamivir is two-compartment with first-order absorption and",
    "first-order elimination, written in micro-constants (ka, ke, kcp, kpc,",
    "Vc/F); everything eliminated from the oseltamivir central compartment",
    "passes through two sequential conversion-delay compartments (shared",
    "transit rate constant tau) into a one-compartment oseltamivir",
    "carboxylate (OC) model with apparent clearance CLm/fm and volume",
    "Vm/fm. Body weight enters CLm/fm and Vm/fm, and postnatal age the",
    "transit rate constant, all as additive linear effects centred on the",
    "cohort medians (8 kg, 8 months). The plasma OC concentration drives an",
    "Imax (fixed at 1) signal that passes through three first-order delay",
    "compartments with a common transit time T1; the delayed effect epsilon",
    "multiplies a maximal virion death rate c in a logistic viral-load",
    "model, dv/dt = (c1 * (1 - v / vmax) - c * epsilon) * v, with vmax",
    "fixed at 1e9 copies/mL. Fitted in ADAPT 5 (MLEM) with log-normal",
    "inter-individual variability; residual-error magnitudes were not",
    "reported and are carried as fixed(0)."
  )
  reference <- "Wang K, Kimberlin DW, Whitley RJ, Acosta EP. Oseltamivir drug-disease modelling in infants with influenza virus infection. Br J Clin Pharmacol. 2026;92(9):3214. doi:10.1002/bcp.70623"
  vignette <- "Wang_2026_oseltamivir"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "oseltamivir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "oseltamivir", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "oseltamivir", units = "mg", specimen = "plasma", verified = TRUE),
    transit1 = list(analyte = "oseltamivir carboxylate", units = "mg", specimen = "not applicable", verified = TRUE),
    transit2 = list(analyte = "oseltamivir carboxylate", units = "mg", specimen = "not applicable", verified = TRUE),
    central_oselcarb = list(analyte = "oseltamivir carboxylate", units = "mg", specimen = "plasma", verified = TRUE),
    effect1 = list(
      analyte = "oseltamivir carboxylate antiviral signal",
      units = "(fraction)",
      specimen = "not applicable",
      verified = TRUE
    ),
    effect2 = list(
      analyte = "oseltamivir carboxylate antiviral signal",
      units = "(fraction)",
      specimen = "not applicable",
      verified = TRUE
    ),
    effect3 = list(
      analyte = "oseltamivir carboxylate antiviral signal",
      units = "(fraction)",
      specimen = "not applicable",
      verified = TRUE
    ),
    virus = list(analyte = "influenza virus RNA", units = "copies/mL", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight at the time of PK sampling",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Additive linear effect on the OC apparent clearance and volume, centred on the cohort median of 8 kg: CLm/fm = CLm1 + CLm2 * (WT - 8) and Vm/fm = Vm1 + Vm2 * (WT - 8) (Equations 12 and 13). Baseline value per subject (sampled on treatment Day 3).",
      source_name = "weight"
    ),
    PNA = list(
      description = "Postnatal (chronological) age at the time of PK sampling",
      units = "months",
      type = "continuous",
      reference_category = NULL,
      notes = "Additive linear effect on the parent-to-OC conversion-delay transit rate constant, centred on the cohort median of 8 months: tau = tau1 + tau2 * (age - 8) (Equation 14). The cohort is birth through 23 months, so chronological age is postnatal age. Mapped to PNA, as in Kamal_2014_oseltamivir for the same CASG 114 infants.",
      source_name = "age"
    ),
    INFLUENZA_VLOAD = list(
      description = "Baseline (pre-treatment, Day 0) influenza viral load in a nasopharyngeal swab",
      units = "copies/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Initial condition of the virus state, virus(0) = INFLUENZA_VLOAD, on the LINEAR copies/mL scale (the paper plots log10 copies/mL; supply 10^log10 value). The paper fits each infant's viral-load course from the observed baseline and does not report an estimated initial-condition parameter. Must be positive: log10_viral_load is log10(virus).",
      source_name = "v (Day 0 viral load)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 78L,
    n_studies = 1L,
    age_range = "0-23 months",
    age_median = "8 months",
    weight_median = "8 kg",
    disease_state = "Infants and toddlers with laboratory-confirmed influenza virus infection",
    dose_range = "Oral oseltamivir twice daily for 5 days: 30 mg (cohort IA, 12-23 months), 3.5 mg/kg (cohorts IB and IIB) or 3.0 mg/kg (cohorts IIA, III, IV and V)",
    regions = "United States (NIAID Collaborative Antiviral Study Group)",
    notes = "CASG 114 (Kimberlin 2013). 373 oseltamivir and 382 OC plasma samples from 78 infants sampled at steady state on treatment Day 3 (predose, 1 h, 2-3 h, 5-7 h, 10-12 h); 239 viral-load measurements from 62 infants over 12 days. Sex and race are not reported in this paper (see Kimberlin 2013). Methods and Results, Wang 2026."
  )

  ini({
    # Oseltamivir (parent), Table 1. ADAPT reports micro-constants.
    lka <- log(0.878); label("Oseltamivir absorption rate constant Ka (1/h)") # Table 1 Ka = 0.878 1/h (RSE 26.5%)
    lkel <- log(0.653); label("Oseltamivir elimination (conversion) rate constant Ke (1/h)") # Table 1 Ke = 0.653 1/h (RSE 21.7%)
    lvc <- log(86.8); label("Oseltamivir apparent central volume Vc/F (L)") # Table 1 Vc = 86.8 L (RSE 22.5%)
    lk12 <- log(0.226); label("Oseltamivir central-to-peripheral rate constant Kcp (1/h)") # Table 1 Kcp = 0.226 1/h (RSE 29.3%)
    lk21 <- log(0.0805); label("Oseltamivir peripheral-to-central rate constant Kpc (1/h)") # Table 1 Kpc = 0.0805 1/h (RSE 51.6%)

    # Oseltamivir carboxylate (OC) and conversion delay, Table 1 and Equations 12-14
    lktr <- log(6.15); label("Conversion-delay transit rate constant tau at 8 months (1/h)") # Table 1 tau1 = 6.15 (RSE 22.0%)
    e_pna_ktr <- 0.0949; label("Additive slope of (PNA - 8 months) on tau (1/h per month)") # Table 1 tau2 = 0.0949 (RSE 21.1%); Equation 14
    lcl_oselcarb <- log(5.93); label("OC apparent clearance CLm/fm at 8 kg (L/h)") # Table 1 CLm1 = 5.93 (RSE 6.94%)
    e_wt_cl_oselcarb <- 0.151; label("Additive slope of (WT - 8 kg) on CLm/fm (L/h per kg)") # Table 1 CLm2 = 0.151 (RSE 6.85%); Equation 12
    lvc_oselcarb <- log(68.2); label("OC apparent volume Vm/fm at 8 kg (L)") # Table 1 Vm1 = 68.2 (RSE 11.7%)
    e_wt_vc_oselcarb <- 0.182; label("Additive slope of (WT - 8 kg) on Vm/fm (L per kg)") # Table 1 Vm2 = 0.182 (RSE 5.51%); Equation 13

    # Viral dynamics, Table 1 and Equations 7-10. Rates are per day in the paper;
    # model() converts them to the model's hour time base.
    lic50 <- log(7.02); label("OC concentration at half-maximal antiviral signal IC50 (ng/mL)") # Table 1 IC50 = 7.02 ng/mL (RSE 89.0%)
    limax <- fixed(log(1)); label("Maximal antiviral signal Imax (fraction)") # Results: 'The Imax and Vmax were fixed to 1 and 10^9 copies/mL'
    ltau <- log(0.274); label("Transit time T1 of each of the three effect-delay compartments (day)") # Table 1 T1 = 0.274 d (RSE 23.3%)
    lkkill <- log(3.11); label("Maximal virion death rate c (1/day)") # Table 1 c = 3.11 (RSE 17.5%); unit printed as mL/d
    lp <- log(0.515); label("Net virion production rate over natural death c1 (1/day)") # Table 1 c1 = 0.515 (RSE 42.0%); unit printed as mL/d
    lvirus_max <- fixed(log(1e9)); label("Maximal viral load vmax (copies/mL)") # Table 1 Vmax = 10^9 copies/mL, held constant; Results 'fixed to 1 and 10^9 copies/mL'

    # Inter-individual variability: Table 1 'IIV CV%' of log-normal parameters,
    # omega^2 = log(1 + CV^2); ADAPT prints CVs only, so the matrix is diagonal.
    etalkel ~ 0.014062 # Table 1 Ke IIV CV 11.9%
    etalvc ~ 0.240151 # Table 1 Vc IIV CV 52.1%
    etalka ~ 0.194225 # Table 1 Ka IIV CV 46.3%
    etalk12 ~ 0.629193 # Table 1 Kcp IIV CV 93.6%
    etalk21 ~ 0.051112 # Table 1 Kpc IIV CV 22.9%
    etalcl_oselcarb ~ 0.137547 # Table 1 CLm IIV CV 38.4%
    etalvc_oselcarb ~ 0.188159 # Table 1 Vm IIV CV 45.5%
    etalktr ~ 0.002752 # Table 1 tau IIV CV 5.25%
    etalic50 ~ 0.338773 # Table 1 IC50 IIV CV 63.5%
    etaltau ~ 0.891998 # Table 1 T1 IIV CV 120%
    etalkkill ~ 0.363423 # Table 1 c IIV CV 66.2%
    etalp ~ 0.578424 # Table 1 c1 IIV CV 88.5%

    # Residual error: Results state 'A proportional model with fixed additive
    # error' for the concentrations but report no magnitudes; the viral-load
    # error model is not described. Carried as fixed(0).
    propSd <- fixed(0); label("Proportional residual SD, oseltamivir (fraction; not reported)") # Results: proportional with fixed additive error, values not reported
    addSd <- fixed(0); label("Additive residual SD, oseltamivir (ng/mL; not reported)") # Results: proportional with fixed additive error, values not reported
    propSd_oselcarb <- fixed(0); label("Proportional residual SD, OC (fraction; not reported)") # Results: proportional with fixed additive error, values not reported
    addSd_oselcarb <- fixed(0); label("Additive residual SD, OC (ng/mL; not reported)") # Results: proportional with fixed additive error, values not reported
    addSd_log10_viral_load <- fixed(0); label("Additive residual SD, log10 viral load (log10 copies/mL; not reported)") # Methods: viral-load error model not described
  })

  model({
    ref_wt <- 8 # Equations 12-13: median weight 8 kg
    ref_pna <- 8 # Equation 14: median age 8 months

    # Oseltamivir
    ka <- exp(lka + etalka)
    kel <- exp(lkel + etalkel)
    vc <- exp(lvc + etalvc)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)

    # Conversion delay and OC: additive linear covariate effects (Equations 12-14)
    ktr <- (exp(lktr) + e_pna_ktr * (PNA - ref_pna)) * exp(etalktr)
    cl_oselcarb <- (exp(lcl_oselcarb) + e_wt_cl_oselcarb * (WT - ref_wt)) * exp(etalcl_oselcarb)
    vc_oselcarb <- (exp(lvc_oselcarb) + e_wt_vc_oselcarb * (WT - ref_wt)) * exp(etalvc_oselcarb)
    kel_oselcarb <- cl_oselcarb / vc_oselcarb

    # Viral dynamics; per-day parameters converted to the hour time base
    ic50 <- exp(lic50 + etalic50)
    imax <- exp(limax)
    tau <- exp(ltau + etaltau) * 24
    kkill <- exp(lkkill + etalkkill) / 24
    p <- exp(lp + etalp) / 24
    virus_max <- exp(lvirus_max)

    # Equations 1-6 (A2 = depot, A1 = central, A3 = peripheral1,
    # A4-A5 = transit1-transit2, A6 = central_oselcarb)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot + k21 * peripheral1 - k12 * central - kel * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(transit1) <- kel * central - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(central_oselcarb) <- ktr * transit2 - kel_oselcarb * central_oselcarb

    # mg / L = ug/mL; x 1000 gives ng/mL, the unit of IC50
    Cc <- central / vc * 1000
    Cc_oselcarb <- central_oselcarb / vc_oselcarb * 1000

    # Equations 7-10 (D1, D2, epsilon = effect1-effect3; v = virus)
    d/dt(effect1) <- (imax * Cc_oselcarb / (ic50 + Cc_oselcarb) - effect1) / tau
    d/dt(effect2) <- (effect1 - effect2) / tau
    d/dt(effect3) <- (effect2 - effect3) / tau
    d/dt(virus) <- (p * (1 - virus / virus_max) - kkill * effect3) * virus
    virus(0) <- INFLUENZA_VLOAD

    log10_viral_load <- log10(virus)

    Cc ~ add(addSd) + prop(propSd)
    Cc_oselcarb ~ add(addSd_oselcarb) + prop(propSd_oselcarb)
    log10_viral_load ~ add(addSd_log10_viral_load)
  })
}
