HerediaOrtiz_2014_benzopyrene <- function() {
  description <- paste(
    "Toxicokinetic (single-compartment, whole-body) model for occupational",
    "benzo(a)pyrene (BaP) exposure, tracking urinary excretion of the",
    "3-hydroxybenzo(a)pyrene (3-OHBaP) biomarker. Absorbed BaP enters a whole-body",
    "pool through route-specific first-order absorption (inhalation and dermal",
    "routes were parameterised for PAH-exposed workers; the published general",
    "equations also carry oral and IV routes, which the paper did not fit); the",
    "body pool is eliminated first-order at rate kb, with a route-specific",
    "fraction alpha recovered in urine as 3-OHBaP and the remainder lost by other",
    "routes. Human parameters were heuristically extrapolated from an in-vivo rat",
    "toxicokinetic fit. Deterministic: the publication reports no inter-individual",
    "or residual variability. This is the simple compartmental toxicokinetic model",
    "of the paper; the companion whole-body PBPK model (Figure 1) inherits human",
    "organ volumes and blood flows from Brown 1997 and Davies & Morris 1993 and",
    "its structure from the rat PBPK model of Heredia-Ortiz & Bouchard 2013, none",
    "of which are reproduced here, so that model is not packaged.",
    sep = " "
  )
  reference <- paste(
    "Heredia-Ortiz R, Maitre A, Barbeau D, Lafontaine M, Bouchard M. Use of",
    "Physiologically-Based Pharmacokinetic Modeling to Simulate the Profiles of",
    "3-Hydroxybenzo(a)pyrene in Workers Exposed to Polycyclic Aromatic",
    "Hydrocarbons. PLoS One. 2014;9(7):e102570.",
    "doi:10.1371/journal.pone.0102570. PMCID: PMC4102510. PMID: 25032692.",
    "Simple compartmental toxicokinetic model: Figure 2 and Appendix S1",
    "Eqs. 27-32; human parameter values in Results ('Single compartment",
    "toxicokinetic model').",
    sep = " "
  )
  vignette <- "HerediaOrtiz_2014_benzopyrene"
  units <- list(
    time = "h",
    dosing = "nmol (absorbed BaP dose; the model is linear, so the absolute dose scale is arbitrary and all outputs scale proportionally with it)",
    concentration = "nmol/h (urinary 3-OHBaP appearance rate, ohbap_urine_rate)"
  )

  # Absorbed BaP is administered into the two route-specific absorption
  # compartments. In the paper's exposure scenarios, D_route(t) is a time-varying
  # input over the work shift; in an event table it is a bolus or a zero-order
  # infusion (rate) into the matching depot.
  dosing <- c("depot_inh", "depot_der")

  # Paper-mechanistic states. The single whole-body pool B of Eq. 31 is split by
  # route of entry (body_inh / body_der) so that the route-specific urinary
  # recovery fractions (alpha_inh, alpha_der; both estimated) can each act on the
  # body burden that arose from their own route. The system is linear and the
  # paper simulates one route at a time, so for a single-route exposure the split
  # is identical to Eq. 31; their sum is the paper's B. Neither a route depot nor
  # a route body pool is an instance of a canonical PK compartment role, so they
  # are whitelisted here (urine is canonical).
  paper_specific_compartments <- c("depot_inh", "depot_der", "body_inh", "body_der")

  compartmentData <- list(
    depot_inh = list(
      analyte = "benzo(a)pyrene (BaP)",
      units = "nmol",
      specimen = "administration site",
      verified = TRUE
    ),
    depot_der = list(
      analyte = "benzo(a)pyrene (BaP)",
      units = "nmol",
      specimen = "administration site",
      verified = TRUE
    ),
    body_inh = list(
      analyte = "benzo(a)pyrene (BaP) whole-body burden from the inhalation route",
      units = "nmol",
      specimen = "not applicable",
      verified = TRUE
    ),
    body_der = list(
      analyte = "benzo(a)pyrene (BaP) whole-body burden from the dermal route",
      units = "nmol",
      specimen = "not applicable",
      verified = TRUE
    ),
    urine = list(
      analyte = "3-hydroxybenzo(a)pyrene (3-OHBaP)",
      units = "nmol",
      specimen = "urine",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 14L,
    n_studies = 3L,
    age_range = "28-54 years (third cohort; full range not reported)",
    sex_female_pct = NA_real_,
    disease_state = "occupationally PAH-exposed workers (healthy)",
    dose_range = "ambient BaP air concentrations ~8-9300 ng/m3 over work shifts; dermal contact with coal-tar pitch",
    regions = "France",
    notes = paste(
      "Human metabolism/absorption/elimination constants were fit to the urinary",
      "3-OHBaP elimination time course of one worker (worker 4 of an artificial",
      "shooting-target factory, over the 35-45 h after onset of a weekly",
      "exposure) and then evaluated against three further industrial cohorts:",
      "five subjects in an artificial shooting-target factory (two occupational",
      "health practitioners + three workers), five workers in a carbon disk-brake",
      "production plant, and four workers repairing metallurgical furnaces in a",
      "silicon production plant (aged 28-54 y, two of four smokers). Workers wore",
      "varying respiratory and glove protection. Urinary 3-OHBaP was measured in",
      "complete voids over a work week and the following days off and corrected",
      "for creatinine. Worker ventilation rate taken as 7.98 L/min. Direct oral",
      "exposure was considered insignificant for these workers."
    )
  )

  ini({
    # ---- Route-specific first-order absorption rate constants (Results) -------
    # 'the absorption rate (ka) of BaP by the different routes-of-entry
    # (inhalation, kainh; dermal, kader)'. Oral (kaoral) and IV (kaiv) appear in
    # the general equations (Eqs. 27-30) but were not parameterised: direct oral
    # exposure was deemed insignificant for workers and IV is not an occupational
    # route, so only the inhalation and dermal arms are carried.
    ka_inh <- fixed(0.04835)
    label("Inhalation absorption rate constant of BaP (1/h)") # Results: kainh = 48.35e-3 /h
    ka_der <- fixed(0.06181)
    label("Dermal absorption rate constant of BaP (1/h)") # Results: kader = 61.81e-3 /h

    # ---- Whole-body elimination rate constant (Results) ----------------------
    # 'The elimination rate from the body (kb) was equal to 198.26e-3 h-1'. This
    # lumps metabolism of BaP to 3-OHBaP and onward clearance; body half-life is
    # log(2)/kb ~ 3.5 h.
    kb <- fixed(0.19826)
    label("Whole-body elimination rate constant of BaP as 3-OHBaP (1/h)") # Results: kb = 198.26e-3 /h

    # ---- Route-specific urinary recovery fractions (Results) -----------------
    # 'the fractions of dose (alpha) recovered in urine as 3-OHBaP following
    # inhalation and dermal exposure were found to be 3.99e-2 and 1.64e-2'. In
    # Eq. 32, alpha multiplies kb*B, so at complete elimination the cumulative
    # urinary 3-OHBaP equals alpha times the absorbed dose.
    falpha_inh <- fixed(0.0399)
    label("Fraction of absorbed inhalation dose recovered in urine as 3-OHBaP (unitless)") # Results: alpha_inh = 3.99e-2
    falpha_der <- fixed(0.0164)
    label("Fraction of absorbed dermal dose recovered in urine as 3-OHBaP (unitless)") # Results: alpha_der = 1.64e-2
  })

  model({
    # Appendix S1 Eqs. 27-32 (simple compartmental toxicokinetic model). Each
    # route depot releases absorbed BaP into its body pool at the route's ka; the
    # body pool is eliminated at kb; a route-specific fraction alpha of that
    # elimination appears in urine as 3-OHBaP (Eq. 32), the rest as 'other
    # excretion' ((1 - alpha) * kb * B in Figure 2).
    d/dt(depot_inh) <- -ka_inh * depot_inh # Eq. 30: dA_inh/dt = -ka_inh*A_inh + D_inh(t)
    d/dt(depot_der) <- -ka_der * depot_der # Eq. 28: dA_der/dt = -ka_der*A_der + D_der(t)

    d/dt(body_inh) <- ka_inh * depot_inh - kb * body_inh # inhalation share of Eq. 31
    d/dt(body_der) <- ka_der * depot_der - kb * body_der # dermal share of Eq. 31

    # Urine accumulates 3-OHBaP from each route's body burden at its own alpha.
    d/dt(urine) <-
      falpha_inh * kb * body_inh + falpha_der * kb * body_der # Eq. 32: dU/dt = alpha*kb*B

    # Observables. Total whole-body BaP burden (sum of the route shares, i.e.
    # Eq. 31's B) and the instantaneous urinary 3-OHBaP appearance rate, which is
    # proportional to the measured urinary 3-OHBaP concentration in a void.
    body <- body_inh + body_der
    ohbap_urine_rate <- falpha_inh * kb * body_inh + falpha_der * kb * body_der
  })
}
