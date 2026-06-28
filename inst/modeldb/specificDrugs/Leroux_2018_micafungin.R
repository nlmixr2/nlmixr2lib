Leroux_2018_micafungin <- function() {
  description <- paste(
    "One-compartment population PK model of intravenous micafungin in",
    "preterm and term neonates with suspected or proven systemic",
    "candidiasis (Leroux 2018), with linear current-weight scaling of",
    "CL and V. Typical-value structural model only: the source paper",
    "and Data S1 supplement (goodness-of-fit plots only) do not",
    "report inter-individual variability magnitudes, residual error",
    "structure, or the functional form / coefficient of the",
    "corrected-gestational-age effect on CL that the paper mentions,",
    "so IIV and RUV are encoded as fixed(0) and the CGA covariate is",
    "omitted. See vignette Errata."
  )
  reference <- paste(
    "Leroux S, Jacqz-Aigrain E, Elie V, Legrand F, Barin-Le Guellec C,",
    "Aurich B, Biran V, Dusang B, Goudjil S, Coopman S,",
    "Garcia Sanchez R, Zhao W, Manzoni P;",
    "FP7 TINN (Treat Infections in NeoNates) consortium.",
    "Pharmacokinetics and safety of fluconazole and micafungin in",
    "neonates with systemic candidiasis: a randomized, open-label",
    "clinical trial. Br J Clin Pharmacol. 2018;84(9):1989-1999.",
    "doi:10.1111/bcp.13628.",
    sep = " "
  )
  vignette <- "Leroux_2018_fluconazole_micafungin"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Current body weight at the time of dosing",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Median 1.090 kg (range 0.640-4.615) at randomization in the",
        "micafungin arm (Table 2). Enters CL and V as a linear per-kg",
        "scalar: CL_i = TVCL * WT_i and V_i = TVV * WT_i, matching the",
        "L/h/kg and L/kg units in which the population estimates are",
        "reported in Table 3. The paper text 'These models included the",
        "impact of current weight on fluconazole and micafungin",
        "clearance' (Results 'Population pharmacokinetic results') does",
        "not state an allometric exponent; the per-kg units in Table 3",
        "imply linear (exponent 1) scaling."
      ),
      source_name        = "Current weight"
    )
  )

  covariatesDataExcluded <- list(
    PAGE = list(
      description = "Postmenstrual age (gestational age at birth + postnatal age)",
      units       = "months",
      type        = "continuous",
      notes       = paste(
        "Synonymous with the paper's 'corrected gestational age' (CGA)",
        "in the neonatal context. The Results section states 'These",
        "models included ... the impact of CGA on micafungin clearance'",
        "but the functional form (linear, power, Hill maturation, ...)",
        "and the coefficient are not reported anywhere in the main",
        "paper or the Data S1 supplement (which contains only",
        "goodness-of-fit plots). The Discussion contradicts the",
        "Results: 'our data confirmed that weight (for both antifungal",
        "agents) and postmenstrual age (only for fluconazole) influence",
        "PK in preterm neonates'. The CGA / PAGE covariate is therefore",
        "omitted from this library implementation; see vignette Errata."
      )
    ),
    GA = list(
      description = "Gestational age at birth",
      units       = "weeks",
      type        = "continuous",
      notes       = paste(
        "Reported in Table 2 (median 26 + 6 weeks, range 23 + 4 to",
        "40 + 0 weeks). Used by the original authors to derive the",
        "corrected gestational age (CGA = GA + postnatal age); not",
        "encoded here because the functional form on CL is unreported",
        "(see PAGE notes)."
      )
    ),
    PNA = list(
      description = "Postnatal age",
      units       = "months",
      type        = "continuous",
      notes       = paste(
        "Reported in Table 2 in days (median 12.5, range 3.0-115.0).",
        "Used by the original authors to derive CGA; not encoded here",
        "because the maturation functional form on CL is unreported."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 18L,
    n_studies      = 1L,
    n_observations = 81L,
    age_range      = "postnatal 3.0-115.0 days (median 12.5)",
    ga_range       = "23 + 4 to 40 + 0 weeks (median 26 + 6)",
    pma_range      = "25.4-41.9 weeks at randomization (median 29.9)",
    weight_range   = "0.640-4.615 kg current weight at randomization (median 1.090); 0.500-3.630 kg birth weight (median 0.885)",
    sex_female_pct = 38.9,
    disease_state  = paste(
      "Preterm and term neonates and young infants (24-42 weeks",
      "corrected gestational age; postnatal age between 48 h and day",
      "of life 120 at culture acquisition) with suspected or",
      "microbiologically documented systemic Candida infection,",
      "treated in five French and one Spanish neonatal intensive care",
      "units between 2013 and 2015."
    ),
    dose_range     = paste(
      "Intravenous micafungin, 2-hour infusion. Loading dose 15",
      "mg/kg/day on day 1 (cohort median 15.1 mg/kg/dose, range",
      "12.1-18.8); maintenance dose 10 mg/kg/day regardless of CGA",
      "(cohort median 9.9 mg/kg/day through day 5, range 8.1-11.4;",
      "reduced after day 6 to median 4.76 mg/kg/day, range 3.48-10.4).",
      "Treatment duration median 4 days (range 1-35)."
    ),
    regions        = "France (Amiens, Lille, Paris, Saint Pierre de la Reunion) and Spain (Salamanca); FP7 TINN consortium, 2013-2015.",
    notes          = paste(
      "Randomized 1:1 to fluconazole vs micafungin; PK samples on",
      "treatment days 1 and 5 per the limited-PK schedule in Table 1",
      "(2-3 samples per occasion, 4 alternative schedules depending",
      "on CGA). Estimation: NONMEM v7.2 with Monte Carlo simulation.",
      "Baseline demographics in Table 2; population PK results in",
      "Table 3."
    )
  )

  ini({
    # Structural population PK parameters (Leroux 2018 Table 3 'Micafungin
    # group' and Abstract 'and for micafungin were: 0.020 (95% CI 0.010,
    # 0.023) and 0.354 (95% CI 0.225, 0.482)'). Reported per kg of
    # current body weight; encoded as linear-in-weight typical values.
    lcl <- fixed(log(0.020));  label("Typical clearance (L/h/kg)")           # Leroux 2018 Table 3 micafungin CL = 0.020 L/h/kg
    lvc <- fixed(log(0.354));  label("Typical central volume (L/kg)")        # Leroux 2018 Table 3 micafungin V  = 0.354 L/kg

    # Inter-individual variability is not reported numerically in the
    # main paper or the Data S1 supplement (goodness-of-fit plots only).
    # The text notes individual CL ranged 0.010-0.024 L/h/kg
    # ('Fluconazole and micafungin clearances were highly variable in
    # the study population, with respective ranges of 0.008-0.042 L/h/kg
    # and 0.010-0.024 L/h/kg', Results 'Population pharmacokinetic
    # results') and individual V ranged 0.225-0.482 L/kg, but no omega
    # magnitudes are given. Per the skill policy for 'unreported IIV/RUV
    # (structural values present)', encode IIVs as fixed(0) and document
    # in vignette Errata.
    etalcl ~ fixed(0)                                                        # Leroux 2018: IIV magnitude on CL not reported -> fixed(0) per policy
    etalvc ~ fixed(0)                                                        # Leroux 2018: IIV magnitude on V not reported -> fixed(0) per policy

    # Residual error structure is not reported in the main paper or the
    # supplement. Per the policy for ambiguous / unreported residual
    # error with structural values present, encode a combined
    # proportional + additive residual with both fixed at zero. See
    # vignette Errata.
    propSd <- fixed(0);  label("Proportional residual error (fraction)")    # Leroux 2018: residual error structure not reported -> fixed(0) per policy
    addSd  <- fixed(0);  label("Additive residual error (mg/L)")            # Leroux 2018: residual error magnitude not reported -> fixed(0) per policy
  })

  model({
    # Individual PK parameters. Linear weight scaling per the L/h/kg and
    # L/kg units in Table 3.
    cl <- exp(lcl + etalcl) * WT
    vc <- exp(lvc + etalvc) * WT

    kel <- cl / vc

    # One-compartment IV (Leroux 2018 Results 'Population pharmacokinetic
    # results': 'Two one-compartment population PK models were developed
    # and validated'). Doses are 2-hour intravenous infusions per the
    # protocol (Procedures); the library model does not hard-code the
    # infusion duration so users supply rate or dur per dose in the
    # event table.
    d/dt(central) <- -kel * central

    # Plasma concentration: dose in mg, vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
