Leroux_2018_fluconazole <- function() {
  description <- paste(
    "One-compartment population PK model of intravenous fluconazole in",
    "preterm and term neonates with suspected or proven systemic",
    "candidiasis (Leroux 2018), with linear current-weight scaling of",
    "CL and V. Typical-value structural model only: the source paper",
    "and Data S1 supplement (goodness-of-fit plots only) do not",
    "report inter-individual variability magnitudes, residual error",
    "structure, or a maturation covariate functional form, so IIV and",
    "RUV are encoded as fixed(0) and no postmenstrual / corrected",
    "gestational age effect is encoded. See vignette Errata."
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
        "Median 1.255 kg (range 0.750-4.255) at randomization in the",
        "fluconazole arm (Table 2). Enters CL and V as a linear per-kg",
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
        "The Discussion ('our data confirmed that weight (for both",
        "antifungal agents) and postmenstrual age (only for fluconazole)",
        "influence PK in preterm neonates') indicates a postmenstrual",
        "age effect on fluconazole CL was retained in the final model.",
        "The Results section in contrast states 'These models included",
        "the impact of current weight on fluconazole and micafungin",
        "clearance and the impact of CGA on micafungin clearance' (no",
        "maturation effect on fluconazole). The two passages contradict",
        "each other. The functional form and coefficient of the",
        "maturation effect are not reported anywhere in the main paper",
        "or the Data S1 supplement (which contains only goodness-of-fit",
        "plots). The maturation covariate is therefore omitted from",
        "this library implementation; see vignette Errata."
      )
    ),
    GA = list(
      description = "Gestational age at birth",
      units       = "weeks",
      type        = "continuous",
      notes       = paste(
        "Reported in Table 2 (median 28 + 2 weeks, range 24 + 1 to",
        "40 + 1 weeks). Used by the original authors to derive PAGE",
        "(via PAGE = GA + postnatal age); not retained as a covariate",
        "on fluconazole CL in the final model (see PAGE notes)."
      )
    ),
    PNA = list(
      description = "Postnatal age",
      units       = "months",
      type        = "continuous",
      notes       = paste(
        "Reported in Table 2 in days (median 13.5, range 2.0-101.0).",
        "Used by the original authors to derive PAGE; not retained as",
        "a covariate on fluconazole CL in the final model."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 18L,
    n_studies      = 1L,
    n_observations = 82L,
    age_range      = "postnatal 2.0-101.0 days (median 13.5)",
    ga_range       = "24 + 1 to 40 + 1 weeks (median 28 + 2)",
    pma_range      = "25.6-49.1 weeks at randomization (median 30.4)",
    weight_range   = "0.750-4.255 kg current weight at randomization (median 1.255); 0.640-3.960 kg birth weight (median 0.995)",
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
      "Intravenous fluconazole, 2-hour infusion. Loading dose 25",
      "mg/kg on day 1 (cohort median 25.0 mg/kg/dose, range 21.4-27.3);",
      "maintenance dose 12 mg/kg/day for CGA < 30 weeks or 20",
      "mg/kg/day for CGA >= 30 weeks (cohort median 18.6 mg/kg/day",
      "through day 5, range 9.3-21.5; reduced after day 6 to median",
      "10.52 mg/kg/day, range 4.78-20.43). Treatment duration median",
      "6 days (range 1-20)."
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
    # Structural population PK parameters (Leroux 2018 Table 3 'Fluconazole
    # group' and Abstract 'Based on 163 PK samples, the median population
    # clearance ... and volume of distribution ... for fluconazole were:
    # 0.015 [95% CI 0.008, 0.039] and 0.913'). Reported per kg of current
    # body weight; encoded as linear-in-weight typical values.
    lcl <- fixed(log(0.015));  label("Typical clearance (L/h/kg)")           # Leroux 2018 Table 3 fluconazole CL = 0.015 L/h/kg
    lvc <- fixed(log(0.913));  label("Typical central volume (L/kg)")        # Leroux 2018 Table 3 fluconazole V  = 0.913 L/kg (95% CI 0.913, 0.913 -> effectively fixed in the original NONMEM run)

    # Inter-individual variability is not reported numerically in the main
    # paper or the Data S1 supplement (which contains only goodness-of-fit
    # plots). The text notes individual CL ranged 0.008-0.042 L/h/kg
    # ('Fluconazole and micafungin clearances were highly variable in the
    # study population, with respective ranges of 0.008-0.042 L/h/kg and
    # 0.010-0.024 L/h/kg', Results 'Population pharmacokinetic results'),
    # but no omega magnitude is given. The fluconazole V 95% CI in Table
    # 3 is degenerate (0.913, 0.913), indicating V had no inter-individual
    # variability or was fixed during estimation. Per the skill policy
    # for 'unreported IIV/RUV (structural values present)', encode IIVs
    # as fixed(0) and document in vignette Errata.
    etalcl ~ fixed(0)                                                        # Leroux 2018: IIV magnitude on CL not reported -> fixed(0) per policy
    etalvc ~ fixed(0)                                                        # Leroux 2018: IIV on V absent (Table 3 95% CI degenerate) -> fixed(0)

    # Residual error structure is not reported in the main paper or the
    # supplement (the supplement shows DV vs PRED, DV vs IPRED, WRES vs
    # TIME, WRES vs PRED, and NPDE diagnostics, but does not report the
    # sigma magnitudes or the residual model form). Per the policy for
    # ambiguous / unreported residual error with structural values
    # present, encode a combined proportional + additive residual with
    # both fixed at zero. See vignette Errata.
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
