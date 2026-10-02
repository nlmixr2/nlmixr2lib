FergusonSells_2022_tadalafil <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption for oral",
    "tadalafil in adult and pediatric patients (2 to < 18 years) with pulmonary",
    "arterial hypertension, fitted to pooled PHIRST-1 (adults) and H6D-MC-LVIG",
    "(children) data. Apparent clearance is higher in patients taking",
    "concomitant bosentan (a CYP3A inducer), apparent volume scales linearly",
    "with body weight (exponent fixed to 1, reference 70 kg), and relative",
    "bioavailability is a power function of dose (falling with increasing dose)",
    "and of age (falling with decreasing age). Clearance has no weight effect.",
    "Residual error is combined proportional plus additive."
  )

  reference <- paste(
    "Ferguson-Sells L, Velez de Mendizabal N, Li B, Small D.",
    "Population Pharmacokinetics of Tadalafil in Pediatric Patients with",
    "Pulmonary Arterial Hypertension: A Combined Adult/Pediatric Model.",
    "Clin Pharmacokinet. 2022;61(2):249-262.",
    "doi:10.1007/s40262-021-01052-8"
  )

  vignette <- "FergusonSells_2022_tadalafil"

  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "tadalafil", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tadalafil", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight at study entry",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline (study-entry) weight; the paper states within-patient weight",
        "changes were negligible over the pediatric PK visits (Results).",
        "Linear effect on V/F, V = TVV * (WT/70)^1, with the exponent fixed to",
        "the allometric value of 1 (Table 2 footnote d). An allometric exponent",
        "on CL/F was estimated at about 0 and removed, so CL/F carries no",
        "weight effect. Pooled range 10-140 kg (Discussion 4.1)."
      ),
      source_name = "WT"
    ),
    AGE = list(
      description = "Age at study entry",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on bioavailability, F x (AGE/52.4)^0.100 (Table 2",
        "footnote e; the paper's column is AGEE, age at study entry). F falls",
        "with decreasing age. The reference 52.4 years is printed in the",
        "equation; the paper does not say how it was chosen. Pooled range",
        "2.5-90.3 years (Table 1). The paper cautions the model must not be",
        "extrapolated below 2 years (no patient < 2 years enrolled)."
      ),
      source_name = "AGEE"
    ),
    DOSE = list(
      description = "Administered tadalafil dose at the dose record",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Use case (a) of the DOSE canonical: the patient's reported dose drives",
        "a power effect on bioavailability, F x (DOSE/16.27)^-0.227 (Table 2",
        "footnote e), so F FALLS as dose rises (less-than-dose-proportional",
        "exposure). Time-varying within a pediatric patient, who received a",
        "low dose for 5 weeks then a high dose for 5 weeks. The reference",
        "16.27 mg is printed in the equation; the paper does not say how it was",
        "chosen. Doses studied 2-40 mg once daily."
      ),
      source_name = "DOSE"
    ),
    CONMED_BOSENTAN = list(
      description = "Concomitant bosentan (CYP3A inducer) treatment indicator, 1 = yes",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (bosentan-taking) in this paper - see notes",
      notes = paste(
        "The canonical orientation is kept (1 = taking bosentan, 0 = not),",
        "which is also the paper's BOS coding. The paper's typical CL/F, 3.23",
        "L/h, is for patients TAKING bosentan, and the estimated coefficient is",
        "the fractional change WITHOUT bosentan: CL = TVCL * BOS + TVCL * (1 +",
        "EffNoBos) * (1 - BOS), EffNoBos = -0.418 (Table 2 footnote c). The",
        "published coefficient is kept verbatim and applied to (1 -",
        "CONMED_BOSENTAN). Patients taking ambrisentan were grouped with",
        "non-bosentan patients (Methods 2.2), so ambrisentan users take 0.",
        "About 50% of adults and children took bosentan (Discussion)."
      ),
      source_name = "BOS"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      notes = "Evaluated and not significant (Conclusion)."
    ),
    FORM_TABLET = list(
      description = "Tablet (1) versus oral suspension (0) formulation indicator",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Formulation was tested and not found significant (Discussion 4.1).",
        "Light-weight (< 25 kg) children received a 2 mg/mL oral suspension;",
        "all other patients received the commercial tablet, so formulation is",
        "confounded with age and dose."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 324L,
    n_studies = 2L,
    n_observations = 1430L,
    age_range = "2.5-90.3 years (adults 14.7-90.3; children 2.5-18.0)",
    age_median = "53.9 years (adults); 14.6 / 11.0 / 5.0 years in the heavy / middle / light-weight pediatric cohorts",
    weight_range = "10.0-140 kg",
    weight_median = "73.0 kg (adults); 49.0 / 30.1 / 14.7 kg in the heavy / middle / light-weight pediatric cohorts",
    sex_female_pct = 76.9,
    disease_state = "Pulmonary arterial hypertension (PAH)",
    dose_range = paste(
      "Adults: 2.5, 10, 20 or 40 mg tablet once daily (PHIRST-1). Children:",
      "low dose 5 weeks then high dose 5 weeks once daily -- heavy-weight",
      "(>= 40 kg) 10 then 20-40 mg tablet, middle-weight (25 to < 40 kg) 5",
      "then 10-20 mg tablet, light-weight (< 25 kg) 2-4 then 8-20 mg",
      "suspension (LVIG)"
    ),
    co_medication = paste(
      "About 50% of adults and children received concomitant bosentan;",
      "ambrisentan users were grouped with non-bosentan patients"
    ),
    regions = "Multinational (PHIRST-1 and LVIG were multicenter international studies)",
    notes = paste(
      "305 adults from PHIRST-1 (NCT00125918; 69 male, 236 female; 1102",
      "observations; sparse sampling at weeks 4, 8, 12 and 16) and 19 children",
      "from LVIG (NCT01484431; 6 male, 13 female; 328 observations; serial",
      "sampling at predose and 2, 4, 8, 12 and 24 h on day 1, day 14 and day",
      "49, plus one sparse Period 2 sample). Demographics from Table 1 and",
      "Results. PHIRST-1 enrolled one patient under 18 years (14 years), who",
      "is counted with the adults. LLOQ 0.500 ng/mL; BLQ samples excluded."
    )
  )

  ini({
    # Absorption
    lka <- log(0.860); label("First-order absorption rate constant Ka (1/h)") # Table 2: Ka = 0.860 1/h (%SEE 10.9)

    # Apparent clearance. The typical value is for patients TAKING bosentan
    # (Table 2 footnote c); patients not taking bosentan have
    # CL = 3.23 * (1 - 0.418) = 1.88 L/h, the value Table 2 prints as
    # 'Not taking bosentan (calculated)'.
    lcl <- log(3.23); label("Apparent clearance CL/F in patients taking bosentan (L/h)") # Table 2: 'Patients taking bosentan' 3.23 L/h (%SEE 4.03)
    e_conmed_bosentan_cl <- -0.418; label("Fractional change in CL/F in patients NOT taking bosentan (unitless)") # Table 2: 'Effect of non-bosentan' -0.418 (%SEE 7.87)

    # Apparent volume, linear in weight with a 70 kg reference (Table 2
    # footnote d).
    lvc <- log(88.1); label("Apparent volume of distribution V/F for a 70 kg patient (L)") # Table 2: 'V/F 70 kg patient' 88.1 L (%SEE 4.97)
    e_wt_vc <- fixed(1); label("Power exponent of body weight (/70 kg) on V/F (unitless)") # Table 2: 'Effect of weight' 1 fixed

    # Bioavailability: F = TVF * (DOSE/16.27)^EffDoseF * (AGEE/52.4)^EffAgeF
    # with TVF fixed to 1 (Table 2 footnote e).
    lfdepot <- fixed(log(1)); label("Typical relative bioavailability F at 16.27 mg and 52.4 years (fraction)") # Table 2: 'F' 1 fixed
    e_dose_fdepot <- -0.227; label("Power exponent of dose (/16.27 mg) on F (unitless)") # Table 2: 'Effect of dose (continuous) on F' -0.227 (%SEE 12.6)
    e_age_fdepot <- 0.100; label("Power exponent of age (/52.4 years) on F (unitless)") # Table 2: 'Effect of age on F' 0.100 (%SEE 45.6)

    # IIV. Table 2 footnote b reports IIV as the log-normal CV,
    # 100% * sqrt(exp(omega^2) - 1), so omega^2 = log(1 + CV^2). The CL/F-V/F
    # correlation of the base model was removed from the final model
    # (Results), so OMEGA is diagonal. No IIV on F.
    etalka ~ 1.617426 # Table 2: Ka IIV 201% CV -> log(1 + 2.01^2)
    etalcl ~ 0.211253 # Table 2: CL/F IIV 48.5% CV -> log(1 + 0.485^2)
    etalvc ~ 0.098071 # Table 2: V/F IIV 32.1% CV -> log(1 + 0.321^2)

    # Residual error. Table 2 footnote f: proportional CV = 100% * sqrt(sigma1)
    # and additive SD = sqrt(x^2 * sigma1), both on one epsilon, i.e. the
    # variance is sigma1 * (IPRED^2 + x^2) = propSd^2 * IPRED^2 + addSd^2.
    propSd <- 0.258; label("Proportional residual error SD (fraction)") # Table 2: Proportional 25.8% (%SEE 10.4)
    addSd <- 11.6; label("Additive residual error SD (ng/mL)") # Table 2: Additive 11.6 ng/mL (%SEE 31.3)
  })

  model({
    # 1. Derived covariate terms
    # Relative bioavailability (Table 2 footnote e)
    fdepot <- exp(lfdepot) * (DOSE / 16.27)^e_dose_fdepot * (AGE / 52.4)^e_age_fdepot

    # 2. Individual PK parameters
    ka <- exp(lka + etalka)
    # CL = TVCL * BOS + TVCL * (1 + EffNoBos) * (1 - BOS) (Table 2 footnote c)
    cl <- exp(lcl + etalcl) * (1 + e_conmed_bosentan_cl * (1 - CONMED_BOSENTAN))
    # V = TVV * (WT / 70)^EffWt (Table 2 footnote d)
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc

    # 3. Micro-constants
    kel <- cl / vc

    # 4. ODE system
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # 5. Bioavailability
    f(depot) <- fdepot

    # 6. Observation and error model; dose mg / volume L = mg/L = ug/mL,
    # times 1000 for ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
