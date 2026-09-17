Ambery_2015_batefenterol <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order absorption for",
    "inhaled GSK961081 (batefenterol) DISKUS, a bifunctional muscarinic",
    "antagonist / beta2-agonist (MABA), in patients with moderate-to-severe",
    "chronic obstructive pulmonary disease (COPD) (Ambery 2015, substudy of",
    "GSK MAB115032 / NCT01319019). All disposition parameters are apparent",
    "(CL/F, V2/F, Q, V3/F) because the inhaled bioavailability is not",
    "identifiable from inhaled-only data. Plasma concentrations were",
    "log-transformed and fitted in NONMEM 7 by Monte Carlo importance",
    "sampling with M3 handling of the below-LLOQ data; the model-building",
    "data set was restricted to the two arms with less than half their data",
    "below the 25 pg/mL LLOQ (800 ug once daily and 400 ug twice daily, day",
    "28 only). Inter-individual variability could be estimated reliably only",
    "for apparent clearance, and no covariates were retained (the authors",
    "cite Ribbing and Jonsson against covariate selection in data sets of",
    "fewer than 50 to 100 subjects)."
  )

  reference <- paste(
    "Ambery CL, Wielders P, Ludwig-Sengpiel A, Chan R, Riley JH.",
    "Population Pharmacokinetics and Pharmacodynamics of GSK961081",
    "(Batefenterol), a Muscarinic Antagonist and beta2-Agonist, in",
    "Moderate-to-Severe COPD Patients: Substudy of a Randomized Trial.",
    "Drugs R D. 2015;15(3):281-291. doi:10.1007/s40268-015-0104-x.",
    "PMID 26286203; PMCID PMC4561049.",
    sep = " "
  )

  vignette <- "Ambery_2015_batefenterol"

  # Doses are administered in ug and the apparent volumes are in L, so the
  # amount/volume ratio Cc is in ug/L, which is the same as ng/mL. Ambery
  # 2015 reports plasma concentrations in pg/mL (LLOQ 25 pg/mL); 1 ug/L =
  # 1000 pg/mL. No parameter has been rescaled -- the conversion is applied
  # in the vignette only where the simulation is compared with Fig. 2.
  units <- list(time = "h", dosing = "ug", concentration = "ug/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Ambery 2015 Sect. 3.1 ("A
  # two-compartment disposition model (ADVAN4 TRANS4) with first-order
  # absorption") and Sect. 2.3 (plasma sampling, HPLC-MS/MS assay).
  compartmentData <- list(
    depot = list(analyte = "batefenterol", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "batefenterol", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "batefenterol", units = "ug", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at screening",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Reported for the PK model-building subset in Ambery 2015 Table 1",
        "(63 +/- 8.9 years). Covariates were screened only for the",
        "pharmacodynamic model; Sect. 3.1 states 'Covariates were not",
        "included in the PK model, because of the limited data set', and the",
        "Discussion cites Ribbing and Jonsson against covariate selection in",
        "data sets of fewer than 50 to 100 subjects (the PK set had 47)."
      ),
      source_name = "Age (years)"
    ),
    WT = list(
      description = "Body weight at screening",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Reported for the PK model-building subset in Ambery 2015 Table 1",
        "(72 +/- 13 kg). Not screened on the PK model; see the AGE entry."
      ),
      source_name = "Weight (kg)"
    ),
    HT = list(
      description = "Height at screening",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Reported for the PK model-building subset in Ambery 2015 Table 1",
        "(169 +/- 7.7 cm). Not screened on the PK model; see the AGE entry."
      ),
      source_name = "Height (cm)"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = male",
      notes = paste(
        "Ambery 2015 Table 1 reports 62% male in the PK model-building",
        "subset, i.e. 38% female. Not screened on the PK model; see the",
        "AGE entry."
      ),
      source_name = "Male (%)"
    ),
    SMOKE = list(
      description = "Current-smoker indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = former smoker",
      notes = paste(
        "Ambery 2015 Table 1 reports 49% current smokers in the PK",
        "model-building subset; all patients were current or former smokers",
        "with at least a 10 pack-year history (Sect. 2.2). Not screened on",
        "the PK model; see the AGE entry."
      ),
      source_name = "Current smoker (%)"
    ),
    CONMED_ICS = list(
      description = "Concurrent inhaled-corticosteroid use",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no concurrent ICS",
      notes = paste(
        "Ambery 2015 Table 1 reports 53% concurrent ICS use in the PK",
        "model-building subset. Not screened on the PK model; see the AGE",
        "entry."
      ),
      source_name = "Concurrent ICS use (%)"
    ),
    FEV1_BL = list(
      description = "Baseline (day 1) trough FEV1",
      units = "L",
      type = "continuous",
      notes = paste(
        "Ambery 2015 Table 1 reports 1.36 +/- 0.42 L for the PK",
        "model-building subset. It is a retained covariate of the companion",
        "pharmacodynamic model (see modellib('Ambery_2015_batefenterol_fev1'))",
        "but was not screened on the PK model; see the AGE entry."
      ),
      source_name = "Baseline FEV1 on day 1 (L)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 47L,
    n_studies = 1L,
    n_observations = 405L,
    age_mean_sd = "63 +/- 8.2 years (PD / PK-PD set); 63 +/- 8.9 years (PK model-building set)",
    weight_mean_sd = "76 +/- 14 kg (PD / PK-PD set); 72 +/- 13 kg (PK model-building set)",
    height_mean_sd = "171 +/- 8.5 cm (PD / PK-PD set); 169 +/- 7.7 cm (PK model-building set)",
    sex_female_pct = 38,
    disease_state = paste(
      "Moderate-to-severe stable COPD: post-bronchodilator FEV1/FVC below",
      "70% and FEV1 30-70% of predicted by NHANES III normative values;",
      "current or former smokers aged 40 years or above with at least a",
      "10 pack-year history; diagnosed asthma excluded"
    ),
    dose_range = paste(
      "GSK961081 DISKUS 100, 400 and 800 ug once daily and 100, 200 and",
      "400 ug twice daily for 4 weeks; the PK model-building data set used",
      "only the 800 ug once-daily and 400 ug twice-daily arms on day 28"
    ),
    regimens = "Inhaled dry powder (DISKUS), once or twice daily for 28 days",
    regions = "Multicenter international (GSK MAB115032 / NCT01319019)",
    notes = paste(
      "Ambery 2015 Table 1 gives demographics for two overlapping analysis",
      "sets: the PD / PK-PD set (n = 347) and the PK model-building set",
      "(n = 47). The population block above describes the PK model-building",
      "set, which is the set this model was fitted to. More than half of the",
      "PK observations were below the 25 pg/mL LLOQ for every treatment",
      "except 800 ug once daily (days 1 and 28) and 400 ug twice daily (day",
      "28), so the model-building data set was restricted to the 800 ug",
      "once-daily and 400 ug twice-daily arms on day 28, where 30% and 27%",
      "of the data respectively were below the LLOQ and were handled by the",
      "M3 likelihood method. The full PK data set used for the visual",
      "predictive check comprised 1112 samples from 128 patients across all",
      "six GSK961081 arms. PK sampling on day 28 was after both the morning",
      "and the evening dose at 1 h to 0 min pre-dose and then 0-30 min,",
      "30 min-2 h, 2-6 h and 6-11 h post-dose."
    )
  )

  ini({
    # ================================================================
    # Structural disposition parameters. Ambery 2015 Table 2 reports
    # NONMEM THETAs directly on the LOG scale (the column header is
    # "Log estimate"), so the values below are transcribed verbatim
    # rather than back-transformed and re-logged. Back-transformed
    # values are given in each label for readability.
    #
    # The log scale is confirmed independently by the paper's own
    # below-LLOQ percentages: exp(6.85) = 945 L/h gives a steady-state
    # AUC(0-24) of 800 ug / 945 L/h = 847 pg*h/mL at a total daily dose
    # of 800 ug, i.e. an average concentration of ~35 pg/mL against a
    # 25 pg/mL LLOQ -- which reproduces Sect. 2.4's statement that the
    # 800 ug once-daily and 400 ug twice-daily day-28 arms were 30% and
    # 27% below LLOQ while every other arm exceeded 50%.
    # ================================================================
    lka <- -0.89
    label("Absorption rate constant KA (1/h)")  # Table 2, log estimate -0.890 (95% CI -1.16 to -0.625); KA = exp(-0.89) = 0.411 /h

    lcl <- 6.85
    label("Apparent clearance CL/F (L/h)")  # Table 2, log estimate 6.85 (95% CI 6.62 to 7.08); CL/F = exp(6.85) = 945 L/h

    lvc <- 6.26
    label("Apparent central volume of distribution V2/F (L)")  # Table 2, log estimate 6.26 (95% CI 5.82 to 6.72); V2/F = exp(6.26) = 523 L

    lq <- 7.25
    label("Apparent intercompartmental clearance Q (L/h)")  # Table 2, log estimate 7.25 (95% CI 6.94 to 7.56); Q = exp(7.25) = 1408 L/h

    lvp <- 9.97
    label("Apparent peripheral volume of distribution V3/F (L)")  # Table 2, log estimate 9.97 (95% CI 9.40 to 10.5); V3/F = exp(9.97) = 21,402 L

    # ================================================================
    # Inter-individual variability. Ambery 2015 Sect. 3.1: 'Inter-
    # individual variability could be estimated reliably for the
    # apparent elimination clearance' -- CL/F is the only parameter
    # carrying an eta.
    #
    # Table 2's value is the NONMEM $OMEGA VARIANCE, not an SD or a
    # %CV. In the last two rows of Table 2 the column headed 'RSE (%)'
    # actually holds the standard error rather than a percentage
    # (0.137 reproduces the 0.325-0.863 CI as 0.594 +/- 1.96 * 0.137,
    # whereas 0.137% does not), which is the signature of a block
    # pasted verbatim from a NONMEM run summary. The variance reading
    # is confirmed against Fig. 2: in panel f (800 ug once daily) the
    # 95% prediction interval at 20 h post-dose spans roughly 2.3 log
    # units either side of the population mean, which needs
    # var(eta) + var(eps) of order 1 -- reproduced by 0.594 + 0.365
    # (SD ~0.98) and not by the SD reading (0.353 + 0.133, SD ~0.70).
    # ================================================================
    etalcl ~ 0.594  # Table 2, row 'Interindividual variability CL/F' = 0.594 (SE 0.137, 95% CI 0.325 to 0.863); omega^2 on log CL/F, i.e. 90.1% CV

    # ================================================================
    # Residual error. Plasma concentrations were 'expressed as natural
    # logarithms' and modelled on that scale (Sect. 2.4), so the
    # reported proportional error is additive on the natural-log scale
    # -- exactly nlmixr2's lnorm() endpoint. Table 2's 0.365 is the
    # $SIGMA VARIANCE (see the etalcl note above for the evidence),
    # so the SD passed to lnorm() is sqrt(0.365) = 0.6042 log units,
    # equivalent to a 66.4% CV on the linear scale.
    # ================================================================
    expSd <- 0.6042
    label("Residual SD on natural-log plasma concentration (log units)")  # Table 2, row 'Proportional residual error' = 0.365 (SE 0.0175, 95% CI 0.331 to 0.399); sigma^2, so SD = sqrt(0.365) = 0.6042
  })

  model({
    # ----------------------------------------------------------------
    # Individual PK parameters. No covariates were retained (Ambery
    # 2015 Sect. 3.1). Only CL/F carries an eta.
    # ----------------------------------------------------------------
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    # ----------------------------------------------------------------
    # Two-compartment disposition with first-order absorption
    # (Ambery 2015 Sect. 3.1: 'A two-compartment disposition model
    # (ADVAN4 TRANS4) with first-order absorption was adequate to
    # describe the plasma GSK961081 concentration-time data'). NONMEM
    # ADVAN4 TRANS4 parameterises on CL, V2, Q, V3 and KA, which is the
    # parameterisation transcribed above; the micro-constants below are
    # the standard TRANS4 reparameterisation.
    #
    # Relative bioavailability is not identifiable from inhaled-only
    # data, so every disposition parameter is apparent (/F) and no
    # f(depot) term is applied.
    # ----------------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka  * depot
    d/dt(central)     <-  ka  * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ----------------------------------------------------------------
    # Observation. Doses enter in ug and vc is in L, so Cc is in ug/L
    # (= ng/mL); the paper tabulates pg/mL (multiply by 1000).
    # ----------------------------------------------------------------
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
