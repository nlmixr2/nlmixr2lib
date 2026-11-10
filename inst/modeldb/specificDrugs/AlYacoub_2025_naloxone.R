AlYacoub_2025_naloxone <- function() {
  description <- paste(
    "Two-compartment population PK model for intranasal naloxone (4 mg",
    "Narcan nasal spray, four doses at 0, 45, 300 and 345 min) in adults",
    "with type 1 diabetes mellitus during recurrent exercise (Al Yacoub",
    "2025): sequential zero-order input into the depot over D1 followed by",
    "first-order absorption into the central compartment, linear",
    "elimination, apparent clearances and volumes allometrically scaled by",
    "body weight (fixed exponents 0.75 and 1, reference 70 kg),",
    "between-subject variability on CL/F, V1/F and D1, and combined",
    "proportional and additive residual error. Time unit is minutes.",
    sep = " "
  )
  reference <- paste(
    "Al Yacoub ON, Cheng S, Fayed MS, Fisher J, Brooks J, Seaquist E, Kumar",
    "A, Moheet A, Eberly LE, Coles LD. Intranasal Naloxone During Recurrent",
    "Exercise in Individuals with Type-1 Diabetes Mellitus: Evaluation of",
    "the Clinical Predictors of Pharmacokinetics and Exposure-Response. J",
    "Clin Pharmacol. 2025;65:1577-1586. doi:10.1002/jcph.70067",
    sep = " "
  )
  vignette <- "AlYacoub_2025_naloxone"
  units <- list(time = "min", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "naloxone", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "naloxone", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "naloxone", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling of CL/F and Q/F (exponent 0.75) and of V1/F and",
        "V2/F (exponent 1), reference weight 70 kg, exponents fixed",
        "(Results 'Population PK Analysis'; supplement NONMEM code TVCL =",
        "THETA(2)*(WT/70)**0.75, TVV2 = THETA(3)*(WT/70)). Body weight of",
        "the 27 participants is not tabulated; Table 1 reports BMI 26.5 +/-",
        "4.3 kg/m^2."
      ),
      source_name = "WT"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 27L,
    n_studies = 1L,
    n_observations = 271L,
    age_range = "18-65 years (inclusion criterion)",
    age_mean_sd = "29.8 +/- 9.2 years",
    bmi_mean_sd = "26.5 +/- 4.3 kg/m^2",
    weight_range = "not reported",
    sex_female_pct = 40.7,
    race_ethnicity = "not reported",
    disease_state = paste(
      "Type 1 diabetes mellitus with intact awareness of hypoglycemia;",
      "HbA1c 6.9 +/- 0.8%, diabetes duration 16 +/- 7.9 years (Table 1).",
      "Euglycemia maintained with insulin and variable-rate glucose",
      "infusion during the PK day."
    ),
    dose_range = paste(
      "Four 4 mg intranasal naloxone doses on one day: at the start of a",
      "90-min morning exercise bout (0 min), 45 min later, and again at",
      "the start of and 45 min into an afternoon bout (300 and 345 min)."
    ),
    regions = "United States (single centre, University of Minnesota).",
    notes = paste(
      "Single-centre, single-blinded, placebo-controlled crossover study;",
      "only the naloxone period contributes PK data. Samples at -5, 5, 15,",
      "30, 60, 120, 295, 310, 325, 355, 385 and 415 min. One participant",
      "was excluded as an outlier because the model failed to converge",
      "with their data. 11 female / 16 male."
    )
  )

  ini({
    # Structural parameters -- Table 2 'Final PopPK Model Parameter Estimates'
    # (apparent values, reference weight 70 kg; time in minutes).
    lka <- log(0.0272) ; label("First-order absorption rate constant from depot KA (1/min)") # Table 2 'KA (per min)' = 0.0272; THETA(1)
    lcl <- log(6.82) ; label("Apparent clearance CL/F at 70 kg (L/min)") # Table 2 'CL/F (L/min/70 kg)' = 6.82; THETA(2)
    lvc <- log(171) ; label("Apparent central volume V1/F at 70 kg (L)") # Table 2 'V1/F (L/70 kg)' = 171; THETA(3)
    lq <- log(2.97) ; label("Apparent intercompartmental clearance Q/F at 70 kg (L/min)") # Table 2 'CL2/F (L/min/70 kg)' = 2.97; THETA(4)
    lvp <- log(278) ; label("Apparent peripheral volume V2/F at 70 kg (L)") # Table 2 'V2/F (L/70 kg)' = 278; THETA(5)
    ld1 <- log(7.33) ; label("Duration of zero-order input into the depot D1 (min)") # Table 2 'D1 (min)' = 7.33; THETA(6)

    # Allometric exponents, fixed (Results: 'The exponents for clearance
    # (0.75) and volume of distribution (1.0) were fixed in the model').
    e_wt_cl <- fixed(0.75) ; label("Allometric exponent of body weight on CL/F and Q/F (unitless)") # Results equations; supplement code (WT/70)**0.75
    e_wt_vc <- fixed(1) ; label("Allometric exponent of body weight on V1/F and V2/F (unitless)") # Results equations; supplement code (WT/70)

    # Between-subject variability -- Table 2 'BSV' column holds the omega^2
    # variance with CV% = sqrt(exp(omega^2) - 1) in parentheses (footnote b).
    etalcl ~ 0.103 # Table 2 CL/F BSV = 0.103 (32.9 %CV); ETA(1)
    etalvc ~ 0.268 # Table 2 V1/F BSV = 0.268 (55.4 %CV); ETA(2)
    etald1 ~ 0.108 # Table 2 D1 BSV = 0.108 (33.8 %CV); ETA(3)

    # Residual error. Supplement $ERROR Y = F + EPS(1)*IPRED + EPS(2), so the
    # Table 2 RUV estimates are SIGMA variances; the 'RUV CV%' column is
    # sqrt(sigma^2) (footnote c), which is the SD carried here.
    propSd <- 0.237 ; label("Proportional residual error (fraction)") # Table 2 'RUV% (proportional)' = 5.62 (sigma^2 = 0.0562), 'RUV CV%' = 23.7
    addSd <- 0.348 ; label("Additive residual error (ng/mL)") # Table 2 'RUV (ng/mL)' = 0.121 (sigma^2), 'RUV CV%' = 34.8 = sqrt(0.121)
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q <- exp(lq) * (WT / 70)^e_wt_cl
    vp <- exp(lvp) * (WT / 70)^e_wt_vc
    d1 <- exp(ld1 + etald1)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ADVAN4 TRANS4: each intranasal dose enters the depot as a zero-order
    # input of duration D1 and is then absorbed first-order into the central
    # compartment. Dose records need rate = -2 so rxode2 uses dur(depot).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    dur(depot) <- d1

    # Dose in mg and V1/F in L give mg/L; x1000 gives ng/mL.
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
