Okubo_2021_apremilast <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption and an",
    "absorption lag time for oral apremilast in Japanese and non-Japanese",
    "adults with moderate to severe plaque psoriasis and non-Japanese healthy",
    "adults (Okubo 2021). Apparent clearance depends on psoriasis disease",
    "status, sex, Japanese race and age; apparent central volume depends on",
    "body weight. The individual steady-state AUC over the 12 h dosing",
    "interval drives the companion exposure-response models",
    "Okubo_2021_apremilast_pasi75, Okubo_2021_apremilast_pasi50 and",
    "Okubo_2021_apremilast_spga."
  )
  reference <- paste(
    "Okubo Y, Ohtsuki M, Komine M, Imafuku S, Kassir N, Petric R, Nemoto O.",
    "Population pharmacokinetic and exposure-response analysis of apremilast",
    "in Japanese subjects with moderate to severe psoriasis.",
    "J Dermatol. 2021;48(11):1652-1664. doi:10.1111/1346-8138.16068"
  )
  vignette <- "Okubo_2021_apremilast"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on Vc/F normalised to 82 kg (Okubo 2021 Table 4 footnote e and Results 3.2).",
      source_name = "weight"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL/F normalised to 45 years (Okubo 2021 Table 4 footnote d and Results 3.2).",
      source_name = "age"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (female); the paper's factor 1.25 applies to males",
      notes = paste(
        "Okubo 2021 Appendix S2 codes sex as male = 1, female = 0 and Table 4",
        "prints 'Sex: male 1.25' (CL/F = 1.25 * 9.25 L/h for males, footnote b).",
        "Encoded on the canonical SEXF as a factor e_sexm_cl^(1 - SEXF), so the",
        "typical value 9.25 L/h is the female value."
      ),
      source_name = "sex (male = 1, female = 0)"
    ),
    RACE_JAPANESE = list(
      description = "Japanese race indicator, 1 = Japanese, 0 = non-Japanese",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Japanese: Caucasian, black or African American, or other)",
      notes = "Okubo 2021 Methods 2.4: race recorded as Caucasian, Japanese, black or African American, or other, and tested as Japanese vs non-Japanese.",
      source_name = "race (Japanese vs non-Japanese)"
    ),
    DIS_PSORIASIS = list(
      description = "Plaque psoriasis disease-state indicator, 1 = psoriasis patient, 0 = otherwise",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy subject, or 'other' disease: psoriatic arthritis or rheumatoid arthritis in study PK-010)",
      notes = paste(
        "Okubo 2021 Appendix S2 codes disease as psoriasis = 1, healthy = 0,",
        "other = 2 (psoriatic arthritis or rheumatoid arthritis). Only the",
        "psoriasis effect is retained in the final model (Table 4 row",
        "'Disease status: psoriasis' 0.834, footnote a), so 'other' shares",
        "the healthy reference. Time-fixed per subject."
      ),
      source_name = "disease"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "apremilast", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "apremilast", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 517L,
    n_studies = 9L,
    age_range = "Study means 28.0-52.4 years; PK-024 enrolled healthy adults 65-85 years (Okubo 2021 Tables 1 and 3)",
    weight_range = "Study means 70.2-96.3 kg (Okubo 2021 Table 3)",
    sex_female_pct = 30.6,
    race_ethnicity = c(Japanese = 20.1, Caucasian = 66.7, Black = 10.4, Other = 2.7),
    disease_state = paste(
      "Moderate to severe plaque psoriasis (104 Japanese patients in PSOR-011;",
      "233 non-Japanese patients in PSOR-005 and PSOR-008/ESTEEM 1) and 180",
      "non-Japanese healthy subjects or patients with psoriasis, psoriatic",
      "arthritis or rheumatoid arthritis in six phase 1 studies."
    ),
    dose_range = "Oral apremilast 10, 20 or 30 mg twice daily (phase 2b/3); 20-50 mg single or 30-50 mg twice-daily doses (phase 1).",
    regions = "Japan, USA, Canada, UK, Europe, Australia",
    notes = paste(
      "Okubo 2021 Results 3.2: 5752 samples from 517 subjects, 48 below the",
      "LOQ of 1.00 ng/mL excluded, leaving 5704 (651 from 104 Japanese",
      "psoriasis patients; 5053 from 233 non-Japanese psoriasis patients and",
      "180 healthy subjects). The Table 4 caption prints n = 533. Female",
      "percentage and race percentages are computed from the per-study counts",
      "of Okubo 2021 Table 3 (359 of 517 male; 104 Japanese; 345 Caucasian;",
      "54 black or African American; the remainder other). NONMEM 7.3, FOCE",
      "with interaction (Appendix S2)."
    )
  )

  ini({
    # Structural parameters -- Okubo 2021 Table 4 (final PPK model). Typical
    # values refer to a healthy, female, non-Japanese, 45-year-old, 82-kg subject.
    lka <- log(1.83); label("Absorption rate constant Ka (1/h)") # Table 4: Ka 1.83 1/h (RSE 11.7%)
    lcl <- log(9.25); label("Apparent clearance CL/F for the reference subject (L/h)") # Table 4: CL/F 9.25 L/h (RSE 3.4%)
    lvc <- log(115); label("Apparent central volume Vc/F at 82 kg (L)") # Table 4: Vc/F 115 L (RSE 2.2%)
    ltlag <- log(0.290); label("Absorption lag time (h)") # Table 4: Lag 0.290 h (RSE 11.7%)

    # Covariate effects -- Okubo 2021 Table 4 and Results 3.2
    e_dis_psoriasis_cl <- 0.834; label("Multiplicative factor on CL/F for psoriasis patients vs healthy (unitless)") # Table 4: Disease status psoriasis 0.834 (RSE 4.0%)
    e_sexm_cl <- 1.25; label("Multiplicative factor on CL/F for males vs females (unitless)") # Table 4: Sex male 1.25 (RSE 2.9%)
    e_race_japanese_cl <- 1.17; label("Multiplicative factor on CL/F for Japanese vs non-Japanese (unitless)") # Table 4: Race Japanese 1.17 (RSE 4.3%)
    e_age_cl <- -0.148; label("Power exponent of age (normalised to 45 years) on CL/F (unitless)") # Table 4: Age (age/45)^-0.148 (RSE 26.9%)
    e_wt_vc <- 0.591; label("Power exponent of body weight (normalised to 82 kg) on Vc/F (unitless)") # Table 4: Bodyweight (weight/82)^0.591 (RSE 13.7%)

    # Between-subject variability -- Okubo 2021 Table 4 prints CV%; converted
    # with omega^2 = log(1 + CV^2). The paper estimated a CL/F-Vc/F block but
    # does not print the covariance, so the block is carried as diagonal.
    # The lag time carried a BSV fixed to 0 and is encoded without an eta.
    etalcl ~ 0.1349 # Table 4: BSV CL/F 38.0% (RSE 4.3%); log(1 + 0.380^2)
    etalvc ~ 0.07087 # Table 4: BSV Vc/F 27.1% (RSE 6.1%); log(1 + 0.271^2)
    etalka ~ 0.5280 # Table 4: BSV Ka 83.4% (RSE 8.7%); log(1 + 0.834^2)

    # Residual error -- Okubo 2021 Table 4 'mixed residual error model'
    propSd <- 0.365; label("Proportional residual error (fraction)") # Table 4: Proportional error 36.5% (RSE 2.6%)
    addSd <- 0.658; label("Additive residual error (ng/mL)") # Table 4: Additive error 0.658 ng/mL (RSE 30.0%)
  })

  model({
    # Individual parameters (Okubo 2021 Table 4 footnotes a-e)
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) *
      e_dis_psoriasis_cl^DIS_PSORIASIS *
      e_sexm_cl^(1 - SEXF) *
      e_race_japanese_cl^RACE_JAPANESE *
      (AGE / 45)^e_age_cl
    vc <- exp(lvc + etalvc) * (WT / 82)^e_wt_vc
    tlag <- exp(ltlag)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    alag(depot) <- tlag

    # Dose in mg and volume in L give mg/L; x 1000 gives ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
