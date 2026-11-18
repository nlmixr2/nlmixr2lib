Miura_2026_tadalafil <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption for oral",
    "tadalafil in Japanese paediatric (< 15 years) and adult (>= 15 years)",
    "patients with pulmonary arterial hypertension, fitted to routine clinical",
    "(real-world) plasma concentrations. Apparent clearance scales with body",
    "weight (power, centred on the cohort median 14.2 kg) and is 40% lower",
    "with concomitant spironolactone. The absorption rate constant was fixed",
    "to a literature value (0.84 1/h). Residual error is proportional."
  )

  reference <- paste(
    "Miura M, Tanaka S, Aihara R, Yamamoto S, Nakashima Y, Kashiwagura Y,",
    "Hakamata A, Miyakawa S, Tatsuguchi M, Odagiri K, Namiki N, Inui N,",
    "Watanabe H, Uchida S. Development of a population pharmacokinetic model",
    "using combined paediatric and adult data for four pulmonary arterial",
    "hypertension drugs. Br J Clin Pharmacol. 2026;92(10):3611.",
    "doi:10.1002/bcp.70676"
  )

  vignette <- "Miura_2026_pulmonary_arterial_hypertension"

  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "tadalafil", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tadalafil", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on CL/F, (WT / BWmedian)^0.602 (Table 2 footnote).",
        "BWmedian is not printed as a number; 14.2 kg is the median body",
        "weight of the 31 tadalafil patients in Table 1 and is used as the",
        "centring value. Range 5.4-73.5 kg."
      ),
      source_name = "BW"
    ),
    CONMED_SPIRON = list(
      description = "Concomitant spironolactone indicator, 1 = yes",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant spironolactone)",
      notes = paste(
        "Table 2 footnote: CL/F = 1.76 x (BW/BWmedian)^0.602 x (1 - 0.400)^Spironolactone;",
        "for a 0/1 indicator this equals (1 + (-0.400) x CONMED_SPIRON), a 40%",
        "lower CL/F. The authors attribute it to right-heart failure",
        "(congestive hepatopathy) in spironolactone-treated patients rather",
        "than to a drug interaction (Discussion). 13 of 31 tadalafil patients",
        "took spironolactone (Table 1)."
      ),
      source_name = "Spironolactone (SPI)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Candidate covariate in the stepwise search (Methods); not retained."
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      notes = "Candidate covariate in the stepwise search (Methods); not retained."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Candidate covariate in the stepwise search (Methods); not retained for tadalafil."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 31L,
    n_studies = 1L,
    n_observations = 157L,
    age_range = "0.1-74.2 years (paediatric 0.1-8.3, adult 15.5-74.2)",
    age_median = "4.6 years (paediatric 1.9, adult 37.8)",
    weight_range = "5.4-73.5 kg",
    weight_median = "14.2 kg (paediatric 10.7, adult 52.3)",
    sex_female_pct = 61.3,
    race_ethnicity = "Japanese",
    disease_state = "Pulmonary arterial hypertension (PAH)",
    dose_range = "4-40 mg per dose (median 10 mg); paediatric median 9 mg, adult median 40 mg",
    co_medication = paste(
      "Bosentan 4, ambrisentan 15, beraprost 4, furosemide 16, spironolactone 13,",
      "warfarin 8 of 31 patients (Table 1)"
    ),
    regions = "Japan (Hamamatsu University Hospital and Seirei Hamamatsu Hospital)",
    notes = paste(
      "20 paediatric (< 15 years, 62 concentrations) and 11 adult (>= 15 years,",
      "95 concentrations) patients; routine clinical sampling. Demographics",
      "from Table 1 (median (range)). Part of a four-drug analysis of 35",
      "paediatric and 16 adult patients overall."
    )
  )

  ini({
    lka <- fixed(log(0.84)); label("First-order absorption rate constant Ka (1/h)") # Table 2: Ka 0.84 1/h, not estimated (literature value, ref 26)
    lcl <- log(1.76); label("Apparent clearance CL/F at the median weight without spironolactone (L/h)") # Table 2: CL/F 1.76 L/h (95% CI 1.38-2.14)
    lvc <- log(80.6); label("Apparent volume of distribution Vd/F (L)") # Table 2: Vd/F 80.6 L (95% CI 38.3-123)

    e_wt_cl <- 0.602; label("Power exponent of body weight (/14.2 kg) on CL/F (unitless)") # Table 2: Effect of BW 0.602 (95% CI 0.408-0.796)
    e_conmed_spiron_cl <- -0.400; label("Fractional change in CL/F with concomitant spironolactone (unitless)") # Table 2: Effect of SPI -0.400 (95% CI -0.542 to -0.258)

    # Tables 2-3 report IIV as 100 x sqrt(omega^2) (the footnote attributes
    # the undefined lower CI bound to the square root of a negative variance),
    # so omega^2 = (IIV / 100)^2. IIV on Ka is 0 in Table 2, so Ka has no eta.
    etalcl ~ 0.138384 # Table 2: IIV CL/F 37.2 percent -> 0.372^2
    etalvc ~ 0.456976 # Table 2: IIV Vd/F 67.6 percent -> 0.676^2

    propSd <- 0.333; label("Proportional residual error (fraction)") # Table 2: residual error epsilon 1 = 33.3 percent
  })

  model({
    ka <- exp(lka)
    # Table 2 footnote: CL/F = 1.76 x (BW/BWmedian)^0.602 x (1 - 0.400)^Spironolactone
    cl <- exp(lcl + etalcl) * (WT / 14.2)^e_wt_cl * (1 + e_conmed_spiron_cl * CONMED_SPIRON)
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # dose mg / volume L = mg/L; x 1000 for ng/mL
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
