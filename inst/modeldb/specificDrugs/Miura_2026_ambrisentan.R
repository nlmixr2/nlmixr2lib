Miura_2026_ambrisentan <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption for oral",
    "ambrisentan in Japanese paediatric (< 15 years) and adult (>= 15 years)",
    "patients with pulmonary arterial hypertension, fitted to routine clinical",
    "(real-world) plasma concentrations. Apparent clearance scales with body",
    "weight (power, centred on the cohort median 13.4 kg). Residual error is",
    "proportional."
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
    depot = list(analyte = "ambrisentan", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ambrisentan", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on CL/F, (WT / BWmedian)^0.561 (Table 3 footnote).",
        "BWmedian is not printed as a number; 13.4 kg is the median body",
        "weight of the 29 ambrisentan patients in Table 1 and is used as the",
        "centring value. Range 7.1-62 kg."
      ),
      source_name = "BW"
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
      notes = "Candidate covariate in the stepwise search (Methods); not retained for ambrisentan."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 29L,
    n_studies = 1L,
    n_observations = 177L,
    age_range = "0.4-69 years (paediatric 0.4-8.8, adult 29.4-69)",
    age_median = "4.8 years (paediatric 1.8, adult 44.6)",
    weight_range = "7.1-62 kg",
    weight_median = "13.4 kg (paediatric 10.1, adult 51.4)",
    sex_female_pct = 62.1,
    race_ethnicity = "Japanese",
    disease_state = "Pulmonary arterial hypertension (PAH)",
    dose_range = "0.6-10 mg per dose (median 2 mg); paediatric median 1 mg, adult median 10 mg",
    co_medication = paste(
      "Sildenafil 12, tadalafil 13, beraprost 8, furosemide 14, spironolactone 11,",
      "warfarin 6 of 29 patients (Table 1)"
    ),
    regions = "Japan (Hamamatsu University Hospital and Seirei Hamamatsu Hospital)",
    notes = paste(
      "18 paediatric (< 15 years, 66 concentrations) and 11 adult (>= 15 years,",
      "111 concentrations) patients; routine clinical sampling. Demographics",
      "from Table 1 (median (range)). Part of a four-drug analysis of 35",
      "paediatric and 16 adult patients overall."
    )
  )

  ini({
    lka <- log(1.03); label("First-order absorption rate constant Ka (1/h)") # Table 3: Ka 1.03 1/h (95% CI 0.638-1.46)
    lcl <- log(0.594); label("Apparent clearance CL/F at the median weight (L/h)") # Table 3: CL/F 0.594 L/h (95% CI 0.432-0.756)
    lvc <- log(12.4); label("Apparent volume of distribution Vd/F (L)") # Table 3: Vd/F 12.4 L (95% CI 9.26-15.5)

    e_wt_cl <- 0.561; label("Power exponent of body weight (/13.4 kg) on CL/F (unitless)") # Table 3: Effect of BW 0.561 (95% CI 0.247-0.875)

    # Tables 2-3 report IIV as 100 x sqrt(omega^2) (the footnote attributes
    # the undefined lower CI bound to the square root of a negative variance),
    # so omega^2 = (IIV / 100)^2. IIV on Ka is 0 in Table 3, so Ka has no eta.
    etalcl ~ 0.369664 # Table 3: IIV CL/F 60.8 percent -> 0.608^2
    etalvc ~ 0.153664 # Table 3: IIV Vd/F 39.2 percent -> 0.392^2

    propSd <- 0.358; label("Proportional residual error (fraction)") # Table 3: residual error epsilon 1 = 35.8 percent
  })

  model({
    ka <- exp(lka)
    # Table 3 footnote: CL/F = 0.594 x (BW/BWmedian)^0.561
    cl <- exp(lcl + etalcl) * (WT / 13.4)^e_wt_cl
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # dose mg / volume L = mg/L; x 1000 for ng/mL
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
