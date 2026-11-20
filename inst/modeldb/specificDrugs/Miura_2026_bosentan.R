Miura_2026_bosentan <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption and an",
    "absorption lag time for oral bosentan in Japanese paediatric (< 15",
    "years) and adult (>= 15 years) patients with pulmonary arterial",
    "hypertension, fitted to routine clinical (real-world) plasma",
    "concentrations. Apparent clearance scales with body weight (power,",
    "centred on the cohort median 42.2 kg) and decreases with alanine",
    "aminotransferase (power, centred on the cohort median 16.5 IU/L).",
    "Residual error is proportional."
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
    depot = list(analyte = "bosentan", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "bosentan", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on CL/F, (WT / BWmedian)^0.472 (Table 3 footnote).",
        "BWmedian is not printed as a number; 42.2 kg is the median body",
        "weight of the 16 bosentan patients in Table 1 and is used as the",
        "centring value. Range 6.9-58.7 kg."
      ),
      source_name = "BW"
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on CL/F, (ALT / ALTmedian)^-0.203 (Table 3 footnote).",
        "ALTmedian is not printed as a number; 16.5 is the median ALT of the",
        "16 bosentan patients in Table 1 and is used as the centring value.",
        "Table 1 prints the unit as IU/mL, which cannot be right for ALT",
        "(values 7-50); IU/L is assumed. Range 7-50 IU/L. Whether ALT was",
        "baseline or time-varying is not stated."
      ),
      source_name = "ALT"
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
    )
  )

  population <- list(
    species = "human",
    n_subjects = 16L,
    n_studies = 1L,
    n_observations = 189L,
    age_range = "2.2-69 years (paediatric 2.2-9.2, adult 15.6-69)",
    age_median = "26.8 years (paediatric 4.8, adult 47)",
    weight_range = "6.9-58.7 kg",
    weight_median = "42.2 kg (paediatric 12, adult 48.5)",
    sex_female_pct = 68.8,
    race_ethnicity = "Japanese",
    disease_state = "Pulmonary arterial hypertension (PAH)",
    dose_range = "12-62.5 mg per dose (median 62.5 mg); paediatric median 25 mg, adult 62.5 mg",
    co_medication = paste(
      "Sildenafil 11, tadalafil 4, beraprost 7, furosemide 7, spironolactone 6,",
      "warfarin 2 of 16 patients (Table 1)"
    ),
    regions = "Japan (Hamamatsu University Hospital and Seirei Hamamatsu Hospital)",
    notes = paste(
      "7 paediatric (< 15 years, 19 concentrations) and 9 adult (>= 15 years,",
      "170 concentrations) patients; routine clinical sampling. Demographics",
      "from Table 1 (median (range)); median ALT 16.5 (range 7-50). Part of a",
      "four-drug analysis of 35 paediatric and 16 adult patients overall."
    )
  )

  ini({
    lka <- log(2.47); label("First-order absorption rate constant Ka (1/h)") # Table 3: Ka 2.47 1/h (95% CI 1.25-3.70)
    lcl <- log(8.85); label("Apparent clearance CL/F at the median weight and ALT (L/h)") # Table 3: CL/F 8.85 L/h (95% CI 6.81-10.9)
    lvc <- log(41.0); label("Apparent volume of distribution Vd/F (L)") # Table 3: Vd/F 41.0 L (95% CI 30.5-51.5)
    ltlag <- log(0.379); label("Absorption lag time (h)") # Table 3: Absorption lag time 0.379 h (95% CI 0.309-0.449)

    e_wt_cl <- 0.472; label("Power exponent of body weight (/42.2 kg) on CL/F (unitless)") # Table 3: Effect of BW 0.472 (95% CI 0.192-0.752)
    e_alt_cl <- -0.203; label("Power exponent of ALT (/16.5 IU/L) on CL/F (unitless)") # Table 3: Effect of ALT -0.203 (95% CI -0.341 to -0.065)

    # Tables 2-3 report IIV as 100 x sqrt(omega^2) (the footnote attributes
    # the undefined lower CI bound to the square root of a negative variance),
    # so omega^2 = (IIV / 100)^2. IIV on Ka is 0 in Table 3, so Ka has no eta;
    # Table 3 lists no IIV on the lag time.
    etalcl ~ 0.108241 # Table 3: IIV CL/F 32.9 percent -> 0.329^2
    etalvc ~ 0.088209 # Table 3: IIV Vd/F 29.7 percent -> 0.297^2

    propSd <- 0.590; label("Proportional residual error (fraction)") # Table 3: residual error epsilon 1 = 59.0 percent
  })

  model({
    ka <- exp(lka)
    tlag <- exp(ltlag)
    # Table 3 footnote: CL/F = 8.85 x (BW/BWmedian)^0.472 x (ALT/ALTmedian)^-0.203
    cl <- exp(lcl + etalcl) * (WT / 42.2)^e_wt_cl * (ALT / 16.5)^e_alt_cl
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    alag(depot) <- tlag

    # dose mg / volume L = mg/L; x 1000 for ng/mL
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
