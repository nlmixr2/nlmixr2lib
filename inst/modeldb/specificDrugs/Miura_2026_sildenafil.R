Miura_2026_sildenafil <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption for oral",
    "sildenafil in Japanese paediatric (< 15 years) and adult (>= 15 years)",
    "patients with pulmonary arterial hypertension, fitted to routine clinical",
    "(real-world) plasma concentrations. Apparent clearance scales with body",
    "weight (power, centred on the cohort median 32.3 kg) and is 1.633-fold",
    "higher with concomitant bosentan (a CYP3A4 inducer). The absorption rate",
    "constant was fixed to a literature value (4.51 1/h). Residual error is",
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
    depot = list(analyte = "sildenafil", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "sildenafil", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on CL/F, (WT / BWmedian)^0.537 (Table 2 footnote).",
        "BWmedian is not printed as a number; 32.3 kg is the median body",
        "weight of the 24 sildenafil patients in Table 1 and is used as the",
        "centring value. Range 6.5-58.7 kg."
      ),
      source_name = "BW"
    ),
    CONMED_BOSENTAN = list(
      description = "Concomitant bosentan (CYP3A4 inducer) indicator, 1 = yes",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant bosentan)",
      notes = paste(
        "Table 2 footnote: CL/F = 60.6 x (BW/BWmedian)^0.537 x (1 + 0.633)^Bosentan;",
        "for a 0/1 indicator this equals (1 + 0.633 x CONMED_BOSENTAN), a",
        "1.63-fold higher CL/F (Discussion). 11 of 24 sildenafil patients took",
        "bosentan (Table 1)."
      ),
      source_name = "Bosentan (BOS)"
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
      notes = "Candidate covariate in the stepwise search (Methods); not retained for sildenafil."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 24L,
    n_studies = 1L,
    n_observations = 297L,
    age_range = "0.4-69 years (paediatric 0.4-14.8, adult 38-69)",
    age_median = "10.8 years (paediatric 5.0, adult 47)",
    weight_range = "6.5-58.7 kg",
    weight_median = "32.3 kg (paediatric 16.5, adult 48.5)",
    sex_female_pct = 66.7,
    race_ethnicity = "Japanese",
    disease_state = "Pulmonary arterial hypertension (PAH)",
    dose_range = "1.7-20 mg per dose (median 13.3 mg); paediatric median 12 mg, adult median 20 mg",
    co_medication = paste(
      "Bosentan 11, ambrisentan 5, beraprost 7, furosemide 9, spironolactone 7,",
      "warfarin 4 of 24 patients (Table 1)"
    ),
    regions = "Japan (Hamamatsu University Hospital and Seirei Hamamatsu Hospital)",
    notes = paste(
      "15 paediatric (< 15 years, 37 concentrations) and 9 adult (>= 15 years,",
      "260 concentrations) patients; routine clinical sampling at steady state.",
      "Demographics from Table 1 (median (range)). Part of a four-drug analysis",
      "of 35 paediatric and 16 adult patients overall."
    )
  )

  ini({
    lka <- fixed(log(4.51)); label("First-order absorption rate constant Ka (1/h)") # Table 2: Ka 4.51 1/h, not estimated (literature value, ref 25)
    lcl <- log(60.6); label("Apparent clearance CL/F at the median weight without bosentan (L/h)") # Table 2: CL/F 60.6 L/h (95% CI 44.4-76.8)
    lvc <- log(322); label("Apparent volume of distribution Vd/F (L)") # Table 2: Vd/F 322 L (95% CI 258-387)

    e_wt_cl <- 0.537; label("Power exponent of body weight (/32.3 kg) on CL/F (unitless)") # Table 2: Effect of BW 0.537 (95% CI 0.196-0.878)
    e_conmed_bosentan_cl <- 0.633; label("Fractional increase in CL/F with concomitant bosentan (unitless)") # Table 2: Effect of BOS 0.633 (95% CI 0.455-0.811)

    # Tables 2-3 report IIV as 100 x sqrt(omega^2) (the footnote attributes
    # the undefined lower CI bound to the square root of a negative variance),
    # so omega^2 = (IIV / 100)^2. IIV on Ka is 0 in Table 2, so Ka has no eta.
    etalcl ~ 0.316969 # Table 2: IIV CL/F 56.3 percent -> 0.563^2
    etalvc ~ 0.065025 # Table 2: IIV Vd/F 25.5 percent -> 0.255^2

    propSd <- 0.585; label("Proportional residual error (fraction)") # Table 2: residual error epsilon 1 = 58.5 percent
  })

  model({
    ka <- exp(lka)
    # Table 2 footnote: CL/F = 60.6 x (BW/BWmedian)^0.537 x (1 + 0.633)^Bosentan
    cl <- exp(lcl + etalcl) * (WT / 32.3)^e_wt_cl * (1 + e_conmed_bosentan_cl * CONMED_BOSENTAN)
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # dose mg / volume L = mg/L; x 1000 for ng/mL
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
