Xu_2024_linezolid <- function() {
  description <- paste(
    "Population PK model for intravenous linezolid in hospital-acquired",
    "pneumonia patients with renal insufficiency (Xu 2024). One-compartment",
    "disposition with linear elimination; clearance carries power-model",
    "effects of age (exponent -0.56, reference 82 years) and Cockcroft-Gault",
    "creatinine clearance (exponent 0.50, reference 42 mL/min). Volume of",
    "distribution is a typical value with no inter-individual variability;",
    "inter-individual variability was estimated on clearance only, and",
    "residual variability is additive."
  )
  reference <- paste(
    "Xu J, Chen X, Zhang Q, Zhuang Z, Yuan Y, Duan L, Shi L, Zhu C, Li J,",
    "Lu J, Yu Y, Tang L. (2024).",
    "Population Pharmacokinetic/Pharmacodynamic Study of Linezolid in",
    "Hospital-Acquired Pneumonia Patients with Renal Insufficiency.",
    "Drug Design, Development and Therapy 18:5073-5086.",
    "doi:10.2147/DDDT.S474470"
  )
  vignette <- "Xu_2024_linezolid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: one entry per d/dt() state. Linezolid was given as a 1 h
  # intravenous infusion and only serum concentrations were measured
  # (Xu 2024, "Dosing Regimen" and "Concentration Determination"), so the
  # single disposition state holds linezolid amount in plasma.
  compartmentData <- list(
    central = list(
      analyte = "linezolid", units = "mg",
      specimen = "plasma", verified = TRUE
    )
  )

  covariateData <- list(
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power-model effect on CL centred at 82 years, per the final-model",
        "equation on Xu 2024 p. 5078: CL = 2.7 x (age/82)^-0.56 x",
        "(CrCL/42.0)^0.50 x exp(eta_CL). The negative exponent makes CL fall",
        "with increasing age. The centering value is exactly the PPK-cohort",
        "median: supplementary Table S1 reports age median (IQR) 82.0",
        "(75.0, 89.0) years for the n = 166 PPK cohort, which independently",
        "corroborates the transcription of the printed equation. Age of 30,",
        "50, 70, 90 and 100 years was simulated in the Monte Carlo dosing",
        "analysis (Table 4)."
      ),
      source_name        = "age"
    ),
    CRCL = list(
      description        = paste(
        "Creatinine clearance estimated by the Cockcroft-Gault formula",
        "(raw mL/min, NOT BSA-normalized)."
      ),
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Carried as the raw Cockcroft-Gault value in mL/min rather than the",
        "canonical BSA-normalized mL/min/1.73 m^2 form - the same deviation",
        "documented in Tsuji_2017_linezolid.R, Delattre_2010_amikacin.R and",
        "Wada_2023_sparsentan.R. Xu 2024 Methods states 'The CrCL value was",
        "estimated using the Cockroft-Gault formula' with no BSA step.",
        "Power-model effect on CL centred at 42.0 mL/min, per the final-model",
        "equation on Xu 2024 p. 5078. As with age, the centering value is",
        "exactly the PPK-cohort median: supplementary Table S1 reports CrCL",
        "median (IQR) 42.0 (29.2, 54.0) mL/min for the n = 166 PPK cohort.",
        "Enrolment required CrCL < 90 mL/min",
        "(renal insufficiency per the 2010 FDA guidance), so the model is",
        "calibrated over roughly 15-90 mL/min; the authors explicitly caution",
        "(Discussion, limitation five) that predictions at CrCL = 90 mL/min",
        "sit at the edge of the observed range. CrCL of 15, 30, 45, 60 and 90",
        "mL/min was simulated in the Monte Carlo dosing analysis (Table 4)."
      ),
      source_name        = "CrCL"
    )
  )

  # Screened during covariate model building but not retained in the final
  # model, so it is documentation only and is deliberately not referenced in
  # model(). Xu 2024 Discussion: "Body weight was not identified as a
  # significant covariate affecting linezolid clearance in our study, which
  # may be attributed to the narrow range of body weights evaluated here."
  covariatesDataExcluded <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Tested as a covariate on clearance and rejected by the stepwise",
        "forward-inclusion / backward-elimination procedure (Xu 2024",
        "Discussion). The authors attribute the null result to the narrow",
        "weight range of the cohort - median (IQR) 60.0 (55.0, 65.0) kg in",
        "the PPK cohort (supplementary Table S1) and identically 60.0",
        "(55.0, 65.0) kg in the pharmacodynamic-study linezolid arm",
        "(Table 1). No allometric term appears in the final model:",
        "supplementary Table S2 tabulates the covariate hypothesis tests and",
        "the full covariate model contains CL-CrCL and CL-age only (basic",
        "model OFV 1327.8; full covariate model 1266.5; backward elimination",
        "removing CL-age +6.9 and removing CL-CrCL +50.7, both P < 0.01), so",
        "weight is absent from the reported model-building path and its",
        "rejection is documented only in the Discussion prose."
      ),
      source_name        = "weight"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 166L,
    n_studies      = 1L,
    age_range      = paste(
      "Adults >= 18 years; PPK cohort (n = 166) median (IQR) 82.0",
      "(75.0, 89.0) years per supplementary Table S1 - the same 82 years",
      "used as the centering value of the final CL equation."
    ),
    age_median     = "82.0 years (supplementary Table S1, n = 166)",
    weight_range   = paste(
      "PPK cohort (n = 166) median (IQR) 60.0 (55.0, 65.0) kg per",
      "supplementary Table S1."
    ),
    weight_median  = "60.0 kg (supplementary Table S1, n = 166)",
    sex_female_pct = 26.5,
    race_ethnicity = "Chinese (single-centre cohort)",
    disease_state  = paste(
      "Adults with hospital-acquired pneumonia diagnosed per the 2016",
      "ATS/IDSA guidelines and renal insufficiency (Cockcroft-Gault",
      "CrCL < 90 mL/min per the 2010 FDA guidance). Patients on renal",
      "replacement therapy, with baseline platelet count < 75000/uL, with",
      "severe hepatic impairment, or with known bleeding disorders were",
      "excluded. More than half had a concomitant non-respiratory infection",
      "site and MRSA was the most frequent respiratory isolate."
    ),
    renal_function = paste(
      "Renal insufficiency by design (CrCL < 90 mL/min); mild to severe.",
      "PPK cohort (n = 166) median (IQR) Cockcroft-Gault CrCL 42.0",
      "(29.2, 54.0) mL/min and serum creatinine 106.2 (70.2, 181.7) umol/L",
      "per supplementary Table S1 - the CrCL median is the same 42.0 mL/min",
      "used as the centering value of the final CL equation."
    ),
    dose_range     = paste(
      "Linezolid 600 mg every 12 h by intravenous infusion (package label",
      "standard dose); clinicians were free to adjust the regimen on the",
      "basis of steady-state trough concentrations. Monte Carlo simulations",
      "explored 600, 400, 300 and 200 mg every 12 h as 1 h infusions."
    ),
    regions        = paste(
      "Single centre: Suzhou Municipal Hospital (The Affiliated Suzhou",
      "Hospital of Nanjing Medical University), Suzhou, China.",
      "July 2018 - August 2023."
    ),
    notes          = paste(
      "Retrospective study. 207 linezolid serum concentrations from 166",
      "patients with renal insufficiency entered the PPK analysis; samples",
      "were predominantly steady-state troughs drawn 30 min or immediately",
      "before the next dose and at least 48 h after therapy started, which",
      "the authors note (limitation four) limits resolution of the",
      "distribution phase. Assay: UPLC-MS/MS (Acquity BEH C18, deuterated",
      "levofloxacin internal standard, MRM 338.6 -> 296.2), linear",
      "0.5-50.0 mg/L. All demographic figures above are the PPK cohort's",
      "own, from supplementary Table S1 (n = 166): 73.5% male, respiratory",
      "failure 47.6%, septic shock 29.5%, MODS 14.5%, albumin 31.8",
      "(28.2, 34.5) g/L, treatment duration 9.0 (8.0, 12.0) days.",
      "Table S1 reports no significant difference on any characteristic",
      "versus the companion pharmacodynamic study's linezolid arm",
      "(n = 81; all P >= 0.105)."
    )
  )

  ini({
    # Structural parameters. Reference age 82 years and reference CrCL
    # 42.0 mL/min are the centering values printed in the final-model
    # equation (Xu 2024 p. 5078). Table 3 reports the same point estimates
    # rounded to one decimal place; the equation is the more precise form
    # for the two covariate exponents and is used here.
    lcl <- log(2.7)
    label("Clearance at age 82 y and CrCL 42.0 mL/min (L/h)")           # Xu 2024 Table 3 (TVCL = 2.7) and final-model equation p. 5078
    lvc <- log(57.1)
    label("Central volume of distribution (L)")                         # Xu 2024 Table 3 (TVv = 57.1) and final-model equation p. 5078

    # Covariate effects on CL, both power-model exponents.
    e_age_cl <- -0.56
    label("Power exponent of age on CL (unitless)")                     # Xu 2024 final-model equation p. 5078 ((age/82)^-0.56); Table 3 theta_Age = -0.6
    e_crcl_cl <- 0.50
    label("Power exponent of creatinine clearance on CL (unitless)")    # Xu 2024 final-model equation p. 5078 ((CrCL/42.0)^0.50); Table 3 theta_CrCL = 0.5

    # Inter-individual variability. Table 3 reports omega^2_CL = 0.1 under
    # the heading "variance of inter-individual variability of clearance"
    # (table footnote), so the value is already on the variance scale that
    # nlmixr2 ini() expects: omega = 0.316, i.e. roughly 32% CV. Eta
    # shrinkage was 23.3%. The paper estimated no IIV on volume of
    # distribution - the printed equation carries exp(eta_CL) on CL only
    # and V is a bare typical value.
    etalcl ~ 0.1                                                        # Xu 2024 Table 3 (omega^2 CL = 0.1, shrinkage 23.3%)

    # Residual variability. Xu 2024 Results: "The additive error model was
    # found to best fit the data", and the Table 3 footnote defines sigma as
    # the "square root of residual variability", i.e. an SD on the
    # concentration scale (mg/L). The stray "(%)" in the Table 3 row label
    # is inconsistent with both statements and is read as a template
    # artefact; see vignette Errata.
    addSd <- 3.4
    label("Additive residual SD (mg/L)")                                # Xu 2024 Table 3 (sigma = 3.4)
  })

  model({
    # Individual clearance: power model in age and Cockcroft-Gault CrCL,
    # centred at the PPK-cohort medians (Xu 2024 p. 5078)
    #   CL (L/h) = 2.7 * (age/82)^-0.56 * (CrCL/42.0)^0.50 * exp(eta_CL)
    cl <- exp(lcl + etalcl) * (AGE / 82)^e_age_cl * (CRCL / 42.0)^e_crcl_cl

    # V (L) = 57.1, a typical value with no inter-individual variability
    vc <- exp(lvc)

    kel <- cl / vc

    # One-compartment disposition with linear elimination; linezolid was
    # administered as a 1 h intravenous infusion directly into `central`
    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
