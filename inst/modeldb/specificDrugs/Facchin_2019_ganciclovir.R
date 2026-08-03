Facchin_2019_ganciclovir <- function() {
  description <- paste(
    "Two-compartment population PK model for ganciclovir after oral valganciclovir",
    "in children with a renal transplant (Facchin 2019), parameterized in apparent",
    "oral clearances and volumes, with first-order absorption after a lag time and",
    "serum-creatinine, body-surface-area and sex effects on CL/F and Vc/F.",
    "Parameters transcribed from the Yang 2023 ganciclovir / valganciclovir",
    "population-PK model repository (Table 3), not from the primary publication;",
    "re-verify against Facchin 2019 when the primary is obtained.",
    sep = " "
  )
  reference <- paste(
    "Facchin A, Elie V, Benyoub N, Magreault S, Maisin A, Storme T, Zhao W,",
    "Deschenes G, Jacqz-Aigrain E. Population pharmacokinetics of ganciclovir after",
    "valganciclovir treatment in children with renal transplant.",
    "Antimicrob Agents Chemother. 2019;63(12):e01192-19. doi:10.1128/AAC.01192-19.",
    "Parameters transcribed from Yang W, Mak W, Gwee A, Gu M, Wu Y, Shi Y, He Q,",
    "Xiang X, Han B, Zhu X. Establishment and Evaluation of a Parametric Population",
    "Pharmacokinetic Model Repository for Ganciclovir and Valganciclovir.",
    "Pharmaceutics. 2023;15(7):1801. doi:10.3390/pharmaceutics15071801 (Table 3).",
    sep = " "
  )
  vignette <- "Yang_2023_ganciclovir_model_repository"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters as a RAW (un-normalized) power term, BSA^1.31 on CL/F and BSA^1.28",
        "on Vc/F, so the structural reference subject is BSA = 1 m^2 rather than a",
        "cohort median or the conventional 1.73 m^2. Body surface area rather than",
        "body weight was the retained body-size descriptor in this model",
        "(Yang 2023 Table 4). Cohort weight median 30.35 kg, range 11.9-83.0 kg;",
        "age median 12.2 years, range 2.1-20.5 (Yang 2023 Table 2).",
        sep = " "
      ),
      source_name        = "BSA"
    ),
    CREAT = list(
      description        = "Serum creatinine concentration",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Yang 2023 Table 3 footnote defines SCR as serum creatinine concentration",
        "in umol/L. Inverse power effect on CL/F referenced to 72.5 umol/L:",
        "(SCR/72.5)^-0.768, so higher serum creatinine (worse renal function)",
        "lowers apparent oral clearance. Serum creatinine, not a derived creatinine",
        "clearance, is the renal-function descriptor retained in this model.",
        sep = " "
      ),
      source_name        = "SCR"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = male",
      notes              = paste(
        "VALUE INVERSION: the source paper codes GENDER = 1 for male and 0 for",
        "female (Yang 2023 Table 3 footnote: 'GENDER: gender, 1 for male and 0 for",
        "female'), and applies multiplicative factors 1.15^GENDER on CL/F and",
        "1.14^GENDER on Vc/F. Stored under the canonical SEXF orientation",
        "(1 = female), so the model applies the factors to (1 - SEXF): males have",
        "15% higher CL/F and 14% higher Vc/F than females of the same body surface",
        "area and serum creatinine. Yang 2023 Section 4.3 notes that although sex",
        "reached statistical significance here, the subsequent effect evaluation",
        "found it was not a clinically meaningful covariate on CL",
        "(less than 20% difference).",
        sep = " "
      ),
      source_name        = "GENDER"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 104L,
    n_studies      = 1L,
    n_observations = 1212L,
    age_median     = "12.2 years (range 2.1-20.5)",
    weight_median  = "30.35 kg (range 11.9-83.0)",
    sex_female_pct = 36.5,
    race_ethnicity = "Not reported.",
    disease_state  = "Children with a renal transplant receiving valganciclovir.",
    dose_range     = "Oral valganciclovir 18.5 mg/kg once or twice daily (range 5.0-70.2 mg/kg).",
    regions        = "France (retrospective study).",
    bioassay       = paste(
      "HPLC. Yang 2023 Table 2 reports the LLOQ as '0.25 mg/mL', which is",
      "implausible for a plasma ganciclovir assay and is most likely a typographic",
      "error for 0.25 mg/L (0.25 ug/mL, matching the LLOQ that the same review",
      "reports for Zhao 2009); flagged for verification against the primary.",
      sep = " "
    ),
    notes          = paste(
      "Demographics and dosing from Yang 2023 Table 2. Intensive sampling at",
      "0 (pre-dose), 0.5, 0.75, 1, 1.5, 2, 6 and 12 h post dose. Covariates tested:",
      "weight, gender, age, height, BSA, underlying disease, serum creatinine,",
      "serum uremia, proteinuria, CrCL, and the co-medications mycophenolate",
      "mofetil, tacrolimus, cyclosporine and azathioprine; BSA, SCR and GENDER were",
      "retained (Yang 2023 Table 4). This study did NOT convert valganciclovir",
      "doses to ganciclovir equivalents; oral doses are supplied as valganciclovir",
      "milligrams and the molecular-weight ratio plus the true bioavailability are",
      "absorbed into the apparent parameters CL/F, Vc/F, Q/F and Vp/F.",
      "The primary study also reported inter-occasion variability (IOV) of 14.4% on",
      "CL, 77.2% on Vp and 111.4% on ka (Yang 2023 Section 3.2.1); per registry",
      "convention IOV is omitted from this packaged model -- see the validation",
      "vignette Assumptions and deviations section.",
      sep = " "
    )
  )

  ini({
    # Structural PK -- Yang 2023 Table 3, Facchin et al. (2019) row. Reference
    # subject: BSA = 1 m^2, SCR = 72.5 umol/L, female (GENDER = 0, i.e. SEXF = 1).
    # Apparent oral clearances in L/h, apparent volumes in L, ka in 1/h, lag in h.
    lcl   <- log(9.07); label("Apparent oral clearance CL/F at BSA = 1 m^2, SCR = 72.5 umol/L, female (L/h)") # Yang 2023 Table 3 (Facchin 2019): CL/F = 9.07 * (SCR/72.5)^-0.768 * BSA^1.31 * 1.15^GENDER
    lvc   <- log(45)  ; label("Apparent central volume Vc/F at BSA = 1 m^2, female (L)")                      # Yang 2023 Table 3 (Facchin 2019): Vc/F = 45 * BSA^1.28 * 1.14^GENDER
    lq    <- log(1.46); label("Apparent inter-compartmental clearance Q/F (L/h)")                             # Yang 2023 Table 3 (Facchin 2019): Q/F = 1.46
    lvp   <- log(18.5); label("Apparent peripheral volume Vp/F (L)")                                          # Yang 2023 Table 3 (Facchin 2019): Vp/F = 18.5
    lka   <- log(6.96); label("First-order oral absorption rate constant (ka, 1/h)")                          # Yang 2023 Table 3 (Facchin 2019): Ka = 6.96
    ltlag <- log(0.86); label("Absorption lag time (Tlag, h)")                                                # Yang 2023 Table 3 (Facchin 2019): Tlag = 0.86

    # Covariate effects. All are non-canonical estimated values.
    e_creat_cl <- -0.768; label("Power exponent of serum creatinine on CL/F (unitless; reference 72.5 umol/L)") # Yang 2023 Table 3 (Facchin 2019): (SCR/72.5)^-0.768
    e_bsa_cl   <- 1.31  ; label("Power exponent of BSA on CL/F (unitless; raw BSA, reference 1 m^2)")           # Yang 2023 Table 3 (Facchin 2019): BSA^1.31
    e_bsa_vc   <- 1.28  ; label("Power exponent of BSA on Vc/F (unitless; raw BSA, reference 1 m^2)")           # Yang 2023 Table 3 (Facchin 2019): BSA^1.28
    e_male_cl  <- 1.15  ; label("Multiplicative male factor on CL/F (unitless; female reference)")              # Yang 2023 Table 3 (Facchin 2019): 1.15^GENDER, GENDER = 1 for male
    e_male_vc  <- 1.14  ; label("Multiplicative male factor on Vc/F (unitless; female reference)")              # Yang 2023 Table 3 (Facchin 2019): 1.14^GENDER, GENDER = 1 for male

    # Between-subject variability. Yang 2023 Methods: %CV = sqrt(omega^2) * 100%,
    # so variance = (BSV% / 100)^2. Q/F and Tlag carry no BSV in the source table.
    etalcl ~ 0.0256    # Yang 2023 Table 3 (Facchin 2019): BSV CL/F = 16.0% -> 0.160^2
    etalvc ~ 0.008649  # Yang 2023 Table 3 (Facchin 2019): BSV Vc/F = 9.3%  -> 0.093^2
    etalvp ~ 0.298116  # Yang 2023 Table 3 (Facchin 2019): BSV Vp/F = 54.6% -> 0.546^2
    etalka ~ 0.350464  # Yang 2023 Table 3 (Facchin 2019): BSV ka   = 59.2% -> 0.592^2

    # Residual unexplained variability: proportional.
    propSd <- 0.235; label("Proportional residual error (fraction)")  # Yang 2023 Table 3 (Facchin 2019): 23.5% proportional error
  })

  model({
    # Sex indicator conversion: the source GENDER column is 1 for male, 0 for
    # female; the canonical SEXF column is 1 for female, so GENDER = 1 - SEXF.
    male <- 1 - SEXF

    cl <- exp(lcl + etalcl) * (CREAT / 72.5)^e_creat_cl * BSA^e_bsa_cl * e_male_cl^male
    vc <- exp(lvc + etalvc) * BSA^e_bsa_vc * e_male_vc^male
    q  <- exp(lq)
    vp <- exp(lvp + etalvp)
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Bioavailability is absorbed into the apparent parameters (CL/F, Vc/F, Q/F,
    # Vp/F), so no explicit f(depot) is applied. Oral valganciclovir doses enter
    # `depot` as valganciclovir milligrams.
    alag(depot) <- tlag

    # Dose in mg, apparent volume in L -> concentration in mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
