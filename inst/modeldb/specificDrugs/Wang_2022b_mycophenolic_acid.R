Wang_2022b_mycophenolic_acid <- function() {
  description <- "Population PK model for mycophenolic acid (MPA) after oral enteric-coated mycophenolate sodium (EC-MPS) in children after renal transplantation (Wang 2022b). One-compartment disposition with linear elimination and a double extravascular absorption: a first-order arm (rate ka, fraction F of the MPA dose, lag time Tlag1) running in parallel with a zero-order arm (duration Tk02, fraction 1 - F) whose lag time Tlag2 = Tlag1 + diffTlag2 is longer than the first. The EC-MPS dose is converted to MPA by the factor 0.936. Body surface area enters the apparent clearance as a power term (BSA / 1.23)^1.30. Log-normal IIV on all seven absorption and disposition parameters; combined additive plus proportional residual error. Each oral dose is given as two event records, one to depot and one to central with rate = -2."
  reference <- paste(
    "Wang G, Ye Q, Huang Y, Xu H, Li Z.",
    "Population Pharmacokinetics of Enteric-Coated Mycophenolate Sodium in",
    "Children after Renal Transplantation and Initial Dosage Recommendation",
    "Based on Body Surface Area.",
    "Comput Math Methods Med. 2022;2022:1881176.",
    "doi:10.1155/2022/1881176.",
    sep = " "
  )
  vignette <- "Wang_2022b_mycophenolic_acid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "mycophenolic acid", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "mycophenolic acid", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    BSA = list(
      description = "Body surface area (Mosteller formula, sqrt(weight [kg] * height [cm] / 3600)).",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on the apparent clearance, CL/F = 4.28 * (BSA / 1.23)^1.30 L/h (Wang 2022 Results 3.2 final-model equations; Table 2). The Methods covariate template divides each continuous covariate by its mean, so 1.23 m^2 is read as the modelling-cohort mean BSA (Table 1 reports the median 1.26 m^2, range 0.66-1.78 m^2). BSA was computed with the Mosteller formula (Methods 2.1).",
      source_name = "BSA"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 96L,
    n_studies = 1L,
    n_observations = 384L,
    age_range = "4.3-18.0 years",
    age_median = "13.3 years",
    weight_range = "15.0-67.0 kg",
    weight_median = "39.0 kg",
    height_range = "102-171 cm",
    height_median = "147 cm",
    bsa_range = "0.66-1.78 m^2",
    bsa_median = "1.26 m^2",
    sex_female_pct = 45.8,
    race_ethnicity = "Not tabulated; single-centre Chinese paediatric cohort (Children's Hospital of Fudan University, Shanghai).",
    disease_state = "Paediatric renal transplant recipients more than ten days after transplantation (median 396 days, range 13-2258) on a stable oral EC-MPS regimen for at least five days (steady state). Median creatinine clearance 78.8 mL/min/1.73 m^2 (20.1-164.6).",
    dose_range = "Oral enteric-coated mycophenolate sodium (Myfortic, 180 mg tablets); median 10.5 mg/kg/day (3.4-24.0), 325 mg/m^2/day (117-605). Doses are EC-MPS mass; the model converts to MPA with the factor 0.936.",
    regions = "China (Shanghai).",
    co_medication = "Tacrolimus in 94 of 96 patients (97.9%), cyclosporin in 2 (2.08%).",
    notes = "Retrospective TDM data, June 2018 - August 2019 (modelling group, Table 1); MPA sampled 30 min before a dose and 20 min, 1 h and 3 h after it, by EMIT on a Viva-E analyser (LLOQ 0.1 mg/L, BLOQ data retained). Monolix 2019R1. An external validation group of 32 children (September 2019 - July 2020) gave MPE -0.880%, MAPE 10.0% and RMSE 0.810 mg/L (Table 3)."
  )

  ini({
    # Structural parameters: Wang 2022 Table 2, final-model column, and the
    # Results 3.2 final-model equations. CL and V are apparent (CL/F, V/F).
    lka <- log(0.123); label("First-order absorption rate constant ka1 (1/h)") # Table 2: theta_ka1 = 0.123 (RSE 20.5%)
    ld2 <- log(2.90); label("Duration of the zero-order absorption arm Tk02 (h)") # Table 2: theta_Tk02 = 2.90 (RSE 42.1%)
    lffo <- log(0.553); label("Fraction F of the MPA dose absorbed by the first-order arm (fraction)") # Table 2: theta_F = 0.553 (RSE 6.30%)
    ltlag <- log(8.45); label("Lag time of the first-order arm Tlag1 (h)") # Table 2: theta_Tlag1 = 8.45 (RSE 9.30%)
    ldtlag2 <- log(5.78); label("Increment of the zero-order arm lag time over Tlag1, diffTlag2 (h)") # Table 2: theta_diffTlag2 = 5.78 (RSE 18.0%)
    lvc <- log(3.73); label("Apparent volume of distribution V/F (L)") # Table 2: theta_V = 3.73 (RSE 13.0%)
    lcl <- log(4.28); label("Apparent clearance CL/F at BSA = 1.23 m^2 (L/h)") # Table 2: theta_CL = 4.28 (RSE 6.53%)
    e_bsa_cl <- 1.30; label("Power exponent of BSA / 1.23 on CL/F (unitless)") # Table 2: theta_CL-BSA = 1.30 (RSE 20.4%)

    # IIV: Table 2 reports the SD (omega) of each log-normal random effect
    # (Methods 2.2, ln(Pi) = ln(Ppop) + eta_pi with SD omega); the variance
    # below is omega^2.
    etalka ~ 0.783225 # Table 2: omega_ka1 = 0.885 (SD)
    etald2 ~ 5.3361 # Table 2: omega_Tk02 = 2.31 (SD)
    etalffo ~ 0.263169 # Table 2: omega_F1 = 0.513 (SD)
    etaltlag ~ 0.331776 # Table 2: omega_Tlag1 = 0.576 (SD)
    etaldtlag2 ~ 0.885481 # Table 2: omega_diffTlag2 = 0.941 (SD)
    etalvc ~ 0.113569 # Table 2: omega_V = 0.337 (SD)
    etalcl ~ 0.231361 # Table 2: omega_CL = 0.481 (SD)

    # Residual error Cobs = Cpred + sqrt(a^2 + (b * Cpred)^2) * eps
    # (Methods 2.2; Table 2 footnote), i.e. combined error on the SD scale.
    addSd <- 0.0631; label("Additive residual error a (mg/L)") # Table 2: a = 0.0631 (RSE 25.2%)
    propSd <- 0.199; label("Proportional residual error b (fraction)") # Table 2: b = 0.199 (RSE 7.5%)
  })

  model({
    # EC-MPS dose to MPA amount (Methods 2.2, citing de Winter 2009).
    fmpa <- 0.936

    ka <- exp(lka + etalka)
    d2 <- exp(ld2 + etald2)
    # F is log-normal as stated in Methods 2.2, so it can exceed 1 in the
    # upper tail (about 12% of subjects); the zero-order arm 1 - F is then
    # negative. Kept as published; the total absorbed fraction is fmpa
    # for every subject.
    ffo <- exp(lffo + etalffo)
    tlag <- exp(ltlag + etaltlag)
    dtlag2 <- exp(ldtlag2 + etaldtlag2)
    tlag2 <- tlag + dtlag2
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl) * (BSA / 1.23)^e_bsa_cl

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # First-order arm: dose record to depot. Zero-order arm: a parallel
    # dose record to central with rate = -2 so dur(central) applies.
    f(depot) <- fmpa * ffo
    alag(depot) <- tlag
    f(central) <- fmpa * (1 - ffo)
    alag(central) <- tlag2
    dur(central) <- d2

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
