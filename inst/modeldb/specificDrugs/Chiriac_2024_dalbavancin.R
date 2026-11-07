Chiriac_2024_dalbavancin <- function() {
  description <- paste(
    "Two-compartment intravenous population PK model for total serum dalbavancin in adult men",
    "with a persistent ventricular assist device (VAD) driveline infection receiving long-term",
    "suppressive outpatient dalbavancin (1500 mg on days 1 and 8, repeated every 6 weeks).",
    "Intercompartmental clearance is fixed to a previously published value; no covariate was",
    "retained.",
    sep = " "
  )
  reference <- paste(
    "Chiriac U, Liebchen U, Frey OR, Lanzinger H, Klein S, Hoppe-Tichy T, Karck M, Meyer A,",
    "Morath B, Monchi M. Model-Based Dose Identification of Dalbavancin for",
    "Long-Term Suppressive Outpatient Treatment of Ventricular Assist Device Infections.",
    "Antibiotics (Basel). 2024;13(11):1103. doi:10.3390/antibiotics13111103",
    sep = " "
  )
  vignette <- "Chiriac_2024_dalbavancin"

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "dalbavancin", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "dalbavancin", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list()

  # Methods 4.2: 'The physiologically plausible covariates BSA, albumin, body
  # weight, age, creatinine, and CrCl were tested'; Results 2.2: 'no covariate
  # relationships for any of the model parameters reached statistical
  # significance'. Supplement Table S3 lists the best candidates (log-weight
  # and log-BSA on V1, log-age on V2), none of which was retained.
  covariatesDataExcluded <- list(
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      notes = "Screened on V1 (Supplement Table S3, log(BSA/mean) power form) but not retained. Table 1 mean 1.95 +/- 0.17 m^2."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = paste(
        "Screened (Methods 4.2) but not retained. Table 1 prints 42.5 +/- 5.7 'mg/dL', which is",
        "physiologically a g/L value (4.25 g/dL); the unit label in the source is a typographical slip."
      )
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened on V1 (Supplement Table S3, log(WT/mean) power form) but not retained. Table 1 median 76 (IQR 20) kg."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on V2 (Supplement Table S3) but not retained. Table 1 mean 57 +/- 9 years."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "(unit not stated)",
      type = "continuous",
      notes = "Screened (Methods 4.2, Supplement Table S2) but not retained; not summarised in Table 1."
    ),
    CRCL = list(
      description = "Creatinine clearance (estimating equation not stated)",
      units = "mL/min",
      type = "continuous",
      notes = paste(
        "Screened (Methods 4.2, Supplement Table S2) but not retained. Table 1 reports",
        "90.8 +/- 34.73 mL/min, not BSA-normalized; the estimating equation is not stated."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 13L,
    n_studies = 1L,
    n_observations = "38 total dalbavancin serum concentrations (Abstract)",
    age_range = "Mean 57 +/- 9 years (Table 1); range not reported",
    weight_range = "Median 76 kg (IQR 20 kg) (Table 1); range not reported",
    sex_female_pct = 0,
    race_ethnicity = "Not reported; single-centre German cohort.",
    disease_state = paste(
      "Adult men with a HeartMate III left ventricular assist device and a persistent or",
      "relapsing driveline infection (S. aureus 12, MRSA 1; 5 with concomitant bloodstream",
      "infection) on long-term suppressive dalbavancin after surgical source control",
      "(Results 2.1). 12/13 non-ischemic cardiomyopathy; 5/13 type 2 diabetes; none on dialysis."
    ),
    renal_function = "Creatinine clearance 90.8 +/- 34.73 mL/min at first administration (Table 1).",
    dose_range = paste(
      "1500 mg intravenously on day 1 and day 8, cycle repeated every 6 weeks; one patient",
      "received 1000 mg because of CrCl < 30 mL/min (Results 2.1). Infusion duration is not stated."
    ),
    regions = "Germany (Heidelberg University Hospital VAD register)",
    notes = paste(
      "Retrospective register study (ethics S-674/2022). Sampling 30 min after the end of the",
      "first infusion, before the day-8 redose and before the day-42 redose (Methods 4.1).",
      "Total serum dalbavancin by HPLC-UV, LOQ 12.5 mg/L; observed protein binding 96-98%.",
      "Estimation by SAEM in Monolix 2023R1; BLQ data replaced with LOQ/2 (Methods 4.2)."
    )
  )

  ini({
    # Structural parameters -- Chiriac 2024 Table 2 (final model). The paper's
    # V1 / V2 map onto the canonical vc / vp.
    lcl <- log(0.050); label("Clearance CL (L/h)") # Chiriac 2024 Table 2: CL = 0.050 L/h (RSE 6.83%)
    lvc <- log(6.5); label("Central volume of distribution V1 (L)") # Chiriac 2024 Table 2: V1 = 6.5 L (RSE 8.05%)
    lq <- fixed(log(0.476)); label("Intercompartmental clearance Q (L/h)") # Chiriac 2024 Table 2: Q = 0.476 L/h, RSE '-' (Results 2.2: fixed to the value of Buckwalter 2005)
    lvp <- log(15.4); label("Peripheral volume of distribution V2 (L)") # Chiriac 2024 Table 2: V2 = 15.4 L (RSE 12.3%)

    # Inter-individual variability. Monolix reports omega as the standard
    # deviation of the log-normal random effect (Methods 4.2: 'Interindividual
    # variability was assumed to be log-normally distributed'), so the variance
    # is omega^2. No IIV on Q and no correlations are reported.
    etalcl ~ 0.0529 # Chiriac 2024 Table 2: omega CL = 0.230 (RSE 22.5%); 0.230^2
    etalvc ~ 0.0676 # Chiriac 2024 Table 2: omega V1 = 0.260 (RSE 28.4%); 0.260^2
    etalvp ~ 0.1681 # Chiriac 2024 Table 2: omega V2 = 0.410 (RSE 24.3%); 0.410^2

    # Residual variability -- Table 2 footnote: 'b: Residual error (proportional)'.
    propSd <- 0.100; label("Proportional residual error (fraction)") # Chiriac 2024 Table 2: b = 0.100 (RSE 54.6%)
  })

  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Intravenous infusion into central; no absorption compartment.
    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Total (bound + unbound) serum dalbavancin, mg/L (dose mg, volume L).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
