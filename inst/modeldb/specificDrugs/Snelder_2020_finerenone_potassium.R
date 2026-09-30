Snelder_2020_finerenone_potassium <- function() {
  description <- paste(
    "Sequential PK/PD indirect-response model for the effect of oral",
    "finerenone on serum potassium in adults with type 2 diabetes and",
    "chronic kidney disease (ARTS-DN and ARTS-DN Japan phase IIb; Snelder",
    "2020). The finerenone PK layer is the final Snelder_2020_finerenone",
    "model carried at fixed values (the PD stage was fitted to individual",
    "post hoc PK). Finerenone plasma concentration inhibits the first-order",
    "potassium loss rate kout through a log-linear model,",
    "EFF = slope * log(Cc + 1), raising serum potassium in proportion to the",
    "baseline. Baseline serum potassium is 3.13% lower in Japanese subjects;",
    "additive (normal) IIV on baseline; proportional residual error whose",
    "variance is 75.2% as large in Japanese subjects."
  )
  reference <- paste(
    "Snelder N, Heinig R, Drenth HJ, Joseph A, Kolkhof P, Lippert J,",
    "Garmann D, Ploeger B, Eissing T. Population Pharmacokinetic and",
    "Exposure-Response Analysis of Finerenone: Insights Based on Phase IIb",
    "Data and Simulations to Support Dose Selection for Pivotal Trials in",
    "Type 2 Diabetes with Chronic Kidney Disease. Clin Pharmacokinet.",
    "2020;59(3):359-370. doi:10.1007/s40262-019-00820-x.",
    "PK values are ESM Table S1 and serum potassium values are ESM Table S3",
    "(final ARTS-DN+JP estimates) of the Electronic Supplementary Material.",
    sep = " "
  )
  vignette <- "Snelder_2020_finerenone"
  units <- list(
    time = "h",
    dosing = "(oral finerenone, mg)",
    concentration = "mmol/L (serum potassium; finerenone plasma concentration Cc is in ug/L)"
  )

  compartmentData <- list(
    depot = list(analyte = "finerenone", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "finerenone", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "finerenone", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "finerenone", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "finerenone", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "finerenone", units = "mg", specimen = "tissue", verified = TRUE),
    serum_k = list(analyte = "potassium", units = "mmol/L", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "PK layer only: Vc/F factor (1 + 0.449 * (log(WT) - log(88.5))), per ESM Table S1.",
      source_name = "BW"
    ),
    CRCL = list(
      description = "eGFR by the MDRD study equation (BSA-normalised estimated glomerular filtration rate)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "PK layer only: CL/F factor (1 + 0.101 * (log(CRCL) - log(63.53))) and F as its reciprocal, per ESM Table S1.",
      source_name = "eGFR-MDRD"
    ),
    RACE_JAPANESE = list(
      description = "Japanese subject (ARTS-DN Japan) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (global ARTS-DN population, which enrolled no Japanese subjects)",
      notes = paste(
        "1 = Japanese subject. Every ARTS-DN Japan subject was Japanese and",
        "ARTS-DN enrolled none (paper Table 2), so the indicator is also the",
        "study indicator. Lowers the typical baseline serum potassium by",
        "3.13% (4.18 mmol/L in Japanese subjects; ESM Table S3 footnote b)",
        "and multiplies the proportional residual variance by 0.752 (ESM",
        "Table S3 footnote c)."
      ),
      source_name = "ethnicity (Japanese vs global)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 893L,
    n_studies = 2L,
    n_observations = 6968L,
    age_range = "5th-95th percentile 49-78 years (ARTS-DN), 44-78 years (ARTS-DN Japan)",
    weight_range = "5th-95th percentile 64.1-126.3 kg (ARTS-DN), 54-100 kg (ARTS-DN Japan)",
    potassium_range = "baseline 5th-95th percentile 3.6-5.0 mmol/L (ARTS-DN, median 4.3), 3.6-4.8 mmol/L (ARTS-DN Japan, median 4.2)",
    race_ethnicity = c(Caucasian = 75.7, Asian_nonJapanese = 9.1, Japanese = 10.5, African_American = 2.9, Other = 1.8),
    disease_state = "type 2 diabetes with persistent albuminuria (UACR >= 30 mg/g) on a renin-angiotensin system blocker (diabetic kidney disease)",
    dose_range = "finerenone 1.25-20 mg once daily orally for 90 days (plus placebo)",
    regions = "global (ARTS-DN) and Japan (ARTS-DN Japan)",
    notes = "PK/PD dataset of 893 subjects with 6968 serum potassium observations (paper Table 2); covariate percentiles from paper Table 3."
  )

  ini({
    # PK layer: final ESM Table S1 estimates, fixed because the PD stage was
    # fitted sequentially to individual post hoc PK (ESM 'Phase IIb PK/PD
    # Model Development')
    lka <- fixed(log(10.7)); label("Common first-order absorption and transit rate constant Ka (1/h)") # ESM Table S1 final model 'Kad, 1/h' = 10.7
    lcl <- fixed(log(37.3)); label("Apparent clearance CL/F at median eGFR-MDRD (L/h)") # ESM Table S1 final model 'CL/F, L/h' = 37.3
    lvc <- fixed(log(123)); label("Apparent central volume Vc/F at median body weight (L)") # ESM Table S1 final model 'Vc/F, L' = 123
    lq <- fixed(log(0.433)); label("Apparent intercompartmental clearance Q/F (L/h)") # ESM Table S1 final model 'Q/F, L/h' = 0.433
    ltlag <- fixed(log(0.215)); label("Absorption lag time (h)") # ESM Table S1 'Lag time, h' = 0.215 (fixed)
    lfdepot <- fixed(log(1)); label("Relative bioavailability at median eGFR-MDRD (fraction)") # ESM Table S1 'F' = 1 (fixed)
    e_wt_vc <- fixed(0.449); label("Slope of Vc/F on log(WT) - log(88.5 kg) (unitless)") # ESM Table S1 final model 'SLVcBW (BW effect)' = 0.449
    e_crcl_cl <- fixed(0.101); label("Slope of CL/F (and inverse F) on log(eGFR-MDRD) - log(63.53) (unitless)") # ESM Table S1 final model 'SLCLeGFR (eGFR effect)' = 0.101

    # Serum potassium indirect-response model: final ESM Table S3 estimates.
    # The baseline carries an additive (normal) IIV, so it is kept on the
    # natural scale.
    rbase <- 4.31; label("Typical baseline serum potassium, global population (mmol/L)") # ESM Table S3 final model 'Baseline serum potassium concentration, mmol/L' = 4.31
    lkin <- log(0.0815); label("Zero-order serum potassium production rate kin (mmol/L/h)") # ESM Table S3 final model 'kin, mmol/L/h' = 0.0815
    lslope <- log(0.0204); label("Log-linear drug-effect slope on kout per log(Cc + 1), Cc in ug/L (unitless)") # ESM Table S3 final model 'SlopeEFF' = 0.0204 (footnote a: log-linear concentration-effect relationship)
    e_japanese_rbase <- -0.0313; label("Fractional change in baseline serum potassium for Japanese subjects (unitless)") # ESM Table S3 final model 'Ethnic effect on baseline, %' = -3.13 (footnote b: baseline = 4.18 mmol/L in Japanese)
    e_japanese_resvar <- 0.752; label("Multiplicative factor on the serum potassium residual variance for Japanese subjects (unitless)") # ESM Table S3 final model 'Ethnic effect on sigma2prop, %' = 75.2 (footnote c: sigma2prop = 0.00325 in Japanese)

    etalka ~ fixed(0.585) # ESM Table S1 final model 'omega2 Ka (IIV)' = 0.585
    etalcl + etalvc ~ fixed(c(0.2, 0.0928, 0.0927)) # ESM Table S1 final model 'omega2 CL/F' = 0.2, 'omega2 CL/F x V/F' = 0.0928, 'omega2 V/F' = 0.0927
    etarbase ~ 0.102 # ESM Table S3 final model 'omega2 baseline (IIV)' = 0.102, additive on the baseline in (mmol/L)^2

    propSd <- sqrt(0.00432); label("Proportional residual error on serum potassium, global population (fraction)") # ESM Table S3 final model 'sigma1^2 prop' = 0.00432; SD = sqrt(0.00432)
  })
  model({
    # PK layer (Snelder_2020_finerenone)
    cov_crcl <- 1 + e_crcl_cl * (log(CRCL) - log(63.53))
    cov_wt <- 1 + e_wt_vc * (log(WT) - log(88.5))

    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * cov_crcl
    vc <- exp(lvc + etalvc) * cov_wt
    vp <- vc
    q <- exp(lq)
    tlag <- exp(ltlag)
    fdepot <- exp(lfdepot) / cov_crcl

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Individual baseline; kout follows from kin so each subject starts at
    # steady state
    rbase_i <- rbase * (1 + e_japanese_rbase * RACE_JAPANESE) + etarbase
    kin <- exp(lkin)
    kout <- kin / rbase_i
    slope <- exp(lslope)

    d/dt(depot) <- -ka * depot
    d/dt(transit1) <- ka * depot - ka * transit1
    d/dt(transit2) <- ka * transit1 - ka * transit2
    d/dt(transit3) <- ka * transit2 - ka * transit3
    d/dt(central) <- ka * transit3 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot
    alag(depot) <- tlag

    Cc <- central / vc * 1000

    # Log-linear inhibition of the loss rate (ESM Eq. S3 and S5)
    eff <- slope * log(Cc + 1)
    d/dt(serum_k) <- kin - kout * (1 - eff) * serum_k
    serum_k(0) <- rbase_i

    propSd_i <- propSd * sqrt(1 + RACE_JAPANESE * (e_japanese_resvar - 1))
    serum_k ~ prop(propSd_i)
  })
}
