Snelder_2020_finerenone <- function() {
  description <- paste(
    "Two-compartment population PK model for oral finerenone in adults with",
    "type 2 diabetes and chronic kidney disease (ARTS-DN and ARTS-DN Japan",
    "phase IIb; Snelder 2020). Absorption runs through four sequential",
    "first-order steps at a common rate Ka (depot plus three transit",
    "compartments) after a fixed 0.215 h lag time; the peripheral volume is",
    "fixed equal to the central volume. Body weight scales Vc/F and",
    "eGFR-MDRD scales CL/F, both as linear functions of the log-covariate",
    "centred at the pooled median; the eGFR effect also enters the relative",
    "bioavailability inversely, so lower eGFR raises exposure through both",
    "clearance and bioavailability. IIV on Ka and a correlated CL/F-Vc/F",
    "block; proportional residual error."
  )
  reference <- paste(
    "Snelder N, Heinig R, Drenth HJ, Joseph A, Kolkhof P, Lippert J,",
    "Garmann D, Ploeger B, Eissing T. Population Pharmacokinetic and",
    "Exposure-Response Analysis of Finerenone: Insights Based on Phase IIb",
    "Data and Simulations to Support Dose Selection for Pivotal Trials in",
    "Type 2 Diabetes with Chronic Kidney Disease. Clin Pharmacokinet.",
    "2020;59(3):359-370. doi:10.1007/s40262-019-00820-x.",
    "Parameter values are the final (ARTS-DN+JP) estimates of Table S1 in",
    "the Electronic Supplementary Material.",
    sep = " "
  )
  vignette <- "Snelder_2020_finerenone"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  compartmentData <- list(
    depot = list(analyte = "finerenone", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "finerenone", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "finerenone", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "finerenone", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "finerenone", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "finerenone", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters Vc/F (and hence Vp/F, which equals Vc/F) as",
        "(1 + 0.449 * (log(WT) - log(88.5))), where 88.5 kg is the median of",
        "the pooled ARTS-DN + ARTS-DN Japan dataset (ESM Table S1 footnote).",
        "Baseline versus time-varying use is not stated in the source."
      ),
      source_name = "BW"
    ),
    CRCL = list(
      description = "eGFR by the MDRD study equation (BSA-normalised estimated glomerular filtration rate)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "eGFR-MDRD (ESM Eq. S9), used as the PK covariate for consistency",
        "with earlier finerenone popPK analyses even though ARTS-DN measured",
        "eGFR-EPI for the PD endpoints. Enters CL/F as",
        "(1 + 0.101 * (log(CRCL) - log(63.53))) and F as the reciprocal of",
        "the same factor; 63.53 mL/min/1.73 m^2 is the pooled median",
        "(ESM Table S1 footnote). Baseline versus time-varying use is not",
        "stated in the source."
      ),
      source_name = "eGFR-MDRD"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 787L,
    n_studies = 2L,
    n_observations = 4597L,
    age_range = "5th-95th percentile 49-78 years (ARTS-DN), 44-78 years (ARTS-DN Japan); medians 65 and 64 years",
    weight_range = "5th-95th percentile 64.1-126.3 kg (ARTS-DN, median 90.6), 54-100 kg (ARTS-DN Japan, median 71.6)",
    egfr_mdrd_range = "5th-95th percentile 33.5-102.5 mL/min/1.73 m^2 (ARTS-DN, median 63.9), 41.0-87.2 (ARTS-DN Japan, median 61.5)",
    race_ethnicity = c(Caucasian = 75.7, Asian_nonJapanese = 9.1, Japanese = 10.5, African_American = 2.9, Other = 1.8),
    disease_state = "type 2 diabetes with persistent albuminuria (UACR >= 30 mg/g) on a renin-angiotensin system blocker (diabetic kidney disease)",
    dose_range = "finerenone 1.25, 2.5, 5, 7.5, 10, 15 or 20 mg once daily orally for 90 days (plus placebo)",
    regions = "global (ARTS-DN, no sites in Japan) and Japan (ARTS-DN Japan)",
    notes = paste(
      "PK dataset: 705 ARTS-DN plus 82 ARTS-DN Japan subjects with 4597",
      "quantifiable observations (paper Table 2); 607 observations (11.3%)",
      "below the 0.1 ug/L LLOQ were excluded. Covariate percentiles are the",
      "PK/PD analysis dataset (paper Table 3; n = 893). Ethnicity",
      "percentages are the PK/PD dataset (paper Table 2). No ethnicity",
      "effect on PK was retained."
    )
  )

  ini({
    lka <- log(10.7); label("Common first-order absorption and transit rate constant Ka (1/h)") # ESM Table S1 final model 'Kad, 1/h' = 10.7
    lcl <- log(37.3); label("Apparent clearance CL/F at median eGFR-MDRD (L/h)") # ESM Table S1 final model 'CL/F, L/h' = 37.3
    lvc <- log(123); label("Apparent central volume Vc/F at median body weight (L)") # ESM Table S1 final model 'Vc/F, L' = 123
    lq <- log(0.433); label("Apparent intercompartmental clearance Q/F (L/h)") # ESM Table S1 final model 'Q/F, L/h' = 0.433
    ltlag <- fixed(log(0.215)); label("Absorption lag time (h)") # ESM Table S1 'Lag time, h' = 0.215 (fixed)
    lfdepot <- fixed(log(1)); label("Relative bioavailability at median eGFR-MDRD (fraction)") # ESM Table S1 'F' = 1 (fixed)

    e_wt_vc <- 0.449; label("Slope of Vc/F on log(WT) - log(88.5 kg) (unitless)") # ESM Table S1 final model 'SLVcBW (BW effect)' = 0.449
    e_crcl_cl <- 0.101; label("Slope of CL/F (and inverse F) on log(eGFR-MDRD) - log(63.53) (unitless)") # ESM Table S1 final model 'SLCLeGFR (eGFR effect)' = 0.101

    etalka ~ 0.585 # ESM Table S1 final model 'omega2 Ka (IIV)' = 0.585
    etalcl + etalvc ~ c(0.2, 0.0928, 0.0927) # ESM Table S1 final model 'omega2 CL/F' = 0.2, 'omega2 CL/F x V/F' = 0.0928, 'omega2 V/F' = 0.0927

    propSd <- sqrt(0.179); label("Proportional residual error (fraction)") # ESM Table S1 final model 'sigma1^2 prop' = 0.179; SD = sqrt(0.179)
  })
  model({
    # Covariate factors, linear in the log-covariate (ESM Table S1 formulae);
    # medians of the pooled ARTS-DN + ARTS-DN Japan dataset
    cov_crcl <- 1 + e_crcl_cl * (log(CRCL) - log(63.53))
    cov_wt <- 1 + e_wt_vc * (log(WT) - log(88.5))

    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * cov_crcl
    vc <- exp(lvc + etalvc) * cov_wt
    # Peripheral volume assumed equal to the central volume (ESM 'Phase IIa
    # Models, Pharmacokinetics'; paper Discussion 4.1)
    vp <- vc
    q <- exp(lq)
    tlag <- exp(ltlag)
    # F = F / (1 + SLCLeGFR * (ln(eGFR) - ln(median eGFR))) (ESM Table S1)
    fdepot <- exp(lfdepot) / cov_crcl

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # The same Ka applies from the dose compartment through three transit
    # compartments into central (ESM Table S1 footnote d)
    d/dt(depot) <- -ka * depot
    d/dt(transit1) <- ka * depot - ka * transit1
    d/dt(transit2) <- ka * transit1 - ka * transit2
    d/dt(transit3) <- ka * transit2 - ka * transit3
    d/dt(central) <- ka * transit3 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot
    alag(depot) <- tlag

    # Dose in mg and volume in L give mg/L; finerenone is reported in ug/L
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
