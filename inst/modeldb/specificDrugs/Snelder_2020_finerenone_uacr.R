Snelder_2020_finerenone_uacr <- function() {
  description <- paste(
    "Sequential PK/PD indirect-response model for the effect of oral",
    "finerenone on the urinary albumin:creatinine ratio (UACR) in adults with",
    "type 2 diabetes and chronic kidney disease (ARTS-DN and ARTS-DN Japan",
    "phase IIb; Snelder 2020). The finerenone PK layer is the final",
    "Snelder_2020_finerenone model carried at fixed values (the PD stage was",
    "fitted to individual post hoc PK). Finerenone plasma concentration",
    "inhibits the zero-order UACR production rate kin through an Imax model",
    "(Imax fixed to 1, IC50 = 12.9 ug/L); separate typical baselines are",
    "estimated for subjects with screening UACR at or below 300 mg/g and",
    "above 300 mg/g. Log-normal IIV on baseline and log-scale additive",
    "residual error whose variance is 69.2% as large in Japanese subjects."
  )
  reference <- paste(
    "Snelder N, Heinig R, Drenth HJ, Joseph A, Kolkhof P, Lippert J,",
    "Garmann D, Ploeger B, Eissing T. Population Pharmacokinetic and",
    "Exposure-Response Analysis of Finerenone: Insights Based on Phase IIb",
    "Data and Simulations to Support Dose Selection for Pivotal Trials in",
    "Type 2 Diabetes with Chronic Kidney Disease. Clin Pharmacokinet.",
    "2020;59(3):359-370. doi:10.1007/s40262-019-00820-x.",
    "PK values are ESM Table S1 and UACR values are ESM Table S2 (final",
    "ARTS-DN+JP estimates) of the Electronic Supplementary Material.",
    sep = " "
  )
  vignette <- "Snelder_2020_finerenone"
  units <- list(
    time = "h",
    dosing = "(oral finerenone, mg)",
    concentration = "mg/g (UACR; finerenone plasma concentration Cc is in ug/L)"
  )

  compartmentData <- list(
    depot = list(analyte = "finerenone", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "finerenone", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "finerenone", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "finerenone", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "finerenone", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "finerenone", units = "mg", specimen = "tissue", verified = TRUE),
    uacr = list(analyte = "urinary albumin:creatinine ratio", units = "mg/g", specimen = "urine", verified = TRUE)
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
    UACR = list(
      description = "Urinary albumin:creatinine ratio at screening",
      units = "mg/g",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Only used to select the typical UACR baseline: screening UACR",
        "> 300 mg/g (very high albuminuria, KDIGO A3) selects the 633 mg/g",
        "baseline, otherwise the 96.8 mg/g baseline applies (paper Results",
        "3.2.1 and Discussion 4.2; ESM Table S2 'cat3' and 'cat2'). Supply",
        "the screening value, not the modelled time course."
      ),
      source_name = "UACR at screening"
    ),
    RACE_JAPANESE = list(
      description = "Japanese subject (ARTS-DN Japan) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (global ARTS-DN population, which enrolled no Japanese subjects)",
      notes = paste(
        "1 = Japanese subject. Every ARTS-DN Japan subject was Japanese and",
        "ARTS-DN enrolled none (paper Table 2), so the indicator is also the",
        "study indicator. Multiplies the UACR residual variance by 0.692",
        "(ESM Table S2 footnote a)."
      ),
      source_name = "ethnicity (Japanese vs global)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 893L,
    n_studies = 2L,
    n_observations = 5666L,
    age_range = "5th-95th percentile 49-78 years (ARTS-DN), 44-78 years (ARTS-DN Japan)",
    weight_range = "5th-95th percentile 64.1-126.3 kg (ARTS-DN), 54-100 kg (ARTS-DN Japan)",
    uacr_range = "baseline 5th-95th percentile 33.7-1626 mg/g (ARTS-DN, median 192.4), 38.5-1345 mg/g (ARTS-DN Japan, median 216.4)",
    race_ethnicity = c(Caucasian = 75.7, Asian_nonJapanese = 9.1, Japanese = 10.5, African_American = 2.9, Other = 1.8),
    disease_state = "type 2 diabetes with persistent albuminuria (UACR >= 30 mg/g) on a renin-angiotensin system blocker (diabetic kidney disease)",
    dose_range = "finerenone 1.25-20 mg once daily orally for 90 days (plus placebo)",
    regions = "global (ARTS-DN) and Japan (ARTS-DN Japan)",
    notes = "PK/PD dataset of 893 subjects with 5666 UACR observations (paper Table 2); covariate percentiles from paper Table 3."
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

    # UACR indirect-response model: final ESM Table S2 estimates
    lrbase <- log(96.8); label("Typical baseline UACR, screening UACR <= 300 mg/g (mg/g)") # ESM Table S2 final model 'Baseline UACR; cat2, mg/g' = 96.8
    lrbase_uacrgt300 <- log(633); label("Typical baseline UACR, screening UACR > 300 mg/g (mg/g)") # ESM Table S2 final model 'Baseline UACR; cat3, mg/g' = 633
    lkout <- log(0.0014); label("First-order UACR loss rate constant kout (1/h)") # ESM Table S2 final model 'kout, 1/h' = 0.0014
    limax <- fixed(log(1)); label("Maximum fractional inhibition of kin (fraction)") # ESM Table S2 'Imax' = 1 (fixed)
    lic50 <- log(12.9); label("Finerenone concentration giving half-maximal inhibition of kin (ug/L)") # ESM Table S2 final model 'IC50, ug/L' = 12.9
    e_japanese_resvar <- 0.692; label("Multiplicative factor on the UACR residual variance for Japanese subjects (unitless)") # ESM Table S2 final model 'Ethnic effect on sigma2add, %' = 69.2 (footnote a: sigma2add = 0.125 in Japanese)

    etalka ~ fixed(0.585) # ESM Table S1 final model 'omega2 Ka (IIV)' = 0.585
    etalcl + etalvc ~ fixed(c(0.2, 0.0928, 0.0927)) # ESM Table S1 final model 'omega2 CL/F' = 0.2, 'omega2 CL/F x V/F' = 0.0928, 'omega2 V/F' = 0.0927
    etalrbase ~ 0.512 # ESM Table S2 final model 'omega2 baseline (IIV)' = 0.512

    expSd <- sqrt(0.181); label("Additive residual error on log(UACR), global population (SD)") # ESM Table S2 final model 'sigma1^2 add' = 0.181 on log-transformed UACR; SD = sqrt(0.181)
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

    # UACR baseline category: screening UACR > 300 mg/g versus <= 300 mg/g
    uacr_gt300 <- 0
    if (UACR > 300) uacr_gt300 <- 1
    rbase <- (exp(lrbase) * (1 - uacr_gt300) + exp(lrbase_uacrgt300) * uacr_gt300) * exp(etalrbase)
    kout <- exp(lkout)
    kin <- rbase * kout
    imax <- exp(limax)
    ic50 <- exp(lic50)

    d/dt(depot) <- -ka * depot
    d/dt(transit1) <- ka * depot - ka * transit1
    d/dt(transit2) <- ka * transit1 - ka * transit2
    d/dt(transit3) <- ka * transit2 - ka * transit3
    d/dt(central) <- ka * transit3 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot
    alag(depot) <- tlag

    Cc <- central / vc * 1000

    # Imax inhibition of the zero-order production rate (ESM Eq. S2 and S7)
    eff <- imax * Cc / (ic50 + Cc)
    d/dt(uacr) <- kin * (1 - eff) - kout * uacr
    uacr(0) <- rbase

    # UACR was fitted on the log scale; the residual variance is 69.2% as
    # large in Japanese subjects
    expSd_i <- expSd * sqrt(1 + RACE_JAPANESE * (e_japanese_resvar - 1))
    uacr ~ lnorm(expSd_i)
  })
}
